-- Feature detection in the editor: run a search, show proposals as an overlay,
-- review and accept them as real features.
local config = require("plasmid.config")
local detect = require("plasmid.bio.detect")
local featuredb = require("plasmid.featuredb")
local record = require("plasmid.model.record")
local util = require("plasmid.util")

local M = {}

local ID_BASE = -1000000 -- overlay ids for proposals (ORFs use -1, -2, ...)

local function buffer() return require("plasmid.buffer") end

local function location_text(h, L)
  local inner = h.wraps and ("join(%d..%d,1..%d)"):format(h.start, L, h.stop)
    or (h.start == h.stop and tostring(h.start) or ("%d..%d"):format(h.start, h.stop))
  return h.strand == -1 and ("complement(" .. inner .. ")") or inner
end

function M.describe(h)
  return ("%s/%s, %.1f%% identity, %.0f%% of the reference%s"):format(
    h.ref.db or "?", h.ref.id, 100 * h.identity, 100 * h.coverage, h.fragment and " (fragment)" or "")
end

--- Is this hit already annotated? Same name on an overlapping span, or the
--- same type on (nearly) the same span.
local function already_annotated(rec, h)
  local L = rec:len()
  for _, f in ipairs(rec.features) do
    if f.loc then
      for _, e in ipairs(rec:extents(f)) do
        local spans = h.wraps and { { h.start, L }, { 1, h.stop } } or { { h.start, h.stop } }
        for _, s in ipairs(spans) do
          local lo, hi = math.max(e[1], s[1]), math.min(e[2], s[2])
          local ov = hi - lo + 1
          if ov > 0 then
            if f:label():lower() == h.ref.name:lower() and ov >= 0.5 * h.length then return true end
            if f.type == h.ref.type and math.abs(e[1] - s[1]) <= 5 and math.abs(e[2] - s[2]) <= 5 then return true end
          end
        end
      end
    end
  end
  return false
end

--- Overlay features for current proposals (called when overlays are rebuilt).
function M.overlay_features(st, rec)
  if not st.detect or st.detect.rec ~= rec then return {} end
  local out = {}
  for i, h in ipairs(st.detect.hits) do
    if not already_annotated(rec, h) then
      local f = record.Feature.new({
        type = h.ref.type,
        location_text = location_text(h, rec:len()),
        qualifiers = { { key = "label", value = h.ref.name .. (h.fragment and " (fragment)" or ""), quoted = true } },
      })
      f.id = ID_BASE - i
      f.overlay = true
      f.detected = h
      out[#out + 1] = f
    end
  end
  return out
end

--- Run detection on the current record.
function M.run(buf, names, opts)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  if rec:is_protein() then return util.notify("plasmid: detection works on DNA/RNA records", vim.log.levels.WARN) end
  local refs, used = featuredb.references(names)
  if #refs == 0 then
    return util.notify("plasmid: no feature databases yet; :DnaFeatureDbImport {path} (see :help plasmid-detect)",
      vim.log.levels.WARN)
  end
  local t0 = vim.uv.hrtime()
  local o = vim.tbl_extend("force", config.options.detect or {}, opts or {}, { circular = rec:is_circular() })
  local hits = detect.search(rec.seq:tostring(), refs, o)
  st.detect = { rec = rec, hits = hits, dbs = used }
  require("plasmid.orfs").refresh(buf)
  buffer().render(buf)
  local shown = 0
  for _, f in ipairs(rec.overlays or {}) do if f.detected then shown = shown + 1 end end
  util.notify(("plasmid: %d proposed feature%s (%d already annotated) from %s in %.1fs; ga accepts, gO lists")
    :format(shown, shown == 1 and "" or "s", #hits - shown, table.concat(used, ", "), (vim.uv.hrtime() - t0) / 1e9))
  return hits
end

function M.clear(buf)
  local st = buffer().get(buf)
  if not st then return end
  st.detect = nil
  require("plasmid.orfs").refresh(buf)
  buffer().render(buf)
end

--- Accept one proposal (overlay feature id) as a real feature.
function M.accept(buf, fid, batch)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  local o = rec:feature_by_id(fid)
  if not (o and o.detected) then return nil end
  local h = o.detected
  -- keep the "(fragment)" marker the proposal showed: a partial match should
  -- not later read as the full-length element
  local quals = { { key = "label", value = o:label(), quoted = true } }
  if h.ref.blurb and h.ref.blurb ~= "" then quals[#quals + 1] = { key = "note", value = h.ref.blurb, quoted = true } end
  quals[#quals + 1] = { key = "note", value = "detected by plasmid.nvim: " .. M.describe(h), quoted = true }
  local f = record.Feature.new({ type = h.ref.type, location_text = o.location_text, qualifiers = quals })
  f.pinned = true
  if batch then return f end
  require("plasmid.author").commit(buf, { op = "add", feature = f }, "accept " .. h.ref.name)
  util.notify(("plasmid: added %s %s"):format(h.ref.type, h.ref.name))
  return f
end

--- Accept every full-length proposal in one undoable step.
function M.accept_all(buf, include_fragments)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  local ops, left = {}, 0
  for _, f in ipairs(rec.overlays or {}) do
    if f.detected and (include_fragments or not f.detected.fragment) then
      ops[#ops + 1] = { op = "add", feature = M.accept(buf, f.id, true) }
    elseif f.detected then
      left = left + 1
    end
  end
  local tail = left > 0 and ("; %d fragment%s left to review (:DnaDetectAccept! fragments takes them)")
    :format(left, left == 1 and "" or "s") or ""
  if #ops == 0 then return util.notify("plasmid: nothing to accept" .. tail) end
  require("plasmid.author").commit(buf, { op = "batch", ops = ops }, ("accept %d detected features"):format(#ops))
  util.notify(("plasmid: added %d detected features (u undoes all)%s"):format(#ops, tail))
end

--- Picker over current proposals.
function M.review(buf)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  local items = {}
  for _, f in ipairs(rec.overlays or {}) do if f.detected then items[#items + 1] = f end end
  if #items == 0 then return util.notify("plasmid: no proposals; run :DnaDetect") end
  vim.ui.select(items, {
    prompt = "Detected features (accept)",
    format_item = function(f)
      local h = f.detected
      return ("%-24s %-13s %6.1f%%  %4.0f%%  %s  %s"):format(h.ref.name:sub(1, 24), h.ref.type,
        100 * h.identity, 100 * h.coverage, f.location_text, h.ref.db or "")
    end,
  }, function(choice)
    if choice then
      require("plasmid.ui").goto_feature(buf, choice.id)
      M.accept(buf, choice.id)
    end
  end)
end

return M
