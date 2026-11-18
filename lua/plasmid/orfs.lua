-- ORF overlay: transient pseudo-features (negative ids) drawn alongside real
-- features, never written to disk unless promoted to a CDS.
local config = require("plasmid.config")
local orf = require("plasmid.bio.orf")
local record = require("plasmid.model.record")
local translate = require("plasmid.bio.translate")
local util = require("plasmid.util")

local M = {}

local function buffer() return require("plasmid.buffer") end

--- Genetic code for a record: configured, else the most common /transl_table
--- among its CDS, else 11 for circular molecules (plasmids) and 1 otherwise.
function M.table_for(rec)
  local o = config.options.orfs
  if o.table then return o.table end
  local counts, best, best_n = {}, nil, 0
  for _, f in ipairs(rec.features) do
    local t = f.type == "CDS" and tonumber(f:get("transl_table"))
    if t then
      counts[t] = (counts[t] or 0) + 1
      if counts[t] > best_n then best, best_n = t, counts[t] end
    end
  end
  return best or (rec:is_circular() and 11 or 1)
end

--- Build overlay features for a record.
function M.compute(rec, opts)
  local t = M.table_for(rec)
  local found = orf.find(rec.seq:tostring(), {
    circular = rec:is_circular(), min_aa = opts.min_aa, starts = opts.starts, table = t,
  })
  -- stop codons already annotated as CDS 3' ends
  local annotated = {}
  if opts.hide_annotated then
    for _, f in ipairs(rec.features) do
      if f.type == "CDS" and f.strand ~= 0 then
        local ext = rec:extents(f)
        local three
        if #ext == 2 then
          local head, tail
          for _, e in ipairs(ext) do if e[1] == 1 then head = e else tail = e end end
          three = f.strand == 1 and head[2] or tail[1]
        elseif #ext == 1 then
          three = f.strand == 1 and ext[1][2] or ext[1][1]
        end
        if three then annotated[f.strand .. ":" .. three] = true end
      end
    end
  end
  local feats, by_id = {}, {}
  local n = 0
  for _, o in ipairs(found) do
    local three = o.strand == 1 and o.stop or o.start
    if not annotated[o.strand .. ":" .. three] then
      n = n + 1
      local f = record.Feature.new({
        type = "ORF",
        location_text = orf.location(o, rec:len()),
        qualifiers = {
          { key = "label", value = ("ORF %d aa"):format(o.aa_len), quoted = true },
          { key = "codon_start", value = "1", quoted = false },
          { key = "transl_table", value = tostring(t), quoted = false },
        },
      })
      f.id = -n
      f.overlay = true
      f.orf = o
      feats[#feats + 1] = f
      by_id[f.id] = f
    end
  end
  return feats, by_id
end

--- Recompute the overlay for the buffer's current record if ORFs are shown.
function M.refresh(buf)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  for _, r in ipairs(st.records) do
    if r ~= rec then r.overlays, r.overlay_by_id = nil, nil end
  end
  local feats, by_id = {}, {}
  if st.orfs then feats, by_id = M.compute(rec, st.orfs) end
  for _, f in ipairs(require("plasmid.detectui").overlay_features(st, rec)) do
    feats[#feats + 1] = f
    by_id[f.id] = f
  end
  if #feats > 0 then
    rec.overlays, rec.overlay_by_id = feats, by_id
  else
    rec.overlays, rec.overlay_by_id = nil, nil
  end
end

--- Show/hide ORFs. `min_aa` given -> show with that minimum.
function M.toggle(buf, min_aa)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  if rec:is_protein() then return util.notify("plasmid: proteins have no reading frames", vim.log.levels.WARN) end
  if st.orfs and not min_aa then
    st.orfs = nil
  else
    st.orfs = vim.tbl_extend("force", config.options.orfs, { min_aa = min_aa or config.options.orfs.min_aa })
  end
  M.refresh(buf)
  buffer().render(buf)
  if st.orfs then
    util.notify(("plasmid: %d ORFs ≥ %d aa (%s starts, code %d)%s"):format(#(rec.overlays or {}), st.orfs.min_aa,
      st.orfs.starts == "atg" and "ATG" or "alternative", M.table_for(rec),
      st.orfs.hide_annotated and ", annotated CDS hidden" or ""))
  end
end

--- Promote an ORF overlay to a real CDS feature (undoable).
function M.promote(buf, fid, label)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  local o = rec.overlay_by_id and rec.overlay_by_id[fid]
  if not o then return util.notify("plasmid: no ORF under cursor", vim.log.levels.WARN) end
  local function go(name)
    if not name or name == "" then return end
    local t = tonumber(o:get("transl_table")) or 11
    local nt = rec:extract(o)
    local aa = translate.translate(nt, { table = t, cds_start = true }):gsub("%*$", "")
    local f = record.Feature.new({
      type = "CDS", location_text = o.location_text,
      qualifiers = {
        { key = "label", value = name, quoted = true },
        { key = "codon_start", value = "1", quoted = false },
        { key = "transl_table", value = tostring(t), quoted = false },
        { key = "translation", value = aa, quoted = true },
      },
    })
    f.pinned = true
    require("plasmid.author").commit(buf, { op = "add", feature = f }, "add CDS " .. name)
    util.notify(("plasmid: added CDS %s %s (%d aa)"):format(name, f.location_text, #aa))
  end
  if label then return go(label) end
  vim.ui.input({ prompt = "CDS label: ", default = o:label() }, go)
end

function M.promote_under_cursor()
  local buf = vim.api.nvim_get_current_buf()
  local _, fid = buffer().cursor_info(vim.api.nvim_get_current_win())
  if not fid or fid > 0 then return util.notify("plasmid: put the cursor on an ORF (go shows them)", vim.log.levels.WARN) end
  M.promote(buf, fid)
end

return M
