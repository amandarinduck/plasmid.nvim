-- Feature impact of a sequence edit. Pure Lua.
--
-- build(rec, op) reads op.changes (filled in by model/transaction.lua when a
-- sequence op is applied) and classifies each feature whose bases changed:
--   destroyed  all bases deleted
--   truncated  fewer bases
--   extended   more bases
--   modified   same length, different bases (substitution or partial inversion)
-- For CDS (and ORF) features the old and new translations are compared:
--   frameshift (length change not a multiple of 3), premature stop, lost start
--   or stop codon, and the number of amino-acid changes (0 = silent).
local translate = require("plasmid.bio.translate")

local M = {}

local CODING = { CDS = true, ORF = true }

local function translation(f, res, default_table)
  local cs = tonumber(f:get("codon_start")) or 1
  local t = tonumber(f:get("transl_table")) or default_table or 11
  return translate.translate(res:sub(cs), { table = t, cds_start = true }), translate.get_table(t), res:sub(cs, cs + 2)
end

local function coding(entry, f, old_res, new_res, default_table)
  local oldp, t, old_start = translation(f, old_res, default_table)
  local newp, _, new_start = translation(f, new_res, default_table)
  local dlen = #new_res - #old_res
  entry.frameshift = dlen % 3 ~= 0
  entry.lost_start = translate.is_start(t, old_start) and not translate.is_start(t, new_start)
  local old_stop = oldp:sub(-1) == "*"
  entry.lost_stop = old_stop and newp:sub(-1) ~= "*" and not entry.frameshift
  -- first stop before the last codon that was not there before
  local body_new = old_stop and newp:sub(1, -2) or newp
  local body_old = old_stop and oldp:sub(1, -2) or oldp
  local s_new = body_new:find("*", 1, true)
  if s_new and not body_old:find("*", 1, true) then entry.premature_stop = s_new end
  if not entry.frameshift then
    local diff = 0
    if #oldp == #newp then
      for i = 1, #oldp do if oldp:sub(i, i) ~= newp:sub(i, i) then diff = diff + 1 end end
    else
      diff = math.abs(#newp - #oldp)
    end
    entry.aa_changes = diff
  end
  entry.old_aa, entry.new_aa = #oldp, #newp
end

--- @return { entries: table[], moved: integer, cut: boolean }
function M.build(rec, op)
  local report = { entries = {}, moved = op.moved or 0, cut = false }
  for _, c in ipairs(op.changes or {}) do
    local f = c.new or c.old
    -- `source` spans the whole record; every edit "changes" it
    if f.type == "source" then goto continue end
    local e = {
      id = f.id, type = f.type, label = f:label(), old_len = #(c.old_res or ""),
      new_len = c.new and #(c.new_res or "") or 0, inverted = c.inverted,
    }
    if not c.new then e.kind = "destroyed"
    elseif e.new_len < e.old_len then e.kind = "truncated"
    elseif e.new_len > e.old_len then e.kind = "extended"
    else e.kind = "modified" end
    if e.kind == "destroyed" or e.kind == "truncated" then report.cut = true end
    if CODING[f.type] and c.new and (c.old_res or "") ~= "" then
      coding(e, f, c.old_res, c.new_res, rec and (rec:is_circular() and 11 or 1) or 11)
    end
    report.entries[#report.entries + 1] = e
    ::continue::
  end
  return report
end

local function signed(n) return (n >= 0 and "+" or "−") .. math.abs(n) end

--- One line per entry, e.g. "tet  CDS  extended +3 bp, frame kept, 1 aa inserted".
function M.describe_entry(e, unit)
  unit = unit or "bp"
  local parts = {}
  if e.kind == "destroyed" then
    parts[1] = ("destroyed (%d %s removed)"):format(e.old_len, unit)
  elseif e.kind == "modified" then
    parts[1] = e.inverted and "partly inverted" or "bases changed"
  else
    parts[1] = ("%s %s %s"):format(e.kind, signed(e.new_len - e.old_len), unit)
  end
  if e.frameshift then parts[#parts + 1] = "FRAMESHIFT" end
  if e.premature_stop then parts[#parts + 1] = ("stop codon at aa %d"):format(e.premature_stop) end
  if e.lost_start then parts[#parts + 1] = "start codon lost" end
  if e.lost_stop then parts[#parts + 1] = "stop codon lost" end
  if e.translation == "updated" then parts[#parts + 1] = "/translation updated"
  elseif e.translation == "removed" then parts[#parts + 1] = "/translation removed" end
  if e.aa_changes and not e.frameshift then
    if e.aa_changes == 0 then parts[#parts + 1] = "silent"
    elseif e.new_aa ~= e.old_aa then parts[#parts + 1] = ("in frame, %s aa"):format(signed(e.new_aa - e.old_aa))
    else parts[#parts + 1] = ("%d aa changed"):format(e.aa_changes) end
  end
  -- unnamed features fall back to their type as label: show it once
  local name = (e.label == e.type or e.label:gsub(" ", "_") == e.type) and e.type or (e.label .. "  " .. e.type)
  return ("%s  %s"):format(name, table.concat(parts, ", "))
end

local function is_alarming(e)
  return e.kind == "destroyed" or e.kind == "truncated" or e.frameshift or e.premature_stop or e.lost_start or e.lost_stop
end

-- Most important entry: problems first, then by feature-type priority.
local function lead(report)
  local ok, cfg = pcall(function() return require("plasmid.config").options.priority end)
  local types = ok and cfg and cfg.types or {}
  local best, score
  for _, e in ipairs(report.entries) do
    local sc = (is_alarming(e) and 10000 or 0) + (types[e.type] or 40)
    if not best or sc > score then best, score = e, sc end
  end
  return best
end

--- Short summary for a notification.
function M.summary(report, unit)
  local n = #report.entries
  if n == 0 then
    if report.moved > 0 then return ("%d feature%s shifted"):format(report.moved, report.moved == 1 and "" or "s") end
    return "no features affected"
  end
  local first = M.describe_entry(lead(report), unit)
  if n == 1 then return first end
  return ("%s; +%d more (:DnaEditReport)"):format(first, n - 1)
end

--- Is there something the user must look at (cut feature or broken CDS)?
function M.alarming(report)
  for _, e in ipairs(report.entries) do
    if e.kind == "destroyed" or e.kind == "truncated" or e.frameshift or e.premature_stop or e.lost_start or e.lost_stop then
      return true
    end
  end
  return false
end

return M
