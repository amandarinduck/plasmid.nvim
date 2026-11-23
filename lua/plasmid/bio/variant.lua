-- What alignment differences do to the features they fall in.
--
-- The differences come from bio/align.differences: substitutions, insertions
-- and deletions in reference coordinates. Here each one is placed in the
-- features it touches, and coding features are re-translated from the
-- observed sequence so a substitution can be called silent, missense or
-- nonsense and an indel in frame or a frameshift.
local alphabet = require("plasmid.model.alphabet")
local translate = require("plasmid.bio.translate")
local impact = require("plasmid.bio.impact")

local M = {}

local CODING = { CDS = true, ORF = true }

local function default_table(rec, f)
  return tonumber(f:get("transl_table")) or (rec:is_circular() and 11 or 1)
end

--- Reference positions of a feature, 5'->3' along the feature's own strand.
--- Returns the list and a lookup from position to its index.
local function feature_positions(rec, f)
  local L = rec:len()
  local pos = {}
  for _, e in ipairs(rec:extents(f)) do
    local a, b = e[1], e[2]
    if a <= b then
      for p = a, b do pos[#pos + 1] = p end
    else
      for p = a, L do pos[#pos + 1] = p end
      for p = 1, b do pos[#pos + 1] = p end
    end
  end
  if f.strand == -1 then
    local rev = {}
    for i = #pos, 1, -1 do rev[#rev + 1] = pos[i] end
    pos = rev
  end
  local index = {}
  for i, p in ipairs(pos) do index[p] = i end
  return pos, index
end

--- The feature's residues with the differences applied, in the feature's own
--- orientation. Returns the observed residues and how many bases changed.
local function observed(rec, f, diffs)
  local seq = rec.seq:tostring()
  local pos, index = feature_positions(rec, f)
  local minus = f.strand == -1
  local sub_at, ins_after, deleted = {}, {}, {}
  local touched = 0
  for _, d in ipairs(diffs) do
    if d.kind == "sub" then
      for k = 0, d.len - 1 do
        if index[d.pos + k] then sub_at[d.pos + k] = d.query:sub(k + 1, k + 1); touched = touched + 1 end
      end
    elseif d.kind == "del" then
      for k = 0, d.len - 1 do
        if index[d.pos + k] then deleted[d.pos + k] = true; touched = touched + 1 end
      end
    else
      if index[d.pos] then ins_after[d.pos] = (ins_after[d.pos] or "") .. d.query; touched = touched + #d.query end
    end
  end
  if touched == 0 then return nil end
  local out = {}
  -- walk the feature's bases in reference order so insertions land correctly,
  -- then flip once for a minus-strand feature
  local ordered = pos
  if minus then
    ordered = {}
    for i = #pos, 1, -1 do ordered[#ordered + 1] = pos[i] end
  end
  for _, p in ipairs(ordered) do
    if not deleted[p] then out[#out + 1] = sub_at[p] or seq:sub(p, p) end
    if ins_after[p] then out[#out + 1] = ins_after[p] end
  end
  local res = table.concat(out)
  if minus then res = alphabet.revcomp(res) end
  return res, touched
end

--- Features that any difference falls in, with the effect on each.
--- Returns a list of
--- { feature, id, type, label, diffs, kind, aa_changes, frameshift, ... }.
function M.by_feature(rec, diffs)
  local out = {}
  for _, f in ipairs(rec.features) do
    if f.type ~= "source" and f.loc then
      local _, index = feature_positions(rec, f)
      local mine = {}
      for _, d in ipairs(diffs) do
        local hit = false
        if d.kind == "ins" then hit = index[d.pos] ~= nil
        else
          for k = 0, d.len - 1 do if index[d.pos + k] then hit = true break end end
        end
        if hit then mine[#mine + 1] = d end
      end
      if #mine > 0 then
        local entry = { feature = f, id = f.id, type = f.type, label = f:label(), diffs = mine }
        local old_res = rec:extract(f)
        local new_res = observed(rec, f, mine)
        entry.old_res, entry.new_res = old_res, new_res
        if CODING[f.type] and new_res then
          impact.coding(entry, f, old_res, new_res, default_table(rec, f))
          if entry.frameshift then entry.kind = "frameshift"
          elseif entry.premature_stop then entry.kind = "nonsense"
          elseif (entry.aa_changes or 0) > 0 then entry.kind = "missense"
          elseif #old_res ~= #new_res then entry.kind = "in-frame indel"
          else entry.kind = "silent" end
        else
          entry.kind = (#old_res ~= #(new_res or old_res)) and "length change" or "changed"
        end
        out[#out + 1] = entry
      end
    end
  end
  return out
end

--- Per-difference annotation: which features it falls in and, for a
--- substitution inside a CDS, the codon and amino-acid change.
--- Each difference gains `features` (labels) and, where it applies,
--- `codon`, `aa_from`, `aa_to` and `effect`.
function M.annotate(rec, diffs)
  for _, d in ipairs(diffs) do
    d.features = {}
    for _, f in ipairs(rec.features) do
      if f.type ~= "source" and f.loc then
        local pos, index = feature_positions(rec, f)
        local hit = false
        if d.kind == "ins" then hit = index[d.pos] ~= nil
        else
          for k = 0, d.len - 1 do if index[d.pos + k] then hit = true break end end
        end
        if hit then
          d.features[#d.features + 1] = f:label()
          if CODING[f.type] and d.kind == "sub" and d.len == 1 then
            local cs = (tonumber(f:get("codon_start")) or 1) - 1
            local off = index[d.pos] - cs            -- 1-based along the coding frame
            if off >= 1 then
              local codon_no = math.floor((off - 1) / 3) + 1
              local within = (off - 1) % 3 + 1
              local res = rec:extract(f):sub(cs + 1)
              local from = res:sub((codon_no - 1) * 3 + 1, codon_no * 3)
              if #from == 3 then
                local base = d.query
                if f.strand == -1 then base = alphabet.revcomp(base) end
                local to = from:sub(1, within - 1) .. base .. from:sub(within + 1)
                local t = translate.get_table(default_table(rec, f))
                local aa_from = translate.codon(t, from:upper())
                local aa_to = translate.codon(t, to:upper())
                d.codon, d.codon_from, d.codon_to = codon_no, from:upper(), to:upper()
                d.aa_from, d.aa_to = aa_from, aa_to
                d.effect = (aa_from == aa_to) and "silent"
                  or (aa_to == "*" and "nonsense" or "missense")
                d.feature = f:label()
              end
            end
          end
        end
      end
    end
  end
  return diffs
end

--- One line per difference, for a report: "1734 A>G  tet  codon 21 GCT>GCC silent".
function M.describe(d, unit)
  unit = unit or "bp"
  local where
  if d.kind == "sub" then
    where = ("%d %s>%s"):format(d.pos, d.ref, d.query)
  elseif d.kind == "del" then
    where = ("%d del %d %s (%s)"):format(d.pos, d.len, unit, d.ref:sub(1, 12))
  else
    where = ("%d ins %d %s (%s)"):format(d.pos, d.len, unit, d.query:sub(1, 12))
  end
  local tail = ""
  if d.codon then
    tail = ("  %s codon %d %s>%s %s%s"):format(d.feature or "", d.codon, d.codon_from, d.codon_to,
      d.aa_from == d.aa_to and "" or (d.aa_from .. tostring(d.codon) .. d.aa_to .. " "), d.effect)
  elseif d.features and #d.features > 0 then
    tail = "  " .. table.concat(d.features, ", ")
  end
  return where .. tail
end

return M
