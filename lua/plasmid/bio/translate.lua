-- Codon translation with NCBI genetic codes. Pure Lua.
local tables = require("plasmid.bio.codon_tables")
local alphabet = require("plasmid.model.alphabet")

local M = {}

local idx = { T = 0, C = 1, A = 2, G = 3, U = 0 }

local function codon_index(c1, c2, c3)
  local a, b, c = idx[c1], idx[c2], idx[c3]
  if not (a and b and c) then return nil end
  return a * 16 + b * 4 + c + 1
end

local cache = {}

-- Ambiguous codons resolve to an amino acid only when every expansion agrees
-- (same rule as Biopython); otherwise "X".
local function ambiguous(t, c1, c2, c3)
  local e1, e2, e3 = alphabet.iupac[c1], alphabet.iupac[c2], alphabet.iupac[c3]
  if not (e1 and e2 and e3) then return "X" end
  local seen
  for i = 1, #e1 do
    for j = 1, #e2 do
      for k = 1, #e3 do
        local ci = codon_index(e1:sub(i, i), e2:sub(j, j), e3:sub(k, k))
        local aa = t.aa:sub(ci, ci)
        if seen == nil then
          seen = aa
        elseif seen ~= aa then
          return "X"
        end
      end
    end
  end
  return seen or "X"
end

---@param table_id integer|nil NCBI transl_table (default 1)
function M.get_table(table_id)
  local id = tonumber(table_id) or 1
  local t = cache[id]
  if t then return t end
  local raw = tables[id] or tables[1]
  t = { id = id, name = raw.name, aa = raw.aa, starts = {} }
  for _, s in ipairs(raw.starts) do t.starts[s] = true end
  cache[id] = t
  return t
end

--- Translate one codon (3 chars, any case).
function M.codon(t, codon)
  local up = codon:upper()
  local c1, c2, c3 = up:sub(1, 1), up:sub(2, 2), up:sub(3, 3)
  local i = codon_index(c1, c2, c3)
  if i then return t.aa:sub(i, i) end
  return ambiguous(t, c1, c2, c3)
end

function M.is_start(t, codon)
  return t.starts[codon:upper():gsub("U", "T")] == true
end

--- Translate a coding sequence string.
---@param seq string
---@param opts { table?: integer, cds_start?: boolean, to_stop?: boolean, partial?: boolean }|nil
function M.translate(seq, opts)
  opts = opts or {}
  local t = M.get_table(opts.table)
  local out = {}
  local stopped = false
  for i = 1, #seq - 2, 3 do
    local cod = seq:sub(i, i + 2)
    local aa = M.codon(t, cod)
    if i == 1 and opts.cds_start and M.is_start(t, cod) then aa = "M" end
    if aa == "*" and opts.to_stop then stopped = true break end
    out[#out + 1] = aa
  end
  -- a trailing partial codon (3'-incomplete CDS) counts when its first bases
  -- already fix the amino acid, e.g. "TC" = Ser whatever the third base is
  -- (NCBI's /translation includes it)
  local rest = #seq % 3
  if rest > 0 and not stopped and opts.partial ~= false then
    local aa = M.codon(t, seq:sub(#seq - rest + 1) .. ("N"):rep(3 - rest))
    if aa ~= "X" and aa ~= "*" then out[#out + 1] = aa end
  end
  return table.concat(out)
end

return M
