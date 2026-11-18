-- Motif search. Pure Lua.
--
-- Nucleotide motifs use IUPAC codes on both strands, across the origin of
-- circular records, with an optional number of mismatches. A sequence base
-- matches a motif position when every base it may stand for is allowed there
-- (so an N in the sequence only matches N, R matches R/D/V/N, ...).
-- Amino-acid motifs are searched in protein records, or in the six-frame
-- translation of a nucleotide record (hits mapped back to bases); X matches
-- any residue, B = D/N, Z = E/Q, J = I/L.
local M = {}

local alphabet = require("plasmid.model.alphabet")
local bit = require("bit")

-- A=1 C=2 G=4 T=8
local NUC = {
  A = 1, C = 2, G = 4, T = 8, U = 8,
  R = 5, Y = 10, S = 6, W = 9, K = 12, M = 3,
  B = 14, D = 13, H = 11, V = 7, N = 15,
}
M.NUC_MASK = NUC

local AA_SETS = { B = "DN", Z = "EQ", J = "IL" }

--- Normalise a motif: uppercase, spaces dropped. Returns motif, kind
--- ("nuc" | "aa") or nil, error. An "aa:" prefix forces amino acids.
function M.parse_motif(text, protein_record)
  local s = (text or ""):gsub("%s", "")
  local forced = s:match("^[Aa][Aa]:(.*)$")
  if forced then s = forced end
  s = s:upper()
  if s == "" then return nil, "empty motif" end
  if protein_record or forced then
    if s:find("[^ACDEFGHIKLMNPQRSTVWYBZJXUO%*]") then return nil, "not an amino-acid motif: " .. text end
    return s, "aa"
  end
  for i = 1, #s do
    if not NUC[s:sub(i, i)] then
      return nil, ("'%s' is not an IUPAC nucleotide code (use aa: for a protein motif)"):format(s:sub(i, i))
    end
  end
  return s, "nuc"
end

local function masks(s)
  local t = {}
  for i = 1, #s do t[i] = NUC[s:sub(i, i)] or 0 end
  return t
end

-- Lua pattern for an exact IUPAC motif (matched against an uppercase sequence).
local CLASS = {}
for code, m in pairs(NUC) do
  local letters = {}
  for _, b in ipairs({ { "A", 1 }, { "C", 2 }, { "G", 4 }, { "T", 8 } }) do
    if bit.band(m, b[2]) ~= 0 then letters[#letters + 1] = b[1] end
  end
  if code == "U" then letters = { "T", "U" } end
  CLASS[code] = letters
end

-- Sequence letters allowed at a motif position: those whose mask is a subset.
local function allowed_letters(code)
  local m = NUC[code]
  local out = {}
  for letter, lm in pairs(NUC) do
    if bit.band(lm, bit.bnot(m)) == 0 then out[#out + 1] = letter end
  end
  table.sort(out)
  return out
end

local function lua_pattern(motif)
  local parts = {}
  for i = 1, #motif do
    local c = motif:sub(i, i)
    local al = allowed_letters(c)
    parts[i] = #al == 1 and al[1] or ("[" .. table.concat(al) .. "]")
  end
  return table.concat(parts)
end

--- IUPAC reverse complement of a motif.
M.lua_pattern = function(m) return lua_pattern(m) end

function M.revcomp_motif(motif)
  return alphabet.revcomp(motif):upper()
end

-- Positions (1-based starts) where motif matches `hay` with at most k
-- mismatches, considering starts 1..last.
local function scan(hay, motif, k, last)
  local out = {}
  local m = #motif
  if k == 0 then
    local pat = lua_pattern(motif)
    local init = 1
    while true do
      local s = hay:find(pat, init)
      if not s or s > last then break end
      out[#out + 1] = { s, 0 }
      init = s + 1
    end
    return out
  end
  local pm = masks(motif)
  local hb = {}
  for i = 1, #hay do hb[i] = NUC[hay:sub(i, i)] or 0 end
  local band, bnot = bit.band, bit.bnot
  for s = 1, last do
    local mm = 0
    for j = 1, m do
      local b = hb[s + j - 1]
      if b == 0 or band(b, bnot(pm[j])) ~= 0 then
        mm = mm + 1
        if mm > k then break end
      end
    end
    if mm <= k then out[#out + 1] = { s, mm } end
  end
  return out
end

--- Find a nucleotide motif. seq: string; opts { circular, mismatches, strands
--- = "both"|"top"|"bottom" }. Returns hits sorted by start:
--- { start, stop, strand = 1|-1, mismatches } (stop < start: wraps the origin).
--- A palindromic motif is reported once per site, as strand 1.
function M.find_nuc(seq, motif, opts)
  opts = opts or {}
  local k = opts.mismatches or 0
  local L, m = #seq, #motif
  if m == 0 or m > L then return {} end
  -- opts.hay: the sequence already upper-cased (and, when circular, extended by
  -- at least m-1 bases) so that many motifs can share one copy.
  local hay = opts.hay
  if not hay then
    hay = seq:upper()
    if opts.circular then hay = hay .. hay:sub(1, m - 1) end
  end
  local last = opts.circular and L or (L - m + 1)
  local hits = {}
  local strands = opts.strands or "both"
  local rc = M.revcomp_motif(motif)
  local seen = {}
  local function add(list, strand)
    for _, h in ipairs(list) do
      local s = h[1]
      local key = s
      if strand == -1 and rc == motif then goto continue end
      local e = s + m - 1
      if e > L then e = e - L end
      hits[#hits + 1] = { start = s, stop = e, strand = strand, mismatches = h[2] }
      seen[key] = true
      ::continue::
    end
  end
  if strands ~= "bottom" then add(scan(hay, motif, k, last), 1) end
  if strands ~= "top" then add(scan(hay, rc, k, last), -1) end
  table.sort(hits, function(a, b)
    if a.start ~= b.start then return a.start < b.start end
    return a.strand > b.strand
  end)
  return hits
end

---------------------------------------------------------------------------
-- Amino-acid motifs
---------------------------------------------------------------------------

local function aa_match(c, p)
  if p == "X" or c == p then return true end
  local set = AA_SETS[p]
  return set ~= nil and set:find(c, 1, true) ~= nil
end

local function scan_aa(prot, motif, k)
  local out = {}
  local m = #motif
  for s = 1, #prot - m + 1 do
    local mm = 0
    for j = 1, m do
      if not aa_match(prot:sub(s + j - 1, s + j - 1), motif:sub(j, j)) then
        mm = mm + 1
        if mm > k then break end
      end
    end
    if mm <= k then out[#out + 1] = { s, mm } end
  end
  return out
end

--- Find an amino-acid motif in a protein sequence: hits { start, stop, strand = 1 }.
function M.find_protein(seq, motif, opts)
  local hits = {}
  for _, h in ipairs(scan_aa(seq:upper(), motif, (opts or {}).mismatches or 0)) do
    hits[#hits + 1] = { start = h[1], stop = h[1] + #motif - 1, strand = 1, mismatches = h[2] }
  end
  return hits
end

--- Find an amino-acid motif in the six-frame translation of a nucleotide
--- sequence (genetic code opts.table, default 11). Hits are in bases:
--- { start, stop, strand, frame = 1..3, mismatches } (stop < start: wraps).
function M.find_translated(seq, motif, opts)
  opts = opts or {}
  local translate = require("plasmid.bio.translate")
  local t = translate.get_table(opts.table or 11)
  local L = #seq
  local up = seq:upper()
  local hits = {}
  local nt = #motif * 3
  for strand = 1, -1, -2 do
    local s = strand == 1 and up or alphabet.revcomp(up)
    -- circular: extend so codons crossing the origin are translated
    local ext = opts.circular and (s .. s:sub(1, nt + 2)) or s
    for frame = 1, 3 do
      local aa = {}
      for i = frame, #ext - 2, 3 do aa[#aa + 1] = translate.codon(t, ext:sub(i, i + 2)) end
      local prot = table.concat(aa)
      for _, h in ipairs(scan_aa(prot, motif, opts.mismatches or 0)) do
        local a = frame + (h[1] - 1) * 3            -- start in this strand's coordinates
        if a <= L then                             -- skip duplicates from the extension
          local b = a + nt - 1
          local start, stop
          if strand == 1 then
            start, stop = a, b
          else
            start, stop = L - b + 1, L - a + 1
          end
          if start < 1 then start = start + L end
          if stop > L then stop = stop - L end
          if stop < 1 then stop = stop + L end
          if start > L then start = start - L end
          hits[#hits + 1] = { start = start, stop = stop, strand = strand, frame = frame, mismatches = h[2] }
        end
      end
    end
  end
  table.sort(hits, function(x, y) return x.start < y.start end)
  return hits
end

--- Search a record: dispatches on motif kind. Returns hits, kind or nil, err.
function M.find(rec, text, mismatches)
  local motif, kind = M.parse_motif(text, rec:is_protein())
  if not motif then return nil, kind end
  local seq = rec.seq:tostring()
  local opts = { circular = rec:is_circular(), mismatches = mismatches or 0 }
  if mismatches and mismatches >= #motif then return nil, "more mismatches than motif letters" end
  if rec:is_protein() then return M.find_protein(seq, motif, opts), "aa", motif end
  if kind == "aa" then
    opts.table = tonumber(rec.transl_table) or (rec:is_circular() and 11 or 1)
    return M.find_translated(seq, motif, opts), "aa", motif
  end
  return M.find_nuc(seq, motif, opts), "nuc", motif
end

return M
