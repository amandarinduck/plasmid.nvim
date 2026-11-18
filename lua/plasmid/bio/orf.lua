-- Six-frame open reading frame finder. Pure Lua.
--
-- An ORF runs from the first start codon after the previous in-frame stop
-- (longest ORF per stop) up to and including the stop codon. ORFs that run off
-- the end of a linear sequence without a stop are not reported.
--
-- Circular records: each frame is scanned over three copies of the sequence
-- and ORFs that start in the middle copy are kept, so every ORF sees its full
-- upstream context and ORFs crossing the origin come out naturally, whatever
-- the length modulo 3.
local alphabet = require("plasmid.model.alphabet")
local translate = require("plasmid.bio.translate")

local M = {}

-- byte -> 0..3 for T C A G (U as T); others -> nil
local NUC = {}
for ch, v in pairs({ T = 0, C = 1, A = 2, G = 3, U = 0 }) do
  NUC[ch:byte()] = v
  NUC[ch:lower():byte()] = v
end

--- Per-table codon classes: class[index] = "*" stop, "M" start, true other.
local function classes(table_id, start_mode)
  local t = translate.get_table(table_id)
  local cls = {}
  local bases = { "T", "C", "A", "G" }
  for a = 0, 3 do
    for b = 0, 3 do
      for c = 0, 3 do
        local codon = bases[a + 1] .. bases[b + 1] .. bases[c + 1]
        local i = a * 16 + b * 4 + c
        local aa = t.aa:sub(i + 1, i + 1)
        if aa == "*" then
          cls[i] = "*"
        elseif (start_mode == "atg" and codon == "ATG") or (start_mode ~= "atg" and t.starts[codon]) then
          cls[i] = "M"
        else
          cls[i] = true
        end
      end
    end
  end
  return cls
end

-- Scan one strand string `s` (positions 1..#s) in all frames.
-- Calls emit(start_pos, stop_end_pos) for each ORF (1-based, inclusive, stop codon included).
local function scan(s, cls, min_codons, emit)
  local n = #s
  local byte = string.byte
  for frame = 0, 2 do
    local orf_start = nil
    local i = 1 + frame
    while i + 2 <= n do
      local a, b, c = NUC[byte(s, i)], NUC[byte(s, i + 1)], NUC[byte(s, i + 2)]
      local k = (a and b and c) and cls[a * 16 + b * 4 + c] or true
      if k == "*" then
        if orf_start and (i - orf_start) / 3 >= min_codons then emit(orf_start, i + 2) end
        orf_start = nil
      elseif k == "M" and not orf_start then
        orf_start = i
      end
      i = i + 3
    end
  end
end

--- Find ORFs.
---@param seq string residues
---@param opts { circular?: boolean, min_aa?: integer, table?: integer, starts?: "atg"|"alternative" }
---@return table[] orfs { start, stop, strand, aa_len, wraps } on the forward axis; for
---   wrapping ORFs start > stop. Sorted by start, then strand.
function M.find(seq, opts)
  opts = opts or {}
  local L = #seq
  local min_aa = opts.min_aa or 75
  local cls = classes(opts.table, opts.starts or "atg")
  local circ = opts.circular and L >= 3
  local out, seen = {}, {}

  local function add(strand, s1, e1)
    -- s1, e1: 1-based inclusive on the scanned strand; map into 1..L on the forward axis
    local len = e1 - s1 + 1
    if len > L then return end -- a frame with no stop all the way round
    local a, b
    if strand == 1 then
      a, b = (s1 - 1) % L + 1, (e1 - 1) % L + 1
    else
      -- position p on the reverse complement corresponds to L - p + 1 (mod L)
      a, b = (L - e1) % L + 1, (L - s1) % L + 1
    end
    local key = strand .. ":" .. a .. ":" .. b
    if seen[key] then return end
    seen[key] = true
    out[#out + 1] = { start = a, stop = b, strand = strand, aa_len = len / 3 - 1, wraps = a > b }
  end

  local rc = alphabet.revcomp(seq)
  for _, strand in ipairs({ 1, -1 }) do
    local s = strand == 1 and seq or rc
    if circ then
      scan(s .. s .. s, cls, min_aa, function(a, b)
        if a > L and a <= 2 * L then add(strand, a, b) end
      end)
    else
      scan(s, cls, min_aa, function(a, b) add(strand, a, b) end)
    end
  end
  table.sort(out, function(x, y)
    if x.start ~= y.start then return x.start < y.start end
    return x.strand > y.strand
  end)
  return out
end

--- GenBank location text for an ORF on a record of length L.
function M.location(orf, L)
  local inner
  if orf.wraps then
    inner = ("join(%d..%d,1..%d)"):format(orf.start, L, orf.stop)
  else
    inner = ("%d..%d"):format(orf.start, orf.stop)
  end
  return orf.strand == -1 and ("complement(" .. inner .. ")") or inner
end

return M
