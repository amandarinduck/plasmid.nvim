-- Nucleotide alphabet helpers. Pure Lua: no Neovim API.
local M = {}

-- IUPAC complement, case preserving. U complements to A (RNA input -> DNA output).
local pairs_ = {
  A = "T", T = "A", U = "A", G = "C", C = "G",
  R = "Y", Y = "R", K = "M", M = "K", S = "S", W = "W",
  B = "V", V = "B", D = "H", H = "D", N = "N",
}

local comp = {}
for i = 0, 255 do
  comp[i] = string.char(i)
end
for k, v in pairs(pairs_) do
  comp[k:byte()] = v
  comp[k:lower():byte()] = v:lower()
end
M.complement_map = comp

-- IUPAC expansions (uppercase).
M.iupac = {
  A = "A", C = "C", G = "G", T = "T", U = "T",
  R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC",
  B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT",
}

--- Complement without reversing.
function M.complement(s)
  return (s:gsub(".", function(c) return comp[c:byte()] end))
end

--- Reverse complement, case preserving, IUPAC aware.
function M.revcomp(s)
  return M.complement(s):reverse()
end

--- Guess alphabet from a residue string.
---@return "dna"|"rna"|"protein"
-- Nucleic acids use only IUPAC codes (ACGTU + RYKMSWBDHVN). Proteins use
-- E F I L P Q (and J O Z X *) routinely, and rarely sit mostly on A/C/G/T/U/N.
-- Degenerate DNA (NNK libraries, primers with R/Y/W...) stays nucleotide.
function M.guess(s)
  -- letters only: gap '-' and stop '*' say nothing about the alphabet
  local sample = s:sub(1, 20000):upper():gsub("[^%a]", "")
  if #sample == 0 then return "dna" end
  local core = select(2, sample:gsub("[ACGTUN]", ""))
  local non_iupac = select(2, sample:gsub("[^ACGTURYKMSWBDHVN]", ""))
  if non_iupac / #sample > 0.01 or core / #sample < 0.6 then return "protein" end
  local u = select(2, sample:gsub("U", ""))
  local t = select(2, sample:gsub("T", ""))
  if u > 0 and t == 0 then return "rna" end
  return "dna"
end

function M.is_nucleotide(alpha)
  return alpha == "dna" or alpha == "rna"
end

--- Complement for display: RNA records complement to U, not T.
function M.complement_for(s, alpha)
  local c = M.complement(s)
  if alpha == "rna" then c = c:gsub("T", "U"):gsub("t", "u") end
  return c
end

--- GC fraction over unambiguous bases.
function M.gc(s)
  local up = s:upper()
  local gc = select(2, up:gsub("[GCS]", ""))
  local at = select(2, up:gsub("[ATUW]", ""))
  if gc + at == 0 then return 0 end
  return gc / (gc + at)
end

return M
