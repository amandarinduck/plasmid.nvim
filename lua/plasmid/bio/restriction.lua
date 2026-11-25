-- Restriction enzymes: cut sites on both strands, overhangs, digests.
--
-- An enzyme is { name, site, cut5, cut3, dam?, dcm?, cpg?, suppliers? } where
-- `site` is the recognition sequence (IUPAC) read 5'->3' on the strand the
-- enzyme binds, and cut5 / cut3 locate the cuts in that frame: the top-strand
-- cut falls just before site offset cut5 (0 = before the first base of the
-- site) and the bottom-strand cut just before offset cut3, both measured along
-- the top strand. EcoRI G^AATTC: cut5 = 1, cut3 = 5. BsaI GGTCTC(1/5): 7 / 11.
-- Values may be negative or exceed #site (Type IIS enzymes cut outside).
--
-- Cut positions are reported as the first base after the cut on the top
-- strand (1-based), as Biopython's Restriction.search does.
local search = require("plasmid.bio.search")

local M = {}

-- Enzymes that cut only methylated DNA: listed, but not counted as cutting
-- unless named explicitly.
M.NEEDS_METHYLATION = { DpnI = "dam" }

local cached

--- The enzyme list: the bundled common set, replaced by an imported REBASE
--- set when one has been imported (see :DnaEnzymeImport).
function M.enzymes()
  if cached then return cached end
  local list
  local ok, rebase = pcall(function() return require("plasmid.rebase").load() end)
  if ok and rebase and #rebase > 0 then list = rebase else list = require("plasmid.data.enzymes") end
  cached = M.index(list)
  return cached
end

function M.reset() cached = nil end

--- Build lookups (by_name, case-insensitive) for a list of enzymes.
function M.index(list)
  local by = {}
  for _, e in ipairs(list) do
    e.motif = e.motif or search.parse_motif(e.site)
    by[e.name:lower()] = e
  end
  list.by_name = by
  return list
end

function M.get(name)
  return M.enzymes().by_name[(name or ""):lower()]
end

--- 5' overhang > 0, 3' overhang < 0, blunt = 0 (length in bases).
function M.overhang(e) return e.cut3 - e.cut5 end

function M.end_kind(e)
  local o = M.overhang(e)
  if o == 0 then return "blunt" end
  return o > 0 and ("5'+%d"):format(o) or ("3'+%d"):format(-o)
end

local function wrap(p, L)
  return (p - 1) % L + 1
end

--- Prepare a sequence for scanning many enzymes. When the sequence is plain
--- ACGT a 4-mer position index is built so that each enzyme only verifies
--- candidate positions instead of scanning the whole sequence.
function M.prepare(rec)
  local seq = type(rec.seq) == "string" and rec.seq or rec.seq:tostring()
  local L = #seq
  local hay = seq:upper():gsub("U", "T")
  local circular = rec.topology == "circular"
  if circular then hay = hay .. hay:sub(1, math.min(L, 64)) end
  local P = { seq = seq, hay = hay, L = L, circular = circular }
  if L >= 2000 and not hay:find("[^ACGT]") then
    local code = { [65] = 0, [67] = 1, [71] = 2, [84] = 3 }
    local idx = {}
    local byte = string.byte
    local c = 0
    for i = 1, #hay do
      c = (c * 4 + code[byte(hay, i)]) % 256
      if i >= 4 then
        local list = idx[c]
        if not list then list = {}; idx[c] = list end
        list[#list + 1] = i - 3
      end
    end
    P.kmers = idx
  end
  return P
end

local BASES = { A = "A", C = "C", G = "G", T = "T", U = "T", R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT",
  M = "AC", B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT" }
local CODE = { A = 0, C = 1, G = 2, T = 3 }

-- Most selective 4-letter window of a motif: offset (0-based) and the 4-mer codes.
local function best_window(motif)
  local best, best_off, best_n
  for off = 0, #motif - 4 do
    local codes = { 0 }
    for j = 1, 4 do
      local set = BASES[motif:sub(off + j, off + j)]
      local nxt = {}
      for _, c in ipairs(codes) do
        for k = 1, #set do nxt[#nxt + 1] = c * 4 + CODE[set:sub(k, k)] end
      end
      codes = nxt
    end
    if not best_n or #codes < best_n then best, best_off, best_n = codes, off, #codes end
  end
  return best_off, best
end

local function indexed_scan(P, motif, strand, out)
  local m = #motif
  local last = P.circular and P.L or (P.L - m + 1)
  local off, codes = best_window(motif)
  local pat = "^" .. search.lua_pattern(motif)
  local hay = P.hay
  local starts = {}
  for _, c in ipairs(codes) do
    for _, pos in ipairs(P.kmers[c] or {}) do
      local s = pos - off
      if s >= 1 and s <= last and hay:find(pat, s) then starts[#starts + 1] = s end
    end
  end
  table.sort(starts)
  for _, s in ipairs(starts) do
    local e = s + m - 1
    if e > P.L then e = e - P.L end
    out[#out + 1] = { start = s, stop = e, strand = strand }
  end
end

-- All matches of an enzyme site: { start, stop, strand } (palindromes once).
local function find_site(P, e)
  local motif = e.motif
  if not P.kmers or #motif < 4 or #motif > 64 then
    return search.find_nuc(P.seq, motif, { circular = P.circular, hay = P.hay })
  end
  local out = {}
  indexed_scan(P, motif, 1, out)
  local rc = search.revcomp_motif(motif)
  if rc ~= motif then indexed_scan(P, rc, -1, out) end
  return out
end

--- Sites of one enzyme: list of
--- { enzyme, start, stop, strand, cut (top strand), cut_bottom, overhang },
--- sorted by top-strand cut. Linear molecules keep only cuts whose top-strand
--- cut falls inside the sequence (2..L).
function M.sites(P, e)
  local m = #e.site
  local hits = find_site(P, e)
  local out = {}
  for _, h in ipairs(hits) do
    local top, bottom
    if h.strand == 1 then
      top, bottom = h.start + e.cut5, h.start + e.cut3
    else
      top, bottom = h.start + m - e.cut3, h.start + m - e.cut5
    end
    local keep = true
    if P.circular then
      top, bottom = wrap(top, P.L), wrap(bottom, P.L)
    else
      -- as Biopython: a linear cut counts when the top-strand cut is inside
      keep = top >= 2 and top <= P.L
    end
    if keep then
      out[#out + 1] = { enzyme = e, start = h.start, stop = h.stop, strand = h.strand,
        cut = top, cut_bottom = bottom, overhang = M.overhang(e) }
    end
  end
  table.sort(out, function(a, b)
    if a.cut ~= b.cut then return a.cut < b.cut end
    if a.cut_bottom ~= b.cut_bottom then return a.cut_bottom < b.cut_bottom end
    return a.strand > b.strand
  end)
  -- the same cut reached from two overlapping sites is one cut
  local dedup = {}
  for _, s in ipairs(out) do
    local last = dedup[#dedup]
    if not (last and last.cut == s.cut and last.cut_bottom == s.cut_bottom) then dedup[#dedup + 1] = s end
  end
  return dedup
end

--- Scan many enzymes: returns { [name] = sites } and the per-name counts.
function M.scan(rec, enzymes)
  local P = M.prepare(rec)
  local res, counts = {}, {}
  for _, e in ipairs(enzymes or M.enzymes()) do
    local s = M.sites(P, e)
    res[e.name] = s
    counts[e.name] = #s
  end
  return res, counts
end

-- Dam (GATC) and Dcm (CCWGG) sites overlapping a recognition site.
local function overlaps(P, s, motif)
  local m = #motif
  local L = P.L
  local from = s.start - (m - 1)
  local span = (s.stop >= s.start and s.stop or s.stop + L) - s.start + 1
  local window_len = span + 2 * (m - 1)
  local text = {}
  for i = 0, window_len - 1 do
    local p = from + i
    if P.circular then p = wrap(p, L) end
    text[#text + 1] = (p >= 1 and p <= L) and P.hay:sub(p, p) or "-"
  end
  local w = table.concat(text)
  return w:find(motif == "GATC" and "GATC" or "CC[AT]GG") ~= nil
end

--- Methylation warnings for a site: list of "dam"/"dcm" whose methylated
--- sequence overlaps a site of an enzyme sensitive to it.
function M.methylation(P, s)
  local e = s.enzyme
  local out = {}
  if e.dam and overlaps(P, s, "GATC") then out[#out + 1] = "dam" end
  if e.dcm and overlaps(P, s, "CCWGG") then out[#out + 1] = "dcm" end
  return out
end

-- Cuts grouped by top-strand position: two enzymes cutting the same phosphate
-- bond are one cut, as in M.digest.
local function cut_positions(cuts)
  local by_pos, pos = {}, {}
  for _, c in ipairs(cuts) do
    if not by_pos[c.cut] then by_pos[c.cut] = {}; pos[#pos + 1] = c.cut end
    table.insert(by_pos[c.cut], c)
  end
  table.sort(pos)
  return by_pos, pos
end

--- Every fragment an incomplete digestion can leave: any span bounded by two
--- cuts, or (on a linear molecule) by a cut and an end. Same shape as
--- M.digest's fragments, plus `partial = true` and `skipped` — the number of
--- the chosen set's sites that lie inside the fragment and were not cut, so
--- the complete-digest fragments are exactly those with `skipped == 0`.
---
--- The answer grows as the square of the site count: n cuts on a circle give
--- n² spans (n of them the whole molecule, linearised at one site) plus the
--- intact circle, and n cuts on a linear molecule give C(n+2, 2) spans over
--- its n+2 boundaries. A set with many sites therefore has more products than
--- anyone can read, so `opts.max_sites` (default 12, which is 145 fragments)
--- refuses to enumerate rather than returning a list nobody asked for, and
--- `opts.max_fragments` caps what is returned. Both are reported in `info`.
---
--- Returns the fragments, sorted by `skipped` and then position — the
--- complete digest first, then the one-site read-throughs, and so on — and
--- `info` = { sites, total, shown, capped, too_many, max_sites, max_fragments }.
function M.partial(rec, cuts, opts)
  opts = opts or {}
  local L = rec.seq:len()
  local circular = rec.topology == "circular"
  local by_pos, pos = cut_positions(cuts)
  local n = #pos
  local max_sites = opts.max_sites or 12
  local max_fragments = opts.max_fragments or 500
  local info = { sites = n, max_sites = max_sites, max_fragments = max_fragments }
  if n > max_sites then
    info.too_many = true
    info.total = circular and (n * n + 1) or ((n + 2) * (n + 1) / 2)
    info.shown = 0
    return {}, info
  end
  local frags = {}
  if circular then
    -- no site cut at all: the plasmid comes back intact
    frags[1] = { start = 1, stop = L, length = L, uncut = true, partial = true, skipped = n }
    for i = 1, n do
      for j = 1, n do
        local a, b = pos[i], pos[j]
        local len = (b - a) % L
        if len == 0 then len = L end   -- cut at one site only: the whole circle, linearised
        local stop = b - 1
        if stop < 1 then stop = L end
        frags[#frags + 1] = { start = a, stop = stop, length = len,
          left = by_pos[a], right = by_pos[b], partial = true, skipped = (j - i - 1) % n }
      end
    end
  else
    -- boundaries of a linear molecule: its two ends and every cut between them
    local b = { 1 }
    for _, p in ipairs(pos) do b[#b + 1] = p end
    b[#b + 1] = L + 1
    for i = 1, #b - 1 do
      for j = i + 1, #b do
        frags[#frags + 1] = { start = b[i], stop = b[j] - 1, length = b[j] - b[i],
          left = (i > 1) and by_pos[b[i]] or nil,
          right = (j < #b) and by_pos[b[j]] or nil,
          -- the span between the molecule's own two ends is the molecule: no
          -- site was cut, and it is reported the way a digest with no cut is
          uncut = (i == 1 and j == #b) or nil,
          partial = true, skipped = j - i - 1 }
      end
    end
  end
  table.sort(frags, function(x, y)
    if x.skipped ~= y.skipped then return x.skipped < y.skipped end
    if x.start ~= y.start then return x.start < y.start end
    return x.length < y.length
  end)
  info.total = #frags
  if #frags > max_fragments then
    for k = #frags, max_fragments + 1, -1 do frags[k] = nil end
    info.capped = true
  end
  info.shown = #frags
  return frags, info
end

--- Digest with a set of cuts (list of sites from M.sites, any enzymes).
--- Returns fragments sorted by position: { start, stop (stop < start wraps),
--- length, left = site|nil, right = site|nil }. Fragment boundaries follow the
--- top-strand cut; the overhangs are described by the sites at each end.
function M.digest(rec, cuts)
  local L = rec.seq:len()
  local circular = rec.topology == "circular"
  local by_pos, pos = cut_positions(cuts)
  local frags = {}
  local function ends(p) return by_pos[p] end
  if #pos == 0 then
    return { { start = 1, stop = L, length = L, uncut = true } }
  end
  if circular then
    for i = 1, #pos do
      local a, b = pos[i], pos[i % #pos + 1]
      local stop = b - 1
      if stop < 1 then stop = L end
      local len = (b - a) % L
      if len == 0 then len = L end
      frags[#frags + 1] = { start = a, stop = stop, length = len, left = ends(a), right = ends(b) }
    end
  else
    local prev = 1
    for _, p in ipairs(pos) do
      frags[#frags + 1] = { start = prev, stop = p - 1, length = p - prev, left = prev > 1 and ends(prev) or nil, right = ends(p) }
      prev = p
    end
    frags[#frags + 1] = { start = prev, stop = L, length = L - prev + 1, left = ends(prev) }
  end
  return frags
end

return M
