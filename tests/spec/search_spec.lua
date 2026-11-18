local H = require("helpers")
local pio = require("plasmid.io")
local search = require("plasmid.bio.search")
local translate = require("plasmid.bio.translate")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

-- Independent oracle: IUPAC codes as explicit letter sets, own complement table.
local EXPAND = {
  A = "A", C = "C", G = "G", T = "TU", U = "TU",
  R = "AG", Y = "CTU", S = "CG", W = "ATU", K = "GTU", M = "AC",
  B = "CGTU", D = "AGTU", H = "ACTU", V = "ACG", N = "ACGTU",
}
local MEANS = { A = "A", C = "C", G = "G", T = "T", U = "T", R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC",
  B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT" }
local COMP = { A = "T", C = "G", G = "C", T = "A", U = "A", R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
  B = "V", V = "B", D = "H", H = "D", N = "N" }
local function rc(s)
  local t = {}
  for i = #s, 1, -1 do t[#t + 1] = COMP[s:sub(i, i)] end
  return table.concat(t)
end
-- sequence letter x matches motif letter p when every base x may be is allowed by p
local function ok(x, p)
  local mx, mp = MEANS[x], MEANS[p]
  if not mx then return false end
  for i = 1, #mx do if not mp:find(mx:sub(i, i), 1, true) then return false end end
  return true
end
local function brute(seq, motif, circular, k)
  seq = seq:upper()
  local L, m = #seq, #motif
  local out = {}
  local function at(i) if i > L then i = i - L end return seq:sub(i, i) end
  local function try(pat, strand)
    for s = 1, circular and L or (L - m + 1) do
      local mm = 0
      for j = 1, m do if not ok(at(s + j - 1), pat:sub(j, j)) then mm = mm + 1 end end
      if mm <= k then out[#out + 1] = ("%d:%d"):format(s, strand) end
    end
  end
  try(motif, 1)
  if rc(motif) ~= motif then try(rc(motif), -1) end
  table.sort(out)
  return out
end
local function keys(hits)
  local out = {}
  for _, h in ipairs(hits) do out[#out + 1] = ("%d:%d"):format(h.start, h.strand) end
  table.sort(out)
  return out
end

case("nucleotide_search_matches_brute_force", function()
  math.randomseed(5)
  local codes = "RYSWKMBDHVN"
  local n = 0
  for _, rel in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001422.gb", "ncbi/L09137.gb" }) do
    local rec = assert(pio.read(H.fixture(rel))).records[1]
    local seq = rec.seq:tostring()
    for trial = 1, 40 do
      -- motif taken from the sequence (so it occurs), with some positions made degenerate
      local m = math.random(4, 9)
      local s = math.random(1, #seq - m)
      local t = {}
      for i = 1, m do
        local c = seq:sub(s + i - 1, s + i - 1):upper()
        if math.random() < 0.25 then local r = math.random(#codes); c = codes:sub(r, r) end
        t[i] = c
      end
      local motif = table.concat(t)
      local k = trial % 4 == 0 and 1 or 0
      for _, circ in ipairs({ false, true }) do
        local got = keys(search.find_nuc(seq, motif, { circular = circ, mismatches = k }))
        local want = brute(seq, motif, circ, k)
        H.eq(table.concat(want, " "), table.concat(got, " "), ("%s %s k=%d circular=%s"):format(rel, motif, k, tostring(circ)))
        n = n + #want
      end
    end
  end
  H.log(("    [search] 120 motifs x 2 topologies (240 searches): %d hits identical to brute force"):format(n))
end)

case("specific_cases", function()
  -- EcoRI is palindromic: one hit per site, strand 1; pUC19 has one site at 396
  local rec = assert(pio.read(H.fixture("ncbi/L09137.gb"))).records[1]
  local hits = search.find(rec, "gaattc")
  H.eq(1, #hits)
  H.eq(396, hits[1].start)
  H.eq(401, hits[1].stop)
  H.eq(1, hits[1].strand)
  -- across the origin of a circular record
  local seq = rec.seq:tostring()
  local L = #seq
  local motif = (seq:sub(L - 2) .. seq:sub(1, 4)):upper()
  local wrap
  for _, h in ipairs(search.find_nuc(seq, motif, { circular = true })) do if h.start == L - 2 then wrap = h end end
  H.ok(wrap, "hit across the origin")
  H.eq(4, wrap.stop)
  H.eq(0, #vim.tbl_filter(function(h) return h.start == L - 2 end, search.find_nuc(seq, motif, { circular = false })))
  -- N in the sequence only matches N
  H.eq(0, #search.find_nuc("AANAA", "AAAAA"))
  H.eq(1, #search.find_nuc("AANAA", "AANAA"))
  H.eq(1, #search.find_nuc("AARAA", "AARAA"))
  H.eq(1, #search.find_nuc("AARAA", "AANAA"), "R in the sequence fits N in the motif")
  -- errors
  H.eq(nil, (search.parse_motif("GAXTC")))
  H.eq("aa", select(2, search.parse_motif("aa:MKV")))
end)

case("translated_search_maps_back_to_codons", function()
  local rec = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  local seq = rec.seq:tostring()
  -- residues 10..17 of the tet protein (forward) and of bla (reverse strand)
  local tet, bla
  for _, f in ipairs(rec.features) do
    if f.type == "CDS" and f:get("gene") == "tet" then tet = f end
    if f.type == "CDS" and f:get("gene") == "bla" then bla = f end
  end
  for _, f in ipairs({ tet, bla }) do
    local prot = f:get("translation")
    local motif = prot:sub(10, 17)
    local hits = search.find(rec, "aa:" .. motif)
    local found
    for _, h in ipairs(hits) do
      local nt = h.strand == 1 and seq:sub(h.start, h.stop) or require("plasmid.model.alphabet").revcomp(seq:sub(h.start, h.stop))
      H.eq(motif, translate.translate(nt, { table = 11 }), "hit translates to the motif")
      local p = require("plasmid.model.location").positions(f.parts)
      if (h.strand == 1 and h.start == p[28]) or (h.strand == -1 and h.stop == p[28]) then found = true end
    end
    H.ok(found, f:get("gene") .. ": hit at codon 10 of the CDS")
  end
  -- protein records search directly; X is a wildcard
  local prec = require("plasmid.model.record").Record.new({ seq = "MKVLAAGIVGMKVL", alphabet = "protein" })
  local hits = search.find(prec, "MKV")
  H.eq(2, #hits)
  H.eq(2, #search.find(prec, "MXV"))
  H.eq(11, hits[2].start)
end)

case("performance", function()
  local rec = assert(pio.read(H.fixture("ncbi/NC_001133.gb"))).records[1]
  local t0 = vim.uv.hrtime()
  local hits = search.find(rec, "GATC")
  local t1 = vim.uv.hrtime()
  local hits2 = search.find(rec, "GGTCTC", 1)
  local t2 = vim.uv.hrtime()
  H.ok(#hits > 100 and #hits2 > 0)
  H.log(("    [perf] 230 kb: GATC exact %d hits in %.0f ms; GGTCTC with 1 mismatch %d hits in %.0f ms"):format(
    #hits, (t1 - t0) / 1e6, #hits2, (t2 - t1) / 1e6))
  H.ok((t2 - t1) / 1e6 < 2000)
end)

return T
