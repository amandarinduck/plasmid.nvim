-- The pairwise alignment engine: modes, the banded path, strands, circular
-- references, differences and what they do to features. No UI here.
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local A = require("plasmid.bio.align")
local V = require("plasmid.bio.variant")
local alphabet = require("plasmid.model.alphabet")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function load(stem) return assert(pio.read(H.fixture("ncbi/" .. stem .. ".gb"))).records[1] end
local function rnd(n, seed)
  if seed then math.randomseed(seed) end
  local b, t = { "A", "C", "G", "T" }, {}
  for i = 1, n do t[i] = b[math.random(4)] end
  return table.concat(t)
end
local function mutate(s, subs, dels, ins)
  local t = {}
  for i = 1, #s do t[i] = s:sub(i, i) end
  for _ = 1, subs or 0 do
    local i = math.random(#t)
    t[i] = ({ A = "G", G = "A", C = "T", T = "C" })[t[i]:upper()] or "A"
  end
  for _ = 1, dels or 0 do table.remove(t, math.random(#t)) end
  for _ = 1, ins or 0 do table.insert(t, math.random(#t), "A") end
  return table.concat(t)
end

---------------------------------------------------------------------------
-- Modes
---------------------------------------------------------------------------

case("global_local_and_semi_modes", function()
  -- identical sequences: every column a match
  local g = A.align("ACGTACGTAC", "ACGTACGTAC", { mode = "global" })
  H.eq(20, g.score, "10 matches at +2")
  H.eq("10M", A.cigar(g))
  H.eq(1, g.identity)
  -- one substitution costs match + mismatch
  local s = A.align("ACGTACGTAC", "ACGTTCGTAC", { mode = "global" })
  H.eq(15, s.score)
  H.eq(9, s.matches)
  H.eq(1, s.mismatches)
  -- a deletion opens a gap: 11 matches, one gap of one base
  local d = A.align("ACGTACGTACGT", "ACGTCGTACGT", { mode = "global" })
  H.eq(22 - 7, d.score, "11 matches, gap open + extend")
  H.eq("4M1D7M", A.cigar(d))
  -- semi: the reference's ends are free, so the read sits inside it
  local read = A.align("TTTACGTACGTAAA", "ACGTACGT", { mode = "semi" })
  H.eq(16, read.score)
  H.eq(4, read.ref_start)
  H.eq(11, read.ref_stop)
  H.eq("8M", A.cigar(read))
  -- local: the best scoring segment of both
  local l = A.align("GGGGACGTACGTGGGG", "TTTTACGTACGTTTTT", { mode = "local" })
  H.eq(16, l.score)
  H.eq(5, l.ref_start)
  H.eq(12, l.ref_stop)
  H.eq(5, l.query_start)
  -- nothing alignable
  H.eq(nil, A.align("AAAAAAAA", "", {}))
end)

case("gap_costs_are_affine", function()
  -- one gap of three bases costs open + 3 x extend, not three gaps
  local one = A.align("ACGTACGTACGTAC", "ACGTACGTAC", { mode = "global" })
  H.eq(20 - 5 - 2 * 4, one.score, "10 matches, one gap of four bases")
  H.eq(1, one.gaps)
  H.eq(4, one.gap_bases)
  local wide = A.align("ACGTACGTACGTAC", "ACGTACGTAC", { mode = "global", gap_open = -1, gap_extend = -1 })
  H.ok(wide.score > one.score, "cheaper gaps score higher")
end)

---------------------------------------------------------------------------
-- The banded path
---------------------------------------------------------------------------

case("banded_matches_the_full_matrix", function()
  math.randomseed(7)
  local ref = load("J01749").seq:tostring():upper():sub(1, 900)
  local query = mutate(ref:sub(50, 800), 12, 2, 2)
  local full = A.align(ref, query, { mode = "semi" })
  local band = A.align(ref, query, { mode = "semi", max_cells = 1000 })
  H.ok(band.banded, "the banded path was taken")
  H.eq(full.score, band.score)
  H.eq(A.cigar(full), A.cigar(band))
  H.eq(full.ref_start, band.ref_start)
  H.eq(full.ref_stop, band.ref_stop)
  -- sequences with no k-mer in common cannot be banded at all
  local _, why = A.align(("A"):rep(3000), ("C"):rep(3000), { mode = "global", max_cells = 1000 })
  H.ok(why and why:find("%-mer"), tostring(why))
  -- two unrelated random sequences do share a few k-mers by chance; they
  -- align, but badly
  local junk = A.align(rnd(3000, 1), rnd(3000, 2), { mode = "local", max_cells = 1000 })
  H.ok(junk == nil or junk.ref_span < 300,
    ("unrelated sequences share only short segments, got %d bp"):format(junk and junk.ref_span or 0))
end)

case("long_pairs_align_in_reasonable_time", function()
  local seq = load("J01749").seq:tostring():upper()
  math.randomseed(3)
  local read = mutate(seq:sub(1200, 2200), 15, 1, 1)
  local t0 = vim.uv.hrtime()
  local aln = A.align(seq, read, { mode = "semi" })
  local ms = (vim.uv.hrtime() - t0) / 1e6
  H.ok(aln ~= nil and aln.identity > 0.95, "the read is placed")
  H.ok(math.abs(aln.ref_start - 1200) <= 5, ("placed at %d"):format(aln.ref_start))
  H.ok(ms < 500, ("1 kb read against a 4.4 kb plasmid took %.0f ms"):format(ms))
  H.log(("[perf] 1 kb read vs 4.4 kb plasmid: %.0f ms"):format(ms))
  local t1 = vim.uv.hrtime()
  local cmp = A.compare(seq, mutate(seq, 40, 0, 0), {})
  local ms2 = (vim.uv.hrtime() - t1) / 1e6
  H.ok(cmp.identity > 0.98)
  H.log(("[perf] 4.4 kb vs 4.4 kb: %.0f ms"):format(ms2))
end)

---------------------------------------------------------------------------
-- Strands and circles
---------------------------------------------------------------------------

case("the_better_strand_is_chosen", function()
  local seq = load("J01749").seq:tostring():upper()
  local rc = alphabet.revcomp(seq:sub(1000, 1600))
  local aln = A.align_strands(seq, rc, { mode = "semi" })
  H.eq(-1, aln.strand)
  H.eq(1000, aln.ref_start)
  H.eq(1600, aln.ref_stop)
  H.eq(1, aln.identity)
  local fwd = A.align_strands(seq, seq:sub(1000, 1600), { mode = "semi" })
  H.eq(1, fwd.strand)
end)

case("a_read_may_cross_the_origin", function()
  local rec = load("J01749")
  local seq = rec.seq:tostring():upper()
  local L = #seq
  local read = seq:sub(L - 400) .. seq:sub(1, 300)
  local aln = A.place(seq, read, { circular = true })
  H.eq(true, aln.wraps)
  H.eq(L - 400, aln.ref_start)
  H.eq(300, aln.ref_stop)
  H.eq(1, aln.identity)
  H.eq(701, aln.ref_span)
  -- on a linear reference the same read only partly matches
  local linear = A.place(seq, read, { circular = false })
  H.ok(linear.ref_span < 701, "a linear reference cannot wrap")
end)

case("circular_comparison_ignores_the_origin", function()
  local rec = load("J01749")
  local seq = rec.seq:tostring():upper()
  math.randomseed(5)
  -- rotate, then change three bases and delete one
  local rot = 1234
  local variant = seq:sub(rot + 1) .. seq:sub(1, rot)
  local t = {}
  for i = 1, #variant do t[i] = variant:sub(i, i) end
  for _, p in ipairs({ 500, 1500, 2500 }) do t[p] = (t[p] == "A") and "G" or "A" end
  table.remove(t, 3000)
  variant = table.concat(t)
  local aln, diffs = A.compare_differences(seq, variant, { circular = true })
  H.eq(1, aln.strand)
  H.eq(4, #diffs, "three substitutions and one deletion, and no origin artefact")
  local kinds = {}
  for _, d in ipairs(diffs) do kinds[#kinds + 1] = d.kind end
  H.eq({ "sub", "sub", "sub", "del" }, kinds)
  for i = 2, #diffs do H.ok(diffs[i].pos > diffs[i - 1].pos, "differences are in order") end
  -- the positions are the reference's own, not the rotated frame's
  H.eq(seq:sub(diffs[1].pos, diffs[1].pos), diffs[1].ref)
  -- the reverse complement of the same molecule is recognised
  local rcaln = A.compare(seq, alphabet.revcomp(variant), { circular = true })
  H.eq(-1, rcaln.strand)
  H.ok(rcaln.identity > 0.999)
end)

---------------------------------------------------------------------------
-- Differences
---------------------------------------------------------------------------

case("differences_are_listed_in_reference_coordinates", function()
  local ref = "ACGTACGTACGTACGTACGT"
  -- substitution at 5, deletion of 2 at 10, insertion of 3 after 15
  local query = "ACGT" .. "T" .. "CGTA" .. "GTACG" .. "AAA" .. "TACGT"
  local aln = A.align(ref, query, { mode = "global" })
  local diffs = A.differences(ref, query, aln)
  local kinds, pos = {}, {}
  for _, d in ipairs(diffs) do kinds[#kinds + 1] = d.kind; pos[#pos + 1] = d.pos end
  H.ok(vim.tbl_contains(kinds, "sub"), vim.inspect(diffs))
  H.eq(5, pos[1])
  H.eq("A", diffs[1].ref)
  H.eq("T", diffs[1].query)
  -- neighbouring substitutions are grouped into one run
  local run = A.differences("ACGTACGT", "ACTTACGT", A.align("ACGTACGT", "ACTTACGT", { mode = "global" }))
  H.eq(1, #run)
  H.eq(1, run[1].len)
  local two = A.align("ACGTACGT", "ACTAACGT", { mode = "global" })
  local tworun = A.differences("ACGTACGT", "ACTAACGT", two)
  H.eq(1, #tworun, "two neighbouring substitutions are one run")
  H.eq(2, tworun[1].len)
  H.eq("GT", tworun[1].ref)
  H.eq("TA", tworun[1].query)
end)

---------------------------------------------------------------------------
-- What the differences hit
---------------------------------------------------------------------------

case("differences_are_placed_in_features_and_codons", function()
  local rec = load("J01749")            -- tet 86..1276 (+), bla complement(3293..4153)
  local seq = rec.seq:tostring():upper()
  local t = {}
  for i = 1, #seq do t[i] = seq:sub(i, i) end
  t[150] = (t[150] == "A") and "G" or "A"       -- inside tet
  t[3400] = (t[3400] == "C") and "T" or "C"     -- inside bla, reverse strand
  local query = table.concat(t)
  local aln, diffs = A.compare_differences(seq, query, { circular = true })
  V.annotate(rec, diffs)
  H.eq(2, #diffs)
  local by_pos = {}
  for _, d in ipairs(diffs) do by_pos[d.pos] = d end
  local tet, bla = by_pos[150], by_pos[3400]
  H.eq("tet", tet.feature)
  H.eq(22, tet.codon)
  H.eq("missense", tet.effect)
  H.eq("G", tet.aa_from)
  H.eq("D", tet.aa_to)
  H.eq("bla", bla.feature)
  H.eq(252, bla.codon, "a reverse-strand CDS is numbered from its own start")
  H.eq("missense", bla.effect)
  H.ok(V.describe(tet, "bp"):find("codon 22"), V.describe(tet, "bp"))
end)

case("frameshifts_and_silent_changes_are_classified", function()
  local rec = load("J01749")
  local seq = rec.seq:tostring():upper()
  -- a deletion inside tet
  local t = {}
  for i = 1, #seq do t[i] = seq:sub(i, i) end
  table.remove(t, 600)
  local _, diffs = A.compare_differences(seq, table.concat(t), { circular = true })
  local by_feature = V.by_feature(rec, diffs)
  local cds
  for _, e in ipairs(by_feature) do if e.type == "CDS" and e.label == "tet" then cds = e end end
  H.ok(cds ~= nil, "the tet CDS is reported")
  H.eq("frameshift", cds.kind)
  H.eq(true, cds.frameshift)
  H.ok(cds.premature_stop ~= nil, "a frameshift brings a stop forward")
  -- a third-position change that keeps the amino acid is silent
  local codon_start = 86
  local silent
  for pos = codon_start + 2, codon_start + 200, 3 do
    local t2 = {}
    for i = 1, #seq do t2[i] = seq:sub(i, i) end
    for _, b in ipairs({ "A", "C", "G", "T" }) do
      if b ~= t2[pos] then
        t2[pos] = b
        local _, ds = A.compare_differences(seq, table.concat(t2), { circular = true })
        V.annotate(rec, ds)
        if #ds == 1 and ds[1].effect == "silent" then silent = ds[1] break end
      end
    end
    if silent then break end
  end
  H.ok(silent ~= nil, "a silent change exists in tet's wobble positions")
  H.eq(silent.aa_from, silent.aa_to)
end)

---------------------------------------------------------------------------
-- Proteins
---------------------------------------------------------------------------

case("protein_alignment_uses_blosum62", function()
  local recs = assert(pio.read(H.fixture("protein/proteins.faa"))).records
  local rec = recs[1]
  for _, r in ipairs(recs) do if (r.id or ""):find("HBA_HUMAN") then rec = r end end
  local p = rec.seq:tostring():upper()
  local self_aln = A.align(p, p, { mode = "global", protein = true })
  -- the score of a sequence against itself is the sum of its diagonal entries
  local data = require("plasmid.data.blosum62")
  local index = {}
  for i = 1, #data.alphabet do index[data.alphabet:sub(i, i)] = i end
  local want = 0
  for i = 1, #p do
    local c = p:sub(i, i)
    want = want + (data.scores[c] and data.scores[c][index[c]] or 0)
  end
  H.eq(want, self_aln.score)
  H.eq(1, self_aln.identity)
  -- a substitution to a similar residue costs less than to a dissimilar one
  local like = p:sub(1, 20):gsub("^(.)(.)", "%1K")
  local unlike = p:sub(1, 20):gsub("^(.)(.)", "%1W")
  local a1 = A.align(p:sub(1, 20), like, { mode = "global", protein = true })
  local a2 = A.align(p:sub(1, 20), unlike, { mode = "global", protein = true })
  H.ok(a1.score >= a2.score, ("K %d vs W %d"):format(a1.score, a2.score))
end)

return T
