local H = require("helpers")
local O = require("plasmid.bio.oligo")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

-- fixture written by tests/tools/oligo_oracle.py (primer3-py values + the
-- agreement/error numbers measured on exactly these cases)
local FX
local function fixture()
  if not FX then
    FX = vim.json.decode(H.read(H.root .. "/tests/oracle/expected/oligo.json"),
      { luanil = { object = true, array = true } })
  end
  return FX
end

local function rc(s)
  local C = { A = "T", C = "G", G = "C", T = "A" }
  return (s:reverse():gsub(".", C))
end

local function near(a, b, tol, msg)
  H.ok(a and b and math.abs(a - b) <= tol,
    ("%s: %s vs %s (tol %s)"):format(msg or "near", tostring(a), tostring(b), tostring(tol)))
end

-- numpy.percentile default (linear interpolation)
local function pct(xs, p)
  local t = vim.deepcopy(xs)
  table.sort(t)
  local x = (#t - 1) * p / 100
  local lo = math.floor(x)
  local f = x - lo
  return t[lo + 1] + ((t[lo + 2] or t[lo + 1]) - t[lo + 1]) * f
end

local function pairset(prs)
  local s = {}
  for _, p in ipairs(prs) do s[#s + 1] = p[1] .. ":" .. p[2] end
  table.sort(s)
  return table.concat(s, " ")
end

---------------------------------------------------------------------------
-- sequence measures on hand cases
---------------------------------------------------------------------------
case("gc_fraction", function()
  H.eq(0.5, O.gc("ACGT"))
  H.eq(1, O.gc("GGCC"))
  H.eq(0, O.gc("ATAT"))
  H.eq(0.25, O.gc("acgtNNaatt")) -- lower case accepted, N ignored: 2/8
  H.eq(0, O.gc(""))
  H.eq(0.5, O.gc("ACGU"))
end)

case("gc_clamp_last_five", function()
  H.eq(0, O.gc_clamp("GGGGGAAAAA"))
  H.eq(5, O.gc_clamp("AAAAAGCGCG"))
  H.eq(2, O.gc_clamp("ATGCGTACGTTAGCCTAGGA")) -- ...TAGGA
  H.eq(4, O.gc_clamp("ACGTACGGC")) -- ...ACGGC
  H.eq(2, O.gc_clamp("GC")) -- shorter than 5: all bases
end)

case("homopolymer_runs", function()
  H.eq(1, O.max_run("ACGT"))
  H.eq(4, O.max_run("ACGTTTTACG"))
  H.eq(6, O.max_run("AAAAAAC"))
  H.eq(5, O.max_run("ACGGGGG"))
  H.eq(0, O.max_run(""))
end)

case("dinucleotide_repeats", function()
  H.eq(4, O.max_repeat("ATATATAT"))
  H.eq(3, O.max_repeat("GATATATG")) -- ATATAT
  H.eq(5, O.max_repeat("CCACACACACAGG")) -- CACACACACA starts at an odd offset
  H.eq(1, O.max_repeat("ACGT"))
  H.eq(0, O.max_repeat("AAAAAA")) -- a homopolymer is not a dinucleotide repeat
  H.eq(2, O.max_repeat("AAGCGCTT"))
end)

---------------------------------------------------------------------------
-- structural cases (exact base pairs)
---------------------------------------------------------------------------
case("palindrome_self_dimer_pairs_end_to_end", function()
  for _, s in ipairs({ "GAATTCGAATTC", "TTTTTTGAATTCAAAAAA", "GCGCATATGCGC" }) do
    local d = O.homodimer(s)
    local want = {}
    for i = 1, #s do want[i] = { i, #s - i + 1 } end
    H.eq(want, d.pairs, s)
    H.eq(true, d.three_prime, s)
    H.ok(d.dg < -5, s .. " dg " .. d.dg)
  end
end)

case("internal_palindrome_is_not_3prime", function()
  -- GCATGCTAGC core only; the 3' TTTT cannot pair with anything opposite it
  local d = O.homodimer("ACGTAGCATGCTAGCATTTT")
  H.eq("10:15 11:14 12:13 13:12 14:11 15:10 16:9 17:8 8:17 9:16", pairset(d.pairs))
  H.eq(false, d.three_prime)
end)

case("heterodimer_with_reverse_complement", function()
  local a = "ATGCGTACGTTAGCCTAGGA"
  local d = O.heterodimer(a, rc(a))
  H.eq(#a, #d.pairs)
  for k, p in ipairs(d.pairs) do H.eq({ k, #a - k + 1 }, p) end
  H.eq(true, d.three_prime)
  -- the same 20-bp duplex scored by primer3.calc_end_stability: dG -21.10, Tm 58.24
  near(d.dg, -21.10, 0.05, "dg")
  near(d.tm, 58.24, 0.2, "tm")
end)

case("three_prime_complementarity", function()
  -- a ends in TTTGCC; b ends in GGCAAA = rc(TTTGCC)
  local a, b = "GTCAGTTCAGTTTGCC", "AAGCTTCGAGGCAAA"
  local d = O.heterodimer(a, b)
  H.eq({ { 11, 15 }, { 12, 14 }, { 13, 13 }, { 14, 12 }, { 15, 11 }, { 16, 10 } }, d.pairs)
  H.eq(true, d.three_prime)
  local e = O.end_dimer(a, b)
  H.eq(d.pairs, e.pairs)
  -- the anchored search pairs a's last base even when a better internal dimer exists
  local x = O.end_dimer("GCGCATATGCGCAAAAACC", "TTTTTTTTTTGGTTT")
  H.eq(19, x.pairs[#x.pairs][1])
end)

case("designed_hairpin_stem_and_loop", function()
  -- stem GCGCAGC, loop TTTT, rc(stem)
  local h = O.hairpin("GCGCAGCTTTTGCTGCGC")
  H.eq({ i = 1, j = 18, len = 7 }, h.stem)
  H.eq(4, h.loop)
  H.eq(1, #h.domains)
  for k = 1, 7 do H.eq({ k, 19 - k }, h.pairs[k]) end
  H.ok(h.dg < -5 and h.tm > 60, ("dg %.2f tm %.1f"):format(h.dg, h.tm))
  -- with non-pairing flanks the stem shifts by the flank length
  local f = O.hairpin("AAAAAGCGCAGCTTTTGCTGCGCAAAAA")
  H.eq({ i = 6, j = 23, len = 7 }, f.stem)
  H.eq(4, f.loop)
  -- dangling A's stabilise the stem
  H.ok(f.dg < h.dg)
end)

case("two_consecutive_hairpins", function()
  local s = "GCGCAGCTTTTGCTGCGC" .. "GGACCAGTTTTCTGGTCC"
  local h = O.hairpin(s)
  H.eq(2, #h.domains)
  H.eq({ i = 1, j = 18, len = 7 }, h.domains[1].stem)
  H.eq({ i = 19, j = 36, len = 7 }, h.domains[2].stem)
  near(h.dg, h.domains[1].dg + h.domains[2].dg, 1e-9, "energies add")
end)

case("poly_a_has_no_structure", function()
  H.eq(nil, O.hairpin(("A"):rep(20)))
  H.eq(nil, O.homodimer(("A"):rep(20)))
  H.eq(nil, O.heterodimer(("A"):rep(20), ("C"):rep(20)))
  H.eq(nil, O.hairpin("ACGT")) -- too short for a stem and a 3-nt loop
end)

case("salt_and_temperature_dependence", function()
  local s = "GAATTCGAATTC"
  local base = O.homodimer(s).dg
  H.ok(O.homodimer(s, { dv_mM = 5 }).dg < base, "more Mg2+ stabilises")
  H.ok(O.homodimer(s, { mv_mM = 200 }).dg < base, "more Na+ stabilises")
  H.ok(O.homodimer(s, { temp_c = 60 }).dg > base, "higher temperature destabilises")
  H.ok(O.homodimer(s, { dna_nM = 1000 }).tm > O.homodimer(s).tm, "Tm rises with concentration")
  -- Tm (°C) is where dG crosses zero for a unimolecular hairpin
  local h = O.hairpin("GCGCAGCTTTTGCTGCGC")
  near(O.hairpin("GCGCAGCTTTTGCTGCGC", { temp_c = h.tm }).dg, 0, 1e-6, "dG(Tm) = 0")
end)

---------------------------------------------------------------------------
-- comparison with primer3 on the fixture
---------------------------------------------------------------------------
case("end_stability_matches_primer3", function()
  local fx = fixture()
  for _, c in ipairs(fx.singles) do
    -- primer3 reports -dG of the 3' pentamer duplex, 2 decimals
    near(-O.end_stability(c.seq), c.p3_end_stability, 0.006, c.seq)
  end
  near(O.end_stability("ATGCGTACGTTAGCCTAGGA"), -2.94, 1e-9, "TAGGA")
end)

local function check_structure(name, list, key, thr, fn)
  local m = fixture().measured[name]
  local errs, agree = {}, 0
  for _, c in ipairs(list) do
    local mine = fn(c)
    local g1 = mine and mine.dg or 0
    local g2 = c[key] and c[key].dg or 0
    errs[#errs + 1] = math.abs(g1 - g2)
    if (g1 < thr) == (g2 < thr) then agree = agree + 1 end
  end
  local rate = agree / #list
  H.log(("  %s: n=%d agree %.1f%% (fixture %.1f%%), median |ddG| %.2f, p95 %.2f, max %.2f"):format(
    name, #list, 100 * rate, 100 * m.agree, pct(errs, 50), pct(errs, 95), pct(errs, 100)))
  H.ok(rate >= m.agree - 1e-9, name .. " threshold agreement dropped")
  H.ok(pct(errs, 50) <= m.median_abs_dg + 1e-6, name .. " median error grew")
  H.ok(pct(errs, 95) <= m.p95_abs_dg + 1e-6, name .. " p95 error grew")
  H.ok(pct(errs, 100) <= m.max_abs_dg + 1e-6, name .. " max error grew")
end

case("hairpin_dg_vs_primer3", function()
  local fx = fixture()
  check_structure("hairpin", fx.singles, "p3_hairpin", fx.thresholds.hairpin_dg, function(c) return O.hairpin(c.seq) end)
end)

case("homodimer_dg_vs_primer3", function()
  local fx = fixture()
  check_structure("homodimer", fx.singles, "p3_homodimer", fx.thresholds.dimer_dg, function(c) return O.homodimer(c.seq) end)
end)

case("heterodimer_dg_vs_primer3", function()
  local fx = fixture()
  check_structure("heterodimer", fx.pairs, "p3_heterodimer", fx.thresholds.dimer_dg,
    function(c) return O.heterodimer(c.a, c.b) end)
end)

case("palindromes_dg_and_tm_vs_primer3", function()
  local n = 0
  for _, c in ipairs(fixture().singles) do
    if c.kind == "palindrome" then
      local d = O.homodimer(c.seq)
      near(d.dg, c.p3_homodimer.dg, 0.30, c.seq .. " dG")
      near(d.tm, c.p3_homodimer.tm, 1.0, c.seq .. " Tm")
      n = n + 1
    end
  end
  H.ok(n >= 10, "palindromes in fixture: " .. n)
end)

---------------------------------------------------------------------------
-- assess
---------------------------------------------------------------------------
case("assess_fields_and_warnings", function()
  local a = O.assess("AGCGTCAGTCACTGGATCCA")
  H.eq(20, a.length)
  H.eq(0.55, a.gc)
  H.eq(2, a.gc_clamp)
  H.eq(2, a.max_run)
  H.eq(1, a.max_repeat)
  near(a.end_stability, -3.41, 1e-9, "end stability")
  H.ok(a.hairpin and a.homodimer, "structures present")
  H.eq({ "3' self-dimer dG -6.7 kcal/mol" }, a.warnings) -- ...GGATCC is palindromic

  H.eq({}, O.assess("CAGGTACCTTGGCGATCAGT").warnings)

  local p = O.assess(("A"):rep(20))
  H.eq({ "GC 0% outside 40-60%", "no G/C in last 5 bases", "homopolymer run of 20" }, p.warnings)

  -- thresholds are overridable
  local q = O.assess(("A"):rep(20), { thresholds = { gc_min = 0, min_gc_clamp = 0, max_run = 30 } })
  H.eq({}, q.warnings)
  H.eq(-9, O.thresholds.dimer_dg)

  local w = table.concat(O.assess("GCGCATATGCGCGCGCATATGCGC").warnings, "|")
  H.ok(w:find("self%-dimer dG") and w:find("hairpin"), w)
end)

---------------------------------------------------------------------------
-- performance
---------------------------------------------------------------------------
case("pair_of_25mers_under_10ms", function()
  local a, b = "ATGCGTACGTTAGCCTAGGATCGAT", "GCTAGCTAGGATCGATTACGGATCC"
  local function run()
    O.assess(a); O.assess(b); O.heterodimer(a, b); O.end_dimer(a, b); O.end_dimer(b, a)
  end
  for _ = 1, 20 do run() end -- JIT warm-up
  local t = {}
  for k = 1, 50 do
    local s = vim.uv.hrtime()
    run()
    t[k] = (vim.uv.hrtime() - s) / 1e6
  end
  table.sort(t)
  H.log(("  25-mer pair (assess x2, heterodimer, 3'-anchored dimers x2): median %.2f ms, max %.2f ms"):format(t[25], t[50]))
  H.ok(t[25] < 10, "median " .. t[25] .. " ms")
end)

return T
