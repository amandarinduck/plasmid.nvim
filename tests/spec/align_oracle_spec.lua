-- The alignment engine against Biopython's PairwiseAligner.
--
-- tests/oracle/dump_alignments.py runs ~1200 pairs through Biopython 1.88 --
-- three modes plus two one-sided variants, eight scoring schemes, nucleotide
-- and protein -- and commits the optimal score and one optimal alignment to
-- tests/oracle/expected/alignments.json. Nothing here touches Python or the
-- network: the JSON is replayed, and for every case the engine's score must
-- equal Biopython's exactly. Where the two spell the alignment differently,
-- the engine's own path is rescored under the same scheme and has to come
-- out at the optimum: a co-optimal path is allowed, a worse one is not.
--
-- The second half cuts simulated reads out of the test plasmids, with known
-- mutations, and checks that A.place puts them back where they came from and
-- A.differences recovers what was done to them.
--
-- docs/alignment_report.md writes up the corpus, the agreement, the four
-- engine bugs this found, and the timings.
local H = require("helpers")

require("plasmid").setup({})
local A = require("plasmid.bio.align")
local O = dofile(H.root .. "/tests/tools/align_oracle.lua")
local R = dofile(H.root .. "/tests/tools/align_reads.lua")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local doc, replay
local function corpus()
  if not doc then doc = O.load(H.root) end
  return doc
end
--- The whole replay, once: every later case reads these counts.
local function results()
  if not replay then
    local counts, bad, by_family = O.run(corpus())
    replay = { counts = counts, bad = bad, by_family = by_family }
  end
  return replay
end

---------------------------------------------------------------------------
-- The corpus itself
---------------------------------------------------------------------------

case("the_committed_corpus_covers_every_mode_and_scheme", function()
  local d = corpus()
  H.ok(#d.cases > 1000, ("%d cases"):format(#d.cases))
  H.eq("1.88", d.biopython)
  local modes, schemes, fams = {}, {}, {}
  for _, c in ipairs(d.cases) do
    modes[c.mode] = (modes[c.mode] or 0) + 1
    schemes[c.scheme] = (schemes[c.scheme] or 0) + 1
    fams[c.family] = (fams[c.family] or 0) + 1
  end
  for _, m in ipairs({ "global", "local", "semi", "semi_q", "overlap" }) do
    H.ok((modes[m] or 0) >= 100, ("mode %s has %d cases"):format(m, modes[m] or 0))
  end
  for _, s in ipairs({ "blastn", "unit", "longgap", "linear", "steep",
                       "blastp", "cheap", "near_linear" }) do
    H.ok((schemes[s] or 0) > 0, "scheme " .. s)
  end
  for _, f in ipairs({ "random_unrelated", "random_related", "fixture_mutated", "revcomp",
                       "hangoff_left", "hangoff_right", "hangoff_both", "short_query",
                       "protein_pair", "protein_mutated", "protein_fragment",
                       "disjoint", "empty_query" }) do
    H.ok((fams[f] or 0) > 0, "family " .. f)
  end
  local nprot = 0
  for _, c in ipairs(d.cases) do if c.protein then nprot = nprot + 1 end end
  H.ok(nprot >= 200, ("%d protein cases"):format(nprot))
  H.log(("[oracle] %d cases, %d protein, Biopython %s"):format(#d.cases, nprot, d.biopython))
end)

---------------------------------------------------------------------------
-- Agreement
---------------------------------------------------------------------------

case("every_score_equals_biopythons", function()
  local r = results()
  local lines = {}
  for i = 1, math.min(#r.bad, 10) do
    local b = r.bad[i]
    lines[#lines + 1] = ("%s %s %s/%s/%s: %s"):format(b.case.id, b.result.status,
      b.case.family, b.case.mode, b.case.scheme, tostring(b.result.detail))
  end
  H.eq(0, #r.bad, table.concat(lines, "\n"))
  H.eq(0, r.counts.score or 0)
  H.eq(0, r.counts.suboptimal or 0)
  H.eq(0, r.counts.nil_engine or 0)
  H.eq(0, r.counts.nil_oracle or 0)
  local total = (r.counts.exact or 0) + (r.counts.cooptimal or 0) + (r.counts.nil_both or 0)
  H.eq(#corpus().cases, total)
  H.log(("[oracle] %d exact, %d co-optimal, %d both decline, 0 disagreements")
    :format(r.counts.exact or 0, r.counts.cooptimal or 0, r.counts.nil_both or 0))
end)

case("a_unique_optimum_is_reproduced_exactly", function()
  -- Where Biopython reports that only one optimal alignment exists, a
  -- co-optimal answer is not an excuse: the engine must produce that one.
  local d = corpus()
  local n, miss = 0, {}
  for _, c in ipairs(d.cases) do
    if c.n_optimal == 1 then
      n = n + 1
      local res = O.check(d, c)
      if res.status ~= "exact" then
        miss[#miss + 1] = ("%s %s/%s: %s"):format(c.id, c.family, c.mode, tostring(res.detail))
      end
    end
  end
  H.ok(n >= 300, ("%d cases with a unique optimum"):format(n))
  H.eq(0, #miss, table.concat(miss, "\n"))
  H.log(("[oracle] %d cases with a unique optimal alignment, all reproduced"):format(n))
end)

case("each_mode_and_scheme_agrees_on_its_own", function()
  -- the same replay, sliced, so a mode-specific regression names itself
  local d = corpus()
  for _, key in ipairs({ "mode", "scheme" }) do
    local seen = {}
    for _, c in ipairs(d.cases) do seen[c[key]] = true end
    for value in pairs(seen) do
      local counts = O.run(d, function(c) return c[key] == value end)
      H.eq(0, counts.score or 0, key .. " " .. value)
      H.eq(0, counts.suboptimal or 0, key .. " " .. value)
      H.eq(0, counts.nil_engine or 0, key .. " " .. value)
    end
  end
end)

case("proteins_score_through_blosum62", function()
  local d = corpus()
  local counts = O.run(d, function(c) return c.protein end)
  H.eq(0, counts.score or 0)
  H.eq(0, counts.suboptimal or 0)
  -- the matrix itself, against the values Biopython reports
  local o = { protein = true }
  local score_of = A._scorer(o)
  H.eq(4, (score_of(("A"):byte(), ("A"):byte())))
  H.eq(11, (score_of(("W"):byte(), ("W"):byte())))
  H.eq(-3, (score_of(("W"):byte(), ("A"):byte())))
  H.eq(0, (score_of(("X"):byte(), ("A"):byte())))
  H.eq(-1, (score_of(("X"):byte(), ("X"):byte())))
  H.log(("[oracle] %d protein cases: %d exact, %d co-optimal")
    :format((counts.exact or 0) + (counts.cooptimal or 0), counts.exact or 0,
      counts.cooptimal or 0))
end)

case("the_better_strand_is_the_one_reported", function()
  -- the revcomp family carries Biopython's score for the query and for its
  -- reverse complement; align_strands must return the larger, with the sign
  local d = corpus()
  local n, bad = 0, {}
  for _, c in ipairs(d.cases) do
    if c.rc_score then
      n = n + 1
      local aln = A.align_strands(c.ref, c.query, O.opts(d, c))
      local want = math.max(c.score, c.rc_score)
      local want_strand = c.rc_score > c.score and -1 or 1
      if not aln or aln.score ~= want or aln.strand ~= want_strand then
        bad[#bad + 1] = ("%s: got %s/%s, want %d/%d"):format(c.id,
          aln and aln.score or "nil", aln and aln.strand or "-", want, want_strand)
      end
    end
  end
  H.ok(n >= 40, ("%d strand cases"):format(n))
  H.eq(0, #bad, table.concat(bad, "\n"))
end)

case("sequences_with_nothing_in_common_are_declined_or_scored_as_biopython_does", function()
  -- Local mode has no positive-scoring segment and returns nil; Biopython
  -- reports no alignments. Global and semi still have an optimum -- for two
  -- unrelated sequences it is the one that gaps past them both.
  local d = corpus()
  local counts = O.run(d, function(c) return c.family == "disjoint" or c.family == "empty_query" end)
  H.eq(0, counts.score or 0)
  H.eq(0, counts.nil_engine or 0)
  H.eq(0, counts.nil_oracle or 0)
  H.ok((counts.nil_both or 0) > 0, "some cases are declined by both")
  -- an empty query is refused outright rather than scored
  local aln, why = A.align("ACGTACGTAC", "", {})
  H.eq(nil, aln)
  H.ok(why:find("empty"), tostring(why))
end)

---------------------------------------------------------------------------
-- The banded path
---------------------------------------------------------------------------

case("the_band_finds_the_full_matrix_answer_on_reads_and_plasmids", function()
  -- Pairs of the size that really does take the band: reads of 300-800 bp
  -- against each test plasmid at four error rates, and each plasmid against
  -- a mutated copy of itself.
  local b = R.band_check(H.root, 8, 4242, 800)
  H.ok(b.n >= 25, ("%d pairs"):format(b.n))
  H.eq(0, b.declined, table.concat(b.cases, "\n"))
  H.eq(0, b.worse, table.concat(b.cases, "\n"))
  H.log(("[oracle] band = full matrix on %d read/plasmid pairs"):format(b.n))
end)

---------------------------------------------------------------------------
-- Simulated reads
---------------------------------------------------------------------------

local reads
local function read_run()
  if not reads then reads = R.run(H.root, 24, 90210) end
  return reads
end

case("simulated_reads_are_placed_on_the_right_strand_and_position", function()
  local r = read_run()
  local t = r.total
  H.ok(t.n >= 250, ("%d reads"):format(t.n))
  H.eq(t.n, t.placed, table.concat(r.failures, "\n"))
  H.eq(t.n, t.strand_ok, "every read's strand")
  H.eq(t.n, t.pos_ok, table.concat(r.failures, "\n"))
  H.eq(t.n, t.wrap_ok, table.concat(r.failures, "\n"))
  H.ok(t.wrapping >= 100, ("%d reads cross the origin"):format(t.wrapping))
  H.eq(t.wrapping, t.wrap_detected)
  H.log(("[oracle] %d simulated reads (%d across the origin): all placed within 2 bp, "
    .. "right strand, right wrap flag"):format(t.n, t.wrapping))
end)

case("the_differences_recovered_are_the_ones_introduced", function()
  local r = read_run()
  -- per error rate, pooled over the three plasmids
  local pool = {}
  for _, row in ipairs(r.per_rate) do
    local p = pool[row.rate] or { n = 0, exact = 0, tp = 0, fp = 0, fn = 0 }
    p.n, p.exact = p.n + row.n, p.exact + row.diffs_ok
    p.tp, p.fp, p.fn = p.tp + row.tp, p.fp + row.fp, p.fn + row.fn
    pool[row.rate] = p
  end
  -- error-free reads must come back with no differences at all
  H.eq(pool[0].n, pool[0].exact, "perfect reads report no differences")
  H.eq(0, pool[0].fp)
  local floors = { [0.005] = 0.98, [0.02] = 0.95, [0.05] = 0.88 }
  for rate, floor in pairs(floors) do
    local p = assert(pool[rate], ("no reads at rate %s"):format(rate))
    local recall = p.tp / (p.tp + p.fn)
    local precision = p.tp / (p.tp + p.fp)
    H.ok(recall >= floor, ("rate %.3f: recall %.4f below %.2f"):format(rate, recall, floor))
    H.ok(precision >= floor, ("rate %.3f: precision %.4f below %.2f"):format(rate, precision, floor))
  end
  local rates = {}
  for rate in pairs(pool) do rates[#rates + 1] = rate end
  table.sort(rates)
  for _, rate in ipairs(rates) do
    local p = pool[rate]
    H.log(("[oracle] reads at %.1f%% substitutions: %d/%d reads exact, %d of %d differences "
      .. "recovered, %d spurious"):format(rate * 100, p.exact, p.n, p.tp, p.tp + p.fn, p.fp))
  end
end)

return T
