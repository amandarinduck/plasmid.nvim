-- Simulated sequencing reads, and what the engine recovers from them.
--
-- Reads are cut from the test plasmids at varied length, mutated at a known
-- error rate, reverse-complemented half the time, and (on circular records)
-- allowed to run through the origin. For each read we know exactly which
-- differences were introduced, so A.place and A.differences can be scored
-- against the truth rather than against themselves.
--
-- The random numbers come from a Lehmer generator written out here rather
-- than math.random, so the corpus is the same on every platform and the
-- thresholds in tests/spec/align_oracle_spec.lua mean something.
--
-- Used by tests/spec/align_oracle_spec.lua and tests/tools/align_report.lua.
local R = {}

local A = require("plasmid.bio.align")
local alphabet = require("plasmid.model.alphabet")

---------------------------------------------------------------------------
-- Deterministic RNG (MINSTD)
---------------------------------------------------------------------------

function R.rng(seed)
  local s = seed % 2147483647
  if s <= 0 then s = s + 2147483646 end
  return function(n)
    s = (s * 48271) % 2147483647
    local u = (s - 1) / 2147483646
    if n then return math.floor(u * n) + 1 end
    return u
  end
end

function R.read(path)
  local fh = assert(io.open(path, "rb"))
  local s = fh:read("*a")
  fh:close()
  return s
end

---------------------------------------------------------------------------
-- Mutation
---------------------------------------------------------------------------

local BASES = { "A", "C", "G", "T" }

--- Mutate `s`, returning the new string. Rates are per base.
function R.mutate(s, sub, ins, del, rnd)
  rnd = rnd or R.rng(1)
  local out = {}
  for i = 1, #s do
    local c = s:sub(i, i)
    local r = rnd()
    if r < del then
      -- dropped
    elseif r < del + ins then
      out[#out + 1] = BASES[rnd(4)]
      out[#out + 1] = c
    elseif r < del + ins + sub then
      local b = BASES[rnd(4)]
      while b == c do b = BASES[rnd(4)] end
      out[#out + 1] = b
    else
      out[#out + 1] = c
    end
  end
  return table.concat(out)
end

--- Cut a read out of `seq` and mutate it, recording the truth.
--- `start` is 1-based; on a circular record the window may pass the origin.
--- Returns read, truth (atomic differences in reference coordinates), wraps.
function R.simulate(seq, start, len, rates, rnd, circular)
  local L = #seq
  local out, truth = {}, {}
  local wraps = circular and (start + len - 1 > L) or false
  for w = 0, len - 1 do
    local p = circular and ((start - 1 + w) % L + 1) or (start + w)
    local c = seq:sub(p, p)
    local r = rnd()
    if r < rates.del then
      truth[#truth + 1] = { kind = "del", pos = p, ref = c }
    elseif r < rates.del + rates.ins then
      local b = BASES[rnd(4)]
      out[#out + 1] = b
      out[#out + 1] = c
      -- A.differences anchors an insertion on the base it follows
      truth[#truth + 1] = { kind = "ins", pos = (p - 2) % L + 1, query = b }
    elseif r < rates.del + rates.ins + rates.sub then
      local b = BASES[rnd(4)]
      while b == c do b = BASES[rnd(4)] end
      out[#out + 1] = b
      truth[#truth + 1] = { kind = "sub", pos = p, ref = c, query = b }
    else
      out[#out + 1] = c
    end
  end
  return table.concat(out), truth, wraps
end

---------------------------------------------------------------------------
-- Comparing recovered differences with the truth
---------------------------------------------------------------------------

--- Split the engine's run-length differences into one entry per base.
function R.atomize(diffs)
  local out = {}
  for _, d in ipairs(diffs) do
    if d.kind == "sub" then
      for k = 0, d.len - 1 do
        out[#out + 1] = { kind = "sub", pos = d.pos + k,
          ref = d.ref:sub(k + 1, k + 1), query = d.query:sub(k + 1, k + 1) }
      end
    elseif d.kind == "del" then
      for k = 0, d.len - 1 do
        out[#out + 1] = { kind = "del", pos = d.pos + k, ref = d.ref:sub(k + 1, k + 1) }
      end
    else
      for k = 0, d.len - 1 do
        out[#out + 1] = { kind = "ins", pos = d.pos, query = d.query:sub(k + 1, k + 1) }
      end
    end
  end
  return out
end

--- An indel inside a homopolymer run has several equally good placements.
--- Slide it to the leftmost one so the two lists can be compared directly.
local function canon(d, seq, L)
  if d.kind == "sub" then return ("sub:%d:%s"):format(d.pos, d.query) end
  local base = (d.kind == "del" and d.ref or d.query):upper()
  local p = d.pos
  if d.kind == "del" then
    -- deleting `base` at p is the same as deleting it at p-1 when the run
    -- extends that way
    local guard = 0
    while guard < L do
      local prev = (p - 2) % L + 1
      if seq:sub(prev, prev):upper() ~= base then break end
      p, guard = prev, guard + 1
    end
    return ("del:%d:%s"):format(p, base)
  end
  -- inserting `base` after p is the same as inserting it after p-1 when
  -- seq[p] is the same base
  local guard = 0
  while guard < L do
    if seq:sub(p, p):upper() ~= base then break end
    p, guard = (p - 2) % L + 1, guard + 1
  end
  return ("ins:%d:%s"):format(p, base)
end

local function multiset(list, seq, L)
  local t = {}
  for _, d in ipairs(list) do
    local k = canon(d, seq, L)
    t[k] = (t[k] or 0) + 1
  end
  return t
end

--- Compare two difference lists up to homopolymer slide. Returns whether
--- they are identical, and the recovered / spurious / missed counts.
function R.compare_diffs(got, truth, seq)
  local L = #seq
  local a, b = multiset(got, seq, L), multiset(truth, seq, L)
  local tp, fp, fn = 0, 0, 0
  for k, v in pairs(a) do
    local w = b[k] or 0
    tp = tp + math.min(v, w)
    if v > w then fp = fp + v - w end
  end
  for k, v in pairs(b) do
    local w = a[k] or 0
    if v > w then fn = fn + v - w end
  end
  return fp == 0 and fn == 0, tp, fp, fn
end

---------------------------------------------------------------------------
-- The experiment
---------------------------------------------------------------------------

R.RATES = {
  { sub = 0, ins = 0, del = 0 },
  { sub = 0.005, ins = 0.001, del = 0.001 },
  { sub = 0.02, ins = 0.004, del = 0.004 },
  { sub = 0.05, ins = 0.01, del = 0.01 },
}

R.PLASMIDS = { "J01749", "L09137", "NC_001422" }

--- Cut, place and score `per` reads for every plasmid and error rate.
--- Returns { per_rate = {...}, total = {...}, failures = {...} }.
function R.run(root, per, seed)
  local pio = require("plasmid.io")
  per = per or 24
  local rnd = R.rng(seed or 90210)
  local rows, failures = {}, {}
  local total = { n = 0, placed = 0, strand_ok = 0, wrap_ok = 0, pos_ok = 0, diffs_ok = 0,
    wrapping = 0, wrap_detected = 0, tp = 0, fp = 0, fn = 0 }
  for _, stem in ipairs(R.PLASMIDS) do
    local rec = assert(pio.read(root .. "/tests/fixtures/ncbi/" .. stem .. ".gb")).records[1]
    local seq = rec.seq:tostring():upper()
    local L = #seq
    local doubled = seq .. seq
    for _, rates in ipairs(R.RATES) do
      local row = { plasmid = stem, rate = rates.sub, n = 0, placed = 0, strand_ok = 0,
        wrap_ok = 0, pos_ok = 0, diffs_ok = 0, errs = {}, tp = 0, fp = 0, fn = 0 }
      for k = 1, per do
        local len = 200 + rnd(1001) - 1
        -- every fourth read is put across the origin
        local start = (k % 4 == 0) and (L - rnd(len - 20)) or rnd(L)
        local read, truth, wraps = R.simulate(seq, start, len, rates, rnd, true)
        local rc = (k % 2 == 0)
        local given = rc and alphabet.revcomp(read) or read
        row.n = row.n + 1
        total.n = total.n + 1
        if wraps then total.wrapping = total.wrapping + 1 end
        local aln = A.place(seq, given, { circular = true, mode = "semi" })
        if not aln then
          failures[#failures + 1] = ("%s rate %.3f read %d (%d bp at %d): not placed")
            :format(stem, rates.sub, k, len, start)
        else
          row.placed = row.placed + 1
          total.placed = total.placed + 1
          local want_strand = rc and -1 or 1
          if aln.strand == want_strand then
            row.strand_ok = row.strand_ok + 1
            total.strand_ok = total.strand_ok + 1
          end
          if aln.wraps == wraps then
            row.wrap_ok = row.wrap_ok + 1
            total.wrap_ok = total.wrap_ok + 1
            if wraps then total.wrap_detected = total.wrap_detected + 1 end
          else
            failures[#failures + 1] = ("%s rate %.3f read %d (%d bp at %d, %s): wraps %s, got %s")
              :format(stem, rates.sub, k, len, start, rc and "rc" or "fwd",
                tostring(wraps), tostring(aln.wraps))
          end
          local err = math.abs(aln.ref_start - start)
          err = math.min(err, L - err)
          row.errs[#row.errs + 1] = err
          if err <= 2 then
            row.pos_ok = row.pos_ok + 1
            total.pos_ok = total.pos_ok + 1
          else
            failures[#failures + 1] = ("%s rate %.3f read %d (%d bp at %d, %s): placed at %d, off by %d")
              :format(stem, rates.sub, k, len, start, rc and "rc" or "fwd", aln.ref_start, err)
          end
          -- differences, in the plasmid's own coordinates
          local q = aln.query_rc or given
          local got = A.differences(doubled, q,
            { ref_start = aln.ref_start_abs or aln.ref_start, query_start = aln.query_start,
              ops = aln.ops })
          for _, d in ipairs(got) do d.pos = (d.pos - 1) % L + 1 end
          local atomic = R.atomize(got)
          local same, tp, fp, fn = R.compare_diffs(atomic, truth, seq)
          row.tp, row.fp, row.fn = row.tp + tp, row.fp + fp, row.fn + fn
          total.tp, total.fp, total.fn = total.tp + tp, total.fp + fp, total.fn + fn
          if same then
            row.diffs_ok = row.diffs_ok + 1
            total.diffs_ok = total.diffs_ok + 1
          end
        end
      end
      table.sort(row.errs)
      row.median_err = row.errs[math.ceil(#row.errs / 2)] or -1
      rows[#rows + 1] = row
    end
  end
  return { per_rate = rows, total = total, failures = failures }
end

---------------------------------------------------------------------------
-- The banded path against the full matrix, on pairs that really do take it
---------------------------------------------------------------------------

--- Reads and whole plasmids aligned twice: once with the full matrix forced,
--- once with the band. Returns { n, same, worse, better, declined, worst_gap,
--- cases }. The band is a heuristic, so "same score" is the figure of merit.
function R.band_check(root, per, seed, maxlen)
  local pio = require("plasmid.io")
  per = per or 12
  maxlen = maxlen or 1500
  local rnd = R.rng(seed or 4242)
  local res = { n = 0, same = 0, worse = 0, declined = 0, worst_gap = 0, cases = {} }
  local function compare(label, ref, query)
    local fullaln = A.align(ref, query, { mode = "semi", max_cells = math.huge })
    local band = A.align(ref, query, { mode = "semi", max_cells = 1 })
    res.n = res.n + 1
    if not band then
      res.declined = res.declined + 1
      res.cases[#res.cases + 1] = label .. ": band declined"
    elseif not fullaln or band.score == fullaln.score then
      res.same = res.same + 1
    else
      res.worse = res.worse + 1
      local gap = fullaln.score - band.score
      if gap > res.worst_gap then res.worst_gap = gap end
      res.cases[#res.cases + 1] = ("%s: band %d, full %d"):format(label, band.score, fullaln.score)
    end
  end
  for _, stem in ipairs(R.PLASMIDS) do
    local seq = assert(pio.read(root .. "/tests/fixtures/ncbi/" .. stem .. ".gb"))
      .records[1].seq:tostring():upper()
    for k = 1, per do
      local len = 300 + rnd(maxlen - 299) - 1
      local at = rnd(#seq - len)
      local rates = R.RATES[(k % #R.RATES) + 1]
      local read = R.simulate(seq, at, len, rates, rnd, false)
      compare(("%s read %d bp @%d rate %.3f"):format(stem, len, at, rates.sub), seq, read)
    end
    -- whole molecule against a mutated copy
    compare(stem .. " vs mutated self", seq,
      R.mutate(seq, 0.01, 0.002, 0.002, rnd))
  end
  return res
end

return R
