-- Pairwise alignment with affine gaps (Gotoh), in three modes:
--
--   global  end gaps on both sequences are charged (Needleman-Wunsch)
--   local   the best-scoring segment of each sequence (Smith-Waterman)
--   semi    end gaps on the query are free: a read placed inside a reference
--           (set free_ref_ends to make the reference's ends free instead)
--
-- Scores follow the usual conventions: `match` / `mismatch` for nucleotides, a
-- substitution matrix (BLOSUM62) for proteins, `gap_open` charged once per gap
-- plus `gap_extend` per base in it, all as numbers where more is better.
--
-- Two code paths produce identical answers where both apply:
--
--   * the full matrix, when #ref * #query is within `max_cells`. Tracebacks
--     are kept as one byte per cell per state in Lua strings, so a 2 kb x 2 kb
--     alignment costs about 12 MB rather than the ~100 MB a table of numbers
--     would take.
--   * a band around the diagonals that shared k-mers vote for, otherwise. The
--     band costs O(#ref * width) and is what makes a read against a plasmid,
--     or two versions of one plasmid, practical.
local M = {}

local alphabet = require("plasmid.model.alphabet")

M.defaults = {
  match = 2,                -- nucleotide scores; blastn's defaults
  mismatch = -3,
  gap_open = -5,            -- charged once per gap, in addition to gap_extend
  gap_extend = -2,
  mode = "global",          -- "global" | "local" | "semi"
  free_ref_ends = true,     -- "semi": the reference's ends are free (a read inside a plasmid)
  free_query_ends = false,  -- "semi": free the query's ends instead (or as well)
  protein = nil,            -- nil: decided from the sequences
  ambiguous = "mismatch",   -- "mismatch" | "match": how N and friends score
  band = 24,                -- half-width of the banded path
  max_cells = 4e6,          -- full-matrix budget (#ref * #query)
  k = 11,                   -- seed length for the banded path
  max_positions = 64,       -- occurrences indexed per seed (repeat guard)
}

local function opts_for(o)
  return vim.tbl_extend("force", vim.deepcopy(M.defaults), o or {})
end

---------------------------------------------------------------------------
-- Scoring
---------------------------------------------------------------------------

local AMBIG = { A = true, C = true, G = true, T = true, U = true }

--- Scorer for two sequences: returns score(byte_a, byte_b) and whether the
--- pair counts as a match for identity.
local function scorer(o)
  if o.protein then
    local data = require("plasmid.data.blosum62")
    local index, scores = {}, {}
    for i = 1, #data.alphabet do index[data.alphabet:byte(i)] = i end
    for sym, row in pairs(data.scores) do
      local b = sym:byte(1)
      scores[b] = row
      index[sym:lower():byte(1)] = index[b]
    end
    local xrow = scores[string.byte("X")]
    return function(a, b)
      if a >= 97 then a = a - 32 end
      if b >= 97 then b = b - 32 end
      local row = scores[a] or xrow
      local col = index[b] or index[string.byte("X")]
      local s = row[col] or 0
      return s, a == b
    end
  end
  local match, mismatch = o.match, o.mismatch
  local ambiguous_match = o.ambiguous == "match"
  return function(a, b)
    if a >= 97 then a = a - 32 end
    if b >= 97 then b = b - 32 end
    if a == 85 then a = 84 end          -- U counts as T
    if b == 85 then b = 84 end
    if a == b then return match, true end
    if ambiguous_match then
      local ca, cb = string.char(a), string.char(b)
      if not AMBIG[ca] or not AMBIG[cb] then return match, true end
    end
    return mismatch, false
  end
end

local function looks_protein(s)
  local n, odd = 0, 0
  for i = 1, math.min(#s, 200) do
    local c = s:sub(i, i):upper()
    if c:match("%a") then
      n = n + 1
      if not c:match("[ACGTUN]") then odd = odd + 1 end
    end
  end
  return n > 0 and odd / n > 0.15
end

---------------------------------------------------------------------------
-- Full matrix (Gotoh)
---------------------------------------------------------------------------

-- Traceback codes: which state the cell was reached from. One byte per cell
-- and state, packed into a string per row.
local FROM_M, FROM_X, FROM_Y, STOP = 1, 2, 3, 0
local NEG = -1e15

local function pack_row(t, m)
  if m <= 4096 then return string.char(unpack(t, 1, m)) end
  local out = {}
  for a = 1, m, 4096 do out[#out + 1] = string.char(unpack(t, a, math.min(a + 4095, m))) end
  return table.concat(out)
end

--- Align `ref` against `query` over the whole matrix.
--- Returns { score, ops, ref_start, ref_stop, query_start, query_stop } or nil
--- (local mode with no positive-scoring segment).
local function full(ref, query, o)
  local n, m = #ref, #query
  local score_of = scorer(o)
  local go, ge = o.gap_open, o.gap_extend
  local mode = o.mode
  -- Free end gaps: "ref" means the reference may hang off the alignment at no
  -- cost (a read placed inside a plasmid), "query" the other way round.
  local free_ref = mode == "local" or (mode == "semi" and o.free_ref_ends ~= false)
  local free_query = mode == "local" or (mode == "semi" and o.free_query_ends)
  local byte = string.byte

  local prevM, prevX, prevY = {}, {}, {}
  local ptrM, ptrX, ptrY = {}, {}, {}
  prevM[0], prevX[0], prevY[0] = 0, NEG, NEG
  for j = 1, m do
    prevM[j] = NEG
    prevX[j] = NEG
    prevY[j] = free_query and 0 or (go + ge * j)   -- leading gap in ref = query overhang
  end
  -- best end cell, by mode
  local best, bi, bj, bstate = NEG, 0, 0, "M"
  if mode == "local" then best, bi, bj = 0, 0, 0 end
  -- row 0 with a free trailing reference: the whole query inserted ahead of a
  -- reference that then hangs off. Rarely optimal, but it is the answer for
  -- two sequences with nothing in common.
  if free_ref and mode ~= "local" and prevY[m] > best then
    best, bi, bj, bstate = prevY[m], 0, m, "Y"
  end

  local curM, curX, curY = {}, {}, {}
  local rowM, rowX, rowY = {}, {}, {}
  for i = 1, n do
    local rb = byte(ref, i)
    curM[0] = NEG
    curY[0] = NEG
    curX[0] = free_ref and 0 or (go + ge * i)      -- leading gap in query = ref overhang
    for j = 1, m do
      local s = score_of(rb, byte(query, j))
      local dM, dX, dY = prevM[j - 1], prevX[j - 1], prevY[j - 1]
      local bestd, from = dM, FROM_M
      if dX > bestd then bestd, from = dX, FROM_X end
      if dY > bestd then bestd, from = dY, FROM_Y end
      local mv = bestd + s
      if mode == "local" and mv < 0 then mv, from = 0, STOP end
      curM[j], rowM[j] = mv, from
      -- A gap may follow a gap of the other kind: with cheap opens or an
      -- expensive mismatch, "insert then delete" can beat aligning the two
      -- runs against each other, so X and Y see each other as well as M.
      local openx, extx, swx = prevM[j] + go + ge, prevX[j] + ge, prevY[j] + go + ge
      if openx >= extx and openx >= swx then curX[j], rowX[j] = openx, FROM_M
      elseif extx >= swx then curX[j], rowX[j] = extx, FROM_X
      else curX[j], rowX[j] = swx, FROM_Y end
      local openy, exty, swy = curM[j - 1] + go + ge, curY[j - 1] + ge, curX[j - 1] + go + ge
      if openy >= exty and openy >= swy then curY[j], rowY[j] = openy, FROM_M
      elseif exty >= swy then curY[j], rowY[j] = exty, FROM_Y
      else curY[j], rowY[j] = swy, FROM_X end
      if mode == "local" and curM[j] > best then best, bi, bj = curM[j], i, j end
    end
    ptrM[i], ptrX[i], ptrY[i] = pack_row(rowM, m), pack_row(rowX, m), pack_row(rowY, m)
    -- Free trailing reference overhang: the alignment may end in this row.
    -- A trailing query overhang is still charged, so the Y state (the query
    -- hanging off the end of the reference) is an end candidate too.
    if free_ref and mode ~= "local" then
      if curM[m] > best then best, bi, bj, bstate = curM[m], i, m, "M" end
      if curY[m] > best then best, bi, bj, bstate = curY[m], i, m, "Y" end
    end
    prevM, curM = curM, prevM
    prevX, curX = curX, prevX
    prevY, curY = curY, prevY
  end
  if mode == "global" or mode == "semi" then
    if free_query then
      -- Free trailing query overhang: the alignment may end anywhere in row
      -- n. Symmetrically, a trailing reference overhang is still charged, so
      -- the X state is a candidate; at j = 0 nothing of the query has been
      -- consumed, which is the answer when the two share nothing.
      for j = 0, m do
        if j >= 1 and prevM[j] > best then best, bi, bj, bstate = prevM[j], n, j, "M" end
        if prevX[j] > best then best, bi, bj, bstate = prevX[j], n, j, "X" end
      end
    end
    if not free_ref and not free_query then
      best, bi, bj, bstate = prevM[m], n, m, "M"
      if prevX[m] > best then best, bstate = prevX[m], "X" end
      if prevY[m] > best then best, bstate = prevY[m], "Y" end
    elseif best == NEG then
      best, bi, bj, bstate = prevM[m], n, m, "M"
    end
  end
  if mode == "local" and best <= 0 then return nil end

  -- traceback
  local ops = {}
  local function push(op)
    local last = ops[#ops]
    if last and last.op == op then last.n = last.n + 1 else ops[#ops + 1] = { op = op, n = 1 } end
  end
  local i, j, state = bi, bj, bstate
  while i > 0 and j > 0 do
    if state == "M" then
      local from = ptrM[i]:byte(j)
      if mode == "local" and from == STOP then break end
      push("M")
      i, j = i - 1, j - 1
      state = (from == FROM_X and "X") or (from == FROM_Y and "Y") or "M"
    elseif state == "X" then
      local from = ptrX[i]:byte(j)
      push("D")
      i = i - 1
      state = (from == FROM_X and "X") or (from == FROM_Y and "Y") or "M"
    else
      local from = ptrY[i]:byte(j)
      push("I")
      j = j - 1
      state = (from == FROM_Y and "Y") or (from == FROM_X and "X") or "M"
    end
  end
  -- leading overhangs: charged already unless they are free
  if not free_query then while j > 0 do push("I"); j = j - 1 end end
  if not free_ref then while i > 0 do push("D"); i = i - 1 end end
  local rev = {}
  for k = #ops, 1, -1 do rev[#rev + 1] = ops[k] end
  return { score = best, ops = rev, ref_start = i + 1, ref_stop = bi,
    query_start = j + 1, query_stop = bj }
end

---------------------------------------------------------------------------
-- Banded (seeded) alignment
---------------------------------------------------------------------------

-- Shared k-mers vote for diagonals d = j - i; the band is w cells either side
-- of the winner.
local function seed_diagonal(ref, query, o)
  local k, maxpos = o.k, o.max_positions
  local n, m = #ref, #query
  if n < k or m < k then return nil end
  local idx = {}
  local up = query:upper():gsub("U", "T")
  for j = 1, m - k + 1 do
    local key = up:sub(j, j + k - 1)
    if not key:find("[^ACGT]") then
      local t = idx[key]
      if not t then idx[key] = { j }
      elseif #t < maxpos then t[#t + 1] = j end
    end
  end
  local votes, best, bestd = {}, 0, nil
  local rup = ref:upper():gsub("U", "T")
  for i = 1, n - k + 1 do
    local t = idx[rup:sub(i, i + k - 1)]
    if t then
      for _, j in ipairs(t) do
        local d = j - i
        local v = (votes[d] or 0) + 1
        votes[d] = v
        if v > best then best, bestd = v, d end
      end
    end
  end
  if not bestd then return nil end
  -- Group the well-voted diagonals into clusters -- runs of diagonals no
  -- more than a band apart, which is what one alignment with indels looks
  -- like -- and take the cluster with the most votes in total.
  --
  -- The single best diagonal is not enough to choose by. A circular
  -- reference is doubled before placing a read, so a read that crosses the
  -- origin matches the start of the first copy just as well as it matches
  -- the junction, and the two peaks tie: only the bases before the origin,
  -- which vote for the junction alone, break it. Merging the two clusters
  -- instead would make the band span the whole molecule.
  local good = {}
  for d, v in pairs(votes) do
    if v >= best * 0.25 then good[#good + 1] = d end
  end
  table.sort(good)
  local lo, hi, mass, peak = nil, nil, -1, 0
  local i = 1
  while i <= #good do
    local j = i
    while j < #good and good[j + 1] - good[j] <= o.band do j = j + 1 end
    local total, top, topd = 0, 0, good[i]
    for c = i, j do
      local v = votes[good[c]]
      total = total + v
      if v > top then top, topd = v, good[c] end
    end
    if total > mass then lo, hi, mass, peak, bestd = good[i], good[j], total, top, topd end
    i = j + 1
  end
  return bestd, lo, hi, peak
end

--- Gotoh inside a band of diagonals [dlo - w, dhi + w].
local function banded(ref, query, o, dlo, dhi)
  local n, m = #ref, #query
  local w = o.band
  local d0, d1 = dlo - w, dhi + w
  local width = d1 - d0 + 1
  if width >= m then return full(ref, query, o) end
  local score_of = scorer(o)
  local go, ge = o.gap_open, o.gap_extend
  local mode = o.mode
  local free_ref = mode == "local" or (mode == "semi" and o.free_ref_ends ~= false)
  local free_query = mode == "local" or (mode == "semi" and o.free_query_ends)
  local byte = string.byte

  -- cell (i, c) holds column j = i + d0 + c - 1
  local prevM, prevX, prevY = {}, {}, {}
  local ptrM, ptrX, ptrY = {}, {}, {}
  for c = 1, width do
    local j = d0 + c - 1          -- row 0
    if j < 0 then prevM[c], prevX[c], prevY[c] = NEG, NEG, NEG
    elseif j == 0 then prevM[c], prevX[c], prevY[c] = 0, NEG, NEG
    else
      prevM[c], prevX[c] = NEG, NEG
      prevY[c] = free_query and 0 or (go + ge * j)
    end
  end
  local best, bi, bc, bstate = NEG, 0, 0, "M"
  if mode == "local" then best = 0 end
  if free_ref and mode ~= "local" then          -- row 0, as in `full`
    local c = m - d0 + 1
    if c >= 1 and c <= width and prevY[c] > best then
      best, bi, bc, bstate = prevY[c], 0, c, "Y"
    end
  end
  local curM, curX, curY = {}, {}, {}
  local rowM, rowX, rowY = {}, {}, {}
  for i = 1, n do
    local rb = byte(ref, i)
    for c = 1, width do
      local j = i + d0 + c - 1
      local mv, from = NEG, FROM_M
      local xv, xfrom = NEG, FROM_M
      local yv, yfrom = NEG, FROM_M
      if j == 0 then
        -- no query consumed yet: a run of reference-only columns
        xv, xfrom = free_ref and 0 or (go + ge * i), FROM_X
      elseif j >= 1 and j <= m then
        local s = score_of(rb, byte(query, j))
        -- diagonal: same c in the previous row
        local dM, dX, dY = prevM[c], prevX[c], prevY[c]
        local bestd = dM
        if dX > bestd then bestd, from = dX, FROM_X end
        if dY > bestd then bestd, from = dY, FROM_Y end
        if bestd > NEG / 2 then mv = bestd + s end
        if mode == "local" and mv < 0 then mv, from = 0, STOP end
        -- gap in query (consume ref): previous row, c + 1
        if c < width then
          local openx, extx = prevM[c + 1] + go + ge, prevX[c + 1] + ge
          local swx = prevY[c + 1] + go + ge
          if openx >= extx and openx >= swx then xv, xfrom = openx, FROM_M
          elseif extx >= swx then xv, xfrom = extx, FROM_X
          else xv, xfrom = swx, FROM_Y end
        end
        -- gap in ref (consume query): same row, c - 1
        if c > 1 then
          local openy, exty = curM[c - 1] + go + ge, curY[c - 1] + ge
          local swy = curX[c - 1] + go + ge
          if openy >= exty and openy >= swy then yv, yfrom = openy, FROM_M
          elseif exty >= swy then yv, yfrom = exty, FROM_Y
          else yv, yfrom = swy, FROM_X end
        end
      end
      curM[c], rowM[c] = mv, from
      curX[c], rowX[c] = xv, xfrom
      curY[c], rowY[c] = yv, yfrom
      if mode == "local" and mv > best then best, bi, bc, bstate = mv, i, c, "M" end
    end
    ptrM[i], ptrX[i], ptrY[i] = pack_row(rowM, width), pack_row(rowX, width), pack_row(rowY, width)
    if free_ref and mode ~= "local" then
      local c = m - i - d0 + 1
      if c >= 1 and c <= width then
        if curM[c] > best then best, bi, bc, bstate = curM[c], i, c, "M" end
        if curY[c] > best then best, bi, bc, bstate = curY[c], i, c, "Y" end
      end
    end
    prevM, curM = curM, prevM
    prevX, curX = curX, prevX
    prevY, curY = curY, prevY
  end
  if mode ~= "local" then
    if free_query then
      for c = 1, width do
        local j = n + d0 + c - 1
        if j >= 1 and j <= m and prevM[c] > best then best, bi, bc, bstate = prevM[c], n, c, "M" end
        if j >= 0 and j <= m and prevX[c] > best then best, bi, bc, bstate = prevX[c], n, c, "X" end
      end
    end
    if not free_ref and not free_query then
      local c = m - n - d0 + 1
      if c < 1 or c > width then return nil end   -- the band does not reach the corner
      best, bi, bc, bstate = prevM[c], n, c, "M"
      if prevX[c] > best then best, bstate = prevX[c], "X" end
      if prevY[c] > best then best, bstate = prevY[c], "Y" end
    end
  end
  if best <= NEG / 2 or (mode == "local" and best <= 0) then return nil end

  local ops = {}
  local function push(op)
    local last = ops[#ops]
    if last and last.op == op then last.n = last.n + 1 else ops[#ops + 1] = { op = op, n = 1 } end
  end
  local i, c, state = bi, bc, bstate
  local j = i + d0 + c - 1
  while i > 0 and j > 0 do
    if c < 1 or c > width then break end
    if state == "M" then
      local from = ptrM[i]:byte(c)
      if from == STOP then break end
      push("M")
      i, j = i - 1, j - 1
      state = (from == FROM_X and "X") or (from == FROM_Y and "Y") or "M"
    elseif state == "X" then
      local from = ptrX[i]:byte(c)
      push("D")
      i, c = i - 1, c + 1
      state = (from == FROM_X and "X") or (from == FROM_Y and "Y") or "M"
    else
      local from = ptrY[i]:byte(c)
      push("I")
      j, c = j - 1, c - 1
      state = (from == FROM_Y and "Y") or (from == FROM_X and "X") or "M"
    end
  end
  if not free_query then while j > 0 do push("I"); j = j - 1 end end
  if not free_ref then while i > 0 do push("D"); i = i - 1 end end
  local rev = {}
  for k = #ops, 1, -1 do rev[#rev + 1] = ops[k] end
  return { score = best, ops = rev, ref_start = i + 1, ref_stop = bi,
    query_start = j + 1, query_stop = bi + d0 + bc - 1, banded = true }
end

---------------------------------------------------------------------------
-- Public API
---------------------------------------------------------------------------

-- A local alignment never begins or ends with a gap, and a leading or
-- trailing run that contributes nothing (a mismatch cancelled by the match
-- next to it, where the two cost the same) is no part of the segment either.
-- Both are co-optimal with the trimmed form; this is the one every other
-- implementation prints, so trim before reporting.
local function trim_local(ref, query, aln, o)
  local score_of = scorer(o)
  local ops = aln.ops
  local i, j = aln.ref_start, aln.query_start
  -- per column: its score, and how much of each sequence it consumes
  local cols = {}
  for _, op in ipairs(ops) do
    if op.op == "M" then
      for k = 0, op.n - 1 do
        cols[#cols + 1] = { s = (score_of(ref:byte(i + k), query:byte(j + k))), r = 1, q = 1, op = "M" }
      end
      i, j = i + op.n, j + op.n
    elseif op.op == "D" then
      for k = 1, op.n do
        cols[#cols + 1] = { s = (k == 1 and o.gap_open or 0) + o.gap_extend, r = 1, q = 0, op = "D" }
      end
      i = i + op.n
    else
      for k = 1, op.n do
        cols[#cols + 1] = { s = (k == 1 and o.gap_open or 0) + o.gap_extend, r = 0, q = 1, op = "I" }
      end
      j = j + op.n
    end
  end
  local from, to = 1, #cols
  local run, best_cut = 0, 0
  for k = from, to do
    run = run + cols[k].s
    if run <= 0 then best_cut = k end
  end
  from = best_cut + 1
  run, best_cut = 0, to + 1
  for k = to, from, -1 do
    run = run + cols[k].s
    if run <= 0 then best_cut = k end
  end
  to = best_cut - 1
  if from == 1 and to == #cols then return aln end
  if from > to then return nil end
  local dr, dq = 0, 0
  for k = 1, from - 1 do dr, dq = dr + cols[k].r, dq + cols[k].q end
  local er, eq = 0, 0
  for k = to + 1, #cols do er, eq = er + cols[k].r, eq + cols[k].q end
  local trimmed = {}
  for k = from, to do
    local last = trimmed[#trimmed]
    if last and last.op == cols[k].op then last.n = last.n + 1
    else trimmed[#trimmed + 1] = { op = cols[k].op, n = 1 } end
  end
  aln.ops = trimmed
  aln.ref_start, aln.query_start = aln.ref_start + dr, aln.query_start + dq
  aln.ref_stop, aln.query_stop = aln.ref_stop - er, aln.query_stop - eq
  return aln
end

--- Statistics and the aligned strings of an alignment.
function M.annotate(ref, query, aln, o)
  o = opts_for(o)
  o.protein = o.protein
  if o.protein == nil then o.protein = looks_protein(ref) and looks_protein(query) end
  local score_of = scorer(o)
  local rparts, qparts, mid = {}, {}, {}
  local i, j = aln.ref_start, aln.query_start
  local matches, mismatches, gap_bases, gaps, columns = 0, 0, 0, 0, 0
  for _, op in ipairs(aln.ops) do
    if op.op == "M" then
      for _ = 1, op.n do
        local a, b = ref:sub(i, i), query:sub(j, j)
        local _, same = score_of(a:byte(), b:byte())
        rparts[#rparts + 1], qparts[#qparts + 1] = a, b
        mid[#mid + 1] = same and "|" or " "
        if same then matches = matches + 1 else mismatches = mismatches + 1 end
        i, j, columns = i + 1, j + 1, columns + 1
      end
    elseif op.op == "D" then
      gaps = gaps + 1
      for _ = 1, op.n do
        rparts[#rparts + 1], qparts[#qparts + 1], mid[#mid + 1] = ref:sub(i, i), "-", " "
        i, columns, gap_bases = i + 1, columns + 1, gap_bases + 1
      end
    else
      gaps = gaps + 1
      for _ = 1, op.n do
        rparts[#rparts + 1], qparts[#qparts + 1], mid[#mid + 1] = "-", query:sub(j, j), " "
        j, columns, gap_bases = j + 1, columns + 1, gap_bases + 1
      end
    end
  end
  aln.ref_aln, aln.query_aln, aln.mid = table.concat(rparts), table.concat(qparts), table.concat(mid)
  aln.matches, aln.mismatches, aln.gaps, aln.gap_bases, aln.columns =
    matches, mismatches, gaps, gap_bases, columns
  aln.identity = columns > 0 and matches / columns or 0
  aln.ref_span = aln.ref_stop - aln.ref_start + 1
  aln.query_span = aln.query_stop - aln.query_start + 1
  aln.ref_coverage = #ref > 0 and aln.ref_span / #ref or 0
  aln.query_coverage = #query > 0 and aln.query_span / #query or 0
  return aln
end

--- Align `query` to `ref`. opts as M.defaults; `mode` is "global", "local" or
--- "semi" (free reference ends, i.e. the query placed inside the reference;
--- free_query_ends = true frees the query's ends instead).
--- Returns the alignment, or nil plus a reason.
function M.align(ref, query, opts)
  local o = opts_for(opts)
  if #ref == 0 or #query == 0 then return nil, "empty sequence" end
  if o.protein == nil then o.protein = looks_protein(ref) and looks_protein(query) end
  local aln, why
  if #ref * #query <= o.max_cells then
    aln = full(ref, query, o)
    if not aln then return nil, "no alignment" end
  else
    local _, lo, hi, votes = seed_diagonal(ref, query, o)
    if not lo then
      return nil, ("the sequences are too long for the full matrix (%d x %d) and share no %d-mer")
        :format(#ref, #query, o.k)
    end
    if (hi - lo) + 2 * o.band > #query then
      return nil, "the shared k-mers are spread over too many diagonals to band"
    end
    aln = banded(ref, query, o, lo, hi)
    if not aln then
      return nil, ("no alignment within %d bases of the seeded diagonals"):format(o.band)
    end
    aln.seed_votes = votes
  end
  if o.mode == "local" then
    aln = trim_local(ref, query, aln, o)
    if not aln then return nil, "no alignment" end
  end
  aln.mode = o.mode
  aln.strand = 1
  return M.annotate(ref, query, aln, o)
end

--- The better of the two strands. Returns the alignment with `strand` set to
--- 1 or -1; a -1 alignment's query coordinates refer to the reverse
--- complement of the query that was passed in.
function M.align_strands(ref, query, opts)
  local fwd, why = M.align(ref, query, opts)
  local rc = alphabet.revcomp(query)
  local rev = M.align(ref, rc, opts)
  if not fwd and not rev then return nil, why end
  if rev and (not fwd or rev.score > fwd.score) then
    rev.strand = -1
    rev.query_rc = rc
    return rev
  end
  return fwd
end

---------------------------------------------------------------------------
-- Circular references
---------------------------------------------------------------------------

local function wrap_pos(p, L) return (p - 1) % L + 1 end

--- Place a query on a reference: the better strand, and on a circular
--- reference an alignment that may cross the origin. opts.circular selects
--- that; opts.mode defaults to "semi" (the query sits inside the reference).
--- The result carries strand, wraps, and ref_start/ref_stop in 1..#ref
--- (ref_start > ref_stop when the alignment crosses the origin).
function M.place(ref, query, opts)
  local mode = opts and opts.mode or "semi"
  local o = opts_for(opts)
  o.mode = mode
  if not o.circular then
    local aln, why = M.align_strands(ref, query, o)
    if not aln then return nil, why end
    aln.wraps = false
    return aln
  end
  local L = #ref
  -- Doubling the reference lets an alignment run through the origin; only as
  -- much as the query can reach is added.
  local extra = math.min(L, #query + o.band)
  local doubled = ref .. ref:sub(1, extra)
  local aln, why = M.align_strands(doubled, query, o)
  if not aln then return nil, why end
  if aln.ref_span > L then return nil, "the query is longer than the circular reference" end
  aln.wraps = aln.ref_stop > L
  aln.ref_start_abs, aln.ref_stop_abs = aln.ref_start, aln.ref_stop
  aln.ref_start, aln.ref_stop = wrap_pos(aln.ref_start, L), wrap_pos(aln.ref_stop, L)
  aln.ref_length = L
  return aln
end

--- Compare two sequences end to end. On circular records the second sequence
--- is rotated into phase with the first (and reverse-complemented if that
--- scores better) so a different origin is not reported as a difference.
--- Returns the alignment with `rotation` (bases the query was rotated by) and
--- `strand`.
function M.compare(a, b, opts)
  local mode = opts and opts.mode or "global"
  local o = opts_for(opts)
  o.mode = mode
  if not o.circular then
    local aln, why = M.align_strands(a, b, o)
    if not aln then return nil, why end
    aln.rotation = 0
    return aln
  end
  local best
  for _, cand in ipairs({ { b, 1 }, { alphabet.revcomp(b), -1 } }) do
    local seq, strand = cand[1], cand[2]
    local d = seed_diagonal(a, seq, o)
    local rot = d and ((d % #seq) + #seq) % #seq or 0
    local rotated = rot == 0 and seq or (seq:sub(rot + 1) .. seq:sub(1, rot))
    local aln = M.align(a, rotated, o)
    if aln and (not best or aln.score > best.score) then
      aln.strand = strand
      aln.rotation = rot
      aln.query_rotated = rotated
      aln.query_oriented = seq
      best = aln
    end
  end
  if not best then return nil, "no alignment" end
  -- Both molecules are circles, so cutting them open anywhere can leave an
  -- artificial gap at the ends. Cut again in the middle of the longest run of
  -- matched columns, where the two agree, and align from there.
  local run, run_ref, run_query = 0, nil, nil
  local i, j = best.ref_start, best.query_start
  for _, op in ipairs(best.ops) do
    if op.op == "M" and op.n > run then run, run_ref, run_query = op.n, i, j end
    if op.op == "M" then i, j = i + op.n, j + op.n
    elseif op.op == "D" then i = i + op.n
    else j = j + op.n end
  end
  if run >= 2 * o.k and run_ref then
    local La, Lb = #a, #best.query_rotated
    local ra = (run_ref + math.floor(run / 2) - 1) % La
    local rb = (run_query + math.floor(run / 2) - 1) % Lb
    local a2 = a:sub(ra + 1) .. a:sub(1, ra)
    local b2 = best.query_rotated:sub(rb + 1) .. best.query_rotated:sub(1, rb)
    local aln2 = M.align(a2, b2, o)
    if aln2 and aln2.score >= best.score then
      aln2.strand = best.strand
      aln2.rotation = (best.rotation + rb) % Lb
      aln2.query_rotated = b2
      aln2.query_oriented = best.query_oriented
      aln2.ref_rotation = ra
      aln2.ref_rotated = a2
      aln2.ref_length = La
      best = aln2
    end
  end
  best.ref_rotation = best.ref_rotation or 0
  best.ref_rotated = best.ref_rotated or a
  best.ref_length = best.ref_length or #a
  return best
end

--- Compare two sequences and report the differences in the first sequence's
--- own coordinates (M.compare may cut a circular reference elsewhere).
--- Returns the alignment and the difference list.
function M.compare_differences(a, b, opts)
  local aln, why = M.compare(a, b, opts)
  if not aln then return nil, why end
  local diffs = M.differences(aln.ref_rotated, aln.query_rotated, aln)
  local L, rot = aln.ref_length, aln.ref_rotation
  if rot ~= 0 then
    for _, d in ipairs(diffs) do d.pos = (d.pos - 1 + rot) % L + 1 end
    table.sort(diffs, function(x, y) return x.pos < y.pos end)
  end
  return aln, diffs
end

--- Differences of an aligned query from the reference, in reference
--- coordinates: { kind = "sub"|"ins"|"del", pos, ref, query, len }.
--- `pos` is the first reference base involved; for an insertion it is the base
--- the inserted sequence follows.
function M.differences(ref, query, aln)
  local out = {}
  local i, j = aln.ref_start, aln.query_start
  for _, op in ipairs(aln.ops) do
    if op.op == "M" then
      local run
      for k = 0, op.n - 1 do
        local a, b = ref:sub(i + k, i + k), query:sub(j + k, j + k)
        if a:upper() ~= b:upper() then
          if run and run.pos + #run.ref == i + k then
            run.ref, run.query = run.ref .. a, run.query .. b
            run.len = run.len + 1
          else
            run = { kind = "sub", pos = i + k, ref = a, query = b, len = 1 }
            out[#out + 1] = run
          end
        else
          run = nil
        end
      end
      i, j = i + op.n, j + op.n
    elseif op.op == "D" then
      out[#out + 1] = { kind = "del", pos = i, ref = ref:sub(i, i + op.n - 1), query = "", len = op.n }
      i = i + op.n
    else
      out[#out + 1] = { kind = "ins", pos = i - 1, ref = "", query = query:sub(j, j + op.n - 1), len = op.n }
      j = j + op.n
    end
  end
  return out
end

--- CIGAR string of an alignment ("120M2I35M1D900M").
function M.cigar(aln)
  local parts = {}
  for _, op in ipairs(aln.ops) do parts[#parts + 1] = ("%d%s"):format(op.n, op.op) end
  return table.concat(parts)
end

M._full = full
M._banded = banded
M._scorer = scorer

return M
