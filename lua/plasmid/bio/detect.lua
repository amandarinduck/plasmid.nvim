-- Feature detection: find reference features (a library of named sequences)
-- in a query sequence. Pure Lua.
--
--   1. seeds: every 12-mer of the query (both strands; circular queries are
--      extended past the origin) is indexed as a 24-bit integer
--   2. candidates: each reference's 12-mers vote for diagonals; nearby
--      diagonals are clustered
--   3. extension: banded Smith-Waterman around the best diagonal of each
--      cluster gives identity, reference coverage and exact coordinates,
--      tolerating small indels
--   4. filtering: identity / coverage thresholds, then redundant hits on the
--      same span (database variants of one part) are resolved by score while
--      genuinely nested features (an operator inside a promoter) are kept
local alphabet = require("plasmid.model.alphabet")

local M = {}

local K = 12
local CODE = {}
for ch, v in pairs({ A = 0, C = 1, G = 2, T = 3, U = 3 }) do
  CODE[ch:byte()] = v
  CODE[ch:lower():byte()] = v
end
local MASK = 4 ^ K

M.defaults = {
  min_identity = 0.95,     -- fraction of aligned columns that match
  min_coverage = 0.95,     -- fraction of the reference covered for a full-length hit
  fragments = true,        -- also report partial matches...
  min_fragment_coverage = 0.5,
  min_fragment_length = 50,
  band = 16,               -- indel tolerance (bp) of the banded alignment
  max_positions = 32,      -- k-mer occurrences indexed per query k-mer (repeat guard)
  -- Drop a hit lying inside a higher-scoring hit of the same feature type
  -- (CMVd1 promoter inside CMV promoter, ISS inside AmpR): such hits are
  -- sub-parts of the larger feature. Primers, operators etc. inside other
  -- types are kept.
  nested_same_type = false,
}

--- Index k-mers of s: kmer -> { positions }.
local function index_kmers(s, max_pos)
  local idx = {}
  local v, valid = 0, 0
  local byte = string.byte
  for i = 1, #s do
    local c = CODE[byte(s, i)]
    if c then
      v = (v * 4 + c) % MASK
      valid = valid + 1
    else
      valid = 0
    end
    if valid >= K then
      local p = i - K + 1
      local t = idx[v]
      if not t then
        idx[v] = { p }
      elseif #t < max_pos then
        t[#t + 1] = p
      end
    end
  end
  return idx
end

--- Diagonal votes of a reference against an indexed query: diag = qpos - rpos.
--- Also returns, per diagonal, the ref span covered by its seeds.
local function diagonal_votes(ref, qidx)
  local votes, lo, hi = {}, {}, {}
  local v, valid = 0, 0
  local byte = string.byte
  for i = 1, #ref do
    local c = CODE[byte(ref, i)]
    if c then
      v = (v * 4 + c) % MASK
      valid = valid + 1
    else
      valid = 0
    end
    if valid >= K then
      local t = qidx[v]
      if t then
        local rp = i - K + 1
        for _, qp in ipairs(t) do
          local d = qp - rp
          local c = votes[d]
          if c then
            votes[d] = c + 1
            hi[d] = rp
          else
            votes[d], lo[d], hi[d] = 1, rp, rp
          end
        end
      end
    end
  end
  return votes, lo, hi
end

local MATCH, MISMATCH, GAP = 2, -3, -5

local function eq(a, b)
  if a == b then return true end
  -- case-insensitive, U == T
  if a >= 97 then a = a - 32 end
  if b >= 97 then b = b - 32 end
  if a == 85 then a = 84 end
  if b == 85 then b = 84 end
  return a == b
end

--- Banded local alignment of ref rows r_from..r_to against query positions
--- i + d + k, k in [-w, w]. Returns best { score, r0, r1, q0, q1, matches, columns } or nil.
local function banded_sw(ref, q, d, w, r_from, r_to)
  local n = #q
  r_from, r_to = r_from or 1, r_to or #ref
  local width = 2 * w + 1
  local byte = string.byte
  local prev, cur = {}, {}
  for k = 1, width do prev[k] = 0 end
  local ptr = {}
  local best, bi, bk = 0, 0, 0
  for i = r_from, r_to do
    local rb = byte(ref, i)
    local base = (i - 1) * width
    for kk = 1, width do
      local k = kk - w - 1
      local j = i + d + k
      local h, p = 0, 0
      if j >= 1 and j <= n then
        -- diagonal: (i-1, j-1) has the same k
        local s = prev[kk] + (eq(rb, byte(q, j)) and MATCH or MISMATCH)
        if s > h then h, p = s, 1 end
        -- up: gap in query, (i-1, j) -> k+1
        if kk < width then
          local u = prev[kk + 1] + GAP
          if u > h then h, p = u, 2 end
        end
        -- left: gap in ref, (i, j-1) -> k-1
        if kk > 1 then
          local l = cur[kk - 1] + GAP
          if l > h then h, p = l, 3 end
        end
      end
      cur[kk] = h
      ptr[base + kk] = p
      if h > best then best, bi, bk = h, i, kk end
    end
    prev, cur = cur, prev
  end
  if best <= 0 then return nil end
  -- traceback
  local i, kk = bi, bk
  local matches, cols = 0, 0
  local r1, q1 = bi, bi + d + (bk - w - 1)
  local r0, q0 = r1, q1
  while i >= r_from do
    local p = ptr[(i - 1) * width + kk]
    if not p or p == 0 then break end
    local j = i + d + (kk - w - 1)
    r0, q0 = i, j
    cols = cols + 1
    if p == 1 then
      if eq(byte(ref, i), byte(q, j)) then matches = matches + 1 end
      i = i - 1
    elseif p == 2 then
      i = i - 1
      kk = kk + 1
    else
      kk = kk - 1
    end
  end
  return { score = best, r0 = r0, r1 = r1, q0 = q0, q1 = q1, matches = matches, columns = cols }
end

--- Search one query strand string for one reference.
local function search_ref(ref, q, qidx, opts, out, extra)
  local votes, lo, hi = diagonal_votes(ref, qidx)
  -- A >=50 bp match at >=95% identity has at most 2 mismatches, leaving >= ~15
  -- exact 12-mers; 8 votes is a safe floor that discards chance seed pairs.
  -- Short references (< 36 bp) may legitimately have a single seed.
  local need = #ref < 3 * K and 1 or math.min(8, math.max(2, math.floor((#ref - K) / 4)))
  local diags = {}
  for d, c in pairs(votes) do
    if c >= need then diags[#diags + 1] = { d, c } end
  end
  if #diags == 0 then return end
  table.sort(diags, function(a, b) return a[2] > b[2] end)
  local used = {}
  for _, dc in ipairs(diags) do
    local d = dc[1]
    local near = false
    for _, u in ipairs(used) do
      if math.abs(u - d) <= opts.band then near = true break end
    end
    if not near then
      used[#used + 1] = d
      -- align only the window supported by seeds on nearby diagonals, plus
      -- enough margin for the alignment to reach the reference ends
      local r0, r1 = lo[d], hi[d] + K - 1
      for dd = d - opts.band, d + opts.band do
        if lo[dd] then
          r0 = math.min(r0, lo[dd])
          r1 = math.max(r1, hi[dd] + K - 1)
        end
      end
      local margin = 2 * opts.band + 2 * K
      local a = banded_sw(ref, q, d, opts.band, math.max(1, r0 - margin), math.min(#ref, r1 + margin))
      if a then
        local identity = a.matches / a.columns
        local coverage = (a.r1 - a.r0 + 1) / #ref
        local full = identity >= opts.min_identity and coverage >= opts.min_coverage
        local frag = not full and opts.fragments and identity >= opts.min_identity
          and coverage >= opts.min_fragment_coverage and (a.r1 - a.r0 + 1) >= opts.min_fragment_length
        if full or frag then
          a.identity, a.coverage, a.fragment = identity, coverage, not full
          for k, v in pairs(extra) do a[k] = v end
          out[#out + 1] = a
        end
      end
    end
  end
end

--- Detect library features in a sequence.
---@param seq string
---@param library { id: string, name: string, type: string, blurb?: string, seq: string }[]
---@param opts table|nil see M.defaults; plus circular = boolean
---@return table[] hits { ref, start, stop, strand, wraps, identity, coverage, fragment, score, ref_from, ref_to }
function M.search(seq, library, opts)
  local o = {}
  for k, v in pairs(M.defaults) do o[k] = v end
  for k, v in pairs(opts or {}) do o[k] = v end
  opts = o
  local L = #seq
  local maxref = 0
  for _, r in ipairs(library) do maxref = math.max(maxref, #r.seq) end
  local ext = math.min(L, maxref + opts.band)
  local fwd = seq
  if opts.circular then fwd = seq .. seq:sub(1, ext) end
  local rev = alphabet.revcomp(seq)
  if opts.circular then rev = rev .. rev:sub(1, ext) end
  local fidx, ridx = index_kmers(fwd, opts.max_positions), index_kmers(rev, opts.max_positions)

  local raw = {}
  local fwd_up, rev_up
  for _, r in ipairs(library) do
    if #r.seq < 20 then
      -- below 20 bp a single mismatch already falls under 95% identity, so an
      -- exact search is equivalent (and short parts cannot be seeded anyway)
      fwd_up = fwd_up or fwd:upper():gsub("U", "T")
      rev_up = rev_up or rev:upper():gsub("U", "T")
      local needle = r.seq:upper():gsub("U", "T")
      for strand, hay in pairs({ [1] = fwd_up, [-1] = rev_up }) do
        local from = 1
        while true do
          local a = hay:find(needle, from, true)
          if not a then break end
          raw[#raw + 1] = { ref = r, strand = strand, q0 = a, q1 = a + #needle - 1, r0 = 1, r1 = #needle,
            identity = 1, coverage = 1, fragment = false, score = #needle * MATCH }
          from = a + 1
        end
      end
    else
      search_ref(r.seq, fwd, fidx, opts, raw, { ref = r, strand = 1 })
      search_ref(r.seq, rev, ridx, opts, raw, { ref = r, strand = -1 })
    end
  end

  -- map to forward-axis coordinates on the original sequence
  local hits = {}
  for _, a in ipairs(raw) do
    local q0, q1 = a.q0, a.q1
    local start, stop
    if a.strand == 1 then
      start, stop = q0, q1
    else
      start, stop = L - q1 + 1, L - q0 + 1
      if opts.circular then
        start = (start - 1) % L + 1
        stop = (stop - 1) % L + 1
      end
    end
    local keep = true
    if opts.circular then
      if a.strand == 1 then
        if q0 > L then keep = false end -- duplicate of a hit found at the start
        stop = (stop - 1) % L + 1
      else
        if q0 > L then keep = false end
      end
    end
    if keep and (q1 - q0 + 1) <= L then
      hits[#hits + 1] = {
        ref = a.ref, start = start, stop = stop, strand = a.strand, wraps = start > stop,
        identity = a.identity, coverage = a.coverage, fragment = a.fragment, score = a.score,
        ref_from = a.r0, ref_to = a.r1, length = q1 - q0 + 1,
      }
    end
  end
  M._nested_ok = opts.nested_same_type
  local out = M.resolve(hits, L)
  M._nested_ok = nil
  return out
end

local function span(h, L)
  if h.wraps then return { { h.start, L }, { 1, h.stop } } end
  return { { h.start, h.stop } }
end

local function overlap(a, b, L)
  local n = 0
  for _, x in ipairs(span(a, L)) do
    for _, y in ipairs(span(b, L)) do
      local lo, hi = math.max(x[1], y[1]), math.min(x[2], y[2])
      if hi >= lo then n = n + hi - lo + 1 end
    end
  end
  return n
end

--- Drop redundant hits: two hits on (nearly) the same span keep the better
--- one (full-length beats fragment, then score, then identity). A hit is
--- redundant only when the overlap covers >= 80% of the shorter AND >= 50% of
--- the longer, so a short part nested inside a long one survives.
function M.resolve(hits, L)
  table.sort(hits, function(a, b)
    if a.fragment ~= b.fragment then return not a.fragment end
    if a.score ~= b.score then return a.score > b.score end
    if a.identity ~= b.identity then return a.identity > b.identity end
    return a.ref.id < b.ref.id
  end)
  local kept = {}
  for _, h in ipairs(hits) do
    local redundant = false
    for _, k in ipairs(kept) do
      local ov = overlap(h, k, L)
      local short, long = math.min(h.length, k.length), math.max(h.length, k.length)
      if ov >= 0.8 * short and ov >= 0.5 * long then redundant = true break end
    end
    if not redundant and not M._nested_ok then
      for _, k in ipairs(kept) do
        if k.ref.type == h.ref.type and not k.fragment and k.length > h.length
          and overlap(h, k, L) >= 0.9 * h.length then
          redundant = true
          break
        end
      end
    end
    if not redundant then kept[#kept + 1] = h end
  end
  table.sort(kept, function(a, b) return a.start < b.start end)
  return kept
end

return M
