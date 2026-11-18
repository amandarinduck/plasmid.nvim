-- Primers: where they bind, their sequence, PCR products. Pure Lua.
--
-- A primer binds where its 3' end matches the template exactly for at least
-- `anchor` bases. The annealed region then extends towards the primer's 5'
-- end across a stretch of mismatches (mismatches fewer than `resume` matches
-- apart count as one stretch) only when at least `resume` consecutive
-- matching bases follow it (so a mutagenic codon inside a primer is annealed, while a 5'
-- tail that happens to match a base or two is not), with at most
-- `max_mismatches` in all. Whatever lies 5' of the annealed region is the
-- tail (restriction sites, overlaps, ...). Calibrated on SnapGene's own
-- annealing calls (19 of 20 primer sites in the test files agree).
local search = require("plasmid.bio.search")

local M = {}

M.defaults = { anchor = 10, max_mismatches = 3, resume = 5, min_anneal = 10 }

local COMP = {
  A = "T", C = "G", G = "C", T = "A", U = "A", R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
  B = "V", V = "B", D = "H", H = "D", N = "N",
}
local MEANS = { A = "A", C = "C", G = "G", T = "T", U = "T", R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT",
  M = "AC", B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT" }

function M.clean(s)
  return (s or ""):gsub("[%s%-]", ""):upper():gsub("U", "T")
end

function M.revcomp(s)
  local t = {}
  local u = s:upper()
  for i = #u, 1, -1 do t[#t + 1] = COMP[u:sub(i, i)] or "N" end
  return table.concat(t)
end

function M.valid(s)
  s = M.clean(s)
  if s == "" then return false, "empty primer" end
  local bad = s:match("[^ACGTRYSWKMBDHVN]")
  if bad then return false, ("not a nucleotide: %s"):format(bad) end
  return true
end

-- does template base x satisfy primer base p (subset semantics, as search)
local function ok(x, p)
  local mx, mp = MEANS[x], MEANS[p]
  if not mx or not mp then return false end
  for i = 1, #mx do if not mp:find(mx:sub(i, i), 1, true) then return false end end
  return true
end

--- Binding sites of a primer (5'->3') on a template string.
--- Returns a list sorted by position of
--- { strand = 1|-1, start, stop (annealed region on the top strand; stop <
---   start wraps the origin), anneal = n bases, tail = 5' unannealed part,
---   mismatches = { primer index, ... }, exact3 = contiguous exact 3' bases }.
--- strand 1: the primer anneals to the bottom strand and copies the top one
--- (a forward primer); strand -1: a reverse primer.
function M.bind(template, circular, primer, opts)
  opts = vim.tbl_extend("force", M.defaults, opts or {})
  local P = M.clean(primer)
  local n = #P
  local T = template:upper():gsub("U", "T")
  local L = #T
  local a = math.min(opts.anchor, n)
  if n == 0 or L == 0 or a > L then return {} end
  local anchor = P:sub(n - a + 1)
  local hits = search.find_nuc(T, anchor, { circular = circular })
  local out = {}
  local function at(p)
    if circular then p = (p - 1) % L + 1 elseif p < 1 or p > L then return nil end
    return T:sub(p, p)
  end
  for _, h in ipairs(hits) do
    -- k-th primer base from the 3' end sits at template position pos(k);
    -- forward: pairs with the top strand base equal to it; reverse: with its complement
    local pos, want
    if h.strand == 1 then
      local e = h.start + a - 1
      pos = function(k) return e - (k - 1) end
      want = function(x, pb) return ok(x, pb) end
    else
      local s = h.start
      pos = function(k) return s + (k - 1) end
      want = function(x, pb) return ok(x, COMP[pb]) end
    end
    -- the template base under the k-th primer base from the 3' end (nil past an end)
    local function pairs_at(k)
      if circular and k > L then return nil end
      local x = at(pos(k))
      if not x then return nil end
      return want(x, P:sub(n - k + 1, n - k + 1))
    end
    local k = 1
    while k <= n and pairs_at(k) do k = k + 1 end
    local exact3 = k - 1
    local last_match, mism = exact3, {}
    -- extend through a stretch of mismatches (mismatches closer than `resume`
    -- matches form one stretch, e.g. a changed codon) when at least `resume`
    -- consecutive matches follow it
    while k <= n and pairs_at(k) == false do
      local j, block, run = k, {}, 0
      while j <= n do
        local r = pairs_at(j)
        if r == nil then break end
        if r then
          run = run + 1
          if run >= opts.resume then break end
        else
          block[#block + 1] = n - j + 1
          run = 0
          if #mism + #block > opts.max_mismatches then break end
        end
        j = j + 1
      end
      if run < opts.resume or #mism + #block > opts.max_mismatches then break end
      for _, i in ipairs(block) do mism[#mism + 1] = i end
      k = j + 1
      while k <= n and pairs_at(k) do k = k + 1 end
      last_match = k - 1
    end
    local keep = mism
    local anneal = last_match
    if anneal >= math.min(opts.min_anneal, n) then
      local p3, p5 = pos(1), pos(anneal)
      local s, e
      if h.strand == 1 then s, e = p5, p3 else s, e = p3, p5 end
      if circular then s, e = (s - 1) % L + 1, (e - 1) % L + 1 end
      table.sort(keep)
      out[#out + 1] = { strand = h.strand, start = s, stop = e, anneal = anneal, tail = P:sub(1, n - anneal),
        mismatches = keep, exact3 = exact3, primer = P }
    end
  end
  table.sort(out, function(x, y)
    if x.start ~= y.start then return x.start < y.start end
    return x.strand > y.strand
  end)
  return out
end

--- The full primer sequence of a primer_bind feature: SnapGene's stored
--- sequence, else a "primer sequence: ..." note, else the template bases.
function M.sequence(rec, f)
  if f.primer_seq and f.primer_seq ~= "" then return M.clean(f.primer_seq), "stored" end
  for _, v in ipairs(f:get_all("note")) do
    local s = tostring(v):match("[Pp]rimer sequence:%s*5?'?%-?([ACGTUacgtuRYSWKMBDHVNryswkmbdhvn ]+)")
    if s then return M.clean(s), "note" end
  end
  return M.clean(rec:extract(f)), "template"
end

function M.is_primer(f)
  return f.type == "primer_bind" or f.type == "primer"
end

--- The note that stores a primer's sequence in GenBank files.
function M.note(seq) return "primer sequence: " .. M.clean(seq) end

return M
