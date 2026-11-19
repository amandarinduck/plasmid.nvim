-- Drives bio.assembly over the committed cloning oracle (tests/oracle/
-- expected/cloning.json, written by tests/oracle/dump_cloning.py).
--
-- One case in, one result out: the distinct products the engine finds, keyed
-- the same way the oracle keys them, so that the spec and the development
-- driver (tests/tools/cloning_compare.lua) agree on what "the same product"
-- means.
--
-- A product is identified by its length and the sha256 of a canonical form:
-- for a circular product the lexicographically least rotation of the sequence
-- or of its reverse complement, whichever is smaller; for a linear one
-- min(seq, revcomp(seq)). Rotation and strand therefore do not matter.
local A = require("plasmid.bio.assembly")
local pio = require("plasmid.io")
local Record = require("plasmid.model.record").Record
local alphabet = require("plasmid.model.alphabet")

local M = {}

--- Booth's algorithm: the lexicographically least rotation of `s`, in O(#s).
function M.least_rotation(s)
  local n = #s
  if n == 0 then return s end
  local t = s .. s
  local byte = string.byte
  local f, k = {}, 0
  for j = 1, 2 * n - 1 do
    local sj = byte(t, j + 1)
    local i = f[j - k - 1]
    if i == nil then i = -1 end
    while i ~= -1 and sj ~= byte(t, k + i + 2) do
      if sj < byte(t, k + i + 2) then k = j - i - 1 end
      i = f[i]
      if i == nil then i = -1 end
    end
    if sj ~= byte(t, k + i + 2) then
      if sj < byte(t, k + 1) then k = j end
      f[j - k] = -1
    else
      f[j - k] = i + 1
    end
  end
  return t:sub(k + 1, k + n)
end

function M.canonical(seq, circular)
  local s = seq:upper()
  local r = alphabet.revcomp(s)
  local a, b
  if circular then a, b = M.least_rotation(s), M.least_rotation(r)
  else a, b = s, r end
  if a <= b then return a end
  return b
end

function M.key(seq, circular)
  return { len = #seq, canon = vim.fn.sha256(M.canonical(seq, circular)) }
end

--- Distinct product keys of a list of assemblies, in the oracle's order.
local function products(asms, circular)
  local out, seen = {}, {}
  for _, asm in ipairs(asms) do
    local rec = A.build(asm, { annotate_parts = false, annotate_junctions = false })
    local k = M.key(rec.seq:tostring(), circular)
    local id = k.len .. ":" .. k.canon
    if not seen[id] then seen[id] = true; out[#out + 1] = k end
  end
  table.sort(out, function(x, y)
    if x.len ~= y.len then return x.len < y.len end
    return x.canon < y.canon
  end)
  return out
end

---------------------------------------------------------------------------
-- Inputs
---------------------------------------------------------------------------

local cache = {}
local function fixture(root, rel)
  if not cache[rel] then
    cache[rel] = assert(pio.read(root .. "/tests/fixtures/" .. rel)).records[1]
  end
  return cache[rel]
end

local function linear(name, seq)
  return Record.new({ name = name, seq = seq, topology = "linear", alphabet = "dna" })
end

--- The fragments a case feeds to the assembler, plus any problems the method
--- reported while preparing them.
function M.fragments(case, root)
  if case.method == "ligation" then
    local frags = {}
    for si, src in ipairs(case.sources) do
      local rec = fixture(root, src.file)
      if src.slice then
        local s = rec.seq:tostring():sub(src.slice[1], src.slice[2])
        rec = linear(rec.name, s)
      end
      local got = A.cut(rec, src.enzymes, ("s%d"):format(si))
      if src.take == "all" then
        for _, f in ipairs(got) do frags[#frags + 1] = f end
      else
        for _, i in ipairs(src.take) do frags[#frags + 1] = assert(got[i], "no such fragment") end
      end
    end
    return frags, {}
  elseif case.method == "gibson" then
    local recs = {}
    for _, p in ipairs(case.parts) do recs[#recs + 1] = linear(p.name, p.seq) end
    return A.gibson_parts(recs)
  elseif case.method == "goldengate" then
    local recs = {}
    for _, p in ipairs(case.parts) do recs[#recs + 1] = linear(p.name, p.seq) end
    local frags, problems = A.type_iis_parts(recs, case.enzyme)
    return frags or {}, problems or {}
  end
  error("unknown method " .. tostring(case.method))
end

--- Run one case. Returns
---   { circular = {key...}, linear = {key...} | nil, frag_lens = {...},
---     problems = {...}, nfrags = n, capped = bool }
--- `linear` is nil when the engine hit its solution cap, which makes the
--- linear list incomplete and so not comparable.
function M.run(case, root, opts)
  opts = opts or {}
  local max = opts.max_solutions or 400
  local frags, problems = M.fragments(case, root)
  local lens = {}
  for i, f in ipairs(frags) do lens[i] = f.rec:len() end

  local join = case.method == "gibson" and "homology" or "ends"
  local base = { join = join, min_homology = case.min_homology, max_solutions = max }

  local res = { frag_lens = lens, nfrags = #frags, problems = problems }
  if #frags == 0 then
    res.circular, res.linear = {}, {}
    return res
  end
  for _, f in ipairs(frags) do
    if f.circular then        -- an intact circle cannot take part in a reaction
      res.circular, res.linear = {}, {}
      res.uncut = true
      return res
    end
  end

  local c = vim.tbl_extend("force", base, { circular = true })
  local casms = A.assemble(frags, c)
  res.circular = products(casms, true)
  res.ncircular_assemblies = #casms

  local l = vim.tbl_extend("force", base, { circular = false })
  local lasms = A.assemble(frags, l)
  res.nlinear_assemblies = #lasms
  if #lasms >= max then res.linear = nil else res.linear = products(lasms, false) end
  return res
end

---------------------------------------------------------------------------
-- Comparison
---------------------------------------------------------------------------

local function keyset(list)
  local t = {}
  for _, k in ipairs(list or {}) do t[k.len .. ":" .. k.canon] = true end
  return t
end

local function diff(want, got)
  local a, b = keyset(want), keyset(got)
  local missing, extra = {}, {}
  for k in pairs(a) do if not b[k] then missing[#missing + 1] = k end end
  for k in pairs(b) do if not a[k] then extra[#extra + 1] = k end end
  table.sort(missing); table.sort(extra)
  return missing, extra
end

--- Compare one case against its oracle entry. Returns a list of complaints
--- (empty when they agree) and the result table.
function M.compare(case, root, opts)
  opts = opts or {}
  local res = M.run(case, root, opts)
  local bad = {}

  if case.frag_lens then
    local want = case.frag_lens
    if #want ~= #res.frag_lens then
      bad[#bad + 1] = ("fragment count: oracle %d, engine %d (%s vs %s)")
        :format(#want, #res.frag_lens, table.concat(want, ","), table.concat(res.frag_lens, ","))
    else
      for i = 1, #want do
        if want[i] ~= res.frag_lens[i] then
          bad[#bad + 1] = ("fragment lengths: oracle %s, engine %s")
            :format(table.concat(want, ","), table.concat(res.frag_lens, ","))
          break
        end
      end
    end
  end

  -- pydna reports one product record per assembly it found; the engine
  -- reports one assembly per ordering. They should agree on how many there
  -- are, not only on which molecules come out.
  if case.n_circ_raw and case.n_circ_raw ~= res.ncircular_assemblies then
    bad[#bad + 1] = ("circular assemblies: oracle %d, engine %d")
      :format(case.n_circ_raw, res.ncircular_assemblies)
  end

  local missing, extra = diff(case.expect.circular, res.circular)
  if #missing > 0 or #extra > 0 then
    bad[#bad + 1] = ("circular products: oracle %d, engine %d; missing %d, extra %d")
      :format(#case.expect.circular, #res.circular, #missing, #extra)
  end

  -- Designed failures must be explained, not merely absent: the whole point of
  -- naming the overhangs is that the user is told which one has no partner.
  local SILENT = { lowered_floor = true, flipped_part = true, self_circularisation = true,
    homology_too_short = true, ambiguous_overlap = true }
  if case.method == "goldengate" then
    local want = case.designed ~= nil and not SILENT[case.designed]
    local got = #res.problems > 0
    if want ~= got then
      bad[#bad + 1] = ("problems: expected %s, got {%s}")
        :format(want and "at least one" or "none", table.concat(res.problems, "; "))
    end
  end

  local wl = case.expect.linear
  if wl and res.linear and #case.expect.circular == 0 and #res.circular == 0
    and #wl <= 200 and #res.linear <= 200 then
    local m2, e2 = diff(wl, res.linear)
    if #m2 > 0 or #e2 > 0 then
      bad[#bad + 1] = ("linear products: oracle %d, engine %d; missing %d, extra %d")
        :format(#wl, #res.linear, #m2, #e2)
    end
  end
  res.n_circ_raw = case.n_circ_raw
  return bad, res
end

return M
