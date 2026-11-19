-- How long the cloning engine takes, and how the search grows with the number
-- of parts.
--
--   nvim -l tests/tools/cloning_bench.lua [out.json]
--
-- Reported for each run: wall time, the number of end-compatibility tests the
-- search made (bio.assembly.compatible / .homology are counted by replacing
-- them for the duration of the run), and the solutions found.
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({
  root .. "/lua/?.lua", root .. "/lua/?/init.lua", root .. "/tests/?.lua", package.path,
}, ";")
vim.opt.runtimepath:prepend(root)
vim.notify = function() end
require("plasmid").setup({})

local A = require("plasmid.bio.assembly")
local pio = require("plasmid.io")
local Record = require("plasmid.model.record").Record

local out = { runs = {} }
local function mk(name, seq, topo)
  return Record.new({ name = name, seq = seq, topology = topo or "linear", alphabet = "dna" })
end
local seedstate = 1
local function rnd(n)
  local b, t = { "A", "C", "G", "T" }, {}
  for i = 1, n do
    seedstate = (1103515245 * seedstate + 12345) % 2147483648
    t[i] = b[seedstate % 4 + 1]
  end
  return table.concat(t)
end

--- Run fn with the join functions counted. Returns ms, calls, result.
local function timed(fn)
  local calls = 0
  local realc, realh = A.compatible, A.homology
  A.compatible = function(a, b) calls = calls + 1; return realc(a, b) end
  A.homology = function(a, b, o) calls = calls + 1; return realh(a, b, o) end
  local t0 = vim.uv.hrtime()
  local ok, res = pcall(fn)
  local ms = (vim.uv.hrtime() - t0) / 1e6
  A.compatible, A.homology = realc, realh
  if not ok then error(res) end
  return ms, calls, res
end

local function record(name, ms, calls, extra)
  local r = { name = name, ms = ms, joins = calls }
  for k, v in pairs(extra or {}) do r[k] = v end
  out.runs[#out.runs + 1] = r
  io.stdout:write(("%-54s %9.1f ms %10d joins  %s\n")
    :format(name, ms, calls, vim.inspect(extra or {}):gsub("%s+", " ")))
end

---------------------------------------------------------------------------
-- 1. A realistic 5-part Golden Gate
---------------------------------------------------------------------------
-- none palindromic, and no two the reverse complement of each other
local OV = { "AATG", "GCTT", "CAGA", "ACCT", "TACA", "GGAT", "CCAG", "AAAC",
  "ACTC", "GCAT", "AAGG", "CTGC" }
local function gg_set(n, body)
  local recs = {}
  for i = 1, n do
    local l, r = OV[i], OV[i % n + 1]
    recs[i] = mk("p" .. i, rnd(8) .. "GGTCTC" .. "A" .. l .. rnd(body) .. r .. "A" .. "GAGACC" .. rnd(8))
  end
  return recs
end

do
  local recs = gg_set(5, 700)
  local ms, calls, res = timed(function()
    local parts, problems = A.type_iis_parts(recs, "BsaI")
    local asms = A.assemble(parts, { circular = true })
    local rec = A.build(asms[1], { name = "gg5" })
    return { n = #asms, len = rec:len(), problems = #problems }
  end)
  record("golden gate, 5 parts of ~700 bp (cut+search+build)", ms, calls, res)

  local parts = A.type_iis_parts(gg_set(5, 700), "BsaI")
  local ms2, c2 = timed(function() return A.assemble(parts, { circular = true }) end)
  record("  the search alone", ms2, c2, {})
end

---------------------------------------------------------------------------
-- 2. Religating a large digest
---------------------------------------------------------------------------
local function load(rel) return assert(pio.read(root .. "/tests/fixtures/" .. rel)).records[1] end

for _, t in ipairs({
  { "ncbi/NC_001422.gb", { "HindIII" } }, { "ncbi/NC_001422.gb", { "HaeIII" } },
  { "ncbi/NC_001422.gb", { "AluI" } }, { "ncbi/J01749.gb", { "HinfI" } },
}) do
  local rec = load(t[1])
  local nfrag
  local msc, cc = timed(function() local f = A.cut(rec, t[2]); nfrag = #f; return f end)
  record(("digest %s / %s -> %d fragments"):format(vim.fn.fnamemodify(t[1], ":t"),
    table.concat(t[2], "+"), nfrag), msc, cc, { fragments = nfrag })
end

do  -- religation itself, growing the fragment count
  local rec = load("ncbi/NC_001422.gb")
  for _, spec in ipairs({ { "HaeIII" }, { "MboI" }, { "HinfI" }, { "AluI" },
    { "HindIII", "PstI", "EcoRV" } }) do
    local frags = A.cut(rec, spec)
    if #frags > 1 and not frags[1].uncut then
      local ms, calls, res = timed(function()
        local asms = A.assemble(frags, { circular = true, max_solutions = 50 })
        return { solutions = #asms }
      end)
      record(("religate phiX174 / %s (%d fragments)"):format(table.concat(spec, "+"), #frags),
        ms, calls, res, {})
    end
  end
end

do  -- the 230 kb record: cutting it, and what religation would cost
  local rec = load("ncbi/NC_001133.gb")
  for _, spec in ipairs({ { "NotI" }, { "SbfI" }, { "EcoRI" } }) do
    local frags
    local ms, calls = timed(function() frags = A.cut(rec, spec); return frags end)
    record(("digest NC_001133 230 kb / %s -> %d fragments"):format(spec[1], #frags), ms, calls,
      { fragments = #frags })
    if #frags > 1 and not frags[1].uncut then
      local ms2, c2, res = timed(function()
        return { solutions = #A.assemble(frags, { circular = false, max_solutions = 50 }) }
      end)
      record(("religate NC_001133 / %s (%d fragments)"):format(spec[1], #frags), ms2, c2, res)
    end
  end
end

do  -- pBR322 cut by a frequent sticky cutter, religated
  local rec = load("ncbi/J01749.gb")
  local frags = A.cut(rec, { "HinfI" })
  local ms, calls, res = timed(function()
    return { solutions = #A.assemble(frags, { circular = true, max_solutions = 50 }) }
  end)
  record(("religate pBR322 / HinfI (%d fragments)"):format(#frags), ms, calls, res)
end

---------------------------------------------------------------------------
-- 3. How the search grows
---------------------------------------------------------------------------
-- (a) a defined design: every overhang is used exactly twice, so only one
--     ordering survives and the search is linear in the number of parts.
for _, n in ipairs({ 2, 3, 4, 5, 6, 8, 10, 12 }) do
  local parts = A.type_iis_parts(gg_set(n, 120), "BsaI")
  local ms, calls, res = timed(function()
    return { solutions = #A.assemble(parts, { circular = true }) }
  end)
  record(("defined golden gate, %d parts"):format(n), ms, calls, res)
end

-- (b) the worst case, with the default cap in force: n blunt fragments, every
--     end compatible with every other, so the search stops as soon as it has
--     max_solutions orderings.
local function blunt_set(n)
  local frags = {}
  for i = 1, n do frags[i] = A.from_record(mk("b" .. i, rnd(200))) end
  return frags
end
for _, n in ipairs({ 2, 4, 6, 8, 10, 12, 16, 20 }) do
  local frags = blunt_set(n)
  local ms, calls, res = timed(function()
    return { solutions = #A.assemble(frags, { circular = true, max_solutions = 50 }) }
  end)
  record(("all-blunt, %d parts, cap 50"):format(n), ms, calls, res)
end

-- (c) the same with the cap lifted: what an exhaustive walk actually costs.
--     Stopped as soon as one point passes 20 s.
local stop_b = false
for _, n in ipairs({ 2, 3, 4, 5, 6, 7, 8 }) do
  if not stop_b then
    local frags = blunt_set(n)
    local ms, calls, res = timed(function()
      return { solutions = #A.assemble(frags, { circular = true, max_solutions = 1e9 }) }
    end)
    record(("all-blunt, %d parts, exhaustive"):format(n), ms, calls, res)
    if ms > 20000 then stop_b = true end
  end
end

-- (d) a one-enzyme digest religated: the same combinatorics with sticky ends,
--     where each compatibility test also compares the overhang text.
local function eco_set(n)
  local t = {}
  for i = 1, n do t[i] = "GAATTC" .. rnd(200) end
  return A.cut(mk("multi", table.concat(t), "circular"), { "EcoRI" })
end
for _, n in ipairs({ 2, 4, 6, 8, 10, 12, 16, 20 }) do
  local frags = eco_set(n)
  local ms, calls, res = timed(function()
    return { solutions = #A.assemble(frags, { circular = true, max_solutions = 50 }) }
  end)
  record(("one-enzyme digest religated, %d fragments, cap 50"):format(#frags), ms, calls, res)
end
local stop_e = false
for _, n in ipairs({ 2, 3, 4, 5, 6, 7, 8 }) do
  if not stop_e then
    local frags = eco_set(n)
    local ms, calls, res = timed(function()
      return { solutions = #A.assemble(frags, { circular = true, max_solutions = 1e9 }) }
    end)
    record(("one-enzyme digest religated, %d fragments, exhaustive"):format(#frags), ms, calls, res)
    if ms > 20000 then stop_e = true end
  end
end

local w = assert(io.open(arg[1] or "cloning_bench.json", "w"))
w:write(vim.json.encode(out))
w:close()
