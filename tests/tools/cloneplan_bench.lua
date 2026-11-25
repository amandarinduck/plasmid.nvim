-- Dev tool: what a bench plan costs, stage by stage.
--
--   nvim -l tests/tools/cloneplan_bench.lua [out.json]
--
-- Three assemblies: the Golden Gate fixtures, a vector + insert ligation whose
-- ends are all the same overhang (every by-product family is on the table),
-- and the same ligation with an insert too short to carry a six-base site, so
-- the diagnostic search has to decline. Prints one line per stage.
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({ root .. "/lua/?.lua", root .. "/lua/?/init.lua", package.path }, ";")
vim.opt.runtimepath:prepend(root)
vim.o.columns, vim.o.lines, vim.o.hidden = 150, 60, true
vim.notify = function() end

require("plasmid").setup({})
local C = require("plasmid.clone")
local CP = require("plasmid.cloneplan")
local D = require("plasmid.digest")
local pio = require("plasmid.io")
local Record = require("plasmid.model.record").Record

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function copy_fixture(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(root .. "/tests/fixtures/" .. rel, "rb"))
  local text = fh:read("*a"); fh:close()
  local out = assert(io.open(dst, "wb")); out:write(text); out:close()
  return dst
end
local function ms(fn)
  local t = vim.uv.hrtime()
  local a, b, c = fn()
  return (vim.uv.hrtime() - t) / 1e6, a, b, c
end
local function reset()
  C.staged = {}
  C.last_assembly = nil
  vim.cmd("silent! only")
  vim.cmd("silent! %bwipeout!")
end

local rows = {}
local function bench(name, ctx)
  local t_bp, bps, _, _ = ms(function() return CP.byproducts(ctx, {}) end)
  local t_d, diag, best, search = ms(function() return CP.diagnostic(ctx.product, bps, {}) end)
  local junctions = CP.build(ctx, {}).junctions
  local t_pcr = ms(function() return CP.colony_pcr(ctx.product, junctions, bps, {}) end)
  local t_seq = ms(function() return CP.sequencing(ctx.product, junctions, {}) end)
  local plan = CP.build(ctx, {})
  local t_render = ms(function() return CP.lines(plan) end)
  local t_all = ms(function() return CP.build(ctx, {}) end)
  local row = {
    name = name, product = ctx.product:len(), byproducts = #bps,
    answer = diag and table.concat(diag.names, " + ") or "none",
    told = (diag or best or {}).told or 0,
    sets = search.tried, enzymes = search.candidates,
    byproducts_ms = t_bp, diagnostic_ms = t_d, pcr_ms = t_pcr, sequencing_ms = t_seq,
    render_ms = t_render, total_ms = t_all,
  }
  rows[#rows + 1] = row
  io.stdout:write(("%-28s %6d bp  %2d by-products  %-12s %2d/%2d told  %4d sets/%d enzymes  "
    .. "bp %5.0f  diag %5.0f  pcr %5.0f  seq %5.0f  render %4.1f  total %5.0f ms\n")
    :format(name, row.product, row.byproducts, row.answer, row.told, row.byproducts,
      row.sets, row.enzymes, t_bp, t_d, t_pcr, t_seq, t_render, t_all))
end

-- 1. Golden Gate, the committed fixtures
reset()
local ipath, vpath = copy_fixture("synthetic/gg_insert.gb"), copy_fixture("synthetic/gg_vector.gb")
vim.cmd("silent edit! " .. vim.fn.fnameescape(vpath))
local vbuf = vim.api.nvim_get_current_buf()
C.run(vbuf, "goldengate", { ipath, vpath }, { enzyme = "BsmBI", first = true })
bench("golden gate (BsmBI)", C.last_assembly)

-- 2 and 3. pUC19/EcoRI with a long and a short EcoRI cassette
local function ligation(name, body)
  reset()
  local ins = Record.new({ name = name, seq = "GGGAATTCC" .. body .. "GGAATTCCC",
    alphabet = "dna", topology = "linear" })
  local ipath2 = tmp .. "/" .. name .. ".gb"
  assert(pio.write(ipath2, { ins }, "genbank"))
  local puc = copy_fixture("ncbi/L09137.gb")
  vim.cmd("silent edit! " .. vim.fn.fnameescape(puc))
  local vdb = D.run(vim.api.nvim_get_current_buf(), { "EcoRI" })
  C.add(vdb, "1")
  vim.cmd("silent edit! " .. vim.fn.fnameescape(ipath2))
  local idb = D.run(vim.api.nvim_get_current_buf(), { "EcoRI" })
  C.add(idb, "2")
  C.run(vdb, "ligation", {}, { first = true })
  return C.last_assembly
end

local pbr = assert(pio.read(root .. "/tests/fixtures/ncbi/J01749.gb")).records[1]
bench("ligation, 901 bp cassette", ligation("tetseg", pbr.seq:tostring():sub(1200, 2100):upper()))
math.randomseed(7)
local short = {}
for i = 1, 40 do short[i] = ("ACGT"):sub(math.random(4), math.random(4)) end
bench("ligation, 40 bp cassette", ligation("tiny", table.concat(short):gsub("[^ACGT]", "A")))

if arg[1] then
  local fh = assert(io.open(arg[1], "wb"))
  fh:write(vim.json.encode(rows))
  fh:close()
  io.stdout:write("wrote " .. arg[1] .. "\n")
end
