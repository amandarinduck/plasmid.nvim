-- Dev tool: what the circular map costs and what its rings claim, on records
-- too large to commit. Not part of the suite.
--
--   tests/tools/fetch_big.sh
--   nvim -l tests/tools/ring_bench.lua big/NC_000913.gbk [WIDTH HEIGHT]
--
-- Prints the build times at both scales, the bins and the window the header
-- states, and -- for a circular record -- where the cumulative GC skew turns,
-- which is the claim worth checking against a genome whose replication origin
-- is known independently.
package.path = "lua/?.lua;lua/?/init.lua;" .. package.path
local config = require("plasmid.config")
config.setup({})
require("plasmid.view.highlights").setup()
local circular = require("plasmid.view.circular")
local ring = require("plasmid.view.ring")

local path = assert(arg[1], "usage: ring_bench.lua FILE [W H]")
local W, H = tonumber(arg[2]) or 120, tonumber(arg[3]) or 40

local function say(...) io.stdout:write(string.format(...), "\n") end
local function ms(t0) return (vim.uv.hrtime() - t0) / 1e6 end
local function heap() collectgarbage(); return collectgarbage("count") / 1024 end

local base = heap()
local t0 = vim.uv.hrtime()
local rec = assert(require("plasmid.io").read(path)).records[1]
say("%-28s %8.0f ms  %6.0f MB", "parse", ms(t0), heap() - base)
say("%-28s %8d %s, %d features, %s", "record", rec:len(), rec:unit(), #rec.features, rec.topology)

local after_parse = heap()
t0 = vim.uv.hrtime()
local lay = circular.build(rec, config.options, W, H)
local cold = ms(t0)
local cold_heap = heap() - after_parse
t0 = vim.uv.hrtime()
lay = circular.build(rec, config.options, W, H)
say("%-28s %8.0f ms  %6.0f MB", ("build %dx%d (first)"):format(W, H), cold, cold_heap)
say("%-28s %8.0f ms", "build again (cached)", ms(t0))
say("%-28s %8s", "scale", lay.atlas and "atlas" or "arcs")

if lay.atlas then
  local c = lay.composition
  say("%-28s %8d %s in %d bins", "bin", lay.bin, rec:unit(), c.n)
  say("%-28s %8.4f  max deviation %.4f", "mean gc", c.mean, c.maxdev)
  say("%-28s %8d = %.3f%% of the record", "cumulative skew swing", c.swing, 100 * c.swing / rec:len())
  if c.turning then
    say("%-28s %8d", "skew minimum (origin)", c.origin_bp)
    say("%-28s %8d", "skew maximum (terminus)", c.terminus_bp)
  else
    say("%-28s %8s", "skew turning points", "not claimed")
  end
  say("%-28s %8d of %d features", "in the rings", #lay.items, #lay.items + lay.history + lay.off_seq)
  say("%-28s %8d", "landmarks named", lay.named)
  for _, m in ipairs(lay.landmarks) do say("    %10d  %s", m.pos, m.label) end
end
say("%-28s %8d", "tick interval", ring.tick_step(rec:len(), lay.atlas and "atlas" or "arcs"))
for i = 1, 2 do say("header %d: %s", i, lay.lines[i]) end
