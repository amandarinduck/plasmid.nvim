-- Dev tool: where the windowed open's planning pass spends its time.
--
--   PLASMID_PROFILE="FILE [width=120] [rec=1] [reps=3] [hook=200]" \
--     nvim --headless --clean -u NONE -c 'luafile tests/tools/plan_profile.lua'
--
-- Opening a large record windowed costs two things: the setup every layout
-- does (visible features, global track packing, the two-level block index) and
-- a pass that plans every block purely to learn how many rows it has. This
-- separates them, counts what the pass allocates, and -- with a count hook,
-- which switches the interpreter off the JIT and so inflates the absolute
-- time -- attributes samples to the lines of linear.lua they landed on, so the
-- phases inside plan() can be told apart.
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({ root .. "/lua/?.lua", root .. "/lua/?/init.lua", package.path }, ";")
vim.opt.runtimepath:prepend(root)
vim.notify = function() end

local argv = arg
if not argv or #argv == 0 then
  argv = vim.split(os.getenv("PLASMID_PROFILE") or "", "%s+", { trimempty = true })
end
local path = argv[1]
local A = {}
for i = 2, #argv do
  local k, v = argv[i]:match("^([%w_]+)=(.*)$")
  if k then A[k] = v end
end
local function num(k, d) return tonumber(A[k]) or d end

vim.o.columns, vim.o.lines = num("width", 120), 40
require("plasmid").setup({ linear = { window = { enabled = true } } })
local io_ = require("plasmid.io")
local config = require("plasmid.config")
local linear = require("plasmid.view.linear")

local function ms(fn)
  local t = vim.uv.hrtime()
  local v = fn()
  return (vim.uv.hrtime() - t) / 1e6, v
end

-- the parse is the other half of an open, and the only honest way to say
-- whether the layout is still worth working on
collectgarbage("collect")
local parse_ms, res = ms(function() return io_.read(path) end)
assert(res, "could not read " .. tostring(path))
local rec = res.records[num("rec", 1)]
local o = vim.deepcopy(config.options)
o.linear.window = { enabled = true, margin = 2, cache_blocks = 512 }

collectgarbage("collect")
local build_ms, lay = ms(function() return linear.build(rec, o, vim.o.columns, {}) end)
local nb = lay.nblocks

--- One run of a pass over every block.
local function once(fn)
  collectgarbage("collect")
  return (ms(function()
    local tot = 0
    for b = 0, nb - 1 do tot = tot + fn(b) end
    return tot
  end))
end

--- Medians of several runs of each pass, taken in turn rather than one pass
--- after the other: whatever the machine is doing while this runs, and
--- whatever the tracing compiler has made of the loops by then, both passes
--- meet it.
local function passes(fns)
  local t = {}
  for i in pairs(fns) do t[i] = {} end
  for _ = 1, num("reps", 3) do
    for i, fn in pairs(fns) do table.insert(t[i], once(fn)) end
  end
  local out = {}
  for i, xs in pairs(t) do
    table.sort(xs)
    out[i] = xs[math.ceil(#xs / 2)]
  end
  return out
end
local function planned(b) return lay._plan(b).n end
local med = passes({ plan = planned, count = lay._count })
local plan_ms, count_ms = med.plan, med.count
-- the first build in a process pays for compiling the module; a second one is
-- warm, like the pass measurements above it, so setup and planning can be
-- subtracted from each other without comparing a cold number with a warm one
collectgarbage("collect")
local warm_build_ms = ms(function() return linear.build(rec, o, vim.o.columns, {}) end)

-- what the pass allocates, with the collector off so nothing is reclaimed
local function allocated(fn)
  collectgarbage("collect")
  collectgarbage("stop")
  local k0 = collectgarbage("count")
  for b = 0, nb - 1 do local _ = fn(b) end
  local kb = collectgarbage("count") - k0
  collectgarbage("restart")
  collectgarbage("collect")
  return kb
end
--- The same pass with the collector stopped: the difference against the
--- ordinary pass is what collecting the garbage it makes costs.
local function nogc(fn)
  collectgarbage("collect")
  collectgarbage("stop")
  local d = ms(function()
    local tot = 0
    for b = 0, nb - 1 do tot = tot + fn(b) end
    return tot
  end)
  collectgarbage("restart")
  collectgarbage("collect")
  return d
end
local plan_nogc_ms = nogc(planned)
local alloc_kb = allocated(planned)
local count_alloc_kb = lay._count and allocated(lay._count) or nil

-- the floor: the same loop over the same blocks doing only its arithmetic
local floor_ms = once(function(b)
  local bs = b * lay.W + 1
  return math.min(bs + lay.W - 1, lay.L)
end)

-- Sampled profile of the pass, attributed to the enclosing function rather
-- than to single lines: the line a sample lands on moves whenever the file is
-- edited, the function it is in does not.
local hits, lhits, total = {}, {}, 0
local function hook()
  local i = debug.getinfo(2, "Sl")
  if not i then return end
  local key = vim.fn.fnamemodify(i.short_src, ":t") .. ":" .. (i.linedefined or 0)
  hits[key] = (hits[key] or 0) + 1
  lhits[key .. "/" .. (i.currentline or 0)] = (lhits[key .. "/" .. (i.currentline or 0)] or 0) + 1
  total = total + 1
end
collectgarbage("collect")
local hook_ms = ms(function()
  debug.sethook(hook, "", num("hook", 200))
  for b = 0, nb - 1 do local _ = lay._plan(b).n end
  debug.sethook()
end)

-- What a session that visits a record without scrolling ends up holding: a
-- position query plans the block it lands in and keeps the plan, and only
-- `ensure` sweeps those.
collectgarbage("collect")
local held0 = collectgarbage("count") / 1024
for b = 0, nb - 1 do lay:plan_of(b) end
collectgarbage("collect")
local held_mb = collectgarbage("count") / 1024 - held0
local held_n = 0
for _ in pairs(lay._plans) do held_n = held_n + 1 end

local rows = {}
for k, v in pairs(hits) do rows[#rows + 1] = { k, v } end
table.sort(rows, function(x, y) return x[2] > y[2] end)

--- Source of the modules the samples can land in, so a bucket can be named
--- after the function it is.
local lines = {}
for _, rel in ipairs({ "view/linear.lua", "view/features.lua" }) do
  local fh = io.open(root .. "/lua/plasmid/" .. rel, "r")
  if fh then
    local name, t = vim.fn.fnamemodify(rel, ":t"), {}
    for l in fh:lines() do t[#t + 1] = l end
    fh:close()
    lines[name] = t
  end
end

local w = io.stdout
w:write(("%s  %s  %d bp, %d features, W=%d, %d blocks, %d rows\n")
  :format(vim.fn.fnamemodify(path, ":t"), rec.name, rec:len(), #rec.features, lay.W, nb, #lay.lines))
w:write(("  parsing the file      %8.1f ms\n"):format(parse_ms))
w:write(("  windowed build, cold  %8.1f ms\n"):format(build_ms))
w:write(("  windowed build, warm  %8.1f ms   setup %.1f ms + %.1f ms of %s\n")
  :format(warm_build_ms, warm_build_ms - (count_ms or plan_ms), count_ms or plan_ms,
    count_ms and "counting" or "planning"))
w:write(("  planning pass         %8.1f ms   %.3f us/block\n"):format(plan_ms, plan_ms * 1000 / nb))
w:write(("  it allocates          %8.1f MB   %.0f bytes/block\n"):format(alloc_kb / 1024, alloc_kb * 1024 / nb))
w:write(("  the pass, collector off %6.1f ms   so collecting costs %.1f ms\n"):format(plan_nogc_ms, plan_ms - plan_nogc_ms))
if count_ms then
  w:write(("  counting pass         %8.1f ms   %.3f us/block\n"):format(count_ms, count_ms * 1000 / nb))
  w:write(("  it allocates          %8.1f MB   %.0f bytes/block\n"):format(count_alloc_kb / 1024, count_alloc_kb * 1024 / nb))
end
w:write(("  arithmetic-only floor %8.1f ms   %.3f us/block\n"):format(floor_ms, floor_ms * 1000 / nb))
w:write(("  a plan kept per block %8.1f MB   %d plans held\n"):format(held_mb, held_n))
w:write(("  under the count hook  %8.1f ms   %d samples\n"):format(hook_ms, total))
w:write("  where the samples landed, by function:\n")
for i = 1, math.min(12, #rows) do
  local file, ln = rows[i][1]:match("^(.*):(%d+)$")
  local text = (lines[file] or {})[tonumber(ln)] or ""
  -- the hottest single line inside the bucket, which says what it is doing
  local best, bn = 0, 0
  for k, v in pairs(lhits) do
    local kk, cl = k:match("^(.*)/(%d+)$")
    if kk == rows[i][1] and v > best then best, bn = v, tonumber(cl) end
  end
  w:write(("   %5.1f%%  %-22s %s\n"):format(100 * rows[i][2] / total,
    file .. ":" .. ln, (text:gsub("^%s+", ""):sub(1, 64))))
  w:write(("           %5.1f%% on  %-18s %s\n"):format(100 * best / total, file .. ":" .. bn,
    (((lines[file] or {})[bn] or ""):gsub("^%s+", ""):sub(1, 64))))
end
vim.cmd("qa!")
