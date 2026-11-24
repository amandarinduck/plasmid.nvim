-- What the linear view is allowed to cost on a large record.
--
-- The point of the windowed layout is that opening a record does not lay it
-- out: it plans every block to learn the buffer's length and builds text only
-- for the rows a window is about to show. A change that quietly goes back to
-- laying everything out up front would still be correct, and every other spec
-- would still pass; these limits are what would fail.
--
-- The record is generated here -- 400 kb and 4,000 features, with joins,
-- complements and enough overlap to fill the feature tracks -- so the suite
-- needs no network and no download. It is about a tenth of E. coli K-12, which
-- is the size this milestone is really about; the real genome is measured by
-- tests/tools/big_budget.lua and written up in docs/large_records_report.md.
--
-- Measured by this spec on an idle 12-core x86-64 machine, Neovim 0.12.2 --
-- the median of three runs in each of the sequence layouts the suite runs
-- over, group 10 and group 0 (PLASMID_TEST_GROUP=0). Heap figures are what is
-- still held after a collection, not peak allocation.
--
--                                 group 10        group 0         limit
--   windowed build               47 ms/2.6 MB    48 ms/2.6 MB    400 ms/15 MB
--   eager build, same record    543 ms/43.9 MB  439 ms/43.7 MB   ratio only
--   body rows laid out at build    0               0               0
--   one screen                     2.2 ms          1.8 ms        150 ms
--   200 pages of scrolling       111 ms           83 ms         3000 ms
--   rows still laid out after   2985 of 30668   3388 of 31998    a quarter
--   400 jumps (pos_of only)        3.1 ms          2.7 ms        300 ms
--   200 cold rows through info()  21 ms           19 ms         1500 ms
--   a query in every block      4445 plans/3.4  4000 plans/3.8  one per block
--   ...then one screen           19 plans/0.1    16 plans/0.0   an eighth, half
--   open a 400 kb file, :edit    124 ms/10.8 MB  121 ms/10.4 MB 3000 ms/60 MB
--   200 pages through buffer.fill 62 ms           64 ms         3000 ms
--   40 jumps through goto_bp      87 ms           93 ms         2000 ms
--   re-render after two resizes  113 ms          115 ms         4000 ms
--   process growth over the case   7 MB            7 MB          250 MB
--
-- The wall-clock limits are several times the measured values because a loaded
-- CI runner is several times slower; they are still well below what the eager
-- layout costs here. The ratios asserted below -- the eager layout takes at
-- least 3x the time and 5x the heap -- do not depend on machine speed at all,
-- and are what would catch a return to laying out every row up front. The
-- lowest ratios seen over those six runs were 9.5x on time and 17.2x on heap.
-- The process-growth limit is loose on purpose: resident memory depends on the
-- allocator and on whatever ran before this spec.
--
-- Those are the numbers when this spec is run on its own. Inside the whole
-- suite the buffer cases are two to four times slower -- 351 ms for the 40
-- jumps in the run this comment was last checked against, against 87 ms alone
-- -- because 40 other spec files have been through the editor first. The
-- limits have room for that; a limit set from a solo run would not.
local H = require("helpers")

require("plasmid").setup({})
local config = require("plasmid.config")
local record = require("plasmid.model.record")
local linear = require("plasmid.view.linear")
local buffer = require("plasmid.buffer")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local L, NFEAT, WIDTH = 400000, 4000, 120

--- Base pairs to jump to, from a generator of this spec's own: seeding
--- math.random here would hand a different sequence to every spec that runs
--- after this one, and some of them draw from it without seeding first.
local function jumps(n, seed)
  local x, out = seed, {}
  for i = 1, n do
    x = (1103515245 * x + 12345) % 2147483648
    out[i] = x % L + 1
  end
  return out
end

--- A record the size of a small chromosome, built rather than read: a repeating
--- sequence, and features that overlap in threes so the track packer, the
--- overflow row and the translation rows all have work to do.
local function synthetic()
  local unit = "ATGAGCTTAAGGCATCGATTCCGTAGCTAAGGCTTAACCGGATATCGCGATTTGCAACGTAC"
  local parts = {}
  for i = 1, math.ceil(L / #unit) do parts[i] = unit end
  local feats = {}
  local step = math.floor(L / NFEAT)
  for i = 1, NFEAT do
    local s = (i - 1) * step + 1
    local e = math.min(L, s + step * 2)   -- each feature overlaps the next two
    local loc
    if i % 7 == 0 then loc = ("join(%d..%d,%d..%d)"):format(s, s + 60, s + 200, e)
    elseif i % 5 == 0 then loc = ("complement(%d..%d)"):format(s, e)
    else loc = ("%d..%d"):format(s, e) end
    feats[i] = record.Feature.new({
      type = (i % 3 == 0) and "CDS" or (i % 3 == 1 and "gene" or "misc_feature"),
      location_text = loc,
      qualifiers = { { key = "label", value = ("f%04d"):format(i), quoted = true } },
    })
  end
  return record.Record.new({ name = "SYNTH", seq = table.concat(parts):sub(1, L),
    topology = "circular", features = feats, format = "genbank" })
end

local function opts(windowed)
  local o = vim.deepcopy(config.options)
  -- the suite runs with group 10 and with group 0; both must stay inside these
  -- limits, so the configured one is left alone
  o.linear.window = { enabled = windowed, margin = 2, cache_blocks = 512 }
  return o
end

--- Wall clock and Lua heap growth of one step.
local function cost(fn)
  collectgarbage("collect")
  local h0 = collectgarbage("count")
  local t0 = vim.uv.hrtime()
  local v = fn()
  local ms = (vim.uv.hrtime() - t0) / 1e6
  collectgarbage("collect")
  return v, ms, (collectgarbage("count") - h0) / 1024
end

local function body_rows_laid_out(lay)
  local n = 0
  for i = lay.header_rows + 1, #lay.lines do
    if lay.lines[i] ~= "" then n = n + 1 end
  end
  return n
end

case("opening_a_large_record_lays_out_none_of_it", function()
  local rec = synthetic()
  local lay, ms, mb = cost(function() return linear.build(rec, opts(true), WIDTH, {}) end)
  H.ok(lay.windowed, "a 400 kb record should be windowed")
  H.eq(0, body_rows_laid_out(lay), "rows were laid out before anything asked for them")
  H.ok(ms < 400, ("windowed build took %.0f ms"):format(ms))
  H.ok(mb < 15, ("windowed build kept %.1f MB of heap"):format(mb))
  H.log(("    [budget] %d kb, %d features: windowed layout %.0f ms, %.1f MB, %d rows in %d blocks, none built")
    :format(L / 1000, NFEAT, ms, mb, #lay.lines, lay.nblocks))
end)

-- The ratios are the part that does not care how fast the machine is: whatever
-- a runner costs, laying out every row costs several times laying out none.
case("the_eager_layout_of_the_same_record_costs_much_more", function()
  local rec = synthetic()
  local wlay, wms, wmb = cost(function() return linear.build(rec, opts(true), WIDTH, {}) end)
  local wrows = #wlay.lines
  wlay = nil   -- do not hold both layouts at once
  local elay, ems, emb = cost(function() return linear.build(rec, opts(false), WIDTH, {}) end)
  H.eq(wrows, #elay.lines, "the two layouts disagree about how many rows the record has")
  H.ok(not elay.windowed, "enabled = false should force the eager path")
  H.ok(ems / wms > 3, ("eager %.0f ms vs windowed %.0f ms: only %.1fx"):format(ems, wms, ems / wms))
  H.ok(emb / wmb > 5, ("eager %.1f MB vs windowed %.1f MB: only %.1fx"):format(emb, wmb, emb / wmb))
  H.log(("    [budget] eager layout of the same record: %.0f ms (%.1fx), %.1f MB (%.1fx)")
    :format(ems, ems / wms, emb, emb / wmb))
end)

case("scrolling_and_jumping_stay_cheap", function()
  local rec = synthetic()
  local lay = linear.build(rec, opts(true), WIDTH, {})
  local page = 40
  local _, one = cost(function() return lay:ensure(lay.header_rows, lay.header_rows + page) end)
  H.ok(one < 150, ("one screen took %.0f ms"):format(one))

  local _, scroll = cost(function()
    local row = lay.header_rows
    for _ = 1, 200 do lay:ensure(row, row + page); row = row + page end
  end)
  H.ok(scroll < 3000, ("200 pages took %.0f ms"):format(scroll))
  local live = body_rows_laid_out(lay)
  H.ok(live < #lay.lines / 4, ("%d of %d rows still laid out after scrolling"):format(live, #lay.lines))

  -- a position query is answered from the row plans, so jumping anywhere in
  -- the record costs nothing and builds nothing
  local before = vim.tbl_count(lay._built)
  local where = jumps(400, 7)
  local _, jump_ms = cost(function()
    for _, bp in ipairs(where) do
      local row, col = lay:pos_of(bp)
      H.eq(bp, lay:bp_at(row, col), "bp_at did not invert pos_of at " .. bp)
    end
  end)
  H.ok(jump_ms < 300, ("400 jumps took %.0f ms"):format(jump_ms))
  H.eq(before, vim.tbl_count(lay._built), "position queries built blocks")

  -- asking for a row's metadata does build its block, and that has a cost too
  local cold_at = jumps(200, 11)
  local _, cold = cost(function()
    for _, bp in ipairs(cold_at) do
      local row = lay:pos_of(bp)
      H.ok(lay:info(row) ~= nil, "no row metadata at " .. row)
    end
  end)
  H.ok(cold < 1500, ("200 cold rows took %.0f ms"):format(cold))
  H.log(("    [budget] one screen %.1f ms, 200 pages %.0f ms (%d of %d rows retained), 400 jumps %.1f ms, 200 cold rows %.0f ms")
    :format(one, scroll, live, #lay.lines, jump_ms, cold))
end)

-- A position query plans the block it lands in and keeps the plan, so visiting
-- enough of a record builds up a plan for every block of it even though not one
-- row was laid out. Scrolling has to take that back.
case("jumping_across_the_record_does_not_accumulate_plans", function()
  local rec = synthetic()
  local lay = linear.build(rec, opts(true), WIDTH, {})
  collectgarbage("collect")
  local base = collectgarbage("count") / 1024
  for b = 0, lay.nblocks - 1 do
    H.ok(lay:pos_of(b * lay.W + 1) ~= nil, "no row for block " .. b)
  end
  collectgarbage("collect")
  local plans, grew = vim.tbl_count(lay._plans), collectgarbage("count") / 1024 - base
  H.eq(lay.nblocks, plans, "a position query in every block should plan every block once")
  H.eq(0, vim.tbl_count(lay._built), "position queries built blocks")

  lay:ensure(lay.header_rows, lay.header_rows + 40)
  collectgarbage("collect")
  local after, left = vim.tbl_count(lay._plans), collectgarbage("count") / 1024 - base
  H.ok(after < lay.nblocks / 8, ("%d of %d plans survived a screen of scrolling"):format(after, lay.nblocks))
  H.ok(left < grew / 2, ("%.1f MB of the %.1f MB the jumps cost is still held"):format(left, grew))
  H.log(("    [budget] a position query in every one of %d blocks: %d plans, +%.1f MB; after one screen %d plans, %+.1f MB")
    :format(lay.nblocks, plans, grew, after, left))
end)

-- The same record through the buffer: this is what a user's keypress costs,
-- buffer edits and cursor work included.
case("a_large_record_opens_scrolls_and_resizes_in_a_buffer", function()
  local rec = synthetic()
  local dir = vim.fn.tempname()
  vim.fn.mkdir(dir, "p")
  local path = dir .. "/synth.gb"
  H.ok(require("plasmid.io").write(path, { rec }, "genbank", {}, "genbank"), "could not write the record")
  rec = nil

  local prev = config.options.linear.window
  config.options.linear.window = { enabled = "auto", threshold = 300000, margin = 2, cache_blocks = 512 }
  collectgarbage("collect")
  local rss0 = vim.uv.resident_set_memory() / 1048576
  local _, open_ms, open_mb = cost(function()
    vim.cmd("silent edit " .. vim.fn.fnameescape(path))
  end)
  local buf = vim.api.nvim_get_current_buf()
  local win = vim.api.nvim_get_current_win()
  local st = buffer.get(buf)
  H.ok(st and st.layout.windowed, "a 400 kb file should open windowed")
  H.ok(open_ms < 3000, ("opening took %.0f ms"):format(open_ms))
  H.ok(open_mb < 60, ("opening kept %.1f MB of heap"):format(open_mb))

  local nrows = vim.api.nvim_buf_line_count(buf)
  local height = math.max(5, vim.fn.winheight(win))
  local function scroll_to(top)
    top = math.max(1, math.min(top, math.max(1, nrows - height)))
    vim.fn.winrestview({ topline = top, lnum = math.min(top + 1, nrows), col = 0 })
    buffer.fill(buf, win)
  end
  local _, scroll_ms = cost(function()
    local top = math.floor(nrows / 4)
    for i = 1, 200 do scroll_to(top + i * height) end
  end)
  H.ok(scroll_ms < 3000, ("200 pages through the buffer took %.0f ms"):format(scroll_ms))

  local where = jumps(40, 3)
  local _, jump_ms = cost(function()
    for _, bp in ipairs(where) do
      buffer.goto_bp(win, bp)
      buffer.fill(buf, win)
      H.eq(bp, (buffer.cursor_info(win)), "cursor base pair after jumping to " .. bp)
    end
  end)
  H.ok(jump_ms < 2000, ("40 jumps through the buffer took %.0f ms"):format(jump_ms))

  local _, resize_ms = cost(function()
    vim.o.columns = vim.o.columns - 20
    buffer.render(buf)
    vim.o.columns = vim.o.columns + 20
    buffer.render(buf)
  end)
  H.ok(resize_ms < 4000, ("two resizes took %.0f ms"):format(resize_ms))

  collectgarbage("collect")
  local rss = vim.uv.resident_set_memory() / 1048576 - rss0
  H.ok(rss < 250, ("the process grew by %.0f MB"):format(rss))
  H.log(("    [budget] 400 kb buffer: open %.0f ms / %.1f MB, 200 pages %.0f ms, 40 jumps %.0f ms, resize %.0f ms, rss +%.0f MB")
    :format(open_ms, open_mb, scroll_ms, jump_ms, resize_ms, rss))
  config.options.linear.window = prev
  vim.cmd("silent bwipeout!")
  os.remove(path)
end)

return T
