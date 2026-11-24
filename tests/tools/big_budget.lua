-- Dev tool: what a large record costs through the real buffer path -- opening
-- it, scrolling it, resizing the window, jumping across it -- in a headless
-- Neovim with the plugin loaded.
--
--   PLASMID_BUDGET="FILE [key=value ...] [OUT.jsonl]" \
--     nvim --headless --clean -u NONE -c 'luafile tests/tools/big_budget.lua'
--
-- keys: mode=window|eager  width=120  height=45  rec=1  pages=300  jumps=20
--
-- The arguments come from the environment because `nvim -l` is the one way of
-- running a script that never attaches a UI: its window stays 80 columns wide
-- whatever 'columns' is set to, and the line width of the view -- which
-- decides how many blocks a record has -- would not be the one asked for.
-- Under `--headless -c` the dummy UI is real enough to resize.
--
-- Every phase goes through plasmid.buffer, not through the layout directly, so
-- the numbers include the buffer edits the fill path makes and the cursor work
-- the autocmds do. Scrolling is emulated with winrestview, which moves topline
-- and botline in `nvim -l` exactly as a real scroll does; nothing here needs a
-- UI, and nothing calls :redraw, which is not safe headless on 0.12.2.
--
-- One record and one mode per process: the eager layout of a bacterial genome
-- is hundreds of megabytes and must not be live while the windowed one is
-- measured.
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({ root .. "/lua/?.lua", root .. "/lua/?/init.lua", package.path }, ";")
vim.opt.runtimepath:prepend(root)
vim.notify = function() end

local argv = arg
if not argv or #argv == 0 then
  argv = vim.split(os.getenv("PLASMID_BUDGET") or "", "%s+", { trimempty = true })
end
local path = argv[1]
local A, out_path = {}, nil
for i = 2, #argv do
  local k, v = argv[i]:match("^([%w_]+)=(.*)$")
  if k then A[k] = v else out_path = argv[i] end
end
local function num(k, d) return tonumber(A[k]) or d end

local windowed = (A.mode or "window") == "window"
vim.o.columns, vim.o.lines = num("width", 120), num("height", 40)
require("plasmid").setup({
  linear = { window = { enabled = windowed, threshold = 300000, margin = 2, cache_blocks = 512 } },
})
local buffer = require("plasmid.buffer")

local peak_heap, peak_rss = 0, 0
local function mark()
  peak_heap = math.max(peak_heap, collectgarbage("count") / 1024)
  peak_rss = math.max(peak_rss, vim.uv.resident_set_memory() / 1048576)
end
local function ms(fn)
  local t = vim.uv.hrtime()
  fn()
  local d = (vim.uv.hrtime() - t) / 1e6
  mark()
  return d
end

--- Median and maximum of a sample, which is what a latency budget is about:
--- the typical keypress and the worst one.
local function stats(xs)
  table.sort(xs)
  local n = #xs
  if n == 0 then return 0, 0 end
  return xs[math.ceil(n / 2)], xs[n]
end

collectgarbage("collect")
local rss0 = vim.uv.resident_set_memory() / 1048576
local out = { file = vim.fn.fnamemodify(path, ":t"), mode = windowed and "window" or "eager",
  width = vim.o.columns, lines = vim.o.lines }

-- opening: read, parse, lay out, fill the first screen
out.open_ms = ms(function() vim.cmd("silent edit " .. vim.fn.fnameescape(path)) end)
local buf = vim.api.nvim_get_current_buf()
local win = vim.api.nvim_get_current_win()
local st = buffer.get(buf)
if not st then error("not a plasmid buffer: " .. tostring(path)) end
if num("rec", 1) > 1 then
  st.current = num("rec", 1)
  out.open_ms = out.open_ms + ms(function() buffer.render(buf, { keep_bp = 1 }) end)
end
local rec = buffer.record(st)
local lay = st.layout
out.name, out.len, out.features = rec.name, rec:len(), #rec.features
out.topology, out.rows = rec.topology, #lay.lines
out.windowed = lay.windowed == true
if out.windowed ~= windowed then error("layout mode is not the one asked for") end
out.nblocks = lay.nblocks
-- the width the view actually got, and the residues per line it chose: both
-- are what the row count is made of, so a measurement is only readable with them
out.text_width, out.W = lay.width, lay.W
out.heap_after_open_mb = collectgarbage("count") / 1024
out.rss_after_open_mb = vim.uv.resident_set_memory() / 1048576
mark()

local height = math.max(5, vim.fn.winheight(win))
local nrows = vim.api.nvim_buf_line_count(buf)
out.win_height = height

--- Scroll to a top line the way a key would: move the view, let the cursor
--- follow it, then give the plugin the chance its autocmds would give it.
local function scroll_to(top)
  top = math.max(1, math.min(top, math.max(1, nrows - height)))
  vim.fn.winrestview({ topline = top, lnum = math.min(top + 1, nrows), col = 0 })
  buffer.fill(buf, win)
end

-- one page, measured from a cold window position each time
local page = {}
for i = 1, 10 do
  local top = math.floor(nrows * (0.1 + 0.07 * i))
  scroll_to(top)
  page[#page + 1] = ms(function() scroll_to(top + height) end)
end
out.page_ms, out.page_max_ms = stats(page)

-- a long scroll: consecutive pages, the way holding the key down feels
local pages = num("pages", 300)
local start_top = math.floor(nrows / 3)
scroll_to(start_top)
local each = {}
out.long_scroll_ms = ms(function()
  for i = 1, pages do
    local t = vim.uv.hrtime()
    scroll_to(start_top + i * height)
    each[#each + 1] = (vim.uv.hrtime() - t) / 1e6
  end
end)
out.long_pages = pages
out.long_rows = pages * height
out.long_page_ms, out.long_page_max_ms = stats(each)

-- jumping to a base pair anywhere in the record
local jumps, jms = num("jumps", 20), {}
math.randomseed(20250101)
for _ = 1, jumps do
  local bp = math.random(1, math.max(1, rec:len()))
  jms[#jms + 1] = ms(function()
    buffer.goto_bp(win, bp)
    buffer.fill(buf, win)
  end)
end
out.jump_ms, out.jump_max_ms = stats(jms)
out.jumps = jumps

-- a resize re-lays out the record at the new width
out.resize_ms = ms(function()
  vim.o.columns = vim.o.columns - 20
  buffer.render(buf)
end)
out.resize_back_ms = ms(function()
  vim.o.columns = vim.o.columns + 20
  buffer.render(buf)
end)
out.rows_after_resize = #st.layout.lines

-- how much of the record is laid out after all that
if st.layout.windowed then
  local built = 0
  for _ in pairs(st.layout._built) do built = built + 1 end
  out.blocks_built = built
  local nonblank = 0
  for i = st.layout.header_rows + 1, #st.layout.lines do
    if st.layout.lines[i] ~= "" then nonblank = nonblank + 1 end
  end
  out.rows_live = nonblank
end

mark()
out.peak_heap_mb, out.peak_rss_mb = peak_heap, peak_rss
collectgarbage("collect")
out.heap_mb = collectgarbage("count") / 1024
out.rss_mb = vim.uv.resident_set_memory() / 1048576
out.rss_delta_mb = out.rss_mb - rss0

local js = vim.json.encode(out)
if out_path then local fh = io.open(out_path, "ab"); fh:write(js, "\n"); fh:close() end
io.stdout:write(("%-16s %-7s %9d bp  open %7.0f ms  page %5.1f/%6.1f  long %6.1f/%7.1f  jump %5.1f/%6.1f  resize %7.0f  heap %5.0f  rss %5.0f MB\n")
  :format(out.file, out.mode, out.len, out.open_ms, out.page_ms, out.page_max_ms,
    out.long_page_ms, out.long_page_max_ms, out.jump_ms, out.jump_max_ms, out.resize_ms,
    out.peak_heap_mb, out.peak_rss_mb))
vim.cmd("qa!")
