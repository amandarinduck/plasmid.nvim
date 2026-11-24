-- What the overview costs on a real genome, end to end.
--   nvim -l tests/tools/overview_bench.lua big/NC_000913.gbk
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({ root .. "/lua/?.lua", root .. "/lua/?/init.lua", package.path }, ";")
vim.opt.runtimepath:prepend(root)
require("plasmid").setup({})
local io_ = require("plasmid.io")
local config = require("plasmid.config")
local ov = require("plasmid.view.overview")
local buffer = require("plasmid.buffer")
local overview = require("plasmid.overview")

local path = ({ ... })[1] or (root .. "/big/NC_000913.gbk")
local COLS, ROWS = 140, 32

local function heap() collectgarbage("collect"); return collectgarbage("count") / 1024 end
local function ms(t0) return (vim.uv.hrtime() - t0) / 1e6 end
local function say(...) io.stdout:write(("%s\n"):format(table.concat({ ... }, ""))) end

-- parse
local h0 = heap()
local t0 = vim.uv.hrtime()
local rec = assert(io_.read(path)).records[1]
say(("parse            %7.0f ms   %6.1f MB   %d bp, %d features")
  :format(ms(t0), heap() - h0, rec:len(), #rec.features))

-- the layout on its own, three times (the first pays for the feature list)
for i = 1, 3 do
  local h1 = heap()
  local t1 = vim.uv.hrtime()
  local lay = ov.build(rec, config.options, COLS, ROWS, { cursor_bp = 2000000, visible = { 2000000, 2005000 } })
  say(("overview build %d %7.1f ms   %6.1f MB   %d rows of %d bp, %d cells")
    :format(i, ms(t1), heap() - h1, lay.nrows, lay.span, lay.cells))
end

-- in a window, the way a user gets there: open the record, then :DnaOverview
vim.o.lines, vim.o.columns = ROWS, COLS
local h2 = heap()
local t2 = vim.uv.hrtime()
vim.cmd("silent edit " .. vim.fn.fnameescape(path))
local buf = vim.api.nvim_get_current_buf()
local win = vim.api.nvim_get_current_win()
local open_ms, open_mb = ms(t2), heap() - h2
local st = buffer.get(buf)
say(("open record      %7.0f ms   %6.1f MB   %s layout, %d rows, %d blocks")
  :format(open_ms, open_mb, st.layout.windowed and "windowed" or "eager", #st.layout.lines, st.layout.nblocks or 0))

buffer.goto_bp(win, 2000000)
local built_before = 0
for _ in pairs(st.layout._built) do built_before = built_before + 1 end

local h3 = heap()
local t3 = vim.uv.hrtime()
vim.cmd("DnaOverview")
local obuf = vim.api.nvim_win_get_buf(win)
say((":DnaOverview     %7.1f ms   %6.1f MB   %d rows on screen")
  :format(ms(t3), heap() - h3, #overview.state[obuf].lay.lines))

local lay = overview.state[obuf].lay
local row = lay.header_rows + math.floor(lay.nrows * 0.75)
local t4 = vim.uv.hrtime()
overview.jump(obuf, { row, 0 })
local jump_ms = ms(t4)
local built_after = 0
for _ in pairs(st.layout._built) do built_after = built_after + 1 end
say(("<CR> jump        %7.1f ms   blocks laid out %d → %d of %d, cursor at %d")
  :format(jump_ms, built_before, built_after, st.layout.nblocks, (buffer.cursor_info(win))))

local laid = 0
for i = st.layout.header_rows + 1, #st.layout.lines do
  if st.layout.lines[i] ~= "" then laid = laid + 1 end
end
say(("after the jump   %d of %d rows laid out, heap %.0f MB, rss %.0f MB")
  :format(laid, #st.layout.lines, heap(), (vim.uv.getrusage().maxrss or 0) / 1024))
