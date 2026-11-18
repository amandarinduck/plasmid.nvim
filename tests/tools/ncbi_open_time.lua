-- Dev tool: time :edit (read + parse + first render) of one file in a headless
-- Neovim with the plugin loaded, then the other view, and report memory.
--   nvim --headless --clean -u NONE -l tests/tools/ncbi_open_time.lua FILE [OUT.json]
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
vim.opt.runtimepath:prepend(root)
package.path = root .. "/lua/?.lua;" .. root .. "/lua/?/init.lua;" .. package.path
vim.o.columns, vim.o.lines = 120, 40
local notes = {}
vim.notify = function(msg, level) notes[#notes + 1] = { msg = tostring(msg), level = level } end
require("plasmid").setup({})
local path = arg[1]
collectgarbage("collect")
local rss0, heap0 = vim.uv.resident_set_memory(), collectgarbage("count")
local t0 = vim.uv.hrtime()
local ok, err = pcall(vim.cmd, "silent edit " .. vim.fn.fnameescape(path))
local open_ms = (vim.uv.hrtime() - t0) / 1e6
local out = { file = vim.fn.fnamemodify(path, ":t"), open_ok = ok, open_err = (not ok) and tostring(err) or nil,
  open_ms = open_ms }
local buf = vim.api.nvim_get_current_buf()
local B = require("plasmid.buffer")
local st = B.state[buf]
if st then
  out.view = st.view
  out.lines = vim.api.nvim_buf_line_count(buf)
  local t1 = vim.uv.hrtime()
  local ok2, err2 = pcall(B.toggle_view, buf)
  out.toggle_ms = (vim.uv.hrtime() - t1) / 1e6
  out.toggle_ok, out.toggle_err = ok2, (not ok2) and tostring(err2) or nil
  out.view2 = st.view
  -- force a redraw pass so the decoration provider runs once
  local t2 = vim.uv.hrtime()
  local ok3, err3 = pcall(vim.cmd, "redraw!")
  out.redraw_ms = (vim.uv.hrtime() - t2) / 1e6
  out.redraw_ok, out.redraw_err = ok3, (not ok3) and tostring(err3) or nil
else
  out.no_state = true
end
collectgarbage("collect")
out.heap_mb = (collectgarbage("count") - heap0) / 1024
out.rss_mb = vim.uv.resident_set_memory() / 1048576
out.rss_delta_mb = (vim.uv.resident_set_memory() - rss0) / 1048576
out.notes = notes
local js = vim.json.encode(out)
if arg[2] then local fh = io.open(arg[2], "ab"); fh:write(js, "\n"); fh:close() end
io.stdout:write(js, "\n")
vim.cmd("qa!")
