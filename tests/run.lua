-- Zero-dependency test runner.
--   nvim -l tests/run.lua            run every spec
--   nvim -l tests/run.lua location   run specs whose filename matches
--   PLASMID_TEST_GROUP=0 nvim -l tests/run.lua   same, continuous sequence layout
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h")
package.path = table.concat({
  root .. "/lua/?.lua", root .. "/lua/?/init.lua", root .. "/tests/?.lua", package.path,
}, ";")
vim.opt.runtimepath:prepend(root)

local H = require("helpers")
H.root = root
local config = require("plasmid.config")
-- PLASMID_TEST_GROUP=0 runs every spec with that linear grouping as the default
-- (CI runs the suite with 10 and with 0)
local group = tonumber(os.getenv("PLASMID_TEST_GROUP") or "")
if group then
  config.defaults.linear.group = group
  config.options.linear.group = group
end

-- Capture notifications instead of printing them. Besides keeping output
-- readable, this sidesteps a Neovim 0.12 headless crash: without a UI,
-- successive messages accumulate on one line, and a redraw after that line
-- outgrows 'columns' segfaults. PLASMID_TEST_VERBOSE=1 prints them.
H.messages = {}
vim.notify = function(msg, level)
  H.messages[#H.messages + 1] = { msg = msg, level = level }
  if os.getenv("PLASMID_TEST_VERBOSE") then io.stdout:write("  [notify] ", tostring(msg), "\n") end
end

local filter = arg[1]
local files = vim.fn.globpath(root .. "/tests/spec", "*_spec.lua", false, true)
table.sort(files)

local passed, failed, failures = 0, 0, {}
local t0 = vim.uv.hrtime()
for _, file in ipairs(files) do
  local name = vim.fn.fnamemodify(file, ":t:r")
  if not filter or name:find(filter, 1, true) then
    -- each spec file starts from a clean editor: one window, no buffers, default config
    vim.cmd("silent! only")
    vim.cmd("silent! %bwipeout!")
    package.loaded["plasmid.config"].setup({})
    local cases = dofile(file)
    local order = cases._order or vim.tbl_keys(cases)
    if not cases._order then table.sort(order) end
    for _, case in ipairs(order) do
      if case:sub(1, 1) ~= "_" then
        -- PLASMID_TEST_TRACE=1: name each case before it runs (locates crashes)
        if os.getenv("PLASMID_TEST_TRACE") then io.stderr:write("  > ", name, " :: ", case, "\n") end
        local ok, err = xpcall(cases[case], debug.traceback)
        if ok then
          passed = passed + 1
        else
          failed = failed + 1
          failures[#failures + 1] = ("%s :: %s\n%s"):format(name, case, err)
        end
      end
    end
  end
end
local ms = (vim.uv.hrtime() - t0) / 1e6
for _, f in ipairs(failures) do io.stdout:write("FAIL ", f, "\n\n") end
io.stdout:write(("%d passed, %d failed (%.0f ms)\n"):format(passed, failed, ms))
os.exit(failed == 0 and 0 or 1)
