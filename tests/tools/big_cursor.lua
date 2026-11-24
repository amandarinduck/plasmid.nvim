-- Dev tool: does a jump on a real genome land where it was asked to, with the
-- row actually in the buffer?
--
--   PLASMID_CHECK=FILE nvim --headless --clean -u NONE \
--     -c 'luafile tests/tools/big_cursor.lua'
--
-- tests/tools/big_parity.py compares two layouts; this one goes through
-- plasmid.buffer, where the fill path decides which rows exist at all. For
-- each of 309 base pairs -- the ends of the record, the middle, and a spread
-- of pseudo-random positions -- it asks for the cursor back, checks the row is
-- really in the buffer and matches the layout, and then checks every feature
-- the track rows of that block claim is under the cursor's column really does
-- cover that base. Windowing is forced on so a record under the threshold is
-- exercised too.
--
-- The arguments come from the environment because the window has to be 120
-- columns wide, and `nvim -l` never attaches a UI to resize (see
-- tests/tools/big_budget.lua).
local root = vim.fn.getcwd()
package.path = table.concat({ root .. "/lua/?.lua", root .. "/lua/?/init.lua", package.path }, ";")
vim.opt.runtimepath:prepend(root)
vim.notify = function() end
vim.o.columns, vim.o.lines = 120, 45
require("plasmid").setup({ linear = { window = { enabled = true, margin = 2, cache_blocks = 512 } } })

local path = os.getenv("PLASMID_CHECK")
vim.cmd("silent edit " .. vim.fn.fnameescape(path))
local buf = vim.api.nvim_get_current_buf()
local win = vim.api.nvim_get_current_win()
local B = require("plasmid.buffer")
local st = B.get(buf)
local rec = B.record(st)
local L = rec:len()
if not st.layout.windowed then io.stdout:write("not windowed\n"); vim.cmd("qa!") end

local targets = { 1, 2, L, L - 1, L - 60, 3, L / 2, 1000000, 4000000 }
math.randomseed(2024)
for _ = 1, 300 do targets[#targets + 1] = math.random(1, L) end
local bad, checked, feats = 0, 0, 0
for _, bpf in ipairs(targets) do
  local bp = math.max(1, math.min(math.floor(bpf), L))
  B.goto_bp(win, bp)
  B.fill(buf, win)
  local got, fid = B.cursor_info(win)
  local row = vim.api.nvim_win_get_cursor(win)[1]
  local line = vim.api.nvim_buf_get_lines(buf, row - 1, row, false)[1]
  local problems = {}
  if got ~= bp then problems[#problems + 1] = ("cursor reports %s"):format(tostring(got)) end
  if line == "" then problems[#problems + 1] = "the row is blank in the buffer" end
  if line ~= st.layout.lines[row] then problems[#problems + 1] = "buffer text differs from the layout" end
  -- the cursor sits on the sequence row, where no feature is drawn; the track
  -- rows of the same block are where feature_at has something to say, and each
  -- answer has to be a feature that really covers the base under that column
  local _, vcol = st.layout:pos_of(bp)
  local b = st.layout:block_at_row(row - 1)
  for r = st.layout.offsets[b + 1], st.layout.offsets[b + 1] + st.layout:plan_of(b).n - 1 do
    local inf = st.layout:info(r)
    if inf and inf.kind == "track" then
      local tf = st.layout:feature_at(r, vcol)
      if tf then
        feats = feats + 1
        local f = rec:feature_by_id(tf)
        local covers = false
        for _, e in ipairs(f and rec:extents(f) or {}) do
          if bp >= e[1] and bp <= e[2] then covers = true end
        end
        if not covers then
          problems[#problems + 1] = ("feature %s on row %d does not cover %d"):format(tostring(tf), r, bp)
        end
      end
    end
  end
  if fid then feats = feats + 1 end
  checked = checked + 1
  if #problems > 0 then
    bad = bad + 1
    if bad <= 5 then
      io.stdout:write(("bp %d row %d: %s\n"):format(bp, row, table.concat(problems, "; ")))
    end
  end
end
io.stdout:write(("%s: %d jumps checked, %d landed on a feature, %d problems\n")
  :format(rec.name, checked, feats, bad))
vim.cmd("qa!")
