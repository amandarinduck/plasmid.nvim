local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local config = require("plasmid.config")
local sel = require("plasmid.select")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.get(buf), buffer.record(buffer.get(buf))
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end
local function flush() vim.wait(50, function() return false end) end -- run scheduled callbacks
local function last_msg() local m = H.messages[#H.messages]; return m and m.msg or "", m and m.level end

-- Record clipboard writes instead of touching a real clipboard.
local calls, result = {}, "provider"
local real_write = sel.clipboard_write
local function stub()
  calls = {}
  sel.clipboard_write = function(reg, text) calls[#calls + 1] = { reg = reg, text = text }; return result end
end
local function unstub() sel.clipboard_write = real_write end

vim.o.columns, vim.o.lines = 120, 50

local function with(opts, fn)
  local saved_cb, saved_opt = vim.o.clipboard, config.options.yank_clipboard
  if opts.clipboard then vim.o.clipboard = opts.clipboard end
  if opts.yank_clipboard ~= nil then config.options.yank_clipboard = opts.yank_clipboard end
  local ok, err = pcall(fn)
  vim.o.clipboard, config.options.yank_clipboard = saved_cb, saved_opt
  if not ok then error(err, 0) end
end

case("plain_yank_goes_to_clipboard_by_default", function()
  local _, _, rec = open("ncbi/J01749.gb")
  stub()
  with({ clipboard = "" }, function()
    buffer.goto_bp(0, 101)
    keys("v9ly")
    flush()
    local want = rec.seq:sub(101, 110):upper()
    H.eq(want, vim.fn.getreg('"'))
    H.eq(1, #calls, "one clipboard write")
    H.eq("+", calls[1].reg)
    H.eq(want, calls[1].text)
    H.ok(last_msg():find("and the clipboard", 1, true), last_msg())
    -- reverse complement and picker yanks take the same path
    buffer.goto_bp(0, 101)
    keys("v9lY")
    flush()
    H.eq("+", calls[2].reg)
    H.eq(calls[2].text, vim.fn.getreg('"'))
  end)
  unstub()
end)

case("named_register_stays_in_neovim", function()
  local _, _, rec = open("ncbi/J01749.gb")
  stub()
  with({ clipboard = "" }, function()
    vim.fn.setreg("a", "")
    buffer.goto_bp(0, 201)
    keys('v4l"ay')
    flush()
    H.eq(rec.seq:sub(201, 205):upper(), vim.fn.getreg("a"))
    H.eq(0, #calls, "no clipboard write for \"a")
    -- an explicit "+ does go to the clipboard
    buffer.goto_bp(0, 201)
    keys('v4l"+y')
    flush()
    H.eq(1, #calls)
    H.eq("+", calls[1].reg)
  end)
  unstub()
end)

case("option_values", function()
  open("ncbi/J01749.gb")
  stub()
  local function yank_once()
    calls = {}
    vim.cmd("DnaSelect 10..20")
    vim.cmd("DnaYank top")
    flush()
    return calls[1] and calls[1].reg or nil
  end
  with({ yank_clipboard = false, clipboard = "unnamedplus" }, function() H.eq(nil, yank_once(), "false") end)
  with({ yank_clipboard = "auto", clipboard = "" }, function() H.eq(nil, yank_once(), "auto, clipboard=''") end)
  with({ yank_clipboard = "auto", clipboard = "unnamed" }, function() H.eq("*", yank_once(), "auto, unnamed") end)
  with({ yank_clipboard = "auto", clipboard = "unnamedplus" }, function() H.eq("+", yank_once(), "auto, unnamedplus") end)
  with({ yank_clipboard = true, clipboard = "" }, function() H.eq("+", yank_once(), "true") end)
  -- :DnaYank with a register argument
  with({ yank_clipboard = true, clipboard = "" }, function()
    calls = {}
    vim.cmd("DnaYank bottom a")
    flush()
    H.eq(0, #calls, ":DnaYank ... a")
    vim.cmd("DnaYank bottom +")
    flush()
    H.eq("+", calls[1].reg, ":DnaYank ... +")
  end)
  unstub()
end)

case("no_clipboard_warns_once", function()
  open("ncbi/J01749.gb")
  stub()
  result = nil
  sel._reset_clipboard_warning()
  with({ clipboard = "" }, function()
    vim.cmd("DnaSelect 10..20")
    H.messages = {}
    vim.cmd("DnaYank top")
    flush()
    local msg, level = last_msg()
    H.ok(msg:find("no clipboard", 1, true), msg)
    H.eq(vim.log.levels.WARN, level)
    H.eq(buffer.record(buffer.get(0)).seq:sub(10, 20):upper(), vim.fn.getreg('"'), "unnamed register still set")
    vim.cmd("DnaYank top")
    flush()
    local _, level2 = last_msg()
    H.eq(vim.log.levels.INFO, level2, "second time not a warning")
  end)
  result = "provider"
  unstub()
end)

case("yanks_follow_displayed_case", function()
  stub()
  -- GenBank: stored lower case, shown and yanked upper case
  local _, _, rec = open("ncbi/J01749.gb")
  H.ok(rec.seq:sub(10, 20):find("%l"), "fixture stores lower case")
  vim.cmd("DnaSelect 10..20")
  vim.cmd("DnaYank top")
  H.eq(rec.seq:sub(10, 20):upper(), vim.fn.getreg('"'))
  -- linear.case = "preserve" keeps the file's case
  config.options.linear.case = "preserve"
  vim.cmd("DnaYank top")
  H.eq(rec.seq:sub(10, 20), vim.fn.getreg('"'))
  config.options.linear.case = "auto"
  -- FASTA: soft-masked (lower-case) bases are kept, as displayed
  local path = tmp .. "/masked.fa"
  local fh = assert(io.open(path, "w")); fh:write(">masked\nACGTacgtNNacgtACGT\n"); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(path))
  vim.cmd("DnaSelect 3..14")
  vim.cmd("DnaYank top")
  H.eq("GTacgtNNacgt", vim.fn.getreg('"'))
  vim.cmd("DnaYank bottom")
  H.eq("acgtNNacgtAC", vim.fn.getreg('"'), "reverse complement keeps case per base")
  flush()
  unstub()
end)

case("real_writer_without_provider_or_ui", function()
  -- headless test Neovim: no clipboard provider, no UI attached
  if vim.fn.has("clipboard") == 0 and #vim.api.nvim_list_uis() == 0 then
    H.eq(nil, real_write("+", "ACGT"))
  end
end)

return T
