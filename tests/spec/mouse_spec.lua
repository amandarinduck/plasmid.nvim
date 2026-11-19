-- Mouse: drag to select, double click to select a feature, Ctrl-click two
-- primers to select their product. getmousepos() is stubbed with the screen
-- position of a known base pair, so the tests exercise the real mapping code.
local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local S = require("plasmid.select")
local mouse = require("plasmid.mouse")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

local real_getmousepos = vim.fn.getmousepos
local function open(rel, as)
  local dst = tmp .. "/" .. (as or vim.fn.fnamemodify(rel, ":t"))
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent! only")
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if st.view ~= "linear" then buffer.set_view(buf, "linear") end
  return buf, st, buffer.record(st)
end

-- Put the mouse on a base pair, as the terminal would report it.
local function at(st, win, bp)
  local row0, vcol = st.layout:pos_of(bp)
  local byte = st.layout.byte_col and st.layout:byte_col(row0, vcol) or vcol
  vim.fn.getmousepos = function()
    return { winid = win, line = row0 + 1, column = byte + 1, screenrow = row0 + 1, screencol = byte + 1 }
  end
end
local function messages(fn)
  local msgs = {}
  local notify = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  fn()
  vim.notify = notify
  return table.concat(msgs, " | ")
end

case("dragging_selects_a_region", function()
  local buf, st = open("ncbi/J01749.gb")
  local win = vim.api.nvim_get_current_win()
  at(st, win, 100); mouse.press()
  H.eq(nil, S.get(buf), "pressing clears the previous selection")
  at(st, win, 400); mouse.drag()
  local sel = S.get(buf)
  H.eq(100, sel.start)
  H.eq(400, sel.stop)
  -- dragging backwards selects the same span
  at(st, win, 50); mouse.drag()
  sel = S.get(buf)
  H.eq(50, sel.start)
  H.eq(100, sel.stop)
  local msg = messages(function() mouse.release() end)
  H.ok(msg:find("selected 50%.%.100") and msg:find("51 bp"), msg)
  -- and the cursor followed the mouse to where it was let go
  H.eq(50, (buffer.cursor_info(win)))
end)

case("dragging_works_on_the_circular_map", function()
  local buf, st = open("ncbi/L09137.gb")
  buffer.set_view(buf, "circular")
  st = buffer.get(buf)
  local win = vim.api.nvim_get_current_win()
  at(st, win, 200); mouse.press()
  at(st, win, 1200); mouse.drag()
  local sel = S.get(buf)
  H.ok(sel ~= nil, "a drag on the ring selects")
  H.ok(math.abs(sel.start - 200) < 30 and math.abs(sel.stop - 1200) < 30,
    ("selected %d..%d"):format(sel.start, sel.stop))
  buffer.set_view(buf, "linear")
end)

case("double_click_selects_the_feature", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local win = vim.api.nvim_get_current_win()
  at(st, win, 500)
  local msg = messages(function() mouse.double_click() end)
  local sel = S.get(buf)
  H.eq(86, sel.start)
  H.eq(1276, sel.stop)
  H.ok(msg:find("tet"), msg)
  -- empty space says so
  at(st, win, rec:len())
  local msg2 = messages(function() mouse.double_click() end)
  H.ok(msg2 == "" or msg2:find("selected") or msg2:find("no feature"), msg2)
end)

case("ctrl_clicking_two_primers_selects_the_product", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local win = vim.api.nvim_get_current_win()
  local tmpl = rec.seq:tostring()
  local primers = require("plasmid.primers")
  primers.add_sequence(buf, tmpl:sub(200, 222), "F1")
  primers.add_sequence(buf, require("plasmid.bio.primer").revcomp(tmpl:sub(700, 722)), "R1")
  st = buffer.get(buf)
  at(st, win, 210)
  local first = messages(function() mouse.ctrl_click() end)
  H.ok(first:find("F1 picked"), first)
  at(st, win, 710)
  local second = messages(function() mouse.ctrl_click() end)
  local sel = S.get(buf)
  H.eq(200, sel.start)
  H.eq(722, sel.stop)
  H.ok(second:find("523 bp product") and second:find("Ta"), second)
  -- the order of the two clicks does not matter
  at(st, win, 710); messages(function() mouse.ctrl_click() end)
  at(st, win, 210); messages(function() mouse.ctrl_click() end)
  sel = S.get(buf)
  H.eq(200, sel.start)
  H.eq(722, sel.stop)
  -- clicking something that is not a primer says so and forgets the pick
  at(st, win, 1500)
  local none = messages(function() mouse.ctrl_click() end)
  H.ok(none:find("no primer there"), none)
end)

case("primers_that_make_no_product_are_reported", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local win = vim.api.nvim_get_current_win()
  local tmpl = rec.seq:tostring()
  local P = require("plasmid.bio.primer")
  local primers = require("plasmid.primers")
  -- both on the same strand: nothing to amplify
  primers.add_sequence(buf, tmpl:sub(200, 222), "A1")
  primers.add_sequence(buf, tmpl:sub(700, 722), "A2")
  st = buffer.get(buf)
  S.set(buf, nil)
  at(st, win, 210); messages(function() mouse.ctrl_click() end)
  at(st, win, 710)
  local msg = messages(function() mouse.ctrl_click() end)
  H.ok(msg:find("no product"), msg)
  H.eq(nil, S.get(buf), "nothing is selected when there is no product")
  H.ok(P.valid(tmpl:sub(200, 222)))
end)

case("mouse_maps_follow_the_config", function()
  local buf = open("ncbi/J01749.gb")
  local function lhs_set()
    local out = {}
    for _, m in ipairs(vim.api.nvim_buf_get_keymap(buf, "n")) do out[m.lhs] = true end
    return out
  end
  local maps = lhs_set()
  for _, k in ipairs({ "<LeftMouse>", "<LeftDrag>", "<LeftRelease>", "<2-LeftMouse>", "<C-LeftMouse>" }) do
    H.ok(maps[k], "missing mapping " .. k)
  end
  local cfg = require("plasmid.config")
  local saved = cfg.options.mouse
  cfg.options.mouse = { enabled = false }
  local buf2 = open("ncbi/L09137.gb", "no_mouse.gb")
  local maps2 = {}
  for _, m in ipairs(vim.api.nvim_buf_get_keymap(buf2, "n")) do maps2[m.lhs] = true end
  H.eq(nil, maps2["<LeftDrag>"], "mouse.enabled = false leaves the mouse alone")
  cfg.options.mouse = saved
end)

case("report_buffers_search_without_case", function()
  local buf = open("ncbi/J01749.gb")
  local saved = vim.o.ignorecase
  vim.o.ignorecase = false
  require("plasmid.digest").run(buf, { "EcoRI", "HindIII" })
  local dbuf
  for b, s in pairs(require("plasmid.digest").state) do
    if s.src == buf and vim.api.nvim_buf_is_valid(b) then dbuf = b end
  end
  H.ok(dbuf ~= nil, "the digest table opened")
  local rhs = {}
  for _, m in ipairs(vim.api.nvim_buf_get_keymap(dbuf, "n")) do rhs[m.lhs] = m.rhs end
  H.eq("/\\c", rhs["/"], "/ in a report buffer ignores case")
  H.eq("?\\c", rhs["?"], "? in a report buffer ignores case")
  vim.o.ignorecase = saved
end)

case("teardown", function()
  vim.fn.getmousepos = real_getmousepos
  vim.cmd("silent! only")
end)

return T
