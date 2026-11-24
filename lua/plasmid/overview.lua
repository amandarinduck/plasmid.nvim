-- The overview as a thing you can steer with: open it over the record, pick
-- a row, land there in the sequence view.
--
-- It takes over the window the record is in rather than opening a split,
-- because its whole promise is that the record fits on one screen: half a
-- screen would halve the resolution of every row. The window remembers what
-- it was showing, so closing the overview puts the sequence back exactly
-- where it was, and jumping puts it somewhere else on purpose.
local config = require("plasmid.config")
local buffer = require("plasmid.buffer")
local view = require("plasmid.view.overview")

local M = {}

---@type table<integer, table> overview buffer -> { src, win, lay, restore }
M.state = {}

local ns = vim.api.nvim_create_namespace("plasmid_overview")
M.ns = ns

local function notify(msg, level)
  require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO)
end

local function win_for(buf)
  local cur = vim.api.nvim_get_current_win()
  if vim.api.nvim_win_get_buf(cur) == buf then return cur end
  return vim.fn.win_findbuf(buf)[1]
end

local function text_width(win)
  local info = vim.fn.getwininfo(win)[1]
  return info and (info.width - info.textoff) or vim.api.nvim_win_get_width(win)
end

--- The bases the sequence view is showing in `win`, as { first, last }.
--- Asked of the layout rather than of the record, so a windowed linear view
--- answers for rows it has not laid out yet.
local function visible_range(st, win)
  local lay = st.layout
  if not (lay and lay.bp_at and win and vim.api.nvim_win_is_valid(win)) then return nil end
  local info = vim.fn.getwininfo(win)[1]
  if not info then return nil end
  local lo, hi
  for row = (info.topline or 1) - 1, (info.botline or 1) - 1 do
    local ok, a = pcall(lay.bp_at, lay, row, lay.prefix_w or 0)
    local ok2, b = pcall(lay.bp_at, lay, row, (lay.width or 200) * 2)
    if ok and a then lo = lo and math.min(lo, a) or a end
    if ok2 and b then hi = hi and math.max(hi, b) or b end
  end
  if lo and hi then return { lo, hi } end
  return nil
end

local function paint(obuf, lay)
  vim.bo[obuf].modifiable = true
  vim.api.nvim_buf_set_lines(obuf, 0, -1, false, lay.lines)
  vim.bo[obuf].modifiable = false
  vim.bo[obuf].modified = false
  vim.api.nvim_buf_clear_namespace(obuf, ns, 0, -1)
  for _, h in ipairs(lay.hls) do
    pcall(vim.api.nvim_buf_set_extmark, obuf, ns, h[1], h[2], {
      end_col = h[3], hl_group = h[4], priority = h[5] or 100,
    })
  end
end

--- The overview open on `buf`, if there is one.
function M.for_record(buf)
  for b, s in pairs(M.state) do
    if s.src == buf and vim.api.nvim_buf_is_valid(b) then return b end
  end
  return nil
end

--- Put the record back in its window. Returns the window.
local function restore(s)
  local win = (s.win and vim.api.nvim_win_is_valid(s.win)) and s.win or win_for(s.src)
  if not win or not vim.api.nvim_buf_is_valid(s.src) then return nil end
  vim.api.nvim_win_set_buf(win, s.src)
  if vim.api.nvim_get_current_win() ~= win then pcall(vim.api.nvim_set_current_win, win) end
  return win
end

---------------------------------------------------------------------------
-- Opening and closing
---------------------------------------------------------------------------

---@param buf integer sequence buffer
---@param opts { span?: integer }|nil
function M.open(buf, opts)
  opts = opts or {}
  local st = buffer.get(buf)
  if not st then return notify("current buffer is not a sequence file", vim.log.levels.WARN) end
  local win = win_for(buf)
  if not win then return notify("that record is not in a window", vim.log.levels.WARN) end
  local open = M.for_record(buf)
  if open then M.close(open) end

  local rec = buffer.record(st)
  local o = config.options
  if opts.span then o = vim.tbl_deep_extend("force", o, { overview = { span = opts.span } }) end
  local bp = select(1, buffer.cursor_info(win)) or st.last_bp
  local lay = view.build(rec, o, text_width(win), vim.fn.winheight(win),
    { cursor_bp = bp, visible = visible_range(st, win) })

  local obuf = vim.api.nvim_create_buf(false, true)
  vim.bo[obuf].bufhidden = "wipe"
  vim.bo[obuf].filetype = "plasmid_overview"
  pcall(vim.api.nvim_buf_set_name, obuf, ("plasmid://overview/%s"):format(rec.name))
  paint(obuf, lay)
  M.state[obuf] = {
    src = buf, win = win, lay = lay, rec = rec,
    restore = vim.api.nvim_win_call(win, function() return vim.fn.winsaveview() end),
  }
  vim.api.nvim_create_autocmd("BufWipeout", { buffer = obuf, once = true,
    callback = function() M.state[obuf] = nil end })
  vim.api.nvim_win_set_buf(win, obuf)
  vim.api.nvim_win_call(win, function()
    vim.wo[0][0].winbar = ""
    vim.wo[0][0].wrap = false
    vim.wo[0][0].number = false
    vim.wo[0][0].relativenumber = false
    vim.wo[0][0].cursorline = true
    vim.wo[0][0].list = false
  end)
  M.set_keymaps(obuf)
  local row = bp and lay:row_of_bp(bp) or lay.header_rows
  pcall(vim.api.nvim_win_set_cursor, win, { row + 1, 0 })
  return obuf
end

--- Close the overview and put the sequence view back where it was.
function M.close(obuf)
  local s = M.state[obuf]
  if not s then return end
  local win = restore(s)
  M.state[obuf] = nil
  if not win then return end
  vim.api.nvim_win_call(win, function() vim.fn.winrestview(s.restore) end)
  buffer.fill(s.src, win)
end

--- Open the overview, or close the one already open on this record.
function M.toggle(buf, opts)
  if M.state[buf] then return M.close(buf) end
  local open = M.for_record(buf)
  if open and not (opts and opts.span) then return M.close(open) end
  return M.open(buf, opts)
end

---------------------------------------------------------------------------
-- Jumping
---------------------------------------------------------------------------

--- Go to the region under the cursor (or under `at`, a { row0, vcol0 } pair)
--- in the sequence view, and leave the overview.
function M.jump(obuf, at)
  local s = M.state[obuf]
  if not s then return end
  local win = (s.win and vim.api.nvim_win_is_valid(s.win)) and s.win or win_for(obuf)
  local row0, vcol0
  if at then
    row0, vcol0 = at[1], at[2]
  else
    local c = vim.api.nvim_win_get_cursor(win or 0)
    row0, vcol0 = c[1] - 1, c[2]
  end
  local bs, be = s.lay:region_of_row(row0)
  if not bs then return notify("put the cursor on a region row", vim.log.levels.WARN) end
  local bp = s.lay:bp_at(row0, vcol0) or bs
  M.state[obuf] = nil
  local w = restore(s)
  if not w then return end
  buffer.goto_bp(w, bp)
  buffer.fill(s.src, w)
  vim.api.nvim_win_call(w, function() vim.cmd("normal! zz") end)
  buffer.fill(s.src, w)
  notify(("%d..%d  (%s)"):format(bs, be, view.span_text(be - bs + 1, s.rec:unit())))
  return bp
end

--- Go to a region of the record and mark it, from anywhere.
--- With a windowed layout this lays out only the rows around `s`.
function M.region(buf, s, e)
  local st = buffer.get(buf)
  if not st then return notify("current buffer is not a sequence file", vim.log.levels.WARN) end
  local rec = buffer.record(st)
  local L = math.max(rec:len(), 1)
  if not s then return notify(":DnaRegion {start}..{stop}", vim.log.levels.ERROR) end
  s = math.max(1, math.min(math.floor(s), L))
  e = math.max(1, math.min(math.floor(e or s), L))
  local win = win_for(buf)
  if not win then return end
  local sel = require("plasmid.select")
  sel.set(buf, sel.make(rec, s, e))
  buffer.goto_bp(win, s)
  buffer.fill(buf, win)
  notify(("%d..%d  (%d %s)"):format(s, e, (e >= s and e - s + 1 or L - s + 1 + e), rec:unit()))
end

---------------------------------------------------------------------------
-- Keys
---------------------------------------------------------------------------

local function mouse_jump(obuf)
  local m = vim.fn.getmousepos()
  if not m or m.winid == 0 or vim.api.nvim_win_get_buf(m.winid) ~= obuf then return end
  local s = M.state[obuf]
  local line = s and s.lay.lines[m.line] or ""
  local vcol = vim.fn.charidx(line, math.max(0, (m.column or 1) - 1))
  if vcol < 0 then vcol = vim.fn.strchars(line) end
  M.jump(obuf, { (m.line or 1) - 1, vcol })
end

function M.set_keymaps(obuf)
  local k = (config.options.overview or {}).keymaps or {}
  local o = { buffer = obuf, nowait = true, silent = true }
  local function map(lhs, fn, desc)
    if lhs and lhs ~= "" then
      vim.keymap.set("n", lhs, fn, vim.tbl_extend("force", o, { desc = "plasmid: " .. desc }))
    end
  end
  map(k.jump or "<CR>", function() M.jump(obuf) end, "go to this region in the sequence view")
  map(k.close or "q", function() M.close(obuf) end, "close the overview")
  map(config.options.keymaps.overview, function() M.close(obuf) end, "close the overview")
  if (config.options.mouse or {}).enabled ~= false then
    map(k.click or "<LeftMouse>", function() mouse_jump(obuf) end, "go to the region under the mouse")
    map("<2-LeftMouse>", function() mouse_jump(obuf) end, "go to the region under the mouse")
  end
end

return M
