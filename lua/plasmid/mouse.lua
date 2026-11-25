-- Mouse in a sequence buffer: drag to select a region, double-click a feature
-- to select it, and Ctrl-click two primers to select the amplicon between
-- them. Everything here works in both views; on the circular map the drag
-- follows the ring.
--
-- 'mouse' must include the mode (Neovim's default covers normal mode).
local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level)
  require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO)
end
local function cfg() return require("plasmid.config").options.mouse or {} end

M.state = {}   -- per buffer: { anchor = bp, dragging = bool, primer = index }

local function state(buf)
  M.state[buf] = M.state[buf] or {}
  return M.state[buf]
end

--- Where the mouse is: buffer, window, base pair and feature id (nil outside
--- a sequence buffer or off the sequence).
function M.at_mouse()
  local m = vim.fn.getmousepos()
  if not m or m.winid == 0 then return nil end
  local buf = vim.api.nvim_win_get_buf(m.winid)
  local st = buffer().get(buf)
  if not st or not st.layout then return nil end
  local row0 = (m.line or 1) - 1
  local line = st.layout.lines[row0 + 1]
  if not line then return nil end
  local vcol = vim.fn.charidx(line, math.max(0, (m.column or 1) - 1))
  if vcol < 0 then vcol = vim.fn.strchars(line) end
  return { buf = buf, win = m.winid, row0 = row0, vcol = vcol,
    bp = st.layout:bp_at(row0, vcol), fid = st.layout:feature_at(row0, vcol), line = line }
end

local function move_cursor(p)
  local st = buffer().get(p.buf)
  local line = p.line or ""
  local col = st.layout.byte_col and st.layout:byte_col(p.row0, p.vcol) or p.vcol
  col = math.min(math.max(0, col), math.max(#line - 1, 0))
  pcall(vim.api.nvim_win_set_cursor, p.win, { p.row0 + 1, col })
  if vim.api.nvim_get_current_win() ~= p.win then pcall(vim.api.nvim_set_current_win, p.win) end
end

--- Press: move the cursor there and remember the base pair to drag from.
function M.press()
  local p = M.at_mouse()
  if not p then return end
  move_cursor(p)
  local s = state(p.buf)
  s.anchor, s.dragging = p.bp, false
  require("plasmid.select").set(p.buf, nil)
end

--- Drag: select from the pressed base pair to the one under the mouse.
function M.drag()
  local p = M.at_mouse()
  if not p or not p.bp then return end
  local s = state(p.buf)
  if not s.anchor then s.anchor = p.bp end
  move_cursor(p)
  s.dragging = true
  local rec = buffer().record(buffer().get(p.buf))
  local a, b = s.anchor, p.bp
  if a > b then a, b = b, a end
  require("plasmid.select").set(p.buf, require("plasmid.select").make(rec, a, b))
end

--- Release: report what was selected.
function M.release()
  local p = M.at_mouse()
  if not p then return end
  local s = state(p.buf)
  if not s.dragging then return end
  s.dragging = false
  local S = require("plasmid.select")
  local sel = S.get(p.buf)
  if not sel then return end
  local rec = buffer().record(buffer().get(p.buf))
  notify(("selected %d..%d (%d %s)"):format(sel.start, sel.stop, S.length(rec, sel), rec:unit()))
end

--- Double click: select the feature under the mouse.
function M.double_click()
  local p = M.at_mouse()
  if not p then return end
  move_cursor(p)
  local st = buffer().get(p.buf)
  local rec = buffer().record(st)
  local S = require("plasmid.select")
  local f = p.fid and rec:feature_by_id(p.fid)
  if not f and p.bp then f = S.feature_at_pos(rec, p.bp, require("plasmid.config").options) end
  if not f then return notify("no feature there", vim.log.levels.WARN) end
  S.set(p.buf, S.for_feature(rec, f))
  notify(("selected %s (%d %s)"):format(f:label(), f:length(), rec:unit()))
end

-- The primer under the mouse, as an entry of primers.list.
local function primer_at(rec, p)
  local list = require("plasmid.primers").list(rec)
  if p.fid then
    for _, x in ipairs(list) do if x.f.id == p.fid then return x end end
  end
  if not p.bp then return nil end
  for _, x in ipairs(list) do
    local a, b = x.start, x.stop
    local inside
    if a <= b then inside = p.bp >= a and p.bp <= b     -- plain span
    else inside = p.bp >= a or p.bp <= b end            -- across the origin
    if inside then return x end
  end
  return nil
end

--- Ctrl-click: the first primer is remembered, the second selects the product
--- the pair makes (SnapGene's amplicon highlight).
function M.ctrl_click()
  local p = M.at_mouse()
  if not p then return end
  move_cursor(p)
  local st = buffer().get(p.buf)
  local rec = buffer().record(st)
  local s = state(p.buf)
  local hit = primer_at(rec, p)
  if not hit then
    s.primer = nil
    return notify("no primer there (Ctrl-click a primer, then its partner)", vim.log.levels.WARN)
  end
  if not s.primer or s.primer.f.id == hit.f.id or not rec:feature_by_id(s.primer.f.id) then
    s.primer = hit
    return notify(("%s picked — Ctrl-click the other primer to show the product"):format(hit.name))
  end
  local first, second = s.primer, hit
  s.primer = nil
  local chosen, ta = M.amplicon(rec, first, second)
  if not chosen then
    return notify(("%s and %s make no product (they face away from each other, or do not both bind)")
      :format(first.name, second.name), vim.log.levels.WARN)
  end
  local S = require("plasmid.select")
  S.set(p.buf, S.make(rec, chosen.start, chosen.stop))
  local tails = #chosen.fwd.tail + #chosen.rev.tail
  notify(("%s / %s: %d %s%s (%d..%d)%s"):format(first.name, second.name, chosen.region, rec:unit(),
    tails > 0 and (" of template, %d bp product with the tails"):format(chosen.length) or " product",
    chosen.start, chosen.stop, ta and (", Ta %d °C"):format(ta) or ""))
end

--- The product two primers make on this record, and the annealing temperature
--- for the first cycle. Only the annealed part of a primer binds the template,
--- so a 5' tail is left out of the Ta: a 20 nt Gibson overhang would otherwise
--- raise it by ten degrees and send the reaction in too hot.
---@return table|nil product, integer|nil ta
function M.amplicon(rec, first, second)
  local PCR = require("plasmid.pcr")
  local prods = PCR.products(rec.seq:tostring(), rec:is_circular(), first.seq, second.seq)
  local fwd = first.strand == 1 and first or second
  local rev = first.strand == 1 and second or first
  local chosen
  for _, pr in ipairs(prods) do
    if pr.start == fwd.start and pr.stop == rev.stop then chosen = pr break end
  end
  chosen = chosen or prods[1]
  if not chosen then return nil end
  local tm = require("plasmid.bio.tm")
  local cond = require("plasmid.primerinfo").conditions()
  local fa = chosen.fwd_primer:sub(#chosen.fwd.tail + 1)
  local ra = chosen.rev_primer:sub(#chosen.rev.tail + 1)
  local ok, ta = pcall(function() return (select(1, tm.ta(fa, ra, cond))) end)
  return chosen, ok and ta or nil
end

--- Buffer-local mouse mappings, following config.mouse.
function M.attach(buf)
  local c = cfg()
  if c.enabled == false then return end
  local function map(lhs, fn, desc)
    vim.keymap.set("n", lhs, fn, { buffer = buf, silent = true, nowait = true,
      desc = "plasmid: " .. desc })
  end
  if c.drag_select ~= false then
    map("<LeftMouse>", M.press, "move the cursor, start a selection")
    map("<LeftDrag>", M.drag, "drag to select a region")
    map("<LeftRelease>", M.release, "finish the selection")
  end
  if c.double_click_feature ~= false then
    map("<2-LeftMouse>", M.double_click, "select the feature under the mouse")
  end
  if c.primer_pair ~= false then
    map("<C-LeftMouse>", M.ctrl_click, "pick two primers to show their product")
  end
end

return M
