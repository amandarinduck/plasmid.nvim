-- Selection: a base-pair range on the model, independent of the rendering.
--
--   sel = { start = pos1, stop = pos1, fid = feature id|nil }
--   On a circular record start > stop means the range wraps the origin.
--
-- The projection buffer shows a rendering, so yanking its text would copy
-- ruler numbers and track glyphs. Everything here yanks residues instead.
local alphabet = require("plasmid.model.alphabet")
local location = require("plasmid.model.location")

local M = {}

local function buffer() return require("plasmid.buffer") end

---------------------------------------------------------------------------
-- Pure helpers (no Neovim state)
---------------------------------------------------------------------------

--- Normalise a range on a record. On linear records start/stop are ordered;
--- on circular records start > stop is kept and means "wraps the origin".
function M.make(rec, a, b, fid)
  local L = rec:len()
  a = math.max(1, math.min(math.floor(a), L))
  b = math.max(1, math.min(math.floor(b), L))
  if a > b and not rec:is_circular() then a, b = b, a end
  return { start = a, stop = b, fid = fid }
end

function M.wraps(sel)
  return sel.start > sel.stop
end

--- Linear intervals covered by the selection.
function M.ranges(rec, sel)
  if M.wraps(sel) then return { { sel.start, rec:len() }, { 1, sel.stop } } end
  return { { sel.start, sel.stop } }
end

function M.length(rec, sel)
  local n = 0
  for _, r in ipairs(M.ranges(rec, sel)) do n = n + r[2] - r[1] + 1 end
  return n
end

function M.contains(rec, sel, pos)
  for _, r in ipairs(M.ranges(rec, sel)) do
    if pos >= r[1] and pos <= r[2] then return true end
  end
  return false
end

--- Top-strand residues of the selection, 5'->3' left to right (wrapping).
function M.top_strand(rec, sel)
  local out = {}
  for _, r in ipairs(M.ranges(rec, sel)) do out[#out + 1] = rec.seq:sub(r[1], r[2]) end
  return table.concat(out)
end

--- Selection covering a feature's whole span (introns included).
function M.for_feature(rec, f)
  local ext = rec:extents(f)
  if #ext == 2 then
    local head, tail
    for _, e in ipairs(ext) do if e[1] == 1 then head = e else tail = e end end
    return { start = tail[1], stop = head[2], fid = f.id }
  end
  return { start = ext[1][1], stop = ext[1][2], fid = f.id }
end

--- Residues for a feature: "a" = whole span, "i" = its own (spliced) parts.
--- Both are returned in the feature's 5'->3' orientation.
function M.feature_residues(rec, f, kind)
  if kind == "i" then return rec:extract(f) end
  local s = M.top_strand(rec, M.for_feature(rec, f))
  if f.strand == -1 and not rec:is_protein() then s = alphabet.revcomp(s) end
  return s
end

--- Innermost visible feature containing pos (smallest span wins).
function M.feature_at_pos(rec, pos, opts)
  local vf = require("plasmid.view.features")
  local best, best_len
  for _, it in ipairs(vf.visible(rec, opts)) do
    for _, e in ipairs(it.extents) do
      if pos >= e[1] and pos <= e[2] then
        local len = location.length(it.f.parts)
        if not best or len < best_len then best, best_len = it.f, len end
      end
    end
  end
  return best
end

---------------------------------------------------------------------------
-- Buffer-bound actions
---------------------------------------------------------------------------

local function redraw(buf) require("plasmid.util").repaint(buf) end

--- Set (or clear with nil) the selection of a buffer and repaint.
function M.set(buf, sel)
  local st = buffer().get(buf)
  if not st then return end
  st.selection = sel
  if st.view == "circular" then
    buffer().render(buf)
  else
    if st.layout then st.layout.selection = sel end
    redraw(buf)
  end
end

function M.get(buf)
  local st = buffer().get(buf)
  return st and st.selection
end

--- Feature for a text object: the one drawn under the cursor, else the
--- innermost feature containing the cursor's base.
function M.feature_under_cursor(win)
  local buf = vim.api.nvim_win_get_buf(win)
  local st = buffer().get(buf)
  if not st then return nil end
  local rec = buffer().record(st)
  local bp, fid = buffer().cursor_info(win)
  if fid then return rec:feature_by_id(fid) end
  if bp then return M.feature_at_pos(rec, bp, require("plasmid.config").options) end
end

---------------------------------------------------------------------------
-- Registers and the system clipboard
---------------------------------------------------------------------------

--- Write `text` to clipboard register "+" or "*". Uses Neovim's clipboard
--- provider when there is one (xclip, wl-copy, pbcopy, win32yank, g:clipboard,
--- ...); otherwise, with a terminal attached, the OSC 52 escape sequence, which
--- most terminals (and tmux with set-clipboard on) turn into a clipboard copy,
--- also over SSH. Returns "provider", "osc52" or nil. Replaceable for tests.
function M.clipboard_write(reg, text)
  if vim.fn.has("clipboard") == 1 then
    vim.fn.setreg(reg, text, "c")
    return "provider"
  end
  if #vim.api.nvim_list_uis() > 0 then
    local ok = pcall(function() require("vim.ui.clipboard.osc52").copy(reg)({ text }) end)
    if ok then return "osc52" end
  end
  return nil
end

--- Which clipboard register a plain yank (no register given) also goes to,
--- per the `yank_clipboard` option: true -> "+", "auto" -> follow 'clipboard',
--- false -> none.
function M.clipboard_target()
  local opt = require("plasmid.config").options.yank_clipboard
  if opt == false then return nil end
  if opt == "auto" then
    local cb = vim.o.clipboard
    if cb:find("unnamedplus", 1, true) then return "+" end
    if cb:find("unnamed", 1, true) then return "*" end
    return nil
  end
  if opt == "*" then return "*" end
  return "+"
end

local warned_no_clipboard = false

--- Put `text` in a register. A plain yank (reg nil, "" or the unnamed
--- register) goes to the unnamed register and, per `yank_clipboard`, the system
--- clipboard; a named register goes there (and to the unnamed register).
--- Residues as the linear view shows them (`linear.case`): "auto" upper-cases
--- GenBank/SnapGene records and keeps FASTA case (soft-masking); "upper";
--- "preserve" keeps the file's case.
function M.display_case(rec, text)
  local case = require("plasmid.config").options.linear.case
  if case == "auto" then case = (rec and rec.format == "fasta") and "preserve" or "upper" end
  return case == "upper" and text:upper() or text
end

local function yank(reg, text, what)
  local st = buffer().get(vim.api.nvim_get_current_buf())
  if st then text = M.display_case(buffer().record(st), text) end
  local explicit = reg ~= nil and reg ~= "" and reg ~= '"'
  -- with 'clipboard' set, Neovim reports "+"/"*" as v:register for plain yanks
  if explicit and (reg == "+" or reg == "*") and vim.o.clipboard:find("unnamed", 1, true) then
    explicit = false
  end
  local clip
  if not explicit then
    clip = M.clipboard_target()
  elseif reg == "+" or reg == "*" then
    clip = reg
  end -- other named registers: Neovim only
  if explicit and reg ~= "+" and reg ~= "*" then vim.fn.setreg(reg, text, "c") end
  vim.fn.setreg('"', text, "c")
  -- The clipboard write (which may start a provider process) and the message
  -- (which may open windows) are deferred: yif/yaf run under textlock.
  vim.schedule(function()
    local via = clip and M.clipboard_write(clip, text)
    local dest = explicit and reg ~= "+" and reg ~= "*" and ("register " .. reg) or 'register "'
    if via then
      dest = dest .. (clip == "+" and " and the clipboard" or " and the selection clipboard (*)")
    end
    local level, tail = vim.log.levels.INFO, ""
    if clip and not via then
      level = (explicit or not warned_no_clipboard) and vim.log.levels.WARN or vim.log.levels.INFO
      tail = "; no clipboard available (:checkhealth plasmid)"
      warned_no_clipboard = true
    end
    require("plasmid.util").notify(("plasmid: yanked %s (%d %s) to %s%s"):format(what, #text,
      #text == 1 and "residue" or "residues", dest, tail), level)
  end)
end

function M._reset_clipboard_warning() warned_no_clipboard = false end

--- Operator-pending text object. Returns keys for an <expr> mapping:
--- for `y` it yanks residues itself and cancels the operator.
function M.textobj_pending(kind)
  local op = vim.v.operator
  local reg = vim.v.register
  local win = vim.api.nvim_get_current_win()
  local f = M.feature_under_cursor(win)
  if not f then
    vim.schedule(function() require("plasmid.util").notify("plasmid: no feature here", vim.log.levels.WARN) end)
    return "<Esc>"
  end
  if op == "d" or op == "c" then
    -- the edit cannot run under textlock: do it right after the mapping
    local buf = vim.api.nvim_get_current_buf()
    vim.schedule(function() require("plasmid.seqedit").feature_op(buf, f, kind, op) end)
    return "<Esc>"
  end
  if op ~= "y" then
    vim.schedule(function()
      require("plasmid.util").notify(("plasmid: '%s' on features is not supported (y, d, c)"):format(op), vim.log.levels.WARN)
    end)
    return "<Esc>"
  end
  local rec = buffer().record(buffer().get(0))
  local text = M.feature_residues(rec, f, kind)
  -- registers may be written under textlock, so the yank is immediate
  yank(reg, text, (kind == "a" and "span of " or "") .. f:label())
  return "<Esc>"
end

--- Visual-mode text object: select the feature and leave visual mode.
function M.textobj_visual(kind)
  local buf = vim.api.nvim_get_current_buf()
  local f = M.feature_under_cursor(vim.api.nvim_get_current_win())
  vim.cmd("normal! \27")
  if not f then return require("plasmid.util").notify("plasmid: no feature here", vim.log.levels.WARN) end
  local rec = buffer().record(buffer().get(buf))
  local sel = M.for_feature(rec, f)
  sel.inner = kind == "i"
  M.set(buf, sel)
end

--- bp range of the current/last visual selection in the projection buffer.
function M.visual_range(buf)
  local st = buffer().get(buf)
  local lay = st.layout
  local p1, p2 = vim.fn.getpos("v"), vim.fn.getpos(".")
  if vim.fn.mode():match("^[vV\22]") == nil then
    p1, p2 = vim.fn.getpos("'<"), vim.fn.getpos("'>")
  end
  local function bp_of(p)
    local row0 = p[2] - 1
    local line = lay.lines[row0 + 1] or ""
    local vcol = vim.fn.charidx(line, math.max(p[3] - 1, 0))
    return lay:bp_at(row0, vcol)
  end
  local a, b = bp_of(p1), bp_of(p2)
  if not a or not b then return nil end
  if a > b then a, b = b, a end -- visual ranges in the linear view never wrap
  return a, b
end

--- Visual y / Y: yank the bases under the visual selection.
function M.yank_visual(revcomp)
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer().get(buf)
  local reg = vim.v.register
  local a, b = M.visual_range(buf)
  vim.cmd("normal! \27")
  if not a then return require("plasmid.util").notify("plasmid: selection is not on the sequence", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  local sel = M.make(rec, a, b)
  M.set(buf, sel)
  local s = M.top_strand(rec, sel)
  if revcomp then s = alphabet.revcomp(s) end
  yank(reg, s, ("%d..%d%s"):format(a, b, revcomp and " (reverse complement)" or ""))
end

--- gy / gY: yank the model selection.
function M.yank_selection(revcomp)
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer().get(buf)
  local sel = st and st.selection
  if not sel then return require("plasmid.util").notify("plasmid: nothing selected", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  local s
  local f = sel.fid and rec:feature_by_id(sel.fid)
  if f and sel.inner then
    s = rec:extract(f)
    if revcomp then s = alphabet.revcomp(s) end
  else
    s = M.top_strand(rec, sel)
    if revcomp then s = alphabet.revcomp(s) end
  end
  yank(vim.v.register, s, ("%d..%d%s"):format(sel.start, sel.stop, revcomp and " (reverse complement)" or ""))
end

---------------------------------------------------------------------------
-- Strand-aware yanks
---------------------------------------------------------------------------

M.STRANDS = {
  { key = "top", label = "top strand 5'->3'" },
  { key = "top35", label = "top strand 3'->5' (reversed)" },
  { key = "bottom", label = "bottom strand 5'->3' (reverse complement)" },
  { key = "bottom35", label = "bottom strand 3'->5' (complement, as displayed)" },
}
local STRAND_ALIASES = {
  top = "top", top53 = "top", forward = "top", ["+"] = "top",
  top35 = "top35", reverse = "top35", rev = "top35",
  bottom = "bottom", bottom53 = "bottom", revcomp = "bottom", rc = "bottom", ["-"] = "bottom",
  bottom35 = "bottom35", complement = "bottom35", comp = "bottom35",
}

function M.strand_key(name) return STRAND_ALIASES[(name or ""):lower()] end

--- `plus` is top-strand text read left to right; return it in orientation `key`.
function M.orient(plus, key)
  if key == "top35" then return plus:reverse() end
  if key == "bottom" then return alphabet.revcomp(plus) end
  if key == "bottom35" then return alphabet.complement(plus) end
  return plus
end

local function label_of(key)
  for _, s in ipairs(M.STRANDS) do if s.key == key then return s.label end end
end

-- Top-strand text of a selection. An inner (spliced) feature selection gives
-- the spliced residues in top-strand orientation; also returns the feature's
-- own strand key so the picker can mark the coding strand.
local function plus_text(rec, sel)
  local f = sel.fid and rec:feature_by_id(sel.fid)
  if f and sel.inner then
    local s = rec:extract(f)
    if f.strand == -1 then return alphabet.revcomp(s), "bottom" end
    return s, "top"
  end
  return M.top_strand(rec, sel), f and (f.strand == -1 and "bottom" or "top") or nil
end

--- Yank `sel` in strand orientation `key`; nil key opens a picker.
function M.yank_strand(buf, sel, key, reg)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  if rec:is_protein() then
    return require("plasmid.util").notify("plasmid: protein records have one strand", vim.log.levels.WARN)
  end
  reg = reg or '"'
  local plus, sense = plus_text(rec, sel)
  local range = ("%d..%d"):format(sel.start, sel.stop)
  local function go(k)
    yank(reg, M.orient(plus, k), ("%s %s"):format(range, label_of(k)))
  end
  if key then return go(key) end
  vim.ui.select(M.STRANDS, {
    prompt = ("Yank %s (%d bp) as"):format(range, #plus),
    format_item = function(s)
      local prev = M.orient(plus, s.key)
      if #prev > 24 then prev = prev:sub(1, 24) .. "…" end
      return ("%-46s %s%s"):format(s.label, prev, s.key == sense and "   (feature's strand)" or "")
    end,
  }, function(choice) if choice then go(choice.key) end end)
end

--- Visual gy: pick the strand for the bases under the visual selection.
function M.yank_visual_pick()
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer().get(buf)
  local reg = vim.v.register
  local a, b = M.visual_range(buf)
  vim.cmd("normal! \27")
  if not a then return require("plasmid.util").notify("plasmid: selection is not on the sequence", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  local sel = M.make(rec, a, b)
  M.set(buf, sel)
  M.yank_strand(buf, sel, nil, reg)
end

--- Winbar fragment.
function M.describe(rec, sel)
  if not sel then return nil end
  return ("sel %d..%d (%d %s)"):format(sel.start, sel.stop, M.length(rec, sel), rec:unit())
end

return M
