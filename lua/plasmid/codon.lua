-- Codon under the cursor: with the cursor on an amino acid of a CDS/ORF
-- translation row, highlight that residue and the three bases of its codon on
-- both strands (codons split by a line break or an intron included).
local M = {}

local ns = vim.api.nvim_create_namespace("plasmid_codon")
M.ns = ns

local function buffer() return require("plasmid.buffer") end

--- Codon at a display position of the linear layout, or nil.
--- Returns the codon_at() table plus f (feature), row0 and vcol of the residue.
function M.at(st, row0, vcol)
  local lay = st and st.layout
  if not lay or lay.view ~= "linear" then return nil end
  local info = lay.rows[row0]
  if not info or info.kind ~= "aa" then return nil end
  local line = lay.lines[row0 + 1] or ""
  local ch = vim.fn.strcharpart(line, vcol, 1)
  if ch == "" or ch == " " then return nil end
  local fid = lay:feature_at(row0, vcol)
  local rec = buffer().record(st)
  local f = fid and rec:feature_by_id(fid)
  if not f then return nil end
  local p2 = lay:bp_at(row0, vcol)
  local c = p2 and require("plasmid.view.features").codon_at(rec, f, p2)
  if not c then return nil end
  c.f, c.row0, c.vcol = f, row0, vcol
  return c
end

local function cursor_vcol(win, row0)
  local pos = vim.api.nvim_win_get_cursor(win)
  local line = vim.api.nvim_buf_get_lines(vim.api.nvim_win_get_buf(win), row0, row0 + 1, false)[1] or ""
  return vim.fn.strchars(line:sub(1, pos[2]))
end

local function mark(buf, row0, vcol, line)
  local b0 = vim.fn.byteidx(line, vcol)
  local b1 = vim.fn.byteidx(line, vcol + 1)
  if b0 < 0 then return end
  if b1 < 0 then b1 = #line end
  pcall(vim.api.nvim_buf_set_extmark, buf, ns, row0, b0, { end_col = b1, hl_group = "PlasmidCodon", priority = 200 })
end

--- Recompute and redraw the codon highlight for the window showing `buf`.
function M.update(buf, win)
  local st = buffer().get(buf)
  win = win or vim.api.nvim_get_current_win()
  if not st or not vim.api.nvim_win_is_valid(win) or vim.api.nvim_win_get_buf(win) ~= buf then return end
  local c
  if st.layout and st.layout.view == "linear" then
    local row0 = vim.api.nvim_win_get_cursor(win)[1] - 1
    c = M.at(st, row0, cursor_vcol(win, row0))
  end
  local key = c and (c.f.id .. ":" .. c.positions[2]) or nil
  if key == st.codon_key and st.codon_layout == st.layout then return end
  st.codon, st.codon_key, st.codon_layout = c, key, st.layout
  vim.api.nvim_buf_clear_namespace(buf, ns, 0, -1)
  if not c then return end
  local lay = st.layout
  mark(buf, c.row0, c.vcol, lay.lines[c.row0 + 1] or "")
  for _, p in ipairs(c.positions) do
    local r, v = lay:pos_of(p)
    if r then
      mark(buf, r, v, lay.lines[r + 1] or "")
      local below = lay.rows[r + 1]
      if below and below.kind == "comp" then mark(buf, r + 1, v, lay.lines[r + 2] or "") end
    end
  end
end

--- Winbar fragment, e.g. "codon 12/287 ATG M 889..891".
function M.describe(st)
  local c = st and st.codon
  if not c then return nil end
  local p = c.positions
  local span
  if math.abs(p[3] - p[1]) == 2 then
    local lo, hi = math.min(p[1], p[3]), math.max(p[1], p[3])
    span = c.f.strand == -1 and ("complement(%d..%d)"):format(lo, hi) or ("%d..%d"):format(lo, hi)
  else -- split by an intron or the origin
    span = ("bp %d,%d,%d"):format(p[1], p[2], p[3])
  end
  return ("codon %d/%d %s→%s %s"):format(c.n, c.total, c.codon, c.aa, span)
end

return M
