-- Buffer integration: one projection buffer per file.
--
-- The Record model is the source of truth. The buffer is a read-only
-- rendering of it (linear or circular); `nomodifiable` is enforced so the
-- view cannot drift from the model. Writing serialises the model.
local config = require("plasmid.config")
local pio = require("plasmid.io")
local linear = require("plasmid.view.linear")
local circular = require("plasmid.view.circular")

local M = {}

---@type table<integer, table>
M.state = {}

local ns = vim.api.nvim_create_namespace("plasmid")
M.ns = ns

local function notify(msg, level)
  require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO)
end

function M.get(buf)
  return M.state[buf == 0 and vim.api.nvim_get_current_buf() or buf]
end

function M.record(st)
  return st.records[st.current]
end

local function win_for(buf)
  local cur = vim.api.nvim_get_current_win()
  if vim.api.nvim_win_get_buf(cur) == buf then return cur end
  local wins = vim.fn.win_findbuf(buf)
  return wins[1]
end

local function text_width(win)
  local info = vim.fn.getwininfo(win)[1]
  return info and (info.width - info.textoff) or vim.api.nvim_win_get_width(win)
end

---------------------------------------------------------------------------
-- Cursor <-> base pair
---------------------------------------------------------------------------

--- Display column (0-based) of the cursor in `win`.
local function cursor_vcol(win, line, bytecol)
  return vim.fn.charidx(line, bytecol)
end

--- Base pair and feature id under the cursor of `win`.
function M.cursor_info(win)
  local buf = vim.api.nvim_win_get_buf(win)
  local st = M.state[buf]
  if not st or not st.layout then return nil end
  local pos = vim.api.nvim_win_get_cursor(win)
  local row0 = pos[1] - 1
  if st.layout.windowed and (st.layout.lines[row0 + 1] or "") == "" then M.fill(buf, win, row0) end
  local line = st.layout.lines[row0 + 1] or ""
  local vcol = cursor_vcol(win, line, pos[2])
  local bp = st.layout:bp_at(row0, vcol)
  local fid = st.layout:feature_at(row0, vcol)
  return bp, fid
end

--- Move the cursor of `win` to a base pair in the current view.
function M.goto_bp(win, bp)
  local st = M.state[vim.api.nvim_win_get_buf(win)]
  if not st or not st.layout then return end
  local rec = M.record(st)
  bp = math.max(1, math.min(math.floor(bp), math.max(rec:len(), 1)))
  local row0, vcol = st.layout:pos_of(bp)
  if not row0 then return end
  local buf = vim.api.nvim_win_get_buf(win)
  if st.layout.windowed then M.fill(buf, win, row0) end
  local line = st.layout.lines[row0 + 1] or ""
  local bytecol = st.layout.byte_col and st.layout:byte_col(row0, vcol) or vcol
  bytecol = math.min(bytecol, math.max(#line - 1, 0))
  vim.api.nvim_win_set_cursor(win, { row0 + 1, bytecol })
  st.last_bp = bp
end

---------------------------------------------------------------------------
-- Windowed layouts: fill the rows a window is about to show
---------------------------------------------------------------------------

--- Write a windowed layout's newly built (or dropped) rows into the buffer.
--- The buffer stays nomodifiable to the user and keeps its modified flag;
--- these edits are the view catching up with itself, not a change to the
--- record, so they are kept out of the undo history.
local function apply_edits(buf, edits)
  if #edits == 0 then return end
  local bo = vim.bo[buf]
  local modified, undo = bo.modified, bo.undolevels
  bo.modifiable = true
  bo.undolevels = -1
  for _, e in ipairs(edits) do
    pcall(vim.api.nvim_buf_set_lines, buf, e[1], e[1] + #e[2], false, e[2])
  end
  bo.undolevels = undo
  bo.modifiable = false
  bo.modified = modified
end

--- Build the rows `win` is showing (plus a margin) if the layout is windowed.
---@param around integer|nil extra row to cover (a jump target, say)
function M.fill(buf, win, around)
  local st = M.state[buf]
  local lay = st and st.layout
  if not lay or not lay.windowed then return end
  win = win or win_for(buf)
  local first, last
  if win and vim.api.nvim_win_is_valid(win) then
    local info = vim.fn.getwininfo(win)[1]
    first, last = (info.topline or 1) - 1, (info.botline or 1) - 1
  else
    first, last = 0, 40
  end
  if around then
    -- a jump lands somewhere else entirely: cover where it lands, not
    -- everything between here and there
    if around < first or around > last then
      local h = math.max(10, last - first + 1)
      first, last = around - h, around + h
    end
  end
  apply_edits(buf, lay:ensure(first, last))
end

---------------------------------------------------------------------------
-- Rendering
---------------------------------------------------------------------------

local function set_lines(buf, lines)
  vim.bo[buf].modifiable = true
  vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
  vim.bo[buf].modifiable = false
  vim.bo[buf].modified = false
end

--- (Re)build the layout for the buffer's current view and window size.
---@param opts { keep_bp?: integer }|nil
function M.render(buf, opts)
  opts = opts or {}
  local st = M.state[buf]
  if not st then return end
  local win = win_for(buf)
  local width = win and text_width(win) or vim.o.columns
  local height = win and vim.fn.winheight(win) or (vim.o.lines - 2)  -- text rows (winbar excluded)
  local keep = opts.keep_bp
  if not keep and win and st.layout then keep = (M.cursor_info(win)) end
  local rec = M.record(st)
  local o = config.options
  if st.track_limit then
    o = vim.tbl_deep_extend("force", o, {
      linear = { max_tracks = st.track_limit }, circular = { max_tracks = st.track_limit },
    })
  end
  if st.group then o = vim.tbl_deep_extend("force", o, { linear = { group = st.group } }) end
  if st.view == "circular" then
    st.layout = circular.build(rec, o, width, height, { selection = st.selection, search = st.search,
      enzymes = st.enzymes and st.enzymes.entries })
  else
    local primers
    if (o.primers or {}).show ~= false and not rec:is_protein() then
      primers = require("plasmid.primers").list(rec)
      if #primers == 0 then primers = nil end
    end
    st.layout = linear.build(rec, o, width, { enzymes = st.enzymes and st.enzymes.entries, primers = primers,
      primer_cells = function(p) return require("plasmid.primers").cells(rec, p) end })
    st.layout.selection = st.selection
    st.layout.search = st.search
  end
  set_lines(buf, st.layout.lines)
  st.modified = st.modified or false
  vim.bo[buf].modified = st.modified
  M.fill(buf, win)
  if win and keep then M.goto_bp(win, keep) end
  M.fill(buf, win)
  if win then require("plasmid.codon").update(buf, win) end
end

--- Decoration provider: highlight only rows that are being drawn.
function M.setup_decorations()
  vim.api.nvim_set_decoration_provider(ns, {
    on_win = function(_, _, buf)
      return M.state[buf] ~= nil and M.state[buf].layout ~= nil
    end,
    on_line = function(_, _, buf, row)
      local st = M.state[buf]
      local lay = st and st.layout
      if not lay then return end
      for _, h in ipairs(lay:row_highlights(row)) do
        local b0, b1 = h[1], h[2]
        local line_len = #(lay.lines[row + 1] or "")
        if b1 == -1 or b1 > line_len then b1 = line_len end
        if b0 < b1 then
          pcall(vim.api.nvim_buf_set_extmark, buf, ns, row, b0, {
            end_col = b1, hl_group = h[3], ephemeral = true, priority = h[4] or 100,
          })
        end
      end
    end,
  })
end

---------------------------------------------------------------------------
-- Views
---------------------------------------------------------------------------

function M.set_view(buf, view, force)
  local st = M.state[buf]
  if not st then return notify("not a plasmid buffer", vim.log.levels.WARN) end
  local rec = M.record(st)
  if view == "circular" and rec:is_protein() then
    return notify(("%s is a protein sequence; only the linear view applies"):format(rec.name), vim.log.levels.WARN)
  end
  if view == "circular" and not rec:is_circular() and not force then
    return notify(("%s is linear; use :DnaView! circular to draw it anyway"):format(rec.name),
      vim.log.levels.WARN)
  end
  local win = win_for(buf)
  local bp = win and (M.cursor_info(win)) or st.last_bp
  st.view = view
  M.render(buf, { keep_bp = bp })
  -- The window keeps the old view's scroll position; when the new layout fits
  -- in the window show all of it, otherwise centre the cursor.
  if win and vim.api.nvim_win_is_valid(win) then
    vim.api.nvim_win_call(win, function()
      if vim.api.nvim_buf_line_count(buf) <= vim.fn.winheight(win) then
        vim.fn.winrestview({ topline = 1 })
      else
        vim.cmd("normal! zz")
      end
    end)
  end
end

--- Override the track limit for this buffer: integer, "all", or nil (config).
function M.set_track_limit(buf, limit)
  local st = M.state[buf]
  if not st then return end
  st.track_limit = limit
  M.render(buf)
end

--- Residues per group in the linear view for this buffer (0 = continuous);
--- nil returns to the configured value.
function M.set_group(buf, n)
  local st = M.state[buf]
  if not st then return end
  st.group = n
  local win = win_for(buf)
  local bp = win and (M.cursor_info(win)) or st.last_bp
  M.render(buf, { keep_bp = bp })
end

function M.group_of(st)
  return st.group or tonumber(config.options.linear.group) or 10
end

function M.toggle_view(buf)
  local st = M.state[buf]
  if not st then return end
  M.set_view(buf, st.view == "linear" and "circular" or "linear", true)
end

---------------------------------------------------------------------------
-- Opening / writing
---------------------------------------------------------------------------

local function initial_view(rec)
  local v = config.options.default_view
  if rec:is_protein() then return "linear" end
  if v == "auto" then return rec:is_circular() and "circular" or "linear" end
  if v == "circular" and not rec:is_circular() then return "linear" end
  return v
end

local function setup_buffer(buf)
  local bo = vim.bo[buf]
  bo.buftype = "acwrite"
  bo.swapfile = false
  bo.undolevels = -1
  bo.modifiable = false
  bo.filetype = "plasmid"
  require("plasmid.ui").attach_keymaps(buf)
end

--- Load a sequence file into `buf` (BufReadCmd handler).
function M.open(buf, path)
  local res, err = pio.read(path)
  if not res then
    if vim.fn.filereadable(path) == 0 then
      -- new file: an empty linear record, written in the extension's format
      res = {
        records = { require("plasmid.model.record").Record.new({
          name = vim.fn.fnamemodify(path, ":t:r"), seq = "", format = pio.format_for_path(path) or "genbank",
        }) },
        file = {}, format = pio.format_for_path(path) or "genbank",
      }
    else
      notify(("could not read %s: %s"):format(path, tostring(err)), vim.log.levels.ERROR)
      return false
    end
  end
  for _, r in ipairs(res.records) do
    if r.no_sequence and r:len() == 0 then
      local what = r.wgs and ("WGS master: %s contig records (%s), no sequence"):format(
        r.record_count and tostring(r.record_count) or "?", r.wgs)
        or ("CON record: %s bp assembled from other records, no sequence in this file"):format(
          tostring(r.declared_length or "?"))
      notify(("%s: %s; fetch it with :DnaFetch %s to get the residues"):format(r.name, what, r.id or r.name),
        vim.log.levels.WARN)
    elseif r.declared_length and r.declared_length ~= r:len() then
      notify(("%s: LOCUS declares %d bp but sequence has %d"):format(r.name, r.declared_length, r:len()),
        vim.log.levels.WARN)
    end
  end
  M.state[buf] = {
    path = path, format = res.format, file = res.file,
    records = res.records, current = 1, modified = false,
  }
  M.state[buf].view = initial_view(res.records[1])
  M.mark_saved(buf)
  setup_buffer(buf)
  M.render(buf, { keep_bp = 1 })
  local win = win_for(buf)
  if win then vim.api.nvim_win_set_cursor(win, { 1, 0 }) ; M.goto_bp(win, 1) end
  return true
end

--- Write the model to `path` (BufWriteCmd handler). Converts format when the
--- target extension implies a different one (`:w copy.fasta` from GenBank).
function M.write(buf, path)
  local st = M.state[buf]
  if not st then return false end
  local target_fmt = pio.format_for_path(path) or st.format
  if target_fmt == "fasta" then
    for _, r in ipairs(st.records) do
      if r.raw and r.raw.features_dirty and #r.features > 0 then
        notify("FASTA cannot store features; your feature changes are not in this file (save as .gb to keep them)",
          vim.log.levels.WARN)
        break
      end
    end
  end
  local ok, err = pio.write(path, st.records, target_fmt, st.file, st.format)
  if not ok then
    notify(("write failed: %s"):format(tostring(err)), vim.log.levels.ERROR)
    return false
  end
  local same = vim.fn.fnamemodify(path, ":p") == vim.fn.fnamemodify(st.path, ":p")
  local renamed = vim.fn.fnamemodify(vim.api.nvim_buf_get_name(buf), ":p") == vim.fn.fnamemodify(path, ":p")
  if same or renamed then
    st.path = path
    if target_fmt ~= st.format then
      -- the buffer now represents the new file; reload so layout metadata matches
      M.open(buf, path)
    else
      M.mark_saved(buf)
    end
  end
  local recs = #st.records
  if target_fmt == "snapgene" then target_fmt = "SnapGene" end
  notify(("wrote %s (%s, %d record%s)"):format(vim.fn.fnamemodify(path, ":~:."), target_fmt, recs,
    recs == 1 and "" or "s"))
  return true
end

--- Remember the current state of every record as "what is on disk".
function M.mark_saved(buf)
  local st = M.state[buf]
  local tx = require("plasmid.model.transaction")
  st.saved = {}
  for i, r in ipairs(st.records) do st.saved[i] = tx.state_id(r) end
  st.modified = false
  vim.bo[buf].modified = false
end

--- Modified = any record differs from its saved state (undo can clear it).
function M.refresh_modified(buf)
  local st = M.state[buf]
  local tx = require("plasmid.model.transaction")
  local mod = false
  for i, r in ipairs(st.records) do
    if tx.state_id(r) ~= (st.saved and st.saved[i] or 0) then mod = true end
  end
  st.modified = mod
  vim.bo[buf].modified = mod
end

function M.mark_modified(buf)
  M.refresh_modified(buf)
end

function M.detach(buf)
  M.state[buf] = nil
end

---------------------------------------------------------------------------
-- Winbar
---------------------------------------------------------------------------

function M.winbar()
  local win = vim.g.statusline_winid or vim.api.nvim_get_current_win()
  local buf = vim.api.nvim_win_get_buf(win)
  local st = M.state[buf]
  if not st or not st.layout then return "" end
  local rec = M.record(st)
  local parts = {
    ("%%#PlasmidTitle# %s%%*"):format(rec.name:gsub("%%", "%%%%")),
    rec:is_protein() and ("%d aa protein"):format(rec:len()) or ("%d %s %s"):format(rec:len(), rec:unit(), rec.topology),
  }
  if #st.records > 1 then parts[#parts + 1] = ("record %d/%d"):format(st.current, #st.records) end
  parts[#parts + 1] = st.view
  if st.modified then parts[#parts + 1] = "%#WarningMsg#unsaved%*" end
  local ok, bp, fid = pcall(M.cursor_info, win)
  if ok and bp then parts[#parts + 1] = ("%s %d"):format(rec:unit(), bp) end
  local seltxt = require("plasmid.select").describe(rec, st.selection)
  if seltxt then parts[#parts + 1] = "%#PlasmidSelectionArc#" .. seltxt .. "%*" end
  local rec_ = M.record(st)
  if rec_ and not rec_:is_protein() and #require("plasmid.primers").list(rec_) > 0 then
    local pol = require("plasmid.primerinfo").conditions().polymerase
    parts[#parts + 1] = "%#PlasmidPrimer#" .. pol .. "%*"
  end
  local enztxt = require("plasmid.enzymes").describe(st)
  if enztxt then parts[#parts + 1] = "%#PlasmidEnzymeUnique#" .. enztxt .. "%*" end
  local searchtxt = require("plasmid.searchui").describe(st)
  if searchtxt then parts[#parts + 1] = "%#PlasmidSearchInfo#" .. searchtxt:gsub("%%", "%%%%") .. "%*" end
  local codontxt = require("plasmid.codon").describe(st)
  if codontxt then parts[#parts + 1] = "%#PlasmidCodonInfo#" .. codontxt .. "%*" end
  if ok and fid then
    local f = rec:feature_by_id(fid)
    if f then
      local lo, hi = f:bounds()
      parts[#parts + 1] = ("%%#%s#%s %s%%* %d..%d"):format(
        require("plasmid.view.highlights").for_type(f.type), f.type,
        f:label(config.options.label_qualifiers):gsub("%%", "%%%%"), lo or 0, hi or 0)
    end
  end
  return table.concat(parts, "  │  ")
end

return M
