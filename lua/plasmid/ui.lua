-- Interactive pieces: keymaps, feature info float, feature picker, navigation.
local config = require("plasmid.config")
local buffer = require("plasmid.buffer")

local M = {}

local function state_and_record(buf)
  local st = buffer.get(buf)
  if not st then return nil end
  return st, buffer.record(st)
end

local function feature_lines(rec, f)
  local lo, hi = f:bounds()
  local strand = ({ [1] = "+ (forward)", [-1] = "− (reverse)", [0] = "mixed" })[f.strand] or "?"
  if rec:is_protein() or f.no_arrow then strand = rec:is_protein() and "n/a" or "none (non-directional)" end
  local lines = {
    ("%s  %s"):format(f.type, f:label(config.options.label_qualifiers)),
    "",
    "location  " .. f.location_text,
    ("span      %s..%s   %d %s   strand %s"):format(lo or "?", hi or "?", f:length(), rec:unit(), strand),
  }
  for _, l in ipairs(require("plasmid.fetch").link_lines(f)) do lines[#lines + 1] = l end
  if f.type == "CDS" and not rec:is_protein() then
    local codons = math.floor(math.max(f:length() - ((tonumber(f:get("codon_start")) or 1) - 1), 0) / 3)
    local tr = f:get("translation")
    if tr and tr ~= "" then
      -- the annotated protein; the CDS span usually also covers the stop codon
      lines[#lines + 1] = ("protein   %d aa   (%d codons)"):format(#tr:gsub("%s", ""), codons)
    else
      lines[#lines + 1] = ("protein   %d codons"):format(codons)
    end
  end
  lines[#lines + 1] = ""
  for _, q in ipairs(f.qualifiers) do
    local v = q.value == true and "" or tostring(q.value)
    if q.key == "translation" then
      lines[#lines + 1] = "/translation"
      for k = 1, #v, 60 do lines[#lines + 1] = "  " .. v:sub(k, k + 59) end
    else
      lines[#lines + 1] = ("/%s %s"):format(q.key, v)
    end
  end
  if f.type == "CDS" and not f:get("translation") and not rec:is_protein() then
    local map = require("plasmid.view.features").codon_map(rec, f)
    local pos = require("plasmid.model.location").positions(f.parts)
    local aa = {}
    for _, p in ipairs(pos) do if map[p] then aa[#aa + 1] = map[p] end end
    local s = table.concat(aa)
    lines[#lines + 1] = "translation (computed)"
    for k = 1, #s, 60 do lines[#lines + 1] = "  " .. s:sub(k, k + 59) end
  end
  return lines
end

--- Floating window describing a feature.
function M.show_feature(buf, fid)
  do
    local _, rec_ = state_and_record(buf)
    local f_ = rec_ and rec_:feature_by_id(fid)
    if f_ and require("plasmid.bio.primer").is_primer(f_) and pcall(require, "plasmid.bio.tm") then
      if require("plasmid.primerinfo").show_primer(buf, f_) then return end
    end
  end
  local _, rec = state_and_record(buf)
  if not rec then return end
  local f = rec:feature_by_id(fid)
  if not f then return end
  local lines = feature_lines(rec, f)
  local fbuf, fwin = vim.lsp.util.open_floating_preview(lines, "", {
    border = "rounded", max_width = math.min(90, vim.o.columns - 4), focus_id = "plasmid_feature",
  })
  local hl = require("plasmid.view.highlights").for_type(f.type)
  vim.api.nvim_buf_set_extmark(fbuf, buffer.ns, 0, 0, { end_col = #lines[1], hl_group = hl })
  return fbuf, fwin
end

local function enzyme_under_cursor(win)
  local buf = vim.api.nvim_win_get_buf(win)
  local st = buffer.get(buf)
  if not st or not st.enzymes then return nil end
  local pos = vim.api.nvim_win_get_cursor(win)
  local line = st.layout.lines[pos[1]] or ""
  local vcol = vim.fn.strdisplaywidth(line:sub(1, pos[2]))
  return require("plasmid.enzymes").entry_at(st, pos[1] - 1, vcol)
end

function M.info_under_cursor()
  local win = vim.api.nvim_get_current_win()
  local entry = enzyme_under_cursor(win)
  if entry then return require("plasmid.enzymes").show_info(vim.api.nvim_get_current_buf(), entry) end
  local _, fid = buffer.cursor_info(win)
  if fid then
    M.show_feature(vim.api.nvim_get_current_buf(), fid)
  else
    require("plasmid.util").notify("plasmid: no feature under cursor", vim.log.levels.INFO)
  end
end

--- Jump to a feature's 5'-most drawn position in the current view.
function M.goto_feature(buf, fid)
  local _, rec = state_and_record(buf)
  local f = rec and rec:feature_by_id(fid)
  if not f then return end
  local win = vim.fn.win_findbuf(buf)[1]
  local lo = f:bounds()
  if not win or not lo then return end
  buffer.goto_bp(win, lo)
  -- on the linear view, land on the feature's own track row
  local st = buffer.get(buf)
  if st.view == "linear" then
    local row0 = vim.api.nvim_win_get_cursor(win)[1] - 1
    for r = row0, row0 + 40 do
      local info = st.layout.rows[r]
      if not info or (info.kind == "blank") then break end
      for _, s in ipairs(info.segs or {}) do
        if s[3] == fid then
          local line = st.layout.lines[r + 1]
          vim.api.nvim_win_set_cursor(win, { r + 1, vim.fn.byteidx(line, s[1]) })
          return
        end
      end
    end
  end
end

--- <CR>: linear -> feature info; circular -> jump into the linear view.
function M.activate()
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if not st then return end
  local win = vim.api.nvim_get_current_win()
  local bp, fid = buffer.cursor_info(win)
  local info = st.layout.rows and st.layout.rows[vim.api.nvim_win_get_cursor(win)[1] - 1]
  if info and info.kind == "overflow" then
    return buffer.set_track_limit(buf, "all")
  end
  local entry = enzyme_under_cursor(win)
  if entry and st.view == "linear" then return require("plasmid.enzymes").show_info(buf, entry) end
  if st.view == "circular" then
    st.view = "linear"
    buffer.render(buf, { keep_bp = bp or 1 })
    if fid then M.goto_feature(buf, fid) end
  elseif fid then
    M.show_feature(buf, fid)
  end
end

local function sorted_starts(rec)
  local vf = require("plasmid.view.features")
  local out = {}
  for _, it in ipairs(vf.visible(rec, config.options)) do
    local lo = it.f:bounds()
    if lo then out[#out + 1] = { lo, it.f.id } end
  end
  table.sort(out, function(a, b) return a[1] < b[1] or (a[1] == b[1] and a[2] < b[2]) end)
  return out
end

function M.next_feature(dir)
  local buf = vim.api.nvim_get_current_buf()
  local _, rec = state_and_record(buf)
  if not rec then return end
  local bp = buffer.cursor_info(vim.api.nvim_get_current_win()) or 1
  local starts = sorted_starts(rec)
  if #starts == 0 then return end
  local pick
  if dir > 0 then
    for _, s in ipairs(starts) do if s[1] > bp then pick = s break end end
    pick = pick or starts[1]
  else
    for i = #starts, 1, -1 do if starts[i][1] < bp then pick = starts[i] break end end
    pick = pick or starts[#starts]
  end
  M.goto_feature(buf, pick[2])
end

--- Picker over all drawn features (uses vim.ui.select, so telescope/fzf-lua
--- take over automatically when they register a ui-select backend).
function M.pick_feature()
  local buf = vim.api.nvim_get_current_buf()
  local _, rec = state_and_record(buf)
  if not rec then return end
  local vf = require("plasmid.view.features")
  local items = vf.visible(rec, config.options)
  table.sort(items, function(a, b) return (a.f:bounds() or 0) < (b.f:bounds() or 0) end)
  vim.ui.select(items, {
    prompt = rec.name .. " features",
    format_item = function(it)
      local lo, hi = it.f:bounds()
      return ("%-14s %-28s %7d..%-7d %s"):format(it.f.type, it.label:sub(1, 28), lo or 0, hi or 0,
        it.f.strand == -1 and "◀" or (it.f.strand == 1 and "▶" or "·"))
    end,
  }, function(choice)
    if choice then
      M.goto_feature(buf, choice.f.id)
      local s = require("plasmid.select")
      s.set(buf, s.for_feature(require("plasmid.buffer").record(require("plasmid.buffer").get(buf)), choice.f))
    end
  end)
end

function M.attach_keymaps(buf)
  local km = config.options.keymaps
  local function map(lhs, fn, desc)
    if lhs and lhs ~= "" then
      vim.keymap.set("n", lhs, fn, { buffer = buf, silent = true, nowait = true, desc = "plasmid: " .. desc })
    end
  end
  map(km.toggle_view, function() buffer.toggle_view(buf) end, "toggle linear/circular view")
  map(km.activate, M.activate, "feature info / jump to linear view")
  map(km.info, M.info_under_cursor, "feature info")
  map(km.next_feature, function() M.next_feature(1) end, "next feature")
  map(km.prev_feature, function() M.next_feature(-1) end, "previous feature")
  map(km.features, M.pick_feature, "feature picker")
  local sel = require("plasmid.select")
  local function xo(mode, lhs, fn, desc, expr)
    if lhs and lhs ~= "" then
      vim.keymap.set(mode, lhs, fn, { buffer = buf, silent = true, nowait = true, expr = expr, desc = "plasmid: " .. desc })
    end
  end
  xo("o", km.feature_outer, function() return sel.textobj_pending("a") end, "feature span (yaf)", true)
  xo("o", km.feature_inner, function() return sel.textobj_pending("i") end, "feature residues (yif)", true)
  xo("x", km.feature_outer, function() sel.textobj_visual("a") end, "select feature span")
  xo("x", km.feature_inner, function() sel.textobj_visual("i") end, "select feature residues")
  xo("x", km.yank_bases, function() sel.yank_visual(false) end, "yank bases (top strand)")
  xo("x", km.yank_bases_rc, function() sel.yank_visual(true) end, "yank bases (reverse complement)")
  xo("x", km.yank_bases_pick, function() sel.yank_visual_pick() end, "yank bases: choose strand and direction")
  map(km.yank_selection, function() sel.yank_selection(false) end, "yank selection")
  map(km.yank_selection_rc, function() sel.yank_selection(true) end, "yank selection, reverse complement")
  map(km.clear_selection, function() sel.set(buf, nil) end, "clear selection")
  local author = require("plasmid.author")
  map(km.toggle_orfs, function() require("plasmid.orfs").toggle(buf) end, "toggle ORFs")
  map(km.add_feature, function() author.quick_add(false) end, "add feature over selection")
  xo("x", km.add_feature, function() author.quick_add(true) end, "add feature over visual range")
  map(km.edit_feature, author.edit_under_cursor, "edit feature")
  map(km.delete_feature, function() author.delete(false) end, "delete feature")
  map(km.undo, author.undo, "undo feature change")
  map(km.redo, author.redo, "redo feature change")
  xo("x", km.add_primer, function() require("plasmid.primers").add_visual(buf) end, "make a primer from the selection")
  map(km.follow_link, function() require("plasmid.fetch").follow(buf) end, "open a linked NCBI record")
  require("plasmid.mouse").attach(buf)
  map(km.clone_add, function() require("plasmid.clone").add(buf) end, "stage this record or selection as a cloning part")
  local enzui = require("plasmid.enzymes")
  map(km.toggle_enzymes, function() enzui.toggle(buf) end, "show/hide restriction sites")
  map(km.pick_enzymes, function() enzui.pick(buf) end, "pick restriction enzymes")
  local su = require("plasmid.searchui")
  map(km.search_forward, function() su.prompt(buf, 1) end, "search a motif forward")
  map(km.search_backward, function() su.prompt(buf, -1) end, "search a motif backward")
  map(km.search_next, function() su.next(buf, true, vim.v.count1) end, "next hit")
  map(km.search_prev, function() su.next(buf, false, vim.v.count1) end, "previous hit")
  map(km.select_hit, function() su.select_hit(buf) end, "select the current hit")
  local se = require("plasmid.seqedit")
  map(km.insert_before, function() se.key_insert(false) end, "insert residues before the cursor")
  map(km.insert_after, function() se.key_insert(true) end, "insert residues after the cursor")
  map(km.delete_base, function() se.key_delete(vim.v.count1) end, "delete bases")
  map(km.replace_base, function() se.key_replace(vim.v.count1) end, "replace bases")
  map(km.paste_after, function() se.key_paste(true) end, "paste residues after the cursor")
  map(km.paste_before, function() se.key_paste(false) end, "paste residues before the cursor")
  xo("x", km.delete_base, se.visual_delete, "delete selected bases")
  xo("x", km.delete_visual, se.visual_delete, "delete selected bases")
  xo("x", km.change_visual, se.visual_change, "replace selected bases")
  xo("x", km.paste_after, se.visual_paste, "replace selected bases with a register")
  xo("x", km.revcomp_visual, se.visual_revcomp, "reverse-complement selected bases")
  map(km.expand_tracks, function() buffer.set_track_limit(buf, "all") end, "show all feature tracks")
  map(km.collapse_tracks, function() buffer.set_track_limit(buf, nil) end, "limit feature tracks")
end

return M
