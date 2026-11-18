local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local E = require("plasmid.enzymes")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent! only")
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if st.view ~= "linear" then buffer.set_view(buf, "linear") end
  return buf, st, buffer.record(st)
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end
local function entry(st, label)
  for _, g in ipairs(st.enzymes.entries) do
    if g.label == label or vim.tbl_contains(g.names, label) then return g end
  end
end
-- display column of the "▏" drawn for an entry on rows of a kind, and the row
local function mark_col(st, g, kind, pos)
  local lay = st.layout
  if kind == "enz" then
    -- "▏" over the base after the cut, or "▕" over the base before it
    local c, r = mark_col(st, g, "enz_right", pos)
    if c then return c, r end
    return mark_col(st, g, "enz_left", pos - 1)
  end
  local glyph = kind == "enz_left" and "▕" or "▏"
  local want_kind = kind
  if kind == "enz_left" or kind == "enz_right" then want_kind = "enz" end
  local brow = lay:pos_of(pos)
  for r = brow - 6, brow + 3 do
    local info = lay.rows[r]
    if info and info.kind == want_kind and pos >= info.bs and pos <= info.be then
      local line = lay.lines[r + 1]
      local want_col = select(2, lay:pos_of(pos))
      local b = vim.fn.byteidx(line, vim.fn.charidx(line, 0)) -- 0
      for c = 0, vim.fn.strchars(line) - 1 do
        if vim.fn.strcharpart(line, c, 1) == glyph and vim.fn.strdisplaywidth(vim.fn.strcharpart(line, 0, c)) == want_col then
          if want_kind ~= "enz" then return want_col, r end
          if E.entry_at(st, r, want_col) == g then return want_col, r end
        end
      end
    end
  end
end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("unique_set_marks_cuts_on_both_strands", function()
  local buf, st = open("ncbi/L09137.gb")
  keys("gz")
  H.eq("unique", st.enzymes.spec.mode)
  local g = assert(entry(st, "EcoRI"), "EcoRI shown")
  H.eq(397, g.cut); H.eq(401, g.cut_bottom); H.eq(1, g.count)
  H.ok(mark_col(st, g, "enz", 397), "top-strand cut mark over bp 397 (cut 396^397)")
  H.ok(mark_col(st, g, "enzcut", 401), "bottom-strand cut mark under bp 401")
  -- isoschizomers with the same cuts share one label
  H.ok(entry(st, "PstI, SbfI"), "PstI/SbfI merged")
  for _, x in ipairs(st.enzymes.entries) do
    H.eq(1, x.count, x.label .. " is unique")
    for _, e in ipairs(x.enzymes) do H.ok(E.specific_len(e.site) >= 6, e.name .. " has a 6+ site") end
  end
  H.ok(buffer.winbar():find("RE unique: ", 1, true))
  -- the site is underlined on the sequence row
  local row = st.layout:pos_of(396)
  local n = 0
  for _, h in ipairs(st.layout:row_highlights(row)) do if h[3] == "PlasmidEnzymeSite" then n = n + 1 end end
  H.ok(n > 0, "site underline")
  keys("gz")
  H.eq(nil, st.enzymes)
end)

case("named_dual_and_all_sets", function()
  local buf, st = open("ncbi/L09137.gb")
  vim.cmd("DnaEnzymes EcoRI BamHI NotI bogus")
  H.eq("names", st.enzymes.spec.mode)
  H.eq(2, #st.enzymes.entries)
  H.eq({ "NotI" }, st.enzymes.summary.zero)
  H.eq({ "bogus" }, st.enzymes.summary.unknown)
  vim.cmd("DnaEnzymes dual")
  local two = false
  for _, x in ipairs(st.enzymes.entries) do
    H.ok(x.count == 1 or x.count == 2, x.label)
    if x.count == 2 then two = true end
  end
  H.ok(two, "some 2-cutters")
  vim.cmd("DnaEnzymes all")
  local many = 0
  for _, x in ipairs(st.enzymes.entries) do if x.count > 2 then many = many + 1 end end
  H.ok(many > 0, "multi-cutters in the all set")
  -- DpnI needs methylated DNA: only shown when named
  for _, x in ipairs(st.enzymes.entries) do H.ok(not vim.tbl_contains(x.names, "DpnI"), "DpnI in all set") end
  vim.cmd("DnaEnzymes DpnI")
  H.ok(#st.enzymes.entries > 0)
  vim.cmd("DnaEnzymes off")
  H.eq(nil, st.enzymes)
end)

case("origin_spanning_site", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  vim.cmd("DnaEnzymes EcoRI")
  H.eq(1, #st.enzymes.entries)
  local g = st.enzymes.entries[1]
  H.eq(4360, g.cut); H.ok(g.stop < g.start, "site wraps the origin")
  -- underline on the last and the first line
  local function underlined(bp)
    local row = st.layout:pos_of(bp)
    local col = select(2, st.layout:pos_of(bp))
    for _, h in ipairs(st.layout:row_highlights(row)) do
      if h[3] == "PlasmidEnzymeSite" and col >= h[1] and col < h[2] then return true end
    end
  end
  H.ok(underlined(rec:len()), "last base underlined")
  H.ok(underlined(1), "first base underlined")
  local lines = E.info_lines(rec, g)
  H.ok(table.concat(lines, "\n"):find("across the origin", 1, true))
end)

case("info_float_and_edit_refresh", function()
  local buf, st, rec = open("ncbi/L09137.gb")
  vim.cmd("DnaEnzymes unique")
  local g = entry(st, "EcoRI")
  local text = table.concat(E.info_lines(rec, g), "\n")
  H.ok(text:find("cuts     top 396^397   bottom 400^401", 1, true), text)
  H.ok(text:find("5'+4", 1, true))
  -- K on the label opens the float
  local c, r = mark_col(st, g, "enz", 397)
  vim.api.nvim_win_set_cursor(0, { r + 1, 0 })
  vim.cmd("normal! 0" .. (c + 1) .. "|")
  local before = #vim.api.nvim_list_wins()
  keys("K")
  H.eq(before + 1, #vim.api.nvim_list_wins(), "info float opened")
  vim.cmd("fclose!")
  -- an edit inside the site removes EcoRI from the unique set; undo restores it
  buffer.goto_bp(0, 398)
  keys("x")
  H.eq(nil, entry(st, "EcoRI"))
  keys("u")
  H.ok(entry(st, "EcoRI"), "EcoRI back after undo")
end)

case("circular_labels_and_picker", function()
  local buf, st = open("ncbi/L09137.gb")
  vim.cmd("DnaEnzymes unique")
  buffer.set_view(buf, "circular")
  local L = st.layout
  H.ok(#L.enz_labels > 0, "labels placed")
  local l = L.enz_labels[1]
  H.eq(l.entry, E.entry_at(st, l.row + L.header_rows, l.col))
  local cut = 0
  for row0 = 0, #L.lines - 1 do
    for _, h in ipairs(L:row_highlights(row0)) do if h[3] == "PlasmidEnzymeCut" then cut = cut + 1 end end
  end
  H.ok(cut > 0, "cut ticks drawn")
  buffer.set_view(buf, "linear")
  -- picker toggles an enzyme in the named set
  local orig = vim.ui.select
  vim.ui.select = function(items, opts, cb)
    for _, it in ipairs(items) do if it.e.name == "HindIII" then
      H.ok(opts.format_item(it):find("1 site", 1, true)); return cb(it) end end
  end
  keys("gZ")
  H.eq({ "HindIII" }, st.enzymes.spec.names)
  keys("gZ")
  vim.ui.select = orig
  H.eq(nil, st.enzymes, "second pick removes it")
end)

case("rebase_import_command_switches_the_set", function()
  local rb = require("plasmid.rebase")
  local R = require("plasmid.bio.restriction")
  local saved = rb.path
  rb.path = tmp .. "/rebase_enzymes.json"
  rb.reset(); R.reset()
  local fx = vim.fn.fnamemodify(H.fixture("rebase/emboss_e.test"), ":h")
  vim.cmd(("DnaEnzymeImport %s/emboss_e.test %s/emboss_r.test %s/emboss_s.test"):format(fx, fx, fx))
  local list = R.enzymes()
  H.ok(#list > 0 and #list < 30, "imported fixture set in use: " .. #list)
  H.ok(R.get("EcoRI"), "EcoRI imported")
  H.eq(1, R.get("EcoRI").cut5)
  local buf, st = open("ncbi/L09137.gb")
  vim.cmd("DnaEnzymes all")
  for _, g in ipairs(st.enzymes.entries) do
    for _, e in ipairs(g.enzymes) do H.ok(e.suppliers ~= "", e.name .. " is commercial") end
  end
  vim.cmd("DnaEnzymeImport clear")
  H.eq(222, #R.enzymes(), "back to the bundled table")
  H.eq(0, vim.fn.filereadable(rb.path))
  rb.path = saved
  rb.reset(); R.reset()
end)

return T
