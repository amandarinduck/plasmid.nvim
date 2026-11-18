local H = require("helpers")
local pio = require("plasmid.io")
local translate = require("plasmid.bio.translate")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local orfs = require("plasmid.orfs")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function open(rel, name)
  local dst = tmp .. "/" .. (name or vim.fn.fnamemodify(rel, ":t"))
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.get(buf), buffer.record(buffer.get(buf)), dst
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end
local function drawn_orfs(st)
  local ids = {}
  for _, info in pairs(st.layout.rows or {}) do
    for _, it in ipairs(info.items or {}) do if it.f.overlay then ids[it.f.id] = true end end
  end
  return vim.tbl_count(ids)
end
local function cursor_on(buf, fid)
  require("plasmid.ui").goto_feature(buf, fid)
  local _, got = buffer.cursor_info(vim.api.nvim_get_current_win())
  H.eq(fid, got)
end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("toggle_draws_orfs_and_hides_annotated", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  keys("go")
  H.ok(rec.overlays and #rec.overlays > 0, "no ORFs")
  H.eq(#rec.overlays, drawn_orfs(st), "every ORF drawn (tracks unlimited for overlays here)")
  -- no ORF shares a stop codon with an annotated CDS
  for _, o in ipairs(rec.overlays) do
    local three = o.orf.strand == 1 and o.orf.stop or o.orf.start
    for _, f in ipairs(rec.features) do
      if f.type == "CDS" and f.strand == o.orf.strand then
        local lo, hi = f:bounds()
        H.ok(three ~= (f.strand == 1 and hi or lo), "ORF duplicates annotated " .. f:label())
      end
    end
    H.ok(o.orf.aa_len >= 75)
  end
  -- translation rows appear under ORFs
  local aa_rows = 0
  for _, info in pairs(st.layout.rows) do
    if info.kind == "aa" then
      for _, s in ipairs(info.segs) do if s[3] < 0 then aa_rows = aa_rows + 1 end end
    end
  end
  H.ok(aa_rows > 0, "no ORF translation drawn")
  vim.cmd("DnaOrfs 150")
  for _, o in ipairs(rec.overlays) do H.ok(o.orf.aa_len >= 150) end
  keys("go")
  H.eq(nil, rec.overlays)
  H.eq(0, drawn_orfs(st))
end)

case("orfs_are_never_saved_and_not_editable", function()
  local buf, st, rec, path = open("ncbi/J01749.gb", "orf_nosave.gb")
  local before = H.read(path)
  vim.cmd("DnaOrfs 75")
  H.eq(false, vim.bo[buf].modified)
  vim.cmd("silent write!")
  H.eq(before, H.read(path))
  local o = rec.overlays[1]
  cursor_on(buf, o.id)
  local wins = #vim.api.nvim_list_wins()
  keys("gE")
  H.eq(wins, #vim.api.nvim_list_wins(), "editor must not open for an ORF")
  vim.cmd("DnaFeatureDelete!")
  H.ok(rec:feature_by_id(o.id), "ORF must not be deletable")
  -- K and yif work on ORFs
  local fbuf = require("plasmid.ui").show_feature(buf, o.id)
  H.ok(fbuf and vim.api.nvim_buf_get_lines(fbuf, 0, 1, false)[1]:find("ORF", 1, true))
  keys("yif")
  local nt = vim.fn.getreg('"')
  H.eq("ATG", nt:sub(1, 3):upper())
  H.eq(o.orf.aa_len + 1, #nt / 3)
end)

case("promote_orf_to_cds_undo_and_save", function()
  local buf, st, rec, path = open("ncbi/J01749.gb", "orf_promote.gb")
  vim.cmd("DnaOrfs 75")
  local n_orfs, n_feats = #rec.overlays, #rec.features
  local o = rec.overlays[1]
  local loc = o.location_text
  cursor_on(buf, o.id)
  local orig = vim.ui.input
  vim.ui.input = function(_, cb) cb("hypothetical") end
  vim.cmd("DnaOrfAdd")
  vim.ui.input = orig
  H.eq(n_feats + 1, #rec.features)
  H.eq(n_orfs - 1, #rec.overlays, "promoted ORF leaves the overlay")
  local cds
  for _, f in ipairs(rec.features) do if f:label() == "hypothetical" then cds = f end end
  H.eq({ "CDS", loc }, { cds.type, cds.location_text })
  local aa = translate.translate(rec:extract(cds), { table = 11, cds_start = true }):gsub("%*$", "")
  H.eq(aa, cds:get("translation"))
  H.eq(o.orf.aa_len, #aa)
  H.eq(true, vim.bo[buf].modified)
  keys("u")
  H.eq(n_feats, #rec.features)
  H.eq(n_orfs, #rec.overlays, "undo brings the ORF back")
  keys("<C-r>")
  vim.cmd("silent write")
  local back = assert(pio.read(path)).records[1]
  local saved
  for _, f in ipairs(back.features) do if f:label() == "hypothetical" then saved = f end end
  H.eq(loc, saved.location_text)
  H.eq(aa, saved:get("translation"))
end)

case("ga_on_orf_promotes_and_snapgene_saves", function()
  local buf, st, rec, path = open("snapgene/pFA-KanMX4.dna", "orf_ga.dna")
  vim.cmd("DnaOrfs 50")
  local o = rec.overlays[1]
  cursor_on(buf, o.id)
  local orig = vim.ui.input
  vim.ui.input = function(_, cb) cb("orfX") end
  keys("ga")
  vim.ui.input = orig
  vim.cmd("silent write")
  local back = assert(pio.read(path)).records[1]
  local saved
  for _, f in ipairs(back.features) do if f:label() == "orfX" then saved = f end end
  H.ok(saved, "promoted ORF not written to SnapGene")
  H.eq(o.location_text, saved.location_text)
end)

case("circular_view_shows_orfs", function()
  local buf, st, rec = open("ncbi/NC_001422.gb", "phix.gb")
  vim.cmd("DnaOrfs 30")
  vim.cmd("DnaView circular")
  local cells = 0
  for _, ow in ipairs(st.layout.owners) do
    for _, fid in pairs(ow) do if fid < 0 then cells = cells + 1 end end
  end
  H.ok(cells > 0, "no ORF drawn on the map")
end)


case("header_counts_overlays_separately", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local n = 0
  for _, f in ipairs(rec.features) do if f.type ~= "source" then n = n + 1 end end
  local before = vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1]
  H.ok(not before:find("ORF"), before)
  keys("go")
  local norf = #rec.overlays
  local after = vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1]
  local feats = tonumber(after:match("(%d+) features"))
  H.eq(tonumber(before:match("(%d+) features")), feats, "annotated count unchanged by ORFs: " .. after)
  H.ok(after:find(("%d ORFs"):format(norf), 1, true), after)
  keys("go")
end)

case("switching_to_circular_shows_whole_map", function()
  local buf = open("ncbi/J01749.gb")
  local win = vim.api.nvim_get_current_win()
  buffer.goto_bp(win, 3500)   -- far down the linear view
  H.ok(vim.fn.line("w0") > 1, "precondition: scrolled")
  buffer.set_view(buf, "circular")
  -- Neovim 0.10 does not keep line("w0") current in headless `nvim -l` runs
  -- (it can even report w$ < w0), so the scroll check needs 0.11+. In a
  -- terminal, 0.10 shows the whole map here too.
  if vim.fn.has("nvim-0.11") == 1 and vim.api.nvim_buf_line_count(buf) <= vim.fn.winheight(win) then
    H.eq(1, vim.fn.line("w0"), "map fits, so it is shown from the top")
  end
  buffer.set_view(buf, "linear")
  local l = vim.fn.line(".")
  H.ok(l >= vim.fn.line("w0") and l <= vim.fn.line("w$"), "cursor visible after switching back")
end)

return T
