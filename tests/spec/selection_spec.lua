local H = require("helpers")
local pio = require("plasmid.io")
local sel = require("plasmid.select")
local alphabet = require("plasmid.model.alphabet")
local translate = require("plasmid.bio.translate")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")

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
local function keys(s)
  vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false)
end
local function feature(rec, pred)
  for _, f in ipairs(rec.features) do if pred(f) then return f end end
end
-- put the cursor on a feature's own track row
local function cursor_on(buf, f)
  require("plasmid.ui").goto_feature(buf, f.id)
  local _, fid = buffer.cursor_info(vim.api.nvim_get_current_win())
  H.eq(f.id, fid, "cursor not on " .. f:label())
end

vim.o.columns, vim.o.lines = 120, 50

case("residue_semantics_across_fixtures", function()
  local n, seen = 0, {}
  for _, path in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001422.gb", "ncbi/U01317.gb", "synthetic/grammar.gb",
    "snapgene/sgffp_origin_spanning_features.dna" }) do
    local rec = assert(pio.read(H.fixture(path))).records[1]
    for _, f in ipairs(rec.features) do
      if f.loc and not require("plasmid.model.location").has_remote(f.parts) and #rec:extents(f) > 0
        and f:length() > 0 and f.strand ~= 0 then
        local inner, outer = sel.feature_residues(rec, f, "i"), sel.feature_residues(rec, f, "a")
        H.eq(rec:extract(f), inner, path .. " " .. f.location_text)
        local s = sel.for_feature(rec, f)
        H.eq(sel.length(rec, s), #outer, path .. " span length " .. f.location_text)
        H.ok(#outer >= #inner, path .. " span shorter than residues " .. f.location_text)
        if #f.parts == 1 then H.eq(inner, outer, path .. " single-part " .. f.location_text) end
        if #f.parts > 1 and f.strand == -1 and #rec:extents(f) == 1 then seen.multi_exon_reverse = true end
        if #rec:extents(f) == 2 then seen.origin_wrap = true end
        if #f.parts > 1 and #rec:extents(f) == 1 and #outer > #inner then seen.intron = true end
        n = n + 1
      end
    end
  end
  -- the hard cases must actually have been exercised
  H.eq({ intron = true, multi_exon_reverse = true, origin_wrap = true }, seen, "checked " .. n .. " features")
end)

case("wrapping_feature_span", function()
  local rec = assert(pio.read(H.fixture("synthetic/grammar.gb"))).records[1]
  local f = feature(rec, function(x) return x:label() == "wrapCDS" end)
  local s = sel.for_feature(rec, f)
  H.eq({ 211, 30 }, { s.start, s.stop })
  H.ok(sel.wraps(s))
  H.eq(60, sel.length(rec, s))
  H.eq(rec.seq:sub(211, 240) .. rec.seq:sub(1, 30), sel.top_strand(rec, s))
end)

case("yif_on_reverse_bla_gives_coding_strand", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local bla = feature(rec, function(f) return f.type == "CDS" and f:get("gene") == "bla" end)
  H.eq(-1, bla.strand)
  cursor_on(buf, bla)
  vim.fn.setreg('"', "")
  keys("yif")
  local got = vim.fn.getreg('"')
  H.eq(rec:extract(bla):upper(), got, "yanked as displayed (upper case for GenBank)")
  H.eq("ATG", got:sub(1, 3):upper())
  H.eq(bla:get("translation"), (translate.translate(got, { table = 11, cds_start = true }):gsub("%*$", "")))
  -- named register
  keys('"ayaf')
  H.eq(got, vim.fn.getreg("a"), "single-part CDS: span == residues")
  H.eq("", vim.bo[buf].modified and "modified" or "", "yank must not modify the buffer")
end)

case("vaf_sets_selection_and_gy_yanks_it", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local tet = feature(rec, function(f) return f.type == "CDS" and f:get("gene") == "tet" end)
  cursor_on(buf, tet)
  keys("vaf")
  H.eq("n", vim.fn.mode())
  H.eq({ 86, 1276, tet.id }, { st.selection.start, st.selection.stop, st.selection.fid })
  H.ok(buffer.winbar():find("sel 86..1276 (1191 bp)", 1, true), buffer.winbar())
  keys("gy")
  H.eq(rec.seq:sub(86, 1276):upper(), vim.fn.getreg('"'))
  keys("gY")
  H.eq(alphabet.revcomp(rec.seq:sub(86, 1276)):upper(), vim.fn.getreg('"'))
  keys("<Esc>")
  H.eq(nil, st.selection)
end)

case("visual_y_and_Y_yank_exact_bases_across_lines", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local W = st.layout.W
  for _, r in ipairs({ { 10, 25 }, { 50, 50 }, { W - 3, W + 12 }, { 5, 3 * W + 7 } }) do
    buffer.goto_bp(0, r[1])
    keys(("v<Cmd>lua require('plasmid.buffer').goto_bp(0, %d)<CR>y"):format(r[2]))
    H.eq(rec.seq:sub(r[1], r[2]):upper(), vim.fn.getreg('"'), ("y %d..%d"):format(r[1], r[2]))
    H.eq({ r[1], r[2] }, { st.selection.start, st.selection.stop })
    buffer.goto_bp(0, r[2])
    keys(("v<Cmd>lua require('plasmid.buffer').goto_bp(0, %d)<CR>Y"):format(r[1]))
    H.eq(alphabet.revcomp(rec.seq:sub(r[1], r[2])):upper(), vim.fn.getreg('"'), ("Y %d..%d (backwards)"):format(r[1], r[2]))
  end
end)

case("selection_highlight_covers_exactly_the_selected_bases", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  vim.cmd("DnaSelect 77 263")
  local covered = {}
  for row0, info in pairs(st.layout.rows) do
    if info.kind == "seq" then
      local line = st.layout.lines[row0 + 1]
      for _, h in ipairs(st.layout:row_highlights(row0)) do
        if h[3] == "PlasmidSelection" then
          covered[#covered + 1] = { info.bs, (line:sub(h[1] + 1, h[2]):gsub(" ", "")) }
        end
      end
    end
  end
  table.sort(covered, function(a, b) return a[1] < b[1] end)
  local text = table.concat(vim.tbl_map(function(c) return c[2] end, covered))
  H.eq(rec.seq:sub(77, 263):upper(), text)
end)

case("origin_wrapping_selection_on_circular_record", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  vim.cmd("DnaSelect 4300 50")
  H.ok(sel.wraps(st.selection))
  H.eq(112, sel.length(rec, st.selection))
  keys("gy")
  H.eq((rec.seq:sub(4300, 4361) .. rec.seq:sub(1, 50)):upper(), vim.fn.getreg('"'))
  vim.cmd("DnaView circular")
  -- adjacent cells merge into one highlight run; count cells (braille = 3 bytes)
  local cells = 0
  for _, runs in ipairs(st.layout.runs) do
    for _, r in ipairs(runs) do
      if r[3] == "PlasmidSelectionArc" then cells = cells + (r[2] - r[1]) / 3 end
    end
  end
  -- 112 of 4361 bp is ~9 degrees of arc: at least a few cells at any window size
  H.ok(cells >= 3, "selection arc covers only " .. cells .. " cells")
  vim.cmd("DnaView linear")
end)

case("select_by_label_and_picker", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  vim.cmd("DnaSelect ROP protein")
  local rop = feature(rec, function(f) return f:label() == "ROP protein" end)
  H.eq(rop.id, st.selection.fid)
  vim.cmd("DnaSelectClear")
  H.eq(nil, st.selection)
  local orig = vim.ui.select
  vim.ui.select = function(items, _, cb)
    for _, it in ipairs(items) do if it.f.id == rop.id then return cb(it) end end
  end
  require("plasmid.ui").pick_feature()
  vim.ui.select = orig
  H.eq(rop.id, st.selection.fid)
end)

case("yaf_without_feature_is_harmless", function()
  local buf, st, rec = open("ncbi/L09137.gb") -- only a source feature
  vim.fn.setreg('"', "sentinel")
  buffer.goto_bp(0, 100)
  keys("yaf")
  H.eq("sentinel", vim.fn.getreg('"'))
end)

return T
