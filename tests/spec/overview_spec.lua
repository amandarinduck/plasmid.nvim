-- The overview claims two things: every base of the record is on exactly one
-- row, and what a row says about its span is true. These check both against
-- the record itself rather than against another part of the view, and then
-- drive the commands the way a user does.
local H = require("helpers")

require("plasmid").setup({})
local io_ = require("plasmid.io")
local config = require("plasmid.config")
local ov = require("plasmid.view.overview")
local overview = require("plasmid.overview")
local buffer = require("plasmid.buffer")
local vf = require("plasmid.view.features")
local Record = require("plasmid.model.record").Record

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function read(rel)
  local res = assert(io_.read(H.fixture(rel)), "could not read " .. rel)
  return res.records[1]
end

local function opts()
  return vim.deepcopy(config.options)
end

--- Features overlapping [bs, be], counted directly off the record: the
--- overview's own binning is what we are checking, so it cannot be used here.
local function direct_counts(rec, bs, be, o)
  local fwd, rev, other = 0, 0, 0
  for _, it in ipairs(vf.visible(rec, o)) do
    local hit = false
    for _, e in ipairs(it.extents) do
      if e[1] <= be and e[2] >= bs then hit = true end
    end
    if hit then
      if it.f.strand == 1 then fwd = fwd + 1
      elseif it.f.strand == -1 then rev = rev + 1
      else other = other + 1 end
    end
  end
  return fwd, rev, other
end

case("the_span_is_a_round_number_that_fits_the_record", function()
  -- the ladder itself: 1, 2 or 5 times a power of ten, never below the ask
  H.eq({ 1, 10, 20, 200, 500, 1000, 200000, 1000000 },
    { ov.ladder(1), ov.ladder(10), ov.ladder(11), ov.ladder(109), ov.ladder(201),
      ov.ladder(501), ov.ladder(100905), ov.ladder(500001) })
  local o = opts()
  local sizes = {
    { 300, 42, 10 },            -- far smaller than a screen: the finest span
    { 4361, 42, 200 },          -- pBR322
    { 230218, 42, 10000 },      -- yeast chromosome I
    { 4641652, 42, 200000 },    -- E. coli K-12
    { 4641652, 10, 500000 },    -- the same genome in a ten-row window
  }
  for _, s in ipairs(sizes) do
    local span = ov.span_for(s[1], s[2], o)
    H.eq(s[3], span, ("span for %d bases in %d rows"):format(s[1], s[2]))
    H.ok(math.ceil(s[1] / span) <= s[2],
      ("%d bases at %d per row needs %d rows, not %d"):format(s[1], span, math.ceil(s[1] / span), s[2]))
  end
  -- a configured span is taken as given, screen or no screen
  o.overview.span = 1000
  H.eq(1000, ov.span_for(4641652, 42, o), "a configured span wins")
  H.log(("    [overview] span ladder: 300 bp → %d, 230 kb → %d, 4.64 Mb → %d bases per row (42 rows)")
    :format(ov.span_for(300, 42, opts()), ov.span_for(230218, 42, opts()), ov.span_for(4641652, 42, opts())))
end)

case("a_tiny_record_still_gets_a_picture", function()
  local rec = Record.new({ name = "tiny", seq = ("ACGT"):rep(75), format = "genbank" })
  local lay = ov.build(rec, opts(), 100, 46, {})
  H.eq(10, lay.span, "300 bp in a 46-row window")
  H.eq(30, lay.nrows)
  H.eq(lay.header_rows + lay.nrows, #lay.lines, "one row per span plus the header")
  H.eq({ 291, 300 }, { lay:region_of_row(#lay.lines - 1) }, "last row")
  -- an empty record must not be a special case anywhere
  local empty = Record.new({ name = "empty", seq = "", format = "genbank" })
  local elay = ov.build(empty, opts(), 100, 46, {})
  H.eq(1, elay.nrows, "an empty record is one empty row")
  H.ok(elay.lines[#elay.lines]:find("no residues", 1, true) ~= nil, "no note on an empty record")
end)

case("rows_and_regions_map_both_ways", function()
  for _, rel in ipairs({ "ncbi/J01749.gb", "ncbi/NC_005816.gb", "ncbi/NC_001133.gb" }) do
    local rec = read(rel)
    local lay = ov.build(rec, opts(), 110, 40, {})
    local covered = 0
    for row = lay.header_rows, #lay.lines - 1 do
      local bs, be = lay:region_of_row(row)
      H.ok(bs ~= nil, rel .. ": no region on row " .. row)
      covered = covered + (be - bs + 1)
      -- the row a base is on, and the bases a row stands for, are inverses
      H.eq(row, lay:row_of_bp(bs), rel .. ": row of first base of row " .. row)
      H.eq(row, lay:row_of_bp(be), rel .. ": row of last base of row " .. row)
      H.eq(row, lay:row_of_bp(math.floor((bs + be) / 2)), rel .. ": row of a middle base")
      if row > lay.header_rows then
        H.eq(bs - 1, select(2, lay:region_of_row(row - 1)), rel .. ": rows must not leave a gap")
      end
      -- a column inside the strip resolves finer than the whole region
      local info = lay.rows[row]
      H.eq(bs, lay:bp_at(row, 0), rel .. ": outside the strip is the region start")
      H.eq(bs, lay:bp_at(row, info.strip_col), rel .. ": first cell")
      local last = lay:bp_at(row, info.strip_col + lay.cells - 1)
      H.ok(last > bs or lay.span <= lay.cells, rel .. ": last cell did not advance")
      H.ok(last <= be, rel .. ": last cell ran past the region")
    end
    H.eq(rec:len(), covered, rel .. ": the rows must cover the record exactly once")
    H.eq(lay.header_rows, lay:row_of_bp(1), rel .. ": base 1")
    H.eq(#lay.lines - 1, lay:row_of_bp(rec:len()), rel .. ": last base")
    H.eq(nil, lay:region_of_row(0), rel .. ": the title is not a region")
  end
  H.log("    [overview] row ↔ region is a bijection on 3 records, cells resolve inside a row")
end)

case("density_counts_match_a_direct_count", function()
  local rows, feats = 0, 0
  for _, rel in ipairs({ "ncbi/J01749.gb", "ncbi/NC_005816.gb", "ncbi/NC_001133.gb", "ncbi/U01317.gb" }) do
    local rec = read(rel)
    local o = opts()
    local lay = ov.build(rec, o, 110, 40, {})
    local total = 0
    for row = lay.header_rows, #lay.lines - 1 do
      local info = lay.rows[row]
      local fwd, rev, other = direct_counts(rec, info.bs, info.be, o)
      H.eq({ fwd, rev, other }, { info.fwd, info.rev, info.other },
        ("%s: counts on %d..%d"):format(rel, info.bs, info.be))
      H.eq(fwd + rev + other, #info.items, ("%s: items listed on %d..%d"):format(rel, info.bs, info.be))
      -- and what the row actually says: a strandless feature is on both strands
      H.ok(lay.lines[row + 1]:find(("%4d+ %4d−"):format(fwd + other, rev + other), 1, true) ~= nil,
        ("%s: row %d does not show %d+ %d−"):format(rel, row, fwd + other, rev + other))
      rows, total = rows + 1, total + fwd + rev + other
    end
    -- every drawn feature reaches at least one row
    H.ok(total >= #vf.visible(rec, o), rel .. ": features were lost between the rows")
    feats = feats + #vf.visible(rec, o)
  end
  H.log(("    [overview] %d rows of 4 records: per-strand counts equal a direct count over %d features")
    :format(rows, feats))
end)

case("a_row_is_never_wider_than_the_window", function()
  for _, w in ipairs({ 80, 100, 130, 200 }) do
    local rec = read("ncbi/NC_001133.gb")
    local lay = ov.build(rec, opts(), w, 40, { cursor_bp = 100000, visible = { 100000, 102000 } })
    for i, line in ipairs(lay.lines) do
      H.ok(vim.fn.strdisplaywidth(line) <= w, ("width %d: row %d is %d columns wide")
        :format(w, i, vim.fn.strdisplaywidth(line)))
    end
    -- highlight runs stay inside their line
    for _, h in ipairs(lay.hls) do
      H.ok(h[3] <= #lay.lines[h[1] + 1] + 1, ("width %d: a highlight runs past row %d"):format(w, h[1]))
      H.ok(h[2] < h[3], "empty highlight run")
    end
  end
  H.log("    [overview] rows fit 80, 100, 130 and 200 columns, highlights inside the text")
end)

-- Protein records, empty ones, SnapGene files, multi-record FASTA: the
-- overview has to have an answer for every record the plugin can open.
case("every_fixture_record_gets_an_overview", function()
  local files, records, rows = 0, 0, 0
  local paths = {}
  for _, pat in ipairs({ "**/*.gb", "**/*.gbk", "**/*.fa", "**/*.fasta", "**/*.dna", "**/*.prot" }) do
    for _, p in ipairs(H.fixtures(pat)) do paths[#paths + 1] = p end
  end
  table.sort(paths)
  for _, path in ipairs(paths) do
    local res = io_.read(path)
    if res then
      files = files + 1
      for ri, rec in ipairs(res.records) do
        records = records + 1
        local tag = ("%s r%d"):format(vim.fn.fnamemodify(path, ":t"), ri)
        local lay = ov.build(rec, opts(), 110, 40, { cursor_bp = 1 })
        H.eq(lay.header_rows + lay.nrows + (rec:len() == 0 and 2 or 0), #lay.lines, tag .. ": row count")
        H.eq(math.max(rec:len(), 1), select(2, lay:region_of_row(#lay.lines - 1 - (rec:len() == 0 and 2 or 0))),
          tag .. ": last row ends at the last residue")
        for i, line in ipairs(lay.lines) do
          H.ok(vim.fn.strdisplaywidth(line) <= 110, ("%s: row %d is too wide"):format(tag, i))
        end
        rows = rows + lay.nrows
      end
    end
  end
  H.ok(files > 40, "expected the fixture corpus, got " .. files .. " files")
  H.log(("    [overview] %d files, %d records (nucleotide, protein and empty), %d rows built")
    :format(files, records, rows))
end)

case("the_viewport_is_marked", function()
  local rec = read("ncbi/NC_001133.gb")
  local lay = ov.build(rec, opts(), 110, 40, { cursor_bp = 100000, visible = { 98000, 112000 } })
  local cursor_row, marked = nil, 0
  for row = lay.header_rows, #lay.lines - 1 do
    local line = lay.lines[row + 1]
    if line:sub(1, 3) == "▶" then cursor_row = row end
    if line:sub(1, 3) == "┃" then marked = marked + 1 end
  end
  H.eq(lay:row_of_bp(100000), cursor_row, "the cursor's row is not marked")
  H.ok(marked >= 1, "no rows marked for a 14 kb viewport at 10 kb per row")
  for row = lay.header_rows, #lay.lines - 1 do
    local info = lay.rows[row]
    local overlaps = info.be >= 98000 and info.bs <= 112000
    local gut = lay.lines[row + 1]:sub(1, 3)
    H.eq(overlaps, gut == "▶" or gut == "┃", ("row %d..%d marked wrongly"):format(info.bs, info.be))
  end
  H.log(("    [overview] a 14 kb viewport marks %d rows plus the cursor's own"):format(marked))
end)

case("a_230_kb_record_fits_one_screen", function()
  local rec = read("ncbi/NC_001133.gb")
  local t0 = vim.uv.hrtime()
  local lay = ov.build(rec, opts(), 130, 46, {})
  local ms = (vim.uv.hrtime() - t0) / 1e6
  H.ok(#lay.lines <= 46, ("the overview is %d rows in a 46-row window"):format(#lay.lines))
  H.eq(10000, lay.span)
  H.ok(ms < 500, ("building the overview took %.0f ms"):format(ms))
  local named = 0
  for row = lay.header_rows, #lay.lines - 1 do
    if lay.lines[row + 1]:find("%a%a%a") then named = named + 1 end
  end
  H.log(("    [overview] NC_001133 (230218 bp, %d features): %d rows of %d kb, %d named, %.1f ms")
    :format(#rec.features, lay.nrows, lay.span / 1000, named, ms))
end)

-- Large records are not committed; when one has been fetched into big/, the
-- same build is measured on it. Nothing here depends on it being there.
case("a_bacterial_genome_if_one_has_been_fetched", function()
  local path = H.root .. "/big/NC_000913.gbk"
  if vim.fn.filereadable(path) == 0 then
    return H.log("    [overview] big/NC_000913.gbk not present, genome timing skipped")
  end
  local rec = assert(io_.read(path)).records[1]
  local t0 = vim.uv.hrtime()
  local lay = ov.build(rec, opts(), 130, 46, { cursor_bp = 1, visible = { 1, 5000 } })
  local ms = (vim.uv.hrtime() - t0) / 1e6
  H.ok(#lay.lines <= 46, ("the genome overview is %d rows"):format(#lay.lines))
  local fwd, rev, other = direct_counts(rec, lay.rows[lay.header_rows].bs, lay.rows[lay.header_rows].be, opts())
  H.eq({ fwd, rev, other },
    { lay.rows[lay.header_rows].fwd, lay.rows[lay.header_rows].rev, lay.rows[lay.header_rows].other },
    "counts on the genome's first row")
  H.log(("    [overview] %s (%d bp, %d features): %d rows of %d kb, %.0f ms")
    :format(rec.name, rec:len(), #rec.features, lay.nrows, lay.span / 1000, ms))
end)

---------------------------------------------------------------------------
-- The commands, driven the way a user drives them
---------------------------------------------------------------------------

local function open_fixture(rel)
  local tmp = vim.fn.tempname()
  vim.fn.mkdir(tmp, "p")
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(dst))
  return vim.api.nvim_get_current_buf(), vim.api.nvim_get_current_win()
end

case("the_command_opens_the_overview_over_the_record", function()
  local buf, win = open_fixture("ncbi/NC_001133.gb")
  buffer.goto_bp(win, 123456)
  vim.cmd("DnaOverview")
  local obuf = vim.api.nvim_win_get_buf(win)
  H.ok(obuf ~= buf, "the overview did not take the window")
  local S = overview.state[obuf]
  H.ok(S ~= nil, "no state for the overview buffer")
  H.eq(buf, S.src)
  H.eq(false, vim.bo[obuf].modifiable, "the overview is editable")
  -- it opened on the row the sequence view was looking at
  H.eq(S.lay:row_of_bp(123456) + 1, vim.api.nvim_win_get_cursor(win)[1], "cursor row on opening")
  H.eq(vim.api.nvim_buf_get_lines(obuf, 0, -1, false), S.lay.lines, "buffer text differs from the layout")
  -- and the command closes it again, putting the record back where it was
  vim.cmd("DnaOverview")
  H.eq(buf, vim.api.nvim_win_get_buf(win), "the record did not come back")
  H.eq(123456, (buffer.cursor_info(win)), "the sequence view moved while the overview was open")
  vim.cmd("silent! %bwipeout!")
end)

case("enter_jumps_to_the_region_under_the_cursor", function()
  local buf, win = open_fixture("ncbi/NC_001133.gb")
  buffer.goto_bp(win, 1)
  vim.cmd("DnaOverview")
  local obuf = vim.api.nvim_win_get_buf(win)
  local lay = overview.state[obuf].lay
  local row = lay.header_rows + math.floor(lay.nrows / 2)
  vim.api.nvim_win_set_cursor(win, { row + 1, 0 })
  local bs, be = lay:region_of_row(row)
  -- the key the user presses, not the function behind it
  local m = vim.fn.maparg(config.options.overview.keymaps.jump, "n", false, true)
  H.ok(m and m.callback, "<CR> is not mapped in the overview")
  m.callback()
  H.eq(buf, vim.api.nvim_win_get_buf(win), "the jump did not restore the record")
  local bp = buffer.cursor_info(win)
  H.ok(bp >= bs and bp <= be, ("jumped to %s, outside %d..%d"):format(tostring(bp), bs, be))
  local line = vim.api.nvim_buf_get_lines(buf, vim.api.nvim_win_get_cursor(win)[1] - 1,
    vim.api.nvim_win_get_cursor(win)[1], false)[1]
  H.ok(line ~= "" and line:match("%u%u%u"), "landed on a row that was not laid out: " .. vim.inspect(line))
  H.eq(nil, overview.state[obuf], "the overview outlived the jump")
  H.log(("    [overview] <CR> on row %d landed on base %d of %d..%d"):format(row, bp, bs, be))
  vim.cmd("silent! %bwipeout!")
end)

case("a_click_jumps_inside_the_row", function()
  local _, win = open_fixture("ncbi/NC_001133.gb")
  vim.cmd("DnaOverview")
  local obuf = vim.api.nvim_win_get_buf(win)
  local lay = overview.state[obuf].lay
  local row = lay.header_rows + 5
  local info = lay.rows[row]
  -- the column the mouse would report, three quarters along the strip
  local col = info.strip_col + math.floor(lay.cells * 3 / 4)
  local want = lay:bp_at(row, col)
  H.ok(want > info.bs + lay.span / 2, "a click late in the row should land late in the region")
  local got = overview.jump(obuf, { row, col })
  H.eq(want, got, "the click did not land where the column points")
  H.eq(want, (buffer.cursor_info(win)), "the cursor is not on the clicked base")
  H.log(("    [overview] a click %d cells along a %d kb row landed on base %d of %d..%d")
    :format(math.floor(lay.cells * 3 / 4), lay.span / 1000, got, info.bs, info.be))
  vim.cmd("silent! %bwipeout!")
end)

case("the_span_can_be_given_to_the_command", function()
  local _, win = open_fixture("ncbi/NC_001133.gb")
  vim.cmd("DnaOverview 1000")
  local obuf = vim.api.nvim_win_get_buf(win)
  local lay = overview.state[obuf].lay
  H.eq(1000, lay.span, "the span argument was ignored")
  H.eq(231, lay.nrows, "230218 bp at 1 kb per row")
  H.ok(#lay.lines > vim.fn.winheight(win), "a forced span is allowed to overflow the screen")
  overview.close(obuf)
  vim.cmd("silent! %bwipeout!")
end)

case("the_keymap_opens_and_closes_it", function()
  local buf, win = open_fixture("ncbi/NC_005816.gb")
  local m = vim.fn.maparg(config.options.keymaps.overview, "n", false, true)
  H.ok(m and m.callback, config.options.keymaps.overview .. " is not mapped in a sequence buffer")
  m.callback()
  local obuf = vim.api.nvim_win_get_buf(win)
  H.ok(overview.state[obuf] ~= nil, "the keymap did not open the overview")
  local back = vim.fn.maparg(config.options.overview.keymaps.close, "n", false, true)
  back.callback()
  H.eq(buf, vim.api.nvim_win_get_buf(win), "q did not put the record back")
  vim.cmd("silent! %bwipeout!")
end)

case("a_region_is_laid_out_without_the_rest_of_the_record", function()
  local prev = config.options.linear.window
  config.options.linear.window = { enabled = true, margin = 1, cache_blocks = 64 }
  local buf, win = open_fixture("ncbi/NC_001133.gb")
  local st = buffer.get(buf)
  H.ok(st.layout.windowed, "the fixture did not open windowed")
  vim.cmd("DnaRegion 150000..151000")
  local bp = buffer.cursor_info(win)
  H.eq(150000, bp, "DnaRegion did not go to the start of the region")
  local sel = require("plasmid.select").get(buf)
  H.eq({ 150000, 151000 }, { sel.start, sel.stop }, "the region was not selected")
  local built = 0
  for _ in pairs(st.layout._built) do built = built + 1 end
  H.ok(built < st.layout.nblocks / 4,
    ("%d of %d blocks laid out to show one region"):format(built, st.layout.nblocks))
  local laid = 0
  for i = st.layout.header_rows + 1, #st.layout.lines do
    if st.layout.lines[i] ~= "" then laid = laid + 1 end
  end
  H.log(("    [overview] :DnaRegion 150000..151000 laid out %d of %d blocks, %d of %d rows")
    :format(built, st.layout.nblocks, laid, #st.layout.lines))
  config.options.linear.window = prev
  vim.cmd("silent! %bwipeout!")
end)

return T
