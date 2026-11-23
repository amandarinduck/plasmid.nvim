-- The alignment UI: :DnaAlign, :DnaAlignReads, :DnaCompare, the view and its
-- keymaps, the difference list and annotating the differences onto a record.
-- The engine itself is tested in align_spec.
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local AL = require("plasmid.alignment")
local C = require("plasmid.clone")
local alphabet = require("plasmid.model.alphabet")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")

vim.o.columns, vim.o.lines, vim.o.hidden = 150, 50, true

local function reset()
  C.staged = {}
  AL.last, AL.last_reads = nil, nil
  vim.cmd("silent! only")
  vim.cmd("silent! %bwipeout!")
end

--- Open a fixture from a copy in the scratch directory.
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if st.view ~= "linear" then buffer.set_view(buf, "linear") end
  return buf, st, buffer.record(st)
end

local function write_file(name, text)
  local p = tmp .. "/" .. name
  local fh = assert(io.open(p, "wb")); fh:write(text); fh:close()
  return p
end

local function seqof(rec) return rec.seq:tostring():upper() end

--- A copy of `seq` with substitutions at the given positions and, optionally,
--- `del` bases removed at `del_at`.
local function variant(seq, subs, del_at, del_n)
  local t = {}
  for i = 1, #seq do t[i] = seq:sub(i, i) end
  for _, p in ipairs(subs or {}) do
    t[p] = ({ A = "G", G = "A", C = "T", T = "C" })[t[p]:upper()] or "A"
  end
  if del_at then for _ = 1, (del_n or 1) do table.remove(t, del_at) end end
  return table.concat(t)
end

local function view_text(b)
  b = b or AL.last
  H.ok(b and vim.api.nvim_buf_is_valid(b), "no alignment view")
  return table.concat(vim.api.nvim_buf_get_lines(b, 0, -1, false), "\n"), b
end

local function msgs_of(fn)
  local out = {}
  local orig = vim.notify
  vim.notify = function(m) out[#out + 1] = tostring(m) end
  local ok, err = pcall(fn)
  vim.notify = orig
  if not ok then error(err, 0) end
  return out, table.concat(out, "\n")
end

local function focus(b)
  local win = vim.fn.bufwinid(b)
  H.ok(win ~= -1, "buffer is not on screen")
  vim.api.nvim_set_current_win(win)
  return win
end

---------------------------------------------------------------------------
-- :DnaAlign
---------------------------------------------------------------------------

case("align_places_an_inline_read_inside_the_record", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local read = variant(seq, { 150, 300, 420 }, 350, 2):sub(100, 700)
  local vb, aln, diffs = AL.align(buf, { read })
  H.ok(vb ~= nil, "a view was opened")
  H.eq("semi", aln.mode)
  H.eq(1, aln.strand)
  H.eq(100, aln.ref_start)
  H.eq(702, aln.ref_stop, "601 read bases cover 603 reference bases: two were deleted")
  H.eq(4, #diffs, "three substitutions and one deletion")
  local text = view_text(vb)
  H.ok(text:find("Alignment: " .. rec.name, 1, true), text:sub(1, 200))
  H.ok(text:find("differences (4)", 1, true), "the difference list is headed")
  H.ok(text:find("plus strand", 1, true), "the strand is in the header")
  H.ok(text:find("identity 99", 1, true), text:sub(1, 300))
  H.log(("[align] 601 bp read on %d bp: identity %.1f%%, %d differences")
    :format(rec:len(), 100 * aln.identity, #diffs))
end)

case("align_takes_a_file_a_buffer_a_part_and_the_selection", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  -- a FASTA file
  local path = write_file("read_file.fa", ">read_from_file\n" .. seq:sub(1000, 1400) .. "\n")
  local _, aln = AL.align(buf, { path })
  H.eq(1000, aln.ref_start)
  H.eq(1400, aln.ref_stop)
  H.ok(view_text():find("read_from_file", 1, true), "the file's record name is used")
  -- an open buffer, by record name
  local other = write_file("piece.gb", pio.format({ (function()
    local R = require("plasmid.model.record").Record
    return R.new({ name = "piece", seq = seq:sub(2000, 2300), topology = "linear", alphabet = "dna" })
  end)() }, "genbank"))
  vim.cmd("silent edit! " .. vim.fn.fnameescape(other))
  vim.cmd("silent buffer " .. buf)
  local _, baln = AL.align(buf, { "piece" })
  H.eq(2000, baln.ref_start)
  H.eq(2300, baln.ref_stop)
  -- a staged cloning part
  C.staged = {}
  C.stage({ name = "insert", rec = (function()
    local R = require("plasmid.model.record").Record
    return R.new({ name = "insert", seq = seq:sub(3000, 3200), topology = "linear", alphabet = "dna" })
  end)(), from = "test" })
  local _, paln = AL.align(buf, { "insert" })
  H.eq(3000, paln.ref_start)
  -- the selection
  local sel = require("plasmid.select")
  focus(buf)
  sel.set(buf, sel.make(rec, 500, 760))
  local _, saln = AL.align(buf, { "selection" })
  H.eq(500, saln.ref_start)
  H.eq(760, saln.ref_stop)
  H.eq(1, saln.identity)
  sel.set(buf, nil)
end)

case("mode_and_strand_options_are_honoured", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local rc = alphabet.revcomp(seq:sub(1200, 1700))
  local _, auto = AL.align(buf, { rc })
  H.eq(-1, auto.strand, "the better strand is found on its own")
  H.eq(1200, auto.ref_start)
  H.ok(view_text():find("reverse%-complemented"), "the view says the query was flipped")
  local _, plus = AL.align(buf, { rc, "strand=plus" })
  H.eq(1, plus.strand, "strand=plus refuses to flip the query")
  H.ok(plus.identity < 0.7, ("forced to the wrong strand: %.2f identity"):format(plus.identity))
  local _, minus = AL.align(buf, { seq:sub(1200, 1700), "strand=minus" })
  H.eq(-1, minus.strand)
  local _, loc = AL.align(buf, { seq:sub(1200, 1700), "mode=local" })
  H.eq("local", loc.mode)
  local _, glob = AL.align(buf, { seq:sub(1200, 1700), "mode=global" })
  H.eq("global", glob.mode)
  H.eq(1, glob.ref_start, "a global alignment charges the record's end gaps")
  -- name= renames the query in the view
  AL.align(buf, { seq:sub(1200, 1300), "name=sanger_01" })
  H.ok(view_text():find("sanger_01", 1, true), "name= is used as the query's name")
end)

case("a_read_may_be_placed_across_the_origin", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local L = #seq
  local read = seq:sub(L - 300) .. seq:sub(1, 200)
  local _, aln = AL.align(buf, { read })
  H.eq(true, aln.wraps)
  H.eq(L - 300, aln.ref_start)
  H.eq(200, aln.ref_stop)
  local text = view_text()
  H.ok(text:find("across the origin", 1, true), "the header says the read crosses the origin")
  -- the reference coordinates wrap back into 1..L part way through the blocks
  local coords = {}
  for n in text:gmatch("\n " .. rec.name .. "%s+(%d+) ") do coords[#coords + 1] = tonumber(n) end
  H.ok(#coords >= 2, ("%d blocks drawn"):format(#coords))
  H.eq(L - 300, coords[1])
  H.ok(coords[#coords] < coords[1], ("the last block starts at %d, after the origin"):format(coords[#coords]))
end)

case("the_user_commands_drive_the_same_thing", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local before = #rec.features
  local read = write_file("cmd_read.fa", ">cmd_read\n" .. variant(seq, { 2500 }):sub(2300, 2800) .. "\n")
  vim.cmd("DnaAlign " .. vim.fn.fnameescape(read))
  local vb = AL.last
  H.ok(vb ~= nil and vim.api.nvim_buf_is_valid(vb), ":DnaAlign opened a view")
  H.eq(1, #AL.state[vb].diffs)
  -- :DnaAlignAnnotate runs from the view as well as from the record
  focus(vb)
  vim.cmd("DnaAlignAnnotate")
  H.eq(before + 1, #rec.features)
  focus(buf)
  require("plasmid.author").undo()
  H.eq(before, #rec.features)
  vim.cmd("DnaAlignReads " .. vim.fn.fnameescape(read))
  H.ok(AL.last_reads ~= nil, ":DnaAlignReads opened a table")
  H.eq(1, #AL.read_state[AL.last_reads].rows)
  focus(buf)
  local other = write_file("cmd_v2.fa", ">cmd_v2\n" .. variant(seq, { 700 }) .. "\n")
  vim.cmd("DnaCompare " .. vim.fn.fnameescape(other))
  H.ok(view_text(AL.last):find("Comparison: ", 1, true), "DnaCompare opened the comparison")
  -- every command is registered with a description
  local cmds = vim.api.nvim_get_commands({})
  for _, name in ipairs({ "DnaAlign", "DnaAlignReads", "DnaCompare", "DnaAlignAnnotate" }) do
    H.ok(cmds[name] ~= nil and cmds[name].definition ~= "", name .. " is registered")
  end
end)

---------------------------------------------------------------------------
-- The view
---------------------------------------------------------------------------

case("the_view_draws_reference_match_and_query_rows", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local read = variant(seq, { 220 }, nil, nil):sub(200, 400)
  local vb = AL.align(buf, { read, "name=read1" })
  local lines = vim.api.nvim_buf_get_lines(vb, 0, -1, false)
  local first
  for i, l in ipairs(lines) do
    if l:find("^ " .. rec.name .. "%s") then first = i break end
  end
  H.ok(first ~= nil, "a reference row was drawn")
  local refrow, midrow, qrow = lines[first], lines[first + 1], lines[first + 2]
  H.ok(refrow:find(" 200 "), refrow)
  H.ok(qrow:find("^ read1%s+1 "), qrow)
  H.ok(midrow:find("|||"), "the match row marks identical columns")
  H.eq(#refrow:gsub("%s+%d+$", ""), #qrow:gsub("%s+%d+$", ""), "the rows line up")
  -- the mismatch is a space in the match row, under both residues
  local col = midrow:find(" ", midrow:find("|") + 1)
  H.ok(col ~= nil and refrow:sub(col, col):match("%a") and qrow:sub(col, col):match("%a"),
    "the mismatched column carries a residue on both rows")
  -- blocks fit the window
  local width = vim.api.nvim_win_get_width(vim.fn.bufwinid(vb))
  H.ok(#refrow <= width, ("block row %d wide, window %d"):format(#refrow, width))
  H.log(("[align] block rows are %d columns wide in a %d column window"):format(#refrow, width))
end)

case("the_view_jumps_between_differences", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local read = variant(seq, { 1050, 1200, 1400 }):sub(1000, 1500)
  local vb, _, diffs = AL.align(buf, { read, "name=read1" })
  H.eq(3, #diffs)
  local win = focus(vb)
  vim.api.nvim_win_set_cursor(win, { 1, 0 })
  local seen = {}
  for _ = 1, 3 do
    local t = AL.jump(vb, 1)
    seen[#seen + 1] = t.diff
  end
  H.eq({ 1, 2, 3 }, seen, "]d walks the difference list in order")
  -- and on into the blocks, then back
  local into = AL.jump(vb, 1)
  H.eq(1, into.diff, "the first block mark comes after the list")
  local back = AL.jump(vb, -1)
  H.eq(3, back.diff, "[d goes back to the last list row")
  -- the keymaps are set on the buffer
  local maps = {}
  for _, m in ipairs(vim.api.nvim_buf_get_keymap(vb, "n")) do maps[m.lhs] = true end
  for _, lhs in ipairs({ "]d", "[d", "<CR>", "ga", "q" }) do
    H.ok(maps[lhs], ("keymap %s is set on the view"):format(lhs))
  end
end)

case("enter_goes_to_that_difference_in_the_record", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local read = variant(seq, { 1234 }):sub(1100, 1400)
  local vb, _, diffs = AL.align(buf, { read })
  H.eq(1234, diffs[1].pos)
  local win = focus(vb)
  -- put the cursor on the difference's row in the list
  local S = AL.state[vb]
  local row
  for line, idx in pairs(S.list_line) do if idx == 1 then row = line end end
  vim.api.nvim_win_set_cursor(win, { row, 0 })
  local pos = AL.goto_difference(vb)
  H.eq(1234, pos)
  local sel = require("plasmid.select").get(buf)
  H.eq(1234, sel.start)
  H.eq(1234, sel.stop)
  H.eq(buf, vim.api.nvim_get_current_buf(), "the record's window is focused")
  require("plasmid.select").set(buf, nil)
end)

case("the_difference_list_names_the_feature_the_codon_and_the_effect", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")     -- tet 86..1276 (+), bla complement(3293..4153)
  local seq = seqof(rec)
  local query = variant(seq, { 150, 3400 })
  local vb, _, diffs = AL.align(buf, { query, "mode=global", "name=resequenced" })
  H.eq(2, #diffs)
  local text = view_text(vb)
  H.ok(text:find("150 ", 1, true), text:sub(1, 400))
  H.ok(text:find("tet codon 22", 1, true), "a tet codon is named")
  H.ok(text:find("bla codon 252", 1, true), "the reverse-strand CDS is numbered from its own start")
  H.ok(text:find("missense", 1, true), "the effect is given")
  H.ok(text:find("features hit (", 1, true), "the features section is drawn")
  H.ok(text:find("2 substitutions, 0 insertions, 0 deletions", 1, true), "the summary counts the kinds")
  H.ok(text:find("0 silent, 2 missense, 0 nonsense", 1, true), "the summary counts the effects")
end)

case("a_frameshift_is_called_out_in_the_features_section", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local query = variant(seq, {}, 600, 1)       -- one base out of tet
  local vb, _, diffs = AL.align(buf, { query, "mode=global" })
  H.eq(1, #diffs)
  H.eq("del", diffs[1].kind)
  local text = view_text(vb)
  H.ok(text:find("frameshift", 1, true), text:match("features hit.-\n.-\n.-\n") or text:sub(1, 400))
  H.ok(text:find("1 frameshifted feature", 1, true), "the summary counts frameshifted features")
  H.ok(text:find("stop codon at aa", 1, true), "a premature stop is reported")
end)

---------------------------------------------------------------------------
-- :DnaAlignReads
---------------------------------------------------------------------------

local function reads_fasta(seq)
  local parts = {
    ">read_mid\n" .. variant(seq, { 2050 }):sub(2000, 2600),
    ">read_start\n" .. seq:sub(100, 700),
    ">read_rc\n" .. alphabet.revcomp(seq:sub(3000, 3600)),
    ">read_junk\n" .. ("ACGT"):rep(100),
  }
  return table.concat(parts, "\n") .. "\n"
end

case("reads_are_placed_in_a_table_sorted_by_position", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local path = write_file("reads.fa", reads_fasta(seqof(rec)))
  local tb, rows, unplaced = AL.align_reads(buf, { path })
  H.ok(tb ~= nil, "a read table was opened")
  H.eq(3, #rows, "three of the four reads place")
  H.eq(1, #unplaced)
  local names = {}
  for i, r in ipairs(rows) do names[i] = r.name end
  H.eq({ "read_start", "read_mid", "read_rc" }, names, "rows are in reference order")
  H.eq(100, rows[1].aln.ref_start)
  H.eq(-1, rows[3].aln.strand, "the reverse read is found on the minus strand")
  H.eq(1, #rows[2].diffs)
  local text = table.concat(vim.api.nvim_buf_get_lines(tb, 0, -1, false), "\n")
  H.ok(text:find("3 of 4 placed", 1, true), text:sub(1, 200))
  H.ok(text:find("read ", 1, true) and text:find("strand", 1, true) and text:find("identity", 1, true),
    "the table has a header row")
  H.ok(text:find("coverage", 1, true), "coverage is a column")
  H.ok(text:find("not placed", 1, true), "unplaced reads get their own section")
  H.ok(text:find("read_junk", 1, true) and text:find("align.min_identity", 1, true),
    "the unplaced read is named with the reason")
  H.log(("[reads] %d of 4 placed; unplaced: %s"):format(#rows, unplaced[1].why))
end)

case("choosing_a_read_opens_its_alignment", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local path = write_file("reads2.fa", reads_fasta(seqof(rec)))
  local tb = AL.align_reads(buf, { path })
  local S = AL.read_state[tb]
  local win = focus(tb)
  vim.api.nvim_win_set_cursor(win, { S.first + 1, 0 })     -- the second row: read_mid
  local vb, row = AL.open_read(tb)
  H.eq("read_mid", row.name)
  local text = view_text(vb)
  H.ok(text:find("read_mid", 1, true), "the view is the read's own alignment")
  H.ok(text:find("differences (1)", 1, true), text:sub(1, 200))
  -- a row that does not exist says so
  local out = msgs_of(function()
    vim.api.nvim_win_set_cursor(vim.fn.bufwinid(tb), { vim.api.nvim_buf_line_count(tb), 0 })
    AL.open_read(tb)
  end)
  H.ok(table.concat(out, "\n"):find("cursor on a read row", 1, true), vim.inspect(out))
end)

case("a_glob_and_several_files_are_both_accepted", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local dir = tmp .. "/sanger"
  vim.fn.mkdir(dir, "p")
  for i, from in ipairs({ 400, 1600, 2800 }) do
    local fh = assert(io.open(("%s/well_%d.fa"):format(dir, i), "wb"))
    fh:write((">well_%d\n%s\n"):format(i, seq:sub(from, from + 500)))
    fh:close()
  end
  local _, rows = AL.align_reads(buf, { dir .. "/*.fa" })
  H.eq(3, #rows, "the glob picked up every file")
  H.eq(400, rows[1].aln.ref_start)
  local _, two = AL.align_reads(buf, { dir .. "/well_1.fa", dir .. "/well_3.fa" })
  H.eq(2, #two)
  H.eq({ 400, 2800 }, { two[1].aln.ref_start, two[2].aln.ref_start })
end)

---------------------------------------------------------------------------
-- :DnaCompare
---------------------------------------------------------------------------

case("compare_reports_the_rotation_the_strand_and_the_differences", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local rot = 1234
  local rotated = seq:sub(rot + 1) .. seq:sub(1, rot)
  local changed = variant(rotated, { 500, 1500 }, 3000, 1)
  local path = write_file("v2.fa", ">pBR322_v2\n" .. changed .. "\n")
  local vb, aln, diffs = AL.compare(buf, { path })
  H.eq(1, aln.strand)
  H.eq(3, #diffs, "two substitutions and one deletion, and no origin artefact")
  for i = 2, #diffs do H.ok(diffs[i].pos > diffs[i - 1].pos, "the differences are in order") end
  local text = view_text(vb)
  H.ok(text:find("Comparison: ", 1, true), text:sub(1, 200))
  H.ok(text:find("pBR322_v2", 1, true), "both record names are in the title")
  H.ok(text:find("rotated ", 1, true), "the rotation that was applied is reported")
  -- the differences are in the record's own coordinates
  H.eq(seq:sub(diffs[1].pos, diffs[1].pos), diffs[1].ref)
  -- a reverse-complemented copy is recognised as the same molecule
  local rcpath = write_file("v3.fa", ">pBR322_rc\n" .. alphabet.revcomp(changed) .. "\n")
  local _, rcaln, rcdiffs = AL.compare(buf, { rcpath })
  H.eq(-1, rcaln.strand)
  H.eq(3, #rcdiffs)
  H.log(("[compare] %d bp vs %d bp: rotation %d, %d differences")
    :format(rec:len(), #changed, aln.rotation or 0, #diffs))
end)

---------------------------------------------------------------------------
-- Annotating
---------------------------------------------------------------------------

case("annotating_writes_variation_features_in_one_undo_step", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local before = #rec.features
  local query = variant(seq, { 150, 3400 }, 2000, 2)
  local vb, _, diffs = AL.align(buf, { query, "mode=global", "name=read_A" })
  H.eq(3, #diffs)
  local n = AL.annotate(vb)
  H.eq(3, n)
  H.eq(before + 3, #rec.features)
  local added = {}
  for _, f in ipairs(rec.features) do
    if f.type == "variation" then added[#added + 1] = f end
  end
  H.eq(3, #added)
  local by_label = {}
  for _, f in ipairs(added) do by_label[f:label()] = f end
  local sub = by_label[("%s150%s"):format(seq:sub(150, 150), diffs[1].query)]
  H.ok(sub ~= nil, vim.inspect(vim.tbl_keys(by_label)))
  H.eq("150..150", sub.location_text)
  H.ok(sub:get("note"):find("codon 22", 1, true), sub:get("note"))
  H.ok(sub:get("note"):find("read_A", 1, true), "the query is named in the note")
  local del = by_label["del 2 bp @2000"]
  H.ok(del ~= nil, "the deletion was annotated")
  H.eq("2000..2001", del.location_text)
  -- one undo takes all three back (undo works on the record's own buffer)
  focus(buf)
  require("plasmid.author").undo()
  H.eq(before, #rec.features)
  require("plasmid.author").redo()
  H.eq(before + 3, #rec.features)
  -- :DnaAlignAnnotate from the record's own buffer finds the same alignment
  require("plasmid.author").undo()
  H.eq(before, #rec.features)
  H.eq(3, AL.annotate(buf))
  H.eq(before + 3, #rec.features)
end)

case("an_insertion_is_annotated_between_two_bases", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local query = seq:sub(1, 1500) .. "GGGG" .. seq:sub(1501)
  local vb, _, diffs = AL.align(buf, { query, "mode=global" })
  H.eq(1, #diffs)
  H.eq("ins", diffs[1].kind)
  AL.annotate(vb)
  local ins
  for _, f in ipairs(rec.features) do if f.type == "variation" then ins = f end end
  H.ok(ins ~= nil, "the insertion was annotated")
  H.eq(("%d^%d"):format(diffs[1].pos, diffs[1].pos + 1), ins.location_text)
  H.ok(ins:get("note"):find("ins 4 bp", 1, true), ins:get("note"))
end)

case("annotated_differences_survive_a_genbank_round_trip", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  local vb = AL.align(buf, { variant(seq, { 150 }, 2000, 3), "mode=global", "name=read_B" })
  H.eq(2, AL.annotate(vb))
  local text = assert(pio.format({ rec }, "genbank"))
  H.ok(text:find("     variation", 1, true), "variation features are written")
  local back = assert(pio.parse(text, "round.gb", "genbank")).records[1]
  local added = {}
  for _, f in ipairs(back.features) do if f.type == "variation" then added[#added + 1] = f end end
  H.eq(2, #added, "both come back")
  local locs, notes = {}, {}
  for _, f in ipairs(added) do
    locs[#locs + 1] = f.location_text
    notes[#notes + 1] = f:get("note")
  end
  table.sort(locs)
  H.eq({ "150..150", "2000..2002" }, locs)
  H.ok(table.concat(notes, " "):find("read_B", 1, true), vim.inspect(notes))
  H.eq(seq, seqof(back), "annotating does not touch the sequence")
end)

case("the_annotate_type_follows_the_config", function()
  reset()
  local config = require("plasmid.config")
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  config.options.align.annotate_type = "misc_difference"
  local vb = AL.align(buf, { variant(seq, { 900 }), "mode=global" })
  AL.annotate(vb)
  local found
  for _, f in ipairs(rec.features) do if f.type == "misc_difference" then found = f end end
  H.ok(found ~= nil, "the configured type is used")
  config.options.align.annotate_type = "variation"
end)

---------------------------------------------------------------------------
-- Proteins
---------------------------------------------------------------------------

case("protein_records_align_too", function()
  reset()
  local recs = assert(pio.read(H.fixture("protein/proteins.faa"))).records
  local rec
  for _, r in ipairs(recs) do if (r.id or ""):find("HBA_HUMAN") then rec = r end end
  rec = rec or recs[1]
  local p = rec.seq:tostring():upper()
  local path = write_file("hba.faa", (">%s\n%s\n"):format(rec.name, p))
  vim.cmd("silent edit! " .. vim.fn.fnameescape(path))
  local buf = vim.api.nvim_get_current_buf()
  local prec = buffer.record(buffer.get(buf))
  H.eq(true, prec:is_protein())
  -- one residue changed, in the middle
  local mutant = p:sub(1, 29) .. (p:sub(30, 30) == "W" and "A" or "W") .. p:sub(31)
  local vb, aln, diffs = AL.align(buf, { mutant, "mode=global", "name=mutant" })
  H.eq(1, aln.strand, "a protein has no strands")
  H.eq(1, #diffs)
  H.eq(30, diffs[1].pos)
  local text = view_text(vb)
  H.ok(not text:find("strand", 1, true), "the header does not talk about strands")
  H.ok(text:find(" aa ", 1, true) or text:find("%d+ aa"), "residues are counted in aa")
  H.ok(text:find("mutant", 1, true), "the query is named")
  H.log(("[align] protein %d aa: score %d, identity %.1f%%"):format(#p, aln.score, 100 * aln.identity))
end)

---------------------------------------------------------------------------
-- Failures
---------------------------------------------------------------------------

case("failure_messages_name_what_went_wrong", function()
  reset()
  local buf, _, rec = open("ncbi/J01749.gb")
  -- no query at all
  local _, m1 = msgs_of(function() AL.align(buf, {}) end)
  H.ok(m1:find("DnaAlign {sequence", 1, true), m1)
  -- a query that is nowhere to be found
  local _, m2 = msgs_of(function() AL.align(buf, { "no_such_thing.gb" }) end)
  H.ok(m2:find("no part called", 1, true), m2)
  -- nothing is selected
  require("plasmid.select").set(buf, nil)
  local _, m3 = msgs_of(function() AL.align(buf, { "selection" }) end)
  H.ok(m3:find("nothing is selected", 1, true), m3)
  -- nothing aligns
  local _, m4 = msgs_of(function() AL.align(buf, { ("N"):rep(60), "mode=local" }) end)
  H.ok(m4:find("does not align", 1, true), m4)
  -- an option that does not exist, and a bad value
  local _, m5 = msgs_of(function() AL.align(buf, { "ACGTACGTACGT", "fuzzy=yes" }) end)
  H.ok(m5:find("no such option", 1, true), m5)
  local _, m6 = msgs_of(function() AL.align(buf, { "ACGTACGTACGT", "mode=fuzzy" }) end)
  H.ok(m6:find("mode must be", 1, true), m6)
  local _, m7 = msgs_of(function() AL.align(buf, { "ACGTACGTACGT", "strand=sideways" }) end)
  H.ok(m7:find("strand must be", 1, true), m7)
  -- reads from a file that is not there
  local _, m8 = msgs_of(function() AL.align_reads(buf, { tmp .. "/absent_*.fa" }) end)
  H.ok(m8:find("no reads in", 1, true), m8)
  -- annotating with nothing aligned yet
  reset()
  local buf2 = open("ncbi/J01749.gb")
  local _, m9 = msgs_of(function() AL.annotate(buf2) end)
  H.ok(m9:find("no alignment to annotate", 1, true), m9)
  -- a command run where there is no record
  vim.cmd("enew")
  local _, m10 = msgs_of(function() AL.compare(vim.api.nvim_get_current_buf(), { "x" }) end)
  H.ok(m10:find("not a sequence file", 1, true), m10)
  local _ = rec
end)

case("a_long_alignment_shows_only_the_regions_around_differences", function()
  reset()
  local config = require("plasmid.config")
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  config.options.align.full_limit = 500
  config.options.align.context = 20
  local vb, _, diffs = AL.align(buf, { variant(seq, { 1000, 3000 }), "mode=global" })
  H.eq(2, #diffs)
  local text = view_text(vb)
  H.ok(text:find("only the 2 regions", 1, true), text:match("%(%d+ columns[^\n]*") or text:sub(-300))
  H.ok(text:find("⋮", 1, true), "the skipped stretch is marked")
  local blocks = select(2, text:gsub("\n " .. rec.name .. "%s+%d+ ", ""))
  H.ok(blocks <= 4, ("%d blocks drawn instead of the whole 4.4 kb"):format(blocks))
  config.options.align.full_limit = 20000
  config.options.align.context = 30
  -- the whole alignment again
  local vb2 = AL.align(buf, { variant(seq, { 1000, 3000 }), "mode=global" })
  local full = select(2, view_text(vb2):gsub("\n " .. rec.name .. "%s+%d+ ", ""))
  H.ok(full > 10, ("the whole alignment is %d blocks"):format(full))
end)

case("the_block_width_follows_the_config", function()
  reset()
  local config = require("plasmid.config")
  local buf, _, rec = open("ncbi/J01749.gb")
  local seq = seqof(rec)
  config.options.align.width = 40
  local vb = AL.align(buf, { seq:sub(1000, 1300), "name=w" })
  local residues
  for _, l in ipairs(vim.api.nvim_buf_get_lines(vb, 0, -1, false)) do
    local body = l:match("^ " .. rec.name .. "%s+%d+ ([ACGTN%-]+) %d+$")
    if body then residues = #body break end
  end
  H.eq(40, residues)
  config.options.align.width = 0
end)

return T
