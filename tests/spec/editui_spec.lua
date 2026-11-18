local H = require("helpers")
local alphabet = require("plasmid.model.alphabet")
local translate = require("plasmid.bio.translate")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local codon = require("plasmid.codon")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local n_open = 0
local function open(rel)
  n_open = n_open + 1
  local dst = tmp .. "/" .. n_open .. "_" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.get(buf), buffer.record(buffer.get(buf))
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end
local function flush() vim.wait(30, function() return false end) end
local function answer(text, fn)
  local orig = vim.ui.input
  local prompt
  vim.ui.input = function(opts, cb) prompt = opts.prompt; cb(text) end
  local ok, err = pcall(fn)
  vim.ui.input = orig
  if not ok then error(err, 0) end
  return prompt
end
local function cds(rec, gene)
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:get("gene") == gene then return f end end
end
local function at() return (buffer.cursor_info(0)) end

vim.o.columns, vim.o.lines = 120, 50

case("insert_keys_and_undo", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  buffer.goto_bp(0, 500)
  local prompt = answer("GCT", function() keys("i") end)
  H.eq("Insert before 500: ", prompt)
  H.eq(orig:sub(1, 499) .. "gct" .. orig:sub(500), rec.seq:tostring(), "inserted before 500 in the record's case")
  H.eq("86..1279", cds(rec, "tet").location_text)
  H.eq(502, at(), "cursor on the last inserted base")
  H.eq(true, vim.bo[buf].modified)
  H.ok(buffer.winbar():find("unsaved", 1, true))
  answer("tt", function() keys("a") end)
  H.eq(orig:sub(1, 499) .. "gcttt" .. orig:sub(500), rec.seq:tostring(), "a inserts after the cursor base")
  keys("u")
  keys("u")
  H.eq(orig, rec.seq:tostring())
  H.eq(false, vim.bo[buf].modified, "back to the saved state")
  keys("<C-r>")
  H.eq(orig:sub(1, 499) .. "gct" .. orig:sub(500), rec.seq:tostring(), "redo")
  -- invalid letters change nothing
  H.messages = {}
  local before = rec.seq:tostring()
  answer("XYZ", function() keys("i") end)
  H.eq(before, rec.seq:tostring())
  H.ok(H.messages[#H.messages].msg:find("not a valid", 1, true))
end)

case("delete_replace_paste_keys", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  buffer.goto_bp(0, 500)
  keys("3x")
  H.eq(orig:sub(1, 499) .. orig:sub(503), rec.seq:tostring(), "3x deletes three bases")
  H.eq(orig:sub(500, 502):upper(), vim.fn.getreg('"'), "deleted bases in the unnamed register")
  H.eq(500, at())
  keys("P")
  H.eq(orig, rec.seq:tostring(), "P puts them back before the cursor")
  -- xp swaps two bases
  buffer.goto_bp(0, 700)
  keys("x")
  keys("p")
  H.eq(orig:sub(1, 699) .. orig:sub(701, 701) .. orig:sub(700, 700) .. orig:sub(702), rec.seq:tostring(), "xp swaps")
  keys("u"); keys("u")
  H.eq(orig, rec.seq:tostring())
  -- r and [count]r
  buffer.goto_bp(0, 800)
  keys("rG")
  H.eq("g", rec.seq:sub(800, 800))
  keys("3rA")
  H.eq("aaa", rec.seq:sub(800, 802))
  H.eq(#orig, rec:len(), "replacing keeps the length")
  keys("u"); keys("u")
  H.eq(orig, rec.seq:tostring())
end)

case("visual_edits", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  buffer.goto_bp(0, 1301)
  keys("v4ld")
  H.eq(orig:sub(1, 1300) .. orig:sub(1306), rec.seq:tostring(), "visual d deletes the selection")
  keys("u")
  buffer.goto_bp(0, 1301)
  answer("ACGTACGT", function() keys("v4lc") end)
  H.eq(orig:sub(1, 1300) .. "acgtacgt" .. orig:sub(1306), rec.seq:tostring(), "visual c replaces")
  H.eq(1301, st.selection.start)
  H.eq(1308, st.selection.stop, "the new bases are selected")
  keys("u")
  vim.fn.setreg("a", "TTTT")
  buffer.goto_bp(0, 1301)
  keys('v1l"ap')
  H.eq(orig:sub(1, 1300) .. "tttt" .. orig:sub(1303), rec.seq:tostring(), "visual p replaces with a register")
  keys("u")
  buffer.goto_bp(0, 1301)
  keys("v9lR")
  H.eq(orig:sub(1, 1300) .. alphabet.revcomp(orig:sub(1301, 1310)) .. orig:sub(1311), rec.seq:tostring(), "visual R")
  keys("u")
  H.eq(orig, rec.seq:tostring())
end)

case("feature_text_object_edits", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  local rop
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:get("product") == "ROP protein" then rop = f end end
  require("plasmid.ui").goto_feature(buf, rop.id)
  keys("daf")
  flush()
  H.eq(#orig - 192, rec:len(), "daf deletes ROP's span (1915..2106)")
  H.eq(orig:sub(1915, 2106):upper(), vim.fn.getreg('"'))
  for _, f in ipairs(rec.features) do H.ok(f:get("product") ~= "ROP protein", "ROP gone") end
  keys("u")
  H.eq(orig, rec.seq:tostring())
  -- caf replaces the span
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:get("product") == "ROP protein" then rop = f end end
  require("plasmid.ui").goto_feature(buf, rop.id)
  answer("ATGTAA", function() keys("caf"); flush() end)
  H.eq(orig:sub(1, 1914) .. "atgtaa" .. orig:sub(2107), rec.seq:tostring())
  keys("u")
  -- dif on a spliced feature deletes only its parts
  local hbuf, hst, hrec = open("ncbi/U01317.gb")
  local spliced
  for _, f in ipairs(hrec.features) do if f.type == "CDS" and #f.parts == 3 then spliced = f break end end
  local horig = hrec.seq:tostring()
  local exons = hrec:extract(spliced)
  require("plasmid.ui").goto_feature(hbuf, spliced.id)
  keys("dif")
  flush()
  H.eq(#horig - #exons, hrec:len(), "only the exons are removed")
  local p = spliced.parts
  H.eq(horig:sub(p[1].stop + 1, p[2].start - 1), hrec.seq:sub(p[1].start, p[1].start + (p[2].start - p[1].stop - 2)),
    "the first intron now starts where exon 1 was")
  keys("u")
  H.eq(horig, hrec.seq:tostring(), "one undo for the whole feature")
end)

case("commands", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  local L = rec:len()
  vim.cmd("DnaInsert 1 ACGT")
  H.eq("acgt" .. orig, rec.seq:tostring())
  vim.cmd("DnaDelete 1..4")
  H.eq(orig, rec.seq:tostring())
  vim.cmd("DnaReplace 10..12 GGG")
  H.eq("ggg", rec.seq:sub(10, 12))
  vim.cmd("DnaSelect 20..29")
  vim.cmd("DnaReplace TT")
  H.eq(orig:sub(1, 9) .. "ggg" .. orig:sub(13, 19) .. "tt" .. orig:sub(30), rec.seq:tostring(), "selection as the range")
  vim.cmd("DnaRevComp 100..200")
  -- positions are current ones: the earlier edits shortened the record by 8
  H.eq(alphabet.revcomp(orig:sub(108, 208)), rec.seq:sub(100, 200), "revcomp of a range")
  while require("plasmid.model.transaction").undo(rec) do end
  buffer.render(buf)
  H.eq(orig, rec.seq:tostring())
  vim.cmd("DnaRevComp all")
  H.eq(alphabet.revcomp(orig), rec.seq:tostring())
  H.eq(rec:extract(cds(rec, "tet")):upper(), orig:sub(86, 1276):upper(), "tet keeps its sequence")
  vim.cmd("DnaRevComp all")
  H.eq(orig, rec.seq:tostring())
  vim.cmd("DnaRotate 1000")
  H.eq(orig:sub(1000) .. orig:sub(1, 999), rec.seq:tostring())
  H.eq("join(3448..4361,1..277)", cds(rec, "tet").location_text, "tet now crosses the origin")
  H.eq(1, at(), "cursor at the new base 1")
  -- errors
  H.messages = {}
  vim.cmd("DnaDelete 0..5")
  vim.cmd("DnaInsert 99999 A")
  H.eq(2, #H.messages)
  H.eq(L, rec:len())
end)

case("view_stays_consistent", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  vim.cmd("DnaOrfs 60")
  local n0 = #rec.overlays
  H.ok(n0 > 0)
  vim.cmd("DnaInsert 500 GCT")
  -- ORFs recomputed on the new sequence: every one translates without an internal stop
  for _, o in ipairs(rec.overlays) do
    local aa = translate.translate(rec:extract(o), { table = 11 })
    H.eq(nil, aa:sub(1, -2):find("*", 1, true), "ORF " .. o.location_text .. " has no internal stop")
    H.ok(select(2, o:bounds()) <= rec:len())
  end
  -- codon highlight on the edited CDS matches the new sequence
  local lay = st.layout
  local r, v = lay:pos_of(501)
  local row
  for k = r + 1, r + 10 do if lay.rows[k] and lay.rows[k].kind == "aa" and lay.lines[k + 1]:sub(v + 1, v + 1):match("%u") then row = k break end end
  H.ok(row, "inserted codon drawn")
  vim.api.nvim_win_set_cursor(0, { row + 1, v })
  codon.update(buf, 0)
  H.eq("GCT", st.codon.codon)
  H.eq("A", st.codon.aa)
  H.eq(139, st.codon.n)
  -- the impact report is available
  vim.cmd("DnaEditReport")
  local fl = vim.api.nvim_get_current_buf()
  H.ok(table.concat(vim.api.nvim_buf_get_lines(fl, 0, -1, false), "\n"):find("tet  CDS  extended +3 bp", 1, true))
  vim.cmd("close")
  -- the view buffer itself cannot be modified by stray keys
  H.eq(false, vim.bo[buf].modifiable)
end)

case("protein_records", function()
  local p = tmp .. "/p.faa"
  local fh = assert(io.open(p, "w")); fh:write(">p\nMKVLAAGIVG\n"); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(p))
  local buf = vim.api.nvim_get_current_buf()
  local rec = buffer.record(buffer.get(buf))
  vim.cmd("DnaInsert 3 WW")
  H.eq("MKWWVLAAGIVG", rec.seq:tostring())
  H.messages = {}
  vim.cmd("DnaRevComp all")
  H.eq("MKWWVLAAGIVG", rec.seq:tostring(), "no reverse complement for proteins")
end)

return T
