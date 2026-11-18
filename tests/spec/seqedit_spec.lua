local H = require("helpers")
local translate = require("plasmid.bio.translate")
local alphabet = require("plasmid.model.alphabet")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local config = require("plasmid.config")
local se = require("plasmid.seqedit")
local tx = require("plasmid.model.transaction")

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
local function entry(rep, label, ftype)
  for _, e in ipairs(rep.entries) do
    if e.label == label and (not ftype or e.type == ftype) then return e end
  end
end
local function cds(rec, gene)
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:get("gene") == gene then return f end end
end

vim.o.columns, vim.o.lines = 120, 50

case("in_frame_insertion_and_frameshift", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  -- 500 = 86 + 3*138: the start of codon 139 of tet
  local rep = assert(se.insert(buf, 500, "GCT"))
  local e = assert(entry(rep, "tet", "CDS"), "tet in report")
  H.eq("extended", e.kind)
  H.eq(1191, e.old_len)
  H.eq(1194, e.new_len)
  H.eq(false, e.frameshift)
  H.eq(e.old_aa + 1, e.new_aa, "one codon inserted")
  H.eq(nil, e.premature_stop)
  local tet = cds(rec, "tet")
  H.eq("86..1279", tet.location_text)
  H.eq("A", translate.translate(rec:extract(tet), { table = 11 }):sub(139, 139), "Ala at codon 139")
  -- the gene feature moved with it; nothing downstream lost bases
  for _, x in ipairs(rep.entries) do H.ok(x.label == "tet", "only tet affected, got " .. x.label) end
  -- one base: frameshift
  local rep2 = assert(se.insert(buf, 600, "A"))
  H.eq(true, entry(rep2, "tet", "CDS").frameshift)
  H.ok(require("plasmid.bio.impact").alarming(rep2))
end)

case("substitutions_stop_and_silent", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local old = rec.seq:sub(500, 502):upper()
  local rep = assert(se.replace(buf, 500, 502, "TAA"))
  local e = entry(rep, "tet", "CDS")
  H.eq("modified", e.kind)
  H.eq(139, e.premature_stop)
  H.eq("86..1276", cds(rec, "tet").location_text, "same length: coordinates unchanged")
  -- back to the original codon, then a synonymous third-position change
  se.replace(buf, 500, 502, old)
  local t = translate.get_table(11)
  local aa = translate.codon(t, old)
  local syn
  for _, b in ipairs({ "A", "C", "G", "T" }) do
    local c = old:sub(1, 2) .. b
    if c ~= old and translate.codon(t, c) == aa then syn = b end
  end
  if syn then
    local rep3 = assert(se.replace(buf, 502, 502, syn))
    local e3 = entry(rep3, "tet", "CDS")
    H.eq(0, e3.aa_changes, "silent")
    H.ok(require("plasmid.bio.impact").describe_entry(e3):find("silent", 1, true))
  end
end)

case("start_codon_loss_and_destruction", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  -- bla is complement(3293..4153): its ATG is 4151..4153 on the bottom strand
  local rep = assert(se.delete(buf, 4151, 4153))
  local e = entry(rep, "bla", "CDS")
  H.eq("truncated", e.kind)
  H.eq(true, e.lost_start)
  H.eq(false, e.frameshift)
  H.eq("complement(3293..4150)", cds(rec, "bla").location_text)
  -- ROP 1915..2106 wholly deleted
  local rep2 = assert(se.delete(buf, 1900, 2120))
  local rop
  for _, x in ipairs(rep2.entries) do if x.label:find("ROP") or x.type == "CDS" then rop = rop or x end end
  H.ok(rop and rop.kind == "destroyed", "ROP destroyed")
  for _, f in ipairs(rec.features) do H.ok(f:get("product") ~= "ROP protein", "ROP feature removed") end
  -- one undo per edit
  local f0 = #rec.features
  require("plasmid.author").undo()
  H.ok(#rec.features > f0, "undo restores ROP")
  require("plasmid.author").undo()
  H.eq("complement(3293..4153)", cds(rec, "bla").location_text)
end)

case("cut_policies", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local before = rec.seq:tostring()
  local depth = #rec.history.undo
  config.options.edit.cut_policy = "refuse"
  H.messages = {}
  local rep, why = se.delete(buf, 500, 502)
  H.eq(nil, rep)
  H.ok(why:find("refused", 1, true), why)
  H.eq(before, rec.seq:tostring(), "nothing changed")
  H.eq(depth, #rec.history.undo, "nothing recorded")
  -- an edit that cuts nothing is allowed
  H.ok(se.insert(buf, 1300, "A"), "insertion outside features allowed")
  -- ! forces it
  H.ok(se.delete(buf, 500, 502, { force = true }))
  H.eq("86..1273", cds(rec, "tet").location_text)
  require("plasmid.author").undo()
  require("plasmid.author").undo()
  H.eq(before, rec.seq:tostring())
  -- flag: the truncated feature gets a note, removed again by one undo
  config.options.edit.cut_policy = "flag"
  H.ok(se.delete(buf, 500, 502))
  local notes = table.concat(cds(rec, "tet"):get_all("note"), " | ")
  H.ok(notes:find("truncated by a sequence edit: 3 of 1191 bases removed", 1, true), notes)
  require("plasmid.author").undo()
  H.ok(not table.concat(cds(rec, "tet"):get_all("note"), " "):find("truncated by", 1, true), "note gone")
  H.eq(before, rec.seq:tostring())
  config.options.edit.cut_policy = "truncate"
end)

case("validation_and_case", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  H.eq(nil, (se.clean(rec, "ACGX")), "X is not a nucleotide")
  H.eq(nil, (se.clean(rec, "  12 ")), "nothing left")
  -- GenBank ORIGIN text pastes cleanly; case follows the stored (lower-case) record
  H.eq("acgtnnrytt", se.clean(rec, "      1 ACGTNN RYTT\n"))
  local fa = tmp .. "/rna.fa"
  local fh = assert(io.open(fa, "w")); fh:write(">r\nACGUACGUacgu\n"); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(fa))
  local frec = buffer.record(buffer.get(0))
  H.eq("UUAcc", se.clean(frec, "TTAcc"), "RNA record: T becomes U; FASTA keeps the typed case")
  local prot = tmp .. "/p.faa"
  fh = assert(io.open(prot, "w")); fh:write(">p\nMKVLAAGIVG\n"); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(prot))
  local prec = buffer.record(buffer.get(0))
  H.eq("WY*", se.clean(prec, "WY*"))
  H.eq(nil, (se.clean(prec, "W1Y2@")), "@ is not an amino acid")
end)

case("wrapped_delete_is_one_step_and_clears_proposals", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local before = rec.seq:tostring()
  local L = rec:len()
  st.detect = { hits = {}, rec = rec } -- stand-in for :DnaDetect proposals
  local depth = #rec.history.undo
  H.ok(se.delete(buf, L - 9, 5))
  H.eq(L - 15, rec:len())
  H.eq(before:sub(6, L - 10), rec.seq:tostring(), "both ends removed")
  H.eq(depth + 1, #rec.history.undo, "one undo step")
  H.eq(nil, st.detect, "proposals refer to old coordinates and are cleared")
  require("plasmid.author").undo()
  H.eq(before, rec.seq:tostring())
  H.ok(buffer.get(buf).modified == false or vim.bo[buf].modified == false, "back to the saved state")
end)

case("report_float", function()
  local buf = open("ncbi/J01749.gb")
  H.eq(nil, se.report_lines(buf))
  se.insert(buf, 600, "A")
  local lines = se.report_lines(buf)
  H.ok(table.concat(lines, "\n"):find("FRAMESHIFT", 1, true), table.concat(lines, "\n"))
  local fbuf, win = se.show_report(buf)
  H.ok(vim.api.nvim_win_is_valid(win))
  vim.api.nvim_win_close(win, true)
end)

return T
