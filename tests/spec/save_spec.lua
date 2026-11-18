local H = require("helpers")
local pio = require("plasmid.io")
local translate = require("plasmid.bio.translate")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local se = require("plasmid.seqedit")

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
  return buf, buffer.get(buf), buffer.record(buffer.get(buf)), dst
end
local function cds(rec, gene)
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:get("gene") == gene then return f end end
end
local function last_msg() return H.messages[#H.messages] and H.messages[#H.messages].msg or "" end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("translations_follow_edits", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local tet0 = cds(rec, "tet")
  local prot0 = tet0:get("translation")
  -- in frame: Ala inserted as residue 139
  local rep = assert(se.insert(buf, 500, "GCT"))
  local tet = cds(rec, "tet")
  H.eq(prot0:sub(1, 138) .. "A" .. prot0:sub(139), tet:get("translation"))
  local e
  for _, x in ipairs(rep.entries) do if x.type == "CDS" then e = x end end
  H.eq("updated", e.translation)
  -- the regenerated value is what translating the CDS gives
  H.eq((translate.translate(rec:extract(tet), { table = 11, cds_start = true }):gsub("%*$", "")), tet:get("translation"))
  -- a frameshift that runs into a stop removes the stale protein
  se.insert(buf, 600, "A")
  tet = cds(rec, "tet")
  local prot = se.translation_of(rec, tet)
  if prot then
    H.eq(prot, tet:get("translation"))
  else
    H.eq(nil, tet:get("translation"), "/translation removed")
    H.ok(table.concat(require("plasmid.seqedit").report_lines(buf), "\n"):find("/translation removed", 1, true))
  end
  -- undo brings back the original feature object and protein
  require("plasmid.author").undo()
  require("plasmid.author").undo()
  H.eq(tet0, cds(rec, "tet"), "original object")
  H.eq(prot0, cds(rec, "tet"):get("translation"))
  -- the option turns it off
  require("plasmid.config").options.edit.update_translations = false
  se.insert(buf, 500, "GCT")
  H.eq(prot0, cds(rec, "tet"):get("translation"), "left as it was")
  require("plasmid.config").options.edit.update_translations = true
end)

case("edited_genbank_writes_and_reads_back", function()
  local buf, st, rec, path = open("ncbi/J01749.gb")
  local text0 = H.read(path)
  se.insert(buf, 500, "GCT")         -- inside tet (translation updated)
  se.delete(buf, 3000, 3010)          -- between features / inside some
  se.replace(buf, 20, 22, "GGGG")     -- near the start: everything downstream shifts
  vim.cmd("write")
  H.eq(false, vim.bo[buf].modified)
  local text = H.read(path)
  local back = assert(pio.read(path)).records[1]
  H.eq(rec.seq:tostring(), back.seq:tostring(), "sequence")
  H.ok(text:find(("LOCUS       SYNPBR322               %d bp"):format(rec:len()), 1, true) or text:find(tostring(rec:len()) .. " bp", 1, true), "LOCUS length updated")
  H.eq(nil, text:find("BASE COUNT", 1, true) and text0:find("BASE COUNT", 1, true) and "stale BASE COUNT" or nil)
  H.eq(#rec.features, #back.features)
  for i, f in ipairs(rec.features) do
    local g = back.features[i]
    H.eq(f.location_text, g.location_text, "location " .. i)
    H.eq(rec:extract(f), back:extract(g), "residues " .. i)
    H.eq(#f.qualifiers, #g.qualifiers, "qualifier count " .. i)
    for k, q in ipairs(f.qualifiers) do
      H.eq(tostring(q.value):gsub("%s+", " "), (tostring(g.qualifiers[k].value):gsub("%s+", " ")), ("qualifier %s of %d"):format(q.key, i))
    end
  end
  -- a moved feature keeps its original qualifier lines verbatim (bla's block)
  local orig = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  local bla0 = cds(orig, "bla")
  local qlines = {}
  for k = 2, #bla0.raw_lines do qlines[#qlines + 1] = bla0.raw_lines[k] end
  local bla_q = table.concat(qlines, "\n")
  H.ok(#qlines > 5, "fixture has bla's qualifier lines")
  H.ok(cds(rec, "bla").location_text ~= bla0.location_text, "bla moved")
  H.ok(text:find(bla_q, 1, true), "bla moved but its qualifier lines are verbatim")
end)

case("snapgene_sources_save_as_genbank", function()
  local buf, st, rec, path = open("snapgene/pFA-KanMX4.dna")
  local bytes0 = H.read(path)
  se.insert(buf, 100, "ACGT")
  H.messages = {}
  pcall(vim.cmd, "write")
  H.eq(bytes0, H.read(path), ".dna file untouched")
  H.ok(table.concat(vim.tbl_map(function(m) return m.msg end, H.messages), "\n"):find("save as .gb", 1, true), last_msg())
  H.eq(true, vim.bo[buf].modified)
  local gb = tmp .. "/pFA_edited.gb"
  vim.cmd("write " .. vim.fn.fnameescape(gb))
  local back = assert(pio.read(gb)).records[1]
  H.eq(rec.seq:tostring():upper(), back.seq:tostring():upper())
  H.eq(rec:len(), back:len())
  local kan
  for _, f in ipairs(back.features) do if f:label() == "KanR" then kan = f end end
  H.eq("463..1272", kan and kan.location_text, "KanR shifted by 4")
end)

case("fasta_keeps_line_width_after_edits", function()
  local p = tmp .. "/w.fa"
  local seq = ("ACGTACGTAC"):rep(13)
  local lines = {}
  for i = 1, #seq, 50 do lines[#lines + 1] = seq:sub(i, i + 49) end
  local fh = assert(io.open(p, "w")); fh:write(">w\n" .. table.concat(lines, "\n") .. "\n"); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(p))
  local buf = vim.api.nvim_get_current_buf()
  vim.cmd("DnaInsert 1 TTTTT")
  vim.cmd("write")
  local out = H.read(p)
  local l2 = out:match("\n([^\n]+)\n")
  H.eq(50, #l2, "50 residues per line kept")
  H.eq("TTTTT" .. seq, (out:gsub("^>[^\n]*\n", ""):gsub("\n", "")))
end)

return T
