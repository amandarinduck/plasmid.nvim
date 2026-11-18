local H = require("helpers")
local pio = require("plasmid.io")
local alphabet = require("plasmid.model.alphabet")
local config = require("plasmid.config")

local T = {}

function T.protein_fasta_detected_by_content_and_extension()
  local text = H.read(H.fixture("protein/proteins.faa"))
  for _, name in ipairs({ "proteins.faa", "proteins.fasta", "proteins.txt" }) do
    local res = assert(pio.parse(text, name))
    H.eq(4, #res.records, name)
    for _, r in ipairs(res.records) do
      H.eq("protein", r.alphabet, name .. " " .. r.id)
      H.eq("aa", r:unit())
      H.eq(false, r:is_circular())
    end
  end
end

function T.alphabet_guess_edge_cases()
  H.eq("dna", alphabet.guess("ACGTACGTTTGACNNNN"))
  H.eq("dna", alphabet.guess("NNKNNKNNKNNKRYSWGATCCATG"), "degenerate library stays DNA")
  H.eq("dna", alphabet.guess("ACGT-ACGT*ACGT"), "gaps/stops are ignored")
  H.eq("rna", alphabet.guess("AUGGCUUACGAUUAG"))
  H.eq("protein", alphabet.guess("MKSNNALIVILGTVTLDAVGIGLVMPVLPGLLRDIVHSD*"))
  H.eq("protein", alphabet.guess("MKKRREREFLWTPIDRQEREKKKEKEKRTERKRHRREEYIDR*N*APARSRS"))
  H.eq("protein", alphabet.guess("GAVLIMFWPSTCYNQDEKRH"))
end

--- Independent cross-check: the TetA protein FASTA from NCBI equals the
--- translation the plugin derives from pBR322's DNA.
function T.teta_protein_equals_dna_translation()
  local prot = assert(pio.read(H.fixture("protein/proteins.faa")))
  local teta
  for _, r in ipairs(prot.records) do if r.id == "AAB59735.1" then teta = r end end
  local dna = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  local cds
  for _, f in ipairs(dna.features) do if f.type == "CDS" and f:get("gene") == "tet" then cds = f end end
  local map = require("plasmid.view.features").codon_map(dna, cds)
  local aa = {}
  for _, p in ipairs(require("plasmid.model.location").positions(cds.parts)) do
    if map[p] then aa[#aa + 1] = map[p] end
  end
  H.eq(teta.seq:tostring(), (table.concat(aa):gsub("%*$", "")))
end

function T.protein_linear_view()
  config.setup({ linear = { color_bases = true } })
  require("plasmid.view.highlights").setup()
  local rec = assert(pio.read(H.fixture("protein/proteins.faa"))).records[1]
  local lay = require("plasmid.view.linear").build(rec, config.options, 100)
  H.ok(lay.lines[1]:find("142 aa", 1, true), lay.lines[1])
  for _, info in pairs(lay.rows) do
    H.ok(info.kind ~= "comp" and info.kind ~= "aa", "protein view drew a " .. info.kind .. " row")
  end
  local row0 = lay.block_row[0]
  local groups = {}
  for _, h in ipairs(lay:row_highlights(row0)) do groups[h[3]] = true end
  H.ok(groups.PlasmidAAHydrophobic and groups.PlasmidAAPositive, vim.inspect(groups))
  config.setup({})
end

function T.protein_genbank_roundtrip_keeps_alphabet()
  local res = assert(pio.read(H.fixture("protein/proteins.faa")))
  local text = assert(pio.format(res.records, "genbank", res.file, "fasta"))
  local back = assert(pio.parse(text, "x.gp"))
  H.eq(4, #back.records)
  for i, r in ipairs(back.records) do
    H.eq("protein", r.alphabet)
    H.eq(res.records[i].seq:tostring(), r.seq:tostring())
  end
end

function T.rna_complement_uses_uracil()
  H.eq("UACG", alphabet.complement_for("AUGC", "rna"))
  H.eq("TACG", alphabet.complement_for("ATGC", "dna"))
end

return T
