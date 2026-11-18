local H = require("helpers")
local pio = require("plasmid.io")
local orf = require("plasmid.bio.orf")
local translate = require("plasmid.bio.translate")

local T = {}

local function load(rel)
  return assert(pio.read(H.fixture(rel))).records[1]
end

function T.matches_independent_oracle_exactly()
  local oracle = vim.json.decode(H.read(H.root .. "/tests/oracle/expected/orfs.json"))
  local compared = 0
  for rel, settings in pairs(oracle) do
    local rec = load(rel)
    for key, want in pairs(settings) do
      local min_aa, starts, tbl = key:match("^(%d+):(%a+):(%d+)$")
      local got = orf.find(rec.seq:tostring():upper(), {
        circular = rec:is_circular(), min_aa = tonumber(min_aa), starts = starts, table = tonumber(tbl),
      })
      local g = vim.tbl_map(function(o) return { o.start, o.stop, o.strand, o.aa_len } end, got)
      H.eq(want, g, rel .. " " .. key)
      compared = compared + #want
    end
  end
  H.ok(compared > 3000, "ORFs compared: " .. compared)
end

--- Biological check: annotated protein-coding genes are found. Each annotated
--- CDS (standard start, no ribosomal slippage/introns) must lie inside an ORF
--- that ends at the same stop codon; the ORF starts at the CDS start unless an
--- in-frame start codon lies further upstream (longest-ORF-per-stop rule).
function T.annotated_cds_are_recovered()
  local exact, contained, total = 0, 0, 0
  local checked = {}
  for _, rel in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001422.gb", "ncbi/NC_005816.gb", "snapgene/pFA-KanMX4.dna" }) do
    local rec = load(rel)
    local L = rec:len()
    local found = orf.find(rec.seq:tostring(), { circular = rec:is_circular(), min_aa = 30, starts = "atg", table = 11 })
    for _, f in ipairs(rec.features) do
      local nt = f.type == "CDS" and rec:extract(f) or ""
      -- contiguous on the genome (no introns): consecutive parts abut, or are split only by the origin
      local single_frame = f.type == "CDS" and #rec:extents(f) >= 1
      for i = 2, #f.parts do
        local a, b = f.parts[i - 1], f.parts[i]
        local abut = (f.strand == 1 and (b.start == a.stop + 1 or (a.stop == L and b.start == 1)))
          or (f.strand == -1 and (b.stop == a.start - 1 or (a.start == 1 and b.stop == L)))
        if not abut then single_frame = false end
      end
      if single_frame and nt:sub(1, 3):upper() == "ATG" and #nt % 3 == 0 and #nt / 3 - 1 >= 30
        and translate.codon(translate.get_table(11), nt:sub(-3)) == "*" then
        total = total + 1
        local ext = rec:extents(f)
        local cds_start = f.strand == 1 and ext[1][1] or nil
        local hit
        for _, o in ipairs(found) do
          if o.strand == f.strand then
            -- same stop codon: forward = ORF.stop == CDS 3' end; reverse = ORF.start == CDS 3' end
            local lo, hi = f:bounds()
            if #ext == 2 then
              lo, hi = ext[1][1], ext[2][2] -- wrapped: 5'..L, 1..3'
              if ext[1][1] == 1 then lo, hi = ext[2][1], ext[1][2] end
            end
            local three = f.strand == 1 and hi or lo
            local five = f.strand == 1 and lo or hi
            local o3 = f.strand == 1 and o.stop or o.start
            local o5 = f.strand == 1 and o.start or o.stop
            if o3 == three then
              hit = true
              if o5 == five then exact = exact + 1 else contained = contained + 1 end
            end
          end
        end
        H.ok(hit, ("%s: CDS %s %s not found as an ORF"):format(rel, f:label(), f.location_text))
        checked[rel .. ":" .. f:label()] = true
      end
    end
  end
  io.stdout:write(("    [orf] %d annotated CDS: %d exact, %d inside a longer ORF (upstream in-frame ATG)\n"):format(
    total, exact, contained))
  -- the eligibility filter must not have silently skipped the genes we rely on
  for _, k in ipairs({ "ncbi/J01749.gb:bla", "ncbi/J01749.gb:tet", "snapgene/pFA-KanMX4.dna:KanR",
    "snapgene/pFA-KanMX4.dna:AmpR" }) do
    H.ok(checked[k], "not checked: " .. k)
  end
end

function T.phix174_gene_a_crosses_the_origin()
  local rec = load("ncbi/NC_001422.gb")
  local found = orf.find(rec.seq:tostring(), { circular = true, min_aa = 75, table = 11 })
  local wraps = vim.tbl_filter(function(o) return o.wraps end, found)
  H.ok(#wraps >= 1, "no origin-crossing ORF")
  local a
  for _, f in ipairs(rec.features) do
    if f.type == "CDS" and f.location_text:match("^join%(3981") then a = f end
  end
  H.ok(a, "fixture lacks gene A CDS")
  local hit = vim.tbl_filter(function(o) return o.wraps and o.stop == 136 and o.strand == 1 end, found)
  H.eq(1, #hit, "gene A (join(3981..5386,1..136)) should end at 136")
  -- longest ORF per stop: an in-frame ATG at 3954 lies 27 bp upstream of the
  -- annotated start, so the ORF begins there and contains gene A
  H.eq("ATG", rec.seq:sub(3954, 3956):upper())
  H.eq("join(3954..5386,1..136)", orf.location(hit[1], rec:len()))
end

function T.linear_ends_and_edge_cases()
  -- no stop before the end of a linear sequence: not an ORF
  H.eq({}, orf.find("ATG" .. ("GCT"):rep(100), { min_aa = 10 }))
  -- reverse strand coordinates
  local s = "CCC" .. require("plasmid.model.alphabet").revcomp("ATG" .. ("GCT"):rep(40) .. "TAA") .. "GG"
  local r = orf.find(s, { min_aa = 30 })
  H.eq(1, #r)
  -- ATG + 40 GCT codons = 41 residues before the stop
  H.eq({ 4, 129, -1, 41 }, { r[1].start, r[1].stop, r[1].strand, r[1].aa_len })
  -- a circular frame without any stop is not an ORF
  H.eq({}, orf.find(("ATGGCT"):rep(50), { circular = true, min_aa = 10 }))
  -- ambiguous bases neither start nor stop an ORF
  H.eq(1, #orf.find("ATG" .. ("GCN"):rep(40) .. "TAG", { min_aa = 30 }))
end

function T.performance()
  local rec = load("ncbi/NC_001133.gb")
  local t0 = vim.uv.hrtime()
  local found = orf.find(rec.seq:tostring(), { min_aa = 75 })
  local ms = (vim.uv.hrtime() - t0) / 1e6
  io.stdout:write(("    [perf] six-frame ORFs on 230 kb: %d ORFs in %.0f ms\n"):format(#found, ms))
  H.ok(ms < 3000)
end

return T
