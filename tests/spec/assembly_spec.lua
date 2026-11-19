-- The cloning engine: fragment ends, ligation, Golden Gate, Gibson, and the
-- records that come out of them. No UI here (see clone_spec).
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local A = require("plasmid.bio.assembly")
local Record = require("plasmid.model.record").Record

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function load(stem) return assert(pio.read(H.fixture("ncbi/" .. stem .. ".gb"))).records[1] end
local function up(s) return s:upper() end
local function rotation_of(a, b)   -- is b a rotation of a?
  return #a == #b and (up(a) .. up(a)):find(up(b), 1, true) ~= nil
end
local function seqof(rec) return rec.seq:tostring():upper() end
local function labels(rec)
  local out = {}
  for _, f in ipairs(rec.features) do
    if f.type ~= "source" then out[f.type .. "|" .. f:label() .. "|" .. rec:extract(f)] = true end
  end
  return out
end

local function rnd(n)
  local b, t = { "A", "C", "G", "T" }, {}
  for i = 1, n do t[i] = b[math.random(4)] end
  return table.concat(t)
end
local function mk(name, seq, topo)
  return Record.new({ name = name, seq = seq, topology = topo or "linear", alphabet = "dna" })
end

---------------------------------------------------------------------------
-- Ends
---------------------------------------------------------------------------

case("ends_of_a_digest_fragment", function()
  local puc = load("L09137")
  local frags = A.cut(puc, { "EcoRI" })
  H.eq(1, #frags, "one cut in a circle gives one fragment")
  local f = frags[1]
  -- EcoRI G^AATTC: a 4 base 5' overhang at both ends, AATT on the top strand
  H.eq("5′ AATT", A.end_text(f, "left"))
  H.eq("5′ AATT", A.end_text(f, "right"))
  H.eq(1, f.left.strand, "the top strand protrudes on the left")
  H.eq(-1, f.right.strand, "the bottom strand protrudes on the right")
  H.eq(puc:len() + 4, f.rec:len(), "the fragment spans the circle plus one overhang")
  -- PstI CTGCA^G: a 4 base 3' overhang, written on the top strand like the digest table
  local p = A.cut(load("J01749"), { "PstI" })[1]
  H.eq("3′ TGCA", A.end_text(p, "left"))
  H.eq("3′ TGCA", A.end_text(p, "right"))
  H.eq(-1, p.left.strand, "the bottom strand protrudes on the left of a 3' overhang")
  H.eq("TGCA", A.overhang_strand(p, "left"), "the protruding strand reads 5'->3'")
end)

case("flipping_a_fragment_swaps_its_ends", function()
  local f = A.cut(load("J01749"), { "EcoRI", "HindIII" })[1]
  local g = A.flip(f)
  H.eq(require("plasmid.model.alphabet").revcomp(seqof(f.rec)), seqof(g.rec))
  H.eq(A.end_seq(f, "right"), A.end_seq(g, "left"), "the right end becomes the left one")
  H.eq(f.right.len, g.left.len)
  H.eq(-f.right.strand, g.left.strand)
  local back = A.flip(g)
  H.eq(seqof(f.rec), seqof(back.rec))
  H.eq(f.left, back.left)
  H.eq(f.right, back.right)
end)

case("incompatible_ends_report_why", function()
  local pbr = load("J01749")
  local eco = A.cut(pbr, { "EcoRI" })[1]
  local pst = A.cut(pbr, { "PstI" })[1]
  local j, why = A.compatible(eco, pst)
  H.eq(nil, j, "AATT and TGCA do not anneal")
  H.ok(why:find("AATT") and why:find("TGCA"), "the reason names both overhangs: " .. tostring(why))
  H.ok(A.compatible(eco, eco) ~= nil, "EcoRI ends anneal to each other")
end)

---------------------------------------------------------------------------
-- Ligation
---------------------------------------------------------------------------

case("religating_a_single_cut_restores_the_plasmid", function()
  local puc = load("L09137")
  local frags = A.cut(puc, { "EcoRI" })
  local asms = A.assemble(frags, {})
  H.eq(1, #asms)
  H.eq(true, asms[1].circular)
  local rec = A.build(asms[1], { name = "religated" })
  H.eq(puc:len(), rec:len())
  H.eq("circular", rec.topology)
  H.ok(rotation_of(seqof(puc), seqof(rec)), "the sequence is the plasmid again")
end)

case("two_fragments_religate_and_features_survive", function()
  local pbr = load("J01749")
  local frags = A.cut(pbr, { "EcoRI", "HindIII" })
  H.eq(2, #frags)
  local asms = A.assemble(frags, {})
  H.eq(1, #asms, "the two overhangs fix one order and orientation")
  local rec = A.build(asms[1], { name = "religated" })
  H.ok(rotation_of(seqof(pbr), seqof(rec)), "the sequence is pBR322 again")
  -- every original feature is back, including the two the cuts had split
  local want, got = labels(pbr), labels(rec)
  local missing = {}
  for k in pairs(want) do if not got[k] then missing[#missing + 1] = k end end
  H.eq({}, missing, "features carried over and re-joined")
end)

case("blunt_ends_ligate_in_both_orientations", function()
  local pbr = load("J01749")
  local frags = A.cut(pbr, { "EcoRV", "ScaI" })
  H.eq(2, #frags)
  for _, f in ipairs(frags) do
    H.eq("blunt", A.end_text(f, "left"))
    H.eq("blunt", A.end_text(f, "right"))
  end
  local asms = A.assemble(frags, {})
  H.eq(2, #asms, "blunt fragments join either way round")
  local a = A.build(asms[1], { name = "a" })
  local b = A.build(asms[2], { name = "b" })
  H.eq(pbr:len(), a:len())
  H.eq(pbr:len(), b:len())
  local one = rotation_of(seqof(pbr), seqof(a)) or rotation_of(seqof(pbr), seqof(b))
  H.ok(one, "one of the two is the original plasmid")
end)

case("linear_assembly_when_nothing_closes", function()
  -- one fragment with an AGCT end and an AATT end: it cannot close on itself
  local lone = A.cut(load("J01749"), { "EcoRI", "HindIII" })[1]
  local asms, fell_back = A.assemble({ lone }, {})
  H.eq(1, #asms)
  H.eq(false, asms[1].circular, "no circle, so the linear molecule is reported")
  H.eq(true, fell_back, "the circular search fell back to linear")
  -- a blunt fragment, by contrast, does circularise on its own
  local blunt = A.from_record(mk("blunt_part", rnd(200)))
  local self_closed = A.assemble({ blunt }, {})
  H.eq(1, #self_closed)
  H.eq(true, self_closed[1].circular, "blunt ends self-ligate")
  -- an intact circle has no ends to assemble with
  local circ = A.cut(mk("circle", ("ACGT"):rep(125), "circular"), { "EcoRI" })
  H.eq(true, circ[1].circular)
  H.eq(false, pcall(A.assemble, circ, {}))
end)

---------------------------------------------------------------------------
-- Golden Gate
---------------------------------------------------------------------------

local function gg_parts(ovl1, ovl2, body)
  math.randomseed(11)
  local backbone, stuffer = rnd(900), rnd(60)
  local insert = mk("insert", "GGTCTC" .. "a" .. ovl1 .. body .. ovl2 .. "a" .. "GAGACC")
  local vector = mk("vector", ovl2 .. backbone .. ovl1 .. "a" .. "GAGACC" .. stuffer .. "GGTCTC" .. "a",
    "circular")
  return insert, vector, backbone
end

case("golden_gate_assembles_insert_into_vector", function()
  local body = rnd(300)
  local insert, vector, backbone = gg_parts("AATG", "GCTT", body)
  local parts, problems = A.type_iis_parts({ insert, vector }, "BsaI")
  H.eq({}, problems, "clean parts have no problems")
  H.eq(2, #parts)
  H.eq("5′ AATG", A.end_text(parts[1], "left"))
  H.eq("5′ GCTT", A.end_text(parts[1], "right"))
  H.ok(not A.end_seq(parts[1], "left"):find("GGTCTC"), "the enzyme sites are gone")
  local asms = A.assemble(parts, {})
  H.eq(1, #asms, "the overhangs allow exactly one assembly")
  local rec = A.build(asms[1], { name = "GG" })
  H.eq("circular", rec.topology)
  H.ok(rotation_of(("GCTT" .. backbone .. "AATG" .. body), seqof(rec)), "insert joined to backbone")
  H.ok(not seqof(rec):find("GGTCTC"), "no BsaI site is left in the product")
end)

case("golden_gate_reports_part_problems", function()
  local _, vector = gg_parts("AATG", "GCTT", rnd(100))
  local no_site = mk("plain", rnd(300))
  local _, problems = A.type_iis_parts({ no_site, vector }, "BsaI")
  H.ok(#problems > 0 and problems[1]:find("no BsaI site"), vim.inspect(problems))
  -- a palindromic overhang, and one that is used by three ends
  local pal = mk("pal", "GGTCTC" .. "a" .. "AATT" .. rnd(100) .. "GCTT" .. "a" .. "GAGACC")
  local _, p2 = A.type_iis_parts({ pal, vector }, "BsaI")
  local joined = table.concat(p2, " | ")
  H.ok(joined:find("palindromic"), joined)
  local _, p3 = A.type_iis_parts({}, "NotAnEnzyme")
  H.ok(p3[1]:find("unknown enzyme"), p3[1])
  local _, p4 = A.type_iis_parts({ vector }, "EcoRI")
  H.ok(p4[1]:find("Type IIS"), p4[1])
end)

---------------------------------------------------------------------------
-- Gibson
---------------------------------------------------------------------------

case("gibson_joins_parts_by_terminal_homology", function()
  math.randomseed(23)
  local h1, h2, h3 = rnd(25), rnd(30), rnd(22)
  local a1, a2, a3 = rnd(400), rnd(500), rnd(300)
  local parts, problems = A.gibson_parts({
    mk("gA", h1 .. a1 .. h2), mk("gB", h2 .. a2 .. h3), mk("gC", h3 .. a3 .. h1) })
  H.eq({}, problems)
  local asms = A.assemble(parts, { join = "homology", min_homology = 20 })
  H.eq(1, #asms)
  local rec = A.build(asms[1], { name = "gibson" })
  H.eq("circular", rec.topology)
  H.ok(rotation_of(h1 .. a1 .. h2 .. a2 .. h3 .. a3, seqof(rec)), "each overlap is written once")
  -- the junction annotations mark the overlaps
  local notes = {}
  for _, f in ipairs(rec.features) do notes[#notes + 1] = f:label() end
  H.ok(vim.tbl_contains(notes, "overlap 25 bp"), vim.inspect(notes))
end)

case("gibson_needs_long_enough_homology", function()
  math.randomseed(29)
  local h = rnd(12)
  local parts = A.gibson_parts({ mk("x", h .. rnd(300) .. h), mk("y", h .. rnd(200) .. h) })
  H.eq(0, #A.assemble(parts, { join = "homology", min_homology = 20 }), "12 bp is below the floor")
  H.ok(#A.assemble(parts, { join = "homology", min_homology = 10 }) > 0, "lowering the floor assembles")
  local circ = mk("circle", rnd(500), "circular")
  local _, problems = A.gibson_parts({ circ })
  H.ok(problems[1]:find("circular"), problems[1])
end)

---------------------------------------------------------------------------
-- Products
---------------------------------------------------------------------------

case("product_annotations_and_round_trip", function()
  local pbr = load("J01749")
  local asms = A.assemble(A.cut(pbr, { "EcoRI", "HindIII" }), {})
  local rec = A.build(asms[1], { name = "annotated" })
  local parts, junctions = 0, 0
  for _, f in ipairs(rec.features) do
    local note = f:get("note") or ""
    if note:find("^part ") then parts = parts + 1 end
    if note:find("junction") then junctions = junctions + 1 end
  end
  H.eq(2, parts)
  H.eq(2, junctions)
  -- the product is a normal record: it writes and parses back unchanged
  local text = pio.format({ rec }, "genbank")
  local again = assert(pio.parse(text, nil, "genbank")).records[1]
  H.eq(seqof(rec), seqof(again))
  H.eq(#rec.features, #again.features)
  H.eq("circular", again.topology)
end)

return T
