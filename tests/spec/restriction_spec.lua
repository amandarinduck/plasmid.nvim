local H = require("helpers")
local pio = require("plasmid.io")
local R = require("plasmid.bio.restriction")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function load(stem)
  return assert(pio.read(H.fixture("ncbi/" .. stem .. ".gb"))).records[1]
end
local function with_topology(rec, topo)
  local copy = setmetatable({}, { __index = rec })
  copy.topology = topo
  return copy
end

case("table_is_sorted_and_well_formed", function()
  local list = require("plasmid.data.enzymes")
  H.ok(#list >= 200, "enzymes: " .. #list)
  for i, e in ipairs(list) do
    H.ok(type(e.name) == "string" and e.site:match("^[ACGTRYSWKMBDHVN]+$"), e.name)
    if i > 1 then H.ok(list[i - 1].name < e.name, "sorted at " .. e.name) end
  end
end)

case("cuts_match_biopython_oracle", function()
  local expected = vim.json.decode(H.read("tests/oracle/expected/enzymes.json"))
  local enzymes = R.index(vim.deepcopy(require("plasmid.data.enzymes")))
  local checked, cuts = 0, 0
  for stem, by_topo in pairs(expected) do
    local rec = load(stem)
    for topo, by_enz in pairs(by_topo) do
      local res = R.scan(with_topology(rec, topo), enzymes)
      for _, e in ipairs(enzymes) do
        local want = by_enz[e.name]
        H.ok(want, "oracle lacks " .. e.name)
        local got = {}
        for _, s in ipairs(res[e.name]) do got[#got + 1] = s.cut end
        table.sort(want)
        if not vim.deep_equal(want, got) then
          error(("%s %s %s: want %s got %s"):format(stem, topo, e.name, vim.inspect(want), vim.inspect(got)))
        end
        checked = checked + 1
        cuts = cuts + #got
      end
    end
  end
  H.log(("%d enzyme×sequence×topology cases, %d cuts"):format(checked, cuts))
  H.eq(222 * 3 * 2, checked)
end)

case("overhangs_and_strands", function()
  local rec = load("L09137")
  local P = R.prepare(rec)
  local function one(name)
    local e = R.index({ vim.deepcopy(assert(R.get(name))) })[1]
    return R.sites(P, e), e
  end
  local s, e = one("EcoRI")
  H.eq(1, #s); H.eq(397, s[1].cut); H.eq(401, s[1].cut_bottom); H.eq("5'+4", R.end_kind(e))
  H.eq("3'+4", R.end_kind(R.get("KpnI")))
  H.eq("blunt", R.end_kind(R.get("EcoRV")))
  -- BsaI on pUC19: one site, cuts outside the recognition sequence
  s = one("BsaI")
  H.eq(1, #s)
  local site = s[1]
  local L = rec:len()
  if site.strand == 1 then
    H.eq(site.start + 7, site.cut)
  else
    H.eq(site.start - 5, site.cut, "minus-strand BsaI cuts 5 nt upstream on the top strand")
    H.eq(site.start - 1, site.cut_bottom)
  end
end)

case("digest_fragments_cover_the_molecule", function()
  local rec = load("J01749")
  local P = R.prepare(rec)
  local cuts = {}
  for _, n in ipairs({ "EcoRI", "PstI", "HindIII" }) do
    for _, s in ipairs(R.sites(P, R.get(n))) do cuts[#cuts + 1] = s end
  end
  local frags = R.digest(rec, cuts)
  H.eq(3, #frags)
  local total = 0
  for _, f in ipairs(frags) do total = total + f.length end
  H.eq(rec:len(), total)
  -- linear: n cuts give n+1 fragments
  local lin = with_topology(rec, "linear")
  local lc = {}
  for _, n in ipairs({ "PstI", "HindIII" }) do
    for _, s in ipairs(R.sites(R.prepare(lin), R.get(n))) do lc[#lc + 1] = s end
  end
  frags = R.digest(lin, lc)
  H.eq(3, #frags)
  total = 0
  for _, f in ipairs(frags) do total = total + f.length end
  H.eq(rec:len(), total)
  H.eq(1, frags[1].start); H.eq(rec:len(), frags[3].stop)
end)

case("methylation_overlap_warnings", function()
  local function sites(seq, name)
    local P = R.prepare({ seq = seq, topology = "linear" })
    local list = R.sites(P, R.get(name))
    return P, list
  end
  H.ok(R.get("XbaI").dam, "XbaI is Dam-sensitive in the table")
  -- TCTAGA followed by TC: the GATC overlaps the XbaI site
  local P, s = sites("AAAAATCTAGATCAAAAA", "XbaI")
  H.eq(1, #s); H.eq({ "dam" }, R.methylation(P, s[1]))
  -- GATC three bases away does not overlap
  P, s = sites("AAAAATCTAGAGGATCCAAAA", "XbaI")
  H.eq(1, #s); H.eq({}, R.methylation(P, s[1]))
  -- an enzyme without the flag never warns
  P, s = sites("AAAAAGAATTCGATCAAAA", "EcoRI")
  H.eq({}, R.methylation(P, s[1]))
  H.eq("dam", R.NEEDS_METHYLATION.DpnI)
end)

return T
