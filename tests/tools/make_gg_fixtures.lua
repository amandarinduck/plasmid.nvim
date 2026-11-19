-- Build the two Golden Gate demo records in tests/fixtures/synthetic/ from
-- real pBR322 sequence (tests/fixtures/ncbi/J01749.gb), so the cloning
-- screenshots and the help examples are reproducible from the tree:
--
--   gg_insert.gb  the tet CDS as a part, flanked by outward-facing BsmBI sites
--                 and the AATG / GCTT fusion sites
--   gg_vector.gb  a circular destination vector: pBR322's replication origin
--                 and bla, with a stuffer between inward-facing BsmBI sites
--
-- BsmBI rather than BsaI because pBR322 carries a BsaI site inside bla, and a
-- destination vector must not be cut by its own enzyme.
--
-- Run:  nvim -l tests/tools/make_gg_fixtures.lua
vim.opt.rtp:prepend(".")
require("plasmid").setup({})

local pio = require("plasmid.io")
local A = require("plasmid.bio.assembly")
local R = require("plasmid.bio.restriction")
local T = require("plasmid.model.transaction")

local ENZYME = "BsmBI"
local OVL_START, OVL_END = "AATG", "GCTT"   -- fusion sites (MoClo convention)
local INSERT = { 86, 1276 }                 -- tet CDS
local BACKBONE = { 2200, 4361 }             -- rep_origin + bla, no BsmBI site
local STUFFER = 120                         -- bases dropped when the vector is cut

local e = assert(R.get(ENZYME))
local spacer = ("a"):rep(math.max(1, e.cut5 - #e.site))
local site = e.site:lower()
local function revcomp(s) return require("plasmid.model.alphabet").revcomp(s) end

local pbr = assert(pio.read("tests/fixtures/ncbi/J01749.gb")).records[1]
local L = pbr:len()

local function region(a, b)
  local len = (b - a) % L + 1
  local rec = require("plasmid.digest").extract(pbr, { start = a, stop = b, length = len })
  rec.history = nil
  return rec
end

local function prepend(rec, s) T.apply(rec, { op = "splice", start = 1, stop = 0, ins = s }) end
local function append(rec, s)
  local n = rec:len()
  T.apply(rec, { op = "splice", start = n + 1, stop = n, ins = s })
end

-- the part: site -> spacer -> AATG -> tet -> GCTT -> spacer -> site (outward)
local insert = region(INSERT[1], INSERT[2])
prepend(insert, site .. spacer .. OVL_START:lower())
append(insert, OVL_END:lower() .. spacer .. revcomp(site:upper()):lower())
insert.name = "tet_part"
insert.description = ("pBR322 tet CDS as a Golden Gate part (%s, %s / %s fusion sites)")
  :format(ENZYME, OVL_START, OVL_END)
insert.topology = "linear"

-- the destination: GCTT -> backbone -> AATG -> spacer -> site -> stuffer -> site (inward)
local vector = region(BACKBONE[1], BACKBONE[2])
local stuffer = pbr.seq:tostring():sub(1500, 1500 + STUFFER - 1)
prepend(vector, OVL_END:lower())
append(vector, OVL_START:lower() .. spacer .. revcomp(site:upper()):lower()
  .. stuffer .. site .. spacer)
vector.name = "pBR_dest"
vector.description = ("Golden Gate destination vector from pBR322 (%s, %s / %s fusion sites)")
  :format(ENZYME, OVL_END, OVL_START)
vector:set_topology("circular")
vector.history = nil

-- checks: one site each, and the two parts make one circular product
for _, rec in ipairs({ insert, vector }) do
  local n = #R.sites(R.prepare(rec), e)
  assert(n == 2, ("%s has %d %s sites, expected 2"):format(rec.name, n, ENZYME))
end
local parts, problems = A.type_iis_parts({ insert, vector }, ENZYME)
assert(#problems == 0, vim.inspect(problems))
local asms = A.assemble(parts, {})
assert(#asms == 1, ("%d assemblies"):format(#asms))
local product = A.build(asms[1], { name = "product" })
assert(product:is_circular())
assert(not product.seq:tostring():upper():find(e.site:upper()), "a site survived the assembly")

for _, rec in ipairs({ insert, vector }) do
  local path = ("tests/fixtures/synthetic/gg_%s.gb"):format(rec == insert and "insert" or "vector")
  local fh = assert(io.open(path, "wb"))
  fh:write(pio.format({ rec }, "genbank"))
  fh:close()
  io.stdout:write(("%s  %s  %d bp  %d features\n"):format(path, rec.name, rec:len(), #rec.features))
end
io.stdout:write(("product: %d bp circular, %d features\n"):format(product:len(), #product.features))
