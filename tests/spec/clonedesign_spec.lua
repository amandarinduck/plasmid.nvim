-- Assembly-aware primer design: tails that make the amplified region join a
-- destination, checked by amplifying and assembling what was designed.
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local PD = require("plasmid.primerdesign")
local A = require("plasmid.bio.assembly")
local PCR = require("plasmid.pcr")
local P = require("plasmid.bio.primer")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent! only")
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.record(buffer.get(buf))
end
local function load(rel) return assert(pio.read(H.fixture(rel))).records[1] end
local function seqof(rec) return rec.seq:tostring():upper() end

-- The region the examples amplify: pBR322's ROP CDS, which has no BsmBI site.
local SEL = { start = 1915, stop = 2106 }

local function design(buf, rec, args)
  require("plasmid.select").set(buf, SEL)
  local picked
  local select_ = vim.ui.select
  vim.ui.select = function(items, opts, cb) picked = { items = items, opts = opts }; cb(items[1]) end
  local list = PD.run(buf, args)
  vim.ui.select = select_
  return list, picked
end

case("golden_gate_tails_come_from_the_destination", function()
  local vec = load("synthetic/gg_vector.gb")
  local tails, label, warn = PD.assembly_tails(vec, { method = "golden_gate", enzyme = "BsmBI" })
  H.ok(tails ~= nil, tostring(label))
  H.eq({}, warn, "a clean destination has nothing to warn about")
  -- BsmBI CGTCTC(1/5): site, one spacer base, then the overhang the product needs
  H.eq("cgtctca", tails[1]:sub(1, 7))
  H.eq("cgtctca", tails[2]:sub(1, 7))
  H.eq("AATG", tails[1]:sub(8):upper(), "the forward tail carries the destination's right overhang")
  H.eq("GCTT", P.revcomp(tails[2]:sub(8):upper()), "the reverse tail carries its left overhang")
  H.ok(label:find("BsmBI") and label:find("pBR_dest"), label)
  -- overhangs can also be given outright
  local t2 = PD.assembly_tails(vec, { method = "golden_gate", enzyme = "BsmBI", overhangs = { "TTCG", "ACCT" } })
  H.eq("TTCG", t2[1]:sub(8):upper())
  H.eq("ACCT", P.revcomp(t2[2]:sub(8):upper()))
  -- and a non-Type IIS enzyme is refused
  local bad, why = PD.assembly_tails(vec, { method = "golden_gate", enzyme = "EcoRI" })
  H.eq(nil, bad)
  H.ok(why:find("Type IIS"), why)
end)

case("gibson_tails_are_the_destination_ends", function()
  local vec = load("synthetic/gg_vector.gb")
  local lin = A.cut(vec, { "ScaI" })[1].rec     -- linearise it first
  local tails, label = PD.assembly_tails(lin, { method = "gibson", homology = 25 })
  H.ok(tails ~= nil, tostring(label))
  local s = seqof(lin)
  H.eq(s:sub(#s - 24), tails[1]:upper(), "the forward tail continues the destination's 3' end")
  H.eq(s:sub(1, 25), P.revcomp(tails[2]:upper()), "the reverse tail runs back into its 5' end")
  H.ok(label:find("25 bp"), label)
  local bad, why = PD.assembly_tails(vec, { method = "gibson", homology = 25 })
  H.eq(nil, bad)
  H.ok(why:find("circular"), why)
end)

case("a_designed_golden_gate_pair_amplifies_and_assembles", function()
  local buf, rec = open("ncbi/J01749.gb")
  local list, picked = design(buf, rec, { "name=rop", "into=" .. H.fixture("synthetic/gg_vector.gb"),
    "enzyme=BsmBI" })
  H.ok(list and #list > 0, "candidates were offered")
  H.ok(picked.opts.prompt:find("Golden Gate into pBR_dest"), picked.opts.prompt)
  -- the two primers were added with their tails
  local primers = require("plasmid.primers").list(buffer.record(buffer.get(buf)))
  local by = {}
  for _, p in ipairs(primers) do by[p.name] = p.seq:upper() end
  H.ok(by.rop_F and by.rop_F:sub(1, 7) == "CGTCTCA", tostring(by.rop_F))
  H.ok(by.rop_R and by.rop_R:sub(1, 7) == "CGTCTCA", tostring(by.rop_R))
  -- amplify with them, then assemble the product with the destination
  local prods = PCR.products(rec.seq:tostring(), rec:is_circular(), by.rop_F, by.rop_R)
  H.ok(#prods > 0, "the designed primers amplify the template")
  local product = PCR.build(rec, prods[1], { "rop_F", "rop_R" })
  local vec = load("synthetic/gg_vector.gb")
  local parts, problems = A.type_iis_parts({ { rec = product, name = "rop" }, { rec = vec, name = "dest" } },
    "BsmBI")
  H.eq({}, problems)
  local asms = A.assemble(parts, {})
  H.eq(1, #asms, "exactly one assembly")
  local built = A.build(asms[1], { name = "check" })
  H.eq("circular", built.topology)
  H.ok(not seqof(built):find("CGTCTC"), "no BsmBI site survives in the product")
  -- the ROP CDS is inside the product
  local rop = seqof(rec):sub(SEL.start, SEL.stop)
  H.ok((seqof(built) .. seqof(built)):find(rop, 1, true) ~= nil, "the amplified region is in the plasmid")
  -- and the engine agrees with the design-time check
  local ok, msg = PD.verify_assembly(rec, list[1], { by.rop_F:gsub("^(CGTCTCA....).*", "%1"):lower(),
    by.rop_R:gsub("^(CGTCTCA....).*", "%1"):lower() }, vec, { method = "golden_gate", enzyme = "BsmBI" })
  H.ok(ok, tostring(msg))
  H.ok(msg:find("circular"), msg)
end)

case("a_designed_gibson_pair_amplifies_and_assembles", function()
  local buf, rec = open("ncbi/J01749.gb")
  local vec = load("synthetic/gg_vector.gb")
  local lin = A.cut(vec, { "ScaI" })[1].rec
  local path = tmp .. "/dest_linear.gb"
  local fh = assert(io.open(path, "wb")); fh:write(pio.format({ lin }, "genbank")); fh:close()
  local list, picked = design(buf, rec, { "name=rg", "into=" .. path, "homology=25" })
  H.ok(list and #list > 0)
  H.ok(picked.opts.prompt:find("Gibson into"), picked.opts.prompt)
  local primers = require("plasmid.primers").list(buffer.record(buffer.get(buf)))
  local by = {}
  for _, p in ipairs(primers) do by[p.name] = p.seq:upper() end
  local s = seqof(lin)
  H.eq(s:sub(#s - 24), by.rg_F:sub(1, 25), "the forward primer starts with the destination's 3' end")
  local prods = PCR.products(rec.seq:tostring(), rec:is_circular(), by.rg_F, by.rg_R)
  H.ok(#prods > 0)
  local product = PCR.build(rec, prods[1], { "rg_F", "rg_R" })
  local parts = A.gibson_parts({ { rec = product, name = "rg" }, { rec = lin, name = "dest" } })
  local asms = A.assemble(parts, { join = "homology", min_homology = 25 })
  H.eq(1, #asms)
  local built = A.build(asms[1], { name = "check" })
  H.eq("circular", built.topology)
  H.eq(product:len() + lin:len() - 50, built:len(), "each 25 bp overlap is written once")
end)

case("design_refuses_contradictory_options", function()
  local buf, rec = open("ncbi/J01749.gb")
  require("plasmid.select").set(buf, SEL)
  local msgs = {}
  local notify = vim.notify
  vim.notify = function(m, lvl) msgs[#msgs + 1] = m end
  PD.run(buf, { "into=" .. H.fixture("synthetic/gg_vector.gb"), "enzyme=BsmBI", "fwd=EcoRI" })
  PD.run(buf, { "into=nowhere_at_all" })
  PD.run(buf, { "into=" .. H.fixture("synthetic/gg_vector.gb"), "method=teleport" })
  PD.run(buf, { "into=" .. H.fixture("synthetic/gg_vector.gb"), "overhangs=AATG" })
  vim.notify = notify
  local all = table.concat(msgs, " | ")
  H.ok(all:find("into= writes the tails itself"), all)
  H.ok(all:find("no part called nowhere_at_all"), all)
  H.ok(all:find("method wants gibson or golden_gate"), all)
  H.ok(all:find("overhangs wants two sequences"), all)
  H.eq(rec:len(), 4361)
end)

case("a_site_in_the_region_is_flagged", function()
  local buf, rec = open("ncbi/J01749.gb")
  -- 3293..4153 is bla, which carries pBR322's BsaI site
  require("plasmid.select").set(buf, { start = 3293, stop = 4153 })
  local msgs = {}
  local notify = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  local select_ = vim.ui.select
  vim.ui.select = function(items, opts, cb) cb(nil) end
  PD.run(buf, { "into=" .. H.fixture("synthetic/gg_vector.gb"), "enzyme=BsaI", "overhangs=AATG/GCTT" })
  vim.ui.select = select_
  vim.notify = notify
  local all = table.concat(msgs, " | ")
  H.ok(all:find("already contains a BsaI site"), all)
end)

return T
