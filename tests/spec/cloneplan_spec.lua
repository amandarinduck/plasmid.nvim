-- The bench plan: what goes in the tube, what should come out, what else
-- could, and the checks that tell those apart. The assemblers themselves are
-- tested in assembly_spec and clone_spec; here the subject is the plan.
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local C = require("plasmid.clone")
local CP = require("plasmid.cloneplan")
local D = require("plasmid.digest")
local Record = require("plasmid.model.record").Record

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")

vim.o.columns, vim.o.lines, vim.o.hidden = 150, 50, true

local function reset()
  C.staged = {}
  C.last = nil
  C.last_assembly = nil
  CP.last = nil
  vim.cmd("silent! only")
  vim.cmd("silent! %bwipeout!")
end

local function rnd(n)
  local b, t = { "A", "C", "G", "T" }, {}
  for i = 1, n do t[i] = b[math.random(4)] end
  return table.concat(t)
end
local function mk(name, seq, topo)
  return Record.new({ name = name, seq = seq, topology = topo or "linear", alphabet = "dna" })
end
local function write_rec(rec, stem)
  local p = ("%s/%s.gb"):format(tmp, stem)
  assert(pio.write(p, { rec }, "genbank"))
  return p
end
local function open_path(path)
  vim.cmd("silent edit! " .. vim.fn.fnameescape(path))
  return vim.api.nvim_get_current_buf()
end
local function open_fixture(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  return open_path(dst), dst
end
local function load_fixture(rel) return assert(pio.read(H.fixture(rel))).records[1] end

--- The lines of the plan report buffer last shown.
local function report()
  local b = CP.last
  H.ok(b and vim.api.nvim_buf_is_valid(b), "no plan report buffer")
  return table.concat(vim.api.nvim_buf_get_lines(b, 0, -1, false), "\n"), b
end

local function labelled(plan, text)
  for _, b in ipairs(plan.byproducts) do
    if b.label:find(text, 1, true) then return b end
  end
  return nil
end
--- The by-product labels, for a readable failure message (vim.inspect of a
--- by-product would print its whole record).
local function labels(plan)
  local out = {}
  for i, b in ipairs(plan.byproducts) do out[i] = ("%s %s (%d bp)"):format(b.tag, b.label, b.length) end
  return table.concat(out, "; ")
end

---------------------------------------------------------------------------
-- The three reactions
---------------------------------------------------------------------------

--- pUC19 cut with EcoRI, plus a pBR322 segment given EcoRI ends: one overhang
--- on all four ends, so every classic by-product is on the table.
local function ecori_ligation(body_from, body_to, stem)
  local pbr = load_fixture("ncbi/J01749.gb")
  local body = pbr.seq:tostring():sub(body_from, body_to):upper()
  local ins = mk(stem, "GGGAATTCC" .. body .. "GGAATTCCC")
  local ipath = write_rec(ins, stem)
  local _, vpath = open_fixture("ncbi/L09137.gb")
  local vbuf = open_path(vpath)
  local vdb = D.run(vbuf, { "EcoRI" })
  C.add(vdb, "1")
  local ibuf = open_path(ipath)
  local idb = D.run(ibuf, { "EcoRI" })
  C.add(idb, "2")   -- the middle piece, the insert proper
  local _, rec, ctx = C.run(vdb, "ligation", {}, { first = true })
  H.ok(rec, "the ligation produced nothing")
  return rec, ctx, vdb
end

case("golden_gate_plan_names_the_parts_the_reaction_and_the_product", function()
  reset()
  local _, ipath = open_fixture("synthetic/gg_insert.gb")
  local vbuf, vpath = open_fixture("synthetic/gg_vector.gb")
  local _, rec = C.run(vbuf, "goldengate", { ipath, vpath }, { enzyme = "BsmBI", first = true })
  H.ok(rec, "the Golden Gate produced nothing")
  local plan = CP.build(C.last_assembly, {})
  H.eq("goldengate", plan.method)
  H.eq("BsmBI", plan.enzyme)
  H.eq(rec:len(), plan.length)
  H.eq("circular", plan.topology)
  -- both parts, where each came from, how it was prepared, and its ends
  H.eq(2, #plan.parts)
  local by = {}
  for _, p in ipairs(plan.parts) do by[p.name] = p end
  H.ok(by.tet_part and by.pBR_dest, vim.inspect(vim.tbl_keys(by)))
  H.eq("BsmBI digest", by.tet_part.prepared)
  H.eq("BsmBI digest", by.pBR_dest.prepared)
  H.ok(by.tet_part.from:find("gg_insert", 1, true), by.tet_part.from)
  -- the two fusion sites are the junctions, with their position in the product
  H.eq(2, #plan.junctions)
  local seqs = {}
  for _, j in ipairs(plan.junctions) do seqs[j.seq] = j.at end
  H.ok(seqs.AATG and seqs.GCTT, vim.inspect(seqs))
  for _, j in ipairs(plan.junctions) do
    H.eq(4, j.len)
    H.eq("sticky", j.kind)
    H.eq(j.seq, rec.seq:tostring():upper():sub(j.at - 3, j.at), "the junction is where the plan says it is")
  end
  -- the destination went in as a circle, and an uncut circle is the background
  -- a Golden Gate is screened against
  H.ok(labelled(plan, "pBR_dest carried through uncut"), labels(plan))
  local text = CP.text(plan)
  for _, want in ipairs({ "Golden Gate (BsmBI)", " reaction", " parts", " expected product",
    "by-products", " diagnostic digest", " colony PCR", " sequencing primers" }) do
    H.ok(text:find(want, 1, true), ("plan has no %q"):format(want))
  end
  H.log(("golden gate plan: %d bp product, %d parts, %d junctions, by-products: %s")
    :format(plan.length, #plan.parts, #plan.junctions, labels(plan)))
end)

case("gibson_plan_reports_the_overlaps_and_the_vector_closing_on_itself", function()
  reset()
  math.randomseed(31)
  -- One arm on both ends of the backbone and both ends of the insert: the
  -- insert goes in either way round and the backbone can close on itself,
  -- which is the background a Gibson reaction actually has.
  local arm, backbone, body = rnd(25), rnd(1200), rnd(500)
  local vec = mk("vecA", arm .. backbone .. arm)
  local ins = mk("insA", arm .. body .. arm)
  local vpath, ipath = write_rec(vec, "gib_vec"), write_rec(ins, "gib_ins")
  local buf = open_path(vpath)
  local _, rec = C.run(buf, "gibson", { vpath, ipath }, { first = true })
  H.ok(rec, "the Gibson produced nothing")
  local plan = CP.build(C.last_assembly, {})
  H.eq("gibson", plan.method)
  H.eq(2, #plan.junctions)
  for _, j in ipairs(plan.junctions) do
    H.eq("homology", j.kind)
    H.eq(25, j.len)
    H.eq(arm, j.seq)
  end
  -- each junction writes its arm once, so both copies of the arm survive
  H.eq(2 * #arm + #backbone + #body, plan.length)
  for _, p in ipairs(plan.parts) do H.eq("whole record", p.prepared) end
  -- the backgrounds a Gibson reaction really has: either part closing on
  -- itself across its own two arms, and either part taken twice
  H.ok(labelled(plan, "vecA alone (self-closure)"), labels(plan))
  H.ok(labelled(plan, "insA alone (self-closure)"), labels(plan))
  H.ok(labelled(plan, "vecA twice"), labels(plan))
  -- and the one a ligation has and this does not: a homology arm anneals in
  -- one direction only, so the insert cannot go in the other way round
  H.eq(nil, labelled(plan, "reverse orientation"), labels(plan))
  for _, b in ipairs(plan.byproducts) do H.ok(b.kind ~= "rearranged", b.label) end
  H.log(("gibson plan: %d bp product, %d bp arms, %d by-products (%s)")
    :format(plan.length, plan.junctions[1].len, #plan.byproducts, labels(plan)))
end)

case("restriction_ligation_plan_finds_the_classic_by_products", function()
  reset()
  local rec, ctx = ecori_ligation(1200, 2100, "insert")
  local plan = CP.build(ctx, {})
  H.eq("ligation", plan.method)
  H.eq(rec:len(), plan.length)
  -- the four families, each a molecule the engine really built
  local self_vec = labelled(plan, "alone (self-closure)")
  local reverse = labelled(plan, "in the reverse orientation")
  local tandem = labelled(plan, "twice:")
  H.ok(self_vec and reverse and tandem, labels(plan))
  H.eq(2686, labelled(plan, "SYNPUC19CV_1 alone").length, "pUC19 religated is pUC19")
  H.eq(plan.length, reverse.length, "the reverse-orientation product is the same size")
  H.ok(tandem.length > plan.length, "a tandem is longer than the product")
  -- every by-product is a distinct molecule, and none of them is the product
  local seen = {}
  for _, b in ipairs(plan.byproducts) do
    local key = b.rec.seq:tostring():upper()
    H.ok(not seen[key], "duplicate by-product " .. b.label)
    seen[key] = true
  end
  H.log(("ligation plan: %d bp product, %d by-products (%d self, %d reverse, %d tandem)"):format(
    plan.length, #plan.byproducts,
    #vim.tbl_filter(function(b) return b.kind == "self" end, plan.byproducts),
    #vim.tbl_filter(function(b) return b.kind == "rearranged" end, plan.byproducts),
    #vim.tbl_filter(function(b) return b.kind == "tandem" end, plan.byproducts)))
end)

---------------------------------------------------------------------------
-- The diagnostic digest
---------------------------------------------------------------------------

case("the_diagnostic_search_finds_a_set_that_tells_the_product_apart", function()
  reset()
  local t0 = vim.uv.hrtime()
  local _, ctx = ecori_ligation(1200, 2100, "insert")
  local plan = CP.build(ctx, {})
  local ms = (vim.uv.hrtime() - t0) / 1e6
  local d = plan.diagnostic
  H.ok(d, "no discriminating set was found")
  H.ok(#d.names >= 1 and #d.names <= 2, vim.inspect(d.names))
  H.eq(#plan.byproducts, d.told, "every by-product is told apart")
  H.ok(d.all, "the verdict says so too")
  -- the bands really differ: each by-product has a gel row the product has
  -- not, or the other way round
  local rows = CP.lane_rows(d.lanes)
  for i, b in ipairs(plan.byproducts) do
    local diff = false
    for r in pairs(rows[1]) do if rows[i + 1][r] == nil then diff = true end end
    for r in pairs(rows[i + 1]) do if rows[1][r] == nil then diff = true end end
    H.ok(diff, ("no band tells %s apart, but the plan says one does"):format(b.label))
    H.ok(d.per[i] and #d.per[i].bands > 0, b.label)
  end
  -- the sizes in the lanes are a real digest of each molecule
  local total = 0
  for _, s in ipairs(d.lanes[1].sizes) do total = total + s end
  H.eq(plan.length, total, "the product's bands add up to the product")
  local text = CP.text(plan)
  H.ok(text:find(table.concat(d.names, " + "), 1, true), text:sub(1, 400))
  H.log(("diagnostic search: %s tells all %d by-products apart (%d sets from %d enzymes, plan %.0f ms)")
    :format(table.concat(d.names, " + "), #plan.byproducts, plan.search.tried, plan.search.candidates, ms))
end)

case("the_diagnostic_search_declines_when_no_set_discriminates", function()
  reset()
  math.randomseed(7)
  -- The same ligation with an insert too short to carry a six-base site: the
  -- two orientations put the same block between the same vector positions, so
  -- no enzyme that cuts the vector can tell them apart, and none that cuts the
  -- insert exists.
  local tiny = mk("tiny", "GGGAATTCC" .. rnd(40) .. "GGAATTCCC")
  local tpath = write_rec(tiny, "tiny")
  local _, vpath = open_fixture("ncbi/L09137.gb")
  local vdb = D.run(open_path(vpath), { "EcoRI" })
  C.add(vdb, "1")
  local tdb = D.run(open_path(tpath), { "EcoRI" })
  C.add(tdb, "2")
  local _, rec, ctx = C.run(vdb, "ligation", {}, { first = true })
  H.ok(rec, "the ligation produced nothing")
  local plan = CP.build(ctx, {})
  H.eq(nil, plan.diagnostic, "a set was offered where none discriminates")
  local reverse = labelled(plan, "in the reverse orientation")
  H.ok(reverse, labels(plan))
  H.eq(plan.length, reverse.length)
  -- the one it cannot tell apart is named, and it is the reverse orientation
  local never = {}
  for _, b in ipairs(plan.search.never) do never[#never + 1] = b.label end
  H.eq({ reverse.label }, never)
  -- a partial check exists and is offered as a partial check, not as the answer
  H.ok(plan.best and plan.best.told > 0 and not plan.best.all, vim.inspect(plan.best and plan.best.told))
  local text = CP.text(plan)
  H.ok(text:find("no set of up to 2 enzymes tells the product from every by-product", 1, true), text)
  H.ok(text:find("nothing distinguishes #" .. reverse.n, 1, true), text)
  H.ok(text:find("will be run and believed", 1, true), text)
  H.log(("diagnostic search declines: %d of %d by-products told apart by %s; %s is indistinguishable")
    :format(plan.best.told, #plan.byproducts, table.concat(plan.best.names, " + "), reverse.label))
end)

case("two_sizes_in_one_gel_row_are_one_band_and_not_a_difference", function()
  reset()
  -- The rule the search uses, stated on its own: the gel model decides what
  -- "resolve" means, and sizes that land in the same row do not.
  H.ok(CP.resolvable({ 3000, 1000 }, { 3000, 500 }), "1000 and 500 resolve")
  H.ok(not CP.resolvable({ 3000, 1000 }, { 3000, 1010 }), "1000 and 1010 do not")
  H.ok(CP.resolvable({ 2000 }, {}), "no product at all is a difference")
  H.log("gel resolution: 1000/500 resolve, 1000/1010 do not")
end)

---------------------------------------------------------------------------
-- Primers
---------------------------------------------------------------------------

case("colony_pcr_is_sized_on_the_product_and_on_every_by_product", function()
  reset()
  local _, ctx = ecori_ligation(1200, 2100, "insert")
  local plan = CP.build(ctx, {})
  local PCR = require("plasmid.pcr")
  H.eq(#plan.junctions, #plan.pcr)
  local checks = 0
  for _, r in ipairs(plan.pcr) do
    H.ok(not r.error, tostring(r.error))
    H.eq(#plan.byproducts, #r.on)
    -- the pair really spans the junction: one product on the product, and the
    -- template span covers the junction base
    local prods = PCR.products(plan.product.seq:tostring(), true, r.fwd, r.rev,
      { binding = (require("plasmid.config").options.primers or {}).binding })
    H.eq(1, #prods, ("junction %d gives %d products on the product"):format(r.junction, #prods))
    local p = prods[1]
    local L = plan.length
    H.ok((r.at - p.start) % L < p.region, "the junction is inside the product span")
    H.eq({ p.length }, r.product_sizes)
    for i, b in ipairs(plan.byproducts) do
      local on = r.on[i]
      -- the verdict and the sizes agree with each other
      H.eq(CP.resolvable(r.product_sizes, on.sizes), on.tells, b.label)
      if #on.sizes == 0 then H.ok(on.tells, "no product at all is a check") end
      if on.tells then checks = checks + 1 end
    end
  end
  -- the vector closing on itself gives no product with a pair across a
  -- junction, which is the whole point of a colony PCR
  local vec = labelled(plan, "SYNPUC19CV_1 alone")
  H.eq({}, plan.pcr[1].on[vec.n].sizes)
  local text = CP.text(plan)
  H.ok(text:find("same size — not a check", 1, true) or text:find("tells it apart", 1, true), text)
  H.log(("colony PCR: junction 1 gives %d bp on the product; %d of %d by-product sizings are checks")
    :format(plan.pcr[1].product_sizes[1], checks, #plan.junctions * #plan.byproducts))
end)

case("sequencing_primers_read_across_each_junction_from_far_enough_back", function()
  reset()
  local _, ctx = ecori_ligation(1200, 2100, "insert")
  local plan = CP.build(ctx, {})
  local P = require("plasmid.bio.primer")
  local seq = plan.product.seq:tostring():upper()
  H.eq(#plan.junctions, #plan.sequencing)
  for _, r in ipairs(plan.sequencing) do
    H.ok(r.fwd or r.rev, tostring(r.error))
    -- each primer binds exactly where the plan says, and the junction falls
    -- past the unreadable first bases of the read
    if r.fwd then
      local sites = P.bind(seq, true, r.fwd.seq)
      H.eq(1, #sites, "a sequencing primer that binds twice is not a sequencing primer")
      H.eq(r.fwd.stop, sites[1].stop)
      H.eq(1, sites[1].strand)
      H.ok(r.fwd.into >= 50, ("the junction is only %d bases into the read"):format(r.fwd.into))
      H.ok(r.fwd.into <= 250, ("the junction is %d bases into the read"):format(r.fwd.into))
    end
    if r.rev then
      local sites = P.bind(seq, true, r.rev.seq)
      H.eq(1, #sites)
      H.eq(-1, sites[1].strand)
      H.eq(r.rev.start, sites[1].start)
      H.ok(r.rev.into >= 50 and r.rev.into <= 250, tostring(r.rev.into))
    end
    -- the junction the primers are named for is the one the plan recorded
    H.eq(plan.junctions[r.junction].at, r.at)
  end
  H.log(("sequencing primers: junction 1 read starts %d bases before the junction, junction 2 %d")
    :format(plan.sequencing[1].fwd.into, plan.sequencing[2].fwd.into))
end)

---------------------------------------------------------------------------
-- The command, the report buffer and the file
---------------------------------------------------------------------------

case("the_command_writes_the_plan_to_a_file_and_shows_the_same_content", function()
  reset()
  local _, ipath = open_fixture("synthetic/gg_insert.gb")
  local vbuf, vpath = open_fixture("synthetic/gg_vector.gb")
  C.run(vbuf, "goldengate", { ipath, vpath }, { enzyme = "BsmBI", first = true })
  local out = tmp .. "/gg_plan.txt"
  vim.cmd("DnaClonePlan " .. vim.fn.fnameescape(out))
  H.eq(1, vim.fn.filereadable(out), "no file was written")
  local text = H.read(out)
  local shown = report()
  -- the buffer and the file say the same thing, bar the buffer's key hint
  H.eq(nil, H.first_diff(text, shown:gsub("\n*%s*q close$", "") .. "\n"),
    H.first_diff(text, shown:gsub("\n*%s*q close$", "") .. "\n"))
  H.ok(text:find("Assembly plan — Golden Gate (BsmBI)", 1, true), text:sub(1, 200))
  -- the gel of the diagnostic digest is in both
  H.ok(text:find("ladder", 1, true) and text:find("━", 1, true), "no gel in the file")
  H.ok(shown:find("━", 1, true), "no gel in the report buffer")
  -- and the written plan is a plain text file that reads back line for line
  local back = vim.split(text, "\n", { plain = true })
  H.eq("", back[#back], "the file ends with a newline")
  H.log(("plan file: %s, %d lines, %d bytes"):format(vim.fn.fnamemodify(out, ":t"), #back - 1, #text))
end)

case("the_plan_command_says_so_when_there_is_no_assembly_to_plan", function()
  reset()
  local msgs = {}
  local orig = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  local ok = pcall(function() vim.cmd("DnaClonePlan") end)
  vim.notify = orig
  H.ok(ok, "the command errored instead of reporting")
  H.ok(msgs[1] and msgs[1]:find("no finished assembly to plan", 1, true), vim.inspect(msgs))
  H.log("no assembly: " .. (msgs[1] or ""))
end)

return T
