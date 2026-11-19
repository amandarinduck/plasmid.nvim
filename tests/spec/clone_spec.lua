-- The cloning UI: staging parts, the four commands, the product record and
-- the assembly report (including what it says when nothing assembles).
-- The engine itself is tested in assembly_spec.
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local C = require("plasmid.clone")
local A = require("plasmid.bio.assembly")
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
  vim.cmd("silent! only")
  vim.cmd("silent! %bwipeout!")
end

--- Open a fixture from a copy in the scratch directory, so products land
--- beside it instead of in the source tree.
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if st.view ~= "linear" then buffer.set_view(buf, "linear") end
  return buf, st, buffer.record(st)
end

local function write_rec(rec, stem)
  local p = ("%s/%s.gb"):format(tmp, stem)
  assert(pio.write(p, { rec }, "genbank"))
  return p
end
local function open_rec(rec, stem)
  local p = write_rec(rec, stem)
  vim.cmd("silent edit! " .. vim.fn.fnameescape(p))
  return vim.api.nvim_get_current_buf(), p
end

local function rnd(n)
  local b, t = { "A", "C", "G", "T" }, {}
  for i = 1, n do t[i] = b[math.random(4)] end
  return table.concat(t)
end
local function mk(name, seq, topo)
  return Record.new({ name = name, seq = seq, topology = topo or "linear", alphabet = "dna" })
end
local function seqof(rec) return rec.seq:tostring():upper() end
local function rotation_of(a, b)
  return #a == #b and (a:upper() .. a:upper()):find(b:upper(), 1, true) ~= nil
end

--- The text of the report buffer last shown.
local function report()
  local b = C.last
  H.ok(b and vim.api.nvim_buf_is_valid(b), "no report buffer")
  return table.concat(vim.api.nvim_buf_get_lines(b, 0, -1, false), "\n"), b
end

local function msgs_of(fn)
  local out = {}
  local orig = vim.notify
  vim.notify = function(m) out[#out + 1] = m end
  local ok, err = pcall(fn)
  vim.notify = orig
  if not ok then error(err, 0) end
  return out
end

--- A record digested in a visible digest table, cursor on fragment n.
local function digest_at(rel, names, n)
  local buf = open(rel)
  local db = D.run(buf, names)
  vim.api.nvim_set_current_win(vim.fn.bufwinid(db))
  vim.api.nvim_win_set_cursor(0, { D.state[db].first + (n or 1) - 1, 0 })
  return db, buf
end

-- Golden Gate parts: an insert flanked by outward-facing BsaI sites, and a
-- vector whose BsaI sites excise a stuffer and leave the matching overhangs.
local function gg_parts(ovl1, ovl2, body, seed)
  math.randomseed(seed or 11)
  local backbone, stuffer = rnd(900), rnd(60)
  local insert = mk("insert", "GGTCTC" .. "a" .. ovl1 .. body .. ovl2 .. "a" .. "GAGACC")
  local vector = mk("vec", ovl2 .. backbone .. ovl1 .. "a" .. "GAGACC" .. stuffer .. "GGTCTC" .. "a", "circular")
  return insert, vector, backbone
end

---------------------------------------------------------------------------
-- Staging
---------------------------------------------------------------------------

case("staging_the_record_the_selection_and_a_digest_fragment", function()
  reset()
  local lin = mk("linpart", rnd(400))
  local buf = open_rec(lin, "linpart")
  C.add(buf)
  H.eq(1, #C.staged)
  H.eq("linpart", C.staged[1].name)
  H.eq("blunt … blunt", C.ends_text(C.staged[1]))
  -- the selection, under a name of its own
  local sel = require("plasmid.select")
  sel.set(buf, sel.make(buffer.record(buffer.get(buf)), 10, 109))
  C.add(buf, "window")
  H.eq(2, #C.staged)
  H.eq("window", C.staged[2].name)
  H.eq(100, C.staged[2].rec:len())
  H.eq(lin.seq:sub(10, 109):upper(), seqof(C.staged[2].rec))
  sel.set(buf, nil)
  -- a fragment of a digest table, by cursor and by number
  local db = digest_at("ncbi/J01749.gb", { "EcoRI", "HindIII" }, 1)
  C.add(db)
  C.add(db, "2")
  H.eq(4, #C.staged)
  H.eq(2, #D.state[db].frags)
  H.ok(C.staged[3].frag and C.staged[3].frag.left.len == 4, "a digest fragment keeps its sticky ends")
  H.ok(C.ends_text(C.staged[3]):find("5′", 1, true), C.ends_text(C.staged[3]))
  -- a circular record has no ends until it is cut
  local circ = open_rec(mk("circpart", rnd(600), "circular"), "circpart")
  C.add(circ)
  H.eq(nil, C.staged[5].frag)
  H.eq("circular (no ends)", C.ends_text(C.staged[5]))
  local m = msgs_of(function() vim.cmd("DnaCloneClear") end)
  H.eq(0, #C.staged)
  H.ok(m[1]:find("cleared 5 staged parts", 1, true), m[1])
end)

case("the_staged_list_shows_name_length_and_ends", function()
  reset()
  local m = msgs_of(function() vim.cmd("DnaCloneList") end)
  H.ok(m[1]:find("no parts staged", 1, true), m[1])
  digest_at("ncbi/L09137.gb", { "EcoRI" }, 1)
  vim.cmd("DnaCloneAdd vector")
  vim.cmd("DnaCloneList")
  local text = report()
  H.ok(text:find("Staged parts (1)", 1, true), text)
  H.ok(text:find("vector", 1, true), text)
  H.ok(text:find("2690 bp", 1, true), text)       -- pUC19 2686 bp + the 4 base overhang
  H.ok(text:find("5′ AATT", 1, true), text)
  H.ok(text:find("digest of", 1, true), text)
end)

case("part_arguments_name_staged_parts_buffers_files_and_fragments", function()
  reset()
  local a = mk("alpha", rnd(300))
  local abuf, apath = open_rec(a, "alpha")
  C.add(abuf, "alpha")
  H.eq("alpha", C.resolve(abuf, "ALPHA").name)
  -- by open buffer: its record name, its file name
  local b = mk("beta", rnd(250))
  local bbuf = open_rec(b, "beta")
  H.eq("beta", C.resolve(bbuf, "beta").name)
  H.eq("beta", C.resolve(bbuf, "beta.gb").name)
  -- by file path
  H.eq(300, C.resolve(bbuf, apath).rec:len())
  -- the current selection
  local sel = require("plasmid.select")
  sel.set(bbuf, sel.make(buffer.record(buffer.get(bbuf)), 5, 54))
  H.eq(50, C.resolve(bbuf, "selection").rec:len())
  sel.set(bbuf, nil)
  H.eq(nil, C.resolve(bbuf, "selection"))
  -- a digest fragment number, only inside a digest table
  local db = digest_at("ncbi/J01749.gb", { "EcoRI", "HindIII" }, 1)
  H.ok(C.resolve(db, "2") ~= nil, "fragment 2")
  local _, why = C.resolve(db, "9")
  H.ok(why:find("not in this digest", 1, true), why)
  local _, why2 = C.resolve(abuf, "nosuchthing")
  H.ok(why2:find("no part called nosuchthing", 1, true), why2)
end)

---------------------------------------------------------------------------
-- Ligation
---------------------------------------------------------------------------

case("ligating_a_two_fragment_digest_restores_the_plasmid", function()
  reset()
  local db, buf = digest_at("ncbi/J01749.gb", { "EcoRI", "HindIII" }, 1)
  local pbr = buffer.record(buffer.get(buf))
  local nb, rec, ctx = C.run(db, "ligation", {})
  H.eq("ligation", ctx.method)
  H.eq(1, ctx.nsolutions, "the two different overhangs fix one order")
  H.eq(pbr:len(), rec:len())
  H.eq("circular", rec.topology)
  H.ok(rotation_of(seqof(pbr), seqof(rec)), "the product is pBR322 again")
  H.ok(nb and buffer.get(nb).modified, "the product opens as an unsaved record")
  local text = report()
  H.ok(text:find("Ligation: 2 parts →", 1, true), text)
  H.ok(text:find("4361 bp circular", 1, true), text)
  for _, want in ipairs({ " parts", " junctions", " features carried", "5′ AATT", "5′ AGCT", "forward" }) do
    H.ok(text:find(want, 1, true), ("report has no %q:\n%s"):format(want, text))
  end
  -- every feature of the plasmid is in the product, the cut ones re-joined
  local want = {}
  for _, f in ipairs(pbr.features) do
    if f.type ~= "source" then want[f.type .. "|" .. f:label() .. "|" .. pbr:extract(f)] = true end
  end
  local missing = 0
  for k in pairs(want) do
    local found = false
    for _, g in ipairs(rec.features) do
      if k == g.type .. "|" .. g:label() .. "|" .. rec:extract(g) then found = true break end
    end
    if not found then missing = missing + 1 end
  end
  H.eq(0, missing, "every feature carried over")
  H.log(("religated pBR322: %d bp, %d features in the product"):format(rec:len(), #rec.features))
end)

case("the_product_is_unsaved_and_round_trips_through_genbank", function()
  reset()
  local db = digest_at("ncbi/L09137.gb", { "EcoRI" }, 1)
  local nb, rec = C.run(db, "ligation", {})
  H.ok(buffer.get(nb).modified, "not on disk yet")
  vim.api.nvim_set_current_buf(nb)
  vim.cmd("silent write")
  local path = vim.api.nvim_buf_get_name(nb)
  local back = assert(pio.read(path)).records[1]
  H.eq(seqof(rec), seqof(back))
  H.eq("circular", back.topology)
  H.eq(#rec.features, #back.features)
  H.ok(not buffer.get(nb).modified, "saved")
  H.log(("product round trip: %s, %d bp, %d features"):format(vim.fn.fnamemodify(path, ":t"), back:len(), #back.features))
end)

case("parts_and_junctions_are_annotated_in_the_product", function()
  reset()
  local db = digest_at("ncbi/J01749.gb", { "EcoRI", "HindIII" }, 1)
  local _, rec = C.run(db, "ligation", {})
  local parts, junctions = {}, {}
  for _, f in ipairs(rec.features) do
    local note = f:get("note") or ""
    if note:find("^part %d+ of ") then parts[#parts + 1] = f:label() end
    if note:find("junction: ") then junctions[#junctions + 1] = f:label() end
  end
  H.eq(2, #parts)
  H.eq(2, #junctions)
  H.ok(vim.tbl_contains(junctions, "AATT") and vim.tbl_contains(junctions, "AGCT"), vim.inspect(junctions))
  -- the annotations can be switched off
  require("plasmid.config").options.cloning.annotate_parts = false
  require("plasmid.config").options.cloning.annotate_junctions = false
  local _, plain = C.run(db, "ligation", {})
  require("plasmid.config").setup({})
  local extra = 0
  for _, f in ipairs(plain.features) do
    local note = f:get("note") or ""
    if note:find("^part %d+ of ") or note:find("junction: ") then extra = extra + 1 end
  end
  H.eq(0, extra, "annotate_parts/annotate_junctions = false leaves only the carried features")
end)

case("a_snapgene_plasmid_religates_with_its_features", function()
  reset()
  local buf = open("snapgene/pFA-KanMX4.dna")
  local rec = buffer.record(buffer.get(buf))
  local single
  for _, e in ipairs(require("plasmid.data.enzymes")) do
    local frags = D.compute(rec, { e.name })
    if #frags == 1 and not frags[1].uncut then single = e.name break end
  end
  H.ok(single, "a unique cutter in pFA-KanMX4")
  local db = D.run(buf, { single })
  local _, prod = C.run(db, "ligation", {})
  H.eq(rec:len(), prod:len())
  H.eq("circular", prod.topology)
  H.ok(rotation_of(seqof(rec), seqof(prod)), "the religated plasmid is the original")
  local before, after = 0, 0
  for _, f in ipairs(rec.features) do if f.type ~= "source" then before = before + 1 end end
  for _, f in ipairs(prod.features) do
    local note = f:get("note") or ""
    if f.type ~= "source" and not note:find("^part %d+ of ") and not note:find("junction: ") then after = after + 1 end
  end
  H.eq(before, after, "SnapGene features survive the round trip through the assembler")
  H.log(("pFA-KanMX4 cut with %s and religated: %d bp, %d features"):format(single, prod:len(), after))
end)

---------------------------------------------------------------------------
-- Golden Gate
---------------------------------------------------------------------------

case("golden_gate_from_files_builds_and_names_the_product", function()
  reset()
  local body = rnd(300)
  local insert, vector, backbone = gg_parts("AATG", "GCTT", body, 11)
  local ipath, vpath = write_rec(insert, "gg_insert"), write_rec(vector, "gg_vector")
  local buf = open_rec(mk("bench", rnd(100)), "gg_bench")
  local nb, rec, ctx = C.run(buf, "goldengate", { ipath, vpath }, { enzyme = "BsaI" })
  H.eq("goldengate", ctx.method)
  H.eq("circular", rec.topology)
  H.eq("insert_BsaI_vec", rec.name)
  H.ok(rotation_of("GCTT" .. backbone .. "AATG" .. body, seqof(rec)), "insert joined to backbone")
  H.ok(not seqof(rec):find("GGTCTC"), "no BsaI site is left")
  H.ok(rec.description:find("Golden Gate (BsaI) of insert + vec", 1, true), rec.description)
  H.ok(nb ~= nil, "a product buffer")
  local text = report()
  H.ok(text:find("Golden Gate (BsaI): 2 parts →", 1, true), text)
  H.ok(text:find("5′ AATG", 1, true) and text:find("5′ GCTT", 1, true), text)
  H.ok(text:find("insert → vec", 1, true), text)
  -- the same thing through the command, with the enzyme as the first argument
  reset()
  open_rec(mk("bench2", rnd(100)), "gg_bench2")
  vim.cmd(("DnaGoldenGate BsaI %s %s"):format(vim.fn.fnameescape(ipath), vim.fn.fnameescape(vpath)))
  H.ok(report():find("Golden Gate (BsaI)", 1, true), report())
end)

case("golden_gate_part_problems_reach_the_report", function()
  reset()
  local _, vector = gg_parts("AATG", "GCTT", rnd(120), 11)
  local vpath = write_rec(vector, "gg_vector2")
  local buf = open_rec(mk("plain", rnd(300)), "gg_plain")
  local msgs = msgs_of(function() C.run(buf, "goldengate", { "plain", vpath }, { enzyme = "BsaI" }) end)
  local text = report()
  H.ok(text:find("nothing assembled", 1, true), text)
  H.ok(text:find("plain has no BsaI site", 1, true), text)
  H.ok(text:find("warnings", 1, true), text)
  H.ok(msgs[#msgs]:find("1 of 2 parts could not be prepared (plain)", 1, true), msgs[#msgs])
  H.ok(text:find("plain contributed no fragment", 1, true), text)
  -- a palindromic overhang is reported too
  reset()
  local pbuf = open_rec(mk("pal", "GGTCTC" .. "a" .. "AATT" .. rnd(100) .. "GCTT" .. "a" .. "GAGACC"), "gg_pal")
  C.run(pbuf, "goldengate", { "pal", vpath }, { enzyme = "BsaI" })
  H.ok(report():find("palindromic", 1, true), report())
  -- an enzyme that is not Type IIS
  reset()
  local b2 = open_rec(mk("bench3", rnd(100)), "gg_bench3")
  msgs_of(function() C.run(b2, "goldengate", { vpath }, { enzyme = "EcoRI" }) end)
  H.ok(report():find("Type IIS", 1, true), report())
end)

---------------------------------------------------------------------------
-- Gibson
---------------------------------------------------------------------------

case("gibson_joins_staged_parts_by_homology", function()
  reset()
  math.randomseed(23)
  local h1, h2, h3 = rnd(25), rnd(30), rnd(22)
  local a1, a2, a3 = rnd(400), rnd(500), rnd(300)
  for _, p in ipairs({ { "gA", h1 .. a1 .. h2 }, { "gB", h2 .. a2 .. h3 }, { "gC", h3 .. a3 .. h1 } }) do
    open_rec(mk(p[1], p[2]), p[1])
    vim.cmd("DnaCloneAdd")
  end
  H.eq(3, #C.staged)
  local nb, rec, ctx = C.run(vim.api.nvim_get_current_buf(), "gibson", {})
  H.eq("gibson", ctx.method)
  H.eq("circular", rec.topology)
  H.eq("gA_gB_gC_gibson", rec.name)
  H.ok(rotation_of(h1 .. a1 .. h2 .. a2 .. h3 .. a3, seqof(rec)), "each overlap is written once")
  H.ok(nb ~= nil)
  local text = report()
  H.ok(text:find("Gibson/HiFi: 3 parts →", 1, true), text)
  H.ok(text:find("overlap 25 bp", 1, true), text)
  H.ok(text:find("overlap 30 bp", 1, true), text)
  H.ok(text:find("gA → gB", 1, true), text)
  local overlaps = 0
  for _, f in ipairs(rec.features) do
    if (f:get("note") or ""):find("Gibson junction", 1, true) then overlaps = overlaps + 1 end
  end
  H.eq(3, overlaps)
  H.log(("Gibson of 3 staged parts: %d bp, overlaps 25/30/22 bp"):format(rec:len()))
end)

case("gibson_failure_names_the_homology_floor", function()
  reset()
  math.randomseed(29)
  local h = rnd(12)
  local x = open_rec(mk("x", h .. rnd(300) .. h), "gib_x")
  vim.cmd("DnaCloneAdd")
  open_rec(mk("y", h .. rnd(200) .. h), "gib_y")
  vim.cmd("DnaCloneAdd")
  local msgs = msgs_of(function() vim.cmd("DnaGibson") end)
  local text = report()
  H.ok(text:find("Gibson/HiFi: 2 parts — nothing assembled", 1, true), text)
  H.ok(text:find("no terminal homology of 20 bp or more", 1, true), text)
  H.ok(text:find("lower cloning.min_homology", 1, true), text)
  H.ok(msgs[#msgs]:find("do not assemble", 1, true), msgs[#msgs])
  -- lowering the floor makes the same parts assemble
  require("plasmid.config").options.cloning.min_homology = 10
  local _, rec = C.run(x, "gibson", {})
  require("plasmid.config").setup({})
  H.ok(rec ~= nil and rec:len() > 0, "12 bp overlaps join once the floor is 10")
end)

---------------------------------------------------------------------------
-- :DnaClone, pickers, failures
---------------------------------------------------------------------------

case("dna_clone_infers_the_method", function()
  reset()
  local insert, vector = gg_parts("AATG", "GCTT", rnd(200), 11)
  H.eq("goldengate", C.infer({ { name = "i", rec = insert }, { name = "v", rec = vector } }))
  math.randomseed(31)
  local h = rnd(28)
  local g1, g2 = mk("h1", h .. rnd(300) .. h), mk("h2", h .. rnd(200) .. h)
  H.eq("gibson", C.infer({
    { name = "h1", rec = g1, frag = A.from_record(g1, "h1") },
    { name = "h2", rec = g2, frag = A.from_record(g2, "h2") },
  }))
  local pbr = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  local frags = A.cut(pbr, { "EcoRI", "HindIII" })
  H.eq("ligation", C.infer({
    { name = frags[1].name, rec = frags[1].rec, frag = frags[1] },
    { name = frags[2].name, rec = frags[2].rec, frag = frags[2] },
  }))
  -- end to end: :DnaClone on a digest picks ligation and rebuilds the plasmid
  local _, buf = digest_at("ncbi/J01749.gb", { "EcoRI", "HindIII" }, 1)
  local want = buffer.record(buffer.get(buf)):len()
  vim.cmd("DnaClone")
  local text = report()
  H.ok(text:find("Ligation: 2 parts →", 1, true), text)
  H.eq(want, buffer.record(buffer.get(vim.api.nvim_get_current_buf())):len())
end)

case("several_assemblies_offer_a_picker", function()
  reset()
  digest_at("ncbi/J01749.gb", { "EcoRV", "ScaI" }, 1)
  local orig = vim.ui.select
  local shown
  vim.ui.select = function(items, o, cb) shown = { items = items, o = o }; cb(items[1]) end
  local ok, err = pcall(function() vim.cmd("DnaLigate") end)
  vim.ui.select = orig
  H.ok(ok, tostring(err))
  H.eq(2, #shown.items, "two blunt fragments join either way round")
  H.ok(shown.o.prompt:find("2 possible assemblies", 1, true), shown.o.prompt)
  local line = shown.o.format_item(shown.items[1])
  H.ok(line:find("circular", 1, true) and line:find("bp", 1, true), line)
  H.ok(line:find("▸", 1, true) or line:find("◂", 1, true), line)
  local text = report()
  H.ok(text:find("2 assemblies are possible", 1, true), text)
  H.ok(text:find("warnings", 1, true), text)
  -- cancelling the picker builds nothing
  reset()
  local db2 = digest_at("ncbi/J01749.gb", { "EcoRV", "ScaI" }, 1)
  local before = #vim.api.nvim_list_bufs()
  vim.ui.select = function(_, _, cb) cb(nil) end
  C.run(db2, "ligation", {})
  vim.ui.select = orig
  H.eq(before, #vim.api.nvim_list_bufs(), "nothing opened")
end)

case("a_ligation_that_cannot_close_lists_every_end", function()
  reset()
  local db = digest_at("ncbi/J01749.gb", { "EcoRI" }, 1)
  C.add(db, "eco")
  local db2 = digest_at("ncbi/L09137.gb", { "PstI" }, 1)
  C.add(db2, "pst")
  H.eq(2, #C.staged)
  local msgs = msgs_of(function() vim.cmd("DnaLigate") end)
  local text = report()
  H.ok(text:find("Ligation: 2 parts — nothing assembled", 1, true), text)
  H.ok(text:find(" joins", 1, true), text)
  H.ok(text:find("5′ AATT", 1, true) and text:find("3′ TGCA", 1, true), text)
  H.ok(text:find("eco → pst", 1, true) and text:find("pst → eco", 1, true), text)
  H.ok(text:find("(flipped)", 1, true), text)
  H.ok(text:find("a 5′ overhang cannot anneal to a 3′ overhang", 1, true), text)
  H.ok(text:find("no two parts have ends that anneal", 1, true), text)
  H.ok(msgs[#msgs]:find("the report lists every end and why each join failed", 1, true), msgs[#msgs])
  -- ends of different lengths are named as such
  local blunt = A.from_record(mk("blunty", rnd(200)), "blunty")
  local sticky = A.cut(assert(pio.read(H.fixture("ncbi/L09137.gb"))).records[1], { "EcoRI" }, "sticky")[1]
  local joined = {}
  for _, r in ipairs(C.join_rows({ blunt, sticky }, "ligation")) do joined[#joined + 1] = r[1] .. " = " .. r[2] end
  local all = table.concat(joined, "\n")
  H.ok(all:find("blunt end and a 4 base end", 1, true), all)
end)

case("a_circular_part_cannot_be_ligated", function()
  reset()
  open_rec(mk("intact", rnd(800), "circular"), "intact")
  vim.cmd("DnaCloneAdd")
  local lin = open_rec(mk("stick", rnd(300)), "stick")
  vim.cmd("DnaCloneAdd")
  local msgs = msgs_of(function() C.run(lin, "ligation", {}) end)
  local text = report()
  H.ok(text:find("Ligation: 2 parts — nothing assembled", 1, true), text)
  H.ok(text:find("intact is circular; cut it first", 1, true), text)
  H.ok(msgs[#msgs]:find("could not be prepared (intact)", 1, true), msgs[#msgs])
  -- the staged circle is still listed, with no ends of its own
  H.ok(text:find("intact", 1, true) and text:find("  —  ", 1, true), text)
end)

case("unknown_parts_and_an_empty_stage_are_refused", function()
  reset()
  open_rec(mk("lonely", rnd(300)), "lonely")
  local m1 = msgs_of(function() vim.cmd("DnaClone nosuchpart") end)
  H.ok(m1[#m1]:find("no part called nosuchpart", 1, true), m1[#m1])
  H.eq(nil, C.last, "a typo does not open a report")
  -- no staged parts, and not a digest table
  vim.cmd("silent! only")
  vim.cmd("silent enew")
  local m2 = msgs_of(function() vim.cmd("DnaLigate") end)
  H.ok(m2[#m2]:find("stage them with :DnaCloneAdd", 1, true), m2[#m2])
  -- :DnaCloneAdd somewhere that holds no sequence
  local m3 = msgs_of(function() vim.cmd("DnaCloneAdd") end)
  H.ok(m3[#m3]:find("nothing to stage here", 1, true), m3[#m3])
end)

case("config_defaults_and_the_staging_keymap", function()
  reset()
  local c = require("plasmid.config").options.cloning
  H.eq(20, c.min_homology)
  H.eq(200, c.max_homology)
  H.eq("BsaI", c.enzyme)
  H.eq(true, c.annotate_parts)
  H.eq(true, c.annotate_junctions)
  H.eq(50, c.max_solutions)
  H.eq("gA", require("plasmid.config").options.keymaps.clone_add)
  local buf = open_rec(mk("mapped", rnd(200)), "mapped")
  local found
  for _, m in ipairs(vim.api.nvim_buf_get_keymap(buf, "n")) do if m.lhs == "gA" then found = m end end
  H.ok(found, "gA is mapped in a plasmid buffer")
  vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes("gA", true, false, true), "x", false)
  H.eq(1, #C.staged)
  H.eq("mapped", C.staged[1].name)
  local cmds = vim.api.nvim_get_commands({})
  for _, n in ipairs({ "DnaClone", "DnaLigate", "DnaGoldenGate", "DnaGibson",
    "DnaCloneAdd", "DnaCloneList", "DnaCloneClear" }) do
    H.ok(cmds[n], n .. " is not registered")
  end
  reset()
end)

return T
