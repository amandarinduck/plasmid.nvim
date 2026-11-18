local H = require("helpers")
local detect = require("plasmid.bio.detect")
local alphabet = require("plasmid.model.alphabet")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local fdb = require("plasmid.featuredb")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

-- deterministic pseudo-random DNA: LuaJIT's math.random (Tausworthe, period
-- ~2^223) reseeded per stream. A power-of-two LCG is NOT usable here: its low
-- bits cycle quickly, so different seeds produce overlapping sequences.
local function random_dna(_, n)
  local b, out = { "A", "C", "G", "T" }, {}
  for i = 1, n do out[i] = b[math.random(4)] end
  return table.concat(out)
end
local function stream(seed)
  math.randomseed(seed)
  return function() return nil, seed end
end

local r = stream(42)
local L = 6000
local bg = random_dna(r, L)
local refs = {
  R1 = random_dna(r, 800), R2 = random_dna(r, 600), R3 = random_dna(r, 500), R4 = random_dna(r, 700),
  R5 = random_dna(r, 400), R6 = random_dna(r, 1000), R7 = random_dna(r, 15), R8 = random_dna(r, 900),
}
local function put(s, pos, piece) return s:sub(1, pos - 1) .. piece .. s:sub(pos + #piece) end
local q = bg
q = put(q, 1000, refs.R1)                                   -- exact, forward: 1000..1799
q = put(q, 2500, alphabet.revcomp(refs.R2))                 -- reverse: 2500..3099
local r4 = refs.R4:sub(1, 349) .. (refs.R4:sub(350, 350) == "A" and "C" or "A") .. refs.R4:sub(351)
q = put(q, 3300, r4)                                        -- 1 mismatch: 3300..3999
local r5 = refs.R5:sub(1, 199) .. refs.R5:sub(202)          -- 2 bp deletion (398 bp)
q = put(q, 4100, r5)                                        -- 4100..4497
q = put(q, 4600, refs.R6:sub(201, 800))                     -- 60% fragment: 4600..5199
q = put(q, 600, refs.R7)                                    -- 15 bp exact: 600..614
q = put(q, 5500, alphabet.revcomp(refs.R7))                 -- and reversed: 5500..5514
-- R3 across the origin: last 200 bp of the plasmid + first 300 bp
q = refs.R3:sub(201) .. q:sub(301, L - 200) .. refs.R3:sub(1, 200)
assert(#q == L)
-- R3 overwrites 1..300 and L-199..L; everything planted above lies outside those ranges

local lib = {
  { id = "R1", name = "alpha", type = "CDS", seq = refs.R1 },
  { id = "R2", name = "beta", type = "promoter", seq = refs.R2 },
  { id = "R3", name = "gamma", type = "rep_origin", seq = refs.R3 },
  { id = "R4", name = "delta", type = "terminator", seq = refs.R4 },
  { id = "R5", name = "epsilon", type = "misc_feature", seq = refs.R5 },
  { id = "R6", name = "zeta", type = "CDS", seq = refs.R6 },
  { id = "R7", name = "tag", type = "protein_bind", seq = refs.R7 },
  { id = "R8", name = "absent", type = "CDS", seq = refs.R8 },
  { id = "R9", name = "alpha-domain", type = "CDS", seq = refs.R1:sub(301, 450) }, -- same type, nested
  { id = "R10", name = "alpha-primer", type = "primer_bind", seq = refs.R1:sub(500, 524) }, -- other type, nested
}

local function by_name(hits)
  local out = {}
  for _, h in ipairs(hits) do
    out[h.ref.name] = out[h.ref.name] or {}
    table.insert(out[h.ref.name], h)
  end
  return out
end

case("planted_features_found_exactly", function()
  local hits = detect.search(q, lib, { circular = true })
  local n = by_name(hits)
  local function one(name, start, stop, strand, frag)
    H.ok(n[name], name .. " not found")
    H.eq(1, #n[name], name .. " count")
    local h = n[name][1]
    H.eq({ start, stop, strand, frag }, { h.start, h.stop, h.strand, h.fragment }, name)
    return h
  end
  H.eq(1, one("alpha", 1000, 1799, 1, false).identity)
  one("beta", 2500, 3099, -1, false)
  local g = one("gamma", L - 199, 300, 1, false)
  H.ok(g.wraps)
  local d = one("delta", 3300, 3999, 1, false)
  H.eq(699 / 700, d.identity)
  local e = one("epsilon", 4100, 4497, 1, false)
  H.ok(e.identity >= 0.99 and e.identity < 1, "gapped identity " .. e.identity)
  local z = one("zeta", 4600, 5199, 1, true)
  H.eq(0.6, z.coverage)
  H.eq(nil, n.absent)
  H.eq(2, #n.tag)
  table.sort(n.tag, function(a, b) return a.start < b.start end)
  H.eq({ 600, 614, 1 }, { n.tag[1].start, n.tag[1].stop, n.tag[1].strand })
  H.eq({ 5500, 5514, -1 }, { n.tag[2].start, n.tag[2].stop, n.tag[2].strand })
  -- nesting: same-type sub-part hidden, other-type part inside kept
  H.eq(nil, n["alpha-domain"])
  H.eq({ 1499, 1523 }, { n["alpha-primer"][1].start, n["alpha-primer"][1].stop })
  local with_nested = by_name(detect.search(q, lib, { circular = true, nested_same_type = true }))
  H.eq({ 1300, 1449 }, { with_nested["alpha-domain"][1].start, with_nested["alpha-domain"][1].stop })
end)

case("linear_record_does_not_wrap", function()
  local hits = by_name(detect.search(q, lib, { circular = false, fragments = false }))
  H.eq(nil, hits.gamma, "origin-spanning hit must not appear on a linear record")
  H.ok(hits.alpha and hits.beta)
end)

case("no_hits_on_unrelated_sequence", function()
  math.randomseed(7)
  local neg = random_dna(function() end, 20000)
  H.eq({}, detect.search(neg, lib, { circular = true }))
end)

case("import_fasta_with_descriptions_csv", function()
  local tmp = vim.fn.tempname()
  vim.fn.mkdir(tmp, "p")
  require("plasmid.config").options.featuredb_dir = tmp .. "/db"
  local fa = { }
  for _, x in ipairs(lib) do fa[#fa + 1] = ">" .. x.id .. " ignored header"; fa[#fa + 1] = x.seq end
  H.write = H.write or function(p, s) local f = assert(io.open(p, "wb")); f:write(s); f:close() end
  H.write(tmp .. "/parts.fasta", table.concat(fa, "\n") .. "\n")
  local csv = { "Accession,Feature,Type,Description" }
  for _, x in ipairs(lib) do csv[#csv + 1] = ('%s,"%s",%s,"a ""quoted"", note"'):format(x.id, x.name, x.type) end
  H.write(tmp .. "/parts.csv", table.concat(csv, "\n") .. "\n")
  local n, name = fdb.import(tmp .. "/parts.fasta")
  H.eq(#lib, n)
  H.eq("parts", name)
  local loaded = fdb.load("parts")
  H.eq({ "R1", "alpha", "CDS", 'a "quoted", note' }, { loaded[1].id, loaded[1].name, loaded[1].type, loaded[1].blurb })
  H.eq(refs.R1, loaded[1].seq)
  -- GenBank folder import: each annotated feature becomes a reference
  vim.fn.mkdir(tmp .. "/gbs", "p")
  H.write(tmp .. "/gbs/pbr.gb", H.read(H.fixture("ncbi/J01749.gb")))
  local n2, name2 = fdb.import(tmp .. "/gbs", "mine")
  -- every annotated feature of >= 8 bp (except source) becomes a reference
  local rec = assert(require("plasmid.io").read(H.fixture("ncbi/J01749.gb"))).records[1]
  local expected = 0
  for _, f in ipairs(rec.features) do
    if f.type ~= "source" and f.loc and f:length() >= 8 then expected = expected + 1 end
  end
  H.eq(expected, n2, tostring(name2))
  H.eq("mine", name2)
  local bla
  for _, x in ipairs(fdb.load("mine")) do if x.name == "bla" and x.type == "CDS" then bla = x end end
  H.eq("ATG", bla.seq:sub(1, 3))
  H.eq({ "mine", "parts" }, fdb.list())
  -- errors
  local bad, err = fdb.import(tmp .. "/nothing-here")
  H.eq(nil, bad)
  H.ok(err:match("cannot read") or err:match("no GenBank"), err)
  -- a real pLannotate BLAST export needs blastdbcmd; exercised when available
  if vim.fn.executable("blastdbcmd") == 1 and vim.env.PLASMID_PLANNOTATE_DIR then
    local n3 = fdb.import(vim.env.PLASMID_PLANNOTATE_DIR)
    H.ok(n3 and n3 > 1000)
  end
end)

case("review_flow_in_buffer", function()
  local tmp = vim.fn.tempname()
  vim.fn.mkdir(tmp, "p")
  require("plasmid.config").options.featuredb_dir = tmp .. "/db"
  require("plasmid.config").options.library_path = tmp .. "/none.gb"
  fdb.save("synthetic", lib)
  local fh = assert(io.open(tmp .. "/plasmid.fasta", "wb"))
  fh:write(">synthetic test plasmid [topology=circular]\n" .. q .. "\n")
  fh:close()
  vim.o.hidden = true
  vim.cmd("silent edit " .. tmp .. "/plasmid.fasta")
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  local rec = buffer.record(st)
  H.eq("circular", rec.topology)
  vim.cmd("DnaDetect")
  local proposals = vim.tbl_filter(function(f) return f.detected end, rec.overlays or {})
  H.eq(9, #proposals) -- alpha beta gamma delta epsilon zeta(fragment) tag x2 alpha-primer
  -- accept one with ga
  local alpha
  for _, f in ipairs(proposals) do if f.detected.ref.name == "alpha" then alpha = f end end
  require("plasmid.ui").goto_feature(buf, alpha.id)
  vim.api.nvim_feedkeys("ga", "x", false)
  local added
  for _, f in ipairs(rec.features) do if f:label() == "alpha" then added = f end end
  H.ok(added, "ga did not accept")
  H.eq("1000..1799", added.location_text)
  local notes = table.concat(added:get_all("note"), " | ")
  H.ok(notes:match("synthetic/R1, 100.0%% identity"), notes)
  local after = vim.tbl_filter(function(f) return f.detected end, rec.overlays or {})
  H.eq(8, #after, "accepted proposal leaves the overlay")
  vim.api.nvim_feedkeys("u", "x", false)
  H.eq(9, #vim.tbl_filter(function(f) return f.detected end, rec.overlays or {}), "undo restores the proposal")
  -- accept all full-length in one undoable step
  local n0 = #rec.features
  vim.cmd("DnaDetectAccept!")
  H.eq(n0 + 8, #rec.features, "all but the fragment accepted")
  vim.api.nvim_feedkeys("u", "x", false)
  H.eq(n0, #rec.features, "one undo reverts the whole batch")
  -- fragments on request; their label and note keep saying so
  vim.cmd("DnaDetectAccept! fragments")
  H.eq(n0 + 9, #rec.features, "fragments included")
  local frag
  for _, f in ipairs(rec.features) do if f:label():find("zeta", 1, true) then frag = f end end
  H.ok(frag, "zeta accepted")
  H.eq("zeta (fragment)", frag:label())
  local notes = {}
  for _, q in ipairs(frag.qualifiers) do if q.key == "note" then notes[#notes + 1] = q.value end end
  H.ok(table.concat(notes, " "):find("% of the reference", 1, true), "note records coverage")
  vim.api.nvim_feedkeys("u", "x", false)
  vim.cmd("DnaDetect off")
  H.eq(nil, rec.overlays)
  require("plasmid.config").options.featuredb_dir = nil
  require("plasmid.config").options.library_path = nil
end)

case("performance", function()
  local big = {}
  math.randomseed(99)
  for i = 1, 1500 do big[i] = { id = "x" .. i, name = "x" .. i, type = "misc_feature", seq = random_dna(function() end, math.random(100, 1600)) } end
  vim.list_extend(big, lib)
  local t0 = vim.uv.hrtime()
  local hits = detect.search(q, big, { circular = true })
  local ms = (vim.uv.hrtime() - t0) / 1e6
  io.stdout:write(("    [perf] detection: 6 kb plasmid vs %d references: %d hits in %.0f ms\n"):format(#big, #hits, ms))
  H.ok(ms < 5000)
end)

return T
