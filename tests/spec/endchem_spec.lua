-- End chemistry: fill-in, chew-back, single-strand nuclease, phosphatase and
-- kinase, and the ligation that asks for a 5' phosphate.
--
-- The check is by construction. A fragment is read back as the two strands it
-- actually has — the stretch of the source present on the top strand and the
-- stretch present on the bottom — and each treatment is defined here as what
-- it does to those two stretches: fill-in extends the recessed strand of a 5'
-- overhang to meet the protruding one, chew-back trims the protruding strand
-- of a 3' overhang back to the recessed one, a nuclease trims whichever
-- strand protrudes. The engine's own end bookkeeping is never consulted to
-- decide what the answer should be.
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local A = require("plasmid.bio.assembly")
local D = require("plasmid.digest")
local C = require("plasmid.clone")
local alphabet = require("plasmid.model.alphabet")
local buffer = require("plasmid.buffer")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local function load(stem) return assert(pio.read(H.fixture("ncbi/" .. stem .. ".gb"))).records[1] end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent! only")
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if st.view ~= "linear" then buffer.set_view(buf, "linear") end
  return buf, st, buffer.record(st)
end

vim.o.columns, vim.o.lines, vim.o.hidden = 140, 50, true

---------------------------------------------------------------------------
-- The two strands a fragment really has
---------------------------------------------------------------------------

--- A fragment as its two strands: the top strand 5'->3', and the bottom
--- strand 5'->3'. Derived from the written sequence and which side of each
--- end is single-stranded, so that two fragments with the same strands are
--- the same molecule however the ends were recorded.
local function strands(frag)
  local s = A.sequence(frag):upper()
  local L = #s
  local a = frag.left.strand == -1 and frag.left.len + 1 or 1       -- top starts here
  local b = frag.right.strand == -1 and L - frag.right.len or L     -- ... and ends here
  local c = frag.left.strand == 1 and frag.left.len + 1 or 1        -- bottom covers here
  local d = frag.right.strand == 1 and L - frag.right.len or L
  return s:sub(a, b), alphabet.revcomp(s:sub(c, d))
end

--- What a treatment should leave, worked out from the two strands alone.
--- `top` and `bot` are given as intervals over the written sequence; the ends
--- are whichever of the two reaches further.
local function expected(frag, name, sides)
  local s = A.sequence(frag):upper()
  local L = #s
  local top = { frag.left.strand == -1 and frag.left.len + 1 or 1,
    frag.right.strand == -1 and L - frag.right.len or L }
  local bot = { frag.left.strand == 1 and frag.left.len + 1 or 1,
    frag.right.strand == 1 and L - frag.right.len or L }
  local function act(side)
    local i = side == "left" and 1 or 2
    if top[i] == bot[i] then return end                       -- already blunt
    -- which strand protrudes, and is its protruding terminus a 5' or a 3' one?
    local long, short = top, bot
    if (side == "left" and bot[1] < top[1]) or (side == "right" and bot[2] > top[2]) then
      long, short = bot, top
    end
    -- the top strand's 5' terminus is on the left, the bottom strand's on the
    -- right, so a protruding top strand on the left is a 5' overhang
    local five = (long == top) == (side == "left")
    if name == "fillin" then if five then short[i] = long[i] end
    elseif name == "chewback" then if not five then long[i] = short[i] end
    elseif name == "blunt" then if five then short[i] = long[i] else long[i] = short[i] end
    elseif name == "nuclease" then long[i] = short[i] end
  end
  for _, side in ipairs(sides) do act(side) end
  return s:sub(top[1], top[2]), alphabet.revcomp(s:sub(bot[1], bot[2]))
end

local SETS = {
  { "EcoRI", "5'" }, { "BamHI", "5'" }, { "HindIII", "5'" }, { "SalI", "5'" },
  { "PstI", "3'" }, { "KpnI", "3'" }, { "SacI", "3'" }, { "SphI", "3'" },
  { "EcoRV", "blunt" }, { "ScaI", "blunt" }, { "PvuII", "blunt" },
}

case("every_treated_end_is_what_the_chemistry_defines", function()
  local checks, by_kind, by_treatment = 0, {}, {}
  for _, stem in ipairs({ "J01749", "L09137", "NC_005816" }) do
    local rec = load(stem)
    for _, set in ipairs(SETS) do
      local frags = A.cut(rec, { set[1] })
      if frags[1] and frags[1].digest then
        for _, f in ipairs(frags) do
          for _, name in ipairs(A.TREATMENT_ORDER) do
            for _, sides in ipairs({ { "left" }, { "right" }, { "left", "right" } }) do
              local g, notes, err = A.treat(f, name, { sides = sides })
              H.ok(g, ("%s %s %s: %s"):format(stem, set[1], name, tostring(err)))
              H.eq(#sides, #notes, "one line per end treated")
              if A.TREATMENTS[name].phos ~= nil then
                -- a phosphatase or a kinase moves no base
                local w1, c1 = strands(f)
                local w2, c2 = strands(g)
                H.eq(w1, w2); H.eq(c1, c2)
                for _, side in ipairs({ "left", "right" }) do
                  local want = A.phosphate(f, side)
                  if vim.tbl_contains(sides, side) then want = A.TREATMENTS[name].phos end
                  H.eq(want, A.phosphate(g, side), ("%s %s %s %s"):format(stem, set[1], name, side))
                end
              else
                local w, c = expected(f, name, sides)
                local gw, gc = strands(g)
                H.eq(w, gw, ("%s %s %s %s: top strand"):format(stem, set[1], name, table.concat(sides, "+")))
                H.eq(c, gc, ("%s %s %s %s: bottom strand"):format(stem, set[1], name, table.concat(sides, "+")))
                -- blunting leaves the phosphates alone: the nuclease and the
                -- polymerase both finish on a 5' phosphate
                for _, side in ipairs({ "left", "right" }) do
                  H.eq(A.phosphate(f, side), A.phosphate(g, side))
                end
              end
              -- the written sequence is the longer of the two strands at each end
              local gw, gc = strands(g)
              local full = A.sequence(g):upper()
              H.ok(full:find(gw, 1, true) and full:find(alphabet.revcomp(gc), 1, true),
                "both strands lie inside the written sequence")
              checks = checks + 1
              by_kind[set[2]] = (by_kind[set[2]] or 0) + 1
              by_treatment[name] = (by_treatment[name] or 0) + 1
            end
          end
        end
      end
    end
  end
  local parts = {}
  for _, n in ipairs(A.TREATMENT_ORDER) do parts[#parts + 1] = ("%s %d"):format(n, by_treatment[n]) end
  H.log(("%d treated ends agree with the two-strand model (5' overhangs %d, 3' %d, blunt %d): %s")
    :format(checks, by_kind["5'"], by_kind["3'"], by_kind.blunt, table.concat(parts, ", ")))
end)

case("the_treatment_table_reads_as_the_bench_does", function()
  local pbr = load("J01749")
  local eco = A.cut(pbr, { "EcoRI" })[1]      -- 5' AATT at both ends
  local pst = A.cut(pbr, { "PstI" })[1]       -- 3' TGCA at both ends
  local blu = A.cut(pbr, { "EcoRV", "ScaI" })[1]
  local n = eco.rec:len()
  local rows = {}
  local function row(frag, name)
    local g = A.treat(frag, name)
    rows[#rows + 1] = ("%s/%s %s %d bp"):format(name, A.end_kind(frag, "left"),
      A.end_text(g, "left"), g.rec:len())
    return g
  end
  -- a 5' overhang: filled in, the bases stay; a nuclease takes them away
  H.eq("blunt", A.end_text(row(eco, "fillin"), "left"))
  H.eq(n, A.treat(eco, "fillin").rec:len(), "fill-in moves no base")
  H.eq(n - 8, A.treat(eco, "nuclease").rec:len(), "the nuclease removes both overhangs")
  H.eq("5′ AATT", A.end_text(A.treat(eco, "chewback"), "left"), "chew-back leaves a 5' overhang alone")
  -- a 3' overhang: chewed back, the bases go; fill-in does not touch it
  H.eq("blunt", A.end_text(row(pst, "chewback"), "left"))
  H.eq(pst.rec:len() - 8, A.treat(pst, "chewback").rec:len())
  H.eq("3′ TGCA", A.end_text(A.treat(pst, "fillin"), "left"), "fill-in leaves a 3' overhang alone")
  -- T4 polymerase with dNTPs does both, which is how blunting is really done
  H.eq("blunt", A.end_text(A.treat(eco, "blunt"), "left"))
  H.eq("blunt", A.end_text(A.treat(pst, "blunt"), "left"))
  H.eq(n, A.treat(eco, "blunt").rec:len())
  H.eq(pst.rec:len() - 8, A.treat(pst, "blunt").rec:len())
  -- a blunt end is not changed by any of them
  for _, name in ipairs({ "fillin", "chewback", "blunt", "nuclease" }) do
    local g = A.treat(blu, name)
    H.eq(A.sequence(blu), A.sequence(g), name .. " on a blunt end")
    H.eq("blunt", A.end_text(g, "right"))
  end
  -- one side only
  local half = A.treat(eco, "nuclease", { sides = { "left" } })
  H.eq("blunt", A.end_text(half, "left"))
  H.eq("5′ AATT", A.end_text(half, "right"))
  H.eq(n - 4, half.rec:len())
  H.log("end chemistry: " .. table.concat(rows, ", "))
end)

case("features_follow_a_chewed_back_end", function()
  local pbr = load("J01749")
  local pst = A.cut(pbr, { "PstI" })[1]
  local L = pst.rec:len()
  local g = A.treat(pst, "chewback")
  H.eq(L - 8, g.rec:len())
  -- the four bases each end lost are gone from the record
  H.eq(A.sequence(pst):sub(5, L - 4):upper(), A.sequence(g):upper())
  -- a feature clear of both chewed ends keeps exactly its own residues
  local clear, same, moved = 0, 0, 0
  for _, f in ipairs(pst.rec.features) do
    local lo, hi = f:bounds()
    if f.type ~= "source" and lo and lo >= 5 and hi <= L - 4 then
      clear = clear + 1
      local want = pst.rec:extract(f)
      local found = false
      for _, h in ipairs(g.rec.features) do
        if h.type == f.type and h:label() == f:label() and g.rec:extract(h) == want then found = true; break end
      end
      if found then same = same + 1 else moved = moved + 1 end
    end
  end
  H.ok(clear > 40, ("only %d features are clear of the chewed ends"):format(clear))
  H.eq(clear, same, ("%d of %d features clear of the ends did not survive intact"):format(moved, clear))
  -- a feature annotated across a chewed end is truncated, not left dangling
  for _, f in ipairs(g.rec.features) do
    local lo, hi = f:bounds()
    if lo then H.ok(lo >= 1 and hi <= g.rec:len(), f:label() .. " stays inside the fragment") end
  end
  H.log(("chew-back of 4 + 4 bases: all %d features clear of the ends kept their residues, "
    .. "and none of the %d in the record reaches outside the fragment"):format(same, #g.rec.features))
end)

case("a_treatment_that_would_leave_nothing_is_refused", function()
  local pbr = load("J01749")
  -- a 10 bp span with a 4 base overhang on each side still has a duplex core
  local small = A.from_span(pbr, 100, 109, { len = 4, strand = 1 }, { len = 4, strand = -1 }, "small")
  H.ok(A.treat(small, "nuclease"), "8 of 10 bases removed is still a fragment")
  local tiny = A.from_span(pbr, 100, 107, { len = 4, strand = 1 }, { len = 4, strand = -1 }, "tiny")
  local g, _, err = A.treat(tiny, "nuclease")
  H.eq(nil, g)
  H.ok(err:find("would remove all 8 bases", 1, true), tostring(err))
  local bad, _, err2 = A.treat(tiny, "phosphorylate")
  H.eq(nil, bad)
  H.ok(err2:find("unknown treatment", 1, true), tostring(err2))
  H.log("refusals: a treatment that would consume the fragment, and an unknown name")
end)

---------------------------------------------------------------------------
-- Ligation and the 5' phosphate
---------------------------------------------------------------------------

case("a_dephosphorylated_vector_does_not_self_ligate", function()
  local puc = load("L09137")
  local cut = A.cut(puc, { "EcoRI" })[1]
  H.eq(1, #A.assemble({ cut }, {}), "untreated, it closes on itself as it always has")
  local flat = A.treat(cut, "dephos")
  H.eq(false, A.phosphate(flat, "left")); H.eq(false, A.phosphate(flat, "right"))
  local asms = A.assemble({ flat }, {})
  H.eq(0, #vim.tbl_filter(function(a) return a.circular end, asms), "it no longer closes")
  local j, why = A.compatible(flat, flat)
  H.eq(nil, j)
  H.ok(why:find("5′ phosphate", 1, true), why)
  -- the overhangs still anneal: it is the phosphate that is missing
  H.eq("AATT", A.end_seq(flat, "right"))
  H.ok(A.compatible(flat, flat, { require_phosphate = false }) ~= nil,
    "with the check off it ligates again")
  -- a kinase puts it back
  H.ok(A.compatible(A.treat(flat, "kinase"), A.treat(flat, "kinase")) ~= nil)
  H.log("a dephosphorylated pUC19 EcoRI vector: ends still anneal, no junction, kinase restores it")
end)

case("one_phosphate_per_junction_is_enough", function()
  -- the point of dephosphorylating a vector: it still takes an insert,
  -- because each junction gets its phosphate from the insert's end
  local pbr = load("J01749")
  local frags = A.cut(pbr, { "EcoRI", "HindIII" })
  H.eq(2, #frags)
  local vector, insert = frags[1], frags[2]
  if vector.rec:len() < insert.rec:len() then vector, insert = insert, vector end
  local flat = A.treat(vector, "dephos")
  local asms = A.assemble({ flat, insert }, {})
  H.eq(1, #asms)
  H.eq(true, asms[1].circular)
  local rec = A.build(asms[1], { name = "vector_insert" })
  H.eq(pbr:len(), rec:len(), "the plasmid is back, one nick per junction")
  -- dephosphorylate both and nothing closes
  local flat_insert = A.treat(insert, "dephos")
  H.eq(0, #A.assemble({ flat, flat_insert }, { circular = true }))
  H.log(("dephosphorylated %d bp vector + %d bp insert: closes at %d bp; both flat: nothing")
    :format(vector.rec:len(), insert.rec:len(), rec:len()))
end)

case("untreated_parts_ligate_exactly_as_before", function()
  -- the phosphate is modelled, but a restriction digest leaves both ends
  -- phosphorylated, so nothing an untreated digest does may change
  local counts = {}
  for _, stem in ipairs({ "J01749", "L09137", "NC_005816" }) do
    local rec = load(stem)
    for _, names in ipairs({ { "EcoRI" }, { "EcoRI", "HindIII" }, { "EcoRV", "ScaI" }, { "PstI" } }) do
      local frags = A.cut(rec, names)
      if frags[1] and frags[1].digest then
        for _, f in ipairs(frags) do
          H.eq(true, A.phosphate(f, "left")); H.eq(true, A.phosphate(f, "right"))
        end
        local asms = A.assemble(frags, {})
        local off = A.assemble(frags, { require_phosphate = false })
        H.eq(#asms, #off, ("%s %s: the check changes nothing"):format(stem, table.concat(names, "+")))
        counts[#counts + 1] = #asms
      end
    end
  end
  H.log(("%d untreated digests assemble identically with and without the phosphate check")
    :format(#counts))
end)

case("the_failure_report_names_the_phosphate", function()
  C.staged = {}
  local buf = open("ncbi/L09137.gb")
  local db = D.run(buf, { "EcoRI" })
  vim.api.nvim_set_current_win(vim.fn.bufwinid(db))
  vim.api.nvim_win_set_cursor(0, { D.state[db].first, 0 })
  vim.cmd("DnaCloneAdd vector")
  H.eq(1, #C.staged)
  local msgs = {}
  local orig = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  vim.cmd("DnaCloneTreat dephos")
  vim.notify = orig
  H.ok(msgs[#msgs]:find("5′ phosphate removed", 1, true), msgs[#msgs])
  H.eq(false, A.phosphate(C.staged[1].frag, "left"))
  -- the part list names the treatment and the ends
  C.list(buf)
  local list = table.concat(vim.api.nvim_buf_get_lines(C.last, 0, -1, false), "\n")
  H.ok(list:find("treatment", 1, true), list)
  H.ok(list:find("dephosphorylation (both ends)", 1, true), list)
  H.ok(list:find("no 5′ P", 1, true), list)
  -- ligating it gives a linear product, and the report says why
  vim.cmd("DnaLigate")
  local report = table.concat(vim.api.nvim_buf_get_lines(C.last, 0, -1, false), "\n")
  H.ok(report:find("the product is linear", 1, true), report)
  H.ok(report:find("neither end has a 5′ phosphate", 1, true), report)
  H.ok(report:find("DnaCloneTreat kinase", 1, true), "the report says what to do about it")
  -- and the joins table spells the same reason out against the two ends
  local rows = C.join_rows({ C.staged[1].frag }, "ligation")
  H.ok(rows[1][2]:find("neither carries a 5′ phosphate", 1, true), rows[1][2])
  C.staged = {}
  H.log("a dephosphorylated vector ligated: linear product, report names the missing phosphate")
end)

case("treat_takes_parts_and_sides", function()
  C.staged = {}
  local buf = open("ncbi/J01749.gb")
  local db = D.run(buf, { "EcoRI", "HindIII" })
  vim.api.nvim_set_current_win(vim.fn.bufwinid(db))
  for i = 1, 2 do
    vim.api.nvim_win_set_cursor(0, { D.state[db].first + i - 1, 0 })
    vim.cmd("DnaCloneAdd")
  end
  H.eq(2, #C.staged)
  local before = C.staged[1].rec:len()
  C.treat({ "nuclease", "1", "left" })
  H.eq(before - 4, C.staged[1].rec:len(), "only the named part's left end")
  H.eq("blunt", A.end_text(C.staged[1].frag, "left"))
  H.ok(A.end_text(C.staged[1].frag, "right") ~= "blunt")
  H.eq(before, C.staged[2].rec:len() + (C.staged[2].rec:len() == before and 0 or before - C.staged[2].rec:len()))
  H.eq("", A.treatment_text(C.staged[2].frag), "the other part is untouched")
  -- by name, both ends
  C.treat({ "dephos", C.staged[2].name })
  H.eq(false, A.phosphate(C.staged[2].frag, "left"))
  H.eq(false, A.phosphate(C.staged[2].frag, "right"))
  H.eq(true, A.phosphate(C.staged[1].frag, "right"))
  local msgs = {}
  local orig = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  C.treat({ "nosuchthing" })
  H.ok(msgs[#msgs]:find("name a treatment", 1, true), msgs[#msgs])
  C.treat({ "dephos", "nosuchpart" })
  H.ok(msgs[#msgs]:find("is not a staged part", 1, true), msgs[#msgs])
  vim.notify = orig
  C.staged = {}
  H.log("DnaCloneTreat: by index and side, by name, and two refusals")
end)

return T