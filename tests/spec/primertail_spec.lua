-- Primers with 5' tails: a Gibson overhang, a restriction site, an overlap.
--
-- A tail does not anneal, so it has no position of its own on the template.
-- The view draws it anyway, in lower case, on the bases the primer hangs over
-- -- which is what a bench scientist wants to see and what makes the tail the
-- part of a primer most likely to be drawn in the wrong place: it is the part
-- that runs off the end of a linear record, round the origin of a circular
-- one, over the line break, and over its neighbours.
--
-- These cases check that every base the view claims to draw is drawn once,
-- with the right glyph, in the right colour, on the right side of the
-- sequence row, and that nothing downstream mistakes a tail for a base that
-- anneals.
local H = require("helpers")

require("plasmid").setup({})
local io_ = require("plasmid.io")
local config = require("plasmid.config")
local linear = require("plasmid.view.linear")
local buffer = require("plasmid.buffer")
local primers = require("plasmid.primers")
local P = require("plasmid.bio.primer")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

-- tails with no match anywhere in pBR322, so the annealed part is exactly the
-- template stretch asked for
local GIBSON = "GGCCCTTTCGTCTTCAAGAA"
local GIBSON2 = "TAAGCTGTCAAACATGAGAA"

local function load(rel, topology)
  local rec = assert(io_.read(H.fixture(rel))).records[1]
  if topology then rec.topology = topology end
  return rec
end

--- Add `tail` + the template at s..e on `strand` as a primer feature.
local function add_primer(rec, s, e, strand, tail, name)
  local L, tmpl = rec:len(), rec.seq:tostring()
  local core = e >= s and tmpl:sub(s, e) or (tmpl:sub(s, L) .. tmpl:sub(1, e))
  local seq = tail .. (strand == 1 and core or P.revcomp(core))
  for _, x in ipairs(P.bind(tmpl, rec:is_circular(), seq, {})) do
    local three = x.strand == 1 and x.stop or x.start
    if x.strand == strand and three == (strand == 1 and e or s) then
      local f = rec:add_feature(primers.make_feature(rec, x, seq, name))
      rec.revision = (rec.revision or 0) + 1
      return f, x
    end
  end
  error(("no binding site for %s at %d..%d on strand %d"):format(name, s, e, strand))
end

--- The character at display column `vcol` (0-based) of `text`, with its bytes.
local function char_at_col(text, vcol)
  local col, b = 0, 1
  while b <= #text do
    local ch = text:match("^[%z\1-\127\194-\244][\128-\191]*", b)
    if not ch then return nil end
    if col == vcol then return ch, b - 1 end
    col = col + vim.fn.strdisplaywidth(ch)
    b = b + #ch
  end
end

local GROUP = { tail = "PlasmidPrimerTail", anneal = "PlasmidPrimer", mismatch = "PlasmidPrimerMismatch" }

--- Build the view and check every cell of every primer against what is drawn.
local function check_drawing(rec, G, width, tag)
  local o = vim.deepcopy(config.options)
  o.linear.group = G
  o.linear.window = { enabled = false }
  local list = primers.list(rec)
  H.ok(#list > 0, tag .. ": no primers to draw")
  local lay = linear.build(rec, o, width, {
    primers = list,
    primer_cells = function(p) return primers.cells(rec, p) end,
  })
  local by_block = {}
  for row0, info in pairs(lay.rows) do
    if info.kind == "primer" then
      by_block[info.b] = by_block[info.b] or {}
      table.insert(by_block[info.b], row0)
    end
  end
  local drawn = 0
  for _, p in ipairs(list) do
    for _, c in ipairs(primers.cells(rec, p)) do
      local b = math.floor((c.pos - 1) / lay.W)
      local vcol = lay.prefix_w + lay:col0(c.pos - (b * lay.W + 1))
      local hits = 0
      for _, row0 in ipairs(by_block[b] or {}) do
        local info = lay.rows[row0]
        local owner
        for _, sg in ipairs(info.segs or {}) do
          if vcol >= sg[1] and vcol <= sg[2] then owner = sg[3] end
        end
        if owner == p.f.id and info.strand == p.strand then
          local ch, b0 = char_at_col(lay.lines[row0 + 1], vcol)
          if ch == c.ch then
            hits = hits + 1
            local group
            for _, r in ipairs(info.runs or {}) do
              if b0 >= r[1] and b0 < r[2] then group = r[3] end
            end
            H.eq(GROUP[c.kind], group, ("%s %s: base %d (%s) colour"):format(tag, p.name, c.pos, c.kind))
            local seqrow = lay.block_row[b]
            if p.strand == 1 then
              H.ok(row0 < seqrow, ("%s %s: forward primer drawn below the sequence"):format(tag, p.name))
            else
              H.ok(row0 > seqrow, ("%s %s: reverse primer drawn above the sequence"):format(tag, p.name))
            end
          end
        end
      end
      H.eq(1, hits, ("%s %s: base %d (%s, %s) drawn %d times"):format(tag, p.name, c.pos, c.ch, c.kind, hits))
      drawn = drawn + 1
    end
  end
  return drawn, lay
end

case("a_gibson_tail_is_drawn_where_the_view_says", function()
  local total = 0
  for _, topology in ipairs({ "circular", "linear" }) do
    local rec = load("ncbi/J01749.gb", topology)
    add_primer(rec, 1000, 1024, 1, GIBSON, "FWD")
    add_primer(rec, 1400, 1424, -1, GIBSON2, "REV")
    for _, G in ipairs({ 10, 0 }) do
      for _, width in ipairs({ 76, 100, 140 }) do
        local n = check_drawing(rec, G, width, ("%s G=%d w=%d"):format(topology, G, width))
        total = total + n
      end
    end
  end
  H.log(("    [tails] %d primer bases checked against the rows that draw them"):format(total))
end)

-- Where the tail falls off the record, or round it.
case("a_tail_past_the_end_of_a_linear_record_is_clipped", function()
  local rec = load("ncbi/J01749.gb", "linear")
  add_primer(rec, 5, 29, 1, GIBSON, "NEAR0")
  local p = primers.list(rec)[1]
  H.eq(20, #p.tail)
  local cells = primers.cells(rec, p)
  -- 4 of the 20 tail bases fall on positions 1..4; the other 16 are off the end
  H.eq(29, #cells, "cells drawn for a primer whose tail runs off the start")
  for _, c in ipairs(cells) do
    H.ok(c.pos >= 1 and c.pos <= rec:len(), "a cell outside the record at " .. c.pos)
  end
  check_drawing(rec, 10, 100, "linear tail off the start")

  local rec2 = load("ncbi/J01749.gb", "linear")
  local L = rec2:len()
  add_primer(rec2, L - 28, L - 4, -1, GIBSON, "ENDREV")
  local q = primers.list(rec2)[1]
  H.eq(29, #primers.cells(rec2, q), "cells drawn for a reverse tail running off the end")
  check_drawing(rec2, 10, 100, "linear tail off the end")
end)

case("a_tail_wraps_the_origin_of_a_circular_record", function()
  local rec = load("ncbi/J01749.gb")
  local L = rec:len()
  add_primer(rec, 5, 29, 1, GIBSON, "ORI")
  local p = primers.list(rec)[1]
  local cells = primers.cells(rec, p)
  H.eq(#p.seq, #cells, "a circular record should draw every base of the primer")
  local wrapped = 0
  for _, c in ipairs(cells) do if c.pos > L - 100 then wrapped = wrapped + 1 end end
  H.eq(16, wrapped, "tail bases drawn at the far end of the record")
  check_drawing(rec, 10, 100, "tail across the origin")

  -- and the primer whose annealed part crosses the origin
  local rec2 = load("ncbi/J01749.gb")
  add_primer(rec2, rec2:len() - 10, 14, 1, GIBSON, "WRAP")
  local q = primers.list(rec2)[1]
  H.ok(q.start > q.stop, "the annealed region should wrap")
  H.eq(#q.seq, #primers.cells(rec2, q), "every base of an origin-spanning primer is drawn")
  check_drawing(rec2, 10, 100, "primer across the origin")
end)

case("overlapping_tails_stay_distinct", function()
  local rec = load("ncbi/J01749.gb")
  add_primer(rec, 1000, 1024, 1, GIBSON, "A1")
  add_primer(rec, 1010, 1034, 1, GIBSON2, "A2")
  add_primer(rec, 1020, 1044, 1, "TTTGGGACACCCTTTGGGAC", "A3")
  H.eq(3, #primers.list(rec))
  check_drawing(rec, 10, 100, "three overlapping tails")
end)

-- The tail is not part of what anneals, and nothing downstream may treat it as
-- if it were.
case("the_amplicon_leaves_the_tails_out_of_its_annealing_temperature", function()
  local rec = load("ncbi/J01749.gb")
  add_primer(rec, 1000, 1024, 1, GIBSON, "GF")
  add_primer(rec, 1400, 1424, -1, GIBSON2, "GR")
  local list = primers.list(rec)
  local prod, ta = require("plasmid.mouse").amplicon(rec, list[1], list[2])
  H.ok(prod, "a tailed pair makes no product")
  local tm = require("plasmid.bio.tm")
  local cond = require("plasmid.primerinfo").conditions()
  local annealed = select(1, tm.ta(list[1].seq:sub(#list[1].tail + 1), list[2].seq:sub(#list[2].tail + 1), cond))
  local whole = select(1, tm.ta(list[1].seq, list[2].seq, cond))
  H.eq(annealed, ta, "Ta should be the annealed parts' Ta")
  H.ok(whole > ta, ("the whole primers' Ta (%d) should be higher than the annealed one (%d); "
    .. "if they are equal this case no longer tests anything"):format(whole, ta))
  H.log(("    [tails] Ta from the annealed parts %d °C, from the whole primers %d °C"):format(ta, whole))
end)

case("the_product_length_counts_the_tails_and_the_span_does_not", function()
  local rec = load("ncbi/J01749.gb")
  add_primer(rec, 1000, 1024, 1, GIBSON, "GF")
  add_primer(rec, 1400, 1424, -1, GIBSON2, "GR")
  local list = primers.list(rec)
  local prod = require("plasmid.mouse").amplicon(rec, list[1], list[2])
  H.eq(425, prod.region, "template span")
  H.eq(465, prod.length, "product length with both tails")
  H.eq(prod.region, prod.stop - prod.start + 1, "the region is what start..stop covers")
  -- what the product actually is
  local seq = require("plasmid.pcr").product_sequence(rec.seq:tostring(), prod)
  H.eq(prod.length, #seq, "the product sequence is as long as the report says")
  H.eq(GIBSON, seq:sub(1, #GIBSON):upper(), "the product should start with the forward tail")
  H.eq(P.revcomp(GIBSON2), seq:sub(-#GIBSON2):upper(), "the product should end with the reverse tail")
end)

case("a_tail_survives_a_save_and_a_reload", function()
  local tmp = vim.fn.tempname()
  vim.fn.mkdir(tmp, "p")
  local path = tmp .. "/tail.gb"
  local fh = assert(io.open(path, "wb")); fh:write(H.read(H.fixture("ncbi/J01749.gb"))); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(path))
  local buf = vim.api.nvim_get_current_buf()
  local rec = buffer.record(buffer.get(buf))
  primers.add_sequence(buf, GIBSON .. rec.seq:tostring():sub(1000, 1024), "RT")
  local before = primers.list(buffer.record(buffer.get(buf)))[1]
  vim.cmd("silent write")
  vim.cmd("silent bwipeout!")
  vim.cmd("silent edit " .. vim.fn.fnameescape(path))
  local after = primers.list(buffer.record(buffer.get(vim.api.nvim_get_current_buf())))[1]
  H.ok(after, "the primer is gone after a save and reload")
  H.eq(before.seq, after.seq, "primer sequence")
  H.eq(before.tail, after.tail, "tail")
  H.eq(before.anneal, after.anneal, "annealed length")
  H.eq({ before.start, before.stop }, { after.start, after.stop }, "binding site")
  vim.cmd("silent bwipeout!")
end)

-- Real primers from SnapGene files, which is where most tailed primers come
-- from in practice.
case("snapgene_primers_with_tails_draw_every_base", function()
  local n, files = 0, 0
  for _, path in ipairs(H.fixtures("snapgene/*.dna")) do
    local res = io_.read(path)
    for _, rec in ipairs(res and res.records or {}) do
      local tailed = {}
      for _, p in ipairs(primers.list(rec)) do if #p.tail > 0 then tailed[#tailed + 1] = p end end
      if #tailed > 0 then
        files = files + 1
        for _, p in ipairs(tailed) do
          local cells = primers.cells(rec, p)
          if rec:is_circular() then
            H.eq(#p.seq, #cells, ("%s %s: cells drawn"):format(vim.fn.fnamemodify(path, ":t"), p.name))
          end
          local tails = 0
          for _, c in ipairs(cells) do if c.kind == "tail" then tails = tails + 1 end end
          H.ok(tails <= #p.tail, ("%s %s: %d tail cells for a %d nt tail"):format(
            vim.fn.fnamemodify(path, ":t"), p.name, tails, #p.tail))
          n = n + 1
        end
        check_drawing(rec, 10, 100, vim.fn.fnamemodify(path, ":t"))
      end
    end
  end
  H.ok(n >= 8, "expected tailed primers in the SnapGene fixtures, found " .. n)
  H.log(("    [tails] %d tailed primers from %d SnapGene records drawn in full"):format(n, files))
end)

-- A primer that binds in more than one place. The view draws every site it
-- was given, so when one of them is weak the drawing is right and the choice
-- is what matters: which site a reader is offered first, and which one the
-- message describes.
case("a_primer_with_two_sites_offers_the_stronger_one_first", function()
  local rec = load("ncbi/J01749.gb")
  local tmpl = rec.seq:tostring():upper()
  -- pBR322 carries AGCAAAAGGCCA twice, 11 bp apart; a primer ending there
  -- anneals 12 nt at both
  local weak = P.bind(tmpl, true, "GGCCCTTTCGTCTTCAAGAAAGCAAAAGGCCA", {})
  H.eq(2, #weak, "the repeated anchor should give two sites")
  for _, s in ipairs(weak) do H.eq(12, s.anneal) end

  -- a primer whose strong site is later in the record than its weak one
  local strong_at = 3000
  local core = tmpl:sub(strong_at, strong_at + 24)
  local sites = P.bind(tmpl, true, core, {})
  H.eq(1, #sites, "the 25-mer should bind once")
  -- plant the weak site: the same 3' anchor earlier in the record, with the
  -- rest diverging, by searching for the anchor elsewhere
  local anchor = core:sub(-10)
  local first = tmpl:find(anchor, 1, true)
  H.ok(first, "the anchor should occur somewhere")

  -- ranking: strongest first, whatever the order in the record
  local mixed = {
    { start = 100, stop = 111, strand = 1, anneal = 12, exact3 = 12, tail = "", mismatches = {} },
    { start = 3000, stop = 3024, strand = 1, anneal = 25, exact3 = 25, tail = "", mismatches = {} },
    { start = 500, stop = 524, strand = 1, anneal = 25, exact3 = 18, tail = "", mismatches = { 8 } },
  }
  local ranked = primers.rank_sites(mixed)
  H.eq(3000, ranked[1].start, "the fully exact site should come first")
  H.eq(500, ranked[2].start, "then the longer site with a mismatch")
  H.eq(100, ranked[3].start, "then the short one")
  -- ranking must not disturb what it was given
  H.eq(100, mixed[1].start, "rank_sites edited its argument")
end)

case("every_site_of_a_multi_site_primer_is_drawn", function()
  local rec = load("ncbi/J01749.gb")
  local tmpl = rec.seq:tostring():upper()
  local seq = "GGCCCTTTCGTCTTCAAGAAAGCAAAAGGCCA"
  local sites = P.bind(tmpl, true, seq, {})
  H.eq(2, #sites)
  for i, s in ipairs(sites) do
    rec:add_feature(primers.make_feature(rec, s, seq, "REP_" .. i))
  end
  rec.revision = (rec.revision or 0) + 1
  local list = primers.list(rec)
  H.eq(2, #list, "both sites should be listed")
  -- 11 bp apart, so they overlap on the same line and need two lanes
  H.ok(math.abs(list[1].start - list[2].start) < 20, "the two sites should be close together")
  for _, G in ipairs({ 10, 0 }) do
    check_drawing(rec, G, 100, ("two sites G=%d"):format(G))
  end
  H.log("    [tails] a primer binding twice 11 bp apart: both sites drawn, in both layouts")
end)

return T
