-- The genome-scale circular map. Two things have to be true of it: what it
-- draws is what the record says -- the bins, the GC content and the GC skew
-- are checked against direct counts over the same windows, not against
-- another part of the view -- and no feature goes missing without the header
-- saying so. The last case is the other half of the bargain: a plasmid is
-- drawn exactly the way it was before any of this existed.
local H = require("helpers")

require("plasmid").setup({})
local io_ = require("plasmid.io")
local config = require("plasmid.config")
local ring = require("plasmid.view.ring")
local circular = require("plasmid.view.circular")
local vf = require("plasmid.view.features")
local Record = require("plasmid.model.record").Record

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function opts(over)
  local o = vim.deepcopy(config.options)
  if over then o = vim.tbl_deep_extend("force", o, over) end
  return o
end

local function read(rel)
  return assert(io_.read(H.fixture(rel)), "could not read " .. rel).records[1]
end

--- MINSTD, so the generated records are the same on every machine and in
--- both sequence layouts.
local function rng(seed)
  local s = seed
  return function()
    s = (16807 * s) % 2147483647
    return s / 2147483647
  end
end

--- A sequence of `n` bases whose composition is whatever `weights(i)` says:
--- { A, C, G, T } probabilities, which need not be normalised.
local function make_seq(n, weights)
  local rand = rng(20251124)
  local bases = { "a", "c", "g", "t" }
  local out, buf = {}, {}
  for i = 1, n do
    local w = weights(i)
    local total = w[1] + w[2] + w[3] + w[4]
    local r = rand() * total
    local pick = 4
    for k = 1, 4 do
      if r < w[k] then pick = k break end
      r = r - w[k]
    end
    buf[#buf + 1] = bases[pick]
    if #buf == 8192 then out[#out + 1] = table.concat(buf); buf = {} end
  end
  out[#out + 1] = table.concat(buf)
  return table.concat(out)
end

-- A synthetic chromosome built to turn where we say it turns. The forward
-- strand is G-rich from the origin round to the terminus and C-rich on the
-- way back, which is what a real replicating chromosome does, so the
-- cumulative skew has its maximum at TER and its minimum at ORI. One stretch
-- inside the C-rich half is GC-rich without changing the skew's sign, so the
-- two composition rings are not measuring the same thing.
local SYN_LEN, SYN_TER, SYN_ORI = 600000, 150000, 450000
local synthetic
local function synth()
  if synthetic then return synthetic end
  local seq = make_seq(SYN_LEN, function(i)
    local g_rich = i <= SYN_TER or i > SYN_ORI
    if i > 200000 and i <= 250000 then return { 0.125, 0.45, 0.30, 0.125 } end
    if g_rich then return { 0.25, 0.20, 0.30, 0.25 } end
    return { 0.25, 0.30, 0.20, 0.25 }
  end)
  local features = {}
  for i = 1, 400 do
    local start = (i - 1) * 1500 + 100
    features[#features + 1] = {
      type = "CDS",
      location_text = i % 2 == 0 and ("complement(%d..%d)"):format(start, start + 900)
        or ("%d..%d"):format(start, start + 900),
      qualifiers = { { key = "gene", value = ("g%03d"):format(i) } },
    }
  end
  synthetic = Record.new({
    name = "SYNTH", seq = seq, topology = "circular", molecule_type = "DNA",
    features = features,
  })
  return synthetic
end

--- GC content and GC skew of a stretch, counted a byte at a time: a second
--- implementation, so the view's own counting is checked against something
--- rather than against itself.
local function direct_composition(seq, a, z)
  local g, c, at = 0, 0, 0
  for i = a, z do
    local b = seq:byte(i)
    if b == 103 or b == 71 then g = g + 1
    elseif b == 99 or b == 67 then c = c + 1
    elseif b == 97 or b == 65 or b == 116 or b == 84 or b == 117 or b == 85 then at = at + 1 end
  end
  local valid = g + c + at
  return valid > 0 and (g + c) / valid or nil,
    (g + c) > 0 and (g - c) / (g + c) or nil, g, c
end

--- djb2 over the bytes, kept to 32 bits. The multiplier is small enough that
--- every intermediate product is exact in a double, so the number below is
--- the same number everywhere.
local function digest(s)
  local h = 5381
  for i = 1, #s do h = (h * 33 + s:byte(i)) % 4294967296 end
  return h
end

---------------------------------------------------------------------------

case("tick_intervals_are_round_numbers_sized_for_the_record", function()
  -- plasmid scale keeps the intervals it has always had; genome scale asks
  -- for more of them, because at that size the ring carries coordinates
  -- rather than names
  H.eq({ 1000, 2000, 50000, 1000000 }, {
    ring.tick_step(4361, "arcs"),      -- pBR322: every kb
    ring.tick_step(9609, "arcs"),      -- pPCP1
    ring.tick_step(230218, "arcs"),    -- yeast chromosome I
    ring.tick_step(4641652, "arcs"),   -- E. coli, were it drawn as arcs
  })
  H.eq({ 100000, 200000, 500000, 1000000 }, {
    ring.tick_step(600000, "atlas"),
    ring.tick_step(1531933, "atlas"),  -- a eukaryotic chromosome
    ring.tick_step(4641652, "atlas"),  -- E. coli: every 500 kb
    ring.tick_step(8667507, "atlas"),  -- S. coelicolor
  })
  -- every interval is a round number, and between six and twelve of them
  -- cover the record
  for _, L in ipairs({ 1000, 9609, 230218, 600000, 1531933, 4641652, 8667507 }) do
    for _, mode in ipairs({ "arcs", "atlas" }) do
      local step = ring.tick_step(L, mode)
      local lead = step / 10 ^ math.floor(math.log10(step))
      H.ok(({ [1] = true, [2] = true, [2.5] = true, [5] = true })[lead],
        ("%d is not a round interval for %d (%s)"):format(step, L, mode))
      local n = math.ceil(L / step)
      H.ok(n >= 4 and n <= 12, ("%d ticks for %d bp (%s)"):format(n, L, mode))
    end
  end
  H.log("tick step", "pBR322 1 kb", "E. coli 500 kb")
end)

case("bins_are_round_and_map_to_angles_both_ways", function()
  local L, cells = 4641652, 216
  local bin = ring.bin_bp(L, cells)
  H.eq(25000, bin)
  H.eq(186, ring.bin_count(L, bin))
  H.ok(ring.bin_count(L, bin) <= cells, "more bins than cells to draw them in")
  -- bin_of and bin_span are inverses on every bin, and the spans tile the
  -- record exactly once with no gap and no overlap
  local covered = 0
  for b = 1, ring.bin_count(L, bin) do
    local a, z = ring.bin_span(b, bin, L)
    H.eq(b, ring.bin_of(a, bin))
    H.eq(b, ring.bin_of(z, bin))
    H.eq(a, covered + 1)
    covered = z
  end
  H.eq(L, covered)
  -- and the way the drawing goes from an angle back to a bin agrees with it
  local worst = 0
  for p = 1, L, 977 do
    local b = ring.bin_at_angle(ring.angle_of(p, L), L, bin)
    worst = math.max(worst, math.abs(b - ring.bin_of(p, bin)))
  end
  H.eq(0, worst)
  H.log("bins", ("%d of %d bp"):format(ring.bin_count(L, bin), bin), "angle round trip exact")
end)

case("gc_and_skew_match_a_direct_count_over_the_same_window", function()
  local rec = synth()
  local bin = 25000
  local comp = ring.composition(rec, bin)
  H.eq(ring.bin_count(SYN_LEN, bin), comp.n)
  local seq = rec.seq:tostring()
  local worst_gc, worst_skew, sum_g, sum_c = 0, 0, 0, 0
  for b = 1, comp.n do
    local a, z = ring.bin_span(b, bin, SYN_LEN)
    local gc, skew, g, c = direct_composition(seq, a, z)
    worst_gc = math.max(worst_gc, math.abs(comp.gc[b] - gc))
    worst_skew = math.max(worst_skew, math.abs(comp.skew[b] - skew))
    sum_g, sum_c = sum_g + g, sum_c + c
  end
  H.ok(worst_gc < 1e-12, "gc differs from a direct count by " .. worst_gc)
  H.ok(worst_skew < 1e-12, "skew differs from a direct count by " .. worst_skew)
  -- the mean is over the whole record, not the mean of the bins
  local whole_gc = select(1, direct_composition(seq, 1, SYN_LEN))
  H.ok(math.abs(comp.mean - whole_gc) < 1e-12, "mean gc differs from a direct count")
  -- the cumulative skew is the running sum of (G - C)
  H.eq(sum_g - sum_c, comp.cum[comp.n])
  -- the GC-rich stretch is where it was built, and its skew still points the
  -- same way as the half it sits in
  local rich = ring.bin_of(225000, bin)
  H.ok(comp.gc[rich] > 0.7, "the GC-rich stretch is not GC-rich: " .. comp.gc[rich])
  H.ok(comp.skew[rich] < 0, "the GC-rich stretch flipped the skew")
  H.log("composition", ("%d bins"):format(comp.n),
    ("mean gc %.4f"):format(comp.mean), ("max |delta| %.1e"):format(math.max(worst_gc, worst_skew)))
end)

case("the_skew_turns_at_the_origin_and_the_terminus", function()
  local rec = synth()
  local bin = 25000
  local comp = ring.composition(rec, bin)
  H.ok(comp.turning, "a record built with a replication signal claims none")
  H.ok(math.abs(comp.origin_bp - SYN_ORI) <= bin,
    ("origin %d is more than a bin from %d"):format(comp.origin_bp, SYN_ORI))
  H.ok(math.abs(comp.terminus_bp - SYN_TER) <= bin,
    ("terminus %d is more than a bin from %d"):format(comp.terminus_bp, SYN_TER))
  -- the sign of the skew changes across both of them and nowhere else
  local flips = {}
  for b = 2, comp.n do
    if (comp.skew[b] > 0) ~= (comp.skew[b - 1] > 0) then
      flips[#flips + 1] = select(1, ring.bin_span(b, bin, SYN_LEN))
    end
  end
  H.eq(2, #flips)
  H.ok(math.abs(flips[1] - SYN_TER) <= bin and math.abs(flips[2] - SYN_ORI) <= bin,
    "the skew changes sign somewhere else: " .. vim.inspect(flips))
  H.log("skew turns", ("ori %d (built %d)"):format(comp.origin_bp, SYN_ORI),
    ("ter %d (built %d)"):format(comp.terminus_bp, SYN_TER),
    ("swing %.2f%% of L"):format(100 * comp.swing / SYN_LEN))
end)

case("a_record_with_no_replication_signal_claims_no_origin", function()
  -- the same size and the same GC content, but the strands balanced: the
  -- cumulative skew is a random walk and the view must not name its wobbles
  local seq = make_seq(SYN_LEN, function() return { 0.25, 0.25, 0.25, 0.25 } end)
  local rec = Record.new({ name = "FLAT", seq = seq, topology = "circular" })
  local comp = ring.composition(rec, 25000)
  H.ok(not comp.turning, "a balanced record claims an origin")
  local lay = circular.build(rec, opts(), 120, 40)
  H.ok(not table.concat(lay.lines, "\n"):find("origin (skew)", 1, true),
    "a balanced record names an origin on the map")
  H.log("no signal", ("swing %.3f%% of L"):format(100 * comp.swing / SYN_LEN),
    "turning " .. tostring(comp.turning))
end)

case("density_counts_every_feature_in_every_bin_it_reaches", function()
  local rec = synth()
  local bin, L = 25000, SYN_LEN
  local items = vf.visible(rec, opts())
  local dens = ring.density(items, L, bin)
  -- recounted from the record: a feature is in a bin when its extent and the
  -- bin's span overlap at all
  local fwd, rev = {}, {}
  for b = 1, dens.n do
    fwd[b], rev[b] = 0, 0
    local a, z = ring.bin_span(b, bin, L)
    for _, it in ipairs(items) do
      local hit = false
      for _, e in ipairs(it.extents) do
        if e[1] <= z and e[2] >= a then hit = true end
      end
      if hit then
        if it.f.strand == -1 then rev[b] = rev[b] + 1 else fwd[b] = fwd[b] + 1 end
      end
    end
  end
  H.eq(fwd, dens.fwd)
  H.eq(rev, dens.rev)
  H.eq(0, dens.unstranded)
  -- every feature reaches at least one bin, so the rings hold all of them
  local total = 0
  for b = 1, dens.n do total = total + dens.fwd[b] + dens.rev[b] end
  H.ok(total >= #items, "some feature reaches no bin at all")
  H.log("density", ("%d features"):format(#items),
    ("%d bins"):format(dens.n), ("busiest bin %d"):format(dens.max))
end)

case("landmarks_are_the_ones_worth_naming", function()
  -- seven rRNA operons of three genes each, a replication origin, four large
  -- regions and two hundred ordinary genes. What a reader wants named is the
  -- origin, the operons (one label apiece, not three) and the large regions.
  local features = { { type = "rep_origin", location_text = "500000..500200",
    qualifiers = { { key = "label", value = "oriX" } } } }
  for k = 1, 7 do
    local base = 100000 * k
    for j, part in ipairs({ { "rrs", 1500 }, { "rrl", 2900 }, { "rrf", 120 } }) do
      features[#features + 1] = {
        type = "rRNA", location_text = ("%d..%d"):format(base + j * 4000, base + j * 4000 + part[2]),
        qualifiers = { { key = "gene", value = part[1] .. string.char(64 + k) } },
      }
    end
  end
  -- far enough apart not to be read as one region (the merge distance is two
  -- bins, or ten kilobases, whichever is larger)
  for k = 1, 4 do
    features[#features + 1] = {
      type = "misc_feature", location_text = ("%d..%d"):format(150000 * k, 150000 * k + 20000),
      qualifiers = { { key = "note", value = "island " .. k } },
    }
  end
  for k = 1, 200 do
    features[#features + 1] = {
      type = "CDS", location_text = ("%d..%d"):format(1000 + 3000 * k, 1600 + 3000 * k),
      qualifiers = { { key = "gene", value = ("y%03d"):format(k) } },
    }
  end
  local rec = Record.new({ name = "LAND", seq = string.rep("acgt", 200000),
    topology = "circular", features = features })
  local items = vf.visible(rec, opts())
  local marks = ring.landmarks(items, rec:len(), 25000, 12)
  local names = {}
  for _, m in ipairs(marks) do names[#names + 1] = m.label end
  H.eq("oriX", names[1])
  -- one label per operon: the longest gene of each, never three
  local operons = 0
  for _, n in ipairs(names) do if n:match("^rrl") then operons = operons + 1 end end
  H.eq(7, operons)
  for _, n in ipairs(names) do
    H.ok(not n:match("^rrs") and not n:match("^rrf"), "an operon was named more than once: " .. n)
  end
  -- the four large regions get in; the ordinary genes do not
  local islands = 0
  for _, n in ipairs(names) do if n:match("^island") then islands = islands + 1 end end
  H.eq(4, islands)
  for _, n in ipairs(names) do
    H.ok(not n:match("^y%d"), "an ordinary gene was named at genome scale: " .. n)
  end
  H.eq(12, #marks)
  H.log("landmarks", ("%d of %d features"):format(#marks, #items), table.concat(names, " "))
end)

case("the_header_accounts_for_every_feature", function()
  -- one of each kind the map can lose: an ordinary feature, a feature of a
  -- type the map hides as annotation history, and a feature whose location
  -- is past the end of the sequence
  local features = {}
  for k = 1, 50 do
    features[#features + 1] = { type = "CDS", location_text = ("%d..%d"):format(10000 * k, 10000 * k + 800),
      qualifiers = { { key = "gene", value = ("c%02d"):format(k) } } }
  end
  for k = 1, 3 do
    features[#features + 1] = { type = "variation", location_text = ("%d..%d"):format(5000 * k, 5000 * k + 2) }
  end
  for k = 1, 2 do
    features[#features + 1] = { type = "CDS", location_text = ("%d..%d"):format(700000 + k, 700500 + k) }
  end
  local rec = Record.new({ name = "ACCT", seq = string.rep("acgt", 150000),
    topology = "circular", features = features })
  local lay = circular.build(rec, opts(), 120, 40)
  H.ok(lay.atlas, "a 600 kb record was not drawn as an atlas")
  H.eq(3, lay.history)
  H.eq(2, lay.off_seq)
  H.eq(50, #lay.items)
  -- the title counts every feature the linear view would show, and the
  -- legend says what became of each one
  H.ok(lay.lines[1]:find("55 features", 1, true), "title: " .. lay.lines[1])
  H.eq(#lay.items + lay.history + lay.off_seq, 55)
  for _, want in ipairs({ "50 in the rings", "3 annotation history", "2 past the sequence end" }) do
    H.ok(lay.lines[2]:find(want, 1, true), ("legend has no %q: %s"):format(want, lay.lines[2]))
  end
  -- nothing is over a limit at this scale: the rings hold everything
  H.eq(0, lay.hidden)
  H.log("accounting", (lay.lines[1]:gsub("^%s+", "")))
  H.log("accounting", (lay.lines[2]:gsub("^%s+", "")))
end)

case("the_rings_shrink_and_then_drop_when_the_window_is_small", function()
  local wide = ring.bands(54)
  H.eq(4, #wide)
  H.eq({ "fwd", "rev", "gc", "skew" }, { wide[1].kind, wide[2].kind, wide[3].kind, wide[4].kind })
  local depths = {}
  for i, b in ipairs(wide) do depths[i] = b.depth end
  for i = 2, #wide do H.ok(depths[i] > depths[i - 1], "the rings are not nested outside in") end
  -- as the room runs out the bars get shorter first, and only then does the
  -- innermost ring go
  local counts = {}
  for _, avail in ipairs({ 54, 40, 30, 20, 12, 6, 0 }) do
    local bands, used = ring.bands(avail)
    counts[#counts + 1] = #bands
    H.ok(used <= math.max(avail, 0), ("%d dots of rings in %d"):format(used, avail))
  end
  for i = 2, #counts do H.ok(counts[i] <= counts[i - 1], "a narrower window grew a ring") end
  H.eq(0, counts[#counts])

  -- a ring that was not drawn is not claimed: in a window with room for two
  -- of them the key becomes part of the header that cannot be dropped, and
  -- it names those two and no others
  local rec = Record.new({ name = "SMALL", seq = string.rep("acgt", 150000), topology = "circular" })
  local small = circular.build(rec, opts(), 60, 20)
  H.eq(2, #small.bands)
  H.ok(small.lines[2]:find("rings: forward · reverse", 1, true), "legend: " .. small.lines[2])
  H.ok(not small.lines[2]:find("gc", 1, true), "a ring that is not drawn is in the legend")
  local wide_lay = circular.build(rec, opts(), 120, 40)
  H.eq(4, #wide_lay.bands)
  H.ok(wide_lay.lines[2]:find("gc · skew", 1, true), "legend: " .. wide_lay.lines[2])
  H.log("bands", "rings by room " .. table.concat(counts, " "),
    ("60x20 draws %d, 120x40 draws %d"):format(#small.bands, #wide_lay.bands))
end)

case("a_plasmid_is_drawn_exactly_as_it_was", function()
  -- Nothing on a plasmid may change: the atlas is a second scale, not a new
  -- look. These checksums are of the rendered lines at two window sizes. If
  -- one of them fails, the circular map changed at plasmid scale, which is
  -- the thing this track promised not to do.
  for _, c in ipairs({
    { "ncbi/J01749.gb", 120, 40, 2842888957 },
    { "ncbi/J01749.gb", 80, 24, 4071131924 },
    { "ncbi/L09137.gb", 120, 40, 3645292028 },
    { "ncbi/NC_005816.gb", 100, 40, 1979384210 },
  }) do
    local rec = read(c[1])
    H.ok(not ring.wanted(rec, opts()), c[1] .. " would be drawn as an atlas")
    local lay = circular.build(rec, opts(), c[2], c[3])
    H.ok(not lay.atlas, c[1] .. " was drawn as an atlas")
    H.eq(c[4], digest(table.concat(lay.lines, "\n")),
      ("%s at %dx%d no longer renders as it did"):format(c[1], c[2], c[3]))
    -- and no ring of the atlas leaked into it
    for _, runs in ipairs(lay.runs) do
      for _, r in ipairs(runs) do
        H.ok(not r[3]:match("^PlasmidRing"), "an atlas ring on a plasmid: " .. r[3])
      end
    end
  end
  -- every committed fixture is below the threshold, so the suite's own
  -- corpus is all plasmid scale
  local over = {}
  for _, path in ipairs(H.fixtures("**/*.gb")) do
    local res = io_.read(path)
    for _, rec in ipairs((res or {}).records or {}) do
      if ring.wanted(rec, opts()) then over[#over + 1] = path end
    end
  end
  H.eq({}, over)
  H.log("plasmid scale", "4 renders byte-identical", "0 fixtures over the threshold")
end)

case("the_threshold_is_the_only_thing_that_chooses_the_scale", function()
  local rec = read("ncbi/J01749.gb")
  H.ok(not ring.wanted(rec, opts()))
  H.ok(ring.wanted(rec, opts({ circular = { atlas = { enabled = true } } })))
  H.ok(ring.wanted(rec, opts({ circular = { atlas = { threshold = 1000 } } })))
  local big = Record.new({ name = "BIG", seq = string.rep("acgt", 150000), topology = "circular" })
  H.ok(ring.wanted(big, opts()))
  H.ok(not ring.wanted(big, opts({ circular = { atlas = { enabled = false } } })))
  -- forced on, a plasmid gets the atlas and every one of its features is in
  -- the rings
  local lay = circular.build(rec, opts({ circular = { atlas = { enabled = true } } }), 120, 40)
  H.ok(lay.atlas, "forcing the atlas on did not")
  H.eq(#vf.visible(rec, opts()), #lay.items + lay.history + lay.off_seq)
  H.log("threshold", ("%d bp"):format(config.options.circular.atlas.threshold),
    "enabled=auto|true|false honoured")
end)

return T
