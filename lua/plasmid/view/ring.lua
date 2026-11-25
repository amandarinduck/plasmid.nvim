-- What the circular map shows once a record is too large to draw one arc per
-- feature: a genome atlas.
--
-- On a 120x40 terminal the ring is about 870 braille dots around. A plasmid
-- puts a few bases behind each dot, so every feature is still its own arc and
-- the map reads the way a plasmid map should. A 4.6 Mb genome puts five
-- kilobases -- several genes -- behind each dot, and E. coli K-12 codes on
-- both strands almost the whole way round: drawn as arcs that is a solid
-- band, and the track packer has to throw three quarters of the features away
-- to make even that. So above `circular.atlas.threshold` residues the map
-- draws what is still true at that scale instead. How many features fall in
-- each angular bin on each strand, and the two base-composition curves a
-- genome atlas is actually read by: GC content against the record's mean, and
-- GC skew, whose sign changes sit at the replication origin and terminus.
--
-- Everything here is computed from the record, never from the drawing: the
-- bins are round numbers of bases, and the same bin is the window GC and skew
-- are measured over, so the header can state one number for both.
local M = {}

local TAU = 2 * math.pi

-- Rings from the backbone inwards. Feature density first, because it is what
-- the view used to show; composition inside it, the way genome atlases are
-- drawn.
local RING_ORDER = { "fwd", "rev", "gc", "skew" }
M.RING_COUNT = #RING_ORDER

-- A cumulative-skew swing smaller than this fraction of the record is noise
-- rather than a replication signal, and the origin and terminus are not
-- claimed. The floor is set from what a real chromosome and a sequence with
-- no replication signal measure; docs/ring_report.md has both numbers.
local MIN_SWING = 0.004

-- What is worth naming at genome scale, before length is taken into account.
-- A replication origin is the one landmark every reader looks for; the rRNA
-- operons are the next; everything else earns its label by being large.
local LANDMARK_TYPE = {
  rep_origin = 1000, oriT = 1000, rRNA = 600, ncRNA = 150, mobile_element = 120,
}
local LANDMARK_LENGTH = 300

---------------------------------------------------------------------------
-- Round numbers
---------------------------------------------------------------------------

--- The smallest round number -- 1, 2, 2.5 or 5 times a power of ten -- that
--- `target` of them cover L. Tick intervals and bin widths are both read off
--- a ring by eye, so both have to be numbers a reader can divide.
function M.nice_step(L, target)
  local raw = L / target
  local mag = 10 ^ math.floor(math.log10(math.max(raw, 1)))
  for _, m in ipairs({ 1, 2, 2.5, 5, 10 }) do
    if raw <= m * mag then return math.max(1, math.floor(m * mag)) end
  end
  return math.floor(10 * mag)
end

--- Tick interval. The atlas asks for more ticks than the plasmid map because
--- its ring carries coordinates rather than names: ten ticks put a bacterial
--- genome on half-megabase marks, eight put a plasmid on kilobase marks.
function M.tick_step(L, mode)
  return M.nice_step(L, mode == "atlas" and 10 or 8)
end

function M.fmt_bp(n)
  if n >= 1e6 then return ((("%.2f"):format(n / 1e6):gsub("%.?0+$", ""))) .. "M" end
  if n >= 1e4 then return ((("%.1f"):format(n / 1e3):gsub("%.0$", ""))) .. "k" end
  return tostring(n)
end

-- Kilo and mega forms of the record's unit: bases contract to kb and Mb the
-- way everyone writes them, the rest just take a prefix.
local SCALED = { bp = { "kb", "Mb" }, nt = { "knt", "Mnt" }, aa = { "k aa", "M aa" } }

--- A length in prose, with its unit: "10 kb", "500 kb", "1.5 Mb". Ticks are
--- read against the ring and can be terse; the header is read as a sentence.
function M.fmt_span(n, unit)
  unit = unit or "bp"
  local scaled = SCALED[unit] or { "k" .. unit, "M" .. unit }
  if n >= 1e6 then return ((("%.2f"):format(n / 1e6):gsub("%.?0+$", ""))) .. " " .. scaled[2] end
  if n >= 1e3 then return ((("%.1f"):format(n / 1e3):gsub("%.0$", ""))) .. " " .. scaled[1] end
  return ("%d %s"):format(n, unit)
end

---------------------------------------------------------------------------
-- Bins
---------------------------------------------------------------------------

--- Bin width in bases: the finest round number that fits `max_bins` of them
--- around the ring. The caller passes the number of character cells the ring
--- is around, so a bin is about one cell wide. A bin per braille dot was
--- tried first and is worse: at a dot apiece the rings carry two bins per
--- cell, and a ring of 10 kb counts that differ from neighbour to neighbour
--- draws as a band of scribble rather than a profile.
function M.bin_bp(L, max_bins)
  return M.nice_step(L, math.max(1, max_bins))
end

function M.bin_count(L, bin)
  return math.max(1, math.ceil(L / bin))
end

--- The bin a base pair falls in (1-based).
function M.bin_of(pos, bin)
  return math.floor((pos - 1) / bin) + 1
end

--- First and last base pair of a bin; the last bin is short when the record
--- is not a whole number of bins long.
function M.bin_span(b, bin, L)
  return (b - 1) * bin + 1, math.min(b * bin, L)
end

--- The angle a base pair sits at, matching the circular view's own geometry:
--- base 1 at 12 o'clock, increasing clockwise.
function M.angle_of(pos, L)
  return TAU * (pos - 1) / L
end

--- The bin an angle points into: the inverse of `angle_of` followed by
--- `bin_of`, and the only place the drawing turns an angle into data.
function M.bin_at_angle(th, L, bin)
  return M.bin_of(math.floor(th % TAU / TAU * L) + 1, bin)
end

---------------------------------------------------------------------------
-- Base composition
---------------------------------------------------------------------------

--- GC content and GC skew per bin, plus the cumulative skew and where it
--- turns. Ambiguity codes count towards neither: `gc` is (G+C)/(A+C+G+T) and
--- `skew` is (G-C)/(G+C), both over the bin, and a bin with no unambiguous
--- bases is `false` rather than a made-up number.
---
--- Cached on the record (one array of numbers per bin), because a resize
--- rebuilds the layout and the sequence has not changed.
function M.composition(rec, bin)
  local cached = rec._ring_composition
  if cached and cached.bin == bin and cached.revision == rec.revision then return cached end
  local L, s = rec:len(), rec.seq
  local nb = M.bin_count(L, bin)
  local gc, skew, cum = {}, {}, {}
  local tot_gc, tot_valid, acc = 0, 0, 0
  for b = 1, nb do
    local a, z = M.bin_span(b, bin, L)
    local sub = s:sub(a, z)
    local _, ng = sub:gsub("[gG]", "")
    local _, nc = sub:gsub("[cC]", "")
    local _, nat = sub:gsub("[atuATU]", "")
    local valid = ng + nc + nat
    gc[b] = valid > 0 and (ng + nc) / valid or false
    skew[b] = (ng + nc) > 0 and (ng - nc) / (ng + nc) or false
    acc = acc + ng - nc
    cum[b] = acc
    tot_gc = tot_gc + ng + nc
    tot_valid = tot_valid + valid
  end
  local mean = tot_valid > 0 and tot_gc / tot_valid or 0
  local maxdev, maxskew = 0, 0
  for b = 1, nb do
    if gc[b] then maxdev = math.max(maxdev, math.abs(gc[b] - mean)) end
    if skew[b] then maxskew = math.max(maxskew, math.abs(skew[b])) end
  end
  -- Where the cumulative (G-C) turns. On a bacterial chromosome the forward
  -- strand is the leading strand from the origin round to the terminus and
  -- the lagging strand on the way back, so the sum of (G-C) falls to its
  -- minimum at the origin and rises to its maximum at the terminus. The
  -- turning point is the bin boundary the sum is lowest (highest) at.
  local bmin, bmax = 1, 1
  for b = 1, nb do
    if cum[b] < cum[bmin] then bmin = b end
    if cum[b] > cum[bmax] then bmax = b end
  end
  local comp = {
    bin = bin, revision = rec.revision, n = nb, gc = gc, skew = skew, cum = cum,
    mean = mean, maxdev = maxdev, maxskew = maxskew,
    swing = cum[bmax] - cum[bmin],
    origin_bp = select(2, M.bin_span(bmin, bin, L)),
    terminus_bp = select(2, M.bin_span(bmax, bin, L)),
  }
  -- A linear record has no replication origin to find this way, and a protein
  -- has no bases at all.
  comp.turning = rec:is_circular() and not rec:is_protein()
    and comp.swing >= MIN_SWING * math.max(L, 1)
  rec._ring_composition = comp
  return comp
end

---------------------------------------------------------------------------
-- Feature density
---------------------------------------------------------------------------

--- How many features cover each bin, per strand. A feature counts in every
--- bin it reaches, so the two rings read as occupancy rather than as starts;
--- a feature that crosses the origin counts on both sides of it.
---
--- An unstranded feature is counted with the forward strand, and `unstranded`
--- says how many there were so the header can admit it.
function M.density(items, L, bin)
  local nb = M.bin_count(L, bin)
  local fwd, rev = {}, {}
  for b = 1, nb do fwd[b], rev[b] = 0, 0 end
  local unstranded = 0
  for _, it in ipairs(items) do
    local strand = it.f.strand
    if strand ~= 1 and strand ~= -1 then unstranded = unstranded + 1 end
    local into = strand == -1 and rev or fwd
    for _, e in ipairs(it.extents) do
      local b0 = M.bin_of(math.max(1, math.min(e[1], L)), bin)
      local b1 = M.bin_of(math.max(1, math.min(e[2], L)), bin)
      for b = b0, b1 do into[b] = into[b] + 1 end
    end
  end
  local max = 0
  for b = 1, nb do max = math.max(max, fwd[b], rev[b]) end
  return { fwd = fwd, rev = rev, max = max, n = nb, unstranded = unstranded }
end

---------------------------------------------------------------------------
-- Landmarks
---------------------------------------------------------------------------

local function circular_distance(a, b, L)
  local d = math.abs(a - b)
  return math.min(d, L - d)
end

local function midpoint(it, L)
  -- a feature that crosses the origin is named at the origin rather than at
  -- the far side of the record
  if #it.extents > 1 then return 1 end
  local e = it.extents[1]
  return math.floor((e[1] + math.min(e[2], L)) / 2)
end

--- The features worth naming at genome scale, best first, at most `limit` of
--- them. Score is the type's own worth plus its size relative to the largest
--- feature on the record, so a replication origin and the rRNA operons are
--- named whatever their length, and the rest are named for being large.
---
--- Neighbours of the same type within `merge` bases are one landmark: E. coli
--- has twenty-two rRNA genes but seven rRNA operons, and seven is what a
--- reader wants to see.
function M.landmarks(items, L, bin, limit)
  local longest = 1
  for _, it in ipairs(items) do longest = math.max(longest, it.f:length()) end
  local cand = {}
  for _, it in ipairs(items) do
    if it.label and it.label ~= it.f.type then
      cand[#cand + 1] = {
        pos = midpoint(it, L), label = it.label, type = it.f.type, id = it.f.id, hl = it.hl,
        score = (LANDMARK_TYPE[it.f.type] or 0) + LANDMARK_LENGTH * it.f:length() / longest,
      }
    end
  end
  table.sort(cand, function(a, b)
    if a.score ~= b.score then return a.score > b.score end
    if a.pos ~= b.pos then return a.pos < b.pos end
    return a.id < b.id
  end)
  local merge = math.max(2 * bin, 10000)
  local out = {}
  for _, c in ipairs(cand) do
    if #out >= limit then break end
    local near
    for _, k in ipairs(out) do
      if k.type == c.type and circular_distance(k.pos, c.pos, L) <= merge then near = k break end
    end
    if not near then out[#out + 1] = c end
  end
  return out
end

---------------------------------------------------------------------------
-- Geometry
---------------------------------------------------------------------------

--- Radial bands for the rings, outermost first, in vertical dot units:
--- `{ kind, depth (below the first ring's outer edge), h }`. A density ring
--- is a baseline with bars growing inwards; a signed ring is a baseline with
--- bars on both sides of it. When `avail` is too small the bars are made
--- shorter, and when that is not enough the innermost rings are dropped --
--- the caller says in the header which ones it got.
function M.bands(avail, kinds)
  kinds = kinds or RING_ORDER
  local hd = math.max(2, math.min(6, math.floor(avail * 0.13)))
  local hs = math.max(1, math.min(5, math.floor(avail * 0.10)))
  local gap = 3
  local function height(kind)
    if kind == "gc" or kind == "skew" then return 2 * hs + 1 end
    return hd + 1
  end
  local function total(n)
    local t = 0
    for i = 1, n do t = t + height(kinds[i]) + (i > 1 and gap or 0) end
    return t
  end
  local n = #kinds
  while n > 0 and total(n) > avail do
    if hd > 2 or hs > 1 then
      hd, hs = math.max(2, hd - 1), math.max(1, hs - 1)
    else
      n = n - 1
    end
  end
  local bands, depth = {}, 0
  for i = 1, n do
    bands[i] = { kind = kinds[i], depth = depth, h = height(kinds[i]), hd = hd, hs = hs }
    depth = depth + height(kinds[i]) + gap
  end
  return bands, math.max(0, depth - gap)
end

--- Whether this record gets the atlas instead of one arc per feature.
--- `atlas.enabled = true` or `false` forces it either way; "auto" decides by
--- size alone, so a record does not change its mind when a window is resized.
function M.wanted(rec, opts)
  local a = ((opts or {}).circular or {}).atlas or {}
  if a.enabled == true then return true end
  if a.enabled == false then return false end
  return rec:len() >= (a.threshold or 500000)
end

---------------------------------------------------------------------------
-- Drawing
---------------------------------------------------------------------------

--- Paint the rings. `geom` is the circular view's own geometry:
--- `{ pt(theta, r), R, L, steps }`, where `steps` is the number of angular
--- samples around the ring (two per braille dot, so the bars have no gaps).
--- Returns the bands it drew, outermost first.
function M.draw(cv, geom, data)
  local L, R, pt = geom.L, geom.R, geom.pt
  local bands = data.bands
  local comp, dens, bin = data.comp, data.density, data.bin

  local function radial(th, r0, r1, group)
    local step = r1 >= r0 and 0.5 or -0.5
    local r = r0
    while (step > 0 and r <= r1 + 1e-9) or (step < 0 and r >= r1 - 1e-9) do
      local x, y = pt(th, r)
      cv:dot(x, y, group)
      r = r + step
    end
  end

  -- Each ring is a step profile around its baseline rather than a filled
  -- bar: a braille cell is 2x4 dots, so a bar filled at this width comes out
  -- as a solid block and the ring reads as the blob the atlas is here to
  -- replace. The riser between one bin's level and the next is drawn, so the
  -- profile is a continuous line.
  for _, band in ipairs(bands) do
    local base = R - 3 - band.depth
    local signed = band.kind == "gc" or band.kind == "skew"
    local axis = signed and (base - band.hs) or base
    local prev
    for k = 0, geom.steps do
      local th = TAU * (k % geom.steps) / geom.steps
      local b = M.bin_at_angle(th, L, bin)
      local level, group = 0, nil
      if band.kind == "fwd" or band.kind == "rev" then
        local count = (band.kind == "fwd" and dens.fwd or dens.rev)[b] or 0
        if count > 0 and dens.max > 0 then
          level = -math.max(1, math.floor(band.hd * count / dens.max + 0.5))
        end
        group = band.kind == "fwd" and "PlasmidRingForward" or "PlasmidRingReverse"
      elseif band.kind == "gc" then
        local v = comp.gc[b] and comp.maxdev > 0 and (comp.gc[b] - comp.mean) / comp.maxdev or 0
        level = v == 0 and 0 or (v > 0 and 1 or -1) * math.max(1, math.floor(band.hs * math.abs(v) + 0.5))
        group = level >= 0 and "PlasmidRingGCHigh" or "PlasmidRingGCLow"
      else
        local v = comp.skew[b] and comp.maxskew > 0 and comp.skew[b] / comp.maxskew or 0
        level = v == 0 and 0 or (v > 0 and 1 or -1) * math.max(1, math.floor(band.hs * math.abs(v) + 0.5))
        group = level >= 0 and "PlasmidRingSkewPos" or "PlasmidRingSkewNeg"
      end
      if k < geom.steps then
        local bx, by = pt(th, axis)
        cv:dot(bx, by, "PlasmidRingBase")
      end
      if level ~= 0 or (prev or 0) ~= 0 then
        radial(th, axis + (prev or level), axis + level, group)
      end
      prev = level
    end
  end
  return bands
end

---------------------------------------------------------------------------
-- Header
---------------------------------------------------------------------------

local RING_TEXT = {
  fwd = { "forward", "PlasmidRingForward" },
  rev = { "reverse", "PlasmidRingReverse" },
  gc = { "gc", "PlasmidRingGCHigh" },
  skew = { "skew", "PlasmidRingSkewPos" },
}

--- The atlas's own header row: what the bins are, what the rings show, and
--- what became of every feature. Returns the line and its highlight runs.
--- Parts are dropped from the end of `optional` until the line fits, so a
--- narrow window loses the colour key before it loses the accounting.
function M.legend(f, width)
  local sep = "  ·  "
  local function seg(text, group) return { text = text, group = group or "PlasmidHeader" } end
  local bins = seg(("%d bins of %s"):format(f.nbins, M.fmt_span(f.bin, f.unit)))
  local window = seg(("gc and skew over %s"):format(M.fmt_span(f.bin, f.unit)))
  local key = { seg("rings: ") }
  for i, band in ipairs(f.rings) do
    if i > 1 then key[#key + 1] = seg(" · ") end
    local t = RING_TEXT[band.kind]
    key[#key + 1] = seg(t[1], t[2])
  end
  local account = { ("%d in the rings"):format(f.summarised) }
  if f.history > 0 then account[#account + 1] = ("%d annotation history"):format(f.history) end
  if f.off_seq > 0 then account[#account + 1] = ("%d past the sequence end"):format(f.off_seq) end
  if f.unstranded > 0 then account[#account + 1] = ("%d unstranded, counted forward"):format(f.unstranded) end
  account[#account + 1] = ("%d named"):format(f.named)
  if f.unnamed > 0 then account[#account + 1] = ("%d with no room to name"):format(f.unnamed) end

  -- In reading order, each with how readily it goes when the line is too
  -- long for the window. The bins and the accounting never go. The colour
  -- key goes first -- unless the window had no room for every ring, when the
  -- key is the only thing that says which rings are on the map and the
  -- window phrase, which names a ring that may not be there, goes instead.
  local reduced = #f.rings < (f.all_rings or #f.rings)
  key.drop = (not reduced) and 1 or nil
  local parts = {
    { bins },
    { window, drop = reduced and 1 or 2 },
    key,
    { seg(table.concat(account, ", ")) },
  }

  local function assemble(limit)
    local segs = { seg(" ") }
    local first = true
    for _, group in ipairs(parts) do
      if not group.drop or group.drop > limit then
        if not first then segs[#segs + 1] = seg(sep) end
        for _, s in ipairs(group) do segs[#segs + 1] = s end
        first = false
      end
    end
    return segs
  end
  local segs
  for limit = 0, 2 do
    segs = assemble(limit)
    local cols = 0
    for _, s in ipairs(segs) do
      -- columns, not bytes: the separators are multibyte
      local _, n_ch = s.text:gsub("[%z\1-\127\194-\244][\128-\191]*", "")
      cols = cols + n_ch
    end
    if cols <= width then break end
  end
  local line, runs, byte = {}, {}, 0
  for _, s in ipairs(segs) do
    if s.text ~= "" then
      line[#line + 1] = s.text
      runs[#runs + 1] = { byte, byte + #s.text, s.group }
      byte = byte + #s.text
    end
  end
  return table.concat(line), runs
end

return M
