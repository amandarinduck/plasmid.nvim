-- Circular map on a braille canvas, at two scales.
--
-- Geometry: angle 0 at 12 o'clock, increasing clockwise; base `bp` starts at
-- theta = 2*pi*(bp-1)/L. Dot space is ~square (2x4 dots per ~1:2 cell); the
-- configurable `aspect` (dot width / dot height) corrects other fonts.
--
-- At plasmid scale (the default below `circular.atlas.threshold` residues)
-- feature arcs are packed into rings inside the backbone, labels are stacked
-- in left/right columns with leader lines, and tick numbers sit outside.
-- Above it an arc is thinner than a braille dot and there are more features
-- than the rings can hold, so the map draws a density atlas instead;
-- `plasmid.view.ring` is that scale, and says why.
local Canvas = require("plasmid.view.canvas")
local vf = require("plasmid.view.features")
local ring = require("plasmid.view.ring")

local M = {}

local TAU = 2 * math.pi
local HEADER_ROWS = 2

local PRIORITY = {
  CDS = 1, gene = 2, promoter = 3, rep_origin = 3, oriT = 3, terminator = 4,
  primer_bind = 5, protein_bind = 5, RBS = 6, mRNA = 2, misc_feature = 7,
}

local nice_step, fmt_bp = ring.nice_step, ring.fmt_bp

local function utf8_len(s)
  local _, n = s:gsub("[%z\1-\127\194-\244][\128-\191]*", "")
  return n
end

local function utf8_trunc(s, n)
  if utf8_len(s) <= n then return s end
  local out, k = {}, 0
  for ch in s:gmatch("[%z\1-\127\194-\244][\128-\191]*") do
    k = k + 1
    if k >= n then break end
    out[#out + 1] = ch
  end
  return table.concat(out) .. "…"
end

---------------------------------------------------------------------------
-- Label columns
---------------------------------------------------------------------------

-- Labels hug the circle just outside the ring, stacked into the free rows of
-- a left and a right column and tied back to their own angle by a leader
-- line. A label that finds no free row near its angle is dropped rather than
-- routed across the map, and the caller counts what it lost. Both scales
-- label this way: at plasmid scale the names are features and restriction
-- sites, at genome scale they are landmarks.
local Labels = {}
Labels.__index = Labels

function Labels.new(cv, geom, label_w)
  return setmetatable({
    cv = cv, cx = geom.cx, cy = geom.cy, ax = geom.ax, pt = geom.pt,
    W = geom.W, H = geom.H, label_r = geom.R + 8, label_w = label_w,
    used = { right = {}, left = {} }, placed = {},
  }, Labels)
end

--- First column, leader-end dot column and dot row of a label on a row of a
--- side: the text starts where the circle leaves that row.
function Labels:placement(side, row, text)
  local y = row * 4 + 2
  local dy = y - self.cy
  local xoff = math.abs(dy) < self.label_r
    and math.sqrt(self.label_r * self.label_r - dy * dy) / self.ax or 0
  local len = utf8_len(text)
  if side == "right" then
    local col = math.ceil((self.cx + xoff) / 2) + 1
    return col, col * 2 - 1, y
  end
  local col = math.floor((self.cx - xoff) / 2) - len
  return col, (col + len) * 2, y
end

--- Place `text` for the angle `th`, searching out from the row that angle
--- points at. `o.from` is the radius its leader starts at (an inner ring; nil
--- starts it at the hinge), `o.hinge` the radius the leader bends at, and
--- `o.tolerance` how many rows away it may settle for. Returns the placed
--- entry, or nil when no row was free in reach.
function Labels:place(th, text, group, owner, o)
  local side = math.sin(th) >= 0 and "right" or "left"
  local _, wy = self.pt(th, self.label_r)
  local want = math.max(0, math.min(self.H - 1, math.floor(wy / 4)))
  for d = 0, o.tolerance do
    for _, row in ipairs(d == 0 and { want } or { want - d, want + d }) do
      if row >= 0 and row < self.H and not self.used[side][row] then
        local avail = side == "right" and (self.W - self:placement(side, row, ""))
          or (self:placement(side, row, "") - 1)
        if avail >= 4 then
          local t = utf8_trunc(text, math.min(self.label_w, avail))
          local col, edge_x, ey = self:placement(side, row, t)
          if self.cv:text_free(col, row, utf8_len(t)) then
            self.used[side][row] = true
            local x1, y1 = self.pt(th, o.hinge)
            if o.from then
              local x0, y0 = self.pt(th, o.from)
              self.cv:line(x0, y0, x1, y1, "PlasmidLeader")
            end
            self.cv:line(x1, y1, edge_x, ey, "PlasmidLeader")
            local entry = { col = col, row = row, text = t, group = group, owner = owner }
            self.placed[#self.placed + 1] = entry
            self.cv:put(col, row, t, group, owner) -- reserve the cells now
            return entry
          end
        end
      end
    end
  end
  return nil
end

--- Write every label again, over whatever was drawn across it since.
function Labels:repaint()
  for _, l in ipairs(self.placed) do self.cv:put(l.col, l.row, l.text, l.group, l.owner) end
end

---------------------------------------------------------------------------
-- Genome scale
---------------------------------------------------------------------------

--- The map above `circular.atlas.threshold`: density and base composition
--- rings instead of one arc per feature. `g` is the geometry and the shared
--- drawing M.build has already set up; plasmid.view.ring is the rest.
local function build_atlas(rec, opts, g)
  local o = opts.circular
  local a = o.atlas or {}
  local cv, ang, L = g.cv, g.ang, g.L
  local geom = g.geom
  local pt, R = geom.pt, geom.R

  g.draw_backbone()
  if o.show_ticks then g.draw_ticks(ring.tick_step(L, "atlas")) end

  -- `steps` samples the ring twice per braille dot and a cell is two dots
  -- wide, so a quarter of them is one bin per character cell
  local bin = ring.bin_bp(L, math.ceil(g.steps / 4))
  local comp = ring.composition(rec, bin)
  local dens = ring.density(g.items, L, bin)
  -- a protein record has no bases to take a composition from
  local kinds = rec:is_protein() and { "fwd", "rev" } or nil
  local bands, depth = ring.bands(math.max(0, R - 15), kinds)

  -- What gets named: the two turning points of the cumulative skew first,
  -- then the landmarks. Every one of them is marked on the ring whether or
  -- not its name found a row, so the reader can see that something is there.
  local marks = {}
  if comp.turning then
    marks[#marks + 1] = { pos = comp.origin_bp, label = "origin (skew)", hl = "PlasmidRingOri" }
    marks[#marks + 1] = { pos = comp.terminus_bp, label = "terminus (skew)", hl = "PlasmidRingOri" }
  end
  for _, lm in ipairs(ring.landmarks(g.items, L, bin, a.landmarks or 12)) do
    marks[#marks + 1] = { pos = lm.pos, label = lm.label, hl = "PlasmidLandmark", id = lm.id }
  end

  local labels, named = nil, 0
  if g.with_labels then
    labels = Labels.new(cv, geom, g.label_w)
    local tolerance = math.max(2, math.floor(geom.H / 8))
    for _, m in ipairs(marks) do
      if labels:place(ang(m.pos), m.label, m.hl, m.id, { tolerance = tolerance, hinge = R + 3 }) then
        named = named + 1
      end
    end
  end
  for _, m in ipairs(marks) do
    local th = ang(m.pos)
    local x0, y0 = pt(th, R - 3)
    local x1, y1 = pt(th, R + 4)
    cv:line(x0, y0, x1, y1, m.hl)
  end

  ring.draw(cv, { pt = pt, R = R, L = L, steps = g.steps },
    { bands = bands, comp = comp, density = dens, bin = bin })

  g.draw_extras()
  if labels then labels:repaint() end

  local inner = R - 4 - depth
  local centre = {}
  if comp.turning then
    centre[1] = ("ori ~%s"):format(ring.fmt_bp(comp.origin_bp))
    centre[2] = ("ter ~%s"):format(ring.fmt_bp(comp.terminus_bp))
  end
  g.centre_text(inner, centre)

  local clines, cruns, cbytes, cowners = cv:render()
  local legend, legend_runs = ring.legend({
    nbins = dens.n, bin = bin, unit = rec:unit(), rings = bands, summarised = #g.items,
    all_rings = kinds and #kinds or ring.RING_COUNT,
    named = named, unnamed = #marks - named, history = g.history, off_seq = g.off_seq,
    unstranded = dens.unstranded,
  }, geom.W)
  local lines = {
    (" %s  ·  %d %s  ·  %s  ·  %s  ·  %s"):format(
      rec.name, rec:len(), rec:unit(), rec.molecule_type or "DNA", rec.topology,
      vf.count_text(g.items, g.history + g.off_seq)),
    legend,
  }
  for _, l in ipairs(clines) do lines[#lines + 1] = l end

  local layout = {
    view = "circular", atlas = true, lines = lines, header_rows = HEADER_ROWS,
    runs = cruns, cell_bytes = cbytes, owners = cowners, legend_runs = legend_runs,
    cx = geom.cx, cy = geom.cy, R = R, ax = geom.ax, L = L, inner = inner,
    width = g.win_w, height = g.win_h, items = g.items, ntracks = 0, hidden = 0,
    history = g.history, off_seq = g.off_seq, named = named,
    canvas_w = geom.W, enz_labels = {},
    bin = bin, bands = bands, composition = comp, density = dens, landmarks = marks,
  }
  return setmetatable(layout, { __index = M })
end

---@param rec table Record
---@param opts table plasmid config options
---@param win_w integer
---@param win_h integer
---@param extra { selection?: table }|nil
function M.build(rec, opts, win_w, win_h, extra)
  extra = extra or {}
  local o = opts.circular
  local L = math.max(rec:len(), 1)
  local W = math.max(20, win_w)
  local H = math.max(8, win_h - HEADER_ROWS)
  local cv = Canvas.new(W, H)
  local ax = o.aspect or 1.0
  local cx, cy = cv.dw / 2, cv.dh / 2

  -- radius (in vertical-dot units) leaving room for tick text and label columns
  local label_w = math.min(o.label_width, math.floor(W * 0.24))
  local R = math.min(cy - 7, (cx - 2 * (label_w + 3)) * ax)
  local with_labels = true
  if R < 14 then
    with_labels = false
    R = math.min(cy - 7, (cx - 4) * ax)
  end
  R = math.max(R, 6)

  local function pt(theta, r)
    return cx + r * math.sin(theta) / ax, cy - r * math.cos(theta)
  end
  local function ang(pos) return TAU * (pos - 1) / L end

  -- the features the map has to account for, at either scale
  local items = vf.visible(rec, vim.tbl_extend("force", opts, {
    hidden_types = vim.tbl_extend("force", opts.hidden_types or {}, o.hidden_types or {}),
  }))
  -- features the linear view shows but the map omits (annotation history)
  local history = #vf.visible(rec, opts) - #items
  -- a feature that lies beyond the sequence cannot be placed on it (records
  -- without residues, truncated files): leave it off the map, and count it,
  -- because a feature the map does not draw is a feature the header owes the
  -- reader an account of
  local off_seq = 0
  do
    local n = rec:len()
    local keep = {}
    for _, it in ipairs(items) do
      local inside = n > 0
      for _, e in ipairs(it.extents or {}) do
        if e[1] < 1 or e[2] > n then inside = false end
      end
      if inside then keep[#keep + 1] = it else off_seq = off_seq + 1 end
    end
    items = keep
  end
  -- The backbone, the ticks and the overlays drawn on top of them are the
  -- frame at both scales; only what sits inside the backbone differs.
  local steps = math.ceil(TAU * R * 2)
  local function draw_backbone()
    for k = 0, steps - 1 do
      local th = TAU * k / steps
      local x, y = pt(th, R)
      cv:dot(x, y, "PlasmidBackbone")
      x, y = pt(th, R - 1)
      cv:dot(x, y, "PlasmidBackbone")
    end
  end

  local function draw_ticks(step)
    local pos = 0
    while pos < L do
      local th = ang(pos + 1)
      local x0, y0 = pt(th, R + 1)
      local x1, y1 = pt(th, R + 3)
      cv:line(x0, y0, x1, y1, "PlasmidTick")
      local tx, ty = pt(th, R + 6)
      local label = pos == 0 and "1" or fmt_bp(pos)
      local col, row = math.floor(tx / 2), math.floor(ty / 4)
      local s = math.sin(th)
      if s > 0.25 then
        -- right side: text starts at the point
      elseif s < -0.25 then
        col = col - #label + 1
      else
        col = col - math.floor(#label / 2)
      end
      if cv:text_free(col, row, #label) then cv:put(col, row, label, "PlasmidTick") end
      pos = pos + step
    end
  end

  -- what the caller asked to see on top of the map, in the band just outside
  -- the backbone where neither scale draws anything of its own
  local function draw_extras()
    -- selection arc just outside the backbone, painted over everything but text
    if extra.selection then
      local sel = extra.selection
      local a0 = ang(sel.start)
      local a1 = sel.start > sel.stop and (ang(sel.stop + 1) + TAU) or ang(sel.stop + 1)
      local n = math.max(2, math.ceil((a1 - a0) * (R + 2) * 2))
      for k = 0, n do
        local th = a0 + (a1 - a0) * k / n
        for _, r in ipairs({ R + 1, R + 1.5, R + 2 }) do
          local x, y = pt(th, r)
          cv:dot(x, y, "PlasmidSelectionArc")
        end
      end
    end

    -- search hits: short radial marks outside the backbone
    if extra.search and #extra.search.hits > 0 then
      for i, h in ipairs(extra.search.hits) do
        local group = i == extra.search.idx and "PlasmidSearchCurrent" or "PlasmidSearch"
        local a0 = ang(h.start)
        local a1 = h.stop >= h.start and ang(h.stop + 1) or (ang(h.stop + 1) + TAU)
        local n = math.max(1, math.ceil((a1 - a0) * (R + 3) * 2))
        for kk = 0, n do
          local th = a0 + (a1 - a0) * kk / n
          for _, r in ipairs({ R + 1, R + 2, R + 3 }) do
            local x, y = pt(th, r)
            cv:dot(x, y, group)
          end
        end
      end
    end

    -- restriction cuts: radial ticks across the backbone
    for _, g in ipairs(extra.enzymes or {}) do
      local th = ang(g.cut)
      local x0, y0 = pt(th, R - 2)
      local x1, y1 = pt(th, R + 4)
      cv:line(x0, y0, x1, y1, "PlasmidEnzymeCut")
    end
  end

  --- Name, size and a line or two more in the hole in the middle, as much as
  --- a circle of radius `inner` has room for.
  local function centre_text(inner, more)
    local inner_cols = math.floor(2 * inner / ax / 2) - 2
    local crow, ccol = math.floor(cy / 4), math.floor(cx / 2)
    if inner_cols < 6 then return end
    local name = utf8_trunc(rec.name, inner_cols)
    local size = ("%d %s"):format(rec:len(), rec:unit())
    cv:put(ccol - math.floor(utf8_len(name) / 2), crow - 1, name, "PlasmidCenter")
    cv:put(ccol - math.floor(#size / 2), crow, size, "PlasmidHeader")
    -- A further line needs a hole of its own to sit in, but may overhang it
    -- by a column or two: `inner_cols` is the width of a circle, and these
    -- rows are near its middle, where it is widest.
    for i, text in ipairs(more or {}) do
      if inner_cols >= 10 and #text <= inner_cols + 2 then
        cv:put(ccol - math.floor(#text / 2), crow + i, text, "PlasmidHeader")
      end
    end
  end

  local geom = { cx = cx, cy = cy, ax = ax, R = R, W = W, H = H, pt = pt }

  if ring.wanted(rec, opts) then
    return build_atlas(rec, opts, {
      cv = cv, geom = geom, ang = ang, L = L, steps = steps, label_w = label_w,
      with_labels = with_labels, items = items, history = history, off_seq = off_seq,
      win_w = win_w, win_h = win_h, draw_backbone = draw_backbone,
      draw_ticks = draw_ticks, draw_extras = draw_extras, centre_text = centre_text,
    })
  end

  local ring_gap = 4.5
  local ntracks = vf.pack(items, L, math.ceil(L * 0.004), opts)
  local want_rings = o.max_tracks == "all" and math.huge or o.max_tracks
  local max_tracks = math.max(1, math.min(want_rings, math.floor((R * 0.6 - 5) / ring_gap) + 1))
  local function ring_r(t) return R - 5 - ring_gap * (t - 1) end

  draw_backbone()
  if o.show_ticks then draw_ticks(ring.tick_step(L, "arcs")) end

  -- feature arcs (drawn after leaders would hide them, so collect first)
  local drawn, hidden = {}, 0
  for _, it in ipairs(items) do
    if it.track <= max_tracks then drawn[#drawn + 1] = it else hidden = hidden + 1 end
  end

  local function arc(r, a0, a1, group, owner, hw_fn)
    local span = a1 - a0
    local n = math.max(2, math.ceil(span * r * 2))
    for k = 0, n do
      local th = a0 + span * k / n
      local hw = hw_fn and hw_fn(th) or 1
      local d = -hw
      while d <= hw + 1e-9 do
        local x, y = pt(th, r + d)
        cv:dot(x, y, group, owner)
        d = d + 0.5
      end
    end
  end

  local function arc_mid(it)
    -- angular midpoint over the feature's covered arc (handles origin wrap)
    local a_start, a_len
    if #it.extents > 1 then
      local tail, head
      for _, e in ipairs(it.extents) do
        if e[1] == 1 then head = e else tail = e end
      end
      a_start = ang(tail[1])
      a_len = (L - tail[1] + 1 + head[2]) / L * TAU
    else
      local e = it.extents[1]
      a_start = ang(e[1])
      a_len = (e[2] - e[1] + 1) / L * TAU
    end
    return (a_start + a_len / 2) % TAU
  end

  -- labels: hug the circle just outside the ring; placed in priority order
  -- into the nearest free row, and dropped (not routed far away) if none is
  -- close. Leaders are drawn now so feature arcs paint over their inner ends.
  local labels
  local enz_labels = {}
  if with_labels then
    labels = Labels.new(cv, { cx = cx, cy = cy, ax = ax, R = R, W = W, H = H, pt = pt }, label_w)
    local cand = {}
    for _, it in ipairs(drawn) do
      if it.label ~= it.f.type then cand[#cand + 1] = it end
    end
    table.sort(cand, function(a, b)
      local pa, pb = PRIORITY[a.f.type] or 8, PRIORITY[b.f.type] or 8
      if pa ~= pb then return pa < pb end
      local la, lb = a.f:length(), b.f:length()
      if la ~= lb then return la > lb end
      return a.f.id < b.f.id
    end)
    local tolerance = math.max(2, math.floor(H / 8))
    local labelled = {}
    for _, it in ipairs(cand) do
      if not labelled[it.label] then
        local placed = labels:place(arc_mid(it), it.label, it.hl, it.f.id,
          { tolerance = tolerance, from = ring_r(it.track) + 2, hinge = R + 3 })
        if placed then labelled[it.label] = true end
      end
    end

    -- restriction-site labels, after the features: unique cutters first, fanned
    -- out further than feature labels since sites cluster in polylinkers
    if extra.enzymes and #extra.enzymes > 0 then
      local ecand = {}
      for _, g in ipairs(extra.enzymes) do ecand[#ecand + 1] = g end
      table.sort(ecand, function(a, b)
        if a.count ~= b.count then return a.count < b.count end
        return a.cut < b.cut
      end)
      local etol = math.max(3, math.floor(H / 3))
      for _, g in ipairs(ecand) do
        local grp = g.count == 1 and "PlasmidEnzymeUnique" or "PlasmidEnzyme"
        local placed = labels:place(ang(g.cut), g.label, grp, nil, { tolerance = etol, hinge = R + 4 })
        if placed then
          enz_labels[#enz_labels + 1] =
            { col = placed.col, row = placed.row, len = utf8_len(placed.text), entry = g }
        end
      end
    end
  end

  for _, it in ipairs(drawn) do
    local r = ring_r(it.track)
    local f = it.f
    local apos, strand = vf.arrow(f)
    for _, g in ipairs(vf.parts_by_extent(it)) do
      local prev
      for _, p in ipairs(g) do
        local a0, a1
        if p.stop < p.start then
          a0 = ang(p.start)
          a1 = a0
        else
          a0, a1 = ang(p.start), ang(p.stop + 1)
        end
        local min_span = 1.5 / r
        if a1 - a0 < min_span then
          local m = (a0 + a1) / 2
          a0, a1 = m - min_span / 2, m + min_span / 2
        end
        local hw_fn
        local has_arrow = apos and ((strand == 1 and p.stop == apos) or (strand == -1 and p.start == apos))
        if has_arrow then
          local alen = math.min((a1 - a0) * 0.5, 5 / r)
          local tip = strand == 1 and a1 or a0
          hw_fn = function(th)
            local d = math.abs(tip - th)
            if d >= alen then return 1 end
            return 2.5 * (d / alen)
          end
        end
        arc(r, a0, a1, it.hl, f.id, hw_fn)
        if prev and p.start > prev.stop + 1 then
          local b0, b1 = ang(prev.stop + 1), ang(p.start)
          local n = math.max(2, math.ceil((b1 - b0) * r))
          for k = 0, n, 2 do
            local x, y = pt(b0 + (b1 - b0) * k / n, r)
            cv:dot(x, y, it.hl, f.id)
          end
        end
        if not prev or p.stop > prev.stop then prev = p end
      end
    end
  end

  draw_extras()

  if labels then labels:repaint() end

  -- centre text
  local inner = ring_r(math.min(ntracks, max_tracks)) - 3
  centre_text(inner, rec.topology ~= "circular" and { "(linear DNA)" } or nil)

  local clines, cruns, cbytes, cowners = cv:render()
  local lines = {
    (" %s  ·  %d %s  ·  %s  ·  %s  ·  %s%s"):format(
      rec.name, rec:len(), rec:unit(), rec.molecule_type or "DNA", rec.topology,
      vf.count_text(items, history + off_seq),
      (history + off_seq + hidden) > 0 and ("  (%s not on map)"):format(table.concat(vim.tbl_filter(function(x) return x end, {
        history > 0 and ("%d annotation-history"):format(history) or false,
        off_seq > 0 and ("%d past the sequence end"):format(off_seq) or false,
        hidden > 0 and ("%d over ring limit"):format(hidden) or false,
      }), ", ")) or ""),
    "",
  }
  for _, l in ipairs(clines) do lines[#lines + 1] = l end

  local layout = {
    view = "circular", lines = lines, header_rows = HEADER_ROWS,
    runs = cruns, cell_bytes = cbytes, owners = cowners,
    cx = cx, cy = cy, R = R, ax = ax, L = L, inner = inner,
    width = win_w, height = win_h, items = items, ntracks = ntracks, hidden = hidden,
    history = history, off_seq = off_seq,
    canvas_w = W, enz_labels = enz_labels,
  }
  return setmetatable(layout, { __index = M })
end

---------------------------------------------------------------------------
-- Layout queries
---------------------------------------------------------------------------

--- Base pair at a cell, by angle from the centre. Returns bp, radius.
function M.bp_at(self, row0, vcol0)
  local r = row0 - self.header_rows
  if r < 0 then return nil end
  local x, y = vcol0 * 2 + 1, r * 4 + 2
  local dx, dy = (x - self.cx) * self.ax, self.cy - y
  local th = math.atan2(dx, dy)
  if th < 0 then th = th + TAU end
  local bp = math.floor(th / TAU * self.L) + 1
  return math.max(1, math.min(bp, self.L)), math.sqrt(dx * dx + dy * dy)
end

--- Feature id drawn at a cell, if any.
function M.feature_at(self, row0, vcol0)
  local r = row0 - self.header_rows
  local ow = self.owners[r + 1]
  return ow and ow[vcol0] or nil
end

--- Cell on the backbone for a base pair: row0, display col0.
function M.pos_of(self, bp)
  local th = TAU * (bp - 0.5) / self.L
  local x = self.cx + self.R * math.sin(th) / self.ax
  local y = self.cy - self.R * math.cos(th)
  local col = math.max(0, math.min(self.canvas_w - 1, math.floor(x / 2)))
  return self.header_rows + math.max(0, math.floor(y / 4)), col
end

--- Byte offset of a display column on a row (rows contain multibyte braille).
function M.byte_col(self, row0, vcol0)
  local cb = self.cell_bytes[row0 - self.header_rows + 1]
  if not cb then return vcol0 end
  return cb[vcol0] or cb[#cb] or 0
end

function M.row_highlights(self, row0)
  if row0 == 0 then return { { 0, -1, "PlasmidTitle" } } end
  if row0 == 1 and self.legend_runs then return self.legend_runs end
  local r = row0 - self.header_rows
  return self.runs[r + 1] or {}
end

return M
