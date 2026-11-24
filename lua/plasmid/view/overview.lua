-- Whole-record overview: the shape of a record on one screen.
--
-- The linear view answers "what is at this base"; at 250,000 rows for a
-- bacterial genome it cannot answer "what does this record look like" or
-- "where should I go". This view gives up the residues and keeps the shape:
-- one row per fixed span, each row a strip of cells coloured by the feature
-- that matters most in it.
--
--   2300001..2400000 │▀▀█▄▄·····▀▀│  14+   9−  nuoA, nuoB, lrhA +11
--
-- A cell is the upper half block when only the forward strand has features
-- there, the lower half when only the reverse strand does, and a full block
-- when both do, so strand reads as position in the cell rather than as a
-- colour the feature types already use.
--
-- Nothing here touches a buffer: build() returns lines, per-row metadata and
-- highlight runs, which is also what the specs work on.
local vf = require("plasmid.view.features")
local hl = require("plasmid.view.highlights")

local M = {}

local FWD, REV, BOTH, EMPTY = "▀", "▄", "█", "·"
local BAR = "│"

-- Bits in a cell's strand mask. A feature with no strand (strand 0) claims
-- both halves: it is on the molecule, not on one of its strands.
local B_FWD, B_REV = 1, 2

--- Number of display cells in a UTF-8 string (continuation bytes removed).
local function width_of(s)
  return #(s:gsub("[\128-\191]", ""))
end

--- The next round number at or above `n`: 1, 2 or 5 times a power of ten.
--- These are the spans a person reads without converting them, which is the
--- whole point of a row that stands for 200 kb of DNA.
function M.ladder(n)
  n = math.max(1, math.ceil(n))
  local k = 1
  while k * 10 <= n do k = k * 10 end
  for _, m in ipairs({ 1, 2, 5 }) do
    if m * k >= n then return m * k end
  end
  return 10 * k
end

--- Bases per row: enough that the record fits in `rows` of them, rounded up
--- to a round number. A configured `overview.span` overrides it, in which
--- case the view may well be longer than the screen and simply scrolls.
function M.span_for(L, rows, opts)
  local o = (opts and opts.overview) or {}
  local fixed = tonumber(o.span)
  if fixed and fixed >= 1 then return math.floor(fixed) end
  local floor_span = math.max(1, math.floor(tonumber(o.min_span) or 10))
  rows = math.max(1, math.floor(rows or 1))
  return math.max(floor_span, M.ladder(math.max(L, 1) / rows))
end

--- A span as a person would say it ("200 kb", "2 Mb", "500 bp"). Only base
--- pairs get the kb/Mb prefixes; nobody measures a protein in kilo-residues.
function M.span_text(span, unit)
  unit = unit or "bp"
  if unit == "bp" then
    if span >= 1000000 and span % 1000000 == 0 then return ("%d Mb"):format(span / 1000000) end
    if span >= 1000 and span % 1000 == 0 then return ("%d kb"):format(span / 1000) end
  end
  return ("%d %s"):format(span, unit)
end

---------------------------------------------------------------------------
-- Building
---------------------------------------------------------------------------

---@param rec table Record
---@param opts table plasmid config options
---@param width integer columns available
---@param height integer text rows available (the record is made to fit)
---@param extra table|nil { cursor_bp = integer, visible = { lo, hi } }
function M.build(rec, opts, width, height, extra)
  extra = extra or {}
  local o = opts.overview or {}
  local L = rec:len()
  local unit = rec:unit()
  width = math.max(40, math.floor(width or 80))
  local header_rows = 3
  local span = M.span_for(L, math.max(1, (height or 40) - header_rows - 1), opts)
  local nrows = math.max(1, math.ceil(math.max(L, 1) / span))

  -- Column budget. The coordinates and the per-strand counts are fixed; what
  -- is left is split between the strip and the names, with the strip taking
  -- the smaller share because a name is only useful whole.
  local numw = #tostring(math.max(L, 1))
  local coord_w = 2 * numw + 2
  local counts = o.counts ~= false
  local counts_w = counts and 12 or 0
  local gutter_w = 2
  local fixed_w = gutter_w + coord_w + 2 + 2 + counts_w + 1
  local avail = math.max(8, width - fixed_w)
  local cells = tonumber(o.strip) or math.floor(avail * 0.45)
  cells = math.max(8, math.min(cells, avail - 1))
  local labels = o.labels ~= false and (avail - cells) >= 12
  local label_w = labels and (avail - cells - 1) or 0

  local items = vf.visible(rec, opts)
  for _, it in ipairs(items) do vf.priority(it, opts) end

  -- Bin every feature into the rows it reaches. One pass over the features,
  -- and a feature only ever touches the cells it actually covers, so this is
  -- proportional to the picture rather than to the record's length.
  local grid, seen = {}, {}
  local function row_of(pos) return math.floor((math.max(pos, 1) - 1) / span) end
  for _, it in ipairs(items) do
    local strand = it.f.strand
    local bit = (strand == 1 and B_FWD) or (strand == -1 and B_REV) or (B_FWD + B_REV)
    for _, ext in ipairs(it.extents) do
      local s, e = math.max(ext[1], 1), math.min(ext[2], math.max(L, 1))
      if s <= e then
        for r = row_of(s), row_of(e) do
          local g = grid[r]
          if not g then
            g = { fwd = 0, rev = 0, other = 0, items = {}, mask = {}, prio = {}, hl = {} }
            grid[r] = g
          end
          -- an origin-spanning feature reaches a row through two extents;
          -- it is still one feature in that row
          if seen[it] ~= r then
            seen[it] = r
            g.items[#g.items + 1] = it
            if strand == 1 then g.fwd = g.fwd + 1
            elseif strand == -1 then g.rev = g.rev + 1
            else g.other = g.other + 1 end
          end
          local bs = r * span + 1
          local c0 = math.floor((math.max(s, bs) - bs) * cells / span)
          local c1 = math.floor((math.min(e, bs + span - 1) - bs) * cells / span)
          for c = math.max(c0, 0), math.min(c1, cells - 1) do
            local m = g.mask[c]
            g.mask[c] = (m == nil and bit) or ((m == bit) and m) or (B_FWD + B_REV)
            if (g.prio[c] or -math.huge) < it.priority then
              g.prio[c], g.hl[c] = it.priority, it.hl
            end
          end
        end
      end
    end
  end

  local lines, rows, hls = {}, {}, {}
  local function push(text, info)
    lines[#lines + 1] = text
    rows[#lines - 1] = info
  end
  local function mark(row0, b0, b1, group, prio)
    hls[#hls + 1] = { row0, b0, b1, group, prio }
  end

  -- header
  local kind = rec:is_protein() and "protein"
    or ((rec.molecule_type or "DNA") .. "  ·  " .. rec.topology)
  push((" %s  ·  %d %s  ·  %s  ·  %s"):format(rec.name, L, unit, kind, vf.count_text(items)),
    { kind = "title" })
  mark(0, 0, #lines[1], "PlasmidTitle")
  local hint = (" %s per row  ·  %d rows  ·  %s jump to this region, %s back"):format(
    M.span_text(span, unit), nrows,
    (o.keymaps or {}).jump or "<CR>", (o.keymaps or {}).close or "q")
  push(hint, { kind = "header" })
  mark(1, 0, #hint, "PlasmidHeader")
  push("", { kind = "blank" })

  local cursor_bp = tonumber(extra.cursor_bp)
  local vis = extra.visible
  local cursor_row = cursor_bp and row_of(math.min(cursor_bp, math.max(L, 1)))

  for r = 0, nrows - 1 do
    local bs, be = r * span + 1, math.min((r + 1) * span, math.max(L, 1))
    local g = grid[r] or { fwd = 0, rev = 0, other = 0, items = {}, mask = {}, hl = {} }
    local row0 = #lines
    local parts, byte = {}, 0
    local function add(s, group, prio)
      parts[#parts + 1] = s
      if group then mark(row0, byte, byte + #s, group, prio) end
      byte = byte + #s
    end

    -- gutter: where the sequence view is looking, and where its cursor is
    local gut = " "
    if cursor_row == r then gut = "▶"
    elseif vis and be >= vis[1] and bs <= vis[2] then gut = "┃" end
    add(gut ~= " " and (gut .. " ") or "  ", gut ~= " " and "PlasmidOverviewViewport" or nil)

    add(("%" .. numw .. "d..%-" .. numw .. "d"):format(bs, be), "PlasmidOverviewCoord")
    add("  ")
    add(BAR, "PlasmidOverviewAxis")
    local strip0 = byte
    local run_hl, run_b0 = nil, byte
    for c = 0, cells - 1 do
      local m = g.mask[c]
      local ch = (m == B_FWD and FWD) or (m == B_REV and REV) or (m and BOTH) or EMPTY
      local group = m and (g.hl[c] or "PlasmidSequence") or "PlasmidOverviewEmpty"
      if group ~= run_hl then
        if run_hl then mark(row0, run_b0, byte, run_hl) end
        run_hl, run_b0 = group, byte
      end
      parts[#parts + 1] = ch
      byte = byte + #ch
    end
    if run_hl then mark(row0, run_b0, byte, run_hl) end
    local strip1 = byte
    add(BAR, "PlasmidOverviewAxis")

    if counts then
      add(("%4d+ %4d−"):format(math.min(g.fwd + g.other, 9999), math.min(g.rev + g.other, 9999)),
        "PlasmidOverviewCount")
    end

    if labels and #g.items > 0 then
      table.sort(g.items, function(a, b)
        if a.priority ~= b.priority then return a.priority > b.priority end
        return a.f.id < b.f.id
      end)
      -- a gene and its CDS, or a CDS and its exons, carry the same name; the
      -- row has room for one of each name and a count of what is left over
      local distinct, seen_name = {}, {}
      for _, it in ipairs(g.items) do
        local nm = it.label
        if nm and nm ~= "" and not seen_name[nm] then
          seen_name[nm] = true
          distinct[#distinct + 1] = nm
        end
      end
      local shown, used = 0, 0
      for _, nm in ipairs(distinct) do
        local w = width_of(nm) + (shown > 0 and 2 or 0)
        if used + w > label_w - 4 then break end
        shown, used = shown + 1, used + w
      end
      if shown > 0 then
        local text = " " .. table.concat(distinct, ", ", 1, shown)
        if shown < #distinct then text = text .. (" +%d"):format(#distinct - shown) end
        add(text, "PlasmidOverviewLabel")
      end
    end

    -- the cursor's own cell, over the feature colour
    if cursor_row == r and cursor_bp then
      local c = math.max(0, math.min(math.floor((cursor_bp - bs) * cells / span), cells - 1))
      mark(row0, strip0 + c * #EMPTY, strip0 + (c + 1) * #EMPTY, "PlasmidOverviewCursor", 200)
    end

    push(table.concat(parts), {
      kind = "region", r = r, bs = bs, be = be,
      fwd = g.fwd, rev = g.rev, other = g.other, items = g.items,
      strip0 = strip0, strip1 = strip1, strip_col = gutter_w + coord_w + 2 + 1,
    })
  end

  if L == 0 then
    push("", { kind = "blank" })
    push("  (no residues in this record)", { kind = "header" })
  end

  return setmetatable({
    view = "overview", lines = lines, rows = rows, hls = hls,
    span = span, nrows = nrows, cells = cells, header_rows = header_rows,
    L = L, width = width, unit = unit, labels = labels,
  }, { __index = M })
end

---------------------------------------------------------------------------
-- Row <-> region
---------------------------------------------------------------------------

--- First and last base of the region on a buffer row (nil off the strip).
function M.region_of_row(self, row0)
  local info = self.rows[row0]
  if not info or info.kind ~= "region" then return nil end
  return info.bs, info.be
end

--- Buffer row (0-based) whose region contains a base.
function M.row_of_bp(self, bp)
  bp = math.max(1, math.min(math.floor(bp or 1), math.max(self.L, 1)))
  return self.header_rows + math.min(math.floor((bp - 1) / self.span), self.nrows - 1)
end

--- Base at a display column of a row: a column inside the strip is worth
--- span/cells bases, so a click lands nearer than the start of the region.
function M.bp_at(self, row0, vcol0)
  local info = self.rows[row0]
  if not info or info.kind ~= "region" then return nil end
  local c0 = info.strip_col
  local col = math.floor(vcol0 or 0)
  if col < c0 or col >= c0 + self.cells then return info.bs end
  local bp = info.bs + math.floor((col - c0) * self.span / self.cells)
  return math.max(info.bs, math.min(bp, info.be))
end

return M
