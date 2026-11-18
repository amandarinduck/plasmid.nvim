-- View-side feature helpers shared by the linear and circular views:
-- which features to draw, their extents, track packing, and codon maps.
local location = require("plasmid.model.location")
local translate = require("plasmid.bio.translate")
local alphabet = require("plasmid.model.alphabet")

local M = {}

local GENE_PRODUCTS = { CDS = true, mRNA = true, tRNA = true, rRNA = true, ncRNA = true, misc_RNA = true }

--- Features worth drawing, each wrapped as { f, extents, label, hl }.
---@return table[] items sorted by first extent start
function M.visible(rec, opts)
  local hl = require("plasmid.view.highlights")
  local hidden = opts.hidden_types or {}
  local product_locs = {}
  if opts.dedupe_genes then
    for _, f in ipairs(rec.features) do
      if GENE_PRODUCTS[f.type] then product_locs[f.location_text] = true end
    end
  end
  local items = {}
  local all = rec.features
  if rec.overlays and #rec.overlays > 0 then
    all = {}
    for _, f in ipairs(rec.features) do all[#all + 1] = f end
    for _, f in ipairs(rec.overlays) do all[#all + 1] = f end
  end
  for _, f in ipairs(all) do
    if not hidden[f.type] and f.loc
      and not (f.type == "gene" and product_locs[f.location_text]) then
      local ext = rec:extents(f)
      if #ext > 0 then
        items[#items + 1] = {
          f = f, extents = ext, label = f:label(opts.label_qualifiers),
          hl = f.detected and "PlasmidProposed" or hl.for_feature(f),
        }
      end
    end
  end
  table.sort(items, function(a, b)
    local sa, sb = a.extents[#a.extents][1], b.extents[#b.extents][1]
    -- wrapped features sort by their origin-side extent {1, x}
    if #a.extents > 1 then sa = 1 end
    if #b.extents > 1 then sb = 1 end
    if sa ~= sb then return sa < sb end
    return math.abs(a.f.id) < math.abs(b.f.id) or (math.abs(a.f.id) == math.abs(b.f.id) and a.f.id > b.f.id)
  end)
  return items
end

--- Display priority of an item (higher wins space). Sets/returns item.priority.
function M.priority(item, opts)
  if item.priority then return item.priority end
  local p = (opts and opts.priority) or require("plasmid.config").options.priority
  local f = item.f
  local score = p.types[f.type] or p.default
  if item.label and item.label ~= f.type then score = score + p.labelled_bonus end
  if f.pinned then score = score + p.pinned_bonus end
  item.priority = score
  return score
end

--- Sort items by priority (desc), then length (desc), then id.
function M.by_priority(items, opts)
  local out = {}
  for i, it in ipairs(items) do
    M.priority(it, opts)
    it._len = it._len or it.f:length()
    out[i] = it
  end
  table.sort(out, function(a, b)
    if a.priority ~= b.priority then return a.priority > b.priority end
    if a._len ~= b._len then return a._len > b._len end
    return a.f.id < b.f.id
  end)
  return out
end

-- Interval set with overlap query, kept sorted by start. Intervals in a set
-- never overlap each other, so the predecessor-by-start is the only candidate
-- that can reach into a query from the left.
local function overlaps(set, s, e)
  local lo, hi = 1, #set
  local pos = #set + 1
  while lo <= hi do
    local mid = math.floor((lo + hi) / 2)
    if set[mid][1] > e then pos = mid; hi = mid - 1 else lo = mid + 1 end
  end
  -- set[pos] is the first interval starting after e; set[pos-1] starts <= e
  local prev = set[pos - 1]
  return prev ~= nil and prev[2] >= s
end

local function insert(set, s, e)
  local i = #set + 1
  for k = 1, #set do
    if set[k][1] > s then i = k break end
  end
  table.insert(set, i, { s, e })
end

M._overlaps, M._insert = overlaps, insert

--- Assign tracks in priority order: the most important features claim the
--- first tracks, so lower-priority ones are the ones displaced outward.
--- Sets item.track (1-based) and returns the number of tracks.
---@param gap integer minimum bp between features sharing a track
function M.pack(items, seq_len, gap, opts)
  gap = gap or 0
  local tracks = {}
  for _, it in ipairs(M.by_priority(items, opts)) do
    local placed
    for t, set in ipairs(tracks) do
      local free = true
      for _, e in ipairs(it.extents) do
        if overlaps(set, e[1] - gap, e[2] + gap) then free = false break end
      end
      if free then placed = t break end
    end
    if not placed then
      tracks[#tracks + 1] = {}
      placed = #tracks
    end
    for _, e in ipairs(it.extents) do insert(tracks[placed], e[1], e[2]) end
    it.track = placed
  end
  return #tracks
end

--- Arrow position: the 3' end on the forward axis, or nil for mixed strand.
---@return integer|nil pos, integer strand
function M.arrow(f)
  if f.strand == 0 or f.no_arrow then return nil, 0 end
  local last
  for i = #f.parts, 1, -1 do
    if not f.parts[i].ref and f.parts[i].stop >= f.parts[i].start then
      last = f.parts[i]
      break
    end
  end
  if not last then return nil, f.strand end
  if f.strand == -1 then return last.start, -1 end
  return last.stop, 1
end

--- Local parts sorted left-to-right, grouped by extent (so introns are
--- drawn only inside an extent, never across the origin).
function M.parts_by_extent(item)
  local groups = {}
  for gi, ext in ipairs(item.extents) do
    local g = {}
    for _, p in ipairs(item.f.parts) do
      if not p.ref and p.start <= ext[2] and math.max(p.stop, p.start) >= ext[1] then g[#g + 1] = p end
    end
    table.sort(g, function(a, b) return a.start < b.start end)
    groups[gi] = g
  end
  return groups
end

--- For a CDS, map genomic position of each codon's middle base -> amino acid.
--- Cached on the feature; invalidated by record revision.
function M.codon_map(rec, f)
  if f._codon_map and f._codon_rev == rec.revision then return f._codon_map end
  local map, idx = {}, {}
  f._codon_idx, f._codon_pos = idx, nil
  if rec.alphabet == "protein" or location.has_remote(f.parts) then
    f._codon_map, f._codon_rev = map, rec.revision
    return map
  end
  local t = translate.get_table(f:get("transl_table"))
  local pos = location.positions(f.parts)
  local cs = tonumber(f:get("codon_start")) or 1
  local first = f.parts[1]
  local five_partial = first and ((f.strand == -1 and first.stop_fuzzy) or (f.strand ~= -1 and first.start_fuzzy))
  local comp = alphabet.complement_map
  local s = rec.seq
  local k = 0
  for i = cs, #pos - 2, 3 do
    local p1, p2, p3 = pos[i], pos[i + 1], pos[i + 2]
    local b1, b2, b3 = s:sub(p1, p1), s:sub(p2, p2), s:sub(p3, p3)
    if f.strand == -1 then
      b1, b2, b3 = comp[b1:byte()], comp[b2:byte()], comp[b3:byte()]
    end
    local codon = b1 .. b2 .. b3
    local aa = translate.codon(t, codon)
    if k == 0 and not five_partial and cs == 1 and translate.is_start(t, codon) then aa = "M" end
    map[p2] = aa
    idx[p2] = i
    k = k + 1
  end
  f._codon_map, f._codon_rev, f._codon_pos = map, rec.revision, pos
  return map
end

--- The codon whose middle base is at `p2` in CDS/ORF `f`:
--- { n = codon number, total, positions = {p1, p2, p3} (5'->3' on f's strand),
---   codon = "ATG" (sense strand), aa = "M" } or nil.
function M.codon_at(rec, f, p2)
  local map = M.codon_map(rec, f)
  local i = f._codon_idx and f._codon_idx[p2]
  if not i then return nil end
  local pos = f._codon_pos
  local cs = tonumber(f:get("codon_start")) or 1
  local p1, p3 = pos[i], pos[i + 2]
  local s = rec.seq
  local codon = s:sub(p1, p1) .. s:sub(p2, p2) .. s:sub(p3, p3)
  if f.strand == -1 then codon = alphabet.complement(codon) end
  return {
    n = (i - cs) / 3 + 1, total = math.floor((#pos - cs + 1) / 3),
    positions = { p1, p2, p3 }, codon = codon:upper(), aa = map[p2],
  }
end


--- Header text for the item list: "N features" plus overlay counts, e.g.
--- "0 features  ·  11 proposed  ·  3 ORFs". `extra` adds to the annotated count.
function M.count_text(items, extra)
  local n, prop, orfs = extra or 0, 0, 0
  for _, it in ipairs(items) do
    if it.f.detected then prop = prop + 1
    elseif it.f.orf then orfs = orfs + 1
    else n = n + 1 end
  end
  local s = ("%d feature%s"):format(n, n == 1 and "" or "s")
  if prop > 0 then s = s .. ("  ·  %d proposed"):format(prop) end
  if orfs > 0 then s = s .. ("  ·  %d ORF%s"):format(orfs, orfs == 1 and "" or "s") end
  return s
end

return M
