-- Dev tool: dump selected windows of a linear layout as JSON, one layout mode
-- per process.
--
--   nvim -l tests/tools/big_dump.lua FILE OUT.json [key=value ...]
--
-- The eager layout of a bacterial genome is half a gigabyte, so the windowed
-- and the eager layout of the same record are never built in one process:
-- each writes its own dump and tests/tools/big_parity.py compares them. The
-- dump holds every row of the chosen windows -- text, kind, span, highlight
-- runs, feature segments -- plus the position queries answered for those rows,
-- which is everything the view shows for a screenful.
--
-- keys: mode=window|eager  rec=1  G=10  width=120  tracks=4|all
--       enz=none|unique|all  primers=0|1  orfs=0|MIN_AA  search=MOTIF
--       sel=START..STOP  churn=0|1  blocks=12  windows=start,mid,end,dense,origin
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({ root .. "/lua/?.lua", root .. "/lua/?/init.lua", package.path }, ";")
vim.opt.runtimepath:prepend(root)
vim.notify = function() end

local config = require("plasmid.config")
config.setup({})
local io_ = require("plasmid.io")
local linear = require("plasmid.view.linear")

local path, out_path = arg[1], arg[2]
local A = {}
for i = 3, #arg do
  local k, v = arg[i]:match("^([%w_]+)=(.*)$")
  if k then A[k] = v end
end
local function num(k, d) return tonumber(A[k]) or d end

local t_parse = vim.uv.hrtime()
local res = assert(io_.read(path), "could not read " .. tostring(path))
local parse_ms = (vim.uv.hrtime() - t_parse) / 1e6
local rec = res.records[num("rec", 1)]

local o = vim.deepcopy(config.options)
o.linear.group = num("G", 10)
o.linear.max_tracks = A.tracks == "all" and "all" or num("tracks", 4)
local windowed = (A.mode or "window") == "window"
-- margin 0 and a cache big enough for the whole record keep the dump about the
-- layout rather than about the eviction policy; churn=1 exercises eviction.
o.linear.window = { enabled = windowed, margin = 0, cache_blocks = num("cache", 1e9) }

---------------------------------------------------------------------------
-- Overlays
---------------------------------------------------------------------------
local extra = {}
if A.enz and A.enz ~= "none" and not rec:is_protein() then
  local entries = require("plasmid.enzymes").compute(rec, { mode = A.enz })
  if #entries > 0 then extra.enzymes = entries end
end
if A.primers == "1" and not rec:is_protein() then
  local primers = require("plasmid.primers").list(rec)
  if #primers > 0 then
    extra.primers = primers
    extra.primer_cells = function(p) return require("plasmid.primers").cells(rec, p) end
  end
end
if A.wrap == "1" then
  -- Real genomes rarely annotate a feature across the origin, so the seam of a
  -- megabase circular record is only exercised if one is put there: a feature
  -- in two extents, drawn at the end of the last block and again at the start
  -- of the first. A linear record has no seam; the edge worth drawing a
  -- feature on there is its last block, which is the one that is short.
  local L = rec:len()
  rec:add_feature(require("plasmid.model.record").Feature.new({
    type = "CDS",
    location_text = rec:is_circular() and ("join(%d..%d,1..%d)"):format(L - 4000, L, 2600)
      or ("%d..%d"):format(L - 4000, L),
    qualifiers = { { key = "label", value = "edge_spanner", quoted = true } },
  }))
end
if A.orfs and A.orfs ~= "0" then
  rec.overlays, rec.overlay_by_id = require("plasmid.orfs").compute(rec, {
    min_aa = num("orfs", 100), starts = config.options.orfs.starts, hide_annotated = true,
  })
end

local t0 = vim.uv.hrtime()
local lay = linear.build(rec, o, num("width", 120), extra)
local build_ms = (vim.uv.hrtime() - t0) / 1e6
if windowed ~= (lay.windowed == true) then error("layout mode is not the one asked for") end

if A.sel then
  local s, e = A.sel:match("^(%d+)%.%.(%d+)$")
  lay.selection = { start = tonumber(s), stop = tonumber(e) }
end
if A.search then
  local searchui = require("plasmid.searchui")
  lay.search = searchui._compute and searchui._compute(rec, A.search, 0) or nil
  if not lay.search then
    local search = require("plasmid.bio.search")
    local hits = search.find(rec, A.search, 0)
    local wrapped, maxlen = {}, 1
    for i, h in ipairs(hits or {}) do
      if h.stop < h.start then wrapped[#wrapped + 1] = i else maxlen = math.max(maxlen, h.stop - h.start + 1) end
    end
    lay.search = { hits = hits or {}, wrapped = wrapped, maxlen = maxlen, idx = 1 }
  end
end

---------------------------------------------------------------------------
-- Which blocks to dump
---------------------------------------------------------------------------
local W = lay.W
local nblocks = math.max(1, math.ceil(lay.L / W))
local span = num("blocks", 24)

--- The block with the most feature extents starting in it: the part of a
--- record where the track packer has the most to do.
local function densest()
  local hist, best, best_n = {}, 0, -1
  for _, it in ipairs(lay.items) do
    for _, e in ipairs(it.extents) do
      local b = math.floor((math.max(e[1], 1) - 1) / W)
      hist[b] = (hist[b] or 0) + 1
    end
  end
  for b, n in pairs(hist) do
    if n > best_n then best, best_n = b, n end
  end
  return best
end

local wanted = {}
for _, name in ipairs(vim.split(A.windows or "start,mid,end,dense,origin", ",", { plain = true })) do
  local first
  if name == "start" then first = 0
  elseif name == "mid" then first = math.floor(nblocks / 2)
  elseif name == "end" then first = math.max(0, nblocks - span)
  elseif name == "dense" then first = math.max(0, densest() - math.floor(span / 2))
  elseif name == "origin" then
    -- the origin of a circular record is a seam: the last block and the first
    -- one, where a wrapped feature is drawn as two pieces. A linear record has
    -- no seam, so this is its last block instead -- asking for a window and
    -- getting none back would be a case that quietly compares nothing.
    if rec:is_circular() then
      wanted[#wanted + 1] = { name = "origin_end", lo = math.max(0, nblocks - math.ceil(span / 2)), hi = nblocks - 1 }
      first = 0
      name = "origin_start"
    else
      first = math.max(0, nblocks - span)
      name = "last_block"
    end
  elseif name:match("^%d+$") then first = math.min(tonumber(name), nblocks - 1)
  end
  if first then
    wanted[#wanted + 1] = { name = name, lo = first, hi = math.min(nblocks - 1, first + span - 1) }
  end
end

---------------------------------------------------------------------------
-- Position queries over the whole record, before anything is built. The
-- windowed layout answers these from the row plans alone, and a plan that
-- counted a block's rows wrongly shows up here as a row number that does not
-- match the eager layout's, anywhere in the record rather than in the windows
-- that follow.
---------------------------------------------------------------------------
local global = {}
do
  local step = math.max(1, math.floor(lay.L / 4000))
  for bp = 1, lay.L, step do
    local r, c = lay:pos_of(bp)
    global[#global + 1] = { bp, r, c, lay:bp_at(r, c), lay:block_at_row(r) }
  end
end
local built_after_positions = 0
if windowed then
  for _ in pairs(lay._built) do built_after_positions = built_after_positions + 1 end
end

---------------------------------------------------------------------------
-- Building the rows a window needs
---------------------------------------------------------------------------
local ensure_ms, churn_ms = 0, 0
if windowed then
  if A.churn == "1" then
    -- Build the windows, then scroll past them with a cache too small to keep
    -- them, so that what is compared below are rebuilt rows: a block that is
    -- laid out, evicted and laid out again has to come back the same.
    lay.cache_blocks = 64
    lay.margin = 1
    for _, w in ipairs(wanted) do
      for b = w.lo, w.hi do lay:build_block(b) end
    end
    local t = vim.uv.hrtime()
    local row, page = lay.header_rows, 40
    while row + page < #lay.lines do
      lay:ensure(row, row + page)
      row = row + page * 37   -- a coarse pass over the record, not every screen
    end
    churn_ms = (vim.uv.hrtime() - t) / 1e6
    local still = 0
    for _, w in ipairs(wanted) do
      for b = w.lo, w.hi do if lay._built[b] then still = still + 1 end end
    end
    io.stdout:write(("churn: %d of the compared blocks survived eviction\n"):format(still))
  end
  local t = vim.uv.hrtime()
  for _, w in ipairs(wanted) do
    for b = w.lo, w.hi do lay:build_block(b) end
  end
  ensure_ms = (vim.uv.hrtime() - t) / 1e6
end

---------------------------------------------------------------------------
-- Normalising a row for comparison
---------------------------------------------------------------------------
local function ids(list, key)
  if not list then return nil end
  local out = {}
  for i, x in ipairs(list) do out[i] = key(x) end
  return out
end

local function norm_runs(runs)
  return ids(runs, function(r) return { r[1], r[2], r[3], r[4] or 100 } end)
end

local function row_dump(lay_, row0)
  local info = lay_.rows[row0]
  local d = { row = row0, text = lay_.lines[row0 + 1] or "" }
  if not info then d.kind = false; return d end
  d.kind, d.bs, d.be = info.kind, info.bs, info.be
  d.runs = norm_runs(info.runs)
  d.site_runs = norm_runs(info.site_runs)
  d.segs = ids(info.segs, function(s) return { s[1], s[2], s[3] } end)
  -- an enzyme segment's payload is the site group; compare what it identifies
  d.esegs = ids(info.esegs, function(s) return { s[1], s[2], s[3].label, s[3].cut } end)
  d.items = ids(info.items, function(it) return it.f.id end)
  d.hidden = info.hidden
  d.strand = info.strand
  -- what the decoration provider would draw, which is where selection and
  -- search hits appear
  d.hl = norm_runs(lay_:row_highlights(row0))
  return d
end

local windows = {}
for _, w in ipairs(wanted) do
  local lo_row = windowed and lay.offsets[w.lo + 1] or nil
  local hi_row
  if windowed then
    hi_row = lay.offsets[w.hi + 1] + lay:plan_of(w.hi).n - 1
  else
    -- the eager layout keeps no offsets; find the window's rows by their block
    for r = 0, #lay.lines - 1 do
      local b = lay.rows[r] and lay.rows[r].b
      if b then
        if b >= w.lo and lo_row == nil then lo_row = r end
        if b <= w.hi then hi_row = r end
        if b > w.hi then break end
      end
    end
  end
  local rows = {}
  for r = lo_row, hi_row do rows[#rows + 1] = row_dump(lay, r) end
  -- position queries for this window's bases, including the ones the windowed
  -- layout has to answer from a plan rather than from text
  local pos = {}
  local bs, be = w.lo * W + 1, math.min((w.hi + 1) * W, lay.L)
  local stepbp = math.max(1, math.floor((be - bs + 1) / 200))
  for bp = bs, be, stepbp do
    local r, c = lay:pos_of(bp)
    pos[#pos + 1] = { bp, r, c, lay:bp_at(r, c) }
  end
  -- the feature drawn at each column of the window's first track row
  local feats = {}
  for _, d in ipairs(rows) do
    if d.kind == "track" then
      for c = 0, lay.prefix_w + lay.content_w do
        local f = lay:feature_at(d.row, c)
        if f then feats[#feats + 1] = { d.row, c, f } end
      end
    end
  end
  windows[#windows + 1] = { name = w.name, lo = w.lo, hi = w.hi,
    lo_row = lo_row, hi_row = hi_row, rows = rows, pos = pos, feats = feats }
end

collectgarbage("collect")
local dump = {
  file = vim.fn.fnamemodify(path, ":t"), mode = windowed and "window" or "eager",
  name = rec.name, len = lay.L, features = #rec.features, items = #lay.items,
  topology = rec.topology, args = A,
  parse_ms = parse_ms, build_ms = build_ms, ensure_ms = ensure_ms, churn_ms = churn_ms,
  nlines = #lay.lines, header_rows = lay.header_rows, nblocks = nblocks,
  W = lay.W, G = lay.G, prefix_w = lay.prefix_w, content_w = lay.content_w, ntracks = lay.ntracks,
  title = lay.lines[1], heap_mb = collectgarbage("count") / 1024,
  rss_mb = vim.uv.resident_set_memory() / 1048576,
  global = global, built_after_positions = built_after_positions,
  windows = windows,
}
local fh = assert(io.open(out_path, "wb"))
fh:write(vim.json.encode(dump))
fh:close()
io.stdout:write(("%s %s: %d bp, %d features, %d lines, build %.0f ms, heap %.0f MB, rss %.0f MB\n"):format(
  dump.file, dump.mode, dump.len, dump.features, dump.nlines, dump.build_ms, dump.heap_mb, dump.rss_mb))
