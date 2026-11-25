-- Opening a windowed record learns how long the buffer is by counting each
-- block's rows without planning it: no plan table, no track-assignment
-- tables, nothing kept. Everything downstream of the open -- where a block
-- starts, which base pair a row is, where a jump lands -- is read off those
-- counts, so a block that counts as one number and plans as another would
-- misplace every row after it. tests/spec/window_spec.lua proves the built
-- text matches the eager layout row for row; this proves the count a block is
-- opened with is the count it is later built with.
local H = require("helpers")

require("plasmid").setup({})
local io_ = require("plasmid.io")
local config = require("plasmid.config")
local linear = require("plasmid.view.linear")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function opts(G)
  local o = vim.deepcopy(config.options)
  o.linear.group = G
  o.linear.window = { enabled = true, margin = 0, cache_blocks = 1e9 }
  return o
end

--- Count and plan every block of a layout and insist on the same number.
--- Returns the blocks checked and the rows they came to.
local function both_ways(lay, tag)
  local rows = 0
  for b = 0, lay.nblocks - 1 do
    local n = lay._count(b)
    local p = lay._plan(b)
    if n ~= p.n then
      error(("%s: block %d counts as %d rows and plans as %d"):format(tag, b, n, p.n))
    end
    rows = rows + n
  end
  -- and the buffer the open sized is the one those blocks add up to
  H.eq(lay.header_rows + rows, lay.offsets[lay.nblocks] + lay:plan_of(lay.nblocks - 1).n,
    tag .. ": the blocks do not add up to the buffer's length")
  return lay.nblocks, rows
end

local function fixtures()
  local out = {}
  for _, pat in ipairs({ "**/*.gb", "**/*.gbk", "**/*.fa", "**/*.dna" }) do
    for _, p in ipairs(H.fixtures(pat)) do out[#out + 1] = p end
  end
  table.sort(out)
  return out
end

case("every_block_counts_as_it_plans", function()
  local files, blocks, rows = 0, 0, 0
  for _, path in ipairs(fixtures()) do
    local res = io_.read(path)
    if res then
      files = files + 1
      for ri, rec in ipairs(res.records) do
        for _, G in ipairs({ 10, 0 }) do
          local tag = ("%s r%d G=%d"):format(vim.fn.fnamemodify(path, ":t"), ri, G)
          local nb, n = both_ways(linear.build(rec, opts(G), 100, {}), tag)
          blocks, rows = blocks + nb, rows + n
        end
      end
    end
  end
  H.ok(files > 40, "expected the fixture corpus, got " .. files .. " files")
  H.log(("    [rowcount] %d files, %d blocks, %d rows counted and planned alike")
    :format(files, blocks, rows))
end)

-- Enzyme label lanes and primer lanes are assigned per line and are rows of
-- the block like any other, so they are the part of a count most likely to
-- drift from the plan beside it.
case("overlay_rows_count_as_they_plan", function()
  local checked, blocks = 0, 0
  for _, path in ipairs(fixtures()) do
    local res = io_.read(path)
    for _, rec in ipairs(res and res.records or {}) do
      if not rec:is_protein() and rec:len() > 0 and rec:len() < 60000 then
        local extra = {}
        local entries = require("plasmid.enzymes").compute(rec, { mode = "unique" })
        if #entries > 0 then extra.enzymes = entries end
        local primers = require("plasmid.primers").list(rec)
        if #primers > 0 then
          extra.primers = primers
          extra.primer_cells = function(p) return require("plasmid.primers").cells(rec, p) end
        end
        if extra.enzymes or extra.primers then
          checked = checked + 1
          local tag = "overlays " .. vim.fn.fnamemodify(path, ":t")
          local nb = both_ways(linear.build(rec, opts(10), 100, extra), tag)
          blocks = blocks + nb
          nb = both_ways(linear.build(rec, opts(0), 86, extra), tag .. " G=0")
          blocks = blocks + nb
        end
      end
    end
  end
  H.ok(checked >= 10, "expected records with enzymes or primers, got " .. checked)
  H.log(("    [rowcount] %d records with enzyme and primer overlays, %d blocks alike")
    :format(checked, blocks))
end)

-- A feature displaced off the last available track adds a summary row to its
-- block and the note on the title line, and the counting pass has to see both
-- without keeping the list of what was displaced.
case("displaced_features_are_counted_and_announced", function()
  local rec = assert(io_.read(H.fixture("ncbi/NC_005816.gb"))).records[1]
  local seen = {}
  for _, cap in ipairs({ 1, 2, "all" }) do
    local o = opts(10)
    o.linear.max_tracks = cap
    local lay = linear.build(rec, o, 100, {})
    both_ways(lay, "max_tracks=" .. tostring(cap))
    local over = 0
    for b = 0, lay.nblocks - 1 do
      if lay._plan(b).overflow then over = over + 1 end
    end
    seen[#seen + 1] = ("%s:%d"):format(cap, over)
    -- the title row carries the note when, and only when, a block overflowed
    H.eq(over > 0, lay.lines[1]:find("rows per line", 1, true) ~= nil,
      "the title note disagrees with the blocks at max_tracks=" .. tostring(cap))
  end
  H.log(("    [rowcount] blocks with a displaced feature by max_tracks: %s")
    :format(table.concat(seen, ", ")))
end)

-- The interval sets the track assignment decides against are reused between
-- blocks, so a block must leave nothing behind for the next one: counting the
-- blocks in a different order, or twice, has to give the same answers.
case("counting_leaves_nothing_behind", function()
  local rec = assert(io_.read(H.fixture("ncbi/NC_001133.gb"))).records[1]
  local lay = linear.build(rec, opts(10), 100, {})
  local first = {}
  for b = 0, lay.nblocks - 1 do first[b] = lay._count(b) end
  for b = lay.nblocks - 1, 0, -1 do
    H.eq(first[b], lay._count(b), "block " .. b .. " counted differently backwards")
  end
  for _ = 1, 3 do
    for b = 0, lay.nblocks - 1, 7 do
      H.eq(first[b], lay._count(b), "block " .. b .. " counted differently on a repeat")
    end
  end
  -- and a plan taken in between does not disturb the next count either
  for b = 0, lay.nblocks - 1, 13 do
    H.eq(first[b], lay._plan(b).n, "block " .. b .. " planned differently after counting")
    H.eq(first[b], lay._count(b), "block " .. b .. " counted differently after planning")
  end
  H.log(("    [rowcount] %d blocks counted forwards, backwards and repeatedly alike")
    :format(lay.nblocks))
end)

return T
