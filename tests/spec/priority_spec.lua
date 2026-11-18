local H = require("helpers")
local pio = require("plasmid.io")
local config = require("plasmid.config")
local vf = require("plasmid.view.features")
local linear = require("plasmid.view.linear")

local T = {}

local function build(path, cap)
  config.setup({ linear = { max_tracks = cap } })
  require("plasmid.view.highlights").setup()
  local rec = assert(pio.read(H.fixture(path))).records[1]
  return linear.build(rec, config.options, 120), rec
end

local function spans(it, bs, be)
  local out = {}
  for _, e in ipairs(it.extents) do
    local s, x = math.max(e[1], bs), math.min(e[2], be)
    if s <= x then out[#out + 1] = { s, x } end
  end
  return out
end

local function touch(a, b)
  for _, x in ipairs(spans(a, a._bs, a._be)) do
    for _, y in ipairs(spans(b, a._bs, a._be)) do
      if x[1] <= y[2] + 1 and y[1] <= x[2] + 1 then return true end
    end
  end
  return false
end

-- group rows by block
local function blocks(lay)
  local out = {}
  for row0, info in pairs(lay.rows) do
    if info.b and (info.kind == "track" or info.kind == "overflow") then
      out[info.b] = out[info.b] or { bs = info.bs, be = info.be, rows = {}, hidden = {} }
      if info.kind == "track" then table.insert(out[info.b].rows, info.items) end
      if info.kind == "overflow" then out[info.b].hidden = info.hidden_items end
    end
  end
  return out
end

function T.priority_scores_follow_config()
  local rec = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  config.setup({})
  local items = vf.visible(rec, config.options)
  local by = vf.by_priority(items, config.options)
  H.eq("CDS", by[1].f.type)
  for i = 2, #by do H.ok(by[i - 1].priority >= by[i].priority) end
  -- labelled beats unlabelled within a type
  local named = { f = { type = "misc_feature", id = 1, length = function() return 10 end }, label = "lacO" }
  local anon = { f = { type = "misc_feature", id = 2, length = function() return 10 end }, label = "misc_feature" }
  H.ok(vf.priority(named) > vf.priority(anon))
  -- user override
  config.setup({ priority = { types = { misc_feature = 500 } } })
  local fresh = { f = { type = "misc_feature", id = 3 }, label = "x" }
  H.eq(510, vf.priority(fresh, config.options))
  config.setup({})
end

function T.global_tracks_are_priority_justified()
  for _, path in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001133.gb", "ncbi/U01317.gb" }) do
    config.setup({})
    local rec = assert(pio.read(H.fixture(path))).records[1]
    local items = vf.visible(rec, config.options)
    vf.pack(items, rec:len(), 0, config.options)
    for _, it in ipairs(items) do
      -- every lower track holds an overlapping feature that was placed first
      for t = 1, it.track - 1 do
        local blocked = false
        for _, o in ipairs(items) do
          if o.track == t and vf.priority(o) >= vf.priority(it) then
            for _, a in ipairs(it.extents) do
              for _, b in ipairs(o.extents) do
                if a[1] <= b[2] and b[1] <= a[2] then blocked = true end
              end
            end
          end
        end
        H.ok(blocked, ("%s %s on track %d but track %d was free"):format(path, it.f.location_text, it.track, t))
      end
    end
  end
end

function T.hidden_features_only_yield_to_equal_or_higher_priority()
  local overflowed = 0
  for _, cap in ipairs({ 1, 2, 4 }) do
    for _, path in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001133.gb", "ncbi/U01317.gb", "snapgene/pFA-KanMX4.dna" }) do
      local lay = build(path, cap)
      for _, blk in pairs(blocks(lay)) do
        H.ok(#blk.rows <= cap, ("%s: %d rows > cap %d"):format(path, #blk.rows, cap))
        for _, h in ipairs(blk.hidden) do
          overflowed = overflowed + 1
          h._bs, h._be = blk.bs, blk.be
          -- every shown row must be blocked here by a feature ranked at least as high
          H.eq(cap, #blk.rows, path .. ": feature hidden while a row was unused")
          for _, row in ipairs(blk.rows) do
            local blocked = false
            for _, s in ipairs(row) do
              if s.rank < h.rank and touch(h, s) then blocked = true end
            end
            H.ok(blocked, ("%s cap %d: %s hidden but a row had room"):format(path, cap, h.f.location_text))
          end
        end
      end
    end
  end
  H.ok(overflowed > 20, "fixtures produced too little overflow to test: " .. overflowed)
end

function T.overflow_counts_are_exact()
  for _, path in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001133.gb" }) do
    local lay = build(path, 2)
    local rec = assert(pio.read(H.fixture(path))).records[1]
    local items = vf.visible(rec, config.options)
    for _, blk in pairs(blocks(lay)) do
      local expected = 0
      for _, it in ipairs(items) do
        if #spans(it, blk.bs, blk.be) > 0 then expected = expected + 1 end
      end
      local shown = 0
      for _, row in ipairs(blk.rows) do shown = shown + #row end
      H.eq(expected, shown + #blk.hidden, path .. " block " .. blk.bs)
    end
    -- the overflow line's text states the same number
    for row0, info in pairs(lay.rows) do
      if info.kind == "overflow" then
        H.eq(#info.hidden, tonumber(lay.lines[row0 + 1]:match("%+(%d+) more")))
      end
    end
  end
end

function T.expand_all_shows_every_feature()
  local lay, rec = build("ncbi/NC_001133.gb", "all")
  local seen = {}
  for _, info in pairs(lay.rows) do
    H.ok(info.kind ~= "overflow", "overflow row with cap=all")
    for _, it in ipairs(info.items or {}) do seen[it.f.id] = true end
  end
  for _, it in ipairs(vf.visible(rec, config.options)) do
    H.ok(seen[it.f.id], "never drawn: " .. it.f.location_text)
  end
  config.setup({})
end

return T
