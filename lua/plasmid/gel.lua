-- A text agarose gel: lanes of bands placed on a log(size) axis next to a ladder.
local M = {}

-- Band sizes (bp) of generic ladders.
M.LADDERS = {
  ["1kb"] = { 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000 },
  ["100bp"] = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500 },
}
M.LADDER_NAMES = { ["1kb"] = "1 kb ladder", ["100bp"] = "100 bp ladder" }

--- Pick a ladder for a set of sizes.
function M.auto_ladder(sizes)
  local lo, hi = math.huge, 0
  for _, s in ipairs(sizes) do lo = math.min(lo, s); hi = math.max(hi, s) end
  if hi <= 1500 and lo < 500 then return "100bp" end
  return "1kb"
end

local BAND = "━"

--- Render lanes { {label=, sizes={...}}, ... } with a ladder.
--- opts { ladder = "auto"|"1kb"|"100bp", height = rows }.
--- Returns lines and highlights { row0, byte0, byte1, group } (row0 relative
--- to the first returned line), plus a per-lane placement table for tests:
--- placed[lane] = { {size, row}, ... }, off[lane] = { sizes outside the range }.
function M.render(lanes, opts)
  opts = opts or {}
  local all = {}
  for _, ln in ipairs(lanes) do for _, s in ipairs(ln.sizes) do all[#all + 1] = s end end
  local key = opts.ladder
  if not key or key == "auto" or not M.LADDERS[key] then key = M.auto_ladder(all) end
  local ladder = M.LADDERS[key]
  local height = math.max(8, tonumber(opts.height) or 24)
  local top, bottom = ladder[#ladder] * 1.25, ladder[1] * 0.8
  local lt, lb = math.log(top), math.log(bottom)
  local function row_of(size)
    if size > top or size < bottom then return nil end
    return math.floor((lt - math.log(size)) / (lt - lb) * (height - 1) + 0.5)
  end

  local band_w = 7
  local size_w = #tostring(ladder[#ladder])
  -- columns: size labels, ladder lane, then sample lanes (band + size text)
  local cols = {}
  local x = size_w + 3
  cols[0] = x
  x = x + math.max(band_w, #M.LADDER_NAMES[key]) + 3
  local sample_w = {}
  for i, ln in ipairs(lanes) do
    local maxlab = 0
    for _, s in ipairs(ln.sizes) do maxlab = math.max(maxlab, #tostring(s) + 3) end
    sample_w[i] = math.max(band_w + 1 + maxlab, vim.fn.strdisplaywidth(ln.label or ""))
    cols[i] = x
    x = x + sample_w[i] + 2
  end

  local grid = {}
  for r = 0, height - 1 do grid[r] = {} end
  local hls = {}
  local placed, off = {}, {}

  -- rows: put text at a display column (ASCII and "━" only; "━" is 3 bytes)
  local function put(r, col, text, group)
    grid[r][#grid[r] + 1] = { col = col, text = text, group = group }
  end
  for _, s in ipairs(ladder) do
    local r = row_of(s)
    put(r, 1, ("%" .. size_w .. "d"):format(s), "PlasmidGelLabel")
    put(r, cols[0], BAND:rep(band_w), "PlasmidGelLadder")
  end
  for i, ln in ipairs(lanes) do
    placed[i], off[i] = {}, {}
    local by_row = {}
    for _, s in ipairs(ln.sizes) do
      local r = row_of(s)
      if r then
        by_row[r] = by_row[r] or {}
        table.insert(by_row[r], s)
        placed[i][#placed[i] + 1] = { s, r }
      else
        off[i][#off[i] + 1] = s
      end
    end
    for r, list in pairs(by_row) do
      table.sort(list, function(a, b) return a > b end)
      local strong = #list > 1
      put(r, cols[i], BAND:rep(band_w), strong and "PlasmidGelBandStrong" or "PlasmidGelBand")
      local txt = tostring(list[1]) .. (strong and (" ×" .. #list) or "")
      if strong and list[#list] ~= list[1] then txt = list[1] .. "–" .. list[#list] .. " ×" .. #list end
      put(r, cols[i] + band_w + 1, txt, "PlasmidGelLabel")
    end
  end

  local lines, out_hl = {}, {}
  -- header: lane labels and wells
  local function compose(items)
    table.sort(items, function(a, b) return a.col < b.col end)
    local parts, col, byte = {}, 1, 0
    local h = {}
    for _, it in ipairs(items) do
      if it.col > col then
        parts[#parts + 1] = (" "):rep(it.col - col)
        byte = byte + (it.col - col)
        col = it.col
      elseif it.col < col then
        goto skip -- overlapping text: drop
      end
      h[#h + 1] = { byte, byte + #it.text, it.group }
      parts[#parts + 1] = it.text
      byte = byte + #it.text
      col = col + vim.fn.strdisplaywidth(it.text)
      ::skip::
    end
    return table.concat(parts), h
  end
  local head = { { col = cols[0], text = M.LADDER_NAMES[key], group = "PlasmidGelLabel" } }
  for i, ln in ipairs(lanes) do head[#head + 1] = { col = cols[i], text = ln.label or ("lane " .. i), group = "PlasmidGelLabel" } end
  local wells = { { col = cols[0], text = ("▔"):rep(band_w), group = "PlasmidGelWell" } }
  for i = 1, #lanes do wells[#wells + 1] = { col = cols[i], text = ("▔"):rep(band_w), group = "PlasmidGelWell" } end
  for _, items in ipairs({ head, wells }) do
    local l, h = compose(items)
    lines[#lines + 1] = l
    for _, x in ipairs(h) do out_hl[#out_hl + 1] = { #lines - 1, x[1], x[2], x[3] } end
  end
  local first = #lines
  for r = 0, height - 1 do
    local l, h = compose(grid[r])
    lines[#lines + 1] = l
    for _, x in ipairs(h) do out_hl[#out_hl + 1] = { #lines - 1, x[1], x[2], x[3] } end
  end
  local notes = {}
  for i, list in ipairs(off) do
    if #list > 0 then
      table.sort(list, function(a, b) return a > b end)
      local big, small = {}, {}
      for _, s in ipairs(list) do if s > top then big[#big + 1] = s else small[#small + 1] = s end end
      local t = {}
      if #big > 0 then t[#t + 1] = "above the ladder: " .. table.concat(big, ", ") .. " bp" end
      if #small > 0 then t[#t + 1] = "ran off: " .. table.concat(small, ", ") .. " bp" end
      notes[#notes + 1] = ("%s: %s"):format(lanes[i].label or ("lane " .. i), table.concat(t, "; "))
    end
  end
  for _, ln in ipairs(lanes) do
    if ln.note then notes[#notes + 1] = ("%s: %s"):format(ln.label or "lane", ln.note) end
  end
  for _, n in ipairs(notes) do
    lines[#lines + 1] = " " .. n
    out_hl[#out_hl + 1] = { #lines - 1, 0, #lines[#lines], "PlasmidOverflow" }
  end
  return lines, out_hl, { placed = placed, off = off, first = first, ladder = key, row_of = row_of }
end

return M
