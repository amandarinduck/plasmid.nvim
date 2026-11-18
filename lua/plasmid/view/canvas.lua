-- Character canvas with a braille sub-cell layer (2x4 dots per cell) and a
-- text overlay. Each cell remembers a highlight group and an optional owner id
-- (feature id) so the rendered map can be queried under the cursor.
local bit = require("bit")

local Canvas = {}
Canvas.__index = Canvas

-- dot (dx, dy) within a cell -> braille bit
local BITS = {
  [0] = { [0] = 0x01, 0x02, 0x04, 0x40 },
  [1] = { [0] = 0x08, 0x10, 0x20, 0x80 },
}

local function braille(bits)
  -- U+2800 + bits, UTF-8 encoded (always 3 bytes)
  return string.char(0xE2, 0xA0 + bit.rshift(bits, 6), 0x80 + bit.band(bits, 0x3F))
end

function Canvas.new(w, h)
  return setmetatable({
    w = w, h = h,
    dw = 2 * w, dh = 4 * h,
    bits = {}, dot_hl = {}, dot_owner = {},
    text = {}, text_hl = {}, text_owner = {},
  }, Canvas)
end

--- Set a dot in dot coordinates (floats are rounded).
function Canvas:dot(x, y, group, owner)
  x, y = math.floor(x + 0.5), math.floor(y + 0.5)
  if x < 0 or y < 0 or x >= self.dw or y >= self.dh then return end
  local c, r = math.floor(x / 2), math.floor(y / 4)
  local i = r * self.w + c
  self.bits[i] = bit.bor(self.bits[i] or 0, BITS[x % 2][y % 4])
  if group then self.dot_hl[i] = group end
  if owner then self.dot_owner[i] = owner end
end

--- Straight line in dot space.
function Canvas:line(x0, y0, x1, y1, group, owner, dotted)
  local n = math.max(1, math.ceil(math.max(math.abs(x1 - x0), math.abs(y1 - y0))))
  for k = 0, n do
    if not dotted or k % 2 == 0 then
      local t = k / n
      self:dot(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t, group, owner)
    end
  end
end

--- Write text starting at cell (col, row); one codepoint per cell.
function Canvas:put(col, row, s, group, owner)
  if row < 0 or row >= self.h then return end
  local c = col
  for ch in s:gmatch("[%z\1-\127\194-\244][\128-\191]*") do
    if c >= 0 and c < self.w then
      local i = row * self.w + c
      self.text[i], self.text_hl[i], self.text_owner[i] = ch, group, owner
    end
    c = c + 1
  end
end

--- True if any text occupies cells col..col+len-1 on row.
function Canvas:text_free(col, row, len)
  for c = col, col + len - 1 do
    if self.text[row * self.w + c] then return false end
  end
  return true
end

--- Render rows. Returns lines, per-row highlight runs (byte ranges),
--- per-row cell->byte offsets, and per-row cell->owner maps.
function Canvas:render()
  local lines, runs, cell_bytes, owners = {}, {}, {}, {}
  for r = 0, self.h - 1 do
    local parts, rr, cb, ow = {}, {}, {}, {}
    local byte = 0
    local cur_hl, cur_b0
    local last_nonspace = 0
    for c = 0, self.w - 1 do
      local i = r * self.w + c
      local ch, group, owner
      if self.text[i] then
        ch, group, owner = self.text[i], self.text_hl[i], self.text_owner[i]
      elseif self.bits[i] then
        ch, group, owner = braille(self.bits[i]), self.dot_hl[i], self.dot_owner[i]
      else
        ch = " "
      end
      if group ~= cur_hl then
        if cur_hl then rr[#rr + 1] = { cur_b0, byte, cur_hl } end
        cur_hl, cur_b0 = group, byte
      end
      cb[c] = byte
      ow[c] = owner
      parts[#parts + 1] = ch
      byte = byte + #ch
      if ch ~= " " then last_nonspace = byte end
    end
    if cur_hl then rr[#rr + 1] = { cur_b0, byte, cur_hl } end
    cb[self.w] = byte
    lines[r + 1] = table.concat(parts):sub(1, last_nonspace)
    runs[r + 1], cell_bytes[r + 1], owners[r + 1] = rr, cb, ow
  end
  return lines, runs, cell_bytes, owners
end

return Canvas
