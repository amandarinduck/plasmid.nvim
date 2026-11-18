-- ResidueStore: the only way the rest of the code touches residues.
--
-- v1 backend is a plain immutable Lua string. Reads are O(1) slices; each
-- splice is an O(n) copy, which is comfortable for plasmid-sized records
-- (<= ~100 kb) and acceptable for read-only viewing up to a few Mb.
-- A piece table can replace this module later without changing callers.
--
-- Coordinates are 1-based inclusive (pos1/start1/stop1).
local M = {}
M.__index = M

function M.new(s)
  return setmetatable({ s = s or "" }, M)
end

function M:len()
  return #self.s
end

--- Residues start1..stop1 (inclusive, clamped).
function M:sub(start1, stop1)
  return self.s:sub(start1, stop1)
end

--- Residue at pos1 as a byte (fast path for per-base loops).
function M:byte(pos1)
  return self.s:byte(pos1)
end

--- Whole sequence as a string (cheap for this backend).
function M:tostring()
  return self.s
end

--- Replace the whole sequence.
function M:set(s)
  self.s = s
end

--- Replace start1..stop1 with `ins`. stop1 = start1 - 1 means pure insertion.
function M:splice(start1, stop1, ins)
  self.s = self.s:sub(1, start1 - 1) .. (ins or "") .. self.s:sub(stop1 + 1)
end

return M
