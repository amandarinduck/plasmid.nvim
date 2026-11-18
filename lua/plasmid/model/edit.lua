-- Coordinate remapping of feature locations under sequence edits. Pure Lua.
--
-- Works on the location AST (model/location.lua) so a feature keeps its form
-- (join/order/complement, fuzzy marks) and only its numbers change. Every
-- function takes an AST and returns a new one (or nil when nothing is left);
-- inputs are never mutated.
--
-- Splice semantics, edit = replace start..stop (1-based, stop = start - 1 for a
-- pure insertion) with k new residues:
--   * deleted bases leave every feature; a part losing all its bases is dropped
--     and a feature losing all its parts is destroyed
--   * a pure insertion strictly inside a part extends it; at a part boundary it
--     shifts the part (the new bases belong to no feature). On a circular record
--     an insertion at the origin (before base 1 or after base L) extends a
--     feature that crosses the origin
--   * a replacement's new bases belong to the features that covered the bases
--     they replace: a part spanning the region is extended, otherwise the part
--     touching the region on its left (or else its right) takes them, otherwise
--     a part that lay wholly inside the region becomes the new bases
local M = {}

local location = require("plasmid.model.location")

local LEAF = { range = true, site = true, within = true, between = true }

function M.clone(node)
  local c = {}
  for k, v in pairs(node) do c[k] = v end
  if node.child then c.child = M.clone(node.child) end
  if node.parts then
    c.parts = {}
    for i, p in ipairs(node.parts) do c.parts[i] = M.clone(p) end
  end
  return c
end

-- Leaves in tree order (local ones only), with accessors.
local function leaves(node, out)
  out = out or {}
  if LEAF[node.kind] then
    if not node.ref then out[#out + 1] = node end
  elseif node.child then
    leaves(node.child, out)
  else
    for _, p in ipairs(node.parts) do leaves(p, out) end
  end
  return out
end
M.leaves = leaves

local function span(l)
  if l.kind == "site" then return l.pos, l.pos end
  if l.kind == "between" then return l.left, l.right end
  return l.start, l.stop
end

local function set_span(l, a, b)
  if l.kind == "site" then l.pos = a
  elseif l.kind == "between" then l.left, l.right = a, b
  else l.start, l.stop = a, b end
end

-- A site/within/between that must become a range to hold several bases.
local function as_range(l, a, b)
  if l.kind == "range" then l.start, l.stop = a, b; return end
  local mark = l.mark
  for k in pairs(l) do if k ~= "ref" then l[k] = nil end end
  l.kind, l.start, l.stop = "range", a, b
  if mark then l.start_mark, l.stop_mark = mark, mark end
end

-- Remove deleted leaves; collapse single-part join/order; nil if empty.
local function prune(node)
  if LEAF[node.kind] then return (not node._del) and node or nil end
  if node.child then
    local c = prune(node.child)
    if not c then return nil end
    node.child = c
    return node
  end
  local parts = {}
  for _, p in ipairs(node.parts) do
    local q = prune(p)
    if q then parts[#parts + 1] = q end
  end
  if #parts == 0 then return nil end
  if #parts == 1 then return parts[1] end
  node.parts = parts
  return node
end

local function clear_flags(node)
  for _, l in ipairs(leaves(node)) do l._del, l._touched, l._orig_a = nil, nil, nil end
end

--- Splice: replace s..e (e = s - 1: none) by k residues on a record of
--- length L (before the edit). Returns new AST or nil (destroyed), plus
--- info { touched = bool (bases removed from it), grown = bool }.
function M.splice(node, s, e, k, L, circular)
  local t = M.clone(node)
  local ls = leaves(t)
  local n = e - s + 1
  local info = { touched = false, grown = false }
  -- phase A: deletion
  if n > 0 then
    for _, l in ipairs(ls) do
      local a, b = span(l)
      if l.kind == "between" then
        local junction = (a == L and b == 1)
        local function dead(x) return x >= s and x <= e end
        if dead(a) or dead(b) then
          l._del, info.touched = true, true
        else
          local a2 = a > e and a - n or a
          set_span(l, a2, junction and 1 or a2 + 1)
        end
      elseif b < s then
        -- before the deletion
      elseif a > e then
        set_span(l, a - n, b - n)
      else
        l._touched, info.touched = true, true
        l._orig_a = a
        if a >= s and b <= e then
          l._del = true
        else
          local a2 = a < s and a or s
          local b2 = b > e and b - n or s - 1
          if l.kind == "site" then as_range(l, a2, b2) else set_span(l, a2, b2) end
        end
      end
    end
  end
  local L2 = L - math.max(n, 0)
  -- phase B: insertion at s
  if k > 0 then
    local pos = s
    local special -- leaf that takes the new bases instead of the standard rule
    if n > 0 then
      -- replacement: spanning parts extend by the standard rule; otherwise
      -- the touched part on the left, then on the right, then a wholly deleted one
      local spanning = false
      for _, l in ipairs(ls) do
        if not l._del and l.kind ~= "between" then
          local a, b = span(l)
          if a < pos and b >= pos then spanning = true end
        end
      end
      if not spanning then
        for _, l in ipairs(ls) do
          if l._touched and not l._del and select(2, span(l)) == pos - 1 then special = special or { l, "left" } end
        end
        for _, l in ipairs(ls) do
          if not special and l._touched and not l._del and span(l) == pos then special = { l, "right" } end
        end
        -- the leftmost wholly replaced part (order-independent choice)
        local best
        for _, l in ipairs(ls) do
          if l._del and l.kind ~= "between" and (not best or l._orig_a < best._orig_a) then best = l end
        end
        if not special and best then special = { best, "whole" } end
      end
    elseif circular and (pos == 1 or pos == L + 1) then
      -- insertion at the origin extends a feature crossing it
      local at_end, at_start
      for _, l in ipairs(ls) do
        if l.kind == "range" then
          local a, b = span(l)
          if b == L then at_end = l end
          if a == 1 then at_start = l end
        end
      end
      if at_end and at_start and at_end ~= at_start then
        special = pos == 1 and { at_start, "origin_start" } or { at_end, "origin_end" }
      end
    end
    for _, l in ipairs(ls) do
      local a, b = span(l)
      if special and l == special[1] then
        local how = special[2]
        if how == "whole" then
          l._del = nil
          as_range(l, pos, pos + k - 1)
        elseif how == "left" or how == "origin_end" then
          if l.kind ~= "range" then as_range(l, a, b + k) else set_span(l, a, b + k) end
        else -- right / origin_start: starts at pos, takes the new bases
          if l.kind ~= "range" then as_range(l, a, b + k) else set_span(l, a, b + k) end
        end
        info.grown = true
      elseif l._del then
        -- stays deleted
      elseif l.kind == "between" then
        local junction = (a == L2 and b == 1)
        local a2 = a >= pos and a + k or a
        set_span(l, a2, junction and 1 or a2 + 1)
      elseif l.kind == "site" then
        if a >= pos then set_span(l, a + k, a + k) end
      elseif a >= pos then
        set_span(l, a + k, b + k)
      elseif b >= pos then
        set_span(l, a, b + k)
        info.grown = true
      end
    end
  end
  local out = prune(t)
  if out then clear_flags(out) end
  return out, info
end

---------------------------------------------------------------------------
-- Reverse complement of a region
---------------------------------------------------------------------------

--- True if every local leaf lies within s..e.
function M.inside(node, s, e)
  local any = false
  for _, l in ipairs(leaves(node)) do
    local a, b = span(l)
    if l.kind == "between" then a, b = math.min(a, b), math.max(a, b) end
    if a < s or b > e then return false end
    any = true
  end
  return any
end

local SWAP_MARK = { ["<"] = ">", [">"] = "<" }

-- Mirror coordinates in s..e and reverse join order; strands unchanged.
local function mirror(node, s, e)
  local m = function(x) return s + e - x end
  local k = node.kind
  if node.ref then return node end
  if k == "range" then
    node.start, node.stop = m(node.stop), m(node.start)
    node.start_mark, node.stop_mark = SWAP_MARK[node.stop_mark], SWAP_MARK[node.start_mark]
  elseif k == "site" then
    node.pos = m(node.pos)
    node.mark = SWAP_MARK[node.mark] or node.mark
  elseif k == "within" then
    node.start, node.stop = m(node.stop), m(node.start)
  elseif k == "between" then
    node.left, node.right = m(node.right), m(node.left)
  elseif k == "complement" then
    node.child = mirror(node.child, s, e)
  else
    local r = {}
    for i = #node.parts, 1, -1 do r[#r + 1] = mirror(node.parts[i], s, e) end
    node.parts = r
  end
  return node
end

local function all_complement(parts)
  for _, p in ipairs(parts) do if p.kind ~= "complement" then return false end end
  return true
end
local function no_complement(node)
  if node.kind == "complement" then return false end
  if node.parts then for _, p in ipairs(node.parts) do if not no_complement(p) then return false end end end
  return true
end

-- Mirror and flip strand, keeping the conventional form:
-- complement(join(a,b)) <-> join(a',b') rather than nesting complements.
local function flip(node, s, e)
  local k = node.kind
  if k == "complement" then return mirror(node.child, s, e) end
  -- join/order list parts in biological order, which an inversion keeps:
  -- join(complement(a),complement(b)) -> join(a',b'), and mixed joins flip
  -- part by part. Only an all-forward join needs the complement(join(...))
  -- form, whose parts are listed in reverse.
  if (k == "join" or k == "order" or k == "bond") and all_complement(node.parts) then
    local r = {}
    for i, p in ipairs(node.parts) do r[i] = mirror(p.child, s, e) end
    node.parts = r
    return node
  end
  if (k == "join" or k == "order" or k == "bond") and not no_complement(node) then
    local r = {}
    for i, p in ipairs(node.parts) do r[i] = flip(p, s, e) end
    node.parts = r
    return node
  end
  return { kind = "complement", child = mirror(node, s, e) }
end

--- Location after reverse-complementing s..e, for a feature wholly inside it.
function M.flip(node, s, e)
  return flip(M.clone(node), s, e)
end

---------------------------------------------------------------------------
-- Rotate the origin of a circular record
---------------------------------------------------------------------------

--- Location after making base `o` the new base 1 (record length L).
function M.rotate(node, o, L)
  local t = M.clone(node)
  local function r(x) return ((x - o) % L) + 1 end
  local function walk(nd)
    local k = nd.kind
    if nd.ref then return nd end
    if k == "range" then
      local a, b = nd.start, nd.stop
      local junction_end, junction_start = (b == L), (a == 1)
      if a < o and b >= o then
        -- crosses the new origin: split
        -- the pieces keep the original part's origin-junction ends, so a
        -- neighbour that was its junction partner can still merge with it
        local left = { kind = "range", start = r(a), stop = L, start_mark = nd.start_mark,
          _was_start = junction_start }
        local right = { kind = "range", start = 1, stop = r(b), stop_mark = nd.stop_mark,
          _was_end = junction_end }
        return { kind = "join", parts = { left, right }, _split = true }
      end
      nd.start, nd.stop = r(a), r(b)
      nd._was_end, nd._was_start = junction_end, junction_start
      return nd
    elseif k == "site" then nd.pos = r(nd.pos)
    elseif k == "within" then
      if nd.start < o and nd.stop >= o then
        -- a.b cannot cross the origin in GenBank syntax: keep the span as ranges
        return { kind = "join", _split = true, parts = {
          { kind = "range", start = r(nd.start), stop = L },
          { kind = "range", start = 1, stop = r(nd.stop) },
        } }
      end
      nd.start, nd.stop = r(nd.start), r(nd.stop)
    elseif k == "between" then
      nd.left = r(nd.left)
      nd.right = nd.left == L and 1 or nd.left + 1
    elseif k == "complement" then nd.child = walk(nd.child)
    else
      local parts = {}
      for _, p in ipairs(nd.parts) do
        local q = walk(p)
        if q._split and (k == "join" or k == "order" or k == "bond") then
          for _, x in ipairs(q.parts) do parts[#parts + 1] = x end
        else
          parts[#parts + 1] = q
        end
      end
      -- merge a former origin junction (..L)(1..) that is now contiguous
      local merged = {}
      for _, p in ipairs(parts) do
        local prev = merged[#merged]
        if prev and prev.kind == "range" and p.kind == "range" and prev._was_end and p._was_start
          and prev.stop + 1 == p.start then
          prev.stop, prev.stop_mark = p.stop, p.stop_mark
          prev._was_end = p._was_end
        else
          merged[#merged + 1] = p
        end
      end
      nd.parts = merged
      if #merged == 1 then return merged[1] end
    end
    return nd
  end
  local out = walk(t)
  if out._split then out._split = nil end
  -- strip helper flags
  local function strip(nd)
    nd._was_end, nd._was_start, nd._jend, nd._jstart, nd._split = nil, nil, nil, nil, nil
    if nd.child then strip(nd.child) end
    if nd.parts then for _, p in ipairs(nd.parts) do strip(p) end end
  end
  strip(out)
  return out
end

function M.format(node) return location.format(node) end

return M
