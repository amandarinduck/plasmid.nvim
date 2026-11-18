-- Transactions on a Record, with undo/redo. Pure Lua.
--
-- Every change is an op table; apply() performs it and returns the inverse.
--   { op = "add",      feature = Feature, index = i }   -- insert (keeps feature.id if set)
--   { op = "remove",   id = n }
--   { op = "replace",  id = n, feature = Feature }       -- new feature takes the same id
--   { op = "topology", value = "linear"|"circular" }
--   { op = "batch",    ops = { op, ... }, label = "..." }
-- Sequence edits (features are remapped, see model/edit.lua):
--   { op = "splice",   start = s, stop = e, ins = "ACGT" }  -- replace s..e (e = s-1: insert)
--   { op = "revcomp",  start = s, stop = e }                -- reverse-complement s..e
--   { op = "rotate",   origin = o }                         -- base o becomes base 1
-- Their inverses carry exact feature snapshots (`restore`) and the record's
-- dirty flags, so undo gives back the very same objects and file output.
-- After a forward sequence op, op.changes lists { old = Feature|nil,
-- new = Feature|nil, old_res, new_res } for features whose bases changed.
local M = {}

local edit = require("plasmid.model.edit")
local location = require("plasmid.model.location")
local alphabet = require("plasmid.model.alphabet")

local next_change = 0

local function index_of(rec, id)
  for i, f in ipairs(rec.features) do
    if f.id == id then return i end
  end
end

--- Index that keeps features ordered by start position (GenBank/SnapGene
--- convention); `source`-type features stay first.
function M.sorted_index(rec, f)
  local lo = f:bounds() or 0
  for i, g in ipairs(rec.features) do
    local glo = g:bounds() or 0
    if g.type ~= "source" and glo > lo then return i end
  end
  return #rec.features + 1
end

local function mark(rec, f)
  rec.revision = rec.revision + 1
  rec.raw = rec.raw or {}
  rec.raw.features_dirty = true
  if f and f.origin == "primer" then rec.raw.primers_dirty = true end
end

---------------------------------------------------------------------------
-- Sequence edits
---------------------------------------------------------------------------

local FLAGS = { "seq_dirty", "locus_dirty", "features_dirty", "primers_dirty" }

local function get_flags(rec)
  local t = {}
  for _, k in ipairs(FLAGS) do t[k] = rec.raw[k] or false end
  return t
end

local function set_flags(rec, t)
  for _, k in ipairs(FLAGS) do rec.raw[k] = t[k] or nil end
end

-- New Feature object at `node`, otherwise identical to f. Its verbatim source
-- lines are kept as raw_prev (qualifiers can still be written verbatim).
local function derive(f, node)
  local Feature = require("plasmid.model.record").Feature
  local quals = {}
  for i, q in ipairs(f.qualifiers) do quals[i] = { key = q.key, value = q.value, quoted = q.quoted } end
  local g = Feature.new({
    type = f.type, location_text = location.format(node), qualifiers = quals,
    color = f.color, no_arrow = f.no_arrow, origin = f.origin, recent_id = f.recent_id,
    primer_seq = f.primer_seq,
  })
  g.id, g.pinned = f.id, f.pinned
  g.raw_prev = f.raw_lines or f.raw_prev
  g.quals_verbatim = g.raw_prev ~= nil and (f.raw_lines ~= nil or f.quals_verbatim) or nil
  return g
end

-- Put features back exactly: entries { id, obj = Feature|false, index }.
-- Returns the entries that undo this restore.
local function restore_features(rec, entries)
  local undo = {}
  for _, r in ipairs(entries) do
    local i = index_of(rec, r.id)
    undo[#undo + 1] = { id = r.id, obj = i and rec.features[i] or false, index = i or 0 }
  end
  for _, r in ipairs(entries) do
    if not r.obj then
      local i = index_of(rec, r.id)
      if i then table.remove(rec.features, i) end
    end
  end
  local inserts = {}
  for _, r in ipairs(entries) do
    if r.obj then
      local i = index_of(rec, r.id)
      if i then rec.features[i] = r.obj else inserts[#inserts + 1] = r end
    end
  end
  table.sort(inserts, function(a, b) return a.index < b.index end)
  for _, r in ipairs(inserts) do table.insert(rec.features, math.min(r.index, #rec.features + 1), r.obj) end
  return undo
end

-- Replace features by remapped versions. `remap(f)` returns
-- node (moved), false (destroyed) or nil (unchanged). Returns restore entries
-- and the list of { old, new } pairs.
local function remap_features(rec, remap)
  local entries, pairs_, keep = {}, {}, {}
  for i, f in ipairs(rec.features) do
    local node = nil
    if f.loc then node = remap(f) end -- (not `f.loc and remap(f) or nil`: false means destroyed)
    if node == false then
      entries[#entries + 1] = { id = f.id, obj = f, index = i }
      pairs_[#pairs_ + 1] = { old = f }
    elseif node and location.format(node) ~= f.location_text then
      local g = derive(f, node)
      entries[#entries + 1] = { id = f.id, obj = f, index = i }
      pairs_[#pairs_ + 1] = { old = f, new = g }
      keep[#keep + 1] = g
    else
      keep[#keep + 1] = f
    end
  end
  rec.features = keep
  return entries, pairs_
end

local function overlaps(f, s, e)
  for _, p in ipairs(f.parts) do
    if not p.ref and p.stop >= p.start and p.start <= e and p.stop >= s then return true end
  end
  return false
end

local function finish_seq_op(rec, op, flags_before, entries, pairs_, inverse)
  rec.revision = rec.revision + 1
  -- counts residue changes only (undo/redo included): :DnaEdit checks it
  rec.seq_revision = (rec.seq_revision or 0) + 1
  if op.flags then
    set_flags(rec, op.flags)
  else
    rec.raw.seq_dirty = true
    rec.raw.locus_dirty = true
    if #entries > 0 then rec.raw.features_dirty = true end
    for _, p in ipairs(pairs_ or {}) do
      if (p.old and p.old.origin == "primer") then rec.raw.primers_dirty = true end
    end
  end
  inverse.restore = entries
  inverse.flags = flags_before
  inverse.label = op.label
  return inverse
end

local function apply_seq(rec, op)
  local k = op.op
  local flags_before = get_flags(rec)
  local L = rec:len()
  if k == "splice" then
    local s, e, ins = op.start, op.stop, op.ins or ""
    assert(s >= 1 and s <= L + 1 and e >= s - 1 and e <= L, ("splice %d..%d outside 1..%d"):format(s, e, L))
    local removed = rec.seq:sub(s, e)
    local inverse = { op = "splice", start = s, stop = s + #ins - 1, ins = removed }
    if op.restore then
      rec.seq:splice(s, e, ins)
      return finish_seq_op(rec, op, flags_before, restore_features(rec, op.restore), nil, inverse)
    end
    local n, kk, circ = e - s + 1, #ins, rec:is_circular()
    local substitution = n == kk and n > 0
    -- residues of features near the edit, to tell real changes from shifts
    local before = {}
    local lo, hi = math.max(1, s - 1), math.min(L, math.max(e, s))
    if hi >= lo then
      for _, f in ipairs(rec.features) do
        if f.loc and overlaps(f, lo, hi) then before[f] = rec:extract(f) end
      end
    end
    local entries, pairs_ = remap_features(rec, function(f)
      if substitution then return nil end
      local node = edit.splice(f.loc, s, e, kk, L, circ)
      if not node then return false end
      return node
    end)
    rec.seq:splice(s, e, ins)
    local new_of = {}
    for _, p in ipairs(pairs_) do if p.new then new_of[p.new] = p.old end end
    local changes = {}
    for _, f in ipairs(rec.features) do
      local old = new_of[f] or f
      if before[old] then
        local new_res = rec:extract(f)
        if new_res ~= before[old] then
          changes[#changes + 1] = { old = old, new = f, old_res = before[old], new_res = new_res }
        end
      end
    end
    for _, p in ipairs(pairs_) do
      if not p.new then changes[#changes + 1] = { old = p.old, old_res = before[p.old] or "" } end
    end
    op.changes = changes
    op.moved = #pairs_
    return finish_seq_op(rec, op, flags_before, entries, pairs_, inverse)
  elseif k == "revcomp" then
    local s, e = op.start, op.stop
    assert(s >= 1 and e <= L and s <= e, ("revcomp %d..%d outside 1..%d"):format(s, e, L))
    local inverse = { op = "revcomp", start = s, stop = e }
    local seg = rec.seq:sub(s, e)
    if op.restore then
      rec.seq:splice(s, e, alphabet.revcomp(seg))
      return finish_seq_op(rec, op, flags_before, restore_features(rec, op.restore), nil, inverse)
    end
    local before = {}
    for _, f in ipairs(rec.features) do
      if f.loc and overlaps(f, s, e) and not edit.inside(f.loc, s, e) then before[f] = rec:extract(f) end
    end
    local entries, pairs_ = remap_features(rec, function(f)
      if edit.inside(f.loc, s, e) then return edit.flip(f.loc, s, e) end
      return nil
    end)
    rec.seq:splice(s, e, alphabet.revcomp(seg))
    local changes = {}
    for _, f in ipairs(rec.features) do
      if before[f] then
        changes[#changes + 1] = { old = f, new = f, old_res = before[f], new_res = rec:extract(f), inverted = true }
      end
    end
    op.changes, op.moved = changes, #pairs_
    return finish_seq_op(rec, op, flags_before, entries, pairs_, inverse)
  elseif k == "rotate" then
    local o = op.origin
    assert(rec:is_circular(), "only circular records can be rotated")
    assert(o >= 1 and o <= L, ("origin %d outside 1..%d"):format(o, L))
    local inverse = { op = "rotate", origin = (L - o + 1) % L + 1 }
    local str = rec.seq:tostring()
    if op.restore then
      rec.seq:set(str:sub(o) .. str:sub(1, o - 1))
      return finish_seq_op(rec, op, flags_before, restore_features(rec, op.restore), nil, inverse)
    end
    local entries, pairs_ = remap_features(rec, function(f) return edit.rotate(f.loc, o, L) end)
    rec.seq:set(str:sub(o) .. str:sub(1, o - 1))
    op.changes, op.moved = {}, #pairs_
    return finish_seq_op(rec, op, flags_before, entries, pairs_, inverse)
  end
end

M.SEQ_OPS = { splice = true, revcomp = true, rotate = true }

function M.apply(rec, op)
  local k = op.op
  if M.SEQ_OPS[k] then return apply_seq(rec, op) end
  if k == "add" then
    local f = op.feature
    if not f.id then
      rec._next_id = (rec._next_id or 0) + 1
      f.id = rec._next_id
    end
    if rec.alphabet == "protein" then f.no_arrow = true end
    local i = op.index or M.sorted_index(rec, f)
    table.insert(rec.features, i, f)
    mark(rec, f)
    return { op = "remove", id = f.id }
  elseif k == "remove" then
    local i = assert(index_of(rec, op.id), "no feature with id " .. tostring(op.id))
    local f = table.remove(rec.features, i)
    mark(rec, f)
    return { op = "add", feature = f, index = i }
  elseif k == "replace" then
    local i = assert(index_of(rec, op.id), "no feature with id " .. tostring(op.id))
    local old = rec.features[i]
    op.feature.id = old.id
    op.feature.origin = op.feature.origin or old.origin
    rec.features[i] = op.feature
    mark(rec, old)
    return { op = "replace", id = old.id, feature = old }
  elseif k == "topology" then
    local old = rec.topology
    rec:set_topology(op.value)
    return { op = "topology", value = old }
  elseif k == "batch" then
    local inv = {}
    for i, sub in ipairs(op.ops) do inv[#op.ops - i + 1] = M.apply(rec, sub) end
    return { op = "batch", ops = inv, label = op.label }
  end
  error("unknown op " .. tostring(k))
end

local function history(rec)
  rec.history = rec.history or { undo = {}, redo = {} }
  return rec.history
end

--- Apply an op and record it for undo. Returns the change id.
function M.commit(rec, op, label)
  local inv = M.apply(rec, op)
  next_change = next_change + 1
  local h = history(rec)
  h.undo[#h.undo + 1] = { forward = op, inverse = inv, label = label or op.op, id = next_change }
  h.redo = {}
  return next_change
end

--- Record an already-applied change (forward op + its inverse) for undo.
function M.record(rec, forward, inverse, label)
  next_change = next_change + 1
  local h = history(rec)
  h.undo[#h.undo + 1] = { forward = forward, inverse = inverse, label = label or forward.op, id = next_change }
  h.redo = {}
  return next_change
end

function M.undo(rec)
  local h = history(rec)
  local c = table.remove(h.undo)
  if not c then return nil end
  c.forward = M.apply(rec, c.inverse) -- apply() returns the op that re-does it
  h.redo[#h.redo + 1] = c
  return c
end

function M.redo(rec)
  local h = history(rec)
  local c = table.remove(h.redo)
  if not c then return nil end
  c.inverse = M.apply(rec, c.forward)
  h.undo[#h.undo + 1] = c
  return c
end

--- Identifier of the current state (top of the undo stack), for "is this
--- what is on disk?" checks.
function M.state_id(rec)
  local h = history(rec)
  local top = h.undo[#h.undo]
  return top and top.id or 0
end

return M
