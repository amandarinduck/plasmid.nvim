-- Cloning: double-stranded fragments, the joins between them, and the records
-- that come out. Pure Lua; no UI.
--
-- A fragment is everything present on either strand, written as the top strand
-- 5'->3' (`rec`, a linear record that carries the part's features). Each end
-- says how many bases there are single-stranded and which strand carries them:
--
--   left  = { len = n, strand =  1 }  first n bases are top-strand only
--                                     -> the top strand's 5' end protrudes:    5' overhang
--   left  = { len = n, strand = -1 }  first n bases are bottom-strand only
--                                     -> the bottom strand's 3' end protrudes: 3' overhang
--   right = { len = n, strand = -1 }  last n bases are bottom-strand only:     5' overhang
--   right = { len = n, strand =  1 }  last n bases are top-strand only:        3' overhang
--   len = 0                           blunt
--
-- Two ends anneal when they are the same length, lie on opposite strands and
-- read the same along the top strand (their single strands are then
-- complementary). A join writes the shared bases once:
--
--   product = A.seq .. B.seq:sub(n + 1)
--
-- Gibson/HiFi joins use the same arithmetic with n = the length of the
-- terminal homology, so one assembler serves all three methods.
local alphabet = require("plasmid.model.alphabet")
local E = require("plasmid.model.edit")
local T = require("plasmid.model.transaction")

local M = {}

local function wrap(p, L) return (p - 1) % L + 1 end

function M.blunt() return { len = 0, strand = 1 } end

--- Deep copy of a record (features included), via the GenBank writer/parser.
function M.copy(rec)
  local pio = require("plasmid.io")
  local text = pio.format({ rec }, "genbank")
  local out = assert(pio.parse(text, nil, "genbank"), "could not copy the record")
  local copy = out.records[1]
  copy.history = nil
  return copy
end

local function seqof(frag) return frag.rec.seq:tostring() end
M.sequence = seqof

--- Top-strand text of an end ("" when blunt).
function M.end_seq(frag, side)
  local e = frag[side]
  if e.len == 0 then return "" end
  local s = seqof(frag):upper()
  if side == "left" then return s:sub(1, e.len) end
  return s:sub(#s - e.len + 1)
end

--- "blunt", "5'" or "3'" for an end.
function M.end_kind(frag, side)
  local e = frag[side]
  if e.len == 0 then return "blunt" end
  local five = (side == "left") == (e.strand == 1)
  return five and "5'" or "3'"
end

--- Human-readable end: "blunt", "5′ AATT", "3′ TGCA". The bases are written
--- along the top strand, as the digest table does, so that two ends which
--- anneal read the same.
function M.end_text(frag, side)
  local kind = M.end_kind(frag, side)
  if kind == "blunt" then return "blunt" end
  return ("%s %s"):format(kind == "5'" and "5′" or "3′", M.end_seq(frag, side))
end

--- The protruding single strand itself, read 5'->3'.
function M.overhang_strand(frag, side)
  local top = M.end_seq(frag, side)
  if top == "" then return "" end
  return frag[side].strand == 1 and top or alphabet.revcomp(top)
end

--- Can fragment `a` be followed by fragment `b`? Returns a junction
--- { kind = "blunt"|"sticky", len, seq } or nil plus the reason.
function M.compatible(a, b)
  local ar, bl = a.right, b.left
  if ar.len ~= bl.len then
    return nil, ("%s / %s (lengths differ)"):format(M.end_text(a, "right"), M.end_text(b, "left"))
  end
  if ar.len == 0 then return { kind = "blunt", len = 0, seq = "" } end
  if ar.strand == bl.strand then
    return nil, ("%s / %s (same strand)"):format(M.end_text(a, "right"), M.end_text(b, "left"))
  end
  local x, y = M.end_seq(a, "right"), M.end_seq(b, "left")
  if x ~= y then
    return nil, ("%s / %s"):format(M.end_text(a, "right"), M.end_text(b, "left"))
  end
  return { kind = "sticky", len = ar.len, seq = x, five = M.end_kind(a, "right") == "5'" }
end

--- Terminal homology between the end of `a` and the start of `b`: the longest
--- exact match of at least opts.min_homology bases.
function M.homology(a, b, opts)
  opts = opts or {}
  local lo = opts.min_homology or 20
  local hi = opts.max_homology or 200
  local x, y = seqof(a):upper(), seqof(b):upper()
  local lim = math.min(#x, #y, hi)
  for n = lim, lo, -1 do
    if x:sub(#x - n + 1) == y:sub(1, n) then
      return { kind = "homology", len = n, seq = y:sub(1, n) }
    end
  end
  return nil, ("no terminal homology of %d bp or more"):format(lo)
end

---------------------------------------------------------------------------
-- Fragments
---------------------------------------------------------------------------

--- Fragment over `start`..`stop` of `rec` (wrapping when circular), extended
--- by `pad_left` / `pad_right` bases of single-stranded overhang. The record
--- keeps the features of the double-stranded core; the overhang bases are
--- plain sequence, so a feature never appears in two neighbouring fragments.
function M.from_span(rec, start, stop, left, right, name, pad_left, pad_right)
  local L = rec:len()
  local len = (stop - start) % L + 1
  local copy = require("plasmid.digest").extract(rec, { start = start, stop = stop, length = len })
  local seq = rec.seq:tostring()
  local function at(a, b)   -- a..b in source coordinates, wrapping
    local out = {}
    for p = a, b do out[#out + 1] = seq:sub(wrap(p, L), wrap(p, L)) end
    return table.concat(out)
  end
  if (pad_right or 0) > 0 then
    local n = copy:len()
    T.apply(copy, { op = "splice", start = n + 1, stop = n, ins = at(stop + 1, stop + pad_right) })
  end
  if (pad_left or 0) > 0 then
    T.apply(copy, { op = "splice", start = 1, stop = 0, ins = at(start - pad_left, start - 1) })
  end
  copy.history = nil
  copy.description = ""
  copy.name = name or rec.name
  return { name = name or rec.name, source = rec.name, rec = copy,
    left = left or M.blunt(), right = right or M.blunt() }
end

--- A whole linear record as a blunt fragment.
function M.from_record(rec, name)
  assert(not rec:is_circular(), "a circular record has no ends; cut it first")
  local copy = M.copy(rec)
  copy.description = ""
  return { name = name or rec.name, source = rec.name, rec = copy,
    left = M.blunt(), right = M.blunt() }
end

-- Ends of a digest fragment, and how far its single strands reach past the
-- top-strand cuts that bound it. `fr.left` / `fr.right` are the cut sites from
-- bio.restriction (nil at the end of a linear molecule).
--
-- A 5' overhang (o > 0) leaves the top strand protruding on the left of the
-- fragment and the bottom strand protruding on its right; a 3' overhang
-- (o < 0) does the opposite, so the fragment reaches `-o` bases further left.
local function ends_of(rec, fr)
  local left, right, pad_left, pad_right = M.blunt(), M.blunt(), 0, 0
  local circular = rec:is_circular()
  local L = rec:len()
  if fr.left then
    local o = fr.left[1].overhang
    left = { len = math.abs(o), strand = o > 0 and 1 or -1 }
    pad_left = math.max(0, -o)
    if not circular and fr.start - pad_left < 1 then   -- no bases there to protrude
      pad_left = math.max(0, fr.start - 1)
      left = { len = pad_left, strand = -1 }
    end
  end
  if fr.right then
    local o = fr.right[1].overhang
    right = { len = math.abs(o), strand = o > 0 and -1 or 1 }
    pad_right = math.max(0, o)
    if not circular and fr.stop + pad_right > L then
      pad_right = math.max(0, L - fr.stop)
      right = { len = pad_right, strand = -1 }
    end
  end
  return left, right, pad_left, pad_right
end

--- Fragments of a digest, as assembly fragments (same order as the digest).
function M.from_digest(rec, frags, name)
  local out = {}
  for i, fr in ipairs(frags) do
    if fr.uncut then
      if rec:is_circular() then
        -- an intact circle: no ends, so nothing to assemble with
        local copy = M.copy(rec)
        copy.description = ""
        out[#out + 1] = { name = name or rec.name, source = rec.name, rec = copy,
          left = M.blunt(), right = M.blunt(), circular = true, uncut = true }
      else
        out[#out + 1] = M.from_record(rec, name)
      end
    else
      local left, right, pl, pr = ends_of(rec, fr)
      local f = M.from_span(rec, fr.start, fr.stop, left, right,
        ("%s_%d"):format(name or rec.name, i), pl, pr)
      f.source = rec.name
      f.digest = fr
      out[#out + 1] = f
    end
  end
  return out
end

--- Cut a record with named enzymes and return the fragments.
function M.cut(rec, names, name)
  local digest = require("plasmid.digest")
  local frags, cuts, none, unknown = digest.compute(rec, names)
  return M.from_digest(rec, frags, name), cuts, none, unknown
end

--- The ends of the same fragment on the other strand, without copying the
--- record: enough for M.compatible, M.homology and M.end_text, but not for
--- M.assemble or M.build (there are no features on it).
function M.flip_ends(frag)
  local seq = alphabet.revcomp(seqof(frag))
  return { name = frag.name, source = frag.source, ends_only = true,
    rec = { seq = { tostring = function() return seq end } },
    left = { len = frag.right.len, strand = -frag.right.strand },
    right = { len = frag.left.len, strand = -frag.left.strand },
    flipped = not frag.flipped }
end

--- The same fragment on the other strand, record and features included.
function M.flip(frag)
  local copy = M.copy(frag.rec)
  T.apply(copy, { op = "revcomp", start = 1, stop = copy:len() })
  copy.history = nil
  return { name = frag.name, source = frag.source, rec = copy, digest = frag.digest,
    left = { len = frag.right.len, strand = -frag.right.strand },
    right = { len = frag.left.len, strand = -frag.left.strand },
    flipped = not frag.flipped }
end

---------------------------------------------------------------------------
-- Methods
---------------------------------------------------------------------------

M.TYPE_IIS = { "BsaI", "BsmBI", "Esp3I", "BbsI", "BpiI", "SapI", "BsmAI", "AarI" }

local function has_site(rec, e)
  local search = require("plasmid.bio.search")
  local seq = rec.seq:tostring()
  local motif = e.motif or search.parse_motif(e.site)
  local opts = { circular = rec:is_circular() }
  if #search.find_nuc(seq, motif, opts) > 0 then return true end
  local rc = search.revcomp_motif(motif)
  return rc ~= motif and #search.find_nuc(seq, rc, opts) > 0
end

--- Golden Gate parts: cut each record with a Type IIS enzyme and keep the
--- piece that no longer carries a site for it. `recs` is a list of records or
--- { rec = record, name = string }. Returns the fragments and a list of
--- problems (parts with no site, parts that yield several site-free pieces,
--- duplicate or palindromic overhangs, mixed overhang lengths).
function M.type_iis_parts(recs, enzyme_name, opts)
  local R = require("plasmid.bio.restriction")
  local e = R.get(enzyme_name)
  if not e then return nil, { ("unknown enzyme %s"):format(tostring(enzyme_name)) } end
  if e.cut5 >= 0 and e.cut5 <= #e.site then
    return nil, { ("%s cuts inside its recognition site; Golden Gate needs a Type IIS enzyme (%s)")
      :format(e.name, table.concat(M.TYPE_IIS, ", ")) }
  end
  if R.overhang(e) == 0 then
    return nil, { ("%s leaves blunt ends; Golden Gate needs an enzyme that leaves an overhang")
      :format(e.name) }
  end
  local out, problems = {}, {}
  for _, item in ipairs(recs) do
    local rec = item.rec or item
    local name = item.name or rec.name
    local frags = M.cut(rec, { e.name }, name)
    local keep = {}
    for _, f in ipairs(frags) do
      if f.digest and not has_site(f.rec, e) then keep[#keep + 1] = f end
    end
    if #frags == 1 and not frags[1].digest then
      problems[#problems + 1] = ("%s has no %s site"):format(name, e.name)
    elseif #keep == 0 then
      problems[#problems + 1] = ("%s: every piece still carries a %s site"):format(name, e.name)
    else
      if #keep > 1 then
        problems[#problems + 1] = ("%s gives %d pieces without a %s site; using the longest")
          :format(name, #keep, e.name)
        table.sort(keep, function(a, b) return a.rec:len() > b.rec:len() end)
      end
      keep[1].name = name
      out[#out + 1] = keep[1]
    end
  end
  for _, p in ipairs(M.overhang_problems(out)) do problems[#problems + 1] = p end
  return out, problems
end

--- Overhang problems in a set of fragments. Each overhang should appear
--- exactly twice, once on each of the two ends it joins: a third occurrence
--- makes the order ambiguous, a single one has nothing to join. Palindromic
--- overhangs anneal to themselves, and ends of different lengths cannot be
--- part of one one-pot reaction.
---
--- A part may be handed in on either strand, and the assembler is free to
--- flip it. Flipping rewrites an end's top-strand text as its reverse
--- complement and moves it to the other side, so an overhang and its reverse
--- complement are one overhang here, counted under whichever of the two reads
--- first. What flipping cannot change is whether two ends are the two halves
--- of one join: the side an end would take once its text is written in that
--- canonical direction is fixed, and two ends that would both be the left
--- half do not anneal however they are turned.
function M.overhang_problems(frags)
  local problems, where, order, lens = {}, {}, {}, {}
  for _, f in ipairs(frags) do
    for _, side in ipairs({ "left", "right" }) do
      local n = f[side].len
      lens[n] = true
      if n > 0 then
        local top = M.end_seq(f, side)
        local rc = alphabet.revcomp(top)
        if top == rc then
          problems[#problems + 1] = ("%s: %s overhang %s is palindromic"):format(f.name, side, top)
        end
        local key = top <= rc and top or rc
        local role = ((top <= rc) == (side == "left")) and "left" or "right"
        if not where[key] then where[key] = { text = top }; order[#order + 1] = key end
        table.insert(where[key], { at = ("%s %s"):format(f.name, side), role = role })
      end
    end
  end
  for _, key in ipairs(order) do
    local w = where[key]
    local at = {}
    for i, e in ipairs(w) do at[i] = e.at end
    if #w > 2 then
      problems[#problems + 1] = ("overhang %s is used %d times (%s)")
        :format(w.text, #w, table.concat(at, ", "))
    elseif #w == 1 then
      problems[#problems + 1] = ("overhang %s (%s) has no partner"):format(w.text, at[1])
    elseif w[1].role == w[2].role and w.text ~= alphabet.revcomp(w.text) then
      problems[#problems + 1] = ("overhang %s is on two ends of the same kind (%s), which do not anneal")
        :format(w.text, table.concat(at, ", "))
    end
  end
  local n = 0
  for _ in pairs(lens) do n = n + 1 end
  if n > 1 then
    local list = {}
    for k in pairs(lens) do list[#list + 1] = k end
    table.sort(list)
    problems[#problems + 1] = ("ends of different lengths: %s"):format(table.concat(list, ", "))
  end
  return problems
end

--- Gibson / HiFi parts: linear records (PCR products, cut vectors) as blunt
--- fragments. Circular records are rejected: they have no ends to overlap.
function M.gibson_parts(recs)
  local out, problems = {}, {}
  for _, item in ipairs(recs) do
    local rec = item.rec or item
    local name = item.name or rec.name
    if rec:is_circular() then
      problems[#problems + 1] = ("%s is circular; cut or amplify it first"):format(name)
    else
      out[#out + 1] = M.from_record(rec, name)
    end
  end
  return out, problems
end

---------------------------------------------------------------------------
-- Assembly search
---------------------------------------------------------------------------

local function joiner(opts)
  if opts.join == "homology" then
    return function(a, b) return M.homology(a, b, opts) end
  end
  return M.compatible
end

local function chain_key(chain)
  local t = {}
  for _, p in ipairs(chain) do t[#t + 1] = ("%d%s"):format(p.index, p.flipped and "-" or "+") end
  return table.concat(t, ",")
end

local function mirror_key(chain)
  local t = {}
  for i = #chain, 1, -1 do
    local p = chain[i]
    t[#t + 1] = ("%d%s"):format(p.index, p.flipped and "+" or "-")
  end
  return table.concat(t, ",")
end

--- Every assembly that uses each fragment exactly once. A fragment that is an
--- intact circle, or one from M.flip_ends, is a programming error here.
--- opts: { circular = true | false | nil (circular, else linear),
---         join = "ends" | "homology", min_homology, allow_flip = true,
---         max_solutions = 50 }
--- Returns a list of { parts = { { frag, index, flipped, junction } ... },
--- circular = bool } and, when no circular assembly closed and the search fell
--- back to linear products, true as a second value.
function M.assemble(frags, opts)
  opts = opts or {}
  local join = joiner(opts)
  for _, f in ipairs(frags) do
    assert(not f.circular, ("%s is an intact circle; cut it first"):format(f.name))
    assert(not f.ends_only, ("%s carries no record; use M.flip, not M.flip_ends"):format(f.name))
  end
  local allow_flip = opts.allow_flip ~= false
  local n = #frags
  local max_solutions = opts.max_solutions or 50
  local oriented = {}
  for i, f in ipairs(frags) do oriented[i] = { f, allow_flip and M.flip(f) or nil } end

  local circular, linear = {}, {}
  local seen_c, seen_l = {}, {}
  local used, chain = {}, {}
  -- The search walks every ordering of the fragments, so the work is
  -- factorial in their number whenever many ends are interchangeable (a
  -- one-enzyme digest, or blunt fragments, where every end joins every
  -- other). `max_solutions` ends the walk as well as the list: without that
  -- the cap bounds what is reported but not what is computed, and a dozen
  -- same-overhang fragments take hours.
  local full = false

  local function finish()
    if #circular + #linear >= max_solutions then full = true; return end
    local last, first = chain[#chain].frag, chain[1].frag
    local close = join(last, first)
    if close and opts.circular ~= false then
      local key = chain_key(chain)
      if not seen_c[key] then
        seen_c[key] = true
        local parts = {}
        for i, p in ipairs(chain) do parts[i] = { frag = p.frag, index = p.index, flipped = p.flipped,
          junction = i < #chain and chain[i].junction or close } end
        circular[#circular + 1] = { parts = parts, circular = true, closing = close }
      end
    end
    if not opts.circular then
      local key, mkey = chain_key(chain), mirror_key(chain)
      local k = key < mkey and key or mkey
      if not seen_l[k] and (n > 1 or not close) then
        seen_l[k] = true
        local parts = {}
        for i, p in ipairs(chain) do parts[i] = { frag = p.frag, index = p.index, flipped = p.flipped,
          junction = chain[i].junction } end
        linear[#linear + 1] = { parts = parts, circular = false }
      end
    end
  end

  local function dfs()
    if full then return end
    if #chain == n then return finish() end
    for i = 1, n do
      if not used[i] and not full then
        for fl = 1, (allow_flip and 2 or 1) do
          local f = oriented[i][fl]
          local j = join(chain[#chain].frag, f)
          if j then
            chain[#chain].junction = j
            used[i] = true
            chain[#chain + 1] = { frag = f, index = i, flipped = fl == 2 }
            dfs()
            chain[#chain] = nil
            used[i] = false
            chain[#chain].junction = nil
          end
        end
      end
    end
  end

  -- A circular assembly contains every fragment, so fixing the first one (in
  -- one orientation) removes rotations and the mirror image at no cost. A
  -- linear assembly has free ends, so every start has to be tried.
  local starts = {}
  if opts.circular == false then
    for i = 1, n do
      starts[#starts + 1] = { i, 1 }
      if allow_flip then starts[#starts + 1] = { i, 2 } end
    end
  else
    starts[1] = { 1, 1 }
  end
  for _, st in ipairs(starts) do
    if full then break end
    used = {}
    used[st[1]] = true
    chain = { { frag = oriented[st[1]][st[2]], index = st[1], flipped = st[2] == 2 } }
    dfs()
  end

  if opts.circular == false then return linear end
  if #circular > 0 or opts.circular then return circular end

  -- nothing closed: fall back to the linear assemblies
  local keep = opts.circular
  opts = vim.tbl_extend("force", opts, { circular = false })
  local out = M.assemble(frags, opts)
  opts.circular = keep
  return out, true
end

---------------------------------------------------------------------------
-- Building the product
---------------------------------------------------------------------------

local function shifted(node, off, L)
  if off % L == 0 then return E.clone(node) end
  return E.rotate(node, (-off) % L + 1, L)
end

local function copy_qualifiers(f)
  local out = {}
  for _, q in ipairs(f.qualifiers) do out[#out + 1] = { key = q.key, value = q.value, quoted = q.quoted } end
  return out
end

local function same_qualifiers(a, b)
  if #a.qualifiers ~= #b.qualifiers then return false end
  for i, q in ipairs(a.qualifiers) do
    local r = b.qualifiers[i]
    if q.key ~= r.key or q.value ~= r.value then return false end
  end
  return true
end

--- Cutting a record splits any feature that spans a cut. When the assembly
--- puts the two pieces back together (religation, or a part boundary inside a
--- feature), join them into one feature again. Only pieces from neighbouring
--- parts that abut exactly, with the same type, strand and qualifiers, merge.
function M.merge_split(rec, carried, nparts, circular)
  local L = rec:len()
  local function single(f)
    return f.loc and #f.parts == 1 and not f.parts[1].ref and f.parts[1].start and f.parts[1].stop
  end
  local merged = true
  while merged do
    merged = false
    for i = 1, #carried do
      local a = carried[i]
      for j = 1, #carried do
        local b = carried[j]
        local neighbours = (b.part - a.part) % nparts == 1 or nparts == 1
        if i ~= j and a.f and b.f and neighbours and a.f.type == b.f.type
          and a.f.strand == b.f.strand and single(a.f) and single(b.f)
          and same_qualifiers(a.f, b.f) then
          local pa, pb = a.f.parts[1], b.f.parts[1]
          local abut = pa.stop + 1 == pb.start
          local wraps = circular and pa.stop == L and pb.start == 1
          if abut or wraps then
            local text
            if abut then text = ("%d..%d"):format(pa.start, pb.stop)
            else text = ("join(%d..%d,1..%d)"):format(pa.start, L, pb.stop) end
            if a.f.strand == -1 then text = ("complement(%s)"):format(text) end
            local g = require("plasmid.model.record").Feature.new({ type = a.f.type,
              location_text = text, qualifiers = a.f.qualifiers, color = a.f.color })
            T.apply(rec, { op = "remove", id = b.f.id })
            T.apply(rec, { op = "replace", id = a.f.id, feature = g })
            a.f, b.f = g, nil
            a.part = b.part
            merged = true
          end
        end
      end
    end
  end
end

--- The assembled record: part sequences joined once per junction, every
--- part's features carried over, parts and junctions annotated.
--- opts: { name, annotate_parts = true, annotate_junctions = true, method }
function M.build(asm, opts)
  opts = opts or {}
  local model = require("plasmid.model.record")
  local Record, Feature = model.Record, model.Feature
  local parts = asm.parts
  local pieces, off = {}, {}
  local pos = 0
  for i, p in ipairs(parts) do
    local s = seqof(p.frag)
    if i == 1 then
      off[i] = 0
      pieces[1] = s
      pos = #s
    else
      local ov = parts[i - 1].junction.len
      off[i] = pos - ov
      pieces[#pieces + 1] = s:sub(ov + 1)
      pos = pos + #s - ov
    end
  end
  local total = pos
  if asm.circular then total = pos - parts[#parts].junction.len end
  local seq = table.concat(pieces):sub(1, total)
  assert(#seq == total and total > 0, "empty assembly")

  local first = parts[1].frag.rec
  local rec = Record.new({
    name = opts.name or (parts[1].frag.name .. "_assembly"),
    seq = seq,
    alphabet = first.alphabet,
    molecule_type = first.molecule_type,
    topology = asm.circular and "circular" or "linear",
    format = "genbank",
  })

  local function add(f) T.apply(rec, { op = "add", feature = f }) end

  local carried = {}
  for i, p in ipairs(parts) do
    for _, f in ipairs(p.frag.rec.features) do
      if f.type ~= "source" and f.loc then
        local node = shifted(f.loc, off[i], total)
        local g = Feature.new({ type = f.type, location_text = E.format(node),
          qualifiers = copy_qualifiers(f), color = f.color, origin = f.origin,
          primer_seq = f.primer_seq })
        add(g)
        carried[#carried + 1] = { f = g, part = i }
      end
    end
  end
  M.merge_split(rec, carried, #parts, asm.circular)

  if opts.annotate_parts ~= false then
    for i, p in ipairs(parts) do
      local len = #seqof(p.frag)
      local a = wrap(off[i] + 1, total)
      local b = wrap(off[i] + len, total)
      local loc = (a <= b) and ("%d..%d"):format(a, b) or ("join(%d..%d,1..%d)"):format(a, total, b)
      if p.flipped then loc = ("complement(%s)"):format(loc) end
      add(Feature.new({ type = "misc_feature", location_text = loc, qualifiers = {
        { key = "label", value = p.frag.name, quoted = true },
        { key = "note", value = ("part %d of %d%s"):format(i, #parts, p.flipped and ", reverse orientation" or ""),
          quoted = true },
      } }))
    end
  end

  if opts.annotate_junctions ~= false then
    for i, p in ipairs(parts) do
      local j = p.junction
      if j and j.len > 0 then
        local stop = off[i] + #seqof(p.frag)
        local a = wrap(stop - j.len + 1, total)
        local b = wrap(stop, total)
        local loc = (a <= b) and ("%d..%d"):format(a, b) or ("join(%d..%d,1..%d)"):format(a, total, b)
        local nxt = parts[i % #parts + 1]
        local label = j.kind == "homology" and ("overlap %d bp"):format(j.len) or j.seq
        add(Feature.new({ type = "misc_feature", location_text = loc, qualifiers = {
          { key = "label", value = label, quoted = true },
          { key = "note", value = ("%s junction: %s / %s"):format(
            j.kind == "homology" and "Gibson" or (j.kind == "blunt" and "blunt ligation" or "overhang"),
            p.frag.name, nxt and nxt.frag.name or "end"), quoted = true },
        } }))
      end
    end
  end

  rec.history = nil
  rec.description = opts.description or ""
  return rec, off, total
end

return M
