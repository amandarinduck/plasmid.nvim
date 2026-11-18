local H = require("helpers")
local location = require("plasmid.model.location")
local edit = require("plasmid.model.edit")
local tx = require("plasmid.model.transaction")
local pio = require("plasmid.io")
local alphabet = require("plasmid.model.alphabet")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function sp(loc, s, e, k, L, circ)
  local node = edit.splice(assert(location.parse(loc)), s, e, k, L or 100, circ)
  return node and location.format(node) or nil
end

case("splice_rules", function()
  -- pure insertion (k bases before base s)
  H.eq("10..23", sp("10..20", 15, 14, 3), "inside extends")
  H.eq("13..23", sp("10..20", 10, 9, 3), "at the start boundary shifts")
  H.eq("10..20", sp("10..20", 21, 20, 3), "right after the end: unchanged")
  H.eq("13..23", sp("10..20", 5, 4, 3), "before shifts")
  H.eq("<10..>23", sp("<10..>20", 15, 14, 3), "fuzzy marks kept")
  H.eq("complement(join(10..23,33..43))", sp("complement(join(10..20,30..40))", 15, 14, 3), "form kept")
  -- deletion
  H.eq("10..17", sp("10..20", 12, 14, 0), "internal")
  H.eq("5..12", sp("10..20", 5, 12, 0), "start cut")
  H.eq("10..17", sp("10..20", 18, 25, 0), "end cut")
  H.eq(nil, sp("10..20", 5, 25, 0), "whole feature deleted")
  H.eq("10..20", sp("join(10..20,30..40)", 28, 45, 0), "a deleted part leaves the join")
  H.eq("complement(join(10..18,28..38))", sp("complement(join(10..20,30..40))", 12, 13, 0))
  H.eq(nil, sp("15", 15, 15, 0), "site deleted")
  H.eq("12", sp("15", 1, 3, 0), "site shifted")
  -- replacement (k new bases for s..e)
  H.eq("10..22", sp("10..20", 12, 14, 5), "spanning part takes the new bases")
  H.eq("10..21", sp("10..20", 18, 25, 4), "part touching on the left takes them")
  H.eq("5..14", sp("10..20", 5, 12, 2), "part touching on the right takes them")
  H.eq("5..6", sp("10..20", 5, 25, 2), "a wholly replaced feature becomes the new bases")
  H.eq("join(10..21,22..29)", sp("join(10..20,30..40)", 18, 32, 4), "only one part of a join takes them")
  H.eq("join(10..13,14..21)", sp("join(10..14,20..25,30..40)", 12, 32, 2), "the left part takes them")
  H.eq("join(3..4,9..19)", sp("join(5..8,12..14,20..30)", 3, 15, 2), "two parts wiped: the leftmost becomes the new bases")
  -- origin of a circular record (L = 100)
  H.eq("join(95..105,1..15)", sp("join(90..100,1..10)", 1, 0, 5, 100, true), "insert before base 1")
  H.eq("join(90..105,1..10)", sp("join(90..100,1..10)", 101, 100, 5, 100, true), "insert after base L")
  H.eq("join(95..105,6..15)", sp("join(90..100,1..10)", 1, 0, 5, 100, false), "linear: no junction rule")
  -- between
  H.eq("10^11", sp("10^11", 11, 10, 4), "insertion in the gap keeps the site on the left")
  H.eq("14^15", sp("10^11", 5, 4, 4))
  H.eq("100^1", sp("100^1", 50, 49, 0, 100, true), "origin between untouched by an interior insertion")
end)

case("flip_and_rotate_rules", function()
  local function fl(loc, s, e) return location.format(edit.flip(assert(location.parse(loc)), s, e)) end
  H.eq("81..91", fl("complement(10..20)", 1, 100))
  H.eq("complement(join(61..71,81..91))", fl("join(10..20,30..40)", 1, 100))
  H.eq("join(61..71,81..91)", fl("complement(join(10..20,30..40))", 1, 100))
  H.eq("join(61..71,81..91)", fl("join(complement(30..40),complement(10..20))", 1, 100), "biological order kept")
  H.eq("join(complement(81..91),61..71)", fl("join(10..20,complement(30..40))", 1, 100), "mixed: part by part")
  H.eq("complement(81..>91)", fl("<10..20", 1, 100), "fuzzy marks swap ends")
  H.eq("complement(35..45)", fl("15..25", 10, 50), "region-relative mirror")
  local function ro(loc, o, L) return location.format(edit.rotate(assert(location.parse(loc)), o, L)) end
  H.eq("join(96..100,1..6)", ro("10..20", 15, 100), "a feature crossing the new origin splits")
  H.eq("complement(join(96..100,1..6))", ro("complement(10..20)", 15, 100))
  H.eq("41..61", ro("join(90..100,1..10)", 50, 100), "a former origin junction merges")
  H.eq("complement(41..61)", ro("complement(join(90..100,1..10))", 50, 100))
  H.eq("81..91", ro("30..40", 50, 100))
  H.eq("join(81..91,96..100,1..6)", ro("join(30..40,45..55)", 50, 100))
end)


case("whole_record_revcomp_and_rotation_keep_feature_sequences", function()
  for _, rel in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001422.gb", "snapgene/sgffp_origin_spanning_features.dna",
    "snapgene/pFA-KanMX4.dna" }) do
    local file = assert(pio.read(H.fixture(rel)))
    local rec = file.records[1]
    local L = rec:len()
    local res, locs = {}, {}
    for _, f in ipairs(rec.features) do
      if f.loc and not location.has_remote(f.parts) then res[f.id] = rec:extract(f) end
      locs[f.id] = f.location_text
    end
    local function same_residues(what)
      for _, f in ipairs(rec.features) do
        if res[f.id] then H.eq(res[f.id], rec:extract(f), ("%s %s: %s %s"):format(rel, what, f.type, f.location_text)) end
      end
    end
    local orig = rec.seq:tostring()
    tx.commit(rec, { op = "revcomp", start = 1, stop = L })
    H.ok(alphabet.revcomp(orig) == rec.seq:tostring(), rel .. ": residues reverse-complemented")
    same_residues("after whole-record revcomp")
    tx.commit(rec, { op = "revcomp", start = 1, stop = L })
    for _, f in ipairs(rec.features) do H.eq(locs[f.id], f.location_text, rel .. " revcomp twice") end
    if rec:is_circular() then
      for _, o in ipairs({ 2, math.floor(L / 3), L }) do
        tx.commit(rec, { op = "rotate", origin = o })
        same_residues("after rotating to " .. o)
        tx.commit(rec, { op = "rotate", origin = (L - o + 1) % L + 1 })
        for _, f in ipairs(rec.features) do H.eq(locs[f.id], f.location_text, ("%s rotate %d and back"):format(rel, o)) end
      end
    end
  end
end)

---------------------------------------------------------------------------
-- Property test: base tagging oracle
---------------------------------------------------------------------------

-- Oracle: the sequence is an array of { tag, base }; a feature is a list of
-- parts in biological order, each { tags (forward order), strand, kind }.
-- Its bases in biological order are compared with the plugin's.
local function new_oracle(rec)
  local o = { seq = {}, feats = {}, next_tag = 0, circ = rec:is_circular() }
  local s = rec.seq:tostring()
  for i = 1, #s do o.next_tag = o.next_tag + 1; o.seq[i] = { tag = o.next_tag, base = s:sub(i, i) } end
  for _, f in ipairs(rec.features) do
    local ok = f.loc ~= nil
    local parts = {}
    if ok then
      for _, p in ipairs(location.flatten(f.loc)) do
        if p.ref or (p.kind ~= "range" and p.kind ~= "site") or p.stop < p.start then ok = false break end
        local tags = {}
        for x = p.start, p.stop do tags[#tags + 1] = o.seq[x].tag end
        parts[#parts + 1] = { tags = tags, strand = p.strand, kind = p.kind }
      end
    end
    if ok then o.feats[f.id] = { parts = parts } end
  end
  o.tracked = {}
  for id in pairs(o.feats) do o.tracked[id] = true end
  return o
end

local function pos_index(o)
  local idx = {}
  for i, c in ipairs(o.seq) do idx[c.tag] = i end
  return idx
end

local function fresh(o, n)
  local t = {}
  for i = 1, n do o.next_tag = o.next_tag + 1; t[i] = o.next_tag end
  return t
end

-- Oracle version of a splice with the documented rules.
local function oracle_splice(o, s, e, ins)
  local L, n, k = #o.seq, e - s + 1, #ins
  local newtags = fresh(o, k)
  if n == k and n > 0 then -- substitution: new bases take the old ones' places
    local map = {}
    for i = s, e do map[o.seq[i].tag] = newtags[i - s + 1] end
    for i = s, e do o.seq[i] = { tag = newtags[i - s + 1], base = ins:sub(i - s + 1, i - s + 1) } end
    for _, F in pairs(o.feats) do
      for _, p in ipairs(F.parts) do
        for j, t in ipairs(p.tags) do if map[t] then p.tags[j] = map[t] end end
      end
    end
    return
  end
  local idx0 = pos_index(o)
  local dead = {}
  for i = s, e do dead[o.seq[i].tag] = true end
  for _, F in pairs(o.feats) do
    for _, p in ipairs(F.parts) do
      local keep = {}
      p.orig_a = idx0[p.tags[1]]
      for _, t in ipairs(p.tags) do if not dead[t] then keep[#keep + 1] = t end end
      p.touched = #keep < #p.tags
      p.whole = #keep == 0
      p.tags = keep
    end
  end
  for i = e, s, -1 do table.remove(o.seq, i) end
  if k > 0 then
    local idx = pos_index(o)
    local function first(p) return idx[p.tags[1]] end
    local function last(p) return idx[p.tags[#p.tags]] end
    for _, F in pairs(o.feats) do
      local special
      if n > 0 then
        local spanning = false
        for _, p in ipairs(F.parts) do
          if #p.tags > 0 and first(p) < s and last(p) >= s then spanning = true end
        end
        if not spanning then
          for _, p in ipairs(F.parts) do
            if not special and p.touched and #p.tags > 0 and last(p) == s - 1 then special = { p, "append" } end
          end
          for _, p in ipairs(F.parts) do
            if not special and p.touched and #p.tags > 0 and first(p) == s then special = { p, "prepend" } end
          end
          local best
          for _, p in ipairs(F.parts) do
            if p.whole and (not best or p.orig_a < best.orig_a) then best = p end
          end
          if not special and best then special = { best, "become" } end
        end
      elseif o.circ and (s == 1 or s == L + 1) then
        local at_end, at_start
        for _, p in ipairs(F.parts) do
          if p.kind == "range" and #p.tags > 0 then
            if last(p) == L then at_end = p end
            if first(p) == 1 then at_start = p end
          end
        end
        if at_end and at_start and at_end ~= at_start then
          special = s == 1 and { at_start, "prepend" } or { at_end, "append" }
        end
      end
      for _, p in ipairs(F.parts) do
        if special and p == special[1] then
          local how = special[2]
          if how == "become" then p.tags = { unpack(newtags) }
          elseif how == "append" then for _, t in ipairs(newtags) do p.tags[#p.tags + 1] = t end
          else for j = #newtags, 1, -1 do table.insert(p.tags, 1, newtags[j]) end end
          p.kind = "range"
        elseif #p.tags > 0 and p.kind == "range" and first(p) < s and last(p) >= s then
          local at = s - first(p)
          for j = #newtags, 1, -1 do table.insert(p.tags, at + 1, newtags[j]) end
        end
      end
    end
  end
  for i = 1, k do table.insert(o.seq, s + i - 1, { tag = newtags[i], base = ins:sub(i, i) }) end
  for id, F in pairs(o.feats) do
    local parts = {}
    for _, p in ipairs(F.parts) do
      p.touched, p.whole, p.orig_a = nil, nil, nil
      if #p.tags > 0 then parts[#parts + 1] = p end
    end
    if #parts == 0 then o.feats[id] = nil else F.parts = parts end
  end
end

local function oracle_revcomp(o, s, e, inside_ids)
  local idx = pos_index(o)
  local runs = {}
  for id, F in pairs(o.feats) do
    if not inside_ids[id] then
      local r = {}
      for j, p in ipairs(F.parts) do r[j] = { idx[p.tags[1]], idx[p.tags[#p.tags]] } end
      runs[id] = r
    end
  end
  local seg = {}
  for i = s, e do seg[#seg + 1] = o.seq[i] end
  for i = s, e do
    local c = seg[e - i + 1]
    o.seq[i] = { tag = c.tag, base = alphabet.complement(c.base) }
  end
  for id, F in pairs(o.feats) do
    if inside_ids[id] then
      -- bases move with the inversion: same biological order, other strand
      for _, p in ipairs(F.parts) do
        local r = {}
        for q = #p.tags, 1, -1 do r[#r + 1] = p.tags[q] end
        p.tags, p.strand = r, -p.strand
      end
    else
      -- coordinates stay; the bases there are now different ones
      for j, p in ipairs(F.parts) do
        local a, b = runs[id][j][1], runs[id][j][2]
        local tags = {}
        for x = a, b do tags[#tags + 1] = o.seq[x].tag end
        p.tags = tags
      end
    end
  end
end

local function oracle_rotate(o, org)
  local L = #o.seq
  local idx = pos_index(o)
  for _, F in pairs(o.feats) do
    for _, p in ipairs(F.parts) do
      p.was_end = idx[p.tags[#p.tags]] == L
      p.was_start = idx[p.tags[1]] == 1
    end
  end
  local new = {}
  for i = org, L do new[#new + 1] = o.seq[i] end
  for i = 1, org - 1 do new[#new + 1] = o.seq[i] end
  o.seq = new
  idx = pos_index(o)
  for _, F in pairs(o.feats) do
    local parts = {}
    for _, p in ipairs(F.parts) do
      -- split a part that now runs across the origin
      local runs, cur, prev = {}, nil, nil
      for _, t in ipairs(p.tags) do
        if not cur or idx[t] ~= idx[prev] + 1 then cur = { tags = {}, strand = p.strand, kind = p.kind }; runs[#runs + 1] = cur end
        cur.tags[#cur.tags + 1] = t
        prev = t
      end
      runs[1].was_start = p.was_start
      runs[#runs].was_end = p.was_end
      if p.strand == -1 then
        for q = #runs, 1, -1 do parts[#parts + 1] = runs[q] end
      else
        for q = 1, #runs do parts[#parts + 1] = runs[q] end
      end
    end
    -- merge a former origin junction that is now contiguous
    local merged = {}
    for _, p in ipairs(parts) do
      local q = merged[#merged]
      local joined = false
      if q and q.kind == "range" and p.kind == "range" and q.strand == p.strand then
        if p.strand ~= -1 and q.was_end and p.was_start and idx[p.tags[1]] == idx[q.tags[#q.tags]] + 1 then
          for _, t in ipairs(p.tags) do q.tags[#q.tags + 1] = t end
          q.was_end = p.was_end
          joined = true
        elseif p.strand == -1 and q.was_start and p.was_end and idx[q.tags[1]] == idx[p.tags[#p.tags]] + 1 then
          local tags = {}
          for _, t in ipairs(p.tags) do tags[#tags + 1] = t end
          for _, t in ipairs(q.tags) do tags[#tags + 1] = t end
          q.tags, q.was_start = tags, p.was_start
          joined = true
        end
      end
      if not joined then merged[#merged + 1] = p end
    end
    for _, p in ipairs(merged) do p.was_end, p.was_start = nil, nil end
    F.parts = merged
  end
end

-- Compare plugin record with oracle (bases of each feature in biological order).
local function check(rec, o, where)
  local s = rec.seq:tostring()
  H.eq(#o.seq, #s, where .. ": length")
  local bases = {}
  for i, c in ipairs(o.seq) do bases[i] = c.base end
  H.eq(table.concat(bases), s, where .. ": residues")
  local seen = 0
  for _, f in ipairs(rec.features) do
    H.ok(f.loc, where .. ": unparsable location " .. tostring(f.location_text))
    for _, p in ipairs(f.parts) do
      if not p.ref then
        H.ok(p.start >= 1 and p.stop <= #s and p.start <= p.stop + 1,
          ("%s: %s out of range (L=%d)"):format(where, f.location_text, #s))
      end
    end
    local F = o.feats[f.id]
    if F then
      seen = seen + 1
      local got, want = {}, {}
      for _, x in ipairs(location.positions(f.parts)) do got[#got + 1] = o.seq[x].tag end
      for _, p in ipairs(F.parts) do
        if p.strand == -1 then
          for q = #p.tags, 1, -1 do want[#want + 1] = p.tags[q] end
        else
          for _, t in ipairs(p.tags) do want[#want + 1] = t end
        end
      end
      if table.concat(got, ",") ~= table.concat(want, ",") then
        error(("%s: feature %d %s covers the wrong bases"):format(where, f.id, f.location_text), 0)
      end
    end
  end
  local ids = {}
  for _, f in ipairs(rec.features) do ids[f.id] = true end
  for id in pairs(o.feats) do H.ok(ids[id], where .. ": feature " .. id .. " vanished") end
  -- and the other way: features the oracle destroyed are gone from the record
  for id in pairs(o.tracked) do
    if not o.feats[id] then H.ok(not ids[id], where .. ": feature " .. id .. " should have been destroyed") end
  end
  return seen
end

local function rnd_bases(n)
  local t = {}
  for i = 1, n do local r = math.random(4); t[i] = ("ACGT"):sub(r, r) end
  return table.concat(t)
end

case("random_edits_track_tagged_bases_and_undo_exactly", function()
  -- PLASMID_EDIT_SEEDS=n repeats the run with n seeds (default 1)
  local total_ops, total_checks = 0, 0
  local nseeds = tonumber(os.getenv("PLASMID_EDIT_SEEDS") or "1")
  for seed = 20261116, 20261116 + nseeds - 1 do
  math.randomseed(seed)
  for _, rel in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001422.gb", "ncbi/U01317.gb", "synthetic/grammar.gb",
    "snapgene/sgffp_origin_spanning_features.dna", "snapgene/pFA-KanMX4.dna" }) do
    local path = H.fixture(rel)
    local text = H.read(path)
    local file = assert(pio.read(path))
    local rec = file.records[1]
    local original_feats = {}
    for i, f in ipairs(rec.features) do original_feats[i] = f end
    local original_seq = rec.seq:tostring()
    local o = new_oracle(rec)
    local nops = 0
    for step = 1, 60 do
      local L = rec:len()
      local r = math.random()
      local where = ("%s seed %d op %d"):format(rel, seed, step)
      local before_locs = {}
      for _, f in ipairs(rec.features) do before_locs[f.id] = f.location_text end
      local desc
      if r < 0.25 then
        local s = math.random(1, L + 1)
        local ins = rnd_bases(math.random(1, 12))
        desc = ("insert %d bp before %d (L=%d)"):format(#ins, s, L)
        tx.commit(rec, { op = "splice", start = s, stop = s - 1, ins = ins })
        oracle_splice(o, s, s - 1, ins)
      elseif r < 0.45 then
        local s = math.random(1, L)
        local e = math.min(L, s + math.random(0, 25))
        if e - s + 1 < L - 10 then
          desc = ("delete %d..%d (L=%d)"):format(s, e, L)
          tx.commit(rec, { op = "splice", start = s, stop = e, ins = "" })
          oracle_splice(o, s, e, "")
        else
          nops = nops - 1
        end
      elseif r < 0.70 then
        local s = math.random(1, L)
        local e = math.min(L, s + math.random(0, 20))
        local ins = rnd_bases(math.random() < 0.3 and (e - s + 1) or math.random(1, 15))
        desc = ("replace %d..%d with %d bp (L=%d)"):format(s, e, #ins, L)
        tx.commit(rec, { op = "splice", start = s, stop = e, ins = ins })
        oracle_splice(o, s, e, ins)
      elseif r < 0.85 then
        local s = math.random(1, L)
        local e = math.min(L, s + math.random(5, 400))
        local inside = {}
        for _, f in ipairs(rec.features) do if f.loc and edit.inside(f.loc, s, e) then inside[f.id] = true end end
        desc = ("revcomp %d..%d (L=%d)"):format(s, e, L)
        tx.commit(rec, { op = "revcomp", start = s, stop = e })
        oracle_revcomp(o, s, e, inside)
      elseif rec:is_circular() then
        local org = math.random(1, L)
        desc = ("rotate to %d (L=%d)"):format(org, L)
        tx.commit(rec, { op = "rotate", origin = org })
        oracle_rotate(o, org)
      else
        nops = nops - 1
      end
      nops = nops + 1
      if os.getenv("EDIT_DEBUG") and rel:find(os.getenv("EDIT_DEBUG"), 1, true) then
        local idx = pos_index(o)
        local F = o.feats[tonumber(os.getenv("EDIT_FID") or "1")]
        local parts = {}
        for _, p in ipairs(F and F.parts or {}) do
          parts[#parts + 1] = ("%d..%d%s"):format(idx[p.tags[1]], idx[p.tags[#p.tags]], p.strand == -1 and "-" or "+")
        end
        local f
        for _, g in ipairs(rec.features) do if g.id == tonumber(os.getenv("EDIT_FID") or "1") then f = g end end
        H.log(("  %2d %-34s plugin %-40s oracle %s"):format(step, tostring(desc), f and f.location_text or "-", table.concat(parts, " ")))
      end
      local ok, err = pcall(check, rec, o, where .. " [" .. tostring(desc) .. "]")
      if not ok then
        local id = tonumber(tostring(err):match("feature (%d+)"))
        error(tostring(err) .. (id and ("; before: " .. tostring(before_locs[id])) or ""), 0)
      end
      total_checks = total_checks + err
    end
    total_ops = total_ops + nops
    -- undo everything: same objects, same residues, same file bytes
    while tx.undo(rec) do end
    H.eq(original_seq, rec.seq:tostring(), rel .. ": residues after undo")
    H.eq(#original_feats, #rec.features, rel .. ": feature count after undo")
    for i, f in ipairs(original_feats) do
      if rec.features[i] ~= f then error(rel .. ": feature " .. i .. " is not the original object after undo", 0) end
    end
    if rel:match("%.gb$") then
      local out = assert(pio.format(file.records, file.format, file.file))
      H.eq(text, out, rel .. ": written file after undo")
    end
    -- redo everything, then undo again
    local redone = 0
    while tx.redo(rec) do redone = redone + 1 end
    H.eq(nops, redone)
    while tx.undo(rec) do end
    H.eq(original_seq, rec.seq:tostring(), rel .. ": residues after redo+undo")
  end
  end
  H.log(("    [edit] %d random edits on 6 records (%d seed%s); %d feature checks against the tagged-base oracle"):format(
    total_ops, nseeds, nseeds == 1 and "" or "s", total_checks))
end)

return T
