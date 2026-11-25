-- The cloning engine against pydna.
--
-- tests/oracle/dump_cloning.py runs several hundred assemblies through pydna
-- 5.5.16 / Biopython 1.88 and commits the answers to tests/oracle/expected/
-- cloning.json. Nothing here touches Python or the network: the JSON is read,
-- every case is run through bio.assembly, and the products are compared
-- through a canonical form that ignores rotation and strand.
--
-- docs/cloning_report.md writes up what agrees, what does not, and why.
local H = require("helpers")
local D = require("cloning_driver")

require("plasmid").setup({})
local alphabet = require("plasmid.model.alphabet")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local doc
local function oracle()
  if not doc then
    local path = H.root .. "/tests/oracle/expected/cloning.json"
    doc = vim.json.decode(H.read(path), { luanil = { object = true, array = true } })
  end
  return doc
end

local function run_method(method)
  local d = oracle()
  local n, bad = 0, {}
  for _, c in ipairs(d.cases) do
    if c.method == method then
      n = n + 1
      local complaints = D.compare(c, H.root)
      if #complaints > 0 then
        bad[#bad + 1] = ("%s (%s): %s"):format(c.id, c.label, table.concat(complaints, "; "))
      end
    end
  end
  return n, bad
end

---------------------------------------------------------------------------
-- The comparison machinery itself
---------------------------------------------------------------------------

case("canonical_form_ignores_rotation_and_strand", function()
  -- Booth's least rotation, against a naive search over every rotation.
  local function naive(s)
    local best
    for i = 1, #s do
      local r = s:sub(i) .. s:sub(1, i - 1)
      if not best or r < best then best = r end
    end
    return best
  end
  math.randomseed(7)
  local b = { "A", "C", "G", "T" }
  for _ = 1, 200 do
    local t = {}
    for i = 1, math.random(1, 40) do t[i] = b[math.random(4)] end
    local s = table.concat(t)
    H.eq(naive(s), D.least_rotation(s), "least rotation of " .. s)
  end
  H.eq("AAAA", D.least_rotation("AAAA"))
  H.eq("ACGTACGT", D.least_rotation("GTACGTAC"))

  -- a circular sequence has one canonical form however it is written
  local seq = "ATGGCCTAGAACCGTTAGCTAGCATCG"
  local want = D.canonical(seq, true)
  for i = 1, #seq do
    local rot = seq:sub(i) .. seq:sub(1, i - 1)
    H.eq(want, D.canonical(rot, true), "rotation " .. i)
    H.eq(want, D.canonical(alphabet.revcomp(rot), true), "reverse complement of rotation " .. i)
  end
  -- and a linear one reads the same either way, but is not canonical under rotation
  H.eq(D.canonical(seq, false), D.canonical(alphabet.revcomp(seq), false))
  H.ok(D.canonical(seq, false) ~= D.canonical(seq:sub(2) .. seq:sub(1, 1), false),
    "a linear sequence is not canonicalised over rotations")
end)

case("the_committed_oracle_is_complete", function()
  local d = oracle()
  H.ok(d.pydna and d.biopython, "the oracle records which versions wrote it")
  local seen, n = {}, 0
  for _, c in ipairs(d.cases) do
    n = n + 1
    H.ok(not seen[c.id], "duplicate case id " .. tostring(c.id))
    seen[c.id] = true
    H.ok(c.expect and c.expect.circular, tostring(c.id) .. " has no expectation")
  end
  local total = 0
  local methods = vim.tbl_keys(d.counts)
  table.sort(methods)
  for _, method in ipairs(methods) do
    local m = 0
    for _, c in ipairs(d.cases) do if c.method == method then m = m + 1 end end
    H.eq(d.counts[method], m, "case count for " .. method)
    total = total + d.counts[method]
  end
  H.eq(total, n, "every case belongs to a counted method")
  local parts = {}
  for _, method in ipairs(methods) do parts[#parts + 1] = ("%s %d"):format(method, d.counts[method]) end
  H.log(("    [oracle] %d cases from pydna %s / biopython %s: %s; %d hand-built")
    :format(n, d.pydna, d.biopython, table.concat(parts, ", "), d.hand_built))
end)

---------------------------------------------------------------------------
-- The three methods
---------------------------------------------------------------------------

case("restriction_ligation_matches_pydna", function()
  local n, bad = run_method("ligation")
  H.eq({}, bad, ("%d of %d ligation cases disagree"):format(#bad, n))
  H.log(("    [oracle] %d restriction-ligation cases agree with pydna"):format(n))
end)

case("gibson_matches_pydna", function()
  local n, bad = run_method("gibson")
  H.eq({}, bad, ("%d of %d Gibson cases disagree"):format(#bad, n))
  H.log(("    [oracle] %d Gibson cases agree with pydna"):format(n))
end)

case("golden_gate_matches_pydna", function()
  local n, bad = run_method("goldengate")
  H.eq({}, bad, ("%d of %d Golden Gate cases disagree"):format(#bad, n))
  H.log(("    [oracle] %d Golden Gate cases agree with pydna or with the hand-built answer")
    :format(n))
end)

---------------------------------------------------------------------------
-- Single fragments: partial digests and end chemistry
---------------------------------------------------------------------------

-- The oracle writes a fragment as pydna does: the length of everything
-- present on either strand, the two overhangs in pydna's sign convention
-- (negative = 5' protruding, positive = 3'), and a hash of the full
-- sequence.
local A = require("plasmid.bio.assembly")
local DG = require("plasmid.digest")
local pio = require("plasmid.io")

local function frag_key(f)
  local seq = A.sequence(f):upper()
  local l = f.left.len == 0 and 0 or (f.left.strand == 1 and -f.left.len or f.left.len)
  local r = f.right.len == 0 and 0 or (f.right.strand == -1 and -f.right.len or f.right.len)
  return { #seq, l, r, vim.fn.sha256(seq):sub(1, 16) }
end

local function key_text(k) return ("%d/%d/%d/%s"):format(k[1], k[2], k[3], k[4]) end

local function sort_keys(list)
  table.sort(list, function(a, b)
    for i = 1, 3 do if a[i] ~= b[i] then return a[i] < b[i] end end
    return a[4] < b[4]
  end)
  return list
end

local fixture_of = {
  J01749 = "ncbi/J01749.gb", L09137 = "ncbi/L09137.gb",
  NC_001422 = "ncbi/NC_001422.gb", NC_005816 = "ncbi/NC_005816.gb",
  gg_insert = "synthetic/gg_insert.gb", gg_vector = "synthetic/gg_vector.gb",
}
local loaded = {}
--- The fixture as a real record of the asked-for topology. A metatable that
--- only overrides `topology` is not enough here: the fragment machinery
--- writes the record out and reads it back, and the written LOCUS line would
--- carry the fixture's own topology.
local function record(stem, topology)
  local k = stem .. "/" .. topology
  if not loaded[k] then
    local copy = A.copy(assert(pio.read(H.fixture(fixture_of[stem]))).records[1])
    if copy.topology ~= topology then
      require("plasmid.model.transaction").apply(copy, { op = "topology", value = topology })
    end
    copy.history = nil
    loaded[k] = copy
  end
  return loaded[k]
end

local function known(names)
  local R = require("plasmid.bio.restriction")
  for _, n in ipairs(names) do if not R.get(n) then return false, n end end
  return true
end

case("partial_digests_match_pydnas_own_fragment_extraction", function()
  local d = oracle()
  H.ok(d.partial and #d.partial > 0, "the oracle carries a partial-digest corpus")
  local sets, frags, unknown = 0, 0, {}
  for _, c in ipairs(d.partial) do
    local ok, which = known(c.enzymes)
    if not ok then
      unknown[which] = true
    else
      local rec = record(c.stem, c.topology)
      local parts = DG.compute(rec, c.enzymes, { partial = true, max_sites = 20, max_fragments = 100000 })
      local got, extra = {}, 0
      for _, f in ipairs(A.from_digest(rec, parts, c.stem)) do
        -- pydna enumerates pairs of cutsites, so it has no entry for the
        -- circle nobody cut; on a linear molecule that span is the pair of
        -- molecule ends and pydna does list it
        if f.circular then extra = extra + 1 else got[#got + 1] = frag_key(f) end
      end
      H.eq(c.topology == "circular" and 1 or 0, extra,
        ("%s %s %s: intact molecules"):format(c.stem, c.topology, table.concat(c.enzymes, "+")))
      sort_keys(got)
      local want, mine = {}, {}
      for i, k in ipairs(c.fragments) do want[i] = key_text(k) end
      for i, k in ipairs(got) do mine[i] = key_text(k) end
      H.eq(want, mine, ("%s %s %s (%d sites)"):format(c.stem, c.topology,
        table.concat(c.enzymes, "+"), c.sites))
      sets = sets + 1
      frags = frags + #got
    end
  end
  H.eq({}, vim.tbl_keys(unknown), "every enzyme the oracle used is in the bundled set")
  H.log(("    [oracle] %d partial digests (%d fragments) match pydna's apply_cut over every "
    .. "pair of cutsites, sequence and both overhangs"):format(sets, frags))
end)

case("end_chemistry_matches_pydna_where_it_can_model_it", function()
  local d = oracle()
  H.ok(d.ends and #d.ends > 0, "the oracle carries an end-chemistry corpus")
  local checked, chewback, skipped_chewback = 0, 0, 0
  for _, c in ipairs(d.ends) do
    if known(c.enzymes) then
      local rec = record(c.stem, c.topology)
      local by_cut = {}
      for _, p in ipairs(c.pieces) do by_cut[key_text(p.cut)] = p end
      for _, f in ipairs(A.cut(rec, c.enzymes, c.stem)) do
        if f.digest then
          local p = by_cut[key_text(frag_key(f))]
          H.ok(p, ("%s %s %s: pydna has no fragment %s"):format(c.stem, c.topology,
            table.concat(c.enzymes, "+"), key_text(frag_key(f))))
          for _, name in ipairs({ "fillin", "blunt", "nuclease", "chewback" }) do
            if p[name] then
              local g = assert(A.treat(f, name))
              H.eq(key_text(p[name]), key_text(frag_key(g)),
                ("%s %s %s %s"):format(c.stem, c.topology, table.concat(c.enzymes, "+"), name))
              checked = checked + 1
              if name == "chewback" then chewback = chewback + 1 end
            elseif name == "chewback" then
              skipped_chewback = skipped_chewback + 1
            end
          end
        end
      end
    end
  end
  local why = {}
  for _, u in ipairs(d.unmodelled or {}) do why[#why + 1] = u.what end
  H.log(("    [oracle] %d treated ends match pydna (fill_in, t4, mung; %d chew-backs compared, "
    .. "%d fragments carry a 5′ overhang so pydna's t4 is not the chew-back answer). "
    .. "Not modelled by pydna: %s"):format(checked, chewback, skipped_chewback,
    table.concat(why, "; ")))
end)

---------------------------------------------------------------------------
-- What the corpus is made of: a case that always passes is no test
---------------------------------------------------------------------------

case("the_corpus_exercises_both_outcomes", function()
  local d = oracle()
  local with, without, products, designed = 0, 0, 0, 0
  for _, c in ipairs(d.cases) do
    local k = #c.expect.circular
    if k > 0 then with = with + 1 else without = without + 1 end
    products = math.max(products, k)
    if c.designed then designed = designed + 1 end
  end
  H.ok(with > 300, ("only %d cases expect a product"):format(with))
  H.ok(without > 40, ("only %d cases expect no product"):format(without))
  H.ok(products > 1, "no case expects more than one product")
  H.log(("    [oracle] %d cases expect a circular product (up to %d distinct), %d expect none; "
    .. "%d are designed failures"):format(with, products, without, designed))
end)

return T
