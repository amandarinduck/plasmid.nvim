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
