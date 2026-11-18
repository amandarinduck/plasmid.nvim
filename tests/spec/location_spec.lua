local H = require("helpers")
local L = require("plasmid.model.location")

local T = {}

local GRAMMAR = {
  "467", "<1", "340..565", "<345..500", "<1..888", "1..>888", "<1..>888",
  "102.110", "123^124", "complement(34..126)", "complement(<1..>200)",
  "join(12..78,134..202)", "complement(join(2691..4571,4918..5163))",
  "join(complement(4918..5163),complement(2691..4571))",
  "order(1..10,20..30)", "J00194.1:100..202", "join(1..100,J00194.1:100..202)",
  "join(4000..4500,1..200)", "complement(join(complement(1..5),6..9))",
}

function T.grammar_roundtrips_exactly()
  for _, s in ipairs(GRAMMAR) do
    local node, err = L.parse(s)
    H.ok(node, ("parse failed for %s: %s"):format(s, tostring(err)))
    H.eq(s, L.format(node), "round trip of " .. s)
  end
end

function T.rejects_malformed()
  for _, s in ipairs({ "", "join(1..2", "1....2", "complement()", "foo(1..2)", "1..2)", "<<1..2", "1..2,3" }) do
    H.ok(L.parse(s) == nil, "should reject: " .. s)
  end
end

function T.complement_join_flattens_in_biological_order()
  local a = L.flatten(L.parse("complement(join(10..20,30..40))"))
  local b = L.flatten(L.parse("join(complement(30..40),complement(10..20))"))
  H.eq({ { 30, 40, -1 }, { 10, 20, -1 } }, vim.tbl_map(function(p) return { p.start, p.stop, p.strand } end, a))
  H.eq(a, b)
end

function T.nested_complement_restores_strand()
  local p = L.flatten(L.parse("complement(join(complement(1..5),6..9))"))
  H.eq({ { 6, 9, -1 }, { 1, 5, 1 } }, vim.tbl_map(function(x) return { x.start, x.stop, x.strand } end, p))
  H.eq(0, L.strand(p))
end

function T.between_is_empty_interval()
  local p = L.flatten(L.parse("123^124"))[1]
  H.eq({ 124, 123 }, { p.start, p.stop })
  H.eq(0, L.length({ p }))
end

function T.length_skips_remote_parts()
  H.eq(10, L.length(L.flatten(L.parse("join(1..10,J00194.1:100..202)"))))
end

function T.extents_detect_origin_wrap()
  local wrap = L.flatten(L.parse("join(4000..4500,1..200)"))
  H.eq({ { 4000, 5000 }, { 1, 200 } }, L.extents(wrap, 5000, true))
  H.eq({ { 1, 4500 } }, L.extents(wrap, 5000, false)) -- linear: no wrap
  local rwrap = L.flatten(L.parse("complement(join(4000..4500,1..200))"))
  H.eq({ { 4000, 5000 }, { 1, 200 } }, L.extents(rwrap, 5000, true))
  local exons = L.flatten(L.parse("join(1..100,3000..3100)"))
  H.eq({ { 1, 3100 } }, L.extents(exons, 4000, true)) -- forward order: not a wrap
end

function T.extract_reverse_strand_exons()
  local store = require("plasmid.model.store").new("AAAACCCCGGGGTTTT")
  local parts = L.flatten(L.parse("complement(join(1..4,13..16))"))
  -- revcomp(AAAA .. TTTT) = AAAA .. TTTT
  H.eq("AAAATTTT", L.extract(parts, store))
  H.eq({ 16, 15, 14, 13, 4, 3, 2, 1 }, L.positions(parts))
end

--- Every location string in the NCBI corpus must parse and re-format exactly.
function T.corpus_locations_roundtrip()
  local gb = require("plasmid.io.genbank")
  local count = 0
  for _, path in ipairs(H.fixtures("**/*.gb")) do
    local res = assert(gb.parse(H.read(path)))
    for _, rec in ipairs(res.records) do
      for _, f in ipairs(rec.features) do
        H.ok(f.loc, ("%s: %s did not parse: %s"):format(path, f.location_text, tostring(f.loc_error)))
        H.eq(f.location_text, L.format(f.loc), path)
        count = count + 1
      end
    end
  end
  -- every feature key line in the raw files must have produced a feature
  H.eq(H.raw_feature_count(), count, "features parsed vs feature key lines in the files")
end

--- Property test: random location ASTs survive format -> parse -> format.
function T.random_ast_roundtrip()
  math.randomseed(7)
  local function rnd_simple()
    local a = math.random(1, 5000)
    local r = math.random()
    if r < 0.5 then
      return { kind = "range", start = a, stop = a + math.random(0, 500),
        start_mark = math.random() < 0.2 and "<" or nil, stop_mark = math.random() < 0.2 and ">" or nil }
    elseif r < 0.7 then
      return { kind = "site", pos = a }
    elseif r < 0.85 then
      return { kind = "between", left = a, right = a + 1 }
    end
    return { kind = "within", start = a, stop = a + 5 }
  end
  local function rnd(depth)
    local r = math.random()
    if depth > 3 or r < 0.4 then return rnd_simple() end
    if r < 0.6 then return { kind = "complement", child = rnd(depth + 1) } end
    local parts = {}
    for _ = 1, math.random(1, 4) do parts[#parts + 1] = rnd(depth + 1) end
    return { kind = math.random() < 0.8 and "join" or "order", parts = parts }
  end
  for _ = 1, 2000 do
    local s = L.format(rnd(0))
    local node = assert(L.parse(s), s)
    H.eq(s, L.format(node))
  end
end

return T
