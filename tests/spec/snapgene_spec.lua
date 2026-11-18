local H = require("helpers")
local pio = require("plasmid.io")
local sg = require("plasmid.io.snapgene")
local xml = require("plasmid.io.xml")
local gb = require("plasmid.io.genbank")

local T = {}

local function snapgene_files()
  local out = H.fixtures("snapgene/*.dna")
  -- files written by this plugin (authoring_spec), read back by Biopython
  vim.list_extend(out, H.fixtures("written/*.dna"))
  vim.list_extend(out, H.fixtures("snapgene/*.prot"))
  vim.list_extend(out, H.fixtures("snapgene/*.rna"))
  return out
end

local function sorted_parts(parts)
  table.sort(parts, function(a, b)
    if a[1] ~= b[1] then return a[1] < b[1] end
    return a[2] < b[2]
  end)
  return parts
end

local function sorted_quals(list)
  table.sort(list, function(a, b)
    if a[1] ~= b[1] then return a[1] < b[1] end
    return a[2] < b[2]
  end)
  return list
end

function T.corpus_present()
  H.eq(18, #snapgene_files())
end

--- .dna files match Biopython feature-for-feature; .prot/.rna match sgffp.
function T.matches_oracles()
  local features = 0
  for _, path in ipairs(snapgene_files()) do
    local oracle = assert(H.oracle(path), "no oracle for " .. path)
    local res = assert(pio.read(path))
    H.eq("snapgene", res.format, path)
    local rec, want = res.records[1], oracle.records[1]
    H.eq(want.length, rec:len(), path .. " length")
    H.eq(want.sha256, vim.fn.sha256(rec.seq:tostring()), path .. " sequence")
    H.eq(want.topology, rec.topology, path .. " topology")
    if oracle.format == "snapgene_seq" then
      H.eq(want.alphabet, rec.alphabet, path .. " alphabet")
      H.eq(want.nfeatures, #rec.features, path .. " features")
    else
      H.eq("dna", rec.alphabet, path)
      if not want.id:match("^<unknown") then H.eq(want.id, rec.id, path .. " id") end
      if not want.name:match("^<unknown") then H.eq(want.name, rec.name, path .. " name") end
      H.eq(#want.features, #rec.features, path .. " feature count")
      for i, wf in ipairs(want.features) do
        local f = rec.features[i]
        local where = ("%s feature %d %s %s"):format(path, i, f.type, f.location_text)
        H.eq(wf.type, f.type, where)
        -- part order compared as a set: Biopython keeps genomic order for an
        -- origin-wrapping reverse feature; we use biological order (see below)
        local got = vim.tbl_map(function(p) return { p.start, p.stop, p.strand, p.ref } end, f.parts)
        H.eq(sorted_parts(wf.parts), sorted_parts(got), where .. " parts")
        local q = {}
        for _, x in ipairs(f.qualifiers) do q[#q + 1] = { x.key, x.value == true and "" or x.value } end
        -- documented difference: Biopython leaves HTML entities (&lt; &amp;) encoded;
        -- the plugin decodes them. Decode the oracle side before comparing.
        local want_q = vim.tbl_map(function(kv)
          local v = kv[2]:gsub("&lt;", "<"):gsub("&gt;", ">"):gsub("&quot;", '"'):gsub("&amp;", "&")
          return { kv[1], v }
        end, wf.qualifiers)
        H.eq(sorted_quals(want_q), sorted_quals(q), where .. " qualifiers")
        features = features + 1
      end
    end
  end
  H.ok(features >= 60, "only " .. features .. " SnapGene features compared")
end

function T.unmodified_files_rewrite_byte_identical()
  for _, path in ipairs(snapgene_files()) do
    local bytes = H.read(path)
    local res = assert(pio.parse(bytes, path))
    local out = assert(pio.format(res.records, "snapgene", res.file))
    H.ok(out == bytes, path .. " not byte-identical")
  end
end

function T.topology_change_touches_only_the_flag_bit()
  local path = H.fixture("snapgene/sample-d.dna")
  local bytes = H.read(path)
  local res = assert(pio.parse(bytes, path))
  local rec = res.records[1]
  H.eq("linear", rec.topology)
  rec:set_topology("circular")
  local out = assert(pio.format(res.records, "snapgene", res.file))
  H.eq(#bytes, #out)
  local diffs = {}
  for i = 1, #bytes do
    if bytes:byte(i) ~= out:byte(i) then diffs[#diffs + 1] = i end
  end
  H.eq(1, #diffs, "exactly one byte should change")
  H.eq(1, bit.bxor(bytes:byte(diffs[1]), out:byte(diffs[1])), "only bit 0 of the flags byte")
  H.eq("circular", assert(pio.parse(out, path)).records[1].topology)
end

function T.wrapped_reverse_feature_in_biological_order()
  local res = assert(pio.read(H.fixture("snapgene/sgffp_origin_spanning_features.dna")))
  local rec = res.records[1]
  local checked = 0
  for _, f in ipairs(rec.features) do
    local ext = rec:extents(f)
    if #ext == 2 then
      checked = checked + 1
      if f.strand == -1 then
        -- 5' end of a reverse feature that wraps is the part touching base 1
        H.eq(1, f.parts[1].start, f.location_text)
        H.ok(f.location_text:match("^complement%(join%("), f.location_text)
      else
        H.eq(rec:len(), f.parts[1].stop, f.location_text)
      end
    end
  end
  H.ok(checked >= 1, "fixture has no origin-wrapping features")
end

function T.converts_to_genbank_and_fasta()
  for _, path in ipairs(snapgene_files()) do
    local res = assert(pio.read(path))
    local rec = res.records[1]
    local text = assert(pio.format(res.records, "genbank", res.file, "snapgene"))
    local back = assert(gb.parse(text)).records[1]
    H.eq(rec.seq:tostring(), back.seq:tostring(), path)
    H.eq(rec.alphabet, back.alphabet, path .. " alphabet survives GenBank")
    H.eq(#rec.features, #back.features, path)
    for i, f in ipairs(rec.features) do
      H.eq(f.location_text, back.features[i].location_text, path)
      H.eq(f:label(), back.features[i]:label(), path)
    end
    local fa = assert(pio.format(res.records, "fasta", res.file, "snapgene"))
    local fb = assert(pio.parse(fa, "x.fasta")).records[1]
    H.eq(rec.seq:tostring(), fb.seq:tostring(), path)
  end
end

function T.protein_genbank_writes_aa_locus()
  local res = assert(pio.read(H.fixture("snapgene/sgffp_test.prot")))
  local text = assert(pio.format(res.records, "genbank", res.file, "snapgene"))
  local locus = text:match("^[^\n]+")
  H.ok(locus:match(" 51 aa "), locus)
  H.ok(not locus:match("DNA"), locus)
end

function T.refuses_to_write_snapgene_from_other_formats()
  local res = assert(pio.read(H.fixture("ncbi/J01749.gb")))
  local out, err = pio.format(res.records, "snapgene", res.file, "genbank")
  H.eq(nil, out)
  H.ok(err:match("not supported"), err)
end

function T.rejects_non_snapgene_and_truncated_files()
  local _, err = pio.parse("MacVector protein document\0\1\2", "x.prot")
  H.ok(err and err:match("MacVector"), tostring(err))
  local bytes = H.read(H.fixture("snapgene/sample-e.dna"))
  local _, err2 = pio.parse(bytes:sub(1, 200), "cut.dna")
  H.ok(err2 and err2:match("file ends first"), tostring(err2))
end

function T.qualifier_text_cleanup()
  local res = assert(pio.read(H.fixture("snapgene/linebreak_in_qualifier_text.dna")))
  for _, f in ipairs(res.records[1].features) do
    for _, q in ipairs(f.qualifiers) do
      if type(q.value) == "string" then
        H.ok(not q.value:find("[\r\n]"), "newline left in /" .. q.key)
        H.ok(not q.value:find("<%a"), "HTML left in /" .. q.key .. ": " .. q.value)
      end
    end
  end
end

function T.file_colours_and_directionality()
  local res = assert(pio.read(H.fixture("snapgene/sample-d.dna")))
  local by = {}
  for _, f in ipairs(res.records[1].features) do by[f:label()] = f end
  H.eq("#ffffff", by.FeatureA.color)
  H.eq(false, by.FeatureA.no_arrow == true, "directional promoter keeps its arrow")
  H.eq(true, by.FeatureB.no_arrow, "non-directional feature has no arrow")
end

function T.xml_parser()
  local doc = xml.parse([==[<?xml version="1.0"?><!-- c --><A x='1' y="a &amp; b &#65;&#x42;">
    <B z="q>r"/><C>text &lt;here&gt;</C><![CDATA[<raw>]]></A>]==])
  local a = xml.find(doc, "A")
  H.eq("a & b AB", a.attrs.y)
  H.eq("q>r", xml.find(doc, "B").attrs.z)
  H.eq("text <here>", xml.text(xml.find(doc, "C")))
  H.ok(vim.tbl_contains(a.children, "<raw>"))
  H.ok(not pcall(xml.parse, "<a><b></a>"))
end

return T
