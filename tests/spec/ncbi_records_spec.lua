-- Real NCBI records vs the Biopython oracle (tests/oracle/dump_ncbi_real.py).
--
-- Corpus: tests/fixtures/ncbi_real/ (fetched by tests/tools/fetch_ncbi_real.py;
-- why each record was chosen is in manifest.json), plus minimal/ (smallest
-- texts reproducing each bug) and synthetic/ (grammar NCBI no longer serves:
-- old headers, one-of, remote references). Fixtures use .gbk/.gp so that the
-- generic **/*.gb checks in io_spec/location_spec do not pick them up.
--
-- corpus_* cases PASS: they compare everything that already works. Items a
-- bug_* case pins are skipped there, explicitly and counted. Each bug_* case
-- FAILS until its bug is fixed (one case per distinct bug).
local H = require("helpers")
local pio = require("plasmid.io")
local gb = require("plasmid.io.genbank")
local location = require("plasmid.model.location")
local translate = require("plasmid.bio.translate")

require("plasmid").setup({})
vim.o.columns, vim.o.lines = 120, 50

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local BASE = "/tests/fixtures/ncbi_real/"

local oracle_cache
local function oracle()
  if not oracle_cache then
    oracle_cache = vim.json.decode(H.read(H.root .. "/tests/oracle/expected/ncbi_real.json"),
      { luanil = { object = true, array = true } })
  end
  return oracle_cache
end

local function rel(path) return path:sub(#(H.root .. BASE) + 1) end

local function corpus()
  local out = {}
  for _, sub in ipairs({ "", "minimal/", "synthetic/" }) do
    vim.list_extend(out, H.fixtures("ncbi_real/" .. sub .. "*.gbk"))
    vim.list_extend(out, H.fixtures("ncbi_real/" .. sub .. "*.gp"))
  end
  table.sort(out)
  return out
end

local function parsed(path)
  return assert(pio.parse(H.read(path), path))
end

-- Biopython keeps one copy of a repeated value-less qualifier (/pseudo twice
-- -> one); the plugin keeps the file as written, which is the right model.
local function plugin_quals(f)
  local out, seen_flag = {}, {}
  for _, q in ipairs(f.qualifiers) do
    if q.value == true then
      if not seen_flag[q.key] then out[#out + 1] = { q.key, "" } end
      seen_flag[q.key] = true
    else
      out[#out + 1] = { q.key, q.value }
    end
  end
  table.sort(out, function(a, b)
    if a[1] ~= b[1] then return a[1] < b[1] end
    return a[2] < b[2]
  end)
  return out
end

-- Pinned elsewhere: locations the plugin cannot parse yet (bond) or rejects
-- by design (one-of, retired from the INSDC grammar).
local function skip_location(text)
  if text:find("bond(", 1, true) then return "bond" end
  if text:find("one-of(", 1, true) then return "one-of" end
end

local function each_record(fn)
  local o = oracle()
  local nfiles = 0
  for _, path in ipairs(corpus()) do
    local want = o[rel(path)]
    H.ok(want, "no oracle entry for " .. rel(path) .. " (run tests/oracle/dump_ncbi_real.py)")
    H.ok(not want.error, rel(path) .. ": Biopython failed: " .. tostring(want.error))
    local res = parsed(path)
    H.eq(#want.records, #res.records, rel(path) .. " record count")
    for ri, w in ipairs(want.records) do fn(w, res.records[ri], rel(path) .. " #" .. ri, path, res) end
    nfiles = nfiles + 1
  end
  H.ok(nfiles >= 45, "corpus too small: " .. nfiles)
  return nfiles
end

---------------------------------------------------------------------------
-- What already works
---------------------------------------------------------------------------

case("corpus_records_match_biopython", function()
  local n, undefined = 0, 0
  each_record(function(w, got, where)
    H.eq(w.name, got.name, where .. " name")
    H.eq(w.id, got.id, where .. " id")
    H.eq(w.description, got.description, where .. " description")
    H.eq(w.topology, got.topology, where .. " topology")
    H.eq(w.molecule_type, got.molecule_type, where .. " molecule type")
    if w.defined then
      H.eq(w.length, got:len(), where .. " length")
      H.eq(w.sha256, vim.fn.sha256(got.seq:tostring():upper()), where .. " sequence")
    else
      undefined = undefined + 1 -- CON / WGS master: pinned by bug_* cases
    end
    H.eq(#w.features, #got.features, where .. " feature count")
    for fi, wf in ipairs(w.features) do H.eq(wf.type, got.features[fi].type, where .. " feature " .. fi) end
    n = n + 1
  end)
  H.log(("    [ncbi] %d records compared, %d without residues skipped for length/sequence"):format(n, undefined))
end)

case("corpus_locations_match_biopython", function()
  local n, skipped = 0, {}
  each_record(function(w, got, where)
    for fi, wf in ipairs(w.features) do
      local gf = got.features[fi]
      local fw = ("%s feature %d (%s %s)"):format(where, fi, gf.type, gf.location_text)
      local why = skip_location(gf.location_text)
      if why then
        skipped[why] = (skipped[why] or 0) + 1
      else
        H.eq(nil, gf.loc_error, fw .. " parse error")
        local parts = vim.tbl_map(function(p) return { p.start, p.stop, p.strand, p.ref } end, gf.parts)
        H.eq(wf.parts, parts, fw .. " parts")
        n = n + 1
      end
    end
  end)
  H.log(("    [ncbi] %d locations equal; skipped bond=%d one-of=%d"):format(n, skipped.bond or 0, skipped["one-of"] or 0))
  H.ok(n > 3000, "only " .. n .. " locations compared")
end)

case("corpus_qualifiers_match_biopython", function()
  local n, longest = 0, 0
  each_record(function(w, got, where)
    for fi, wf in ipairs(w.features) do
      local gf = got.features[fi]
      H.eq(wf.qualifiers, plugin_quals(gf), ("%s feature %d (%s %s) qualifiers"):format(where, fi, gf.type, gf.location_text))
      for _, q in ipairs(wf.qualifiers) do longest = math.max(longest, #q[2]) end
      n = n + #wf.qualifiers
    end
  end)
  H.log(("    [ncbi] %d qualifier values equal; longest %d chars"):format(n, longest))
  H.ok(longest > 30000, "TTN /translation missing from the corpus")
end)

case("corpus_roundtrip_is_byte_identical", function()
  local n, skipped = 0, 0
  for _, path in ipairs(corpus()) do
    local w = oracle()[rel(path)]
    local all_defined = true
    for _, r in ipairs(w.records) do all_defined = all_defined and r.defined end
    if all_defined then
      local text = H.read(path)
      local res = parsed(path)
      local out = pio.format(res.records, res.format, res.file)
      H.ok(out == text, rel(path) .. " not byte-identical: " .. tostring(H.first_diff(text, out)))
      n = n + 1
    else
      skipped = skipped + 1 -- sequence-less records: bug_con_record_write_adds_origin
    end
  end
  H.log(("    [ncbi] %d files byte-identical, %d sequence-less skipped"):format(n, skipped))
end)

local function three_partial(f)
  local last = f.parts[#f.parts]
  return (f.strand == -1 and last.start_fuzzy) or (f.strand ~= -1 and last.stop_fuzzy)
end

case("corpus_translations_match_annotations", function()
  local checked, gapped, ambiguous = 0, 0, 0
  for _, path in ipairs(corpus()) do
    for _, rec in ipairs(parsed(path).records) do
      for _, f in ipairs(rec.features) do
        local want = f:get("translation")
        if f.type == "CDS" and want and rec:len() > 0 and not f:get("transl_except") and not f:get("exception")
          and #f.parts > 0 and not location.has_remote(f.parts) then
          -- a fuzzy end inside a join marks a gap of unknown size (segmented
          -- sets): the annotated translation then comes from elsewhere
          local internal_fuzzy = false
          for k, p in ipairs(f.parts) do
            if (k > 1 and (p.start_fuzzy or (p.strand == -1 and p.stop_fuzzy)))
              or (k < #f.parts and (p.stop_fuzzy or (p.strand == -1 and p.start_fuzzy))) then
              internal_fuzzy = true
            end
          end
          local cs = tonumber(f:get("codon_start") or 1)
          local nt = rec:extract(f):sub(cs)
          if internal_fuzzy then
            gapped = gapped + 1
          elseif nt:find("[^ACGTacgt]") or (#nt % 3 ~= 0 and three_partial(f)) then
            ambiguous = ambiguous + 1 -- pinned by bug_translate_degenerate_codons
          else
            local five_partial = (f.strand == -1 and f.parts[1].stop_fuzzy) or (f.strand ~= -1 and f.parts[1].start_fuzzy)
            local got = translate.translate(nt, {
              table = f:get("transl_table"), cds_start = not five_partial and cs == 1, to_stop = false,
            }):gsub("%*$", "")
            H.eq(want, got, ("%s %s %s"):format(rel(path), f:label(), f.location_text))
            checked = checked + 1
          end
        end
      end
    end
  end
  H.log(("    [ncbi] %d /translation values reproduced (skipped: %d gap-spanning, %d IUPAC/partial-codon)"):format(
    checked, gapped, ambiguous))
  H.ok(checked >= 350, "only " .. checked .. " translations checked")
end)

case("corpus_opens_and_renders_both_views", function()
  local buffer = require("plasmid.buffer")
  local slow = {}
  for _, path in ipairs(corpus()) do
    local t0 = vim.uv.hrtime()
    vim.cmd("silent edit " .. vim.fn.fnameescape(path))
    local open_ms = (vim.uv.hrtime() - t0) / 1e6
    local buf = vim.api.nvim_get_current_buf()
    local st = buffer.state[buf]
    H.ok(st and st.layout and #st.layout.lines > 0, rel(path) .. " did not render")
    local t1 = vim.uv.hrtime()
    buffer.toggle_view(buf)
    local toggle_ms = (vim.uv.hrtime() - t1) / 1e6
    if open_ms > 250 or toggle_ms > 250 then
      slow[#slow + 1] = ("%s open %.0f ms, other view %.0f ms"):format(rel(path), open_ms, toggle_ms)
    end
    vim.cmd("silent! bwipeout!")
  end
  for _, s in ipairs(slow) do H.log("    [perf] " .. s) end
end)

case("legacy_one_of_degrades_gracefully", function()
  -- one-of() was retired from the INSDC location grammar and NCBI no longer
  -- serves it (U18266's one-of(1888,1901) is now AH005336's misc_feature 1901).
  -- Rejecting it is acceptable; losing the feature or the text is not.
  local path = H.fixture("ncbi_real/synthetic/one_of.gbk")
  local text = H.read(path)
  local res = parsed(path)
  local f = res.records[1].features[1]
  H.eq("one-of(1888,1901)..2200", f.location_text)
  H.ok(f.loc_error, "expected a location error to be reported")
  H.eq("1", f:get("number"))
  H.eq(text, pio.format(res.records, res.format, res.file))
end)

---------------------------------------------------------------------------
-- Bugs (each case fails until fixed)
---------------------------------------------------------------------------

-- GenPept Bond/Het features from UniProtKB and PDB use bond(a,b) / bond(a):
-- 3 of 4 Bond features in P01308 and the Het in 4HHB_A get no parts, so they
-- are not drawn, searchable or navigable.
case("fixed_genpept_bond_location", function()
  local path = H.fixture("ncbi_real/minimal/bond_genpept.gp")
  local want = oracle()["minimal/bond_genpept.gp"].records[1]
  local rec = parsed(path).records[1]
  for fi, wf in ipairs(want.features) do
    local gf = rec.features[fi]
    H.eq(nil, gf.loc_error, gf.location_text .. " parse error")
    local parts = vim.tbl_map(function(p) return { p.start, p.stop, p.strand, p.ref } end, gf.parts)
    H.eq(wf.parts, parts, gf.location_text)
  end
  -- unchanged text is still written back verbatim
  H.eq(H.read(path), pio.format(parsed(path).records, "genbank", parsed(path).file))
end)

-- efetch rettype=gb on a CON/scaffold record returns CONTIG and no ORIGIN.
-- Writing the unmodified record inserts an empty ORIGIN block, so the saved
-- file claims a 0-residue sequence under a 39953 bp LOCUS line.
case("fixed_con_record_write_adds_origin", function()
  local path = H.fixture("ncbi_real/minimal/con_contig.gbk")
  local text = H.read(path)
  local res = parsed(path)
  H.eq(39953, res.records[1].declared_length)
  local out = pio.format(res.records, res.format, res.file)
  H.ok(out == text, "not byte-identical: " .. tostring(H.first_diff(text, out)))
end)

-- A sequence-less CON record keeps its features (1000..2000 on a 0-residue
-- record). The linear view shows "(empty sequence)"; the circular view draws
-- every feature as a multi-turn arc with labels (and takes ~1.5 s for the
-- 45-feature NZ_JCEXQL010000002). Features past the sequence end must not be
-- drawn as if they were on it.
case("fixed_circular_view_draws_features_past_sequence_end", function()
  local config = require("plasmid.config")
  local rec = parsed(H.fixture("ncbi_real/minimal/con_features.gbk")).records[1]
  H.eq(0, rec:len())
  local lin = require("plasmid.view.linear").build(rec, config.options, 100)
  H.ok(not table.concat(lin.lines, "\n"):find("repA", 1, true), "linear view draws repA")
  local circ = require("plasmid.view.circular").build(rec, config.options, 100, 40)
  local drawn = table.concat(circ.lines, "\n")
  H.ok(not drawn:find("repA", 1, true), "circular view draws repA (1000..2000) on a 0 bp record")
end)

-- WGS master records declare "807 rc" (record count, not bp). The LOCUS
-- parser ignores the unit, so declared_length is nil and the record opens as
-- an empty sequence with no warning, unlike CON records.
case("fixed_wgs_master_opens_without_warning", function()
  local path = H.fixture("ncbi_real/minimal/wgs_master.gbk")
  local first = #H.messages + 1
  vim.cmd("silent edit " .. vim.fn.fnameescape(path))
  vim.cmd("silent! bwipeout!")
  local warned = false
  for i = first, #H.messages do
    if (H.messages[i].level or 0) >= vim.log.levels.WARN then warned = true end
  end
  H.ok(warned, "opening a WGS master (no residues, LOCUS '807 rc') gave no warning")
end)

-- Not a parser bug, found by the translation check on two segmented sets: a
-- codon with an ambiguity code at a degenerate position (AH003512: GCN,
-- 4-fold Ala) translates to X, while NCBI's /translation and Biopython give A;
-- and the incomplete last codon of a 3'-partial CDS (AH005336 TMPO: "TC", i.e.
-- TCN = Ser) is dropped, while NCBI translates it. Expected values: Biopython
-- for the codon table, the record's own /translation for the CDSs.
case("fixed_translate_degenerate_codons", function()
  local want = oracle()._codons
  local codons = vim.tbl_keys(want)
  table.sort(codons)
  local bad = {}
  for _, c in ipairs(codons) do
    local got = translate.translate(c, { to_stop = false })
    if got ~= want[c] then bad[#bad + 1] = ("%s: want %s got %s"):format(c, want[c], got) end
  end
  -- the real records: AH003512 CDS join(<456..833,2196..>2398) codon 117 GCN;
  -- AH005336 CDS join(2201..2479,...,10475..>10521) ends in "TC"
  for _, name in ipairs({ "AH003512.gbk", "AH005336.gbk" }) do
    local rec = parsed(H.fixture("ncbi_real/" .. name)).records[1]
    for _, f in ipairs(rec.features) do
      if f.type == "CDS" then
        local nt = rec:extract(f):sub(tonumber(f:get("codon_start") or 1))
        local five = (f.strand == -1 and f.parts[1].stop_fuzzy) or (f.strand ~= -1 and f.parts[1].start_fuzzy)
        local got = translate.translate(nt, { table = f:get("transl_table"), cds_start = not five, to_stop = false })
          :gsub("%*$", "")
        local want_aa = f:get("translation")
        if got ~= want_aa then
          local k = 1
          while k <= #want_aa and got:sub(k, k) == want_aa:sub(k, k) do k = k + 1 end
          bad[#bad + 1] = ("%s %s: residue %d is %q, annotated %q"):format(name, f.location_text, k,
            got:sub(k, k), want_aa:sub(k, k))
        end
      end
    end
  end
  H.eq({}, bad)
end)

-- Multi-line qualifiers are re-joined and re-scanned for balanced quotes on
-- every continuation line: O(lines^2). TTN's 35k-aa /translation costs
-- ~140 ms; doubling the value roughly quadruples the time.
case("fixed_multiline_qualifier_parse_is_quadratic", function()
  local function record(n)
    local aa = ("ACDEFGHIKLMNPQRSTVWY"):rep(math.ceil(n / 20)):sub(1, n)
    local lines = { '                     /translation="' .. aa:sub(1, 44) }
    for i = 45, n, 58 do lines[#lines + 1] = (" "):rep(21) .. aa:sub(i, i + 57) end
    lines[#lines] = lines[#lines] .. '"'
    return "LOCUS       T 60 bp    DNA     linear   SYN 01-JAN-2000\nFEATURES             Location/Qualifiers\n"
      .. "     CDS             1..60\n" .. table.concat(lines, "\n") .. "\nORIGIN\n        1 "
      .. ("acgtacgtac "):rep(6) .. "\n//\n"
  end
  local function best_ms(text)
    local best = math.huge
    for _ = 1, 3 do
      local t0 = vim.uv.hrtime()
      local res = gb.parse(text)
      best = math.min(best, (vim.uv.hrtime() - t0) / 1e6)
      assert(res)
    end
    return best
  end
  local small, big = record(20000), record(80000)
  local a, b = best_ms(small), best_ms(big)
  H.log(("    [perf] /translation 20k aa %.1f ms, 80k aa %.1f ms (x%.1f for 4x input)"):format(a, b, b / a))
  H.ok(b / a < 8, ("parse time grows x%.1f for 4x qualifier length (linear would be ~x4)"):format(b / a))
end)

return T
