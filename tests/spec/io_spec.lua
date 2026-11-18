local H = require("helpers")
local pio = require("plasmid.io")
local gb = require("plasmid.io.genbank")
local fasta = require("plasmid.io.fasta")
local translate = require("plasmid.bio.translate")

local T = {}

local function all_seq_files()
  local out = H.fixtures("**/*.gb")
  vim.list_extend(out, H.fixtures("**/*.fasta"))
  vim.list_extend(out, H.fixtures("**/*.faa"))
  return out
end

--- Unmodified files must be written back byte-for-byte.
function T.lossless_roundtrip_bytes()
  for _, path in ipairs(all_seq_files()) do
    local text = H.read(path)
    local res = assert(pio.parse(text, path))
    local out = pio.format(res.records, res.format, res.file)
    H.ok(out == text, path .. " not byte-identical: " .. tostring(H.first_diff(text, out)))
  end
end

--- Records must match Biopython's interpretation exactly.
function T.matches_biopython_oracle()
  local checked = 0
  for _, path in ipairs(all_seq_files()) do
    local oracle = H.oracle(path)
    if oracle then
      local res = assert(pio.parse(H.read(path), path))
      H.eq(#oracle.records, #res.records, path .. " record count")
      for ri, want in ipairs(oracle.records) do
        local got = res.records[ri]
        local where = path .. " #" .. ri
        H.eq(want.length, got:len(), where .. " length")
        -- Biopython upper-cases GenBank residues; FASTA residues are compared as-is
        local residues = got.seq:tostring()
        if oracle.format == "genbank" then residues = residues:upper() end
        H.eq(want.sha256, vim.fn.sha256(residues), where .. " sequence")
        H.eq(want.id, got.id, where .. " id")
        if oracle.format == "genbank" then
          H.eq(want.name, got.name, where .. " name")
          H.eq(want.topology, got.topology, where .. " topology")
          H.eq(want.description, got.description, where .. " description")
          H.eq(#want.features, #got.features, where .. " feature count")
          for fi, wf in ipairs(want.features) do
            local gf = got.features[fi]
            local fw = ("%s feature %d (%s %s)"):format(where, fi, gf.type, gf.location_text)
            H.eq(wf.type, gf.type, fw)
            if wf.parts ~= nil then
              local parts = vim.tbl_map(function(p)
                return { p.start, p.stop, p.strand, p.ref }
              end, gf.parts)
              H.eq(wf.parts, parts, fw .. " parts")
            end
            local quals = {}
            for _, q in ipairs(gf.qualifiers) do
              quals[#quals + 1] = { q.key, q.value == true and "" or q.value }
            end
            table.sort(quals, function(a, b)
              if a[1] ~= b[1] then return a[1] < b[1] end
              return a[2] < b[2]
            end)
            H.eq(wf.qualifiers, quals, fw .. " qualifiers")
            checked = checked + 1
          end
        else
          H.eq(want.description, got.id .. (got.description ~= "" and (" " .. got.description) or ""))
        end
      end
    end
  end
  H.eq(H.raw_feature_count(), checked, "features compared with the oracle")
end

--- Regenerating every feature canonically must re-parse to the same model.
function T.canonical_writer_is_semantically_stable()
  for _, path in ipairs(H.fixtures("**/*.gb")) do
    local res = assert(gb.parse(H.read(path)))
    for _, rec in ipairs(res.records) do
      rec.raw.locus_dirty = true
      for _, f in ipairs(rec.features) do f:mark_dirty() end
    end
    local text = gb.format(res.records, res.file)
    local again = assert(gb.parse(text))
    for ri, rec in ipairs(res.records) do
      local r2 = again.records[ri]
      H.eq(rec.name, r2.name, path)
      H.eq(rec.topology, r2.topology, path)
      H.eq(rec.seq:tostring(), r2.seq:tostring(), path)
      for fi, f in ipairs(rec.features) do
        local f2 = r2.features[fi]
        H.eq(f.location_text, f2.location_text, path)
        for qi, q in ipairs(f.qualifiers) do
          local q2 = f2.qualifiers[qi]
          -- words longer than a line are hard-wrapped; re-reading inserts a space,
          -- exactly as Biopython does. Compare with whitespace removed for those.
          local a, b = tostring(q.value), tostring(q2.value)
          if a ~= b then a, b = a:gsub("%s", ""), b:gsub("%s", "") end
          H.eq(a, b, ("%s feature %d /%s"):format(path, fi, q.key))
        end
      end
    end
    for line in text:gmatch("[^\n]+") do
      H.ok(#line <= 80, path .. " canonical line too long: " .. line)
    end
  end
end

function T.canonical_locus_matches_ncbi_columns()
  for _, path in ipairs(H.fixtures("ncbi/*.gb")) do
    local rec = gb.parse(H.read(path)).records[1]
    if #rec.name <= 16 and #tostring(rec:len()) <= 11 then
      H.eq(rec.raw.locus, gb.format_locus(rec), path)
    end
  end
end

function T.locus_parsing_is_permissive()
  local cases = {
    { "LOCUS       pET-28a(+)     5369 bp ds-DNA     circular SYN 01-JAN-1980",
      { name = "pET-28a(+)", length = 5369, topology = "circular", molecule_type = "ds-DNA" } },
    { "LOCUS       Exported 3000 bp DNA circular", { name = "Exported", length = 3000, topology = "circular" } },
    { "LOCUS       x 12 bp DNA", { name = "x", length = 12 } },
  }
  for _, c in ipairs(cases) do
    local info = gb.parse_locus(c[1])
    for k, v in pairs(c[2]) do H.eq(v, info[k], c[1] .. " " .. k) end
  end
end

function T.quoted_values_and_flags()
  local res = assert(gb.parse(H.read(H.fixture("synthetic/grammar.gb"))))
  local rec = res.records[1]
  local primer, flagged
  for _, f in ipairs(rec.features) do
    if f:label() == "revPrimer" then primer = f end
    if f:label() == "flagged" then flagged = f end
  end
  H.eq('a quoted "inner" phrase and a /slash at line start /looks_like_a_qualifier but is not', primer:get("note"))
  H.eq(2, #primer.qualifiers)
  H.eq(true, flagged:get("pseudo"))
  H.eq("1", rec.features[3]:get("codon_start"))
end

--- Every /translation in the corpus must be reproduced from the sequence.
function T.translations_match_annotations()
  local checked = 0
  for _, path in ipairs(H.fixtures("ncbi/*.gb")) do
    for _, rec in ipairs(gb.parse(H.read(path)).records) do
      for _, f in ipairs(rec.features) do
        local want = f:get("translation")
        if f.type == "CDS" and want and not f:get("transl_except") and not f:get("exception")
          and not require("plasmid.model.location").has_remote(f.parts) then
          local nt = rec:extract(f)
          local cs = tonumber(f:get("codon_start") or 1)
          nt = nt:sub(cs)
          local five_partial = (f.strand == -1 and f.parts[1].stop_fuzzy) or (f.strand ~= -1 and f.parts[1].start_fuzzy)
          local got = translate.translate(nt, {
            table = f:get("transl_table"), cds_start = not five_partial and cs == 1, to_stop = false,
          }):gsub("%*$", "")
          H.eq(want, got, ("%s %s %s"):format(path, f:label(), f.location_text))
          checked = checked + 1
        end
      end
    end
  end
  H.ok(checked >= 100, "only " .. checked .. " translations checked")
end

function T.fasta_soft_masking_crlf_and_width()
  local text = ">seq1 desc here [topology=circular]\r\nACGTacgtNN\r\nAC-*\r\n>seq2\r\nGGG\r\n"
  local res = assert(fasta.parse(text))
  H.eq(2, #res.records)
  H.eq("ACGTacgtNNAC-*", res.records[1].seq:tostring())
  H.eq("circular", res.records[1].topology)
  H.eq("desc here [topology=circular]", res.records[1].description)
  H.eq(text, fasta.format(res.records, res.file))
end

function T.conversion_genbank_to_fasta_and_back()
  local res = assert(pio.read(H.fixture("ncbi/L09137.gb")))
  local fa = pio.format(res.records, "fasta", res.file, "genbank")
  H.ok(fa:match("^>L09137%.2 "), fa:sub(1, 40))
  local back = assert(pio.parse(fa, "x.fasta"))
  H.eq(res.records[1].seq:tostring(), back.records[1].seq:tostring())
  local gbtext = pio.format(back.records, "genbank", back.file, "fasta")
  local again = assert(gb.parse(gbtext))
  H.eq(res.records[1].seq:tostring(), again.records[1].seq:tostring())
  H.eq(res.records[1]:len(), again.records[1]:len())
end

function T.detect_by_content_then_extension()
  H.eq("genbank", pio.detect("LOCUS  x", "a.fasta"))
  H.eq("fasta", pio.detect(">x\nAC", "a.gb"))
  H.eq("genbank", pio.detect("", "a.gbk"))
  H.eq(nil, pio.detect("hello", "a.txt"))
end

return T
