-- SnapGene binary files: .dna (DNA), .prot (protein), .rna (RNA). Pure Lua.
--
-- A file is a sequence of packets: 1 byte type, uint32 big-endian length,
-- payload. The first packet is the cookie ("SnapGene" + seq type + versions).
--
--   cookie seq type   1 = DNA, 2 = protein, 7 = RNA
--   0x00 / 0x15 / 0x20  DNA / protein / RNA sequence: flags byte + residues
--                       (flags bit 0 = circular)
--   0x0A features (XML)   0x05 primers (XML)   0x06 notes (XML)
--   other packets (history, enzymes, display settings, traces) are kept
--   verbatim and written back untouched.
--
-- Interpretation follows Biopython's SnapGeneIO (feature segments, wrapped
-- ranges, 0-based primer sites, hybridisation filtering, HTML stripping), with
-- one deliberate difference: an origin-wrapping reverse-strand feature is
-- written complement(join(a..L,1..b)) so its parts are in biological order.
local record = require("plasmid.model.record")
local location = require("plasmid.model.location")
local xml = require("plasmid.io.xml")

local M = {}

M.MAGIC = "\9\0\0\0\14SnapGene"

local SEQ_PACKETS = { [0x00] = "dna", [0x15] = "protein", [0x20] = "rna" }
local COOKIE_TYPES = { [1] = "dna", [2] = "protein", [7] = "rna" }
local MOLECULE = { dna = "DNA", protein = "protein", rna = "RNA" }

local function u32(s, i)
  local a, b, c, d = s:byte(i, i + 3)
  return ((a * 256 + b) * 256 + c) * 256 + d
end

local function u16(s, i)
  local a, b = s:byte(i, i + 1)
  return a * 256 + b
end

local function be32(n)
  return string.char(math.floor(n / 16777216) % 256, math.floor(n / 65536) % 256,
    math.floor(n / 256) % 256, n % 256)
end

function M.is_snapgene(text)
  return text:sub(1, #M.MAGIC) == M.MAGIC
end

--- Split into packets { {type=, data=}, ... }.
function M.packets(text)
  local out, i, n = {}, 1, #text
  while i <= n do
    if i + 4 > n then return nil, ("truncated packet header at byte %d"):format(i - 1) end
    local t = text:byte(i)
    local len = u32(text, i + 1)
    if i + 4 + len > n then
      return nil, ("packet 0x%02x at byte %d claims %d bytes, file ends first"):format(t, i - 1, len)
    end
    out[#out + 1] = { type = t, data = text:sub(i + 5, i + 4 + len) }
    i = i + 5 + len
  end
  return out
end

-- SnapGene text values are HTML. As Biopython's _decode: strip tags, then
-- newlines -> space, trim. Unlike Biopython, HTML entities are decoded too,
-- so a literal "<" or "&" written by us (escaped) reads back unchanged.
local HTML_ENT = { amp = "&", lt = "<", gt = ">", quot = '"', apos = "'", nbsp = " " }
local function clean(s)
  s = s:gsub("<[^>]+>", ""):gsub("\r\n", " "):gsub("[\r\n]", " ")
  s = s:gsub("&(#?x?)(%w+);", function(kind, body)
    if kind == "#" then return xml.utf8_char(tonumber(body) or 63) end
    if kind == "#x" then return xml.utf8_char(tonumber(body, 16) or 63) end
    return HTML_ENT[body] or ("&" .. body .. ";")
  end)
  return (s:gsub("^%s+", ""):gsub("%s+$", ""))
end

local function html_escape(s)
  return (tostring(s):gsub("[&<>]", { ["&"] = "&amp;", ["<"] = "&lt;", [">"] = "&gt;" }))
end

--- Build a location node from SnapGene ranges.
---@param ranges { [1]: integer, [2]: integer }[] 1-based inclusive, start > stop = wraps origin
local function ranges_to_location(ranges, strand, seq_len)
  local nodes = {}
  for _, r in ipairs(ranges) do
    local a, b = r[1], r[2]
    if a > b then
      nodes[#nodes + 1] = { kind = "range", start = a, stop = seq_len }
      nodes[#nodes + 1] = { kind = "range", start = 1, stop = b }
    elseif a == b then
      nodes[#nodes + 1] = { kind = "site", pos = a }
    else
      nodes[#nodes + 1] = { kind = "range", start = a, stop = b }
    end
  end
  local node = #nodes == 1 and nodes[1] or { kind = "join", parts = nodes }
  if strand == -1 then node = { kind = "complement", child = node } end
  return location.format(node)
end

local function extend(dst, src)
  for _, v in ipairs(src) do dst[#dst + 1] = v end
  return dst
end

local function parse_range(s)
  local a, b = s:match("^%s*(%d+)%s*%-%s*(%d+)%s*$")
  return tonumber(a), tonumber(b)
end

local function parse_features(data, seq_len, alphabet)
  local doc = xml.parse(data)
  local feats = {}
  for _, fe in ipairs(xml.find_all(doc, "Feature")) do
    local a = fe.attrs
    local dirn = tonumber(a.directionality or "1") or 1
    local strand = (dirn == 2) and -1 or 1
    local ranges, seg_names, color = {}, {}, nil
    local n_parts = 0
    for _, seg in ipairs(xml.find_all(fe, "Segment")) do
      if (seg.attrs.type or "standard") ~= "gap" then
        local s1, s2 = parse_range(seg.attrs.range or "")
        if s1 then
          ranges[#ranges + 1] = { s1, s2 }
          n_parts = n_parts + 1
          if seg.attrs.name and seg.attrs.name ~= "" then
            seg_names[#seg_names + 1] = { n_parts, clean(seg.attrs.name) }
          end
          color = color or seg.attrs.color
        end
      end
    end
    if #ranges > 0 then
      local quals = {}
      if #seg_names > 0 then
        local items = {}
        if strand == -1 then
          for k = #seg_names, 1, -1 do
            items[#items + 1] = (n_parts - seg_names[k][1] + 1) .. ":" .. seg_names[k][2]
          end
        else
          for _, sn in ipairs(seg_names) do items[#items + 1] = sn[1] .. ":" .. sn[2] end
        end
        quals[#quals + 1] = { key = "parts", value = table.concat(items, ";"), quoted = true }
      end
      local has_label, labels = false, {}
      for _, q in ipairs(xml.find_all(fe, "Q")) do
        local qname = q.attrs.name
        for _, v in ipairs(xml.find_all(q, "V")) do
          local val, quoted
          if v.attrs.text then
            val, quoted = clean(v.attrs.text), true
          elseif v.attrs.predef then
            val, quoted = clean(v.attrs.predef), true
          elseif v.attrs.int then
            val, quoted = tostring(tonumber(v.attrs.int) or v.attrs.int), false
          end
          if val then
            quals[#quals + 1] = { key = qname, value = val, quoted = quoted }
            if qname == "label" then
              has_label = true
              labels[val] = true
            end
          end
        end
      end
      local name = a.name and clean(a.name) or nil
      if name and name ~= "" then
        if not has_label then
          table.insert(quals, 1, { key = "label", value = name, quoted = true })
        elseif not labels[name] then
          quals[#quals + 1] = { key = "name", value = name, quoted = true }
        end
      end
      feats[#feats + 1] = {
        type = a.type or "misc_feature",
        location_text = ranges_to_location(ranges, strand, seq_len),
        qualifiers = quals,
        color = color,
        raw_xml = data:sub(fe.s, fe.e),
        recent_id = tonumber(a.recentID),
        no_arrow = alphabet == "protein" or dirn == 0 or dirn == 3 or a.directionality == nil,
      }
    end
  end
  return feats
end

local function parse_primers(data, seq_len)
  local doc = xml.parse(data)
  local min_len, min_tm = 0, 0
  for _, hp in ipairs(xml.find_all(doc, "HybridizationParams")) do
    min_len = tonumber(hp.attrs.minContinuousMatchLen or "0") or 0
    min_tm = tonumber(hp.attrs.minMeltingTemperature or "0") or 0
  end
  local feats = {}
  for _, pr in ipairs(xml.find_all(doc, "Primer")) do
    local name = pr.attrs.name and clean(pr.attrs.name) or nil
    local seen = {}
    for _, site in ipairs(xml.find_all(pr, "BindingSite")) do
      local s1, s2 = parse_range(site.attrs.location or "")
      if s1 then
        -- primer sites are 0-based inclusive
        local a, b = s1 + 1, s2 + 1
        local strand = site.attrs.boundStrand == "1" and -1 or 1
        local loc = ranges_to_location({ { a, b } }, strand, seq_len)
        local simplified = site.attrs.simplified == "1"
        local annealed = site.attrs.annealedBases
        local tm = tonumber(site.attrs.meltingTemperature)
        local keep = not (simplified and seen[loc])
          and not (annealed and #annealed < min_len)
          and not (tm and tm < min_tm)
        if keep then
          seen[loc] = true
          local quals = {}
          if name and name ~= "" then quals[1] = { key = "label", value = name, quoted = true } end
          if site.attrs.annealedBases then
            -- kept off the qualifier list (Biopython parity); available for display
          end
          feats[#feats + 1] = { type = "primer_bind", location_text = loc, qualifiers = quals, origin = "primer",
            primer_seq = pr.attrs.sequence and pr.attrs.sequence:gsub("%s", "") or nil,
            snapgene_site = { annealed = annealed, tm = tm, strand = strand, a = a, b = b } }
        end
      end
    end
  end
  return feats
end

local MONTHS = { "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" }

local function parse_notes(data)
  local doc = xml.parse(data)
  local function child(tag)
    local t = xml.text(xml.find(doc, tag))
    return t and clean(t) or nil
  end
  local info = {}
  info.division = child("Type") == "Synthetic" and "SYN" or "UNC"
  local d = child("LastModified")
  local y, mo, dd = (d or ""):match("^(%d%d%d%d)%.(%d%d?)%.(%d%d?)$")
  if y then info.date = ("%02d-%s-%s"):format(tonumber(dd), MONTHS[tonumber(mo)], y) end
  info.accession = child("AccessionNumber")
  info.comments = child("Comments")
  info.description = child("Description")
  return info
end

--- Parse SnapGene bytes. Returns { records = {Record}, file = {} } or nil, err.
---@param path string|nil used for the record name when the file has none
function M.parse(text, path)
  if not M.is_snapgene(text) then
    return nil, "not a SnapGene file (no SnapGene cookie)"
  end
  local packets, err = M.packets(text)
  if not packets then return nil, err end
  local cookie = packets[1].data
  local seq_type = u16(cookie, 9)
  local alphabet = COOKIE_TYPES[seq_type]

  local seq, flags, seq_index
  for i, p in ipairs(packets) do
    if SEQ_PACKETS[p.type] then
      if seq then return nil, "SnapGene file contains more than one sequence packet" end
      flags = p.data:byte(1) or 0
      seq = p.data:sub(2)
      seq_index = i
      alphabet = alphabet or SEQ_PACKETS[p.type]
    end
  end
  if not seq then
    return nil, ("SnapGene file has no sequence packet (sequence type %d is not supported)"):format(seq_type)
  end
  alphabet = alphabet or "dna"
  -- SnapGene stores RNA with T and displays U; present it as RNA (written back as T)
  if alphabet == "rna" then seq = seq:gsub("T", "U"):gsub("t", "u") end

  local feats, notes = {}, {}
  local ok, perr = pcall(function()
    for _, p in ipairs(packets) do
      if p.type == 0x0A then
        extend(feats, parse_features(p.data, #seq, alphabet))
      elseif p.type == 0x05 then
        extend(feats, parse_primers(p.data, #seq))
      elseif p.type == 0x06 then
        notes = parse_notes(p.data)
      end
    end
  end)
  if not ok then return nil, "SnapGene XML: " .. tostring(perr) end

  local stem = path and (path:match("([^/\\]+)$") or path):gsub("%.[^.]*$", "") or "snapgene"
  local name = notes.comments and notes.comments:match("^(%S+)") or stem
  local rec = record.Record.new({
    name = name,
    id = notes.accession or name,
    description = notes.comments or notes.description or "",
    molecule_type = MOLECULE[alphabet],
    alphabet = alphabet,
    topology = (flags % 2 == 1) and "circular" or "linear",
    division = notes.division,
    date = notes.date,
    seq = seq,
    features = feats,
    format = "snapgene",
    raw = { packets = packets, seq_index = seq_index, flags = flags, seq_type = seq_type },
  })
  return { records = { rec }, file = {} }
end

local SKIP_QUALIFIERS = {
  parts = true, name = true, ApEinfo_fwdcolor = true, ApEinfo_revcolor = true, ApEinfo_graphicformat = true,
}

--- Segments for a feature in SnapGene form: genomic order, 1-based "a-b",
--- an origin-wrapping pair of parts becomes one "a-b" segment with a > b.
local function segments(f, seq_len)
  local parts = {}
  for _, p in ipairs(f.parts) do
    if not p.ref and p.stop >= p.start then parts[#parts + 1] = { p.start, p.stop } end
  end
  if f.strand == -1 then
    local rev = {}
    for i = #parts, 1, -1 do rev[#rev + 1] = parts[i] end
    parts = rev
  end
  -- merge a part ending at the origin with the following part starting at 1
  local out = {}
  local i = 1
  while i <= #parts do
    local p = parts[i]
    local nxt = parts[i + 1]
    if nxt and p[2] == seq_len and nxt[1] == 1 then
      out[#out + 1] = { p[1], nxt[2] }
      i = i + 2
    else
      out[#out + 1] = p
      i = i + 1
    end
  end
  return out
end

--- Generate a <Feature> element for a (new or edited) feature.
function M.feature_xml(f, seq_len, recent_id, default_color)
  local function e(s) return xml.escape(html_escape(s)) end
  local label = f:label()
  local attrs = { ('recentID="%d"'):format(recent_id), ('name="%s"'):format(e(label)) }
  if not f.no_arrow then attrs[#attrs + 1] = ('directionality="%d"'):format(f.strand == -1 and 2 or 1) end
  attrs[#attrs + 1] = ('type="%s"'):format(e(f.type))
  attrs[#attrs + 1] = 'allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"'
  local out = { "<Feature " .. table.concat(attrs, " ") .. ">" }
  local color = f:file_color() or default_color or "#a6acb3"
  for _, seg in ipairs(segments(f, seq_len)) do
    out[#out + 1] = ('<Segment range="%d-%d" color="%s" type="standard"%s/>'):format(
      seg[1], seg[2], xml.escape(color), f.type == "CDS" and ' translated="1"' or "")
  end
  -- group consecutive values of the same qualifier into one <Q>
  local i, q = 1, f.qualifiers
  while i <= #q do
    local key = q[i].key
    local j = i
    while j + 1 <= #q and q[j + 1].key == key do j = j + 1 end
    if not SKIP_QUALIFIERS[key] and not (key == "label" and i == j and q[i].value == label) then
      local vs = {}
      for k = i, j do
        local v = q[k].value
        if v == true then
          vs[#vs + 1] = '<V text=""/>'
        elseif not q[k].quoted and tostring(v):match("^%-?%d+$") then
          vs[#vs + 1] = ('<V int="%s"/>'):format(v)
        else
          vs[#vs + 1] = ('<V text="%s"/>'):format(e(v))
        end
      end
      out[#out + 1] = ('<Q name="%s">%s</Q>'):format(e(key), table.concat(vs))
    end
    i = j + 1
  end
  out[#out + 1] = "</Feature>"
  return table.concat(out)
end

--- Rebuild the Features packet: unchanged features verbatim, others generated.
local function features_packet(rec, old)
  local prolog, suffix = '<?xml version="1.0"?><Features nextValidID="0">', "</Features>"
  if old then
    local first = old:find("<Feature[%s>]")
    local last = old:find("</Features>%s*$")
    if first then
      prolog = old:sub(1, first - 1)
    elseif last then
      prolog = old:sub(1, last - 1)
    end
    if last then suffix = old:sub(last) end
    -- a self-closing <Features .../> root with no children
    if not last and old:match("<Features[^>]*/>") then
      prolog = old:gsub("<Features([^>]*)/>.*$", "<Features%1>")
      suffix = "</Features>"
    end
  end
  local max_id = -1
  for _, f in ipairs(rec.features) do
    if f.recent_id and f.recent_id > max_id then max_id = f.recent_id end
  end
  local next_id = math.max(max_id + 1, tonumber(prolog:match('nextValidID="(%d+)"') or 0) or 0)
  local colors = require("plasmid.config").options.colors or {}
  local body = {}
  for _, f in ipairs(rec.features) do
    if f.origin ~= "primer" then
      if f.raw_xml then
        body[#body + 1] = f.raw_xml
      else
        if not f.recent_id then
          f.recent_id = next_id
          next_id = next_id + 1
        end
        body[#body + 1] = M.feature_xml(f, rec:len(), f.recent_id, colors[f.type])
      end
    end
  end
  if prolog:match('nextValidID="%d+"') then
    prolog = prolog:gsub('nextValidID="%d+"', ('nextValidID="%d"'):format(next_id))
  end
  return prolog .. table.concat(body) .. suffix
end

--- Serialise. Only files read from SnapGene can be written: every original
--- packet is emitted verbatim except the sequence packet (its flags byte
--- carries the topology) and, when features changed, the Features packet.
function M.format(records)
  if #records ~= 1 then return nil, "a SnapGene file holds exactly one sequence" end
  local rec = records[1]
  local raw = rec.raw or {}
  if not raw.packets then
    return nil, "creating SnapGene files from other formats is not supported yet; save as .gb or .fasta"
  end
  if raw.seq_dirty then
    return nil, "writing an edited sequence to SnapGene is not supported yet; save as .gb"
  end
  if raw.primers_dirty then
    return nil, "edited or deleted primer sites cannot be written to SnapGene yet; save as .gb"
  end
  local out = {}
  local wrote_features = false
  for i, p in ipairs(raw.packets) do
    local data = p.data
    if i == raw.seq_index then
      local f = raw.flags or 0
      f = f - (f % 2) + (rec.topology == "circular" and 1 or 0)
      local residues = rec.seq:tostring()
      if rec.alphabet == "rna" then residues = residues:gsub("U", "T"):gsub("u", "t") end
      data = string.char(f) .. residues
    elseif p.type == 0x0A and raw.features_dirty then
      data = features_packet(rec, p.data)
      wrote_features = true
    end
    out[#out + 1] = string.char(p.type) .. be32(#data) .. data
  end
  if raw.features_dirty and not wrote_features then
    local data = features_packet(rec, nil)
    out[#out + 1] = string.char(0x0A) .. be32(#data) .. data
  end
  return table.concat(out)
end

local PACKET_NAMES = {
  [0x00] = "DNA sequence", [0x01] = "compressed DNA", [0x02] = "enzyme visibilities",
  [0x03] = "restriction enzymes", [0x05] = "primers (XML)", [0x06] = "notes (XML)",
  [0x07] = "history tree (LZMA)", [0x08] = "sequence properties (XML)", [0x09] = "cookie",
  [0x0A] = "features (XML)", [0x0B] = "history node", [0x0D] = "display settings",
  [0x0E] = "custom enzyme sets", [0x10] = "trace", [0x11] = "alignable sequences",
  [0x15] = "protein sequence", [0x1C] = "enzyme visibilities", [0x20] = "RNA sequence",
}

--- Human-readable dump for :DnaRaw: packet table, then XML packets verbatim.
function M.describe(text)
  local packets, err = M.packets(text)
  if not packets then return "unreadable SnapGene file: " .. tostring(err) end
  local out = { "<!-- SnapGene file: " .. #packets .. " packets -->" }
  for i, p in ipairs(packets) do
    out[#out + 1] = ("<!-- %2d  0x%02X  %8d bytes  %s -->"):format(i, p.type, #p.data,
      PACKET_NAMES[p.type] or "unknown")
  end
  for _, p in ipairs(packets) do
    if p.data:sub(1, 1) == "<" then
      out[#out + 1] = ""
      out[#out + 1] = ("<!-- packet 0x%02X: %s -->"):format(p.type, PACKET_NAMES[p.type] or "unknown")
      out[#out + 1] = (p.data:gsub("><", ">\n<"))
    end
  end
  return table.concat(out, "\n")
end

return M
