-- GenBank flat file reader/writer. Pure Lua.
--
-- Fidelity policy: anything not modified is written back verbatim (LOCUS line,
-- header blocks, each feature's original lines, BASE COUNT, file trailer).
-- Modified or newly created pieces are regenerated in canonical NCBI layout.
-- The sequence block is always regenerated; NCBI's own layout is canonical,
-- so unmodified NCBI files still round-trip byte-for-byte.
local record = require("plasmid.model.record")

local M = {}

local QUAL_INDENT = 21
local LINE_MAX = 79
local NO_SPACE_JOIN = { translation = true }

---------------------------------------------------------------------------
-- LOCUS
---------------------------------------------------------------------------

--- Permissive LOCUS parsing: token-based, since tools disagree on columns.
function M.parse_locus(line)
  local toks = {}
  for tok in line:gmatch("%S+") do toks[#toks + 1] = tok end
  local info = { name = toks[2] }
  for i = 3, #toks do
    local t = toks[i]
    local lt = t:lower()
    if lt == "rc" and tonumber(toks[i - 1]) then
      -- WGS/TSA master records: "807 rc" counts contig records, not bases
      info.record_count = tonumber(toks[i - 1])
      info.units = "rc"
      local nxt = toks[i + 1]
      if nxt and nxt:lower() ~= "linear" and nxt:lower() ~= "circular" then info.molecule_type = nxt end
    elseif (lt == "bp" or lt == "aa") and tonumber(toks[i - 1]) then
      info.length = tonumber(toks[i - 1])
      info.units = lt
      local nxt = toks[i + 1]
      if nxt and nxt:lower() ~= "linear" and nxt:lower() ~= "circular" then
        info.molecule_type = nxt
      end
    elseif lt == "linear" or lt == "circular" then
      info.topology = lt
    elseif t:match("^%d%d?%-%a%a%a%-%d%d%d%d$") then
      info.date = t
    elseif t:match("^%u%u%u$") and i > 3 then
      info.division = t
    end
  end
  return info
end

function M.format_locus(rec)
  local protein = rec.alphabet == "protein"
  local units = protein and "aa" or "bp"
  local mol = protein and "" or (rec.molecule_type or "DNA")
  local strandness = ""
  local s, m = mol:match("^(%a%a%-)(.+)$")
  if s then strandness, mol = s, m end
  local name = rec.name or "unnamed"
  local n = rec:len()
  if n == 0 and rec.no_sequence then
    if rec.record_count then units = "rc" end
    n = rec.record_count or rec.declared_length or 0
  end
  local len = tostring(n)
  -- NCBI column layout; long names push the length right instead of truncating.
  -- name starts at column 13; length is right-aligned to end at column 40
  local pad = math.max(1, 28 - #name - #len)
  return ("LOCUS       %s%s%s %s %3s%-6s  %-8s %s %s"):format(
    name, (" "):rep(pad), len, units, strandness, mol,
    rec.topology or "linear", rec.division or "UNA",
    rec.date or os.date("%d-%b-%Y"):upper()
  )
end

---------------------------------------------------------------------------
-- Qualifiers
---------------------------------------------------------------------------

local function parse_qualifier(lines)
  local first = lines[1]:sub(2) -- drop '/'
  local key, rest = first:match("^([^=]+)=(.*)$")
  if not key then
    return { key = first, value = true, quoted = false }
  end
  local sep = NO_SPACE_JOIN[key] and "" or " "
  local buf = { rest }
  for i = 2, #lines do buf[#buf + 1] = lines[i] end
  local raw = table.concat(buf, sep)
  if raw:sub(1, 1) == '"' then
    local inner = raw:sub(2, raw:sub(-1) == '"' and -2 or -1)
    return { key = key, value = (inner:gsub('""', '"')), quoted = true }
  end
  return { key = key, value = raw, quoted = false }
end

---------------------------------------------------------------------------
-- Reader
---------------------------------------------------------------------------

local function split_lines(text)
  local crlf = text:find("\r\n", 1, true) ~= nil
  if crlf then text = text:gsub("\r\n", "\n") end
  local lines = {}
  for line in (text .. "\n"):gmatch("(.-)\n") do lines[#lines + 1] = line end
  -- the final split always yields a trailing "" from the appended \n
  if lines[#lines] == "" then lines[#lines] = nil end
  local ends_nl = text:sub(-1) == "\n"
  return lines, crlf, ends_nl
end

--- Parse GenBank text. Returns { records = Record[], file = {...} } or nil, err.
function M.parse(text)
  local lines, crlf, ends_nl = split_lines(text)
  local records = {}
  local file = { crlf = crlf, ends_nl = ends_nl, prefix = {}, trailer = {} }

  local i, n = 1, #lines
  -- anything before the first LOCUS is kept verbatim
  while i <= n and not lines[i]:match("^LOCUS") do
    file.prefix[#file.prefix + 1] = lines[i]
    i = i + 1
  end
  if i > n then return nil, "no LOCUS line found" end

  while i <= n do
    local line = lines[i]
    if not line:match("^LOCUS") then
      -- between/after records: trailer text (e.g. blank line after final //)
      file.trailer[#file.trailer + 1] = line
      i = i + 1
    else
      if #file.trailer > 0 then
        -- text between records; attach to the previous record so order is kept
        records[#records].raw.after = file.trailer
        file.trailer = {}
      end
      local raw = { locus = line, header = {}, pre_origin = {} }
      local info = M.parse_locus(line)
      local feats, seqbuf = {}, {}
      local state = "header"
      local cur_feat, cur_qual
      local qstate = { quoted = false, n = 0, len = 0, last = "" }
      local definition, accession, version
      local last_key
      i = i + 1

      local function close_qual()
        if cur_qual then
          cur_feat.qualifiers[#cur_feat.qualifiers + 1] = parse_qualifier(cur_qual)
          cur_qual = nil
        end
      end
      local function close_feat()
        close_qual()
        if cur_feat then
          cur_feat.location_text = table.concat(cur_feat.loc_lines)
          cur_feat.loc_lines = nil
          feats[#feats + 1] = cur_feat
          cur_feat = nil
        end
      end

      while i <= n do
        line = lines[i]
        if line:match("^//") then
          close_feat()
          raw.terminator = line
          i = i + 1
          break
        end
        if state == "header" then
          if line:match("^FEATURES") then
            raw.features_line = line
            state = "features"
          elseif line:match("^ORIGIN") then
            raw.origin_line = line
            state = "origin"
          else
            raw.header[#raw.header + 1] = line
            local key, val = line:match("^(%u[%u ]-)%s%s+(.*)$")
            if key then
              last_key = key
              if key == "DEFINITION" then definition = { val }
              elseif key == "ACCESSION" then accession = val:match("%S+")
              elseif key == "VERSION" then version = val:match("%S+")
              end
            elseif line:match("^%s") and last_key == "DEFINITION" and definition then
              definition[#definition + 1] = line:match("^%s*(.-)%s*$")
            else
              last_key = line:match("^(%u+)") or last_key
            end
          end
        elseif state == "features" then
          if line:match("^%S") then
            close_feat()
            if line:match("^ORIGIN") then
              raw.origin_line = line
              state = "origin"
            else
              raw.pre_origin[#raw.pre_origin + 1] = line
              state = "post"
            end
          elseif line:match("^%s*$") then
            -- stray blank line inside the table: keep with the feature
            if cur_feat then cur_feat.raw_lines[#cur_feat.raw_lines + 1] = line end
          else
            local indent = #line:match("^( *)")
            local body = line:match("^%s*(.-)%s*$")
            if indent < QUAL_INDENT and indent >= 1 then
              close_feat()
              local key, loc = body:match("^(%S+)%s*(.*)$")
              cur_feat = { type = key, loc_lines = { loc }, qualifiers = {}, raw_lines = { line } }
            elseif cur_feat then
              cur_feat.raw_lines[#cur_feat.raw_lines + 1] = line
              -- is the current qualifier's quoted value still open? Tracked
              -- incrementally (quote count, length, last character) so long
              -- multi-line values parse in linear time; same rule as
              -- the old whole-value check on the joined value.
              local open = cur_qual and qstate.quoted
                and not (qstate.len >= 2 and qstate.last == '"' and qstate.n % 2 == 0)
              if body:sub(1, 1) == "/" and not open then
                close_qual()
                cur_qual = { body }
                local v = body:match("^/[^=]+=(.*)$")
                qstate.quoted = v ~= nil and v:sub(1, 1) == '"'
                qstate.n = v and select(2, v:gsub('"', "")) or 0
                qstate.len = v and #v or 0
                qstate.last = v and v:sub(-1) or ""
              elseif cur_qual then
                cur_qual[#cur_qual + 1] = body
                qstate.n = qstate.n + select(2, body:gsub('"', ""))
                qstate.len = qstate.len + 1 + #body
                qstate.last = #body > 0 and body:sub(-1) or " "
              else
                cur_feat.loc_lines[#cur_feat.loc_lines + 1] = body
              end
            end
          end
        elseif state == "post" then
          if line:match("^ORIGIN") then
            raw.origin_line = line
            state = "origin"
          else
            raw.pre_origin[#raw.pre_origin + 1] = line
          end
        elseif state == "origin" then
          local residues = line:gsub("[%s%d]", "")
          seqbuf[#seqbuf + 1] = residues
        end
        i = i + 1
      end
      close_feat()

      local alpha = "dna"
      local mt = (info.molecule_type or ""):upper()
      if info.units == "aa" then alpha = "protein" elseif mt:find("RNA") then alpha = "rna" end
      local desc = definition and table.concat(definition, " ") or ""
      desc = desc:gsub("%.$", "")
      local seq = table.concat(seqbuf)
      local rec = record.Record.new({
        name = info.name,
        id = version or accession or info.name,
        description = desc,
        molecule_type = info.molecule_type or (alpha == "protein" and "protein" or nil),
        alphabet = alpha,
        topology = info.topology or "linear",
        division = info.division,
        date = info.date,
        seq = seq,
        features = feats,
        format = "genbank",
        raw = raw,
      })
      rec.declared_length = info.length
      rec.record_count = info.record_count
      -- CON / scaffold records and WGS masters carry no residues, only a
      -- CONTIG join of other records or a WGS range of contig records
      if not raw.origin_line and #seq == 0 then
        for _, l in ipairs(raw.pre_origin or {}) do
          local key = l:match("^(%u[%u_]+)%s")
          if key == "CONTIG" then rec.contig = (rec.contig or "") .. vim.trim(l:sub(13))
          elseif key == "WGS" or key == "WGS_SCAFLD" or key == "TSA" or key == "TLS" then
            rec.wgs = rec.wgs or vim.trim(l:sub(13))
          elseif rec.contig and l:match("^%s") then rec.contig = rec.contig .. vim.trim(l) end
        end
        if rec.contig or rec.wgs or info.record_count then rec.no_sequence = true end
      end
      records[#records + 1] = rec
    end
  end
  return { records = records, file = file }
end

---------------------------------------------------------------------------
-- Writer
---------------------------------------------------------------------------

local function wrap_location(text, first_prefix)
  local out, cur = {}, first_prefix
  local width = LINE_MAX
  local pieces = {}
  for piece in text:gmatch("[^,]+,?") do pieces[#pieces + 1] = piece end
  local indent = (" "):rep(QUAL_INDENT)
  for _, piece in ipairs(pieces) do
    if #cur + #piece > width and #cur > QUAL_INDENT then
      out[#out + 1] = cur
      cur = indent
    end
    while #cur + #piece > width do
      local take = width - #cur
      out[#out + 1] = cur .. piece:sub(1, take)
      piece = piece:sub(take + 1)
      cur = indent
    end
    cur = cur .. piece
  end
  out[#out + 1] = cur
  return out
end

local function wrap_qualifier(text, hard)
  local indent = (" "):rep(QUAL_INDENT)
  local width = LINE_MAX - QUAL_INDENT
  local out = {}
  if hard then
    for k = 1, #text, width do out[#out + 1] = indent .. text:sub(k, k + width - 1) end
    return out
  end
  local cur = ""
  for word in text:gmatch("%S+") do
    if cur == "" then
      cur = word
    elseif #cur + 1 + #word <= width then
      cur = cur .. " " .. word
    else
      out[#out + 1] = cur
      cur = word
    end
    while #cur > width do
      out[#out + 1] = cur:sub(1, width)
      cur = cur:sub(width + 1)
    end
  end
  out[#out + 1] = cur
  for k, l in ipairs(out) do out[k] = indent .. l end
  return out
end

--- Canonical lines for one feature.
function M.format_feature(f)
  local key = f.type
  local prefix = "     " .. key .. (" "):rep(math.max(1, QUAL_INDENT - 5 - #key))
  local loc = f.loc and require("plasmid.model.location").format(f.loc) or f.location_text
  local out = wrap_location(loc, prefix)
  for _, q in ipairs(f.qualifiers) do
    local text
    if q.value == true then
      text = "/" .. q.key
    elseif q.quoted then
      text = "/" .. q.key .. '="' .. tostring(q.value):gsub('"', '""') .. '"'
    else
      text = "/" .. q.key .. "=" .. tostring(q.value)
    end
    for _, l in ipairs(wrap_qualifier(text, NO_SPACE_JOIN[q.key])) do out[#out + 1] = l end
  end
  return out
end

--- Canonical ORIGIN block lines (without the ORIGIN header).
function M.format_sequence(seq)
  local out = {}
  local n = #seq
  for pos = 1, n, 60 do
    local groups = {}
    for g = pos, math.min(pos + 59, n), 10 do
      groups[#groups + 1] = seq:sub(g, math.min(g + 9, n, pos + 59))
    end
    out[#out + 1] = ("%9d %s"):format(pos, table.concat(groups, " "))
  end
  return out
end

local function wrap_text_field(key, text)
  local out = {}
  local lines = wrap_qualifier(text, false)
  for k, l in ipairs(lines) do
    local body = l:sub(QUAL_INDENT + 1)
    if k == 1 then
      out[#out + 1] = key .. (" "):rep(12 - #key) .. body
    else
      out[#out + 1] = (" "):rep(12) .. body
    end
  end
  return out
end

--- Lines for one record.
function M.format_record(rec)
  local raw = rec.raw or {}
  local out = {}
  local function add(l) out[#out + 1] = l end

  add((raw.locus and not raw.locus_dirty) and raw.locus or M.format_locus(rec))
  if raw.header and #raw.header > 0 then
    for _, l in ipairs(raw.header) do add(l) end
  else
    local desc = rec.description ~= "" and rec.description or "."
    if not desc:match("%.$") then desc = desc .. "." end
    for _, l in ipairs(wrap_text_field("DEFINITION", desc)) do add(l) end
    add("ACCESSION   " .. (rec.id or rec.name))
    add("VERSION     " .. (rec.id or rec.name))
    add("KEYWORDS    .")
    add("SOURCE      .")
    add("  ORGANISM  .")
    add("            .")
  end
  add(raw.features_line or "FEATURES             Location/Qualifiers")
  for _, f in ipairs(rec.features) do
    if f.raw_lines then
      for _, l in ipairs(f.raw_lines) do add(l) end
    elseif f.raw_prev and f.quals_verbatim then
      -- moved by a sequence edit: new location line(s), original qualifier lines
      local fresh = M.format_feature(f)
      local nloc = 0
      for _, l in ipairs(fresh) do
        if l:sub(QUAL_INDENT + 1, QUAL_INDENT + 1) == "/" then break end
        nloc = nloc + 1
      end
      for k = 1, nloc do add(fresh[k]) end
      local in_quals = false
      for k, l in ipairs(f.raw_prev) do
        if k > 1 and l:sub(QUAL_INDENT + 1, QUAL_INDENT + 1) == "/" then in_quals = true end
        if in_quals then add(l) end
      end
    else
      for _, l in ipairs(M.format_feature(f)) do add(l) end
    end
  end
  if raw.pre_origin and not raw.seq_dirty then
    for _, l in ipairs(raw.pre_origin) do add(l) end
  end
  -- records without residues (CON, WGS master) have no ORIGIN block to write
  if raw.origin_line or rec:len() > 0 or not rec.no_sequence then
    add(raw.origin_line or "ORIGIN")
    for _, l in ipairs(M.format_sequence(rec.seq:tostring())) do add(l) end
  end
  add(raw.terminator or "//")
  if raw.after then
    for _, l in ipairs(raw.after) do add(l) end
  end
  return out
end

--- Serialise records. `file` is the metadata returned by parse (optional).
function M.format(records, file)
  file = file or {}
  local out = {}
  for _, l in ipairs(file.prefix or {}) do out[#out + 1] = l end
  for _, rec in ipairs(records) do
    for _, l in ipairs(M.format_record(rec)) do out[#out + 1] = l end
  end
  for _, l in ipairs(file.trailer or {}) do out[#out + 1] = l end
  local nl = file.crlf and "\r\n" or "\n"
  local text = table.concat(out, nl)
  if file.ends_nl ~= false then text = text .. nl end
  return text
end

return M
