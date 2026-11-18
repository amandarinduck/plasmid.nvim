-- FASTA reader/writer. Pure Lua.
--
-- Preserves: header line verbatim, per-record line width, residue case
-- (soft-masking), '*' and '-' characters, CRLF, trailing newline.
local record = require("plasmid.model.record")
local alphabet = require("plasmid.model.alphabet")

local M = {}

local function topology_from_header(desc)
  local d = desc:lower()
  if d:match("%[topology=circular%]") or d:match("topology=circular") or d:match("%(circular%)") then
    return "circular"
  end
  return "linear"
end

--- @param opts { alphabet?: string }|nil force an alphabet (e.g. from a .faa extension)
function M.parse(text, _path, opts)
  opts = opts or {}
  local crlf = text:find("\r\n", 1, true) ~= nil
  if crlf then text = text:gsub("\r\n", "\n") end
  local ends_nl = text:sub(-1) == "\n"
  local records, file = {}, { crlf = crlf, ends_nl = ends_nl, prefix = {} }
  local cur, buf, width, ragged
  local function close()
    if not cur then return end
    local seq = table.concat(buf)
    cur.raw.width = ragged and width or width
    if width and width >= #seq then cur.raw.width = 0 end -- single-line record
    local alpha = opts.alphabet or alphabet.guess(seq)
    local rec = record.Record.new({
      name = cur.id, id = cur.id, description = cur.desc,
      alphabet = alpha,
      molecule_type = ({ dna = "DNA", rna = "RNA", protein = "protein" })[alpha],
      topology = alpha == "protein" and "linear" or topology_from_header(cur.desc),
      seq = seq, format = "fasta", raw = cur.raw,
    })
    records[#records + 1] = rec
  end
  for line in (text .. (ends_nl and "" or "\n")):gmatch("(.-)\n") do
    if line:sub(1, 1) == ">" then
      close()
      local hdr = line:sub(2)
      local id, desc = hdr:match("^(%S*)%s*(.-)%s*$")
      cur = { id = id ~= "" and id or ("record" .. (#records + 1)), desc = desc, raw = { header = line } }
      buf, width, ragged = {}, nil, false
    elseif cur and not line:match("%S") then
      -- blank lines (NCBI separates records with one) are re-emitted after the record
      cur.raw.blank = (cur.raw.blank or 0) + 1
    elseif cur then
      local s = line:gsub("%s", "")
      if #s > 0 then
        if width == nil then width = #s end
        buf[#buf + 1] = s
      end
    elseif line:match("%S") then
      if line:sub(1, 1) == ";" then
        file.prefix[#file.prefix + 1] = line
      else
        return nil, "FASTA: sequence data before first '>' header"
      end
    else
      file.prefix[#file.prefix + 1] = line
    end
  end
  close()
  if #records == 0 then return nil, "FASTA: no records" end
  return { records = records, file = file }
end

function M.format_record(rec, default_width)
  local raw = rec.raw or {}
  local out = {}
  local header = raw.header
  if not header or raw.header_dirty then
    header = ">" .. (rec.id or rec.name)
    if rec.description and rec.description ~= "" then header = header .. " " .. rec.description end
  end
  out[1] = header
  local seq = rec.seq:tostring()
  local w = raw.width
  if w == nil then w = default_width or 60 end
  if w == 0 then
    if #seq > 0 then out[#out + 1] = seq end
  else
    for i = 1, #seq, w do out[#out + 1] = seq:sub(i, i + w - 1) end
  end
  for _ = 1, raw.blank or 0 do out[#out + 1] = "" end
  return out
end

function M.format(records, file)
  file = file or {}
  local out = {}
  for _, l in ipairs(file.prefix or {}) do out[#out + 1] = l end
  for _, rec in ipairs(records) do
    for _, l in ipairs(M.format_record(rec)) do out[#out + 1] = l end
  end
  local nl = file.crlf and "\r\n" or "\n"
  local text = table.concat(out, nl)
  if file.ends_nl ~= false then text = text .. nl end
  return text
end

return M
