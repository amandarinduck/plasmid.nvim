-- Format detection and read/write dispatch. Pure Lua (uses io.open).
local M = {}

M.formats = {
  genbank = require("plasmid.io.genbank"),
  fasta = require("plasmid.io.fasta"),
  snapgene = require("plasmid.io.snapgene"),
}

M.extensions = {
  gb = "genbank", gbk = "genbank", gbff = "genbank", genbank = "genbank", gbx = "genbank", ape = "genbank",
  gp = "genbank", gpff = "genbank", genpept = "genbank",
  fa = "fasta", fasta = "fasta", fna = "fasta", ffn = "fasta", fas = "fasta", seq = "fasta", fsa = "fasta",
  faa = "fasta", pep = "fasta", mpfa = "fasta",
  dna = "snapgene", prot = "snapgene", rna = "snapgene",
}

-- Extensions that pin the alphabet regardless of content.
M.protein_extensions = { faa = true, pep = true, mpfa = true, gp = true, gpff = true, genpept = true }

local function ext_of(path)
  local ext = path and path:match("%.([%w]+)$")
  return ext and ext:lower() or nil
end

--- Detect format from content first, extension second.
function M.detect(text, path)
  if M.formats.snapgene.is_snapgene(text) then return "snapgene" end
  local head = text:sub(1, 4096):gsub("^\239\187\191", "") -- strip UTF-8 BOM
  local first = head:match("^%s*(.)")
  if head:match("^%s*LOCUS") or head:match("\nLOCUS ") then return "genbank" end
  if first == ">" or first == ";" then return "fasta" end
  local ext = ext_of(path)
  if ext and M.extensions[ext] then return M.extensions[ext] end
  return nil
end

function M.format_for_path(path)
  local ext = ext_of(path)
  return ext and M.extensions[ext] or nil
end

--- Parse text. Returns { records, file, format } or nil, err.
function M.parse(text, path, fmt)
  if not M.formats.snapgene.is_snapgene(text) then text = text:gsub("^\239\187\191", "") end
  fmt = fmt or M.detect(text, path)
  if not fmt then return nil, "unrecognised sequence format" end
  if fmt == "snapgene" and not M.formats.snapgene.is_snapgene(text) then
    local ext = ext_of(path)
    if ext == "prot" then
      return nil, "not a SnapGene file: only SnapGene .prot files are supported (MacVector .prot is not)"
    end
    return nil, "not a SnapGene file (missing SnapGene header)"
  end
  local opts = {}
  if M.protein_extensions[ext_of(path) or ""] then opts.alphabet = "protein" end
  local res, err = M.formats[fmt].parse(text, path, opts)
  if not res then return nil, err end
  res.format = fmt
  return res
end

function M.read(path)
  local fh, err = io.open(path, "rb")
  if not fh then return nil, err end
  local text = fh:read("*a")
  fh:close()
  return M.parse(text, path)
end

--- Serialise to text (or bytes). Returns text, or nil, err when the target
--- format cannot represent the records. When converting between formats the
--- source file's layout metadata is dropped (it describes the other format).
function M.format(records, fmt, file, source_fmt)
  if source_fmt and source_fmt ~= fmt then
    file = { crlf = file and file.crlf }
    local copies = {}
    for i, r in ipairs(records) do
      copies[i] = setmetatable({ raw = {} }, { __index = r })
    end
    records = copies
  end
  return M.formats[fmt].format(records, file)
end

function M.write(path, records, fmt, file, source_fmt)
  local text, ferr = M.format(records, fmt, file, source_fmt)
  if not text then return nil, ferr end
  local fh, err = io.open(path, "wb")
  if not fh then return nil, err end
  fh:write(text)
  fh:close()
  return true
end

return M
