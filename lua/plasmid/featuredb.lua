-- Feature databases for detection.
--
-- A database is stored as two plain files under stdpath('data')/plasmid/featuredb:
--   <name>.fasta   reference sequences (>id)
--   <name>.tsv     id <TAB> name <TAB> type <TAB> description
-- Nothing is bundled with the plugin; databases are imported by the user:
--   * a pLannotate BLAST_dbs directory (exported with pLannotate's own tools:
--     blastdbcmd for sequences, sqlite for descriptions)
--   * FASTA plus an optional descriptions CSV/TSV in pLannotate's custom-database
--     schema (sseqid,name,type,blurb; common header aliases accepted)
--   * GenBank / SnapGene files or folders: each annotated feature becomes a reference
-- The personal library (:DnaFeatureSave) is always searched as well.
local config = require("plasmid.config")
local util = require("plasmid.util")

local M = {}

function M.dir()
  return config.options.featuredb_dir or (vim.fn.stdpath("data") .. "/plasmid/featuredb")
end

local function read(path)
  local fh = io.open(path, "rb")
  if not fh then return nil end
  local s = fh:read("*a")
  fh:close()
  return s
end

local function write(path, text)
  vim.fn.mkdir(vim.fn.fnamemodify(path, ":h"), "p")
  local fh = assert(io.open(path, "wb"))
  fh:write(text)
  fh:close()
end

---------------------------------------------------------------------------
-- Parsing helpers
---------------------------------------------------------------------------

local function parse_fasta(text)
  local out, cur = {}, nil
  for line in text:gmatch("[^\r\n]+") do
    if line:sub(1, 1) == ">" then
      local id, rest = line:sub(2):match("^(%S+)%s*(.*)$")
      cur = { id = id or ("seq" .. (#out + 1)), header = rest or "", parts = {} }
      out[#out + 1] = cur
    elseif cur then
      cur.parts[#cur.parts + 1] = line:gsub("%s", "")
    end
  end
  for _, r in ipairs(out) do
    r.seq = table.concat(r.parts):upper()
    r.parts = nil
  end
  return out
end

--- Minimal CSV/TSV reader (RFC 4180 quoting for CSV).
function M.parse_table(text, sep)
  sep = sep or (text:match("^[^\n]*\t") and "\t" or ",")
  local rows, row, field, i, inq = {}, {}, {}, 1, false
  local n = #text
  while i <= n do
    local c = text:sub(i, i)
    if inq then
      if c == '"' then
        if text:sub(i + 1, i + 1) == '"' then field[#field + 1] = '"'; i = i + 1 else inq = false end
      else
        field[#field + 1] = c
      end
    elseif c == '"' and sep == "," then
      inq = true
    elseif c == sep then
      row[#row + 1] = table.concat(field); field = {}
    elseif c == "\n" then
      row[#row + 1] = table.concat(field); rows[#rows + 1] = row; row, field = {}, {}
    elseif c ~= "\r" then
      field[#field + 1] = c
    end
    i = i + 1
  end
  if #field > 0 or #row > 0 then row[#row + 1] = table.concat(field); rows[#rows + 1] = row end
  return rows
end

-- header aliases, as accepted by pLannotate's custom databases
local ALIASES = {
  sseqid = "id", id = "id", accession = "id", qseqid = "id",
  name = "name", feature = "name", gene = "name", label = "name",
  type = "type", blurb = "blurb", description = "blurb", desc = "blurb", note = "blurb",
}

--- Descriptions table -> { [id] = { name, type, blurb } }.
function M.parse_descriptions(text)
  local rows = M.parse_table(text)
  if #rows < 2 then return {} end
  local col = {}
  for i, h in ipairs(rows[1]) do
    local k = ALIASES[h:lower():gsub("^%s+", ""):gsub("%s+$", "")]
    if k and not col[k] then col[k] = i end
  end
  if not col.id then error("descriptions table needs an id column (sseqid/id/accession)") end
  local out = {}
  for r = 2, #rows do
    local row = rows[r]
    local id = row[col.id]
    if id and id ~= "" and not out[id] then
      out[id] = {
        name = col.name and row[col.name] ~= "" and row[col.name] or id,
        type = col.type and row[col.type] ~= "" and row[col.type] or "misc_feature",
        blurb = col.blurb and row[col.blurb] or "",
      }
    end
  end
  return out
end

local function esc(s) return (tostring(s or ""):gsub("[\t\r\n]", " ")) end

--- Write a database from refs { id, name, type, blurb, seq }.
function M.save(name, refs)
  local fa, tsv = {}, { "id\tname\ttype\tdescription" }
  local seen = {}
  for _, r in ipairs(refs) do
    local id = r.id
    local k = 1
    while seen[id] do k = k + 1; id = r.id .. "_" .. k end
    seen[id] = true
    fa[#fa + 1] = ">" .. id
    fa[#fa + 1] = r.seq
    tsv[#tsv + 1] = table.concat({ id, esc(r.name), esc(r.type), esc(r.blurb) }, "\t")
  end
  write(M.dir() .. "/" .. name .. ".fasta", table.concat(fa, "\n") .. "\n")
  write(M.dir() .. "/" .. name .. ".tsv", table.concat(tsv, "\n") .. "\n")
  M._cache = {}
  return #refs
end

---------------------------------------------------------------------------
-- Importers
---------------------------------------------------------------------------

--- Each annotated feature of sequence files becomes a reference.
function M.refs_from_records(paths)
  local pio = require("plasmid.io")
  local sel = require("plasmid.select")
  local refs = {}
  for _, p in ipairs(paths) do
    local res = pio.read(p)
    if res then
      for _, rec in ipairs(res.records) do
        if not rec:is_protein() then
          for _, f in ipairs(rec.features) do
            if f.type ~= "source" and f.loc and f:length() >= 8 then
              local s = sel.feature_residues(rec, f, "i")
              if #s >= 8 then
                refs[#refs + 1] = {
                  id = (f:label():gsub("[^%w%._%-]", "_")),
                  name = f:label(), type = f.type,
                  blurb = f:get("note") or f:get("product") or "",
                  seq = s:upper(),
                  source = vim.fn.fnamemodify(p, ":t"),
                }
              end
            end
          end
        end
      end
    end
  end
  return refs
end

local function run(cmd)
  local ok, res = pcall(function() return vim.system(cmd, { text = true }):wait() end)
  if not ok then return nil, tostring(res) end
  if res.code ~= 0 then return nil, (res.stderr ~= "" and res.stderr or ("exit " .. res.code)) end
  return res.stdout
end

--- Export a pLannotate BLAST_dbs/snapgene database using its own tools.
function M.export_plannotate(dir)
  local db = dir .. "/snapgene"
  if vim.fn.filereadable(db .. ".nsq") == 0 and vim.fn.filereadable(dir .. "/BLAST_dbs/snapgene.nsq") == 1 then
    dir = dir .. "/BLAST_dbs"
    db = dir .. "/snapgene"
  end
  if vim.fn.filereadable(db .. ".nsq") == 0 then return nil, "no snapgene BLAST database in " .. dir end
  if vim.fn.executable("blastdbcmd") == 0 then
    return nil, "blastdbcmd (BLAST+) not found; activate the environment where pLannotate is installed, "
      .. "or export manually: blastdbcmd -db " .. db .. " -entry all -outfmt %f > snapgene.fasta"
  end
  local fasta, err = run({ "blastdbcmd", "-db", db, "-entry", "all", "-outfmt", "%f", "-line_length", "1000000" })
  if not fasta then return nil, "blastdbcmd failed: " .. err end
  local meta_csv
  local sqlite = dir .. "/snapgene.db"
  if vim.fn.filereadable(sqlite) == 1 then
    if vim.fn.executable("sqlite3") == 1 then
      meta_csv = run({ "sqlite3", "-csv", "-header", sqlite, "select * from snapgene" })
    end
    if not meta_csv then
      local py = vim.fn.executable("python3") == 1 and "python3" or (vim.fn.executable("python") == 1 and "python")
      if py then
        meta_csv = run({ py, "-c", table.concat({
          "import csv,sqlite3,sys",
          "c=sqlite3.connect(sys.argv[1]);cur=c.execute('select * from snapgene')",
          "w=csv.writer(sys.stdout);w.writerow([d[0] for d in cur.description]);w.writerows(cur)",
        }, "\n"), sqlite })
      end
    end
  elseif vim.fn.filereadable(dir .. "/snapgene.csv") == 1 then
    meta_csv = read(dir .. "/snapgene.csv")
  end
  local meta = meta_csv and M.parse_descriptions(meta_csv) or {}
  local refs = {}
  for _, r in ipairs(parse_fasta(fasta)) do
    local m = meta[r.id] or {}
    refs[#refs + 1] = { id = r.id, name = m.name or r.id, type = m.type or "misc_feature", blurb = m.blurb or "", seq = r.seq }
  end
  if not meta_csv then
    return refs, "descriptions not found (no sqlite3/python); feature names default to ids"
  end
  return refs, nil
end

local SEQ_EXT = { gb = true, gbk = true, gbff = true, genbank = true, ape = true, dna = true }

--- Import from a path. Returns count, name or nil, err.
---@param path string
---@param name string|nil database name (default from the path)
---@param descriptions string|nil descriptions CSV/TSV for a FASTA import
function M.import(path, name, descriptions)
  path = vim.fn.fnamemodify(vim.fn.expand(path), ":p"):gsub("/$", "")
  local refs, note
  if vim.fn.isdirectory(path) == 1 then
    if vim.fn.filereadable(path .. "/snapgene.nsq") == 1 or vim.fn.filereadable(path .. "/BLAST_dbs/snapgene.nsq") == 1 then
      local err
      refs, err = M.export_plannotate(path)
      if not refs then return nil, err end
      note = err
      name = name or "plannotate"
    else
      local files = {}
      for _, f in ipairs(vim.fn.globpath(path, "**/*", false, true)) do
        local ext = f:match("%.(%w+)$")
        if ext and SEQ_EXT[ext:lower()] then files[#files + 1] = f end
      end
      table.sort(files)
      if #files == 0 then return nil, "no GenBank/SnapGene files or BLAST database in " .. path end
      refs = M.refs_from_records(files)
      name = name or vim.fn.fnamemodify(path, ":t")
    end
  else
    local ext = (path:match("%.(%w+)$") or ""):lower()
    if SEQ_EXT[ext] then
      refs = M.refs_from_records({ path })
    else
      local text = read(path)
      if not text then return nil, "cannot read " .. path end
      local seqs = parse_fasta(text)
      if #seqs == 0 then return nil, path .. " is not a FASTA file" end
      local meta = {}
      if not descriptions then
        for _, cand in ipairs({ ".csv", ".tsv" }) do
          local p = path:gsub("%.%w+$", "") .. cand
          if vim.fn.filereadable(p) == 1 then descriptions = p break end
        end
      end
      if descriptions then
        local ok, res = pcall(M.parse_descriptions, assert(read(descriptions), "cannot read " .. descriptions))
        if not ok then return nil, res end
        meta = res
      end
      refs = {}
      for _, r in ipairs(seqs) do
        local m = meta[r.id] or {}
        refs[#refs + 1] = {
          id = r.id, name = m.name or (r.header ~= "" and r.header:match("^(%S+)") or r.id),
          type = m.type or "misc_feature", blurb = m.blurb or r.header, seq = r.seq,
        }
      end
    end
    name = name or vim.fn.fnamemodify(path, ":t:r")
  end
  if #refs == 0 then return nil, "no reference features found in " .. path end
  name = name:gsub("[^%w%._%-]", "_")
  return M.save(name, refs), name, note
end

---------------------------------------------------------------------------
-- Loading
---------------------------------------------------------------------------

M._cache = {}

function M.list()
  local out = {}
  for _, f in ipairs(vim.fn.globpath(M.dir(), "*.fasta", false, true)) do
    out[#out + 1] = vim.fn.fnamemodify(f, ":t:r")
  end
  table.sort(out)
  return out
end

function M.load(name)
  local fa = M.dir() .. "/" .. name .. ".fasta"
  local mtime = vim.fn.getftime(fa)
  local c = M._cache[name]
  if c and c.mtime == mtime then return c.refs end
  local text = read(fa)
  if not text then return nil end
  local meta = {}
  local tsv = read(M.dir() .. "/" .. name .. ".tsv")
  if tsv then meta = M.parse_descriptions(tsv) end
  local refs = {}
  for _, r in ipairs(parse_fasta(text)) do
    local m = meta[r.id] or {}
    refs[#refs + 1] = { id = r.id, name = m.name or r.id, type = m.type or "misc_feature", blurb = m.blurb or "", seq = r.seq, db = name }
  end
  M._cache[name] = { mtime = mtime, refs = refs }
  return refs
end

--- References to search: chosen (or all) imported databases + personal library.
function M.references(names)
  local refs, used = {}, {}
  for _, n in ipairs(names and #names > 0 and names or M.list()) do
    local r = M.load(n)
    if r then
      vim.list_extend(refs, r)
      used[#used + 1] = n
    else
      util.notify("plasmid: no feature database named " .. n, vim.log.levels.WARN)
    end
  end
  local lib = require("plasmid.author").library_path()
  if vim.fn.filereadable(lib) == 1 then
    local mine = M.refs_from_records({ lib })
    for _, r in ipairs(mine) do r.db = "library" end
    if #mine > 0 then
      vim.list_extend(refs, mine)
      used[#used + 1] = "library"
    end
  end
  return refs, used
end

return M
