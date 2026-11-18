-- REBASE import: read the user's own downloaded REBASE EMBOSS-format files
-- (emboss_e.NNN, optionally emboss_r.NNN and emboss_s.NNN) and store the
-- converted enzyme list under stdpath("data"). Nothing from REBASE is bundled
-- with plasmid.nvim and this module never touches the network.
--
-- Output entries use the convention of plasmid.data.enzymes:
--   { name, site, cut5, cut3, suppliers? }
--   cut5  top-strand cut between site offsets cut5 and cut5+1 (0 = just
--         before the first site base; may be negative or exceed #site)
--   cut3  bottom-strand cut, in the same top-strand coordinates
--
-- EMBOSS cut positions have no 0: -1 is the gap just before the first site
-- base, 1 the gap after it. Biopython's Scripts/Restriction/ranacompiler.py
-- converts them as fst5 = c1 (+1 if negative) and fst3 = c2 - size (+1 if
-- negative); with cut3 = size + fst3 that is the same rule for both strands:
--   cut = c < 0 and c + 1 or c
local M = {}

--- Storage file; tests (or users) may set M.path to override.
M.path = nil

local IUPAC = "^[ACGTRYSWKMBDHVN]+$"

local function store_path()
  return M.path or (vim.fn.stdpath("data") .. "/plasmid/rebase_enzymes.json")
end
M.store_path = function() return store_path() end

local function lines(text)
  return vim.split((text or ""):gsub("\r\n?", "\n"), "\n", { plain = true })
end

local function emboss_cut(c)
  if c < 0 then return c + 1 end
  return c
end

--- Parse emboss_e text.
--- Record: name site length ncuts blunt c1 c2 c3 c4 (whitespace separated).
--- Returns list (file order) and a skipped-count table.
function M.parse_emboss_e(text)
  local out, seen = {}, {}
  local skipped = { two_cuts = 0, no_cut = 0, bad_site = 0, duplicate = 0, malformed = 0 }
  for _, line in ipairs(lines(text)) do
    if line:match("^%s*$") or line:match("^#") then goto continue end
    do
      local f = vim.split(vim.trim(line), "%s+")
      local nums = {}
      for i = 3, 9 do nums[i] = tonumber(f[i]) end
      if #f < 9 or not nums[3] or not nums[4] or not nums[6] or not nums[7]
          or not nums[8] or not nums[9] then
        skipped.malformed = skipped.malformed + 1
        goto continue
      end
      local name, site = f[1], f[2]:upper()
      local ncuts, c1, c2, c3, c4 = nums[4], nums[6], nums[7], nums[8], nums[9]
      if ncuts == 0 or c1 == 0 or c2 == 0 then
        skipped.no_cut = skipped.no_cut + 1
      elseif ncuts > 2 or c3 ~= 0 or c4 ~= 0 then
        skipped.two_cuts = skipped.two_cuts + 1
      elseif not site:match(IUPAC) then
        skipped.bad_site = skipped.bad_site + 1
      elseif seen[name] then
        skipped.duplicate = skipped.duplicate + 1
      else
        seen[name] = true
        out[#out + 1] = { name = name, site = site, cut5 = emboss_cut(c1), cut3 = emboss_cut(c2) }
      end
    end
    ::continue::
  end
  return out, skipped
end

--- Parse emboss_r text -> { [name] = "SUPPLIER_CODES" } ("" = none).
--- Each record is: name, organism, isoschizomers, methylation site, source,
--- supplier codes, reference count, references..., then a line "//".
--- Leading '#' lines are the file header.
function M.parse_emboss_r(text)
  local out = {}
  local block = {}
  local in_header = true
  local function flush()
    -- drop blank lines before the name (e.g. a trailing newline after //)
    local i = 1
    while block[i] and block[i]:match("^%s*$") do i = i + 1 end
    local name = block[i] and vim.trim(block[i])
    if name and name ~= "" and out[name] == nil then
      out[name] = vim.trim(block[i + 5] or "")
    end
    block = {}
  end
  for _, line in ipairs(lines(text)) do
    if in_header and line:match("^#") then
      -- header comment
    else
      in_header = false
      if line:match("^//") then flush() else block[#block + 1] = line end
    end
  end
  flush()
  return out
end

--- Parse emboss_s text -> { [code] = "Supplier name" }.
function M.parse_emboss_s(text)
  local out = {}
  for _, line in ipairs(lines(text)) do
    if not line:match("^#") then
      local code, name = line:match("^%s*(%S+)%s+(.-)%s*$")
      if code then out[code] = name end
    end
  end
  return out
end

local function read(path)
  local fh, err = io.open(vim.fn.expand(path), "rb")
  if not fh then return nil, err end
  local s = fh:read("*a")
  fh:close()
  return s
end

local cache

--- Import REBASE EMBOSS files and store the result.
--- Returns list, summary | nil, err.
--- summary = { count, skipped = { two_cuts, no_cut, bad_site, duplicate, malformed },
---             commercial, path }
function M.import(e_path, r_path, s_path)
  if not e_path or e_path == "" then return nil, "no emboss_e file given" end
  local etext, err = read(e_path)
  if not etext then return nil, "cannot read " .. e_path .. ": " .. tostring(err) end
  local list, skipped = M.parse_emboss_e(etext)
  if #list == 0 then return nil, "no usable enzymes in " .. e_path end

  local commercial
  if r_path and r_path ~= "" then
    local rtext, rerr = read(r_path)
    if not rtext then return nil, "cannot read " .. r_path .. ": " .. tostring(rerr) end
    local sup = M.parse_emboss_r(rtext)
    commercial = 0
    for _, e in ipairs(list) do
      e.suppliers = sup[e.name] or ""
      if e.suppliers ~= "" then commercial = commercial + 1 end
    end
  end

  local supplier_names
  if s_path and s_path ~= "" then
    local stext, serr = read(s_path)
    if not stext then return nil, "cannot read " .. s_path .. ": " .. tostring(serr) end
    supplier_names = M.parse_emboss_s(stext)
  end

  table.sort(list, function(a, b) return a.name < b.name end)

  local path = store_path()
  vim.fn.mkdir(vim.fs.dirname(path), "p")
  local doc = {
    format = "plasmid.nvim/rebase",
    version = 1,
    imported = os.date("!%Y-%m-%dT%H:%M:%SZ"),
    source = vim.fs.basename(e_path),
    supplier_names = supplier_names,
    enzymes = list,
  }
  local fh, werr = io.open(path, "wb")
  if not fh then return nil, "cannot write " .. path .. ": " .. tostring(werr) end
  fh:write(vim.json.encode(doc))
  fh:close()

  cache = list
  local summary = { count = #list, skipped = skipped, commercial = commercial or 0, path = path }
  return list, summary
end

local function read_doc()
  local text = read(store_path())
  if not text or text == "" then return nil end
  local ok, doc = pcall(vim.json.decode, text, { luanil = { object = true, array = true } })
  if not ok or type(doc) ~= "table" then return nil end
  return doc
end

--- The imported enzyme list, or nil when nothing has been imported.
--- Cached after the first read; M.reset() drops the cache.
function M.load()
  if cache then return cache end
  local doc = read_doc()
  if not doc then return nil end
  local list = doc.enzymes or doc -- accept a bare list as well
  if type(list) ~= "table" or #list == 0 then return nil end
  cache = list
  return cache
end

--- Supplier code -> name map from the last import (nil when not imported
--- with an emboss_s file).
function M.supplier_names()
  local doc = read_doc()
  return doc and doc.supplier_names or nil
end

function M.reset() cache = nil end

--- Delete the stored import. Returns true if a file was removed.
function M.clear()
  cache = nil
  local path = store_path()
  if vim.fn.filereadable(path) == 1 then return os.remove(path) and true or false end
  return false
end

return M
