-- Opening NCBI records from inside Neovim: :DnaFetch, :DnaSearch, :DnaGene,
-- and following a feature's cross-references (gF). The network work lives in
-- plasmid.ncbi; this module turns its results into buffers and pickers.
local M = {}

local function ncbi() return require("plasmid.ncbi") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end
local function cfg() return require("plasmid.config").options.fetch or {} end

M.active = {}   -- handle -> description, for :DnaFetch cancel

local function track(handle, what)
  if handle then M.active[handle] = what end
  return handle
end
local function done(handle) if handle then M.active[handle] = nil end end

--- Cancel every request in flight; returns how many were cancelled.
function M.cancel_all()
  local n = 0
  for h in pairs(M.active) do
    pcall(function() h:cancel() end)
    M.active[h] = nil
    n = n + 1
  end
  return n
end

local function fmt_len(n, db)
  if not n then return "?" end
  local unit = db == "protein" and "aa" or "bp"
  if n >= 1e6 then return ("%.2f M%s"):format(n / 1e6, unit) end
  if n >= 1e4 then return ("%.1f k%s"):format(n / 1e3, unit) end
  return ("%d %s"):format(n, unit)
end
M.fmt_len = fmt_len

-- File name for a fetched record: accession.version plus range and strand.
function M.file_stem(res, spec)
  local stem = (res.accession or spec.id):gsub("[^%w%._%-]", "_")
  if spec.from and spec.to then stem = ("%s_%d-%d"):format(stem, spec.from, spec.to) end
  if spec.strand == 2 then stem = stem .. "_rc" end
  return stem
end

--- Turn fetched flat-file text into an open buffer. Returns buffer, record.
function M.open_result(res, spec, src)
  local pio = require("plasmid.io")
  local ext = (res.db or spec.db) == "protein" and "gp" or "gb"
  local parsed, err = pio.parse(res.text, "fetched." .. ext, "genbank")
  if not parsed or not parsed.records or #parsed.records == 0 then
    return nil, "could not read the record NCBI sent: " .. tostring(err)
  end
  local rec = parsed.records[1]
  local stem = M.file_stem(res, spec)
  -- keep the record's own name/accession; the buffer name carries the range
  local name = rec.name
  rec.name = stem
  local msg = ("%s  %s  %s"):format(res.accession or spec.id, fmt_len(rec:len(), res.db or spec.db),
    res.from_cache and "(cached)" or "from NCBI")
  local nb = require("plasmid.digest").open_record(src, rec, msg,
    { dir = cfg().dir or vim.fn.getcwd(), ext = ext, file = parsed.file })
  rec.name = name
  if nb then require("plasmid.buffer").render(nb, { keep_bp = 1 }) end
  return nb, rec
end

--- :DnaFetch[!] {accession}[:from-to][:-]
--- bang: refresh the cache and lift the size limit.
function M.fetch(text, opts, cb)
  opts = opts or {}
  local spec, err = ncbi().parse_accession(text)
  if not spec then return notify(err, vim.log.levels.ERROR) end
  local src = vim.api.nvim_get_current_buf()
  if not require("plasmid.buffer").get(src) then src = nil end
  notify(("fetching %s%s from NCBI…"):format(spec.id,
    spec.from and ((":%d-%d"):format(spec.from, spec.to) .. (spec.strand == 2 and " (minus strand)" or "")) or ""))
  local handle
  handle = ncbi().fetch(spec, { force = opts.force, refresh = opts.refresh }, function(e, res)
    done(handle)
    if e then
      local hint = ""
      if not opts.force and tostring(e):lower():find("limit", 1, true) then
        hint = " (:DnaFetch! " .. text .. " fetches it anyway)"
      end
      notify(tostring(e) .. hint, vim.log.levels.ERROR)
      if cb then cb(e) end
      return
    end
    local nb, rec_or_err = M.open_result(res, spec, src)
    if not nb then notify(rec_or_err, vim.log.levels.ERROR) end
    if cb then cb(nil, nb, rec_or_err, res) end
  end)
  return track(handle, "fetch " .. spec.id)
end

local function item_text(it, db)
  local title = it.title or ""
  if #title > 70 then title = title:sub(1, 67) .. "…" end
  return ("%-16s %11s  %-8s %s"):format(it.accession or it.uid or "?", fmt_len(it.length, db),
    it.topology or "", title)
end
M.item_text = item_text

--- :DnaSearch[!] {query}  (bang: protein database)
function M.search(query, protein, opts, cb)
  opts = opts or {}
  if not query or vim.trim(query) == "" then return notify("give a search term, e.g. :DnaSearch pUC19 complete sequence", vim.log.levels.WARN) end
  local db = protein and "protein" or "nuccore"
  notify(("searching NCBI %s for “%s”…"):format(db, query))
  local handle
  handle = ncbi().search(db, query, { retmax = opts.retmax or cfg().search_results or 20 }, function(e, res)
    done(handle)
    if e then notify(tostring(e), vim.log.levels.ERROR); if cb then cb(e) end return end
    if not res or #res.items == 0 then
      notify(("no %s records match “%s”"):format(db, query), vim.log.levels.WARN)
      if cb then cb(nil, nil) end
      return
    end
    vim.ui.select(res.items, {
      prompt = ("NCBI %s: %d of %d hits for “%s”"):format(db, #res.items, res.count or #res.items, query),
      format_item = function(it) return item_text(it, db) end,
    }, function(choice)
      if not choice then if cb then cb(nil, nil) end return end
      M.fetch((db == "protein" and "protein:" or "nuccore:") .. (choice.accession or choice.uid), opts, function(e2, nb, rec)
        if cb then cb(e2, nb, rec, choice) end
      end)
    end)
  end)
  return track(handle, "search " .. query)
end

--- The region to fetch for a gene: flanks added, strand following the gene
--- unless opts.plus. Returns an accession spec string.
function M.gene_region(g, flank, plus)
  local from = math.max(1, g.start - flank)
  local to = g.stop + flank
  local s = ("%s:%d-%d"):format(g.chr_accession, from, to)
  if g.strand == 2 and not plus then s = s .. ":-" end
  return s
end

--- :DnaGene {symbol} [organism words...] [flank=N] [plus]
function M.gene(args, opts, cb)
  opts = opts or {}
  local symbol, words, flank, plus = nil, {}, cfg().gene_flank or 0, false
  for _, a in ipairs(args or {}) do
    local k, v = a:match("^(%w+)=(.+)$")
    if k == "flank" then flank = tonumber((v:gsub(",", ""))) or flank
    elseif a == "plus" then plus = true
    elseif not symbol then symbol = a
    else words[#words + 1] = a end
  end
  if not symbol then return notify("give a gene symbol, e.g. :DnaGene BRCA1 human", vim.log.levels.WARN) end
  local organism = #words > 0 and table.concat(words, " ") or cfg().organism
  notify(("looking up gene %s%s…"):format(symbol, organism and (" in " .. organism) or ""))
  local handle
  handle = ncbi().gene(symbol, organism, {}, function(e, list)
    done(handle)
    if e then notify(tostring(e), vim.log.levels.ERROR); if cb then cb(e) end return end
    local usable = {}
    for _, g in ipairs(list or {}) do if g.chr_accession and g.start and g.stop then usable[#usable + 1] = g end end
    if #usable == 0 then
      notify(("no annotated genomic location for %s%s"):format(symbol, organism and (" in " .. organism) or ""),
        vim.log.levels.WARN)
      if cb then cb(nil, nil) end
      return
    end
    local function go(g)
      local spec = M.gene_region(g, flank, plus)
      M.fetch(spec, opts, function(e2, nb, rec, res)
        if not e2 and nb then
          notify(("%s (%s): %s %d..%d%s, %d bp%s"):format(g.symbol or symbol, g.organism or "", g.chr_accession,
            g.start, g.stop, g.strand == 2 and " (minus strand)" or "", g.stop - g.start + 1,
            flank > 0 and (" + " .. flank .. " bp flanks") or ""))
        end
        if cb then cb(e2, nb, rec, g) end
      end)
    end
    if #usable == 1 then return go(usable[1]) end
    vim.ui.select(usable, {
      prompt = ("%d genes match %s"):format(#usable, symbol),
      format_item = function(g)
        return ("%-10s %-22s %-14s %s  %s"):format(g.symbol or "", (g.organism or ""):sub(1, 22),
          g.location_text or "", g.chr_accession, (g.description or ""):sub(1, 50))
      end,
    }, function(g) if g then go(g) elseif cb then cb(nil, nil) end end)
  end)
  return track(handle, "gene " .. symbol)
end

---------------------------------------------------------------------------
-- Cross-references from features
---------------------------------------------------------------------------

-- qualifier value -> { kind = "fetch"|"gene", text, label }
local XREF = {
  { "^GeneID:(%d+)$", function(id) return { kind = "gene", text = id .. "[uid]", label = "NCBI Gene " .. id } end },
  { "^UniProtKB/[%w%-]+:([%w_]+)$", function(id) return { kind = "fetch", text = "protein:" .. id, label = "UniProt " .. id } end },
  { "^UniProt:([%w_]+)$", function(id) return { kind = "fetch", text = "protein:" .. id, label = "UniProt " .. id } end },
  { "^GI:(%d+)$", nil },  -- ambiguous database; GIs are retired anyway
}

--- Linked records of a feature: list of { kind, text, label, source }.
function M.links(f)
  local out, seen = {}, {}
  local function add(l, source)
    if l and not seen[l.text] then seen[l.text] = true; l.source = source; out[#out + 1] = l end
  end
  for _, v in ipairs(f:get_all("protein_id")) do
    add({ kind = "fetch", text = "protein:" .. v, label = v .. " (protein)" }, "protein_id")
  end
  for _, v in ipairs(f:get_all("db_xref")) do
    for _, x in ipairs(XREF) do
      local id = tostring(v):match(x[1])
      if id then if x[2] then add(x[2](id), "db_xref") end break end
    end
  end
  for _, v in ipairs(f:get_all("inference")) do
    -- e.g. "similar to AA sequence:UniProtKB:P0A7Y4" or "...:RefSeq:NP_414542.1"
    local db, id = tostring(v):match("sequence:([%w/%-]+):([%w_%.]+)")
    if db and db:find("UniProt") then
      add({ kind = "fetch", text = "protein:" .. id, label = id .. " (protein, inference)" }, "inference")
    elseif db == "RefSeq" or db == "INSD" then
      local acc = ncbi().parse_accession(id)
      if acc then
        add({ kind = "fetch", text = acc.db .. ":" .. id, label = id .. " (" .. (acc.db == "protein" and "protein" or "nucleotide") .. ", inference)" }, "inference")
      end
    end
  end
  return out
end

--- gF: open a record linked from the feature under the cursor.
function M.follow(buf, fid, opts, cb)
  local B = require("plasmid.buffer")
  local st = B.get(buf)
  if not st then return end
  if not fid then _, fid = B.cursor_info(vim.fn.bufwinid(buf)) end
  local rec = B.record(st)
  local f = fid and rec:feature_by_id(fid)
  if not f then return notify("no feature under the cursor", vim.log.levels.WARN) end
  local links = M.links(f)
  if #links == 0 then return notify(("%s has no NCBI cross-references"):format(f:label()), vim.log.levels.WARN) end
  local function go(l)
    if l.kind == "gene" then return M.gene({ l.text }, opts, cb) end
    return M.fetch(l.text, opts, cb)
  end
  if #links == 1 then return go(links[1]) end
  vim.ui.select(links, {
    prompt = ("Records linked from %s"):format(f:label()),
    format_item = function(l) return ("%-34s %s"):format(l.label, l.source) end,
  }, function(l) if l then go(l) end end)
end

--- Lines for the feature float.
function M.link_lines(f)
  local links = M.links(f)
  if #links == 0 then return {} end
  local names = {}
  for _, l in ipairs(links) do names[#names + 1] = l.label end
  local key = (require("plasmid.config").options.keymaps or {}).follow_link or "gF"
  return { ("links     %s  (%s opens)"):format(table.concat(names, ", "), key) }
end

return M
