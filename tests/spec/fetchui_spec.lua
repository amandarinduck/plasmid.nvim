-- The fetch UI against a stand-in NCBI client (the client itself is tested
-- on recorded NCBI responses in ncbi_spec.lua).
local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")

local calls = {}
local stub = {}
local TEXT = {
  L09137 = H.read(H.fixture("ncbi/L09137.gb")),
  J01749 = H.read(H.fixture("ncbi/J01749.gb")),
}
local function later(fn) vim.schedule(fn) end
local function handle() return { cancelled = false, cancel = function(self) self.cancelled = true end } end
function stub.config() return { cache_dir = tmp .. "/cache" } end
function stub.parse_accession(text)
  local db = "nuccore"
  local t = text:gsub("^protein:", function() db = "protein"; return "" end):gsub("^nuccore:", "")
  local id, a, b, minus = t:match("^([%w_%.]+):(%d+)%-(%d+)(:?%-?)$")
  if id then return { id = id, db = db, from = tonumber(a), to = tonumber(b), strand = minus == ":-" and 2 or nil } end
  if t:match("^[%w_%.]+$") then return { id = t, db = db } end
  return nil, "not an accession: " .. text
end
function stub.fetch(spec, opts, cb)
  calls[#calls + 1] = { "fetch", spec, opts }
  local h = handle()
  later(function()
    if h.cancelled then return end
    local key = spec.id:gsub("%.%d+$", "")
    if key == "BIG" then return cb("BIG is 9,000,000 bp, over the 5,000,000 bp limit (max_length)") end
    if not TEXT[key] then return cb(("NCBI has no nuccore record %s"):format(spec.id)) end
    cb(nil, { text = TEXT[key], accession = key .. ".1", length = 0, db = spec.db, from_cache = opts.refresh ~= true and #calls > 3 })
  end)
  return h
end
function stub.search(db, term, opts, cb)
  calls[#calls + 1] = { "search", db, term, opts }
  later(function()
    if term == "nothing" then return cb(nil, { count = 0, items = {} }) end
    cb(nil, { count = 57, items = {
      { uid = "1", accession = "L09137.2", title = "pUC19c, complete sequence", length = 2686, topology = "circular" },
      { uid = "2", accession = "J01749.1", title = "pBR322, complete sequence", length = 4361, topology = "circular" },
    } })
  end)
  return handle()
end
function stub.gene(sym, org, opts, cb)
  calls[#calls + 1] = { "gene", sym, org }
  later(function()
    if sym == "NOPE" then return cb(nil, {}) end
    if sym == "672[uid]" or sym == "BRCA1" then
      return cb(nil, { { geneid = "672", symbol = "BRCA1", organism = "Homo sapiens", chr_accession = "J01749.1",
        start = 43044295, stop = 43170327, strand = 2, location_text = "17q21.31", description = "BRCA1 DNA repair associated" } })
    end
    cb(nil, {
      { geneid = "1", symbol = sym, organism = "Homo sapiens", chr_accession = "L09137.2", start = 100, stop = 200, strand = 1 },
      { geneid = "2", symbol = sym, organism = "Mus musculus", chr_accession = "L09137.2", start = 300, stop = 400, strand = 1 },
    })
  end)
  return handle()
end
function stub.cache_list() return {} end
local real_ncbi = package.loaded["plasmid.ncbi"]
package.loaded["plasmid.ncbi"] = stub
local F = require("plasmid.fetch")
require("plasmid.config").options.fetch = { dir = tmp }

local function wait_for(pred, ms) return vim.wait(ms or 2000, pred, 5) end
local function msgs_of(fn)
  local msgs = {}
  local on = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  fn()
  vim.notify = on
  return msgs
end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("fetch_opens_an_unsaved_record", function()
  calls = {}
  local got
  local msgs = msgs_of(function()
    vim.cmd("DnaFetch L09137")
    wait_for(function() return vim.api.nvim_buf_get_name(0):find("L09137") end)
  end)
  local buf = vim.api.nvim_get_current_buf()
  H.ok(vim.api.nvim_buf_get_name(buf):find(tmp .. "/L09137.1.gb", 1, true), vim.api.nvim_buf_get_name(buf))
  local st = buffer.get(buf)
  local rec = buffer.record(st)
  H.eq(2686, rec:len()); H.eq("circular", rec.topology)
  H.ok(st.modified, "not on disk yet")
  H.ok(msgs[#msgs]:find("L09137.1  2686 bp  from NCBI", 1, true), msgs[#msgs])
  -- :w writes NCBI's text back verbatim
  vim.cmd("silent write")
  H.eq(TEXT.L09137, H.read(vim.api.nvim_buf_get_name(buf)))
  -- fetching again does not reuse the saved file: a second name is proposed
  vim.cmd("DnaFetch L09137")
  wait_for(function() return vim.api.nvim_buf_get_name(0):find("L09137.1_2.gb", 1, true) end)
  H.ok(vim.api.nvim_buf_get_name(0):find("L09137.1_2.gb", 1, true), vim.api.nvim_buf_get_name(0))
end)

case("ranges_bang_errors_and_cancel", function()
  calls = {}
  vim.cmd("DnaFetch! J01749:100-500:-")
  wait_for(function() return vim.api.nvim_buf_get_name(0):find("J01749.1_100-500_rc", 1, true) end)
  H.ok(vim.api.nvim_buf_get_name(0):find("J01749.1_100-500_rc.gb", 1, true), vim.api.nvim_buf_get_name(0))
  local spec, opts = calls[1][2], calls[1][3]
  H.eq(100, spec.from); H.eq(500, spec.to); H.eq(2, spec.strand)
  H.ok(opts.force and opts.refresh, "! forces and refreshes")
  -- size refusal suggests the bang form
  local msgs = msgs_of(function()
    vim.cmd("DnaFetch BIG")
    wait_for(function() return false end, 50)
  end)
  H.ok(msgs[#msgs]:find(":DnaFetch! BIG", 1, true), msgs[#msgs])
  msgs = msgs_of(function()
    vim.cmd("DnaFetch ZZ999999")
    wait_for(function() return false end, 50)
  end)
  H.ok(msgs[#msgs]:find("no nuccore record", 1, true), msgs[#msgs])
  msgs = msgs_of(function() vim.cmd("DnaFetch not an accession!") end)
  H.ok(msgs[1]:find("not an accession", 1, true), msgs[1])
  -- cancel: the request in flight never opens a buffer
  local before = vim.api.nvim_get_current_buf()
  vim.cmd("DnaFetch J01749")
  msgs = msgs_of(function() vim.cmd("DnaFetch cancel") end)
  H.ok(msgs[1]:find("cancelled 1 request", 1, true), msgs[1])
  wait_for(function() return false end, 50)
  H.eq(before, vim.api.nvim_get_current_buf())
end)

case("search_picker", function()
  calls = {}
  local orig = vim.ui.select
  local shown
  vim.ui.select = function(items, o, cb) shown = { items = items, o = o }; cb(items[2]) end
  vim.cmd("DnaSearch pBR322 complete")
  wait_for(function() return #calls >= 2 and vim.api.nvim_buf_get_name(0):find("J01749.1.gb", 1, true) end)
  vim.ui.select = orig
  H.ok(vim.api.nvim_buf_get_name(0):find("J01749.1.gb", 1, true), vim.api.nvim_buf_get_name(0))
  H.eq("nuccore", calls[1][2]); H.eq("pBR322 complete", calls[1][3])
  H.ok(shown.o.prompt:find("2 of 57 hits", 1, true), shown.o.prompt)
  local line = shown.o.format_item(shown.items[1])
  H.ok(line:find("L09137.2", 1, true) and line:find("2686 bp", 1, true) and line:find("circular", 1, true), line)
  H.eq("nuccore:J01749.1", calls[2][2].id and ("nuccore:" .. calls[2][2].id) or "")
  -- protein search with !
  local picked = false
  vim.ui.select = function(_, _, cb) picked = true; cb(nil) end
  vim.cmd("DnaSearch! hemoglobin beta")
  wait_for(function() return picked end)
  vim.ui.select = orig
  H.eq("protein", calls[3][2])
  local msgs = msgs_of(function()
    vim.cmd("DnaSearch nothing")
    wait_for(function() return false end, 50)
  end)
  H.ok(msgs[#msgs]:find("no nuccore records match", 1, true), msgs[#msgs])
end)

case("gene_region_strand_and_flank", function()
  -- region arithmetic
  local g = { chr_accession = "NC_000017.11", start = 43044295, stop = 43170327, strand = 2 }
  H.eq("NC_000017.11:43043295-43171327:-", F.gene_region(g, 1000, false))
  H.eq("NC_000017.11:43044295-43170327", F.gene_region(g, 0, true))
  H.eq("NC_1.1:1-150", F.gene_region({ chr_accession = "NC_1.1", start = 20, stop = 100, strand = 1 }, 50, false))
  -- one match: fetched directly, minus-strand gene read 5'->3'
  calls = {}
  vim.cmd("DnaGene BRCA1 human flank=2,000")
  wait_for(function() return #calls >= 2 end)
  H.eq("BRCA1", calls[1][2]); H.eq("human", calls[1][3])
  local spec = calls[2][2]
  H.eq(43042295, spec.from); H.eq(43172327, spec.to); H.eq(2, spec.strand)
  -- several matches: a picker
  calls = {}
  local orig = vim.ui.select
  local prompt
  vim.ui.select = function(items, o, cb) prompt = o.prompt; cb(items[2]) end
  vim.cmd("DnaGene Tp53")
  wait_for(function() return #calls >= 2 and prompt ~= nil end)
  vim.ui.select = orig
  H.ok(prompt:find("2 genes match Tp53", 1, true), prompt)
  H.eq(300, calls[2][2].from)
  local msgs = msgs_of(function()
    vim.cmd("DnaGene NOPE")
    wait_for(function() return false end, 50)
  end)
  H.ok(msgs[#msgs]:find("no annotated genomic location", 1, true), msgs[#msgs])
end)

case("links_from_features", function()
  local f = require("plasmid.model.record").Feature.new({ type = "CDS", location_text = "1..30", qualifiers = {
    { key = "gene", value = "HBB", quoted = true },
    { key = "protein_id", value = "NP_000509.1", quoted = true },
    { key = "db_xref", value = "GeneID:3043", quoted = true },
    { key = "db_xref", value = "UniProtKB/Swiss-Prot:P68871", quoted = true },
    { key = "db_xref", value = "HGNC:HGNC:4827", quoted = true },
    { key = "db_xref", value = "GI:4504349", quoted = true },
    { key = "inference", value = "similar to AA sequence:UniProtKB:P0A7Y4", quoted = true },
  } })
  local links = F.links(f)
  local texts = {}
  for _, l in ipairs(links) do texts[#texts + 1] = l.kind .. " " .. l.text end
  H.eq({ "fetch protein:NP_000509.1", "gene 3043[uid]", "fetch protein:P68871", "fetch protein:P0A7Y4" }, texts)
  H.ok(F.link_lines(f)[1]:find("NP_000509.1 (protein)", 1, true))
  -- gF on a real feature: pBR322's tet CDS has a protein_id
  calls = {}
  local dst = tmp .. "/pbr.gb"
  local fh = io.open(dst, "wb"); fh:write(TEXT.J01749); fh:close()
  vim.cmd("silent edit! " .. dst)
  local buf = vim.api.nvim_get_current_buf()
  local rec = buffer.record(buffer.get(buf))
  local tet
  for _, x in ipairs(rec.features) do if x.type == "CDS" and x:get("gene") == "tet" then tet = x end end
  H.ok(tet and tet:get("protein_id"), "tet has a protein_id")
  F.follow(buf, tet.id)
  wait_for(function() return #calls >= 1 end)
  H.eq("protein", calls[1][2].db)
  H.eq(tet:get("protein_id"), calls[1][2].id)
  -- a feature without links says so
  local msgs = msgs_of(function()
    for _, x in ipairs(rec.features) do
      if x.type == "rep_origin" then F.follow(buf, x.id) break end
    end
  end)
  H.ok(msgs[1]:find("no NCBI cross-references", 1, true), msgs[1])
end)

-- The same commands end to end on the real client, replaying the NCBI
-- responses recorded in tests/fixtures/ncbi_http/ (no network).
case("end_to_end_on_recorded_ncbi_responses", function()
  package.loaded["plasmid.ncbi"] = nil
  local ncbi = require("plasmid.ncbi")
  ncbi.reset()
  ncbi._forget_rate_clock()
  local saved_cfg = require("plasmid.config").options.ncbi
  require("plasmid.config").options.ncbi = vim.tbl_extend("force", saved_cfg or {}, { cache_dir = tmp .. "/ncbi_cache" })
  ncbi.set_transport(ncbi.transports.replay(H.fixture("ncbi_http")))
  local function fetched(cmd, name)
    vim.cmd(cmd)
    wait_for(function() return vim.api.nvim_buf_get_name(0):find(name, 1, true) end, 20000)
    H.ok(vim.api.nvim_buf_get_name(0):find(name, 1, true), cmd .. " -> " .. vim.api.nvim_buf_get_name(0))
    return buffer.record(buffer.get(vim.api.nvim_get_current_buf()))
  end
  -- a range on the minus strand of E. coli: thrA reads forward
  local r = fetched("DnaFetch NC_000913.3:2800-190", "NC_000913.3_190-2800_rc")
  H.eq(2611, r:len())
  local thrA
  for _, f in ipairs(r.features) do if f.type == "CDS" and f:get("gene") == "thrA" then thrA = f end end
  H.ok(thrA and thrA.strand == -1, "thrA on the minus strand of the reverse extract")
  -- the plasmid: byte-identical to the committed fixture once written
  r = fetched("DnaFetch J01749", "J01749.1")
  vim.cmd("silent write")
  H.eq(H.read(H.fixture("ncbi/J01749.gb")), H.read(vim.api.nvim_buf_get_name(0)))
  -- a protein accession opens as a protein record
  r = fetched("DnaFetch NP_000509.1", "NP_000509.1")
  H.ok(r:is_protein()); H.eq(147, r:len())
  H.ok(vim.api.nvim_buf_get_name(0):find("%.gp$"), "saved as GenPept")
  -- gene lookup: thrA region from NCBI Gene coordinates (plus strand)
  local orig = vim.ui.select
  vim.ui.select = function(items, _, cb) cb(items[1]) end
  r = fetched("DnaGene thrA Escherichia coli str. K-12 substr. MG1655", "NC_000913.3_337-2799")
  vim.ui.select = orig
  H.eq(2463, r:len())
  H.eq("ATG", r.seq:sub(1, 3):upper())
  -- the second fetch of the same accession comes from the cache: no request
  local before = #(ncbi.history or {})
  local first_path = vim.api.nvim_buf_get_name(0)
  vim.cmd("DnaFetch J01749")
  wait_for(function() return vim.api.nvim_buf_get_name(0) ~= first_path end, 20000)
  H.ok(vim.api.nvim_buf_get_name(0):find("J01749.1_%d+%.gb$"), vim.api.nvim_buf_get_name(0))
  H.eq(before, #(ncbi.history or {}), "served from the cache")
  -- the size guard explains how to override it
  require("plasmid.config").options.ncbi.max_length = 1000000
  local msgs = msgs_of(function()
    vim.cmd("DnaFetch NC_000913.3")
  end)
  local on = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  wait_for(function() return #msgs >= 2 end, 20000)
  vim.notify = on
  H.ok(msgs[#msgs]:find(":DnaFetch! NC_000913.3", 1, true), msgs[#msgs])
  require("plasmid.config").options.ncbi = saved_cfg
  ncbi.reset()
  if real_ncbi then package.loaded["plasmid.ncbi"] = real_ncbi end
end)

return T
