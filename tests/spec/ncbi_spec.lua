-- NCBI E-utilities client (lua/plasmid/ncbi.lua).
--
-- Runs offline: recorded responses in tests/fixtures/ncbi_http/ are served by
-- the replay transport; stub transports cover rate limiting, retries, errors
-- and cancel; a local TCP server exercises the real curl transport.
-- PLASMID_LIVE=1 additionally re-runs the fetch/search/gene cases against NCBI.
local H = require("helpers")
local pio = require("plasmid.io")
local ncbi = require("plasmid.ncbi")
local cases = require("ncbi_cases")
local cfgmod = require("plasmid.config")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local HTTP = H.root .. "/tests/fixtures/ncbi_http"
local now = function() return vim.uv.hrtime() / 1e6 end

-- own complement table (independent of the plugin's alphabet module)
local COMP = { A = "T", C = "G", G = "C", T = "A", U = "A", R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
  B = "V", V = "B", D = "H", H = "D", N = "N" }
local function rc(s)
  s = s:upper()
  local t = {}
  for i = #s, 1, -1 do t[#t + 1] = COMP[s:sub(i, i)] or "?" end
  return table.concat(t)
end

--- Fresh client state + ncbi options (api key from the environment masked).
local function setup(opts, keep_clock)
  ncbi.reset()
  if not keep_clock then ncbi._forget_rate_clock() end
  vim.env.NCBI_API_KEY = nil
  cfgmod.options.ncbi = vim.tbl_extend("force", { cache = false, min_interval_ms = 0, retry_base_ms = 1 }, opts or {})
end

--- Replay transport counting calls.
local function replay()
  local inner = ncbi.transports.replay(HTTP)
  local calls = {}
  ncbi.set_transport(function(req, done)
    calls[#calls + 1] = req
    return inner(req, done)
  end)
  return calls
end

--- Stub transport: responder(req, n) -> err, resp ("hang" = never answer).
local function stub(responder)
  local calls = { cancelled = 0 }
  ncbi.set_transport(function(req, done)
    calls[#calls + 1] = { req = req, t = now() }
    local err, resp = responder(req, #calls)
    if err == "hang" then
      return { cancel = function() calls.cancelled = calls.cancelled + 1 end }
    end
    done(err, resp)
  end)
  return calls
end

local function json_resp(tbl, status, headers)
  return nil, { status = status or 200, headers = headers or { ["content-type"] = "application/json" }, body = vim.json.encode(tbl) }
end

local function fake_summary(acc, len)
  return { header = {}, result = { uids = { "1" }, ["1"] = { uid = "1", caption = acc:gsub("%.%d+$", ""),
    accessionversion = acc, title = "test record " .. acc, slen = len, topology = "linear", moltype = "dna",
    organism = "synthetic", updatedate = "2026/01/01" } } }
end

local function fake_gb(acc, seq)
  return ("LOCUS       %s %d bp    DNA     linear   SYN 01-JAN-2026\nDEFINITION  test.\nACCESSION   %s\nVERSION     %s\n"
    .. "FEATURES             Location/Qualifiers\nORIGIN\n        1 %s\n//\n"):format(
    acc:gsub("%.%d+$", ""), #seq, acc:gsub("%.%d+$", ""), acc, seq:lower())
end

--- Stub that answers esummary + efetch for one synthetic record.
local function record_stub(acc, seq)
  return stub(function(req)
    if req.url:find("esummary", 1, true) then return json_resp(fake_summary(acc, #seq)) end
    return nil, { status = 200, headers = {}, body = fake_gb(acc, seq) }
  end)
end

local function run(name, extra)
  local h = cases.run(ncbi, cases.by_name[name], function() end, extra)
  local done, err, res = ncbi.wait(h, 10000)
  H.ok(done, name .. ": no callback")
  return err, res
end

local function parse(text, name)
  local res = assert(pio.parse(text, name .. ".gb", "genbank"))
  return res.records[1], res
end

--------------------------------------------------------------------------------
-- accession parsing
--------------------------------------------------------------------------------

case("parse_accession_accepts", function()
  local want = {
    ["L09137"] = { id = "L09137", db = "nuccore" },
    ["L09137.2"] = { id = "L09137.2", db = "nuccore" },
    ["  NM_000518.5\t"] = { id = "NM_000518.5", db = "nuccore" },
    ["AB000263"] = { id = "AB000263", db = "nuccore" },
    ["l09137"] = { id = "l09137", db = "nuccore" },
    ["NC_000913.3:1000-5000"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000 },
    ["NC_000913.3:5000-1000"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000, strand = 2 },
    ["NC_000913.3:1000-5000:-"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000, strand = 2 },
    ["NC_000913.3:1000-5000:minus"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000, strand = 2 },
    ["NC_000913.3:5000-1000:-"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000, strand = 2 },
    ["NC_000913.3:1000-5000:+"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000, strand = 1 },
    ["NC_000913.3:1000..5000"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000 },
    ["NC_000913.3:5000..1000"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000, strand = 2 },
    ["NC_000913.3:1,000-5,000"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000 },
    ["NC_000913.3: 1,000 .. 5,000 :minus"] = { id = "NC_000913.3", db = "nuccore", from = 1000, to = 5000, strand = 2 },
    ["L09137:-"] = { id = "L09137", db = "nuccore", strand = 2 },
    ["L09137:minus"] = { id = "L09137", db = "nuccore", strand = 2 },
    ["L09137:MINUS"] = { id = "L09137", db = "nuccore", strand = 2 },
    ["NZ_CP009072.1"] = { id = "NZ_CP009072.1", db = "nuccore" },
    ["NZ_JAAAAA010000001.1"] = { id = "NZ_JAAAAA010000001.1", db = "nuccore" },
    ["NW_003315950.2"] = { id = "NW_003315950.2", db = "nuccore" },
    ["XM_005257322.4"] = { id = "XM_005257322.4", db = "nuccore" },
    ["AAAA01000001"] = { id = "AAAA01000001", db = "nuccore" },
    ["CBZL010000100.1"] = { id = "CBZL010000100.1", db = "nuccore" },
    ["LN625376.1"] = { id = "LN625376.1", db = "nuccore" },
    ["20141090"] = { id = "20141090", db = "nuccore" },
    ["NP_000509.1"] = { id = "NP_000509.1", db = "protein" },
    ["XP_011524365.1"] = { id = "XP_011524365.1", db = "protein" },
    ["WP_000184067.1"] = { id = "WP_000184067.1", db = "protein" },
    ["YP_009724390.1"] = { id = "YP_009724390.1", db = "protein" },
    ["AP_000001.1"] = { id = "AP_000001.1", db = "protein" },
    ["AAA12345"] = { id = "AAA12345", db = "protein" },
    ["AAA12345.1"] = { id = "AAA12345.1", db = "protein" },
    ["MBC1234567.1"] = { id = "MBC1234567.1", db = "protein" },
    ["P69905"] = { id = "P69905", db = "protein" },
    ["P69905.2"] = { id = "P69905.2", db = "protein" },
    ["Q9Y6K9"] = { id = "Q9Y6K9", db = "protein" },
    ["A0A024R161"] = { id = "A0A024R161", db = "protein" },
    ["NP_000509.1:1-50"] = { id = "NP_000509.1", db = "protein", from = 1, to = 50 },
    ["protein:L09137"] = { id = "L09137", db = "protein" },
    ["Protein: NP_000509.1"] = { id = "NP_000509.1", db = "protein" },
    ["nuccore:P69905"] = { id = "P69905", db = "nuccore" },
    ["nucleotide:FOO123"] = { id = "FOO123", db = "nuccore" },
    ["nuccore:NC_000913.3:190-2800"] = { id = "NC_000913.3", db = "nuccore", from = 190, to = 2800 },
  }
  local n = 0
  for input, exp in pairs(want) do
    local got, err = ncbi.parse_accession(input)
    H.eq(exp, got, ("parse_accession(%q) %s"):format(input, err or ""))
    n = n + 1
  end
  H.ok(n >= 40)
end)

case("parse_accession_rejects", function()
  local bad = {
    [""] = "empty accession",
    ["   "] = "empty accession",
    ["hello"] = "not an accession",
    ["pUC19"] = "not an accession",
    ["L0913"] = "not an accession",
    ["L091377"] = "not an accession",
    ["XY_123456"] = "unknown RefSeq prefix XY_",
    ["NC_000913.3:abc"] = "cannot read the range",
    ["NC_000913.3:100"] = "cannot read the range",
    ["NC_000913.3:0-100"] = "1-based",
    ["NC_000913.3:5000-1000:+"] = "reversed range",
    ["NP_000509.1:-"] = "reverse strand only applies to nucleotide",
    ["P69905:10-1"] = "reverse strand only applies to nucleotide",
    ["L09137 J01749"] = "not an accession",
    ["L09137;rm -rf"] = "not an accession",
    ["protein:"] = "not an accession",
    ["https://www.ncbi.nlm.nih.gov/nuccore/L09137"] = "not an accession",
  }
  for input, msg in pairs(bad) do
    local got, err = ncbi.parse_accession(input)
    H.ok(got == nil, ("%q should be rejected, got %s"):format(input, vim.inspect(got)))
    H.ok(err and err:find(msg, 1, true), ("%q: error %q lacks %q"):format(input, tostring(err), msg))
  end
  local _, err = ncbi.parse_accession("hello")
  H.ok(err:find("L09137", 1, true) and err:find("NC_000913.3:190-2800", 1, true), "error lists accepted forms")
  H.eq(nil, (ncbi.parse_accession(42)))
end)

--------------------------------------------------------------------------------
-- fetch on replayed responses
--------------------------------------------------------------------------------

case("fetch_overlaps_existing_fixtures", function()
  setup()
  replay()
  for _, id in ipairs({ "L09137", "J01749", "NC_001422" }) do
    local err, res = run(id)
    H.eq(nil, err, id)
    H.eq(false, res.from_cache)
    H.eq("nuccore", res.db)
    local got = parse(res.text, id)
    local want = assert(pio.read(H.fixture("ncbi/" .. id .. ".gb"))).records[1]
    H.eq(want:len(), got:len(), id .. " length")
    H.eq(want.seq:tostring():upper(), got.seq:tostring():upper(), id .. " sequence")
    H.eq(#want.features, #got.features, id .. " feature count")
    H.eq(want.topology, got.topology, id .. " topology")
    H.eq(res.length, got:len(), id .. " summary length = record length")
    H.ok(res.accession:match("^" .. id .. "%.%d+$"), id .. " accession.version " .. tostring(res.accession))
    -- the committed fixtures were downloaded from the same NCBI records
    local same = H.read(H.fixture("ncbi/" .. id .. ".gb")) == res.text
    H.ok(same, id .. " fetched text differs from the fixture: " .. tostring(H.first_diff(H.read(H.fixture("ncbi/" .. id .. ".gb")), res.text)))
  end
end)

--- Plugin parser vs Biopython on every recorded flat file (oracle JSON made by
--- tests/oracle/dump_ncbi_http.py).
case("fetched_records_match_biopython", function()
  local files = vim.fn.globpath(H.root .. "/tests/oracle/expected", "ncbi_http__*.json", false, true)
  table.sort(files)
  H.ok(#files >= 10, "oracle files: " .. #files)
  local checked = 0
  for _, file in ipairs(files) do
    local oracle = vim.json.decode(H.read(file), { luanil = { object = true, array = true } })
    local text = H.read(HTTP .. "/" .. oracle.key .. ".body")
    local res = assert(pio.parse(text, oracle.key .. ".gb", "genbank"))
    local where0 = vim.fn.fnamemodify(file, ":t")
    H.eq(#oracle.records, #res.records, where0 .. " record count")
    for ri, want in ipairs(oracle.records) do
      local got = res.records[ri]
      local where = where0 .. " #" .. ri
      H.eq(want.length, got:len(), where .. " length")
      H.eq(want.sha256, vim.fn.sha256(got.seq:tostring():upper()), where .. " sequence")
      H.eq(want.id, got.id, where .. " id")
      H.eq(want.name, got.name, where .. " name")
      H.eq(want.topology, got.topology, where .. " topology")
      H.eq(want.description, got.description, where .. " description")
      H.eq(#want.features, #got.features, where .. " feature count")
      for fi, wf in ipairs(want.features) do
        local gf = got.features[fi]
        local fw = ("%s feature %d (%s %s)"):format(where, fi, gf.type, gf.location_text)
        H.eq(wf.type, gf.type, fw)
        if wf.parts ~= nil then
          H.eq(wf.parts, vim.tbl_map(function(p) return { p.start, p.stop, p.strand, p.ref } end, gf.parts), fw .. " parts")
        end
        local quals = {}
        for _, q in ipairs(gf.qualifiers) do quals[#quals + 1] = { q.key, q.value == true and "" or q.value } end
        table.sort(quals, function(a, b) if a[1] ~= b[1] then return a[1] < b[1] end return a[2] < b[2] end)
        H.eq(wf.qualifiers, quals, fw .. " qualifiers")
        checked = checked + 1
      end
    end
    -- unmodified records are written back byte-for-byte
    local out = pio.format(res.records, res.format, res.file)
    H.ok(out == text, where0 .. " not byte-identical: " .. tostring(H.first_diff(text, out)))
  end
  H.ok(checked > 100, "features compared: " .. checked)
end)

case("mrna_cds_translates_to_fetched_protein", function()
  setup()
  replay()
  local _, mres = run("NM_000518.5")
  local _, pres = run("NP_000509.1")
  local mrna = parse(mres.text, "NM_000518.5")
  local prot = parse(pres.text, "NP_000509.1")
  H.eq("protein", pres.db)
  H.eq("NP_000509.1", pres.accession)
  H.eq(147, prot:len())
  H.eq(prot:len(), pres.length)
  H.ok(prot:is_protein(), "GenPept record parsed as protein")
  local cds
  for _, f in ipairs(mrna.features) do if f.type == "CDS" then cds = f end end
  H.ok(cds, "CDS in NM_000518.5")
  H.eq("NP_000509.1", cds:get("protein_id"))
  H.eq(prot.seq:tostring():upper(), (cds:get("translation"):gsub("%s", "")))
  H.eq(mres.length, mrna:len())
end)

case("uniprot_accession_fetch", function()
  setup()
  replay()
  local err, res = run("P69905")
  H.eq(nil, err)
  H.eq("P69905.2", res.accession)
  H.eq("protein", res.db)
  H.eq(142, res.length)
  local rec = parse(res.text, "P69905")
  H.eq(142, rec:len())
  H.ok(rec.seq:tostring():upper():match("^MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF"), "HBA_HUMAN N-terminus")
end)

case("range_fetch_both_strands", function()
  setup()
  replay()
  local e1, plus = run("thr_plus")
  local e2, minus = run("thr_minus")
  H.eq(nil, e1); H.eq(nil, e2)
  H.eq(2611, plus.length); H.eq(2611, minus.length)
  H.eq({ 190, 2800, nil }, { plus.from, plus.to, plus.strand })
  H.eq({ 190, 2800, 2 }, { minus.from, minus.to, minus.strand })
  H.eq(4641652, plus.total_length)
  H.eq("NC_000913.3", plus.accession)
  local p, m = parse(plus.text, "plus"), parse(minus.text, "minus")
  H.eq(2611, p:len()); H.eq(2611, m:len())
  H.eq(rc(p.seq:tostring()), m.seq:tostring():upper(), "minus = reverse complement of plus")
  H.ok(p.seq:tostring() ~= m.seq:tostring())
  -- thrA (gene 945803, NC_000913.3 337..2799 plus) sits at 148..2610 of the extract
  local thrA
  for _, f in ipairs(p.features) do if f.type == "CDS" and f:get("gene") == "thrA" then thrA = f end end
  H.ok(thrA, "thrA CDS in the extract")
  H.eq({ 148, 2610, 1 }, { thrA.parts[1].start, thrA.parts[1].stop, thrA.parts[1].strand })
  local mthrA
  for _, f in ipairs(m.features) do if f.type == "CDS" and f:get("gene") == "thrA" then mthrA = f end end
  H.ok(mthrA, "thrA CDS in the reverse extract")
  H.eq({ 2611 - 2610 + 1, 2611 - 148 + 1, -1 }, { mthrA.parts[1].start, mthrA.parts[1].stop, mthrA.parts[1].strand })
end)

case("con_record_comes_with_sequence", function()
  setup()
  replay()
  local err, res = run("con_LN625376")
  H.eq(nil, err)
  H.eq(1765, res.length)
  H.ok(res.text:find("\nORIGIN", 1, true), "gbwithparts gives ORIGIN")
  local rec = parse(res.text, "con")
  H.eq(1765, rec:len())
  -- the same record as rettype=gb: a CONTIG join and no sequence
  local err2, raw = run("raw_con_gb")
  H.eq(nil, err2)
  H.ok(raw.body:find("\nCONTIG      join(CBZL010000100.1:1..1765)", 1, true), "CONTIG line")
  H.ok(not raw.body:find("\nORIGIN", 1, true), "rettype=gb has no ORIGIN")
  H.ok(raw.body:match("^LOCUS%s+LN625376%s+1765 bp.- CON "), "CON division")
end)

case("gi_number_fetch", function()
  setup()
  local calls = replay()
  local err, res = run("gi_number")
  H.eq(nil, err)
  H.eq("L09137.2", res.accession)
  H.eq(2, #calls)
  local _, ref = run("L09137")
  H.eq(ref.text, res.text)
end)

case("unknown_accession_error", function()
  setup()
  local calls = replay()
  local err, res = run("unknown")
  H.eq(nil, res)
  H.eq("NCBI does not know nuccore accession 'XX999999' (NCBI: Invalid uid XX999999 at position= 0)", err)
  H.eq(1, #calls, "no efetch after a failed summary")
  -- efetch's own error format (HTTP 200 + letter-spaced text)
  local err2, resp = run("raw_efetch_unknown")
  H.eq(nil, err2)
  H.eq(200, resp.status)
  H.eq("Failed to understand id: XX999999", ncbi._body_error(resp.body))
end)

case("size_guard_refuses_and_force", function()
  setup()
  local calls = replay()
  local err, res = run("big_refused")
  H.eq(nil, res)
  H.ok(err:find("NC_000913.3 is 4,641,652 bp, above the download limit of 1,000,000 bp", 1, true), err)
  H.ok(err:find("force", 1, true) and err:find("NC_000913.3:1-100000", 1, true), err)
  H.eq(1, #calls)
  H.ok(calls[1].url:find("esummary", 1, true), "only the esummary was made")
  -- small limit: refused, then forced
  cfgmod.options.ncbi.max_length = 1000
  local done, err2 = ncbi.wait(ncbi.fetch("L09137", {}, function() end), 5000)
  H.ok(done)
  H.ok(err2:find("L09137.2 is 2,686 bp, above the download limit of 1,000 bp", 1, true), err2)
  local done3, err3, res3 = ncbi.wait(ncbi.fetch("L09137", { force = true }, function() end), 5000)
  H.ok(done3); H.eq(nil, err3); H.eq(2686, res3.length)
  -- a range under the limit is fine even though the record is not
  cfgmod.options.ncbi.max_length = 3000
  local _, err4, res4 = ncbi.wait(ncbi.fetch("NC_000913.3:190-2800", {}, function() end), 5000)
  H.eq(nil, err4); H.eq(2611, res4.length)
  cfgmod.options.ncbi.max_length = 2000
  local _, err5 = ncbi.wait(ncbi.fetch("NC_000913.3:190-2800", {}, function() end), 5000)
  H.ok(err5:find("NC_000913.3 range 190-2800 is 2,611 bp", 1, true), err5)
end)

case("range_beyond_end_and_clipping", function()
  setup()
  record_stub("TEST1.1", ("ACGT"):rep(25))
  local _, err = ncbi.wait(ncbi.fetch("nuccore:TEST1.1:200-300", {}, function() end), 2000)
  H.ok(err:find("TEST1.1 is 100 bp long; the range 200-300 starts beyond its end", 1, true), err)
  local calls = record_stub("TEST1.1", ("ACGT"):rep(25))
  local _, err2, res = ncbi.wait(ncbi.fetch("nuccore:TEST1.1:51-500", {}, function() end), 2000)
  H.eq(nil, err2)
  H.eq({ 51, 100, 50 }, { res.from, res.to, res.length })
  H.ok(calls[2].req.url:find("seq_start=51&seq_stop=100", 1, true), calls[2].req.url)
end)

--------------------------------------------------------------------------------
-- cache
--------------------------------------------------------------------------------

case("cache_hit_refresh_list_clear", function()
  local dir = vim.fn.tempname()
  setup({ cache = true, cache_dir = dir })
  local calls = replay()
  local function fetch(s, o)
    local done, err, res = ncbi.wait(ncbi.fetch(s, o or {}, function() end), 5000)
    H.ok(done, s); H.eq(nil, err, s)
    return res
  end
  local r1 = fetch("L09137")
  H.eq(2, #calls); H.eq(false, r1.from_cache)
  H.eq(dir .. "/L09137.2.gb", r1.path)
  H.eq(r1.text, H.read(r1.path))
  local r2 = fetch("L09137")
  H.eq(2, #calls, "second fetch made no transport call")
  H.eq(true, r2.from_cache)
  H.eq(r1.text, r2.text)
  H.eq({ "L09137.2", "Cloning vector pUC19c, complete sequence", 2686, "nuccore" },
    { r2.accession, r2.title, r2.length, r2.db })
  fetch("L09137.2"); fetch("l09137.2"); fetch("20141090")
  H.eq(3, #calls, "versioned id hits directly; GI needs one summary, then shares the file")
  fetch("20141090")
  H.eq(3, #calls, "GI alias")
  local r3 = fetch("L09137", { refresh = true })
  H.eq(5, #calls, "refresh bypasses the cache")
  H.eq(false, r3.from_cache)
  -- range and strand get their own files
  local rp = fetch("NC_000913.3:190-2800")
  local rm = fetch("NC_000913.3:2800-190")
  H.eq(dir .. "/NC_000913.3_190-2800.gb", rp.path)
  H.eq(dir .. "/NC_000913.3_190-2800_rc.gb", rm.path)
  local n = #calls
  local rm2 = fetch("NC_000913.3:190..2800:-")
  H.eq(n, #calls); H.eq(true, rm2.from_cache); H.eq(rm.text, rm2.text)
  H.eq({ 190, 2800, 2, 2611 }, { rm2.from, rm2.to, rm2.strand, rm2.length })
  local pp = fetch("NP_000509.1")
  H.eq(dir .. "/NP_000509.1.gp", pp.path)
  H.eq(dir .. "/NP_000509.1.gp", ncbi.cache_path(ncbi.parse_accession("NP_000509.1")))
  H.eq(dir .. "/L09137.2.gb", ncbi.cache_path("L09137"))
  H.eq(nil, ncbi.cache_path("J01749"))
  local list = ncbi.cache_list()
  H.eq(4, #list)
  local accs = vim.tbl_map(function(e) return e.accession end, list)
  H.eq({ "L09137.2", "NC_000913.3", "NC_000913.3", "NP_000509.1" }, accs)
  for _, e in ipairs(list) do H.ok(e.size > 1000 and e.mtime > 0 and e.path:find(dir, 1, true) == 1) end
  -- expiry: an old file is ignored
  vim.uv.fs_utime(r1.path, os.time() - 40 * 86400, os.time() - 40 * 86400)
  n = #calls
  local r4 = fetch("L09137")
  H.eq(n + 2, #calls, "expired file re-fetched"); H.eq(false, r4.from_cache)
  -- cache = false neither reads nor writes
  cfgmod.options.ncbi.cache = false
  n = #calls
  local r5 = fetch("L09137")
  H.eq(n + 2, #calls); H.eq(nil, r5.path)
  cfgmod.options.ncbi.cache = true
  H.eq(4, ncbi.cache_clear())
  H.eq(0, #ncbi.cache_list())
  H.eq(nil, ncbi.cache_path("L09137"))
  vim.fn.delete(dir, "rf")
end)

--------------------------------------------------------------------------------
-- summary / search / gene
--------------------------------------------------------------------------------

case("summary_post_and_order", function()
  setup()
  local calls = replay()
  local err, items = run("summary_post")
  H.eq(nil, err)
  H.eq("POST", calls[1].method)
  H.ok(calls[1].body:find("id=L09137%2CJ01749%2CNC_001422", 1, true), calls[1].body)
  H.ok(not calls[1].url:find("?", 1, true), "POST url carries no query")
  H.eq({ "L09137.2", "J01749.1", "NC_001422.1" }, vim.tbl_map(function(i) return i.accession end, items))
  H.eq({ 2686, 4361, 5386 }, vim.tbl_map(function(i) return i.length end, items))
  local phi = items[3]
  H.eq("circular", phi.topology)
  H.ok(phi.uid and phi.title and phi.organism and phi.updated and phi.moltype, vim.inspect(phi))
  -- order follows the request even when NCBI answers in another order; unknown ids are reported
  setup()
  stub(function()
    return json_resp({ error = "Invalid uid NOPE1 at position= 1", result = { uids = { "7", "3" },
      ["3"] = { uid = "3", caption = "A00001", accessionversion = "A00001.1", slen = 10 },
      ["7"] = { uid = "7", caption = "B00002", accessionversion = "B00002.3", slen = 20 } } })
  end)
  local done, err2, list, missing = nil, nil, nil, nil
  ncbi.summary("nuccore", { "A00001", "NOPE1", "B00002.3" }, function(e, l, m) done, err2, list, missing = true, e, l, m end)
  vim.wait(2000, function() return done end, 5)
  H.eq(nil, err2)
  H.eq({ "A00001.1", "B00002.3" }, vim.tbl_map(function(i) return i.accession end, list))
  H.eq({ "NOPE1" }, missing)
end)

case("search_results", function()
  setup()
  replay()
  local err, res = run("search_puc19")
  H.eq(nil, err)
  H.ok(res.count >= 1 and #res.items == math.min(res.count, 20), vim.inspect(res.count))
  for _, it in ipairs(res.items) do
    H.ok(it.accession:match("%.%d+$") and it.length > 0 and it.title:lower():find("puc19", 1, true), vim.inspect(it))
  end
  local err2, res2 = run("search_protein")
  H.eq(nil, err2)
  H.eq(1, res2.count)
  H.eq("NP_000509.1", res2.items[1].accession)
  H.eq(147, res2.items[1].length)
  local err3, res3 = run("search_nohit")
  H.eq(nil, err3)
  H.eq(0, res3.count); H.eq({}, res3.items)
  H.ok(vim.tbl_contains(res3.warnings, "phrasesnotfound: zzqxqzzq[title]"), vim.inspect(res3.warnings))
end)

case("gene_coordinates_match_ncbi_gene_records", function()
  setup()
  replay()
  local err, genes = run("gene_BRCA1")
  H.eq(nil, err)
  local g = genes[1]
  H.eq({ "672", "BRCA1", "Homo sapiens", "NC_000017.11", 2, "17q21.31" },
    { g.geneid, g.symbol, g.organism, g.chr_accession, g.strand, g.location_text })
  H.ok(g.start <= g.stop)
  -- NCBI's own gene record states the annotation as "NC_000017.11 (start..stop, complement)"
  local _, raw = run("raw_gene_672")
  local acc, a, b, comp = raw.body:match("Annotation:%s+Chromosome %S+ (%S+) %((%d+)%.%.(%d+)(, complement)%)")
  H.ok(acc, "annotation line in " .. raw.body)
  H.eq({ acc, tonumber(a), tonumber(b), 2 }, { g.chr_accession, g.start, g.stop, g.strand })
  local err2, genes2 = run("gene_thrA")
  H.eq(nil, err2)
  local t = genes2[1]
  H.eq({ "945803", "thrA", "NC_000913.3", 337, 2799, 1 }, { t.geneid, t.symbol, t.chr_accession, t.start, t.stop, t.strand })
  local _, raw2 = run("raw_gene_945803")
  local acc2, a2, b2, tail = raw2.body:match("Annotation:%s+(%S+) %((%d+)%.%.(%d+)(.-)%)")
  H.eq({ "NC_000913.3", 337, 2799, "" }, { acc2, tonumber(a2), tonumber(b2), tail })
  _ = comp
end)

case("gene_location_conversion", function()
  -- esummary genomicinfo: 0-based chrstart/chrstop, stop < start on the minus strand
  H.eq({ 337, 2799, 1 }, { ncbi.gene_location({ chrstart = 336, chrstop = 2798 }) })
  H.eq({ 43044295, 43170327, 2 }, { ncbi.gene_location({ chrstart = 43170326, chrstop = 43044294 }) })
  H.eq({ 1, 1, 1 }, { ncbi.gene_location({ chrstart = 0, chrstop = 0 }) })
  H.eq({ 11, 20, 2 }, { ncbi.gene_location({ chrstart = "19", chrstop = "10" }) })
  H.eq(nil, ncbi.gene_location({ chrstart = "", chrstop = "" }))
  H.eq('(BRCA1[sym]) AND "Homo sapiens"[orgn] AND alive[prop]', ncbi._gene_term("BRCA1", "Homo sapiens"))
  H.eq("BRCA1[sym] AND alive[prop]", ncbi._gene_term(" BRCA1 "))
  H.eq("breast cancer 1 AND alive[prop]", ncbi._gene_term("breast cancer 1", ""))
  -- unannotated genes are skipped; no gene at all is an error
  setup()
  stub(function(req)
    if req.url:find("esearch", 1, true) then return json_resp({ esearchresult = { count = "1", idlist = { "9" } } }) end
    return json_resp({ result = { uids = { "9" }, ["9"] = { uid = "9", name = "X", genomicinfo = {} } } })
  end)
  local _, err = ncbi.wait(ncbi.gene("X", nil, {}, function() end), 2000)
  H.eq("no annotated genomic location for X", err)
  setup()
  stub(function() return json_resp({ esearchresult = { count = "0", idlist = {} } }) end)
  local _, err2 = ncbi.wait(ncbi.gene("NOPE", "Homo sapiens", {}, function() end), 2000)
  H.eq("no current gene found for NOPE in Homo sapiens", err2)
end)

--------------------------------------------------------------------------------
-- errors
--------------------------------------------------------------------------------

case("error_messages", function()
  local function fetch_err(responder, opts)
    setup(opts)
    local calls = stub(responder)
    local done, err = ncbi.wait(ncbi.fetch("L09137", {}, function() end), 5000)
    H.ok(done)
    return err, calls
  end
  -- HTTP 500 is retried, then reported
  local err, calls = fetch_err(function() return nil, { status = 500, headers = {}, body = "<html>oops</html>" } end)
  H.eq("NCBI esummary failed: HTTP 500 Internal Server Error", err)
  H.eq(4, #calls, "1 + retries")
  -- NCBI JSON error on HTTP 400 (not retried)
  err, calls = fetch_err(function() return json_resp({ error = "Invalid db name specified: nuccor" }, 400) end)
  H.eq("NCBI esummary failed: HTTP 400 Bad Request: Invalid db name specified: nuccor", err)
  H.eq(1, #calls)
  -- XML error body
  err = fetch_err(function()
    return nil, { status = 200, headers = {}, body = "<?xml version=\"1.0\"?>\n<eSummaryResult><ERROR>Empty id list - nothing todo</ERROR></eSummaryResult>" }
  end)
  H.eq("NCBI esummary: Empty id list - nothing todo", err)
  -- efetch letter-spaced error after a good summary
  err = fetch_err(function(req)
    if req.url:find("esummary", 1, true) then return json_resp(fake_summary("L09137.2", 2686)) end
    return nil, { status = 200, headers = {}, body = "Error: F a i l e d  t o  u n d e r s t a n d  i d :  L 0 9 1 3 7 . 2 \n" }
  end)
  H.eq("NCBI efetch L09137.2: Failed to understand id: L09137.2", err)
  -- garbage from efetch
  err = fetch_err(function(req)
    if req.url:find("esummary", 1, true) then return json_resp(fake_summary("L09137.2", 2686)) end
    return nil, { status = 200, headers = {}, body = "Supplied id parameter is empty." }
  end)
  H.eq("NCBI efetch L09137.2: unexpected reply: Supplied id parameter is empty.", err)
  -- 429 exhausted
  err, calls = fetch_err(function() return json_resp({ error = "API rate limit exceeded", count = "4" }, 429) end)
  H.ok(err:find("HTTP 429 Too Many Requests (rate limit; set an API key", 1, true) and err:find("API rate limit exceeded", 1, true), err)
  H.eq(4, #calls)
  -- timeout: reported, not retried
  err, calls = fetch_err(function() return "timed out after 60 s (timeout_ms = 60000)" end)
  H.eq("NCBI esummary: timed out after 60 s (timeout_ms = 60000)", err)
  H.eq(1, #calls)
  -- network errors are retried
  err, calls = fetch_err(function() return "network error: Could not resolve host: eutils.ncbi.nlm.nih.gov (curl exit 6)" end)
  H.eq("NCBI esummary: network error: Could not resolve host: eutils.ncbi.nlm.nih.gov (curl exit 6)", err)
  H.eq(4, #calls)
  -- a transport that throws
  err = fetch_err(function() error("boom") end)
  H.ok(err:find("transport error", 1, true) and err:find("boom", 1, true), err)
  -- bad accession text: callback still asynchronous
  setup()
  local done, err2 = ncbi.wait(ncbi.fetch("hello", {}, function() end), 1000)
  H.ok(done); H.ok(err2:find("not an accession", 1, true))
  -- replay names the missing key
  setup()
  ncbi.set_transport(ncbi.transports.replay(HTTP))
  local _, err3 = ncbi.wait(ncbi.fetch("J01636", {}, function() end), 2000)
  H.ok(err3:find("replay: no recorded response for key %x+ %(GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi%?db=nuccore&id=J01636"), err3)
end)

case("curl_missing", function()
  setup()
  ncbi.set_transport(ncbi.transports.curl)
  local path = vim.env.PATH
  vim.env.PATH = "/nonexistent"
  local done, err = ncbi.wait(ncbi.fetch("L09137", {}, function() end), 3000)
  vim.env.PATH = path
  H.ok(done)
  H.eq("NCBI esummary: curl is not installed (not found in $PATH); plasmid.nvim needs it to talk to NCBI", err)
end)

--------------------------------------------------------------------------------
-- rate limiting, retries, cancel, scheduling
--------------------------------------------------------------------------------

local function spacing(opts, n)
  setup(opts)
  local calls = stub(function() return json_resp(fake_summary("L09137.2", 2686)) end)
  local hs = {}
  for i = 1, n do hs[i] = ncbi.summary("nuccore", { "L09137" }, function() end) end
  for i = 1, n do H.ok((ncbi.wait(hs[i], 5000)), "summary " .. i) end
  local gaps = {}
  for i = 2, #calls do gaps[#gaps + 1] = calls[i].t - calls[i - 1].t end
  return gaps, calls
end

case("rate_limit_spacing", function()
  -- without a key: >= 334 ms between request starts (default 350)
  local gaps = spacing({ min_interval_ms = false }, 4)
  H.eq(3, #gaps)
  for _, g in ipairs(gaps) do H.ok(g >= 349 and g < 600, ("gap %.1f ms"):format(g)) end
  H.log(("    rate limit, no key: gaps %s ms"):format(table.concat(vim.tbl_map(function(g) return ("%.0f"):format(g) end, gaps), ", ")))
  -- with an API key: >= 100 ms (default 110)
  local gaps2, calls2 = spacing({ min_interval_ms = false, api_key = "SECRETKEY123" }, 4)
  for _, g in ipairs(gaps2) do H.ok(g >= 109 and g < 300, ("gap %.1f ms"):format(g)) end
  H.log(("    rate limit, API key: gaps %s ms"):format(table.concat(vim.tbl_map(function(g) return ("%.0f"):format(g) end, gaps2), ", ")))
  H.ok(calls2[1].req.url:find("api_key=SECRETKEY123", 1, true), "key sent to NCBI")
  H.eq(110, ncbi.config().min_interval_ms)
  -- $NCBI_API_KEY fallback
  setup({ min_interval_ms = false })
  vim.env.NCBI_API_KEY = "ENVKEY"
  local c = ncbi.config()
  vim.env.NCBI_API_KEY = nil
  H.eq({ "ENVKEY", 110 }, { c.api_key, c.min_interval_ms })
  H.eq({ nil, 350 }, { ncbi.config().api_key, ncbi.config().min_interval_ms })
end)

case("api_key_never_in_keys_files_or_history", function()
  local dir = vim.fn.tempname()
  setup({ api_key = "SECRETKEY123", email = "user@example.org", cache = true, cache_dir = dir .. "/cache" })
  local inner = function(req, done)
    if req.url:find("esummary", 1, true) then done(json_resp(fake_summary("TEST1.1", 8))) return end
    done(nil, { status = 200, headers = {}, body = fake_gb("TEST1.1", "ACGTACGT") })
  end
  local seen = {}
  local rec = ncbi.transports.record(dir .. "/rec", inner)
  ncbi.set_transport(function(req, done) seen[#seen + 1] = req; return rec(req, done) end)
  local done, err = ncbi.wait(ncbi.fetch("nuccore:TEST1", {}, function() end), 2000)
  H.ok(done); H.eq(nil, err)
  H.ok(seen[1].url:find("email=user%40example.org", 1, true) and seen[1].url:find("api_key=SECRETKEY123", 1, true))
  local all = {}
  for _, f in ipairs(vim.fn.globpath(dir, "**/*", false, true)) do
    all[#all + 1] = f
    if vim.fn.isdirectory(f) == 0 then H.ok(not H.read(f):find("SECRETKEY123", 1, true), "key in " .. f) end
    H.ok(not f:find("SECRETKEY123", 1, true), f)
  end
  H.ok(#all >= 6, vim.inspect(all))
  H.ok(not vim.inspect(ncbi.history):find("SECRETKEY123", 1, true), "history")
  local plain = vim.tbl_extend("force", {}, seen[1])
  plain.url = plain.url:gsub("&email=[^&]*", ""):gsub("&api_key=[^&]*", "")
  H.eq(ncbi.request_key(plain), ncbi.request_key(seen[1]))
  -- POST body is scrubbed too
  H.eq(ncbi.request_key({ method = "POST", url = "u", body = "db=x&id=1&tool=t" }),
    ncbi.request_key({ method = "POST", url = "u", body = "db=x&id=1&tool=t&api_key=SECRETKEY123" }))
  vim.fn.delete(dir, "rf")
end)

case("retries_with_backoff", function()
  -- 429 with Retry-After: 1 -> waits >= 1 s, then succeeds
  setup({ retry_base_ms = 10 })
  local calls = stub(function(_, n)
    if n == 1 then return json_resp({ error = "API rate limit exceeded" }, 429, { ["retry-after"] = "1" }) end
    return json_resp(fake_summary("L09137.2", 2686))
  end)
  local done, err, items = ncbi.wait(ncbi.summary("nuccore", { "L09137" }, function() end), 5000)
  H.ok(done); H.eq(nil, err); H.eq("L09137.2", items[1].accession)
  H.eq(2, #calls)
  local g = calls[2].t - calls[1].t
  H.ok(g >= 990 and g < 1500, ("Retry-After gap %.0f ms"):format(g))
  -- 503, 503, 200: exponential backoff 40, 80 ms
  setup({ retry_base_ms = 40 })
  calls = stub(function(_, n)
    if n <= 2 then return nil, { status = 503, headers = {}, body = "" } end
    return json_resp(fake_summary("L09137.2", 2686))
  end)
  done, err = ncbi.wait(ncbi.summary("nuccore", { "L09137" }, function() end), 5000)
  H.ok(done); H.eq(nil, err); H.eq(3, #calls)
  local g1, g2 = calls[2].t - calls[1].t, calls[3].t - calls[2].t
  H.ok(g1 >= 35 and g2 >= 75 and g2 > g1 * 1.5, ("backoff gaps %.0f, %.0f ms"):format(g1, g2))
  -- old-style 200 + rate-limit JSON is retried as well
  setup()
  calls = stub(function(_, n)
    if n == 1 then return json_resp({ error = "API rate limit exceeded", count = "11" }, 200) end
    return json_resp(fake_summary("L09137.2", 2686))
  end)
  done, err = ncbi.wait(ncbi.summary("nuccore", { "L09137" }, function() end), 5000)
  H.eq(nil, err); H.eq(2, #calls)
  -- esearch "Search Backend failed" (seen live: HTTP 200 + ERROR) is transient: retried
  setup()
  calls = stub(function(_, n)
    if n == 1 then return json_resp({ esearchresult = { ERROR = "Search Backend failed: \t\t\n" } }) end
    if n == 2 then return json_resp({ esearchresult = { count = "0", idlist = {} } }) end
  end)
  local _, err4, res4 = ncbi.wait(ncbi.search("nuccore", "x", {}, function() end), 5000)
  H.eq(nil, err4); H.eq(0, res4.count); H.eq(2, #calls)
  setup({ retries = 1 })
  stub(function() return json_resp({ esearchresult = { ERROR = "Search Backend failed: \t\t\n" } }) end)
  _, err4 = ncbi.wait(ncbi.search("nuccore", "x", {}, function() end), 5000)
  H.eq("NCBI esearch: Search Backend failed:", err4)
  -- retries = 0
  setup({ retries = 0 })
  calls = stub(function() return nil, { status = 502, headers = {}, body = "" } end)
  _, err = ncbi.wait(ncbi.summary("nuccore", { "L09137" }, function() end), 5000)
  H.eq("NCBI esummary failed: HTTP 502 Bad Gateway", err); H.eq(1, #calls)
end)

case("cancel", function()
  -- in flight: transport cancel is called, callback never runs
  setup()
  local calls = stub(function() return "hang" end)
  local ran = false
  local h = ncbi.fetch("L09137", {}, function() ran = true end)
  vim.wait(200, function() return #calls == 1 end, 5)
  H.eq(1, #calls)
  H.eq(true, h:cancel())
  H.eq(1, calls.cancelled)
  H.eq(false, (ncbi.wait(h, 300)))
  H.eq(false, ran)
  H.eq(false, h:cancel(), "second cancel is a no-op")
  -- queued: never reaches the transport
  setup({ min_interval_ms = 200 })
  calls = stub(function() return json_resp(fake_summary("L09137.2", 2686)) end)
  local r1, r2 = false, false
  local h1 = ncbi.summary("nuccore", { "A" }, function() r1 = true end)
  local h2 = ncbi.summary("nuccore", { "B" }, function() r2 = true end)
  vim.wait(50, function() return false end)
  H.eq(1, #calls)
  H.eq(true, h2:cancel())
  vim.wait(500, function() return false end)
  H.eq(1, #calls); H.eq(true, r1); H.eq(false, r2)
  H.eq(false, h1:cancel(), "cancel after the callback ran")
  -- during a retry backoff
  setup({ retry_base_ms = 300 })
  calls = stub(function() return nil, { status = 503, headers = {}, body = "" } end)
  local r3 = false
  local h3 = ncbi.summary("nuccore", { "A" }, function() r3 = true end)
  vim.wait(100, function() return false end)
  h3:cancel()
  vim.wait(600, function() return false end)
  H.eq(1, #calls); H.eq(false, r3)
  -- between result and scheduled callback
  setup()
  stub(function() return json_resp(fake_summary("L09137.2", 2686)) end)
  local r4 = false
  local h4 = ncbi.summary("nuccore", { "A" }, function() r4 = true end)
  vim.wait(1000, function() return h4.done end, 1)
  if not h4.called then H.eq(true, h4:cancel()) end
  vim.wait(50, function() return false end)
  H.eq(h4.called, r4)
end)

case("callbacks_run_on_main_loop", function()
  local dir = vim.fn.tempname()
  setup({ cache = true, cache_dir = dir })
  replay()
  local fast, sync = {}, {}
  local function probe(name, h_fn)
    local returned = false
    local h = h_fn(function()
      fast[#fast + 1] = vim.in_fast_event() and name or nil
      sync[#sync + 1] = (not returned) and name or nil
      pcall(vim.api.nvim_buf_set_lines, 0, 0, -1, false, { name })  -- UI code may touch buffers
    end)
    returned = true
    H.ok((ncbi.wait(h, 5000)), name)
  end
  probe("fetch", function(cb) return ncbi.fetch("L09137", {}, cb) end)
  probe("cached fetch", function(cb) return ncbi.fetch("L09137", {}, cb) end)
  probe("bad accession", function(cb) return ncbi.fetch("hello", {}, cb) end)
  probe("search", function(cb) return ncbi.search("nuccore", "zzqxqzzq[title]", {}, cb) end)
  probe("gene", function(cb) return ncbi.gene("thrA", "Escherichia coli str. K-12 substr. MG1655", {}, cb) end)
  H.eq({}, fast, "callbacks in fast context")
  H.eq({}, sync, "callbacks run synchronously")
  H.eq({ "gene" }, vim.api.nvim_buf_get_lines(0, 0, -1, false))
  vim.fn.delete(dir, "rf")
end)

case("reset_clears_queue", function()
  setup({ min_interval_ms = 1000 })
  local calls = stub(function() return json_resp(fake_summary("L09137.2", 2686)) end)
  local ran = 0
  for _ = 1, 3 do ncbi.summary("nuccore", { "A" }, function() ran = ran + 1 end) end
  vim.wait(100, function() return false end)
  H.eq(1, #calls)
  ncbi.reset()
  H.eq(ncbi.transports.curl, ncbi.get_transport())
  vim.wait(1200, function() return false end)
  H.eq(1, #calls); H.eq(1, ran)
end)

--------------------------------------------------------------------------------
-- curl transport against a local HTTP server
--------------------------------------------------------------------------------

local function local_server()
  local server = vim.uv.new_tcp()
  server:bind("127.0.0.1", 0)
  local port = server:getsockname().port
  local seen = {}
  local clients = {}
  server:listen(16, function()
    local c = vim.uv.new_tcp()
    server:accept(c)
    clients[#clients + 1] = c
    local buf = ""
    c:read_start(function(rerr, chunk)
      if rerr or not chunk then if not c:is_closing() then c:close() end return end
      buf = buf .. chunk
      local head, rest = buf:match("^(.-)\r\n\r\n(.*)$")
      if not head then return end
      local len = tonumber(head:lower():match("content%-length:%s*(%d+)") or "0")
      if #rest < len then return end
      local method, path = head:match("^(%u+) (%S+)")
      seen[#seen + 1] = { method = method, path = path, head = head, body = rest }
      if path:find("^/slow") then return end
      local status, body, extra = 200, "hello\nworld", "X-Test: yes\r\n"
      if path:find("^/echo") then body = method .. " " .. rest end
      if path:find("^/busy") then status, body, extra = 429, '{"error":"API rate limit exceeded"}', "Retry-After: 0\r\n" end
      c:write(("HTTP/1.1 %d X\r\nContent-Length: %d\r\nContent-Type: text/plain\r\n%sConnection: close\r\n\r\n%s"):format(
        status, #body, extra, body), function() if not c:is_closing() then c:close() end end)
    end)
  end)
  return port, seen, function()
    for _, c in ipairs(clients) do if not c:is_closing() then c:close() end end
    server:close()
  end
end

case("curl_transport_local_server", function()
  if vim.fn.executable("curl") ~= 1 then H.log("    (curl not installed: skipped)") return end
  local port, seen, stop = local_server()
  local saved = { vim.env.no_proxy, vim.env.NO_PROXY }
  vim.env.no_proxy, vim.env.NO_PROXY = "127.0.0.1,localhost", "127.0.0.1,localhost"
  local base = ("http://127.0.0.1:%d/"):format(port)
  local function call(req)
    local out
    ncbi.transports.curl(req, function(err, resp) out = { err = err, resp = resp, fast = vim.in_fast_event() } end)
    vim.wait(8000, function() return out ~= nil end, 5)
    return out
  end
  local ok, e = pcall(function()
    local r = call({ url = base .. "ok?a=1", method = "GET", timeout_ms = 5000 })
    H.eq(nil, r.err)
    H.eq({ 200, "hello\nworld", "yes", false }, { r.resp.status, r.resp.body, r.resp.headers["x-test"], r.fast })
    r = call({ url = base .. "echo", method = "POST", body = "db=nuccore&id=A%2CB", timeout_ms = 5000 })
    H.eq("POST db=nuccore&id=A%2CB", r.resp.body)
    H.ok(seen[#seen].head:lower():find("content-type: application/x-www-form-urlencoded", 1, true))
    H.ok(seen[#seen].head:find("User-Agent: curl/", 1, true), "curl's default User-Agent")
    r = call({ url = base .. "busy", method = "GET", timeout_ms = 5000 })
    H.eq({ 429, "0" }, { r.resp.status, r.resp.headers["retry-after"] })
    local t0 = now()
    r = call({ url = base .. "slow", method = "GET", timeout_ms = 1000 })
    H.eq("timed out after 1 s (timeout_ms = 1000)", r.err)
    H.ok(now() - t0 < 4000)
    r = call({ url = "http://127.0.0.1:1/", method = "GET", timeout_ms = 3000 })
    H.ok(r.err and r.err:find("^network error: .*%(curl exit 7%)$"), tostring(r.err))
    -- cancel kills curl; done is never called
    local got = false
    local c = ncbi.transports.curl({ url = base .. "slow2", method = "GET", timeout_ms = 5000 }, function() got = true end)
    vim.wait(2000, function() return #seen >= 6 end, 5)
    c.cancel()
    vim.wait(500, function() return got end, 5)
    H.eq(false, got)
    -- whole client path through curl with the local server as base_url
    setup({ base_url = base, retry_base_ms = 1 })
    ncbi.set_transport(ncbi.transports.curl)
    local done, err = ncbi.wait(ncbi.summary("nuccore", { "A" }, function() end), 5000)
    H.ok(done)
    H.eq("NCBI esummary: unreadable response (not JSON)", err)
    H.ok(seen[#seen].path:find("^/esummary.fcgi%?db=nuccore&id=A&retmode=json&tool=plasmid.nvim$"), seen[#seen].path)
  end)
  vim.env.no_proxy, vim.env.NO_PROXY = saved[1], saved[2]
  stop()
  if not ok then error(e, 0) end
end)

case("headers_parse_last_block", function()
  local raw = "HTTP/1.1 200 Connection established\r\n\r\nHTTP/1.1 301 Moved\r\nLocation: x\r\n\r\n"
    .. "HTTP/2 429 \r\nRetry-After: 2\r\nContent-Type: application/json\r\n\r\n"
  local h = ncbi._parse_headers(raw)
  H.eq({ 429, "2", "application/json" }, { h.status, h.headers["retry-after"], h.headers["content-type"] })
  H.eq(nil, ncbi._parse_headers(""))
end)

case("config_defaults_and_override", function()
  cfgmod.options.ncbi = nil
  vim.env.NCBI_API_KEY = nil
  local c = ncbi.config()
  H.eq({ "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", nil, nil, "plasmid.nvim", 5000000, 60000, 3, true, 30 },
    { c.base_url, c.api_key, c.email, c.tool, c.max_length, c.timeout_ms, c.retries, c.cache, c.cache_days })
  H.eq(vim.fn.stdpath("cache") .. "/plasmid/ncbi", c.cache_dir)
  cfgmod.options.ncbi = { max_length = 10, base_url = "http://x/y", api_key = "K" }
  c = ncbi.config()
  H.eq({ 10, "http://x/y/", "K", 3 }, { c.max_length, c.base_url, c.api_key, c.retries })
  cfgmod.options.ncbi = nil
end)

--------------------------------------------------------------------------------
-- live mode
--------------------------------------------------------------------------------

case("live_ncbi", function()
  if os.getenv("PLASMID_LIVE") ~= "1" then return end
  local dir = vim.fn.tempname()
  local summary = {}
  for _, c in ipairs(cases.list) do
    if c.kind == "fetch" or c.kind == "search" or c.kind == "gene" then
      setup(nil, true)  -- keep the rate-limit clock between live calls
      replay()
      local e0, want = run(c.name)
      ncbi.reset()
      cfgmod.options.ncbi = { cache = true, cache_dir = dir }
      vim.env.NCBI_API_KEY = os.getenv("NCBI_API_KEY")
      local t0 = now()
      local h = cases.run(ncbi, c, function() end, { cache = true, cache_dir = dir })
      local done, err, got = ncbi.wait(h, 120000)
      H.ok(done, c.name .. " live: no answer")
      local ms = now() - t0
      if err and err:find("Search Backend failed", 1, true) then
        -- NCBI-side outage (persisted for minutes during development); the client
        -- retried and reported it; not a client regression
        H.log(("    %s: skipped, NCBI answered %q after retries"):format(c.name, err))
      elseif c.expect_err then
        H.eq(e0, err, c.name .. " live error")
      elseif c.kind == "fetch" then
        H.eq(nil, err, c.name)
        H.eq({ want.accession, want.length, want.db }, { got.accession, got.length, got.db }, c.name .. " live vs recorded")
        local rec = parse(got.text, c.name)
        H.eq(got.length, rec:len(), c.name .. " parsed length")
        if want.text ~= got.text then H.log(("    %s: live text differs from the recording (%s)"):format(c.name,
          tostring(H.first_diff(want.text, got.text)):gsub("\n", " | "))) end
      elseif c.kind == "search" then
        H.eq(nil, err, c.name)
        H.ok(got.count >= want.count, c.name .. " live count " .. got.count)
      elseif c.kind == "gene" then
        H.eq(nil, err, c.name)
        H.eq({ want[1].geneid, want[1].chr_accession, want[1].start, want[1].stop, want[1].strand },
          { got[1].geneid, got[1].chr_accession, got[1].start, got[1].stop, got[1].strand }, c.name .. " live coordinates")
      end
      local odd = {}
      for _, r in ipairs(ncbi.history) do
        if r.err or r.status ~= 200 then odd[#odd + 1] = tostring(r.status or r.err) end
      end
      summary[#summary + 1] = ("%s %.0f ms %d req%s"):format(c.name, ms, #ncbi.history,
        #odd > 0 and (" (retried: " .. table.concat(odd, ",") .. ")") or "")
    end
  end
  H.log("    live: " .. table.concat(summary, "\n          "))
  vim.fn.delete(dir, "rf")
  ncbi.reset()
end)

case("teardown", function()
  ncbi.reset()
  cfgmod.options.ncbi = nil
end)

return T
