-- Re-record the NCBI E-utilities responses used by tests/spec/ncbi_spec.lua.
--
--   nvim -l tests/tools/ncbi_record.lua            re-record every case (wipes the directory)
--   nvim -l tests/tools/ncbi_record.lua gene_      re-record cases whose name contains "gene_"
--
-- Talks to NCBI live through the curl transport wrapped by the record
-- transport. Requests are rate limited by the client (334 ms without an API
-- key). An API key from $NCBI_API_KEY is used if set but never written to the
-- fixtures (keys, URLs and bodies are scrubbed). No email is sent.
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({ root .. "/lua/?.lua", root .. "/lua/?/init.lua", root .. "/tests/?.lua", package.path }, ";")
vim.opt.runtimepath:prepend(root)

local ncbi = require("plasmid.ncbi")
local cases = require("ncbi_cases")
local dir = root .. "/tests/fixtures/ncbi_http"
local filter = arg[1]

local function out(...) io.stdout:write(table.concat(vim.tbl_map(tostring, { ... }), "\t"), "\n") end

vim.fn.mkdir(dir, "p")
local manifest_path = dir .. "/manifest.json"
local manifest = {}
if filter then
  local fh = io.open(manifest_path, "rb")
  if fh then
    manifest = vim.json.decode(fh:read("*a")).cases or {}
    fh:close()
  end
else
  for _, f in ipairs(vim.fn.globpath(dir, "*", false, true)) do os.remove(f) end
end

local current
local recorder = ncbi.transports.record(dir, ncbi.transports.curl)
ncbi.reset()
ncbi.set_transport(function(req, done)
  local t0 = vim.uv.hrtime()
  return recorder(req, function(err, resp)
    current.requests[#current.requests + 1] = {
      key = ncbi.request_key(req), method = req.method, url = ncbi._scrub_url(req.url),
      status = resp and resp.status, err = err, bytes = resp and #(resp.body or "") or 0,
      ms = math.floor((vim.uv.hrtime() - t0) / 1e6 + 0.5),
      -- error bodies (429/5xx are not stored as fixtures) are kept for the report
      -- NCBI's 429 JSON echoes the client IP as "api-key": redact it
      sample = resp and resp.status ~= 200
        and ((resp.body or ""):sub(1, 200):gsub('"api%-key":"[^"]*"', '"api-key":"<client address redacted>"')) or nil,
    }
    done(err, resp)
  end)
end)

local by = {}
for i, m in ipairs(manifest) do by[m.name] = i end
local t_all = vim.uv.hrtime()
local nreq, nbytes = 0, 0
for _, case in ipairs(cases.list) do
  if not filter or case.name:find(filter, 1, true) then
    current = { name = case.name, kind = case.kind, note = case.note, requests = {} }
    local t0 = vim.uv.hrtime()
    local h = cases.run(ncbi, case, function() end, { cache = false })
    local done, err, res = ncbi.wait(h, 180000)
    current.ms = math.floor((vim.uv.hrtime() - t0) / 1e6 + 0.5)
    current.err = err
    if case.kind == "fetch" and res then current.accession, current.length = res.accession, res.length end
    if case.kind == "search" and res then current.count = res.count end
    if case.kind == "gene" and res then
      current.genes = vim.tbl_map(function(g)
        return ("%s %s %s:%d-%d strand %d"):format(g.geneid, g.symbol, g.chr_accession, g.start, g.stop, g.strand)
      end, res)
    end
    if not done then current.err = "timed out waiting" end
    for _, r in ipairs(current.requests) do nreq = nreq + 1; nbytes = nbytes + r.bytes end
    out(case.name, done and "done" or "TIMEOUT", #current.requests .. " req", current.ms .. " ms", err or (res and res.accession) or "")
    if by[case.name] then manifest[by[case.name]] = current else manifest[#manifest + 1] = current; by[case.name] = #manifest end
  end
end
local fh = assert(io.open(manifest_path, "wb"))
fh:write(vim.json.encode({ recorded = os.date("!%Y-%m-%dT%H:%M:%SZ"), cases = manifest }), "\n")
fh:close()
out(("%d requests, %d bytes, %.1f s"):format(nreq, nbytes, (vim.uv.hrtime() - t_all) / 1e9))
