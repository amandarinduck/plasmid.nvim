--- Asynchronous NCBI E-utilities client (efetch / esummary / esearch).
---
--- Pure Lua; the only external program is `curl`, spawned with vim.system.
--- Every public call returns a handle immediately and reports through a
--- callback that always runs on the main loop (vim.schedule), so callbacks may
--- touch buffers and windows.
---
--- API
---   M.config()                      effective options (see M.defaults)
---   M.parse_accession(text)         -> spec | nil, err
---   M.fetch(spec|text, opts, cb)    -> handle; cb(err, res)
---   M.summary(db, ids, cb)          -> handle; cb(err, items, missing_ids)
---   M.search(db, term, opts, cb)    -> handle; cb(err, { count, items, warnings })
---   M.gene(symbol_or_term, organism, opts, cb) -> handle; cb(err, genes)
---   M.eutil(name, params, opts, cb) -> handle; cb(err, resp)  (low level, raw response)
---   handle:cancel()                 stop the request: curl is killed, queued
---                                   requests are dropped and the callback is
---                                   NEVER called afterwards (returns true when
---                                   it prevented a pending callback)
---   M.wait(handle, timeout_ms)      -> done, err, res  (blocks with vim.wait)
---   M.set_transport(fn), M.transports.{curl, replay(dir), record(dir, inner)}, M.reset()
---   M.cache_path(spec[, accession_version]), M.cache_list(), M.cache_clear()
---   M.request_key(request)          stable fixture key (api_key/email removed)
---
--- Accession rules (M.parse_accession). Input is trimmed; `id` keeps the
--- typed spelling (case included); classification is case-insensitive.
---   [db:]ID[:RANGE][:STRAND]
---   db        "nuccore:" / "nucleotide:" / "protein:" prefix forces the
---             database (any letters/digits/_ id is then accepted as typed)
---   ID        nuccore:
---               GenBank/EMBL/DDBJ  1 letter + 5 digits      L09137 (not O/P/Q: UniProt)
---                                  2 letters + 6 or 8 digits AB000263
---               WGS / TSA          4 letters + 8-9 digits   AAAA01000001
---                                  6 letters + 9-11 digits
---               MGA                5 letters + 7 digits
---               RefSeq             AC_ NC_ NG_ NM_ NR_ NT_ NW_ NZ_ XM_ XR_ + digits
---                                  (NZ_ etc. may wrap a GenBank stem: NZ_CP009072.1,
---                                  NZ_JAAAAA010000001.1)
---               GI number          bare digits              20141090
---             protein:
---               RefSeq             NP_ XP_ WP_ YP_ AP_ + digits
---               GenPept            3 letters + 5 or 7 digits AAA12345
---               UniProt            [OPQ][0-9][A-Z0-9]{3}[0-9]               P69905
---                                  [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2} A0A024R161
---             any of them may carry a .version suffix (NM_000518.5)
---   RANGE     FROM-TO or FROM..TO, 1-based inclusive; commas allowed in
---             numbers (1,000-5,000). FROM > TO means the reverse strand
---             (spec.from/to are always returned ascending, strand = 2).
---   STRAND    ":-" / ":minus" / ":rc" -> strand 2 (reverse complement),
---             ":+" / ":plus" -> strand 1. Allowed with or without a range,
---             nucleotides only. A reversed range with ":+" is an error.
---   Anything else is rejected with a message naming the accepted forms.
---
--- fetch(spec, opts, cb): opts = { force = skip the size guard, refresh = ignore
--- the cache }. Unless served from the cache, an esummary runs first (length and
--- accession.version); a record (or range) longer than max_length is refused
--- without downloading it. nuccore is fetched with rettype=gbwithparts so CON
--- records come back with their sequence; protein with rettype=gp.
--- res = { text, path, from_cache, accession, title, length, db,
---         uid, total_length, from, to, strand, id }
---   length = residues in the returned record (the range length for ranges).
---
--- Cache: cache_dir/<ACCESSION.VERSION>[_FROM-TO][_rc].gb|.gp (raw flat file)
--- plus a <file>.json sidecar with the summary fields; aliases.json maps
--- unversioned / GI input to the versioned file. Files older than cache_days are
--- ignored (re-fetched). The API key never appears in file names, fixture keys
--- or M.history.
---
--- Requests go through one queue: starts are spaced 350 ms apart (110 ms with an
--- API key from config or $NCBI_API_KEY), i.e. under NCBI's 3/s (10/s) limit. HTTP 429, 5xx and network errors are
--- retried `retries` times with exponential backoff (retry_base_ms * 2^n), as are
--- the transient "Search Backend failed" esearch answers (HTTP 200 + ERROR),
--- honouring a numeric Retry-After header. Timeouts are not retried.
--- Long id lists (> post_ids ids or a URL > post_url_len chars) are POSTed.
local M = {}

M.defaults = {
  base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
  api_key = nil,          -- fallback: $NCBI_API_KEY
  email = nil,            -- sent only when the user sets it
  tool = "plasmid.nvim",
  max_length = 5000000,   -- bp/aa above which fetch refuses unless opts.force
  timeout_ms = 60000,
  retries = 3,
  cache = true,
  cache_dir = nil,        -- nil = stdpath("cache") .. "/plasmid/ncbi"
  cache_days = 30,
  -- extras
  min_interval_ms = nil,  -- nil = 350 (no key) / 110 (with key); NCBI allows 3/s and 10/s.
                          -- The margin over 334/100 ms absorbs network jitter: at exactly
                          -- 334 ms NCBI still answered 2 of 38 recording requests with 429.
  retry_base_ms = 1000,
  retry_max_ms = 60000,
  post_ids = 200,
  post_url_len = 2000,
}

--------------------------------------------------------------------------------
-- options
--------------------------------------------------------------------------------

function M.config()
  local ok, cfgmod = pcall(require, "plasmid.config")
  local user = ok and cfgmod.options and cfgmod.options.ncbi or nil
  local o = vim.tbl_deep_extend("force", vim.deepcopy(M.defaults), type(user) == "table" and user or {})
  if o.api_key == nil or o.api_key == "" then
    local env = os.getenv("NCBI_API_KEY")
    o.api_key = (env and env ~= "") and env or nil
  end
  if o.email == "" then o.email = nil end
  if not o.cache_dir then o.cache_dir = vim.fn.stdpath("cache") .. "/plasmid/ncbi" end
  if not o.base_url:match("/$") then o.base_url = o.base_url .. "/" end
  if not o.min_interval_ms then o.min_interval_ms = o.api_key and 110 or 350 end
  return o
end

--------------------------------------------------------------------------------
-- small helpers
--------------------------------------------------------------------------------

local function trim(s) return (s:gsub("^%s+", ""):gsub("%s+$", "")) end

local function commas(n)
  local s = tostring(math.floor(n))
  local out = s:reverse():gsub("(%d%d%d)", "%1,"):reverse()
  return (out:gsub("^,", ""))
end
M._commas = commas

local function enc(s)
  return (tostring(s):gsub("[^%w%-%._~]", function(c) return ("%%%02X"):format(c:byte()) end))
end

local function json_decode(body)
  if type(body) ~= "string" or body == "" then return nil end
  local ok, v = pcall(vim.json.decode, body, { luanil = { object = true, array = true } })
  if ok and type(v) == "table" then return v end
  return nil
end

--- NCBI writes some efetch errors with every letter spaced out
--- ("F a i l e d  t o  u n d e r s t a n d  i d :  X X 9"). Undo that.
local function despace(s)
  s = trim(s)
  if s:match("^%S %S %S ") then
    local words = {}
    for w in (s .. "  "):gmatch("(.-)%s%s+") do
      if w ~= "" then words[#words + 1] = (w:gsub(" ", "")) end
    end
    s = table.concat(words, " "):gsub(" :", ":")
  end
  return s
end
M._despace = despace

--- Remove api_key / email parameters from a query string or form body.
local function scrub(q)
  if not q or q == "" then return q end
  local keep = {}
  for part in q:gmatch("[^&]+") do
    local k = part:match("^([^=]*)")
    if k ~= "api_key" and k ~= "email" then keep[#keep + 1] = part end
  end
  return table.concat(keep, "&")
end

local function scrub_url(url)
  local base, q = url:match("^([^?]*)%?(.*)$")
  if not base then return url end
  local s = scrub(q)
  return s ~= "" and (base .. "?" .. s) or base
end
M._scrub_url = scrub_url

function M.request_key(req)
  local material = (req.method or "GET") .. " " .. scrub_url(req.url) .. "\n" .. (scrub(req.body) or "")
  return vim.fn.sha256(material):sub(1, 24)
end

local function after(ms, fn)
  local t = vim.uv.new_timer()
  local live = true
  t:start(math.max(0, math.floor(ms)), 0, function()
    if not t:is_closing() then t:close() end
    vim.schedule(function()
      if live then live = false; fn() end
    end)
  end)
  return {
    cancel = function()
      live = false
      if not t:is_closing() then t:stop(); t:close() end
    end,
  }
end

local function now_ms() return vim.uv.hrtime() / 1e6 end

--------------------------------------------------------------------------------
-- accession parsing
--------------------------------------------------------------------------------

local REFSEQ_NUC = { AC = true, NC = true, NG = true, NM = true, NR = true, NT = true, NW = true, NZ = true, XM = true, XR = true }
local REFSEQ_PROT = { NP = true, XP = true, WP = true, YP = true, AP = true }

local function classify(idu)
  local acc = idu:match("^(.-)%.%d+$") or idu
  if idu:match("^%d+$") then return "nuccore" end
  local pre, rest = acc:match("^(%u%u)_(%w+)$")
  if pre then
    if REFSEQ_PROT[pre] and rest:match("^%d+$") then return "protein" end
    if REFSEQ_NUC[pre] and (rest:match("^%d+$") or rest:match("^%u%u%d%d%d%d%d%d$") or rest:match("^%u%u%d%d%d%d%d%d%d%d$")
      or rest:match("^%u%u%u%u%d%d%d%d%d%d%d%d%d?$") or rest:match("^%u%u%u%u%u%u%d%d%d%d%d%d%d%d%d%d?%d?$")) then
      return "nuccore"
    end
    return nil, ("unknown RefSeq prefix %s_"):format(pre)
  end
  -- UniProt (checked before the 1+5 GenBank pattern: O/P/Q are UniProt-only)
  if acc:match("^[OPQ]%d%w%w%w%d$") then return "protein" end
  if acc:match("^[A-NR-Z]%d%u%w%w%d$") or acc:match("^[A-NR-Z]%d%u%w%w%d%u%w%w%d$") then return "protein" end
  if acc:match("^%u%d%d%d%d%d$") then return "nuccore" end
  if acc:match("^%u%u%d%d%d%d%d%d$") or acc:match("^%u%u%d%d%d%d%d%d%d%d$") then return "nuccore" end
  if acc:match("^%u%u%u%d%d%d%d%d$") or acc:match("^%u%u%u%d%d%d%d%d%d%d$") then return "protein" end
  if acc:match("^%u%u%u%u%d%d%d%d%d%d%d%d%d?$") then return "nuccore" end           -- WGS 4+8/9
  if acc:match("^%u%u%u%u%u%u%d%d%d%d%d%d%d%d%d%d?%d?$") then return "nuccore" end -- WGS 6+9..11
  if acc:match("^%u%u%u%u%u%d%d%d%d%d%d%d$") then return "nuccore" end               -- MGA 5+7
  return nil
end

local FORMS = "expected e.g. L09137, L09137.2, NM_000518.5, NC_000913.3:190-2800, NC_000913.3:2800-190, "
  .. "NC_000913.3:190..2800:-, NP_000509.1, AAA12345, P69905, protein:ID or a GI number"

local function num(s)
  if not s then return nil end
  local d = s:gsub(",", "")
  if not d:match("^%d+$") then return nil end
  return tonumber(d)
end

function M.parse_accession(text)
  if type(text) ~= "string" then return nil, "accession must be a string" end
  local s = trim(text)
  if s == "" then return nil, "empty accession (" .. FORMS .. ")" end
  local forced
  local p, rest = s:match("^(%a+):%s*(.*)$")
  if p then
    local pl = p:lower()
    if pl == "protein" then forced = "protein"
    elseif pl == "nuccore" or pl == "nucleotide" then forced = "nuccore" end
    if forced then s = rest end
  end
  local id, tail = s:match("^([%w_%.]+)(.*)$")
  if not id then return nil, ("not an accession: %q (%s)"):format(text, FORMS) end
  tail = trim(tail)
  if not forced then
    local okdb, why0 = classify(id:upper())
    if not okdb then return nil, ("not an accession: %q%s (%s)"):format(text, why0 and (": " .. why0) or "", FORMS) end
  end
  local from, to, strand
  if tail ~= "" then
    local t = tail
    -- strand suffix
    local function take_strand(x)
      local body, st = x:match("^(.-)%s*:%s*([%-%+]?)%s*$")
      if body and st ~= "" then return body, st == "-" and 2 or 1 end
      local b2, word = x:match("^(.-)%s*:%s*(%a+)%s*$")
      if b2 then
        local w = word:lower()
        if w == "minus" or w == "rc" then return b2, 2 end
        if w == "plus" then return b2, 1 end
      end
      return x, nil
    end
    t, strand = take_strand(t)
    t = trim(t)
    if t ~= "" then
      local a, b = t:match("^:%s*([%d,]+)%s*%-%s*([%d,]+)$")
      if not a then a, b = t:match("^:%s*([%d,]+)%s*%.%.%s*([%d,]+)$") end
      local na, nb = num(a), num(b)
      if not t:match("^:") then return nil, ("not an accession: %q (%s)"):format(text, FORMS) end
      if not (na and nb) then
        return nil, ("cannot read the range in %q (use FROM-TO or FROM..TO, e.g. NC_000913.3:190-2800)"):format(text)
      end
      if na < 1 or nb < 1 then return nil, ("range positions are 1-based (got %q)"):format(text) end
      if na > nb then
        if strand == 1 then
          return nil, ("%q: a reversed range means the minus strand but ':+' was given"):format(text)
        end
        na, nb, strand = nb, na, 2
      end
      from, to = na, nb
    end
  end
  local idu = id:upper()
  local db, why = classify(idu)
  if forced then
    if not idu:match("^[%w_]+%.?%d*$") then return nil, ("not an accession: %q (%s)"):format(text, FORMS) end
    db = forced
  elseif not db then
    return nil, ("not an accession: %q%s (%s)"):format(text, why and (": " .. why) or "", FORMS)
  end
  if db == "protein" and strand == 2 then
    return nil, ("%q: the reverse strand only applies to nucleotide records"):format(text)
  end
  return { id = id, db = db, from = from, to = to, strand = strand }
end

--------------------------------------------------------------------------------
-- transports
--------------------------------------------------------------------------------

M.transports = {}

local function parse_headers(raw)
  -- the file may hold several blocks (proxy CONNECT, redirects): keep the last
  local blocks, cur = {}, nil
  for line in (raw or ""):gmatch("[^\r\n]*") do
    local status = line:match("^HTTP/[%d%.]+%s+(%d+)")
    if status then
      cur = { status = tonumber(status), headers = {} }
      blocks[#blocks + 1] = cur
    elseif cur and line ~= "" then
      local k, v = line:match("^([^:]+):%s*(.*)$")
      if k then cur.headers[k:lower()] = v end
    end
  end
  return blocks[#blocks]
end
M._parse_headers = parse_headers

local function cfgquote(s) return '"' .. s:gsub("\\", "\\\\"):gsub('"', '\\"') .. '"' end

--- Default transport: curl via vim.system. The URL (and POST body) are passed on
--- stdin as a curl config so an API key does not show up in the process list.
--- `done` is called on the main loop.
function M.transports.curl(req, done)
  if vim.fn.executable("curl") ~= 1 then
    done("curl is not installed (not found in $PATH); plasmid.nvim needs it to talk to NCBI")
    return
  end
  local hdr = vim.fn.tempname()
  local timeout_s = math.max(1, math.ceil((req.timeout_ms or 60000) / 1000))
  local conf = { "url = " .. cfgquote(req.url) }
  if req.method == "POST" then
    conf[#conf + 1] = "data-binary = " .. cfgquote(req.body or "")
    conf[#conf + 1] = "header = " .. cfgquote("Content-Type: application/x-www-form-urlencoded")
  end
  local args = { "curl", "-sS", "-L", "--compressed", "--max-time", tostring(timeout_s), "-D", hdr, "-K", "-" }
  local killed = false
  local fast_done = done
  done = function(...) local a = { n = select("#", ...), ... }; vim.schedule(function() fast_done(unpack(a, 1, a.n)) end) end
  local ok, proc = pcall(vim.system, args, {
    stdin = table.concat(conf, "\n") .. "\n",
    timeout = (timeout_s + 5) * 1000,
  }, function(obj)
    local raw
    local fh = io.open(hdr, "rb")
    if fh then raw = fh:read("*a"); fh:close() end
    os.remove(hdr)
    if killed then return end
    if obj.code == 28 or obj.code == 124 then
      done(("timed out after %d s (timeout_ms = %d)"):format(timeout_s, req.timeout_ms or 60000))
    elseif obj.code ~= 0 then
      local msg = trim((obj.stderr or ""):gsub("^curl: %(%d+%)%s*", ""))
      done(("network error: %s (curl exit %d)"):format(msg ~= "" and msg or "no details", obj.code))
    else
      local h = parse_headers(raw)
      if not h then done("network error: no HTTP response headers from curl"); return end
      done(nil, { status = h.status, headers = h.headers, body = obj.stdout or "" })
    end
  end)
  if not ok then
    os.remove(hdr)
    done("curl could not be started: " .. tostring(proc))
    return
  end
  return {
    cancel = function()
      killed = true
      pcall(function() proc:kill(15) end)
    end,
  }
end

local function read_file(path)
  local fh = io.open(path, "rb")
  if not fh then return nil end
  local s = fh:read("*a")
  fh:close()
  return s
end

local function write_file(path, data)
  vim.fn.mkdir(vim.fn.fnamemodify(path, ":h"), "p")
  local tmp = path .. ".tmp" .. tostring(vim.uv.hrtime())
  local fh = assert(io.open(tmp, "wb"))
  fh:write(data)
  fh:close()
  assert(os.rename(tmp, path))
end

--- Serve recorded responses: dir/<key>.json (meta) + dir/<key>.body (raw body).
function M.transports.replay(dir)
  return function(req, done)
    local key = M.request_key(req)
    local meta = json_decode(read_file(dir .. "/" .. key .. ".json"))
    if not meta then
      done(("replay: no recorded response for key %s (%s %s%s) in %s"):format(
        key, req.method or "GET", scrub_url(req.url), req.body and (" body " .. scrub(req.body)) or "", dir))
      return
    end
    local body = read_file(dir .. "/" .. key .. ".body") or ""
    done(nil, { status = meta.status, headers = meta.headers or {}, body = body })
  end
end

--- Pass requests to `inner` (default curl) and save every HTTP response except
--- 429/5xx (those are retried; the final answer overwrites the same key).
function M.transports.record(dir, inner)
  inner = inner or M.transports.curl
  return function(req, done)
    local t0 = now_ms()
    return inner(req, function(err, resp)
      if not err and resp and resp.status ~= 429 and resp.status < 500 then
        local key = M.request_key(req)
        local keep = {}
        for _, k in ipairs({ "content-type", "retry-after", "x-ratelimit-limit" }) do
          keep[k] = resp.headers and resp.headers[k] or nil
        end
        write_file(dir .. "/" .. key .. ".body", resp.body or "")
        write_file(dir .. "/" .. key .. ".json", vim.json.encode({
          key = key, method = req.method or "GET", url = scrub_url(req.url), body = scrub(req.body),
          status = resp.status, headers = keep, bytes = #(resp.body or ""),
          ms = math.floor(now_ms() - t0 + 0.5), recorded = os.date("!%Y-%m-%dT%H:%M:%SZ"),
        }) .. "\n")
      end
      done(err, resp)
    end)
  end
end

--------------------------------------------------------------------------------
-- request queue (rate limiting) and retries
--------------------------------------------------------------------------------

local state = { transport = nil, queue = {}, last_start = nil, pump_timer = nil }
M.history = {}   -- { url (scrubbed), method, status|err, ms, attempt } most recent last (max 100)

function M.set_transport(fn) state.transport = fn end
function M.get_transport() return state.transport or M.transports.curl end

local pump

local function schedule_pump(ms)
  if state.pump_timer then return end
  state.pump_timer = after(ms, function()
    state.pump_timer = nil
    pump()
  end)
end

pump = function()
  if #state.queue == 0 then return end
  local job = state.queue[1]
  local t = now_ms()
  if state.last_start then
    local wait = state.last_start + job.interval - t
    if wait > 0 then schedule_pump(wait); return end
  end
  table.remove(state.queue, 1)
  state.last_start = t
  job.started = true
  local t0 = t
  local finished = false
  local ok, canceller = pcall(M.get_transport(), job.req, function(err, resp)
    if finished then return end
    finished = true
    vim.schedule(function()
      if job.cancelled then return end
      local hist = M.history
      hist[#hist + 1] = { method = job.req.method, url = scrub_url(job.req.url), status = resp and resp.status,
        err = err, ms = math.floor(now_ms() - t0 + 0.5), bytes = resp and resp.body and #resp.body or 0 }
      if #hist > 100 then table.remove(hist, 1) end
      job.cb(err, resp)
    end)
  end)
  if not ok then
    finished = true
    vim.schedule(function() if not job.cancelled then job.cb("transport error: " .. tostring(canceller)) end end)
  else
    job.canceller = canceller
  end
  if #state.queue > 0 then schedule_pump(job.interval) end
end

local function enqueue(req, interval, cb)
  local job = { req = req, interval = interval, cb = cb }
  function job.cancel()
    job.cancelled = true
    if not job.started then
      for i, j in ipairs(state.queue) do
        if j == job then table.remove(state.queue, i); break end
      end
    elseif job.canceller then
      local c = job.canceller
      if type(c) == "function" then pcall(c)
      elseif type(c) == "table" and c.cancel then pcall(c.cancel, c)
      elseif type(c) == "table" and c.kill then pcall(c.kill, c, 15) end
    end
  end
  state.queue[#state.queue + 1] = job
  vim.schedule(pump)
  return job
end

--- Test hook: forget when the last request started (stub-transport tests only).
function M._forget_rate_clock() state.last_start = nil end

--- Restore the default transport, drop queued requests (their callbacks never
--- run) and clear M.history.
function M.reset()
  for _, j in ipairs(state.queue) do j.cancelled = true end
  state.queue = {}
  if state.pump_timer then state.pump_timer.cancel(); state.pump_timer = nil end
  -- state.last_start is kept on purpose: a reset must not let the next request
  -- jump the rate limit (tests that reset between live calls hit 429 otherwise)
  state.transport = nil
  M.history = {}
end

--------------------------------------------------------------------------------
-- handles
--------------------------------------------------------------------------------

local Handle = {}
Handle.__index = Handle

local function new_handle(cb)
  return setmetatable({ cb = cb, subs = {}, done = false, called = false, cancelled = false }, Handle)
end

function Handle:_track(sub)
  self.subs[#self.subs + 1] = sub
  return sub
end

function Handle:_finish(err, res, extra)
  if self.done or self.cancelled then return end
  self.done = true
  self.err, self.result, self.extra = err, res, extra
  vim.schedule(function()
    if self.cancelled then return end
    self.called = true
    if self.cb then self.cb(err, res, extra) end
  end)
end

--- Cancel: kills curl / drops queued requests / stops backoff timers. The
--- callback is never called after cancel(). Returns true when a callback was
--- still pending.
function Handle:cancel()
  if self.called or self.cancelled then return false end
  self.cancelled = true
  for _, s in ipairs(self.subs) do pcall(s.cancel) end
  return true
end

function M.wait(handle, timeout_ms)
  local ok = vim.wait(timeout_ms or 60000, function() return handle.called or handle.cancelled end, 5)
  return ok and handle.called, handle.err, handle.result
end

--------------------------------------------------------------------------------
-- NCBI error interpretation
--------------------------------------------------------------------------------

local HTTP_REASON = { [400] = "Bad Request", [403] = "Forbidden", [404] = "Not Found", [414] = "URI Too Long",
  [429] = "Too Many Requests", [500] = "Internal Server Error", [502] = "Bad Gateway",
  [503] = "Service Unavailable", [504] = "Gateway Timeout" }

--- Error text embedded in an E-utilities body (JSON, XML or plain text), or nil.
local function body_error(body)
  if type(body) ~= "string" or body == "" then return nil end
  local head = body:sub(1, 4000)
  if head:match("^%s*{") then
    local j = json_decode(body)
    if j then
      if type(j.error) == "string" then return j.error end
      if type(j.ERROR) == "string" then return j.ERROR end
      local es = j.esearchresult
      if es and type(es.ERROR) == "string" then return es.ERROR end
      if type(j.esummaryresult) == "table" and type(j.esummaryresult[1]) == "string" then return j.esummaryresult[1] end
    end
    return nil
  end
  local x = head:match("<ERROR>(.-)</ERROR>") or head:match("<Error>(.-)</Error>")
  if x then return despace(x) end
  local first = head:match("^%s*([^\n]*)")
  if first and first:match("^Error") then
    local msg = despace(first:gsub("^Error%s*:?%s*", ""))
    return msg
  end
  return nil
end
M._body_error = body_error

local function http_error(what, resp)
  local detail = body_error(resp.body)
  if not detail and resp.body and not resp.body:match("^%s*<") then
    local first = trim(resp.body:sub(1, 200):match("^[^\n]*") or "")
    if first ~= "" then detail = first end
  end
  local s = ("NCBI %s failed: HTTP %d%s"):format(what, resp.status, HTTP_REASON[resp.status] and (" " .. HTTP_REASON[resp.status]) or "")
  if resp.status == 429 then s = s .. " (rate limit; set an API key: ncbi.api_key or $NCBI_API_KEY)" end
  if detail then s = s .. ": " .. detail end
  return s
end

local function retryable(err, resp)
  if err then return err:match("^network error") ~= nil end
  if resp.status == 429 or resp.status >= 500 then return true end
  -- older E-utilities answer rate-limit violations with 200 + JSON error
  local head = resp.body and resp.body:sub(1, 600) or ""
  if head:match('"error"%s*:%s*"API rate limit exceeded"') then return true end
  -- transient esearch failure reported as HTTP 200 + {"esearchresult":{"ERROR":"Search Backend failed: ..."}}
  if head:match('"ERROR"%s*:%s*"Search Backend failed') or head:match("<ERROR>Search Backend failed") then return true end
  return false
end

--------------------------------------------------------------------------------
-- E-utility requests
--------------------------------------------------------------------------------

local function build_query(params, cfg)
  local parts = {}
  for _, kv in ipairs(params) do
    if kv[2] ~= nil then parts[#parts + 1] = kv[1] .. "=" .. enc(kv[2]) end
  end
  parts[#parts + 1] = "tool=" .. enc(cfg.tool)
  if cfg.email then parts[#parts + 1] = "email=" .. enc(cfg.email) end
  if cfg.api_key then parts[#parts + 1] = "api_key=" .. enc(cfg.api_key) end
  return table.concat(parts, "&")
end

--- Low-level: run one E-utility (name = "efetch" | "esummary" | "esearch" | ...)
--- with ordered params {{k, v}, ...}. opts = { post = bool }. cb(err, resp) with
--- resp = { status, body, headers }; HTTP errors and NCBI error bodies become err.
local function eutil(h, cfg, name, params, opts, cb)
  local q = build_query(params, cfg)
  local url = cfg.base_url .. name .. ".fcgi"
  local req
  if opts.post or #url + 1 + #q > cfg.post_url_len then
    req = { url = url, method = "POST", body = q, timeout_ms = cfg.timeout_ms }
  else
    req = { url = url .. "?" .. q, method = "GET", timeout_ms = cfg.timeout_ms }
  end
  local attempt = 0
  local function try()
    attempt = attempt + 1
    h:_track(enqueue(req, cfg.min_interval_ms, function(err, resp)
      if h.cancelled then return end
      if attempt <= cfg.retries and retryable(err, resp) then
        local wait = cfg.retry_base_ms * 2 ^ (attempt - 1)
        local ra = resp and resp.headers and tonumber(resp.headers["retry-after"])
        if ra then wait = math.max(wait, ra * 1000) end
        wait = math.min(wait, cfg.retry_max_ms)
        h:_track(after(wait, try))
        return
      end
      if err then cb(("NCBI %s: %s"):format(name, err)); return end
      if resp.status < 200 or resp.status >= 300 then
        cb(http_error(name, resp), resp)
        return
      end
      cb(nil, resp)
    end))
  end
  try()
end

function M.eutil(name, params, opts, cb)
  local h = new_handle(cb)
  eutil(h, M.config(), name, params, opts or {}, function(err, resp) h:_finish(err, resp) end)
  return h
end

--------------------------------------------------------------------------------
-- esummary / esearch
--------------------------------------------------------------------------------

local function seq_item(d)
  return {
    uid = d.uid,
    accession = d.accessionversion ~= "" and d.accessionversion or d.caption,
    title = d.title,
    length = tonumber(d.slen),
    topology = d.topology,
    moltype = d.moltype,
    organism = d.organism,
    updated = d.updatedate,
    created = d.createdate,
    caption = d.caption,
    taxid = d.taxid,
    sourcedb = d.sourcedb,
    biomol = d.biomol,
    completeness = d.completeness,
  }
end

--- Map esummary results back to the requested ids (uid, accession.version or
--- unversioned caption); leftovers are matched in order.
local function order_items(ids, result)
  local items, by = {}, {}
  for _, uid in ipairs(result.uids or {}) do
    local d = result[uid]
    if type(d) == "table" and not d.error then
      local it = seq_item(d)
      items[#items + 1] = it
      by[tostring(uid)] = it
      if it.accession then by[it.accession:upper()] = by[it.accession:upper()] or it end
      if d.caption then by[d.caption:upper()] = by[d.caption:upper()] or it end
    end
  end
  local used, out, missing, pending = {}, {}, {}, {}
  for i, id in ipairs(ids) do
    local key = tostring(id):upper()
    local it = by[key]
    if it and not used[it] then used[it] = true; out[i] = it else pending[#pending + 1] = i end
  end
  local spare = {}
  for _, it in ipairs(items) do if not used[it] then spare[#spare + 1] = it end end
  for _, i in ipairs(pending) do
    local it = table.remove(spare, 1)
    if it then out[i] = it else missing[#missing + 1] = ids[i] end
  end
  local list = {}
  for i = 1, #ids do if out[i] then list[#list + 1] = out[i] end end
  return list, missing
end

local function summary_into(h, cfg, db, ids, cb)
  if #ids == 0 then cb(nil, {}, {}); return end
  local strs = vim.tbl_map(tostring, ids)
  eutil(h, cfg, "esummary", { { "db", db }, { "id", table.concat(strs, ",") }, { "retmode", "json" } },
    { post = #strs > cfg.post_ids }, function(err, resp)
      if err then cb(err); return end
      local j = json_decode(resp.body)
      if not j or type(j.result) ~= "table" then
        local e = body_error(resp.body)
        cb(e and ("NCBI esummary: " .. e) or "NCBI esummary: unreadable response (not JSON)")
        return
      end
      local list, missing = order_items(strs, j.result)
      if #list == 0 then
        local e = body_error(resp.body)
        local what = #strs == 1 and ("%s accession '%s'"):format(db, strs[1]) or (db .. " ids " .. table.concat(strs, ", "))
        cb(("NCBI does not know %s%s"):format(what, e and (" (NCBI: " .. e .. ")") or ""))
        return
      end
      cb(nil, list, missing)
    end)
end

--- cb(err, items, missing): items in the order given; ids NCBI does not know are
--- left out of `items` and listed in `missing` (all unknown -> err).
function M.summary(db, ids, cb)
  local h = new_handle(cb)
  if type(ids) ~= "table" then ids = { ids } end
  summary_into(h, M.config(), db, ids, function(err, list, missing) h:_finish(err, list, missing) end)
  return h
end

function M.search(db, term, opts, cb)
  if type(opts) == "function" then cb, opts = opts, nil end
  opts = opts or {}
  local h = new_handle(cb)
  local cfg = M.config()
  eutil(h, cfg, "esearch", { { "db", db }, { "term", term }, { "retmax", opts.retmax or 20 },
    { "retstart", opts.retstart }, { "sort", opts.sort }, { "retmode", "json" } }, {}, function(err, resp)
      if err then h:_finish(err); return end
      local j = json_decode(resp.body)
      local es = j and j.esearchresult
      if not es then
        local e = body_error(resp.body)
        h:_finish(e and ("NCBI esearch: " .. e) or "NCBI esearch: unreadable response")
        return
      end
      if type(es.ERROR) == "string" then h:_finish("NCBI esearch: " .. trim(es.ERROR)); return end
      local warnings = {}
      for _, lst in ipairs({ es.errorlist or {}, es.warninglist or {} }) do
        for k, v in pairs(lst) do
          if type(v) == "table" and #v > 0 then warnings[#warnings + 1] = k .. ": " .. table.concat(v, ", ") end
        end
      end
      table.sort(warnings)
      local ids = es.idlist or {}
      local count = tonumber(es.count) or 0
      if #ids == 0 then h:_finish(nil, { count = count, items = {}, warnings = warnings }); return end
      summary_into(h, cfg, db, ids, function(err2, list)
        if err2 then h:_finish(err2); return end
        h:_finish(nil, { count = count, items = list, warnings = warnings, ids = ids })
      end)
    end)
  return h
end

--------------------------------------------------------------------------------
-- gene
--------------------------------------------------------------------------------

--- Convert one esummary genomicinfo entry (0-based chrstart/chrstop; stop <
--- start on the minus strand) to 1-based inclusive start <= stop + strand.
function M.gene_location(gi)
  local a, b = tonumber(gi.chrstart), tonumber(gi.chrstop)
  if not (a and b) then return nil end
  if b >= a then return a + 1, b + 1, 1 end
  return b + 1, a + 1, 2
end

local function gene_term(sym, organism)
  local t = trim(sym)
  if t:match("^[%w%-%._/@']+$") then t = t .. "[sym]" end
  if organism and trim(organism) ~= "" then t = ("(%s) AND \"%s\"[orgn]"):format(t, trim(organism)) end
  return t .. " AND alive[prop]"
end
M._gene_term = gene_term

--- cb(err, genes). Genes without a current genomic annotation are skipped.
function M.gene(symbol, organism, opts, cb)
  if type(opts) == "function" then cb, opts = opts, nil end
  opts = opts or {}
  local h = new_handle(cb)
  local cfg = M.config()
  if type(symbol) ~= "string" or trim(symbol) == "" then h:_finish("gene: empty symbol"); return h end
  local term = gene_term(symbol, organism)
  eutil(h, cfg, "esearch", { { "db", "gene" }, { "term", term }, { "retmax", opts.retmax or 20 }, { "retmode", "json" } },
    {}, function(err, resp)
      if err then h:_finish(err); return end
      local j = json_decode(resp.body)
      local es = j and j.esearchresult
      if not es then h:_finish("NCBI esearch (gene): unreadable response"); return end
      if type(es.ERROR) == "string" then h:_finish("NCBI esearch (gene): " .. trim(es.ERROR)); return end
      local ids = es.idlist or {}
      if #ids == 0 then
        h:_finish(("no current gene found for %s%s"):format(symbol, organism and (" in " .. organism) or ""))
        return
      end
      eutil(h, cfg, "esummary", { { "db", "gene" }, { "id", table.concat(ids, ",") }, { "retmode", "json" } },
        { post = #ids > cfg.post_ids }, function(err2, resp2)
          if err2 then h:_finish(err2); return end
          local j2 = json_decode(resp2.body)
          local r = j2 and j2.result
          if not r then h:_finish("NCBI esummary (gene): " .. (body_error(resp2.body) or "unreadable response")); return end
          local out = {}
          for _, uid in ipairs(r.uids or {}) do
            local d = r[uid]
            if type(d) == "table" and not d.error and type(d.genomicinfo) == "table" and d.genomicinfo[1] then
              local gi = d.genomicinfo[1]
              local s, e, strand = M.gene_location(gi)
              if s then
                local locs = {}
                for _, g in ipairs(d.genomicinfo) do
                  local s2, e2, st2 = M.gene_location(g)
                  locs[#locs + 1] = { chr_accession = g.chraccver, start = s2, stop = e2, strand = st2 }
                end
                out[#out + 1] = {
                  geneid = tostring(d.uid),
                  symbol = d.name,
                  description = d.description,
                  organism = type(d.organism) == "table" and d.organism.scientificname or d.organism,
                  taxid = type(d.organism) == "table" and d.organism.taxid or nil,
                  chr_accession = gi.chraccver,
                  start = s, stop = e, strand = strand,
                  location_text = d.maplocation ~= "" and d.maplocation or nil,
                  chromosome = d.chromosome ~= "" and d.chromosome or nil,
                  exon_count = tonumber(gi.exoncount),
                  aliases = d.otheraliases ~= "" and d.otheraliases or nil,
                  locations = locs,
                }
              end
            end
          end
          if #out == 0 then
            h:_finish(("no annotated genomic location for %s%s"):format(symbol, organism and (" in " .. organism) or ""))
            return
          end
          h:_finish(nil, out)
        end)
    end)
  return h
end

--------------------------------------------------------------------------------
-- cache
--------------------------------------------------------------------------------

local function safe(s) return (s:upper():gsub("[^%w%._%-]", "_")) end
local function ext(db) return db == "protein" and ".gp" or ".gb" end
local function suffix(from, to, strand)
  local s = ""
  if from then s = ("_%d-%d"):format(from, to) end
  if strand == 2 then s = s .. "_rc" end
  return s
end
local function versioned(id) return id:match("^.+%.%d+$") ~= nil and not id:match("^%d+$") end

local function alias_key(spec)
  return spec.db .. ":" .. spec.id:upper() .. suffix(spec.from, spec.to, spec.strand)
end

local function read_aliases(cfg)
  return json_decode(read_file(cfg.cache_dir .. "/aliases.json")) or {}
end

--- Path of the cache file for spec. With accession_version (from a summary) or a
--- versioned spec.id the path is computed; for unversioned/GI input it is looked
--- up in aliases.json (nil when that input was never fetched).
function M.cache_path(spec, accession_version, cfg)
  cfg = cfg or M.config()
  if type(spec) == "string" then spec = assert(M.parse_accession(spec)) end
  local acc = accession_version or (versioned(spec.id) and spec.id or nil)
  if acc then
    return cfg.cache_dir .. "/" .. safe(acc) .. suffix(spec.from, spec.to, spec.strand) .. ext(spec.db)
  end
  local name = read_aliases(cfg)[alias_key(spec)]
  return name and (cfg.cache_dir .. "/" .. name) or nil
end

local function fresh(path, cfg)
  local st = vim.uv.fs_stat(path)
  if not st or st.type ~= "file" then return false end
  return (os.time() - st.mtime.sec) <= cfg.cache_days * 86400
end

local function meta_from_text(text)
  local m = {}
  m.accession = text:match("\nVERSION%s+(%S+)")
  local def = text:match("\nDEFINITION%s+(.-)\n%S")
  if def then m.title = (def:gsub("%s*\n%s+", " "):gsub("%.$", "")) end
  m.length = tonumber(text:match("^LOCUS%s+%S+%s+(%d+)"))
  return m
end

local function cache_store(cfg, spec, path, text, meta)
  write_file(path, text)
  write_file(path .. ".json", vim.json.encode(meta) .. "\n")
  if not versioned(spec.id) or spec.id:upper() ~= (meta.accession or ""):upper() then
    local al = read_aliases(cfg)
    al[alias_key(spec)] = vim.fn.fnamemodify(path, ":t")
    write_file(cfg.cache_dir .. "/aliases.json", vim.json.encode(al) .. "\n")
  end
end

local function from_cache(path, spec)
  local text = read_file(path)
  if not text then return nil end
  local meta = json_decode(read_file(path .. ".json")) or meta_from_text(text)
  local res = vim.tbl_extend("force", {}, meta)
  res.text, res.path, res.from_cache = text, path, true
  res.db = res.db or spec.db
  res.id = spec.id
  return res
end

function M.cache_list()
  local cfg = M.config()
  local out = {}
  for _, pat in ipairs({ "*.gb", "*.gp" }) do
    for _, p in ipairs(vim.fn.globpath(cfg.cache_dir, pat, false, true)) do
      local st = vim.uv.fs_stat(p)
      local meta = json_decode(read_file(p .. ".json")) or {}
      out[#out + 1] = { path = p, accession = meta.accession, size = st and st.size or 0,
        mtime = st and st.mtime.sec or 0, db = meta.db, length = meta.length, from = meta.from, to = meta.to,
        strand = meta.strand, title = meta.title }
    end
  end
  table.sort(out, function(a, b) return a.path < b.path end)
  return out
end

--- Remove every cache file; returns the number of cached records removed.
function M.cache_clear()
  local cfg = M.config()
  local n = 0
  for _, it in ipairs(M.cache_list()) do
    if os.remove(it.path) then n = n + 1 end
    os.remove(it.path .. ".json")
  end
  os.remove(cfg.cache_dir .. "/aliases.json")
  return n
end

--------------------------------------------------------------------------------
-- fetch
--------------------------------------------------------------------------------

function M.fetch(spec, opts, cb)
  if type(opts) == "function" then cb, opts = opts, nil end
  opts = opts or {}
  local h = new_handle(cb)
  if type(spec) == "string" then
    local s, perr = M.parse_accession(spec)
    if not s then h:_finish(perr); return h end
    spec = s
  end
  local cfg = M.config()
  if cfg.cache and not opts.refresh then
    local p = M.cache_path(spec, nil, cfg)
    if p and fresh(p, cfg) then
      local res = from_cache(p, spec)
      if res then h:_finish(nil, res); return h end
    end
  end
  summary_into(h, cfg, spec.db, { spec.id }, function(err, items)
    if err then h:_finish(err); return end
    local it = items[1]
    local total = it.length or 0
    local from, to = spec.from, spec.to
    if from then
      if from > total then
        h:_finish(("%s is %s %s long; the range %s-%s starts beyond its end"):format(
          it.accession, commas(total), spec.db == "protein" and "aa" or "bp", commas(from), commas(spec.to)))
        return
      end
      if to > total then to = total end
    end
    local len = from and (to - from + 1) or total
    local unit = spec.db == "protein" and "aa" or "bp"
    if len > cfg.max_length and not opts.force then
      h:_finish(("%s%s is %s %s, above the download limit of %s %s (ncbi.max_length); "
        .. "fetch a range (e.g. %s:1-100000) or force the download"):format(
        it.accession, from and (" range " .. from .. "-" .. to) or "", commas(len), unit,
        commas(cfg.max_length), unit, it.accession))
      return
    end
    local cspec = { id = spec.id, db = spec.db, from = from, to = to, strand = spec.strand }
    local path = cfg.cache and M.cache_path(cspec, it.accession, cfg) or nil
    local meta = { accession = it.accession, title = it.title, length = len, db = spec.db, uid = it.uid,
      total_length = total, from = from, to = to, strand = spec.strand, topology = it.topology,
      organism = it.organism, updated = it.updated }
    if path and not opts.refresh and fresh(path, cfg) then
      local res = from_cache(path, spec)
      if res then
        if not versioned(spec.id) or spec.id:upper() ~= (it.accession or ""):upper() then
          local al = read_aliases(cfg)
          al[alias_key(spec)] = vim.fn.fnamemodify(path, ":t")
          write_file(cfg.cache_dir .. "/aliases.json", vim.json.encode(al) .. "\n")
        end
        h:_finish(nil, res)
        return
      end
    end
    local params = { { "db", spec.db }, { "id", it.accession }, { "rettype", spec.db == "protein" and "gp" or "gbwithparts" },
      { "retmode", "text" } }
    if from then
      params[#params + 1] = { "seq_start", from }
      params[#params + 1] = { "seq_stop", to }
    end
    if spec.strand then params[#params + 1] = { "strand", spec.strand } end
    eutil(h, cfg, "efetch", params, {}, function(err2, resp)
      if err2 then h:_finish(err2); return end
      local text = resp.body or ""
      local e = body_error(text)
      if e or not text:match("^%s*LOCUS") then
        h:_finish(("NCBI efetch %s: %s"):format(it.accession, e or ("unexpected reply: " .. trim(text:sub(1, 120)))))
        return
      end
      if spec.db ~= "protein" and not text:find("\nORIGIN", 1, true) and text:find("\nCONTIG", 1, true) then
        h:_finish(("NCBI efetch %s: CON record returned without sequence"):format(it.accession))
        return
      end
      if path then cache_store(cfg, spec, path, text, meta) end
      local res = vim.tbl_extend("force", {}, meta)
      res.text, res.path, res.from_cache, res.id = text, path, false, spec.id
      h:_finish(nil, res)
    end)
  end)
  return h
end

return M
