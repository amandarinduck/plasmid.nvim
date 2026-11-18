-- Motif search in a plasmid buffer: :DnaFind, / and ?, n/N, gn.
local search = require("plasmid.bio.search")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end

local function index(hits)
  local wrapped, maxlen = {}, 1
  for i, h in ipairs(hits) do
    if h.stop < h.start then wrapped[#wrapped + 1] = i
    else maxlen = math.max(maxlen, h.stop - h.start + 1) end
  end
  return wrapped, maxlen
end

local function compute(rec, text, k)
  local hits, kind, motif = search.find(rec, text, k)
  if not hits then return nil, kind end
  local wrapped, maxlen = index(hits)
  return { text = text, motif = motif, kind = kind, k = k or 0, hits = hits, wrapped = wrapped, maxlen = maxlen, idx = nil }
end

local function hit_text(rec, S, i)
  local h = S.hits[i]
  local strand = h.strand == -1 and "−" or "+"
  local mm = (h.mismatches or 0) > 0 and (", %d mismatch%s"):format(h.mismatches, h.mismatches == 1 and "" or "es") or ""
  return ("%s %d/%d  %d..%d (%s%s)"):format(S.motif, i, #S.hits, h.start, h.stop, strand, mm)
end

-- Index of the next hit from bp in direction dir (wrapping around).
local function nearest(S, bp, dir, strictly)
  local n = #S.hits
  if dir > 0 then
    for i = 1, n do
      local s = S.hits[i].start
      if s > bp or (not strictly and s == bp) then return i end
    end
    return 1
  else
    for i = n, 1, -1 do
      local s = S.hits[i].start
      if s < bp or (not strictly and s == bp) then return i end
    end
    return n
  end
end

local function go(buf, S, i)
  S.idx = i
  local win = vim.fn.bufwinid(buf)
  buffer().render(buf, { keep_bp = S.hits[i].start })
  if win ~= -1 then buffer().goto_bp(win, S.hits[i].start) end
end

--- Search `text` (motif; "aa:" prefix for amino acids) with k mismatches.
--- dir = 1 (/) or -1 (?) sets the direction of the first jump and of n.
function M.run(buf, text, k, dir)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  local S, err = compute(rec, text, k)
  if not S then return notify(err, vim.log.levels.ERROR) end
  S.dir = dir or 1
  st.search = S
  if #S.hits == 0 then
    buffer().render(buf)
    return notify(("%s: no match%s"):format(S.motif, S.k > 0 and (" with ≤%d mismatches"):format(S.k) or ""),
      vim.log.levels.WARN)
  end
  local win = vim.fn.bufwinid(buf)
  local bp = win ~= -1 and (buffer().cursor_info(win)) or 1
  local i = nearest(S, bp or 1, S.dir, false)
  go(buf, S, i)
  notify(("%s: %d hit%s%s"):format(S.motif, #S.hits, #S.hits == 1 and "" or "s",
    S.k > 0 and (" (≤%d mismatches)"):format(S.k) or ""))
  return S
end

--- n / N: next / previous hit, relative to the search direction.
function M.next(buf, forward, count)
  local st = buffer().get(buf)
  local S = st and st.search
  if not S or #S.hits == 0 then return notify("no search; use / or :DnaFind", vim.log.levels.WARN) end
  local dir = (forward and 1 or -1) * (S.dir or 1)
  local win = vim.fn.bufwinid(buf)
  local i = S.idx
  for _ = 1, math.max(1, count or 1) do
    local bp = S.hits[i or 1].start
    if not i then bp = (win ~= -1 and buffer().cursor_info(win)) or 1 end
    local j = nearest(S, bp, dir, true)
    if i and dir > 0 and j <= i and S.hits[j].start == S.hits[i].start then
      -- several hits share a start (both strands): step through them in order
      j = i < #S.hits and i + 1 or 1
    end
    i = j
  end
  go(buf, S, i)
  notify(hit_text(buffer().record(st), S, i))
end

--- gn: select the current hit (or the next one).
function M.select_hit(buf)
  local st = buffer().get(buf)
  local S = st and st.search
  if not S or #S.hits == 0 then return notify("no search; use / or :DnaFind", vim.log.levels.WARN) end
  local i = S.idx or nearest(S, (buffer().cursor_info(vim.fn.bufwinid(buf))) or 1, 1, false)
  local h = S.hits[i]
  local sel = require("plasmid.select")
  sel.set(buf, sel.make(buffer().record(st), h.start, h.stop))
  S.idx = i
end

function M.clear(buf)
  local st = buffer().get(buf)
  if not st then return end
  st.search = nil
  buffer().render(buf)
end

--- After the residues changed: search again, keeping the motif.
function M.refresh(buf)
  local st = buffer().get(buf)
  local S = st and st.search
  if not S then return end
  local fresh = compute(buffer().record(st), S.text, S.k)
  if fresh then fresh.dir = S.dir; st.search = fresh end
end

--- Winbar fragment.
function M.describe(st)
  local S = st and st.search
  if not S then return nil end
  if #S.hits == 0 then return ("/%s: none"):format(S.motif) end
  if not S.idx then return ("/%s: %d"):format(S.motif, #S.hits) end
  return ("/%s %d/%d"):format(S.motif, S.idx, #S.hits)
end

--- The / and ? prompt. Input "MOTIF", "MOTIF 1" (mismatches) or "aa:MOTIF".
function M.prompt(buf, dir)
  vim.ui.input({ prompt = dir > 0 and "/" or "?" }, function(text)
    if not text or text:match("^%s*$") then return end
    local motif, k = text:match("^%s*(%S+)%s+(%d+)%s*$")
    M.run(buf, motif or text, tonumber(k) or 0, dir)
  end)
end

return M
