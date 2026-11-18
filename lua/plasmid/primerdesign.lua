-- Primer design for a selected region.
--
-- The selection is the region to amplify. A forward primer starts at its
-- first base, a reverse primer ends at its last one; both may also start up
-- to `slack` bases further out (the product grows, the selection stays
-- inside), which is what makes a Tm match possible at a GC-rich boundary. Each primer is scored on its distance from the
-- target Tm, GC content, a 3'-terminal G/C, homopolymer runs, how uniquely it
-- binds the template and (when bio/oligo is available) hairpin, self-dimer
-- and 3'-end stability; pairs are then scored on the Tm difference between
-- them. 5' tails (restriction sites, overlaps) are added after the search and
-- do not take part in it.
local P = require("plasmid.bio.primer")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end

M.defaults = {
  len = { 18, 30 },          -- primer length range (annealing part)
  tm = { 57, 68 },           -- acceptable Tm window (°C)
  target_tm = nil,           -- defaults to the middle of the window
  max_dtm = 3,               -- largest Tm difference within a pair (°C)
  gc = { 0.35, 0.65 },       -- acceptable GC fraction
  max_run = 4,               -- longest homopolymer run allowed without penalty
  slack = 8,                 -- bases a primer may start outside the selection (the product grows)
  pairs = 8,                 -- candidates offered
}

local function opts()
  local c = require("plasmid.config").options.primers or {}
  return vim.tbl_deep_extend("force", M.defaults, c.design or {})
end

local function tm_mod() return require("plasmid.bio.tm") end

-- Score one primer. Lower is better; nil means it is unusable.
local function score_primer(seq, template, circular, o, cond, binding, seen)
  local d = tm_mod().tm_detail(seq, cond)
  if not d then return nil end
  local tmv = d.tm
  if seen then seen.lo = math.min(seen.lo, tmv); seen.hi = math.max(seen.hi, tmv) end
  if tmv < o.tm[1] - 3 or tmv > o.tm[2] + 3 then return nil end
  local target = o.target_tm or (o.tm[1] + o.tm[2]) / 2
  local s = math.abs(tmv - target)
  local why = {}
  if tmv < o.tm[1] or tmv > o.tm[2] then s = s + 4; why[#why + 1] = "Tm outside the window" end
  local gc = 0
  for _ in seq:gmatch("[GCgc]") do gc = gc + 1 end
  gc = gc / #seq
  if gc < o.gc[1] or gc > o.gc[2] then s = s + 6; why[#why + 1] = ("GC %d %%"):format(math.floor(gc * 100 + 0.5)) end
  local last = seq:sub(-1)
  if last ~= "G" and last ~= "C" then s = s + 1.5; why[#why + 1] = "no 3' G/C" end
  if seq:sub(-5):find("[GC][GC][GC]") then s = s + 1; why[#why + 1] = "GC-rich 3' end" end
  local run, best_run = 1, 1
  for i = 2, #seq do
    run = (seq:sub(i, i) == seq:sub(i - 1, i - 1)) and run + 1 or 1
    best_run = math.max(best_run, run)
  end
  if best_run > o.max_run then s = s + 2 * (best_run - o.max_run); why[#why + 1] = best_run .. "-base run" end
  local sites = P.bind(template, circular, seq, binding)
  if #sites == 0 then return nil end
  if #sites > 1 then s = s + 8; why[#why + 1] = ("binds %d sites"):format(#sites) end
  local hp, sd, es
  local ok, oligo = pcall(require, "plasmid.bio.oligo")
  if ok then
    local h = oligo.hairpin(seq)
    local dm = oligo.homodimer(seq)
    es = oligo.end_stability(seq)
    hp = h and h.dg
    sd = dm and dm.dg
    if hp and hp < -2 then s = s + (-hp - 2) * 1.5; why[#why + 1] = ("hairpin %.1f"):format(hp) end
    if sd and sd < -6 then s = s + (-sd - 6); why[#why + 1] = ("self-dimer %.1f"):format(sd) end
  end
  return { seq = seq, tm = tmv, tm_rounded = d.tm_rounded, gc = gc, score = s, sites = #sites, why = why,
    hairpin = hp, self_dimer = sd, end_stability = es, len = #seq }
end

--- Rank primer pairs for a region (start..stop, inclusive, may wrap).
--- Returns a list of { fwd, rev, ta, dtm, score, product } sorted by score.
function M.candidates(rec, sel, over)
  local o = vim.tbl_deep_extend("force", opts(), over or {})
  local cfg = require("plasmid.config").options.primers or {}
  local cond = require("plasmid.primerinfo").conditions(o.cond)
  local template = rec.seq:tostring()
  local L = #template
  local circular = rec:is_circular()
  local function sub(a, n)
    local out = {}
    for i = 0, n - 1 do
      local p = a + i
      if circular then p = (p - 1) % L + 1 elseif p < 1 or p > L then return nil end
      out[#out + 1] = template:sub(p, p)
    end
    return table.concat(out):upper()
  end
  local span = (sel.stop >= sel.start) and (sel.stop - sel.start + 1) or (L - sel.start + 1 + sel.stop)
  if span < o.len[1] * 2 then return nil, ("the region is only %d bp; select at least %d"):format(span, o.len[1] * 2) end
  local function wrap(p) return (p - 1) % L + 1 end
  local fwds, revs = {}, {}
  local fseen, rseen = { lo = math.huge, hi = -math.huge }, { lo = math.huge, hi = -math.huge }
  for shift = 0, math.max(0, o.slack) do
    for n = o.len[1], o.len[2] do
      local fstart = circular and wrap(sel.start - shift) or (sel.start - shift)
      if fstart >= 1 then
        local f = sub(fstart, n)
        if f then
          local s = score_primer(f, template, circular, o, cond, cfg.binding, fseen)
          if s then
            s.score = s.score + shift * 0.1
            s.start, s.stop, s.shift, s.strand = fstart, wrap(fstart + n - 1), shift, 1
            fwds[#fwds + 1] = s
          end
        end
      end
      local rstop = circular and wrap(sel.stop + shift) or (sel.stop + shift)
      if rstop <= L or circular then
        local r = sub(circular and wrap(rstop - n + 1) or (rstop - n + 1), n)
        if r then
          local s = score_primer(P.revcomp(r), template, circular, o, cond, cfg.binding, rseen)
          if s then
            s.score = s.score + shift * 0.1
            s.start, s.stop, s.shift, s.strand = circular and wrap(rstop - n + 1) or (rstop - n + 1), rstop, shift, -1
            revs[#revs + 1] = s
          end
        end
      end
    end
  end
  if #fwds == 0 or #revs == 0 then
    local which, seen = "forward", fseen
    if #fwds > 0 then which, seen = "reverse", rseen end
    local reach = seen.hi > -math.huge and (", that side reaches %.1f-%.1f °C"):format(seen.lo, seen.hi) or ""
    return nil, ("no usable %s primer in %d-%d nt at Tm %g-%g °C%s; widen primers.design.tm or len, or raise slack (%d)")
      :format(which, o.len[1], o.len[2], o.tm[1], o.tm[2], reach, o.slack)
  end
  local out = {}
  for _, f in ipairs(fwds) do
    for _, r in ipairs(revs) do
      local dtm = math.abs(f.tm - r.tm)
      if dtm <= o.max_dtm then
        local ta, t1, t2, info = tm_mod().ta(f.seq, r.seq, cond)
        out[#out + 1] = { fwd = f, rev = r, dtm = dtm, ta = ta, tm1 = t1, tm2 = t2, ta_info = info,
          score = f.score + r.score + dtm * 2, product = span + f.shift + r.shift }
      end
    end
  end
  if #out == 0 then
    local best = math.huge
    for _, f in ipairs(fwds) do
      for _, r in ipairs(revs) do best = math.min(best, math.abs(f.tm - r.tm)) end
    end
    return nil, ("no pair within %g °C of each other (closest %.1f °C); raise primers.design.max_dtm or slack (%d)")
      :format(o.max_dtm, best, o.slack)
  end
  table.sort(out, function(a, b)
    if a.score ~= b.score then return a.score < b.score end
    return a.fwd.len + a.rev.len < b.fwd.len + b.rev.len
  end)
  local top = {}
  for i = 1, math.min(o.pairs, #out) do top[i] = out[i] end
  return top
end

local function tail_for(arg)
  if not arg or arg == "" then return "" end
  local e = require("plasmid.bio.restriction").get(arg)
  if e then return e.site:lower(), e.name end
  if P.valid(arg) then return P.clean(arg):lower() end
  return nil
end

--- Add a chosen pair as two primer features.
function M.accept(buf, sel, pair, names, tails)
  local rec = buffer().record(buffer().get(buf))
  local primers = require("plasmid.primers")
  local fseq = (tails and tails[1] or "") .. pair.fwd.seq
  local rseq = (tails and tails[2] or "") .. pair.rev.seq
  local fsite = { strand = 1, start = pair.fwd.start, stop = pair.fwd.stop,
    anneal = pair.fwd.len, tail = tails and tails[1] or "", mismatches = {} }
  local rsite = { strand = -1, start = pair.rev.start, stop = pair.rev.stop,
    anneal = pair.rev.len, tail = tails and tails[2] or "", mismatches = {} }
  local ops = {
    { op = "add", feature = primers.make_feature(rec, fsite, fseq, names[1]) },
    { op = "add", feature = primers.make_feature(rec, rsite, rseq, names[2]) },
  }
  require("plasmid.author").commit(buf, { op = "batch", ops = ops, label = "add primer pair" }, "add primer pair")
  require("plasmid.select").set(buf, nil)
  local total = pair.product + #(tails and tails[1] or "") + #(tails and tails[2] or "")
  notify(("%s / %s: Tm %d / %d °C, Ta %d °C, product %d bp%s"):format(names[1], names[2], pair.tm1, pair.tm2,
    pair.ta, total, (tails and (#tails[1] + #tails[2]) > 0) and " (with tails)" or ""))
end

local function describe(p)
  local flags = #p.fwd.why + #p.rev.why
  local grew = p.fwd.shift + p.rev.shift
  return ("%2d+%2d nt  Tm %4.1f/%4.1f  Ta %2d  ΔTm %3.1f  GC %2d/%2d %%  %d bp%s%s"):format(p.fwd.len, p.rev.len, p.fwd.tm,
    p.rev.tm, p.ta, p.dtm, math.floor(p.fwd.gc * 100 + 0.5), math.floor(p.rev.gc * 100 + 0.5), p.product,
    grew > 0 and ("  (+%d outside)"):format(grew) or "",
    flags > 0 and ("  [" .. table.concat(vim.list_extend(vim.deepcopy(p.fwd.why), p.rev.why), "; ") .. "]") or "")
end

--- :DnaPrimerDesign [name=X] [tm=62] [len=18-30] [fwd=EcoRI|SEQ] [rev=BamHI|SEQ]
function M.run(buf, args)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  if rec:is_protein() then return notify("primer design needs a nucleotide sequence", vim.log.levels.WARN) end
  local sel = require("plasmid.select").get(buf)
  if not sel then return notify("select the region to amplify first", vim.log.levels.WARN) end
  local over, name, tails = {}, nil, { "", "" }
  local tnames = {}
  for _, a in ipairs(args or {}) do
    local k, v = a:match("^(%w+)=(.+)$")
    if k == "tm" then over.target_tm = tonumber(v)
    elseif k == "len" then
      local lo, hi = v:match("^(%d+)%-(%d+)$")
      if not lo then return notify("len wants a range, e.g. len=18-30", vim.log.levels.ERROR) end
      over.len = { tonumber(lo), tonumber(hi) }
    elseif k == "slack" then over.slack = tonumber(v)
    elseif k == "dtm" then over.max_dtm = tonumber(v)
    elseif k == "name" then name = v
    elseif k == "fwd" or k == "rev" then
      local t, enz = tail_for(v)
      if not t then return notify("not an enzyme or sequence: " .. v, vim.log.levels.ERROR) end
      tails[k == "fwd" and 1 or 2] = t
      tnames[k == "fwd" and 1 or 2] = enz or v
    elseif k then return notify("unknown option " .. k, vim.log.levels.ERROR) end
  end
  local list, err = M.candidates(rec, sel, over)
  if not list then return notify(err, vim.log.levels.WARN) end
  name = name or ("P" .. sel.start)
  local names = { name .. "_F", name .. "_R" }
  vim.ui.select(list, {
    prompt = ("Primer pair for %d..%d (%d bp selected)%s"):format(sel.start, sel.stop,
      (sel.stop >= sel.start) and (sel.stop - sel.start + 1) or (rec:len() - sel.start + 1 + sel.stop),
      (#tails[1] + #tails[2]) > 0 and (" + tails " .. (tnames[1] or tails[1]) .. "/" .. (tnames[2] or tails[2])) or ""),
    format_item = describe,
  }, function(choice)
    if not choice then return end
    M.accept(buf, sel, choice, names, tails)
  end)
  return list
end

return M
