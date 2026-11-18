-- PCR simulation: products of a primer pair on a linear or circular template.
local P = require("plasmid.bio.primer")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end
local function cfg() return require("plasmid.config").options.primers or {} end

--- Products of two primers (sequences 5'->3'). Every forward-binding site of
--- either primer is paired with every reverse-binding site of the other one
--- (or of itself when opts.single). Returns products sorted by length:
--- { fwd = site, rev = site, fwd_primer, rev_primer, region = template span
---   length, length = product length incl. tails, start, stop }.
function M.products(template, circular, p1, p2, opts)
  opts = opts or {}
  local L = #template
  local max_len = opts.max_len or (circular and L or math.huge)
  local s1 = P.bind(template, circular, p1, opts.binding)
  local s2 = P.bind(template, circular, p2, opts.binding)
  local fwd, rev = {}, {}
  for _, s in ipairs(s1) do if s.strand == 1 then fwd[#fwd + 1] = { s, 1 } else rev[#rev + 1] = { s, 1 } end end
  for _, s in ipairs(s2) do if s.strand == 1 then fwd[#fwd + 1] = { s, 2 } else rev[#rev + 1] = { s, 2 } end end
  local seqs = { P.clean(p1), P.clean(p2) }
  local out = {}
  for _, f in ipairs(fwd) do
    for _, r in ipairs(rev) do
      if f[2] ~= r[2] or opts.single then
        local a, b = f[1], r[1]
        local region
        if circular then
          region = (b.stop - a.start) % L + 1
        elseif b.stop >= a.start and b.start >= a.start then
          region = b.stop - a.start + 1
        end
        -- the reverse primer's 3' end must not lie upstream of the forward primer's
        local ok = region and region >= math.max(a.anneal, b.anneal)
        if ok and circular then
          local off3 = (b.start - a.start) % L
          ok = off3 < region
        end
        if ok and region <= max_len then
          out[#out + 1] = { fwd = a, rev = b, fwd_primer = seqs[f[2]], rev_primer = seqs[r[2]], fwd_index = f[2],
            rev_index = r[2], region = region, length = region + #a.tail + #b.tail, start = a.start, stop = b.stop }
        end
      end
    end
  end
  table.sort(out, function(x, y)
    if x.length ~= y.length then return x.length < y.length end
    return x.start < y.start
  end)
  return out, s1, s2
end

--- Product sequence (top strand 5'->3').
function M.product_sequence(template, prod)
  local L = #template
  local function sub(a, b)
    if a <= b then return template:sub(a, b) end
    return template:sub(a, L) .. template:sub(1, b)
  end
  local mid = sub(prod.start, prod.stop):upper()
  local fa, ra = prod.fwd.anneal, prod.rev.anneal
  local fp, rp = prod.fwd_primer, P.revcomp(prod.rev_primer)
  if fa + ra <= #mid then
    -- primer bases replace the template under each primer (mismatches end up in the product)
    mid = fp:sub(#fp - fa + 1) .. mid:sub(fa + 1, #mid - ra) .. rp:sub(1, ra)
  end
  return prod.fwd.tail .. mid .. P.revcomp(prod.rev.tail)
end

--- The product as a new linear record with the template's features inside it.
function M.build(rec, prod, names)
  local T = require("plasmid.model.transaction")
  local copy = require("plasmid.digest").extract(rec, { start = prod.start, stop = prod.stop, length = prod.region })
  local want = M.product_sequence(rec.seq:tostring(), prod)
  -- follow the template's letter case
  local lower = rec.seq:sub(prod.start, prod.start):match("%l") ~= nil
  local function cased(x) return lower and x:lower() or x end
  -- substitute primer mismatches (same length: features unaffected)
  local core = want:sub(#prod.fwd.tail + 1, #want - #prod.rev.tail)
  local cur = copy.seq:tostring():upper()
  for i = 1, #core do
    if cur:sub(i, i) ~= core:sub(i, i) then
      T.apply(copy, { op = "splice", start = i, stop = i, ins = cased(core:sub(i, i)) })
    end
  end
  -- drop carried copies of the two primers; they are re-added full length
  for i = #copy.features, 1, -1 do
    local f = copy.features[i]
    local l = f:get("label")
    if P.is_primer(f) and (l == names[1] or l == names[2]) then T.apply(copy, { op = "remove", id = f.id }) end
  end
  if #prod.fwd.tail > 0 then T.apply(copy, { op = "splice", start = 1, stop = 0, ins = cased(prod.fwd.tail) }) end
  if #prod.rev.tail > 0 then
    local n = copy:len()
    T.apply(copy, { op = "splice", start = n + 1, stop = n, ins = cased(P.revcomp(prod.rev.tail)) })
  end
  local n = copy:len()
  local Feature = require("plasmid.model.record").Feature
  for _, f in ipairs(copy.features) do
    if f.type == "source" then
      local g = Feature.new({ type = "source", location_text = "1.." .. n, qualifiers = f.qualifiers })
      T.apply(copy, { op = "replace", id = f.id, feature = g })
    end
  end
  local author = require("plasmid.author")
  local function add(seq, name, s, e, strand)
    local f = author.make_feature(copy, { type = "primer_bind", label = name, sel = { start = s, stop = e },
      strand = strand, note = P.note(seq) })
    f.primer_seq = seq
    T.apply(copy, { op = "add", feature = f })
  end
  add(prod.fwd_primer, names[1], 1, #prod.fwd_primer, 1)
  add(prod.rev_primer, names[2], n - #prod.rev_primer + 1, n, -1)
  copy.history = nil
  copy.name = ("%s_PCR_%s_%s"):format(rec.name, names[1], names[2]):gsub("[^%w%._%-]", "_"):sub(1, 40)
  copy.description = ("PCR product of %s with %s and %s (%d bp; template %d..%d)"):format(rec.name, names[1],
    names[2], n, prod.start, prod.stop)
  assert(copy.seq:tostring():upper() == want, "product sequence mismatch")
  return copy
end

-- A primer argument: the label of a primer feature, or a sequence.
local function resolve(rec, arg)
  for _, p in ipairs(require("plasmid.primers").list(rec)) do
    if p.name == arg then return p.seq, p.name end
  end
  if P.valid(arg) then return P.clean(arg), arg:sub(1, 12) end
  return nil
end

local function ta_text(f, r)
  local ok, tm = pcall(require, "plasmid.bio.tm")
  if not ok or not tm.ta then return "" end
  local o = require("plasmid.config").options.primers or {}
  local ta, t1, t2 = tm.ta(f, r, { polymerase = o.polymerase, dmso = o.dmso, conc_nM = o.conc_nM, salt_mM = o.salt_mM })
  if not ta then return "" end
  return (", Tm %d/%d °C, Ta %d °C (%s)"):format(t1, t2, ta, o.polymerase or "q5")
end

--- :DnaPCR [primer1] [primer2]
function M.run(buf, a1, a2)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  local function go(s1, n1, s2, n2)
    local prods = M.products(rec.seq:tostring(), rec:is_circular(), s1, s2, { binding = cfg().binding })
    if #prods == 0 then
      return notify(("no product: %s and %s do not bind facing each other"):format(n1, n2), vim.log.levels.WARN)
    end
    local function make(prod)
      local names = prod.fwd_index == 1 and { n1, n2 } or { n2, n1 }
      local newrec = M.build(rec, prod, names)
      -- Tm/Ta on the annealed parts
      local fa = prod.fwd_primer:sub(#prod.fwd.tail + 1)
      local ra = prod.rev_primer:sub(#prod.rev.tail + 1)
      require("plasmid.digest").open_record(buf, newrec, ("PCR product %d bp (template %d..%d)%s%s"):format(
        newrec:len(), prod.start, prod.stop, ta_text(fa, ra),
        #prods > 1 and (", 1 of " .. #prods .. " products") or ""))
      return newrec
    end
    if #prods == 1 then return make(prods[1]) end
    vim.ui.select(prods, {
      prompt = ("%d possible products"):format(#prods),
      format_item = function(p) return ("%6d bp  template %d..%d"):format(p.length, p.start, p.stop) end,
    }, function(choice) if choice then make(choice) end end)
  end
  if a1 and a2 then
    local s1, n1 = resolve(rec, a1)
    local s2, n2 = resolve(rec, a2)
    if not s1 then return notify("not a primer name or sequence: " .. a1, vim.log.levels.ERROR) end
    if not s2 then return notify("not a primer name or sequence: " .. a2, vim.log.levels.ERROR) end
    return go(s1, n1, s2, n2)
  end
  local list = require("plasmid.primers").list(rec)
  if #list < 2 then return notify("give two primers (:DnaPCR fwd rev) or add primers first", vim.log.levels.WARN) end
  if #list == 2 then return go(list[1].seq, list[1].name, list[2].seq, list[2].name) end
  local fmt = function(p) return ("%-14s %s %d..%d"):format(p.name, p.strand == 1 and "→" or "←", p.start, p.stop) end
  vim.ui.select(list, { prompt = "First primer", format_item = fmt }, function(p1)
    if not p1 then return end
    vim.ui.select(list, { prompt = "Second primer", format_item = fmt }, function(p2)
      if p2 then go(p1.seq, p1.name, p2.seq, p2.name) end
    end)
  end)
end

return M
