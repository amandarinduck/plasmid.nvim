-- Primers in a plasmid buffer: adding, listing and laying them out.
local P = require("plasmid.bio.primer")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end
local function cfg() return require("plasmid.config").options.primers or {} end

local function next_name(rec)
  local used = {}
  for _, f in ipairs(rec.features) do
    local l = f:get("label")
    if l then used[l] = true end
  end
  local i = 1
  while used["P" .. i] do i = i + 1 end
  return "P" .. i
end

--- A primer feature for a binding site.
function M.make_feature(rec, site, seq, name)
  local author = require("plasmid.author")
  local f = author.make_feature(rec, {
    type = "primer_bind", label = name, sel = { start = site.start, stop = site.stop }, strand = site.strand,
    note = P.note(seq),
  })
  f.primer_seq = P.clean(seq)
  return f
end

local function commit_sites(buf, rec, sites, seq, name)
  local ops = {}
  for i, s in ipairs(sites) do
    local f = M.make_feature(rec, s, seq, #sites > 1 and (name .. "_" .. i) or name)
    ops[#ops + 1] = { op = "add", feature = f }
  end
  local op = #ops == 1 and ops[1] or { op = "batch", ops = ops, label = "add primer " .. name }
  require("plasmid.author").commit(buf, op, "add primer " .. name)
  -- with several sites annotated, report the one that primes best: the others
  -- are drawn too, and the strongest is not always the first in the record
  local s = #sites > 1 and M.rank_sites(sites)[1] or sites[1]
  notify(("primer %s: %d nt, anneals %d nt at %d..%d (%s)%s%s%s"):format(name, #P.clean(seq), s.anneal, s.start, s.stop,
    s.strand == 1 and "forward" or "reverse",
    #s.tail > 0 and (", 5' tail " .. #s.tail .. " nt") or "",
    #s.mismatches > 0 and (", " .. #s.mismatches .. " mismatch" .. (#s.mismatches > 1 and "es" or "")) or "",
    #sites > 1 and ("; %d sites annotated, this is the strongest"):format(#sites) or ""))
  return ops
end

local function site_text(s)
  return ("%d..%d %s  anneals %d nt%s%s  (%d exact from the 3' end)"):format(s.start, s.stop,
    s.strand == 1 and "forward" or "reverse", s.anneal,
    #s.tail > 0 and (", tail " .. #s.tail) or "", #s.mismatches > 0 and (", " .. #s.mismatches .. " mm") or "",
    s.exact3 or s.anneal)
end

--- Binding sites strongest first. What decides whether a site primes is the
--- run of bases that match exactly from the 3' end, so a site that matches 25
--- of 25 outranks one that matches 12 even where the weaker one comes first in
--- the record -- which is the order bind() returns them in, and the order a
--- picker would otherwise offer them in.
---@return table[] a sorted copy
function M.rank_sites(sites)
  local out = vim.deepcopy(sites)
  table.sort(out, function(a, b)
    local ea, eb = a.exact3 or a.anneal, b.exact3 or b.anneal
    if ea ~= eb then return ea > eb end
    if a.anneal ~= b.anneal then return a.anneal > b.anneal end
    if #a.mismatches ~= #b.mismatches then return #a.mismatches < #b.mismatches end
    if a.start ~= b.start then return a.start < b.start end
    return a.strand > b.strand
  end)
  return out
end

--- :DnaPrimerAdd {seq} [name]: annotate where a primer binds.
function M.add_sequence(buf, seq, name)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  local ok, err = P.valid(seq)
  if not ok then return notify(err, vim.log.levels.ERROR) end
  name = name or next_name(rec)
  local sites = P.bind(rec.seq:tostring(), rec:is_circular(), seq, cfg().binding)
  if #sites == 0 then
    return notify(("%s does not bind (its 3'-terminal %d bases must match)"):format(name,
      (cfg().binding or P.defaults).anchor or P.defaults.anchor), vim.log.levels.WARN)
  end
  if #sites == 1 then return commit_sites(buf, rec, sites, seq, name) end
  local items = M.rank_sites(sites)
  local best = items[1]
  table.insert(items, 1, { all = true })
  vim.ui.select(items, {
    prompt = ("%s binds at %d sites, strongest first"):format(name, #sites),
    format_item = function(s)
      if s.all then return ("all %d sites"):format(#sites) end
      return site_text(s) .. (s == best and "   <- strongest" or "")
    end,
  }, function(choice)
    if not choice then return end
    commit_sites(buf, rec, choice.all and sites or { choice }, seq, name)
  end)
end

--- From the selection: strand "fwd" (top strand) or "rev" (reverse complement).
function M.add_selection(buf, strand, name, sel)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  sel = sel or require("plasmid.select").get(buf)
  if not sel then return notify("select the primer's binding region first", vim.log.levels.WARN) end
  local top = require("plasmid.select").top_strand(rec, sel)
  local seq = strand == "rev" and P.revcomp(top) or top
  name = name or next_name(rec)
  local site = { strand = strand == "rev" and -1 or 1, start = sel.start, stop = sel.stop, anneal = #seq, tail = "",
    mismatches = {} }
  commit_sites(buf, rec, { site }, seq, name)
  require("plasmid.select").set(buf, nil)
end

--- Visual gp: pick the strand, then a name.
function M.add_visual(buf)
  local sel_ = require("plasmid.select")
  local a, b = sel_.visual_range(buf)
  vim.cmd("normal! \27")
  if not a then return end
  local rec = buffer().record(buffer().get(buf))
  local sel = sel_.make(rec, a, b)
  local top = sel_.top_strand(rec, sel)
  vim.ui.select({ "fwd", "rev" }, {
    prompt = ("Primer from %d..%d (%d nt)"):format(sel.start, sel.stop, #top),
    format_item = function(x)
      return x == "fwd" and ("forward  5'-%s-3'"):format(top) or ("reverse  5'-%s-3'"):format(P.revcomp(top))
    end,
  }, function(strand)
    if not strand then return end
    vim.ui.input({ prompt = "Primer name: ", default = next_name(rec) }, function(name)
      if not name or name == "" then return end
      M.add_selection(buf, strand, name, sel)
    end)
  end)
end

--- Primers of a record with their alignment to the template:
--- { f, name, seq, strand, start, stop, anneal, tail, mism = { [primer index] = true } }.
function M.list(rec)
  local out = {}
  local L = rec:len()
  local tmpl = rec.seq
  for _, f in ipairs(rec.features) do
    if P.is_primer(f) and f.loc and #f.parts > 0 then
      local seq = P.sequence(rec, f)
      local ext = rec:extents(f)
      local s, e
      if #ext == 1 then s, e = ext[1][1], ext[1][2]
      else
        -- across the origin: the extent that starts at 1 is the tail end
        for _, x in ipairs(ext) do if x[1] == 1 then e = x[2] else s = x[1] end end
        s, e = s or ext[1][1], e or ext[#ext][2]
      end
      local span = (e >= s) and (e - s + 1) or (L - s + 1 + e)
      local n = #seq
      local anneal = math.min(span, n)
      local strand = f.strand == -1 and -1 or 1
      local mism = {}
      -- compare the annealed part with the template
      for k = 1, anneal do
        local pos = strand == 1 and (e - (k - 1)) or (s + (k - 1))
        pos = (pos - 1) % L + 1
        local x = tmpl:sub(pos, pos):upper()
        local pb = seq:sub(n - k + 1, n - k + 1)
        local want = strand == 1 and pb or P.revcomp(pb)
        if x ~= want and want ~= "N" then mism[n - k + 1] = true end
      end
      out[#out + 1] = { f = f, name = f:get("label") or f:label(), seq = seq, strand = strand, start = s, stop = e,
        anneal = anneal, tail = seq:sub(1, n - anneal), mism = mism }
    end
  end
  return out
end

--- Display cells of a primer: list of { pos, ch, kind = "tail"|"anneal"|"mismatch" }
--- in template coordinates (wrapped on circular records, clipped on linear).
function M.cells(rec, p)
  local L = rec:len()
  local circ = rec:is_circular()
  local n = #p.seq
  local out = {}
  for i = 1, n do
    local k = n - i + 1 -- position from the 3' end
    local pos = p.strand == 1 and (p.stop - (k - 1)) or (p.start + (k - 1))
    if circ then pos = (pos - 1) % L + 1 end
    if pos >= 1 and pos <= L then
      local kind = (k > p.anneal) and "tail" or (p.mism[i] and "mismatch" or "anneal")
      local ch = p.seq:sub(i, i)
      if kind == "tail" then ch = ch:lower() end
      out[#out + 1] = { pos = pos, ch = ch, kind = kind }
    end
  end
  return out
end

return M
