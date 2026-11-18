-- Restriction sites shown in a plasmid buffer: :DnaEnzymes, the picker, info.
local R = require("plasmid.bio.restriction")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function cfg() return require("plasmid.config").options.enzymes or {} end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end

--- Specific letters in a site (N does not count): BglI GCCNNNNNGGC = 6.
function M.specific_len(site)
  local n = 0
  for c in site:gmatch("[^N]") do n = n + 1 end
  return n
end

local MODES = { unique = true, dual = true, all = true, names = true }

--- Sites to display for a spec { mode = "unique"|"dual"|"all"|"names", names = {...} }.
--- Returns entries sorted by cut and a summary. An entry groups enzymes that
--- cut at the same place (isoschizomers): { cut, cut_bottom, start, stop,
--- strand, label, names, enzymes, count (cuts per enzyme), overhang, methyl }.
function M.compute(rec, spec)
  local o = cfg()
  local min_site = tonumber(o.min_site) or 6
  local list, unknown = {}, {}
  if spec.mode == "names" then
    for _, n in ipairs(spec.names or {}) do
      local e = R.get(n)
      if e then list[#list + 1] = e else unknown[#unknown + 1] = n end
    end
  else
    for _, e in ipairs(R.enzymes()) do
      -- imported REBASE sets: only commercially available enzymes, unless named
      local available = e.suppliers == nil or e.suppliers ~= "" or o.commercial_only == false
      if available and not R.NEEDS_METHYLATION[e.name] and M.specific_len(e.site) >= min_site then
        list[#list + 1] = e
      end
    end
  end
  local P = R.prepare(rec)
  local groups, entries, zero, cutting = {}, {}, {}, 0
  for _, e in ipairs(list) do
    local sites = R.sites(P, e)
    local n = #sites
    local show = (spec.mode == "names" and n > 0) or (spec.mode == "unique" and n == 1)
      or (spec.mode == "dual" and (n == 1 or n == 2)) or (spec.mode == "all" and n > 0)
    if n == 0 and spec.mode == "names" then zero[#zero + 1] = e.name end
    if show then
      cutting = cutting + 1
      for _, s in ipairs(sites) do
        local key = s.cut .. ":" .. s.cut_bottom
        local g = groups[key]
        if not g then
          g = { cut = s.cut, cut_bottom = s.cut_bottom, start = s.start, stop = s.stop, strand = s.strand,
            names = {}, enzymes = {}, count = n, overhang = s.overhang, methyl = {} }
          groups[key] = g
          entries[#entries + 1] = g
        end
        g.names[#g.names + 1] = e.name
        g.enzymes[#g.enzymes + 1] = e
        g.count = math.min(g.count, n)
        for _, m in ipairs(R.methylation(P, s)) do g.methyl[m] = true end
      end
    end
  end
  for _, g in ipairs(entries) do g.label = table.concat(g.names, ", ") end
  table.sort(entries, function(a, b)
    if a.cut ~= b.cut then return a.cut < b.cut end
    return a.label < b.label
  end)
  return entries, { enzymes = cutting, sites = #entries, zero = zero, unknown = unknown, considered = #list }
end

local function summary_text(spec, sum)
  local what = ({ unique = "unique cutters", dual = "1–2 cutters", all = "cutting enzymes" })[spec.mode]
    or "of the named enzymes cut"
  local s = ("%d %s, %d site%s"):format(sum.enzymes, what, sum.sites, sum.sites == 1 and "" or "s")
  if #sum.zero > 0 then s = s .. "; no site: " .. table.concat(sum.zero, ", ") end
  if #sum.unknown > 0 then s = s .. "; unknown: " .. table.concat(sum.unknown, ", ") end
  return s
end

--- Show a set: "unique" | "dual" | "all" | list of names | "off".
function M.set(buf, what, quiet)
  local st = buffer().get(buf)
  if not st then return end
  if buffer().record(st):is_protein() then return notify("restriction sites need a nucleotide sequence", vim.log.levels.WARN) end
  if what == "off" or what == nil then
    st.enzymes = nil
    buffer().render(buf)
    return
  end
  local spec
  if type(what) == "table" then spec = { mode = "names", names = what }
  elseif MODES[what] then spec = { mode = what }
  else return notify("unknown enzyme set " .. tostring(what), vim.log.levels.ERROR) end
  local entries, sum = M.compute(buffer().record(st), spec)
  st.enzymes = { spec = spec, entries = entries, summary = sum }
  buffer().render(buf)
  if not quiet then
    notify(summary_text(spec, sum), (#sum.unknown > 0) and vim.log.levels.WARN or vim.log.levels.INFO)
  end
  return st.enzymes
end

function M.toggle(buf)
  local st = buffer().get(buf)
  if not st then return end
  if st.enzymes then return M.set(buf, "off") end
  return M.set(buf, cfg().default_set or "unique")
end

--- After an edit: recompute the shown set.
function M.refresh(buf)
  local st = buffer().get(buf)
  local E = st and st.enzymes
  if not E then return end
  local entries, sum = M.compute(buffer().record(st), E.spec)
  E.entries, E.summary = entries, sum
end

function M.describe(st)
  local E = st and st.enzymes
  if not E then return nil end
  local tag = E.spec.mode == "names" and "named" or E.spec.mode
  return ("RE %s: %d"):format(tag, E.summary.enzymes)
end

--- Entry under a display position (linear label rows, circular labels).
function M.entry_at(st, row0, vcol)
  local lay = st and st.layout
  if not lay then return nil end
  if lay.view == "linear" then
    local info = lay.rows[row0]
    for _, s in ipairs(info and info.esegs or {}) do
      if vcol >= s[1] and vcol <= s[2] then return s[3] end
    end
  else
    for _, l in ipairs(lay.enz_labels or {}) do
      if row0 == l.row + (lay.header_rows or 0) and vcol >= l.col and vcol < l.col + l.len then return l.entry end
    end
  end
end

local function cut_text(L, p)
  -- cut between p-1 and p
  local before = p - 1
  if before < 1 then before = L end
  return ("%d^%d"):format(before, p)
end

function M.info_lines(rec, entry)
  local L = rec:len()
  local lines = { entry.label }
  for _, e in ipairs(entry.enzymes) do
    lines[#lines + 1] = ("  %-10s %-14s %s%s"):format(e.name, e.site, R.end_kind(e),
      e.suppliers and ("  suppliers " .. e.suppliers) or "")
  end
  local where = entry.stop >= entry.start and ("%d..%d"):format(entry.start, entry.stop)
    or ("%d..%d (across the origin)"):format(entry.start, entry.stop)
  lines[#lines + 1] = ("site     %s on the %s strand"):format(where, entry.strand == 1 and "top" or "bottom")
  lines[#lines + 1] = ("cuts     top %s   bottom %s"):format(cut_text(L, entry.cut), cut_text(L, entry.cut_bottom))
  lines[#lines + 1] = ("sites    %d in this molecule"):format(entry.count)
  local m = {}
  for k in pairs(entry.methyl) do m[#m + 1] = k == "dam" and "Dam (GATC)" or "Dcm (CCWGG)" end
  table.sort(m)
  if #m > 0 then
    lines[#lines + 1] = "warning  overlaps " .. table.concat(m, " and ") .. " methylation: may not cut DNA from dam+/dcm+ E. coli"
  end
  return lines
end

function M.show_info(buf, entry)
  local st = buffer().get(buf)
  local lines = M.info_lines(buffer().record(st), entry)
  local fbuf = vim.lsp.util.open_floating_preview(lines, "", {
    border = "rounded", max_width = math.min(90, vim.o.columns - 4), focus_id = "plasmid_enzyme",
  })
  vim.api.nvim_buf_set_extmark(fbuf, buffer().ns, 0, 0, { end_col = #lines[1], hl_group = "PlasmidEnzymeUnique" })
end

--- Picker: every enzyme with its number of sites; choosing one toggles it in
--- the named set.
function M.pick(buf)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  local P = R.prepare(rec)
  local items = {}
  for _, e in ipairs(R.enzymes()) do
    items[#items + 1] = { e = e, n = #R.sites(P, e) }
  end
  local shown = {}
  if st.enzymes and st.enzymes.spec.mode == "names" then
    for _, n in ipairs(st.enzymes.spec.names) do shown[n:lower()] = true end
  end
  vim.ui.select(items, {
    prompt = "Restriction enzyme (toggles it in the named set)",
    format_item = function(it)
      local n = it.n == 0 and "none" or (it.n == 1 and "1 site" or (it.n .. " sites"))
      return ("%s %-10s %-14s %-8s %s"):format(shown[it.e.name:lower()] and "●" or " ", it.e.name, it.e.site, n, R.end_kind(it.e))
    end,
  }, function(it)
    if not it then return end
    local names = {}
    for _, n in ipairs((st.enzymes and st.enzymes.spec.mode == "names") and st.enzymes.spec.names or {}) do
      if n:lower() ~= it.e.name:lower() then names[#names + 1] = n end
    end
    if not shown[it.e.name:lower()] then names[#names + 1] = it.e.name end
    if #names == 0 then M.set(buf, "off") else M.set(buf, names) end
  end)
end

--- Entries of the shown set as a list of cuts (for digests).
function M.cuts(st)
  local out = {}
  for _, g in ipairs(st and st.enzymes and st.enzymes.entries or {}) do out[#out + 1] = g end
  return out
end

return M
