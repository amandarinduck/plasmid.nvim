-- Tm calculator and primer details: :DnaTm, :DnaPolymerase, K on a primer.
local P = require("plasmid.bio.primer")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end
local function cfg() return require("plasmid.config").options.primers or {} end
local function tm() return require("plasmid.bio.tm") end
local function oligo() return require("plasmid.bio.oligo") end

--- Reaction settings from the config (and :DnaPolymerase).
function M.conditions(over)
  local o = cfg()
  local c = { polymerase = o.polymerase or "q5", conc_nM = o.conc_nM, salt_mM = o.salt_mM, dmso = o.dmso,
    buffer = o.buffer }
  for k, v in pairs(over or {}) do c[k] = v end
  return c
end

function M.polymerase_name(key)
  local p = tm().polymerase(key or M.conditions().polymerase)
  return p and p.name or tostring(key)
end

local function fmt_dg(x) return x and ("%.1f kcal/mol"):format(x) or "—" end

--- Tm line for a sequence (annealed part, A/C/G/T only).
function M.tm_text(seq, cond)
  cond = cond or M.conditions()
  local d, err = tm().tm_detail(seq, cond)
  if not d then return nil, err end
  local extra = {}
  extra[#extra + 1] = ("%d nM"):format(d.conc_nM)
  extra[#extra + 1] = ("%s, %d mM salt"):format(d.method, d.salt_mM)
  if cond.dmso and cond.dmso > 0 then extra[#extra + 1] = ("%g%% DMSO"):format(cond.dmso) end
  return ("%d °C  (%s; %s)"):format(d.tm_rounded, M.polymerase_name(d.polymerase), table.concat(extra, ", ")), d
end

--- A primer feature's details as lines for a float.
function M.primer_lines(rec, p)
  local cond = M.conditions()
  local lines = {}
  local strand = p.strand == 1 and "forward" or "reverse"
  lines[1] = ("%s  %s primer, %d nt"):format(p.name, strand, #p.seq)
  lines[#lines + 1] = ("5'-%s-3'"):format(#p.tail > 0 and (p.tail:lower() .. p.seq:sub(#p.tail + 1)) or p.seq)
  local mm = vim.tbl_count(p.mism)
  lines[#lines + 1] = ("anneals   %d..%d (%d nt)%s%s"):format(p.start, p.stop, p.anneal,
    #p.tail > 0 and (", 5' tail " .. #p.tail .. " nt") or "", mm > 0 and (", " .. mm .. " mismatch" .. (mm > 1 and "es" or "")) or "")
  local ann = p.seq:sub(#p.tail + 1)
  local ok_tm, txt = pcall(M.tm_text, ann, cond)
  if ok_tm and txt then
    lines[#lines + 1] = "Tm        " .. txt .. (mm > 0 and "  [mismatches ignored]" or "")
    if #p.tail > 0 then
      local full = M.tm_text(p.seq, cond)
      if full then lines[#lines + 1] = "Tm full   " .. full:match("^%d+ °C") .. "  (whole primer, once the tail is copied)" end
    end
  end
  local ok_o, a = pcall(function() return oligo().assess(p.seq) end)
  if ok_o and a then
    lines[#lines + 1] = ("GC        %d %%   3' clamp %d G/C in last 5   longest run %d"):format(math.floor(a.gc * 100 + 0.5),
      a.gc_clamp, a.max_run)
    lines[#lines + 1] = "3' end    " .. fmt_dg(a.end_stability)
    -- only folds that are actually stabilising are worth reporting
    if a.hairpin and a.hairpin.dg < 0 then lines[#lines + 1] = "hairpin   " .. fmt_dg(a.hairpin.dg) end
    if a.homodimer and a.homodimer.dg < 0 then
      lines[#lines + 1] = "dimer     " .. fmt_dg(a.homodimer.dg) .. " with itself" ..
        (a.homodimer.three_prime and ", 3' end paired" or "")
    end
    for _, w in ipairs(a.warnings or {}) do lines[#lines + 1] = "warning   " .. w end
  end
  local sites = P.bind(rec.seq:tostring(), rec:is_circular(), p.seq, cfg().binding)
  lines[#lines + 1] = ("binds     %d site%s (3'-terminal %d bases exact)"):format(#sites, #sites == 1 and "" or "s",
    (cfg().binding or P.defaults).anchor or P.defaults.anchor)
  -- the facing primer that makes the shortest product, for Ta
  local best
  local pcr = require("plasmid.pcr")
  for _, q in ipairs(require("plasmid.primers").list(rec)) do
    if q.f ~= p.f and q.strand ~= p.strand then
      local prods = pcr.products(rec.seq:tostring(), rec:is_circular(), p.seq, q.seq, { binding = cfg().binding })
      if prods[1] and (not best or prods[1].length < best.len) then best = { q = q, len = prods[1].length } end
    end
  end
  if best then
    local ok, ta, t1, t2, info = pcall(tm().ta, ann, best.q.seq:sub(#best.q.tail + 1), cond)
    if ok and ta then
      lines[#lines + 1] = ("pair      with %s: product %d bp, Ta %d °C (Tm %d / %d)%s"):format(best.q.name, best.len, ta, t1, t2,
        info and info.note and ("  " .. info.note) or "")
      local ok2, het = pcall(function() return oligo().heterodimer(p.seq, best.q.seq) end)
      if ok2 and het and het.dg < 0 then
        lines[#lines + 1] = ("dimer     %s with %s%s"):format(fmt_dg(het.dg), best.q.name,
          het.three_prime and ", 3' end paired" or "")
      end
    end
  end
  return lines
end

function M.show_primer(buf, f)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  local p
  for _, q in ipairs(require("plasmid.primers").list(rec)) do if q.f == f then p = q end end
  if not p then return false end
  local lines = M.primer_lines(rec, p)
  local fbuf = vim.lsp.util.open_floating_preview(lines, "", {
    border = "rounded", max_width = math.min(100, vim.o.columns - 4), focus_id = "plasmid_primer",
  })
  vim.api.nvim_buf_set_extmark(fbuf, buffer().ns, 0, 0, { end_col = #lines[1], hl_group = "PlasmidPrimer" })
  for i, l in ipairs(lines) do
    if l:find("^warning") then vim.api.nvim_buf_set_extmark(fbuf, buffer().ns, i - 1, 0, { end_col = #l, hl_group = "WarningMsg" }) end
  end
  return true
end

--- :DnaTm [seq [seq2]] [polymerase=key] [dmso=n] [conc=nM] [salt=mM]; no
--- sequence: the selection's top strand.
function M.command(buf, args)
  local seqs, over = {}, {}
  for _, a in ipairs(args) do
    local k, v = a:match("^(%w+)=(.+)$")
    if k == "polymerase" or k == "pol" then over.polymerase = v
    elseif k == "dmso" then over.dmso = tonumber(v)
    elseif k == "conc" then over.conc_nM = tonumber(v)
    elseif k == "salt" then over.salt_mM = tonumber(v)
    elseif k then return notify("unknown option " .. k, vim.log.levels.ERROR)
    else seqs[#seqs + 1] = P.clean(a) end
  end
  if over.polymerase and not tm().polymerase(over.polymerase) then
    return notify("unknown polymerase " .. over.polymerase .. " (see :DnaPolymerase)", vim.log.levels.ERROR)
  end
  if #seqs == 0 and buf then
    local st = buffer().get(buf)
    local sel = st and require("plasmid.select").get(buf)
    if not sel then return notify("give a sequence or select one", vim.log.levels.WARN) end
    seqs[1] = P.clean(require("plasmid.select").top_strand(buffer().record(st), sel))
  end
  local cond = M.conditions(over)
  local out = {}
  for i, s in ipairs(seqs) do
    local ok, err = tm().valid(s)
    if not ok then return notify(("sequence %d: %s"):format(i, err), vim.log.levels.ERROR) end
    out[#out + 1] = ("Tm %s  %d nt  %s"):format(M.tm_text(s, cond), #s, #s > 30 and (s:sub(1, 27) .. "…") or s)
  end
  if #seqs == 2 then
    local ta, _, _, info = tm().ta(seqs[1], seqs[2], cond)
    out[#out + 1] = ("Ta %d °C%s"):format(ta, info and info.note and ("  (" .. info.note .. ")") or "")
  end
  notify(table.concat(out, "\n"))
  return out
end

--- :DnaPolymerase [key]: choose the polymerase used for Tm/Ta.
function M.choose(key)
  local o = require("plasmid.config").options
  o.primers = o.primers or {}
  local function set(k)
    o.primers.polymerase = k
    notify(("polymerase: %s"):format(M.polymerase_name(k)))
    for _, b in ipairs(vim.api.nvim_list_bufs()) do if buffer().get(b) then pcall(buffer().render, b) end end
  end
  if key then
    local p = tm().polymerase(key)
    if not p then return notify("unknown polymerase " .. key, vim.log.levels.ERROR) end
    return set(p.key)
  end
  local list = tm().POLYMERASES
  vim.ui.select(list, {
    prompt = "Polymerase for Tm / Ta",
    format_item = function(p)
      local cur = p.key == (o.primers.polymerase or "q5") and "●" or " "
      return ("%s %-18s %-48s %3d mM  %d nM"):format(cur, p.key, p.name, tm().BUFFERS[p.buffer] or 0, p.conc_nM)
    end,
  }, function(p) if p then set(p.key) end end)
end

function M.complete_polymerase(lead)
  local out = {}
  for _, p in ipairs(tm().POLYMERASES) do
    if p.key:find(lead:lower(), 1, true) == 1 then out[#out + 1] = p.key end
  end
  return out
end

return M
