-- Cloning simulation: :DnaClone, :DnaLigate, :DnaGoldenGate, :DnaGibson and
-- the staging list they draw their parts from.
--
-- The biology lives in plasmid.bio.assembly; this module only decides which
-- records become parts, runs the assembler, opens the product as a new record
-- and writes the assembly report. The report is the point of the exercise:
-- when an assembly fails it lists every part's ends and the reason each
-- candidate join was rejected, so the next step is obvious.
local A = require("plasmid.bio.assembly")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function digest() return require("plasmid.digest") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end
local function cfg() return require("plasmid.config").options.cloning or {} end

--- Parts staged with :DnaCloneAdd. Session-wide, not buffer-local: the parts
--- of one assembly usually come from several buffers.
--- Each entry: { name, rec, frag (nil for an intact circle), from }.
M.staged = {}

M.state = {}    -- report buffer -> { src, method, parts, product }

M.METHODS = { ligation = "Ligation", goldengate = "Golden Gate", gibson = "Gibson/HiFi" }

local ns = vim.api.nvim_create_namespace("plasmid_clone")

---------------------------------------------------------------------------
-- Parts
---------------------------------------------------------------------------

local function homology_opts()
  local c = cfg()
  return { min_homology = c.min_homology or 20, max_homology = c.max_homology or 200 }
end

--- Make the blunt fragment of a linear part; an intact circle keeps frag = nil
--- (it has no ends, so only Golden Gate can use it).
local function ensure_frag(p)
  if not p.frag and p.rec and not p.rec:is_circular() then
    p.frag = A.from_record(p.rec, p.name)
  end
  return p.frag
end

--- "5′ AATT … blunt", or "circular" for a part that still has to be cut.
function M.ends_text(p)
  if p.frag then return ("%s … %s"):format(A.end_text(p.frag, "left"), A.end_text(p.frag, "right")) end
  return "circular (no ends)"
end

--- The assembly fragment of a digest buffer's fragment `n` (default: the row
--- under the cursor). Returns the fragment, or nil and why not.
function M.digest_fragment(b, n)
  local D = digest().state[b]
  if not D then return nil, "not a digest table" end
  if not n then
    local win = vim.fn.bufwinid(b)
    if win == -1 then return nil, "the digest table is not on screen" end
    n = vim.api.nvim_win_get_cursor(win)[1] - D.first + 1
  end
  if not D.frags[n] then
    return nil, ("fragment %s is not in this digest (1..%d)"):format(tostring(n), #D.frags)
  end
  local frags = A.from_digest(D.rec, D.frags, D.rec.name)
  return frags[n], nil, D
end

--- :DnaCloneAdd [name] — stage the fragment under the cursor (digest table),
--- the selection, or the whole record of the current buffer.
function M.add(buf, name)
  if digest().state[buf] then
    local n = name and tonumber(name)
    local f, err, D = M.digest_fragment(buf, n)
    if not f then return notify(err .. " (put the cursor on a fragment row, or :DnaCloneAdd {n})", vim.log.levels.WARN) end
    if name and not n then f.name = name end
    return M.stage({ name = f.name, rec = f.rec, frag = f, from = ("digest of %s"):format(D.rec.name) })
  end
  local st = buffer().get(buf)
  if not st then
    return notify("nothing to stage here: open a sequence file, or run this from a digest table", vim.log.levels.WARN)
  end
  local rec = buffer().record(st)
  if rec:is_protein() then return notify("cloning needs a nucleotide sequence", vim.log.levels.WARN) end
  local sel = require("plasmid.select").get(buf)
  if sel then
    local S = require("plasmid.select")
    local copy = digest().extract(rec, { start = sel.start, stop = sel.stop, length = S.length(rec, sel) })
    copy.name = name or ("%s_%d-%d"):format(rec.name, sel.start, sel.stop)
    copy.description = ""
    return M.stage({ name = copy.name, rec = copy,
      from = ("%s %d..%d"):format(rec.name, sel.start, sel.stop) })
  end
  local copy = A.copy(rec)
  copy.name = name or rec.name
  copy.description = ""
  return M.stage({ name = copy.name, rec = copy, from = ("%s (%s, %d bp)"):format(rec.name, rec.topology, rec:len()) })
end

--- Append a part to the staging list.
function M.stage(p)
  ensure_frag(p)
  M.staged[#M.staged + 1] = p
  notify(("staged %s (%d bp, %s) — %d part%s staged"):format(p.name, p.rec:len(), M.ends_text(p),
    #M.staged, #M.staged == 1 and "" or "s"))
  return p
end

--- :DnaCloneClear
function M.clear()
  local n = #M.staged
  M.staged = {}
  notify(("cleared %d staged part%s"):format(n, n == 1 and "" or "s"))
  return n
end

--- A part argument: a staged part, "selection", a digest fragment number (in
--- a digest table), an open buffer, or a file path.
function M.resolve(buf, arg)
  for _, p in ipairs(M.staged) do
    if p.name:lower() == arg:lower() then return p end
  end
  if arg:lower() == "selection" then
    local st = buffer().get(buf)
    local sel = st and require("plasmid.select").get(buf)
    if not sel then return nil, "nothing is selected in this buffer (vaf, visual mode, or :DnaSelect)" end
    local rec = buffer().record(st)
    local S = require("plasmid.select")
    local copy = digest().extract(rec, { start = sel.start, stop = sel.stop, length = S.length(rec, sel) })
    copy.name = ("%s_%d-%d"):format(rec.name, sel.start, sel.stop)
    copy.description = ""
    return { name = copy.name, rec = copy, from = ("%s %d..%d"):format(rec.name, sel.start, sel.stop) }
  end
  local n = tonumber(arg)
  if n and digest().state[buf] then
    local f, err = M.digest_fragment(buf, n)
    if not f then return nil, err end
    return { name = f.name, rec = f.rec, frag = f, from = ("fragment %d of the digest"):format(n) }
  end
  local want = arg:lower()
  for b, st in pairs(buffer().state) do
    if vim.api.nvim_buf_is_valid(b) then
      local rec = buffer().record(st)
      local base = vim.fn.fnamemodify(st.path or "", ":t")
      if rec.name:lower() == want or base:lower() == want or vim.fn.fnamemodify(base, ":r"):lower() == want then
        local copy = A.copy(rec)
        copy.description = ""
        return { name = rec.name, rec = copy, from = ("buffer %s"):format(base ~= "" and base or rec.name) }
      end
    end
  end
  local path = vim.fn.expand(arg)
  if vim.fn.filereadable(path) == 1 then
    local res, err = require("plasmid.io").read(path)
    if not res or not res.records or #res.records == 0 then
      return nil, ("could not read %s: %s"):format(arg, tostring(err))
    end
    local rec = res.records[1]
    return { name = rec.name, rec = rec, from = vim.fn.fnamemodify(path, ":t") }
  end
  return nil, ("no part called %s: not staged, not an open buffer, not a file"):format(arg)
end

--- The parts a command works on: its arguments, else the staged list, else
--- the fragments of the current digest table. Returns parts and problems.
function M.parts(buf, args)
  local out, problems = {}, {}
  for _, a in ipairs(args or {}) do
    local p, err = M.resolve(buf, a)
    if p then ensure_frag(p); out[#out + 1] = p else problems[#problems + 1] = err end
  end
  if #args > 0 then return out, problems end
  if #M.staged > 0 then
    for _, p in ipairs(M.staged) do out[#out + 1] = p end
    return out, problems
  end
  local D = digest().state[buf]
  if D then
    for i, f in ipairs(A.from_digest(D.rec, D.frags, D.rec.name)) do
      out[i] = { name = f.name, rec = f.rec, frag = f, from = ("digest of %s"):format(D.rec.name) }
    end
    return out, problems
  end
  return out, { "no parts: stage them with :DnaCloneAdd, name them as arguments, or run this from a digest table" }
end

---------------------------------------------------------------------------
-- Method
---------------------------------------------------------------------------

local function has_site(rec, ename)
  local R = require("plasmid.bio.restriction")
  local e = R.get(ename)
  if not e then return false end
  return #R.sites(R.prepare(rec), e) > 0
end

--- The method :DnaClone picks: Golden Gate when every part carries a site for
--- one Type IIS enzyme, Gibson when two parts share terminal homology,
--- otherwise ligation. Returns method and (for Golden Gate) the enzyme.
function M.infer(parts)
  if #parts == 0 then return "ligation" end
  local names = { cfg().enzyme or "BsaI" }
  for _, n in ipairs(A.TYPE_IIS) do if n ~= names[1] then names[#names + 1] = n end end
  for _, n in ipairs(names) do
    local all = true
    for _, p in ipairs(parts) do
      if not has_site(p.rec, n) then all = false break end
    end
    if all then return "goldengate", n end
  end
  local frags, o = {}, homology_opts()
  for _, p in ipairs(parts) do if p.frag then frags[#frags + 1] = p.frag end end
  for i = 1, #frags do
    for j = 1, #frags do
      if i ~= j and A.homology(frags[i], frags[j], o) then return "gibson" end
    end
  end
  return "ligation"
end

--- Fragments for a method. Returns fragments and problems.
function M.fragments(parts, method, enzyme)
  if method == "goldengate" then
    local recs = {}
    for i, p in ipairs(parts) do recs[i] = { rec = p.rec, name = p.name } end
    local frags, problems = A.type_iis_parts(recs, enzyme or cfg().enzyme or "BsaI")
    return frags or {}, problems or {}
  end
  if method == "gibson" then
    local recs = {}
    for i, p in ipairs(parts) do recs[i] = { rec = p.rec, name = p.name } end
    return A.gibson_parts(recs)
  end
  local frags, problems = {}, {}
  for _, p in ipairs(parts) do
    local f = p.frag or ensure_frag(p)
    if not f then
      problems[#problems + 1] = ("%s is circular; cut it first (:DnaDigest, then :DnaCloneAdd) or use :DnaGoldenGate")
        :format(p.name)
    elseif f.circular then
      problems[#problems + 1] = ("%s is an intact circle; no enzyme cut it"):format(p.name)
    else
      frags[#frags + 1] = f
    end
  end
  return frags, problems
end

---------------------------------------------------------------------------
-- Report
---------------------------------------------------------------------------

local function pad(s, w) return s .. (" "):rep(math.max(0, w - vim.fn.strdisplaywidth(s))) end
local function rpad(s, w) return (" "):rep(math.max(0, w - vim.fn.strdisplaywidth(s))) .. s end

--- Rows of a small table: `right` lists the columns to right-align.
local function table_rows(header, rows, right)
  local w = {}
  for c = 1, #header do
    w[c] = vim.fn.strdisplaywidth(header[c])
    for _, r in ipairs(rows) do w[c] = math.max(w[c], vim.fn.strdisplaywidth(tostring(r[c] or ""))) end
  end
  local function fmt(r)
    local out = {}
    for c = 1, #header do
      local cell = tostring(r[c] or "")
      out[c] = (right and right[c]) and rpad(cell, w[c]) or (c < #header and pad(cell, w[c]) or cell)
    end
    return "  " .. table.concat(out, "  ")
  end
  local lines = { fmt(header) }
  for _, r in ipairs(rows) do lines[#lines + 1] = fmt(r) end
  return lines
end

local function short(seq, n)
  n = n or 24
  if #seq <= n then return seq end
  return seq:sub(1, n - 9) .. "…" .. seq:sub(#seq - 7)
end

--- One junction as text: "5′ AATG (4 bp)", "blunt", "overlap 25 bp ACGT…".
function M.junction_text(j)
  if not j then return "—" end
  if j.kind == "homology" then return ("overlap %d bp  %s"):format(j.len, short(j.seq)) end
  if j.kind == "blunt" or j.len == 0 then return "blunt" end
  return ("%s %s  (%d bp)"):format(j.five and "5′" or "3′", j.seq, j.len)
end

--- Features of the product that came from the parts (not our own part and
--- junction annotations).
function M.carried(rec)
  local out, seen, total = {}, {}, 0
  for _, f in ipairs(rec.features) do
    local note = f:get("note") or ""
    if f.type ~= "source" and not note:find("^part %d+ of ") and not note:find("junction: ") then
      total = total + 1
      local label = f:label()
      if vim.fn.strdisplaywidth(label) > 26 then label = label:sub(1, 25) .. "…" end
      local text = ("%s (%s)"):format(label, f.type)
      if not seen[text] then seen[text] = true; out[#out + 1] = text end
    end
  end
  return out, total
end

--- Why a's right end does not join b's left end, in the terms a cloner
--- thinks in. The engine's own reason is terse ("same strand"); here the two
--- ends are named and the mismatch is spelled out.
function M.why_text(a, b, method, why)
  if method == "gibson" then
    return tostring(why)
  end
  local ka, kb = A.end_kind(a, "right"), A.end_kind(b, "left")
  local la, lb = a.right.len, b.left.len
  local ta, tb = A.end_text(a, "right"), A.end_text(b, "left")
  if la ~= lb then
    return ("%s / %s — %s end and %s end"):format(ta, tb,
      la == 0 and "a blunt" or ("a %d base"):format(la), lb == 0 and "a blunt" or ("a %d base"):format(lb))
  end
  if ka ~= kb then
    return ("%s / %s — a %s overhang cannot anneal to a %s overhang"):format(ta, tb,
      ka == "5'" and "5′" or "3′", kb == "5'" and "5′" or "3′")
  end
  if A.end_seq(a, "right") ~= A.end_seq(b, "left") then
    return ("%s / %s — the overhangs are not complementary"):format(ta, tb)
  end
  return tostring(why)
end

--- Every ordered pair of ends, with the reason it does or does not join.
--- This is what a failed assembly reports. Returns the rows and how many
--- pairs of different parts do join.
function M.join_rows(frags, method)
  local o = homology_opts()
  local join = method == "gibson" and function(a, b) return A.homology(a, b, o) end or A.compatible
  local flip = {}
  for i, f in ipairs(frags) do flip[i] = A.flip_ends(f) end
  local rows, joining = {}, 0
  for i, a in ipairs(frags) do
    for j, b in ipairs(frags) do
      local tries = { { b, "" } }
      if i ~= j then tries[2] = { flip[j], " (flipped)" } end
      for _, t in ipairs(tries) do
        local jn, why = join(a, t[1])
        if jn and i ~= j then joining = joining + 1 end
        rows[#rows + 1] = {
          ("%s → %s%s"):format(a.name, b.name, t[2]),
          jn and ("joins: " .. M.junction_text(jn)) or M.why_text(a, t[1], method, why),
        }
      end
    end
  end
  return rows, joining
end

local function title_of(ctx)
  local label = M.METHODS[ctx.method] or ctx.method
  if ctx.method == "goldengate" then label = ("%s (%s)"):format(label, ctx.enzyme) end
  if ctx.product then
    return (" %s: %d part%s → %s, %d bp %s"):format(label, #ctx.asm.parts, #ctx.asm.parts == 1 and "" or "s",
      ctx.product.name, ctx.product:len(), ctx.product.topology)
  end
  return (" %s: %d part%s — nothing assembled"):format(label, #ctx.parts, #ctx.parts == 1 and "" or "s")
end

--- The report: title, parts, junctions, features, warnings — and, when the
--- assembly failed, every candidate join with the reason it was rejected.
--- Returns the lines and the highlights.
function M.report_lines(ctx)
  local lines, hl = { title_of(ctx) }, {}
  local function head(text)
    lines[#lines + 1] = ""
    lines[#lines + 1] = " " .. text
    hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "PlasmidHeader" }
  end
  hl[#hl + 1] = { 0, 0, #lines[1], "PlasmidTitle" }

  head("parts")
  local rows, right
  if ctx.asm then
    right = { true, false, true, false, false, false }
    local body = {}
    for i, p in ipairs(ctx.asm.parts) do
      body[i] = { tostring(i), p.frag.name, ("%d bp"):format(p.frag.rec:len()),
        p.flipped and "reverse" or "forward",
        A.end_text(p.frag, "left"), A.end_text(p.frag, "right") }
    end
    rows = table_rows({ "#", "part", "length", "strand", "left end", "right end" }, body, right)
  else
    right = { true, false, true, false, false, false }
    local body = {}
    for i, p in ipairs(ctx.parts) do
      local f = p.frag
      body[i] = { tostring(i), p.name, ("%d bp"):format(p.rec:len()), p.from or "",
        f and A.end_text(f, "left") or "—", f and A.end_text(f, "right") or "—" }
    end
    rows = table_rows({ "#", "part", "length", "from", "left end", "right end" }, body, right)
  end
  hl[#hl + 1] = { #lines, 0, #rows[1], "PlasmidHeader" }
  for _, l in ipairs(rows) do lines[#lines + 1] = l end

  if ctx.asm then
    head("junctions")
    local parts = ctx.asm.parts
    local body = {}
    for i, p in ipairs(parts) do
      local nxt = parts[i % #parts + 1]
      if p.junction and (i < #parts or ctx.asm.circular) then
        body[#body + 1] = { ("%s → %s"):format(p.frag.name, nxt.frag.name), M.junction_text(p.junction) }
      end
    end
    if #body == 0 then body[1] = { "—", "a single part, no junction" } end
    for _, l in ipairs(table_rows({ "from → to", "join" }, body)) do lines[#lines + 1] = l end

    head("features carried")
    local carried, total = M.carried(ctx.product)
    if total == 0 then
      lines[#lines + 1] = "  (none: the parts had no features)"
    else
      local shown = {}
      for i = 1, math.min(6, #carried) do shown[i] = carried[i] end
      lines[#lines + 1] = ("  %d feature%s: %s%s"):format(total, total == 1 and "" or "s",
        table.concat(shown, ", "), #carried > 6 and (", +%d more"):format(#carried - 6) or "")
    end
  else
    head("joins")
    local jrows, joining = M.join_rows(ctx.frags, ctx.method)
    local body, extra = {}, 0
    for i, r in ipairs(jrows) do
      if i <= 40 then body[#body + 1] = r else extra = extra + 1 end
    end
    for _, l in ipairs(table_rows({ "ends", "result" }, body)) do lines[#lines + 1] = l end
    if extra > 0 then lines[#lines + 1] = ("  +%d more pairs"):format(extra) end
    lines[#lines + 1] = ""
    local hint
    if ctx.method == "gibson" then
      hint = ("no pair of parts overlaps by %d bp or more; lower cloning.min_homology, or give the primers "
        .. "longer homology tails (:DnaPrimerDesign)"):format(homology_opts().min_homology)
    elseif joining == 0 then
      hint = "no two parts have ends that anneal; cut them with enzymes that leave matching overhangs, "
        .. "or blunt both ends"
    else
      hint = ("%d pair%s of ends anneal, but no order uses every part%s"):format(joining,
        joining == 1 and "" or "s", #ctx.frags > 1 and " and closes" or "")
    end
    lines[#lines + 1] = "  " .. hint
    hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "ErrorMsg" }
  end

  if #(ctx.warnings or {}) > 0 then
    head("warnings")
    for _, w in ipairs(ctx.warnings) do
      lines[#lines + 1] = "  " .. w
      hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], ctx.product and "WarningMsg" or "ErrorMsg" }
    end
  end
  lines[#lines + 1] = ""
  lines[#lines + 1] = " q close"
  hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "Comment" }
  return lines, hl
end

local function report_buffer(src)
  for b, s in pairs(M.state) do
    if s.src == src and vim.api.nvim_buf_is_valid(b) then return b end
  end
  local b = vim.api.nvim_create_buf(false, true)
  vim.bo[b].bufhidden = "wipe"
  vim.bo[b].filetype = "plasmid_clone"
  pcall(vim.api.nvim_buf_set_name, b, "plasmid://clone/" .. tostring(src or 0))
  vim.api.nvim_create_autocmd("BufWipeout", { buffer = b, once = true, callback = function() M.state[b] = nil end })
  return b
end

--- Show lines in the cloning report split for `src`.
function M.show(src, lines, hl, state)
  local b = report_buffer(src)
  vim.bo[b].modifiable = true
  vim.api.nvim_buf_set_lines(b, 0, -1, false, lines)
  vim.bo[b].modifiable = false
  vim.api.nvim_buf_clear_namespace(b, ns, 0, -1)
  for _, h in ipairs(hl or {}) do
    pcall(vim.api.nvim_buf_set_extmark, b, ns, h[1], h[2], { end_col = h[3], hl_group = h[4] })
  end
  M.state[b] = vim.tbl_extend("force", { src = src }, state or {})
  M.last = b
  local win = vim.fn.bufwinid(b)
  if win == -1 then
    vim.cmd("botright " .. math.min(#lines + 1, math.floor(vim.o.lines * 0.55)) .. "split")
    vim.api.nvim_win_set_buf(0, b)
    win = vim.api.nvim_get_current_win()
    vim.wo[win].wrap = false
    vim.wo[win].number = false
    vim.wo[win].relativenumber = false
  end
  vim.api.nvim_win_set_cursor(win, { 1, 0 })
  vim.keymap.set("n", "q", "<cmd>close<cr>", { buffer = b, nowait = true, silent = true })
  return b
end

--- :DnaCloneList
function M.list(src)
  if #M.staged == 0 then
    notify("no parts staged (:DnaCloneAdd stages the record, the selection or a digest fragment)")
    return nil
  end
  local body = {}
  for i, p in ipairs(M.staged) do
    body[i] = { tostring(i), p.name, ("%d bp"):format(p.rec:len()),
      p.frag and A.end_text(p.frag, "left") or "—", p.frag and A.end_text(p.frag, "right") or "—",
      p.from or "" }
  end
  local lines = { (" Staged parts (%d) — :DnaClone assembles them, :DnaCloneClear empties the list"):format(#M.staged), "" }
  local hl = { { 0, 0, #lines[1], "PlasmidTitle" } }
  local rows = table_rows({ "#", "part", "length", "left end", "right end", "source" }, body,
    { true, false, true, false, false, false })
  hl[#hl + 1] = { #lines, 0, #rows[1], "PlasmidHeader" }
  for _, l in ipairs(rows) do lines[#lines + 1] = l end
  lines[#lines + 1] = ""
  lines[#lines + 1] = " q close"
  hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "Comment" }
  return M.show(src, lines, hl, { method = "staged" }), lines
end

---------------------------------------------------------------------------
-- Running an assembly
---------------------------------------------------------------------------

local function sanitize(s) return (s:gsub("[^%w%._%-]", "_")):sub(1, 40) end

--- Product name: the parts in assembly order, with the enzyme after the first
--- one for Golden Gate (pUC19_BsaI_insert).
function M.product_name(asm, method, enzyme)
  local names = {}
  for i, p in ipairs(asm.parts) do names[i] = p.frag.name end
  if method == "goldengate" then
    table.insert(names, math.min(2, #names + 1), enzyme or cfg().enzyme or "BsaI")
  elseif method == "gibson" then
    names[#names + 1] = "gibson"
  elseif #names == 1 then
    names[2] = "ligation"
  end
  return sanitize(table.concat(names, "_"))
end

local function describe(asm, method, enzyme, rec)
  local names = {}
  for i, p in ipairs(asm.parts) do names[i] = p.frag.name .. (p.flipped and " (reverse)" or "") end
  local label = M.METHODS[method] or method
  if method == "goldengate" then label = ("%s (%s)"):format(label, enzyme) end
  return ("%s of %s — %d bp %s"):format(label, table.concat(names, " + "), rec:len(), rec.topology)
end

local function assembly_item(asm)
  local names = {}
  for i, p in ipairs(asm.parts) do names[i] = p.frag.name .. (p.flipped and "◂" or "▸") end
  local len = 0
  for _, p in ipairs(asm.parts) do len = len + p.frag.rec:len() - (p.junction and p.junction.len or 0) end
  return ("%7d bp %-9s %s"):format(len, asm.circular and "circular" or "linear", table.concat(names, " "))
end

--- Build the product of one assembly, open it and write the report.
function M.finish(src, ctx, asm)
  local c = cfg()
  ctx.asm = asm
  local name = M.product_name(asm, ctx.method, ctx.enzyme)
  local ok, rec = pcall(A.build, asm, {
    name = name,
    annotate_parts = c.annotate_parts ~= false,
    annotate_junctions = c.annotate_junctions ~= false,
  })
  if not ok then
    notify("could not build the product: " .. tostring(rec), vim.log.levels.ERROR)
    return nil
  end
  rec.description = describe(asm, ctx.method, ctx.enzyme, rec)
  ctx.product = rec
  if not asm.circular then
    table.insert(ctx.warnings, 1, "the parts do not close into a circle; the product is linear")
  end
  if ctx.nsolutions and ctx.nsolutions > 1 then
    table.insert(ctx.warnings, 1, ("%d assemblies are possible with these parts; this is the one you picked")
      :format(ctx.nsolutions))
  end
  local nb = digest().open_record(src, rec, ("%s — %d bp %s"):format(M.METHODS[ctx.method] or ctx.method,
    rec:len(), rec.topology))
  local lines, hl = M.report_lines(ctx)
  local pwin = nb and vim.fn.bufwinid(nb) or -1
  M.show(src, lines, hl, { method = ctx.method, parts = ctx.parts, product = rec })
  if pwin ~= -1 then vim.api.nvim_set_current_win(pwin) end
  return nb, rec, ctx
end

--- Run an assembly. `method` is "ligation" | "goldengate" | "gibson" | nil
--- (inferred). Returns the product buffer and record, or nil.
function M.run(buf, method, args, opts)
  opts = opts or {}
  local src = buffer().get(buf) and buf or (digest().state[buf] and digest().state[buf].src) or nil
  local parts, problems = M.parts(buf, args or {})
  if #(args or {}) > 0 and #problems > 0 then
    -- a named part that is not there is a typo, not a biological result
    return notify(table.concat(problems, "; "), vim.log.levels.ERROR)
  end
  if #parts == 0 then
    return notify(problems[1] or "no parts to assemble", vim.log.levels.WARN)
  end
  local enzyme = opts.enzyme
  if not method then
    method, enzyme = M.infer(parts)
    enzyme = opts.enzyme or enzyme
  end
  if method == "goldengate" then enzyme = enzyme or cfg().enzyme or "BsaI" end
  local frags, fproblems = M.fragments(parts, method, enzyme)
  for _, p in ipairs(fproblems) do problems[#problems + 1] = p end
  local ctx = { method = method, enzyme = enzyme, parts = parts, frags = frags, warnings = problems }
  local label = M.METHODS[method] or method
  local function fail(msg)
    local lines, hl = M.report_lines(ctx)
    M.show(src, lines, hl, { method = method, parts = parts })
    notify(("%s: %s — see the report"):format(label, msg), vim.log.levels.WARN)
    return nil
  end
  -- A part that yields no fragment is a dead end, not a reason to assemble
  -- what is left: assembling fewer parts than asked for would be a lie.
  local used, dropped = {}, {}
  for _, f in ipairs(frags) do used[f.name] = true end
  for _, p in ipairs(parts) do
    if not used[p.name] then
      dropped[#dropped + 1] = p.name
      problems[#problems + 1] = ("%s contributed no fragment and was not used"):format(p.name)
    end
  end
  if #frags == 0 then return fail("no part could be prepared") end
  if #dropped > 0 then
    return fail(("%d of %d parts could not be prepared (%s)"):format(#dropped, #parts, table.concat(dropped, ", ")))
  end
  local c = cfg()
  local aopts = vim.tbl_extend("force", homology_opts(), {
    join = method == "gibson" and "homology" or "ends",
    max_solutions = c.max_solutions or 50,
    circular = opts.circular,
  })
  local ok, asms = pcall(A.assemble, frags, aopts)
  if not ok then
    notify("could not assemble: " .. tostring(asms), vim.log.levels.ERROR)
    return nil
  end
  if #asms == 0 then
    local lines, hl = M.report_lines(ctx)
    M.show(src, lines, hl, { method = method, parts = parts })
    notify(("%s: the %d parts do not assemble — the report lists every end and why each join failed"):format(
      M.METHODS[method] or method, #frags), vim.log.levels.WARN)
    return nil
  end
  ctx.nsolutions = #asms
  if #asms == 1 or opts.first then return M.finish(src, ctx, asms[1]) end
  vim.ui.select(asms, {
    prompt = ("%d possible assemblies"):format(#asms),
    format_item = assembly_item,
  }, function(choice)
    if choice then M.finish(src, ctx, choice) end
  end)
  return nil, nil, ctx
end

return M
