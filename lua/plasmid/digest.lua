-- Simulated restriction digests: :DnaDigest, a fragment table (and gel) in a
-- split, fragment selection and extraction as a new record.
local R = require("plasmid.bio.restriction")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end

M.state = {}   -- digest buffer -> { src, rec, frags, names, revision, rows }

local function seq_at(rec, a, b)
  local L = rec:len()
  if a <= b then return rec.seq:sub(a, b) end
  return rec.seq:sub(a, L) .. rec.seq:sub(1, b)
end

--- Overhang text of a cut: "blunt", "5′ AATT", "3′ TGCA" (top-strand bases).
function M.end_text(rec, site)
  local o = site.overhang
  if o == 0 then return "blunt" end
  local L = rec:len()
  local function w(p) return (p - 1) % L + 1 end
  if o > 0 then return "5′ " .. seq_at(rec, w(site.cut), w(site.cut + o - 1)):upper() end
  return "3′ " .. seq_at(rec, w(site.cut_bottom), w(site.cut_bottom - o - 1)):upper()
end

local function end_label(rec, sites)
  if not sites then return "(end)" end
  local s = sites[1]
  local names = {}
  for _, x in ipairs(sites) do names[#names + 1] = x.enzyme.name end
  return ("%s %s"):format(table.concat(names, "/"), M.end_text(rec, s))
end

--- Features in a fragment: whole ones by label, cut ones with a trailing "…".
function M.fragment_features(rec, fr)
  local L = rec:len()
  local function off(p) return (p - fr.start) % L end
  local whole, cut = {}, {}
  local vf = require("plasmid.view.features")
  local opts = require("plasmid.config").options
  local seen = {}
  for _, it in ipairs(vf.by_priority(vf.visible(rec, opts), opts)) do
    local f = it.f
    if f.type ~= "source" and f.loc then
      local any, all = false, true
      for _, e in ipairs(rec:extents(f)) do
        local a, b = e[1], e[2]
        -- does the extent [a, b] fall inside the fragment?
        local oa, ob = off(a), off(b)
        local inside = oa < fr.length and ob < fr.length and oa <= ob
        if inside then any = true else
          all = false
          -- partial overlap: either end inside, or the extent covers the fragment start
          if oa < fr.length or ob < fr.length then any = true end
          local covers = (a <= b and fr.start >= a and fr.start <= b) or (a > b and (fr.start >= a or fr.start <= b))
          if covers then any = true end
        end
      end
      local name = it.label or f:label()
      if any and not seen[name] then
        seen[name] = true
        if all then whole[#whole + 1] = name else cut[#cut + 1] = name .. "…" end
      end
    end
  end
  local out = {}
  for _, x in ipairs(whole) do out[#out + 1] = x end
  for _, x in ipairs(cut) do out[#out + 1] = x end
  return out
end

--- Digest a record with enzymes (names). Returns fragments (position order),
--- the cuts and the enzymes with no site.
function M.compute(rec, names)
  local P = R.prepare(rec)
  local cuts, none, unknown = {}, {}, {}
  for _, n in ipairs(names) do
    local e = R.get(n)
    if not e then unknown[#unknown + 1] = n
    else
      local s = R.sites(P, e)
      if #s == 0 then none[#none + 1] = e.name end
      for _, x in ipairs(s) do cuts[#cuts + 1] = x end
    end
  end
  return R.digest(rec, cuts), cuts, none, unknown
end

local HEADER = { "#", "from", "to", "length", "left end", "right end", "features" }

local function table_lines(rec, frags)
  local rows = {}
  for i, fr in ipairs(frags) do
    rows[i] = {
      tostring(i), tostring(fr.start), tostring(fr.stop), tostring(fr.length) .. " bp",
      fr.uncut and "—" or end_label(rec, fr.left), fr.uncut and "—" or end_label(rec, fr.right),
      (function()
        local list = M.fragment_features(rec, fr)
        local shown = {}
        for k = 1, math.min(5, #list) do shown[k] = list[k] end
        local t = table.concat(shown, ", ")
        if #list > 5 then t = t .. (", +%d more"):format(#list - 5) end
        return t
      end)(),
    }
  end
  local w = {}
  for c = 1, #HEADER do
    w[c] = vim.fn.strdisplaywidth(HEADER[c])
    for _, r in ipairs(rows) do w[c] = math.max(w[c], vim.fn.strdisplaywidth(r[c])) end
  end
  local function fmt(r)
    local parts = {}
    for c = 1, #r do
      local cell = r[c]
      local pad = w[c] - vim.fn.strdisplaywidth(cell)
      if c <= 4 then parts[c] = (" "):rep(pad) .. cell else parts[c] = cell .. (c < #r and (" "):rep(pad) or "") end
    end
    return " " .. table.concat(parts, "  ")
  end
  local out = { fmt(HEADER) }
  for _, r in ipairs(rows) do out[#out + 1] = fmt(r) end
  return out
end

local function digest_buffer(src)
  for b, s in pairs(M.state) do
    if s.src == src and vim.api.nvim_buf_is_valid(b) then return b end
  end
  local b = vim.api.nvim_create_buf(false, true)
  vim.bo[b].bufhidden = "wipe"
  vim.bo[b].filetype = "plasmid_digest"
  pcall(vim.api.nvim_buf_set_name, b, "plasmid://digest/" .. src)
  vim.api.nvim_create_autocmd("BufWipeout", { buffer = b, once = true, callback = function() M.state[b] = nil end })
  return b
end

local ns = vim.api.nvim_create_namespace("plasmid_digest")

--- :DnaDigest — names default to the shown named set.
function M.run(src, names)
  local st = buffer().get(src)
  if not st then return end
  local rec = buffer().record(st)
  if rec:is_protein() then return notify("digests need a nucleotide sequence", vim.log.levels.WARN) end
  if not names or #names == 0 then
    if st.enzymes and st.enzymes.spec.mode == "names" then names = st.enzymes.spec.names
    else return notify("name the enzymes, e.g. :DnaDigest EcoRI HindIII (or pick them with gZ first)", vim.log.levels.WARN) end
  end
  local frags, cuts, none, unknown = M.compute(rec, names)
  if #unknown > 0 then return notify("unknown enzyme: " .. table.concat(unknown, ", "), vim.log.levels.ERROR) end
  local b = digest_buffer(src)
  local title = (" Digest of %s (%d bp, %s) with %s — %d cut%s, %d fragment%s"):format(
    rec.name, rec:len(), rec.topology, table.concat(names, " + "), #cuts, #cuts == 1 and "" or "s",
    #frags, #frags == 1 and "" or "s")
  local lines = { title }
  if #none > 0 then lines[#lines + 1] = " no site: " .. table.concat(none, ", ") end
  lines[#lines + 1] = ""
  local first_row = #lines + 1   -- the column header row
  for _, l in ipairs(table_lines(rec, frags)) do lines[#lines + 1] = l end
  lines[#lines + 1] = ""
  lines[#lines + 1] = " <CR> select fragment   e extract as new record   q close"
  local gel_at = #lines
  local gel = require("plasmid.gel")
  local glines, ghl = gel.render({ { label = table.concat(names, "+"), sizes = vim.tbl_map(function(f) return f.length end, frags) } },
    require("plasmid.config").options.gel or {})
  lines[#lines + 1] = ""
  for _, l in ipairs(glines) do lines[#lines + 1] = l end

  vim.bo[b].modifiable = true
  vim.api.nvim_buf_set_lines(b, 0, -1, false, lines)
  vim.bo[b].modifiable = false
  vim.api.nvim_buf_clear_namespace(b, ns, 0, -1)
  vim.api.nvim_buf_set_extmark(b, ns, 0, 0, { end_col = #lines[1], hl_group = "PlasmidTitle" })
  vim.api.nvim_buf_set_extmark(b, ns, first_row - 1, 0, { end_col = #lines[first_row], hl_group = "PlasmidHeader" })
  vim.api.nvim_buf_set_extmark(b, ns, gel_at - 1, 0, { end_col = #lines[gel_at], hl_group = "Comment" })
  for _, h in ipairs(ghl) do
    vim.api.nvim_buf_set_extmark(b, ns, gel_at + 1 + h[1], h[2], { end_col = h[3], hl_group = h[4] })
  end
  M.state[b] = { src = src, names = names, frags = frags, first = first_row + 1, revision = rec.seq_revision or 0, rec = rec }

  local win = vim.fn.bufwinid(b)
  if win == -1 then
    local srcwin = vim.fn.bufwinid(src)
    if srcwin ~= -1 then vim.api.nvim_set_current_win(srcwin) end
    vim.cmd("botright " .. math.min(#lines + 1, math.floor(vim.o.lines * 0.6)) .. "split")
    vim.api.nvim_win_set_buf(0, b)
    win = vim.api.nvim_get_current_win()
    vim.wo[win].wrap = false
    vim.wo[win].number = false
    vim.wo[win].relativenumber = false
  end
  vim.api.nvim_win_set_cursor(win, { first_row + 1, 0 })
  local o = { buffer = b, nowait = true, silent = true }
  vim.keymap.set("n", "<CR>", function() M.select(b) end, vim.tbl_extend("force", o, { desc = "select fragment" }))
  vim.keymap.set("n", "e", function() M.extract_under_cursor(b) end, vim.tbl_extend("force", o, { desc = "extract fragment" }))
  vim.keymap.set("n", "q", "<cmd>close<cr>", o)
  return b, frags
end

local function fragment_under_cursor(b)
  local D = M.state[b]
  if not D then return end
  local i = vim.api.nvim_win_get_cursor(0)[1] - D.first + 1
  local fr = D.frags[i]
  if not fr then return notify("cursor is not on a fragment row", vim.log.levels.WARN) end
  local st = buffer().get(D.src)
  if not st then return notify("the digested buffer is gone", vim.log.levels.WARN) end
  if buffer().record(st) ~= D.rec or (D.rec.seq_revision or 0) ~= D.revision then
    return notify("the sequence changed since this digest; run :DnaDigest again", vim.log.levels.WARN)
  end
  return fr, D, st
end

--- <CR>: select the fragment in the plasmid view.
function M.select(b)
  local fr, D = fragment_under_cursor(b)
  if not fr then return end
  local sel = require("plasmid.select")
  local srcwin = vim.fn.bufwinid(D.src)
  if srcwin == -1 then return end
  vim.api.nvim_set_current_win(srcwin)
  local st = buffer().get(D.src)
  sel.set(D.src, sel.make(buffer().record(st), fr.start, fr.stop))
  buffer().goto_bp(srcwin, fr.start)
end

--- A new linear record holding the fragment (top strand, from the left
--- top-strand cut to the base before the right one). Features are carried
--- over by the editing machinery: whole ones intact, cut ones truncated.
function M.extract(rec, fr)
  local pio = require("plasmid.io")
  local T = require("plasmid.model.transaction")
  local text = pio.format({ rec }, "genbank")
  local copy = assert(pio.parse(text, nil, "genbank")).records[1]
  local L = copy:len()
  if fr.start ~= 1 then
    if copy:is_circular() then T.apply(copy, { op = "rotate", origin = fr.start })
    else T.apply(copy, { op = "splice", start = 1, stop = fr.start - 1, ins = "" }) end
  end
  if fr.length < copy:len() then
    T.apply(copy, { op = "splice", start = fr.length + 1, stop = copy:len(), ins = "" })
  end
  if copy.topology ~= "linear" then T.apply(copy, { op = "topology", value = "linear" }) end
  copy.history = nil
  local function ends(sites) return sites and sites[1].enzyme.name or "end" end
  copy.name = (rec.name .. "_" .. fr.start .. "-" .. fr.stop):sub(1, 40)
  copy.description = ("%s–%s fragment (%d..%d, %d bp) of %s"):format(ends(fr.left), ends(fr.right),
    fr.start, fr.stop, fr.length, rec.name)
  return copy, L
end

--- Open a new record (not yet on disk) in a plasmid buffer next to `src`.
--- opts: { dir = directory for the proposed file name (default: next to src,
--- else the working directory), ext = "gb" | "gp", file = file-level data
--- from the parser (kept so :w can write the text back verbatim) }.
--- An existing file is never reused: NAME_2.gb, NAME_3.gb, ... are tried.
function M.open_record(src, newrec, quiet_msg, opts)
  opts = opts or {}
  local st = src and buffer().get(src)
  local dir = opts.dir or (st and st.path and vim.fn.fnamemodify(st.path, ":h")) or vim.fn.getcwd()
  local base = newrec.name:gsub("[^%w%._%-]", "_")
  local ext = opts.ext or "gb"
  local path = dir .. "/" .. base .. "." .. ext
  local k = 1
  while vim.fn.filereadable(path) == 1 or vim.fn.bufexists(path) == 1 do
    k = k + 1
    path = ("%s/%s_%d.%s"):format(dir, base, k, ext)
  end
  local srcwin = src and vim.fn.bufwinid(src) or -1
  if srcwin ~= -1 then vim.api.nvim_set_current_win(srcwin) end
  vim.cmd("edit " .. vim.fn.fnameescape(path))
  local nb = vim.api.nvim_get_current_buf()
  local ns_ = buffer().get(nb)
  if not ns_ then notify("could not open a buffer for the new record", vim.log.levels.ERROR) return nil end
  newrec.format = newrec.format or "genbank"
  ns_.records = { newrec }
  ns_.current = 1
  if opts.file then ns_.file = opts.file end
  ns_.format = "genbank"
  ns_.saved = { -1 }   -- nothing on disk yet
  ns_.view = newrec:is_protein() and "linear" or ns_.view
  buffer().mark_modified(nb)
  buffer().render(nb, { keep_bp = 1 })
  if quiet_msg then notify(quiet_msg .. " — opened as " .. vim.fn.fnamemodify(path, ":~:.") .. "; :w to save") end
  return nb, path
end

function M.extract_under_cursor(b)
  local fr, D, st = fragment_under_cursor(b)
  if not fr then return end
  local rec = buffer().record(st)
  local frag = M.extract(rec, fr)
  local nb = M.open_record(D.src, frag, ("fragment %d..%d (%d bp)"):format(fr.start, fr.stop, fr.length))
  return nb, frag
end

--- :DnaGel "EcoRI | EcoRI HindIII | uncut [ladder=100bp]": one lane per digest.
function M.gel(src, spec)
  local st = buffer().get(src)
  if not st then return end
  local rec = buffer().record(st)
  local opts = vim.deepcopy(require("plasmid.config").options.gel or {})
  spec = spec:gsub("ladder=(%S+)", function(l) opts.ladder = l; return "" end)
  local lanes = {}
  for part in (spec .. "|"):gmatch("([^|]*)|") do
    local names = vim.split(vim.trim(part), "%s+", { trimempty = true })
    if #names > 0 then
      if #names == 1 and names[1]:lower() == "uncut" then
        lanes[#lanes + 1] = { label = "uncut", sizes = { rec:len() },
          note = rec:is_circular() and "drawn at its linear size; supercoiled and nicked forms run differently" or nil }
      else
        local frags, _, _, unknown = M.compute(rec, names)
        if #unknown > 0 then return notify("unknown enzyme: " .. table.concat(unknown, ", "), vim.log.levels.ERROR) end
        lanes[#lanes + 1] = { label = table.concat(names, "+"), sizes = vim.tbl_map(function(f) return f.length end, frags) }
      end
    end
  end
  if #lanes == 0 then return notify("no lanes; e.g. :DnaGel EcoRI | EcoRI HindIII", vim.log.levels.WARN) end
  local glines, ghl, info = require("plasmid.gel").render(lanes, opts)
  local b = vim.api.nvim_create_buf(false, true)
  vim.bo[b].bufhidden = "wipe"
  vim.bo[b].filetype = "plasmid_gel"
  local lines = { (" Gel: %s (%d bp, %s)"):format(rec.name, rec:len(), rec.topology), "" }
  for _, l in ipairs(glines) do lines[#lines + 1] = l end
  vim.api.nvim_buf_set_lines(b, 0, -1, false, lines)
  vim.bo[b].modifiable = false
  vim.api.nvim_buf_set_extmark(b, ns, 0, 0, { end_col = #lines[1], hl_group = "PlasmidTitle" })
  for _, h in ipairs(ghl) do vim.api.nvim_buf_set_extmark(b, ns, 2 + h[1], h[2], { end_col = h[3], hl_group = h[4] }) end
  vim.cmd("botright " .. math.min(#lines + 1, vim.o.lines - 4) .. "split")
  vim.api.nvim_win_set_buf(0, b)
  vim.wo.wrap = false
  vim.keymap.set("n", "q", "<cmd>close<cr>", { buffer = b, nowait = true, silent = true })
  return b, lanes, info
end

return M
