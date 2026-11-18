-- Adding, editing and deleting features; undo/redo; the personal library.
local config = require("plasmid.config")
local location = require("plasmid.model.location")
local record = require("plasmid.model.record")
local tx = require("plasmid.model.transaction")
local select_ = require("plasmid.select")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end

-- Common GenBank feature keys, for completion.
M.FEATURE_TYPES = {
  "CDS", "gene", "promoter", "terminator", "rep_origin", "RBS", "primer_bind", "protein_bind",
  "misc_feature", "misc_binding", "regulatory", "enhancer", "polyA_signal", "sig_peptide",
  "mat_peptide", "mRNA", "exon", "intron", "ncRNA", "tRNA", "rRNA", "misc_RNA", "repeat_region",
  "mobile_element", "LTR", "oriT", "5'UTR", "3'UTR", "stem_loop", "variation", "misc_difference",
}

function M.complete_types(arglead)
  local out = {}
  for _, t in ipairs(M.FEATURE_TYPES) do
    if t:lower():find(arglead:lower(), 1, true) == 1 then out[#out + 1] = t end
  end
  return out
end

---------------------------------------------------------------------------
-- Pure helpers
---------------------------------------------------------------------------

--- Location text for a selection on a record with a strand (1, -1, or 0).
function M.location_for(rec, sel, strand)
  local nodes = {}
  for _, r in ipairs(select_.ranges(rec, sel)) do
    if r[1] == r[2] then
      nodes[#nodes + 1] = { kind = "site", pos = r[1] }
    else
      nodes[#nodes + 1] = { kind = "range", start = r[1], stop = r[2] }
    end
  end
  local node = #nodes == 1 and nodes[1] or { kind = "join", parts = nodes }
  if strand == -1 then node = { kind = "complement", child = node } end
  return location.format(node)
end

--- Build a new Feature from simple fields.
---@param t { type: string, label: string, sel: table, strand: integer, color?: string, note?: string }
function M.make_feature(rec, t)
  local quals = { { key = "label", value = t.label, quoted = true } }
  if t.note and t.note ~= "" then quals[#quals + 1] = { key = "note", value = t.note, quoted = true } end
  if t.color then
    quals[#quals + 1] = { key = "ApEinfo_fwdcolor", value = t.color, quoted = true }
    quals[#quals + 1] = { key = "ApEinfo_revcolor", value = t.color, quoted = true }
  end
  local f = record.Feature.new({
    type = t.type, location_text = M.location_for(rec, t.sel, t.strand),
    qualifiers = quals, color = t.color, no_arrow = t.strand == 0,
  })
  f.pinned = true
  return f
end

--- Text shown in the feature editor.
function M.editor_text(rec, f, header)
  local out = {
    "# " .. header,
    ("# Record %s, %d %s, %s.  :w applies, :q discards."):format(rec.name, rec:len(), rec:unit(), rec.topology),
    "# First line: TYPE LOCATION (GenBank syntax, e.g. complement(100..200) or join(1..50,90..120)).",
    "# Then one /qualifier=value per line. Indentation is free.",
    "",
  }
  out[#out + 1] = (f and f.type or "misc_feature") .. " " .. (f and f.location_text or "")
  for _, q in ipairs(f and f.qualifiers or { { key = "label", value = "", quoted = true } }) do
    if q.value == true then
      out[#out + 1] = "/" .. q.key
    elseif q.quoted then
      out[#out + 1] = ("/%s=\"%s\""):format(q.key, tostring(q.value):gsub('"', '""'))
    else
      out[#out + 1] = ("/%s=%s"):format(q.key, tostring(q.value))
    end
  end
  return out
end

--- Parse editor text into a Feature, validating against the record.
---@return table|nil feature, string|nil err
function M.parse_editor(rec, lines)
  local body, first = {}, true
  for _, raw in ipairs(lines) do
    local l = raw:gsub("%s+$", "")
    if not l:match("^%s*#") and l:match("%S") then
      local t = l:gsub("^%s+", "")
      if first then
        local key, loc = t:match("^(%S+)%s*(.*)$")
        body[#body + 1] = "     " .. key .. (" "):rep(math.max(1, 16 - #key)) .. loc
        first = false
      else
        body[#body + 1] = (" "):rep(21) .. t
      end
    end
  end
  if first then return nil, "empty feature" end
  local text = table.concat({
    "LOCUS       editor 1 bp DNA linear UNA 01-JAN-1980",
    "FEATURES             Location/Qualifiers",
    table.concat(body, "\n"),
    "ORIGIN",
    "//",
  }, "\n") .. "\n"
  local res, err = require("plasmid.io.genbank").parse(text)
  if not res then return nil, err end
  local feats = res.records[1].features
  if #feats ~= 1 then return nil, ("expected one feature, found %d"):format(#feats) end
  local f = feats[1]
  if not f.loc then return nil, ("invalid location %q: %s"):format(f.location_text, tostring(f.loc_error)) end
  local L = rec:len()
  for _, p in ipairs(f.parts) do
    if not p.ref and (p.start < 1 or p.stop > L or (p.stop < p.start and p.kind ~= "between")) then
      return nil, ("location %s is outside 1..%d"):format(f.location_text, L)
    end
  end
  f.raw_lines = nil
  f.id = nil
  local c = f:get("ApEinfo_fwdcolor")
  if type(c) == "string" and c:match("^#%x%x%x%x%x%x$") then f.color = c end
  return f
end

---------------------------------------------------------------------------
-- Buffer-bound actions
---------------------------------------------------------------------------

local function current()
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer().get(buf)
  if not st then notify("not a sequence buffer", vim.log.levels.WARN) return nil end
  return buf, st, buffer().record(st)
end

--- Commit an op on the current record, refresh view and modified state.
function M.commit(buf, op, label)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  tx.commit(rec, op, label)
  buffer().refresh_modified(buf)
  require("plasmid.orfs").refresh(buf)
  buffer().render(buf)
end

--- Base pair where a sequence op happened (for placing the cursor).
function M.edit_site(op)
  if not op then return nil end
  if op.op == "batch" then return M.edit_site(op.ops[1]) end
  if op.op == "rotate" then return 1 end
  return op.start
end

function M.undo()
  local buf, st, rec = current()
  if not buf then return end
  local c = tx.undo(rec)
  if not c then return notify("already at oldest change") end
  buffer().refresh_modified(buf)
  local seq = require("plasmid.seqedit").is_seq_op(c.forward)
  if seq then require("plasmid.seqedit").after_change(buf) else require("plasmid.orfs").refresh(buf) end
  buffer().render(buf, { keep_bp = seq and M.edit_site(c.forward) or nil })
  notify("undid: " .. c.label)
end

function M.redo()
  local buf, st, rec = current()
  if not buf then return end
  local c = tx.redo(rec)
  if not c then return notify("already at newest change") end
  buffer().refresh_modified(buf)
  local seq = require("plasmid.seqedit").is_seq_op(c.forward)
  if seq then require("plasmid.seqedit").after_change(buf) else require("plasmid.orfs").refresh(buf) end
  buffer().render(buf, { keep_bp = seq and M.edit_site(c.forward) or nil })
  notify("redid: " .. c.label)
end

--- Selection to act on: explicit selection, else the visual range.
local function target_selection(buf, from_visual)
  if from_visual then
    local a, b = select_.visual_range(buf)
    vim.cmd("normal! \27")
    if not a then return nil end
    local sel = select_.make(buffer().record(buffer().get(buf)), a, b)
    select_.set(buf, sel)
    return sel
  end
  return select_.get(buf)
end

--- Quick add: three prompts (label, type, strand) over the selection.
function M.quick_add(from_visual)
  local buf, st, rec = current()
  if not buf then return end
  local sel = target_selection(buf, from_visual)
  if not sel then
    -- on an ORF with nothing selected: promote it
    local _, fid = buffer().cursor_info(vim.api.nvim_get_current_win())
    if fid and fid < 0 then
      local o = rec:feature_by_id(fid)
      if o and o.detected then return require("plasmid.detectui").accept(buf, fid) end
      return require("plasmid.orfs").promote(buf, fid)
    end
    return notify("select a range first (visual mode, vaf, or :DnaSelect)", vim.log.levels.WARN)
  end
  vim.ui.input({ prompt = ("Label for %d..%d: "):format(sel.start, sel.stop) }, function(label)
    if not label or label == "" then return end
    vim.ui.input({
      prompt = "Type: ", default = "misc_feature",
      completion = "customlist,v:lua.require'plasmid.author'.complete_types",
    }, function(ftype)
      if not ftype or ftype == "" then return end
      local default_strand = rec:is_protein() and "." or "+"
      vim.ui.input({ prompt = "Strand (+ / - / . for none): ", default = default_strand }, function(s)
        if not s then return end
        local strand = (s == "-" and -1) or ((s == "." or s == "0" or rec:is_protein()) and 0) or 1
        local f = M.make_feature(rec, { type = ftype, label = label, sel = sel, strand = strand })
        M.commit(buf, { op = "add", feature = f }, "add " .. label)
        select_.set(buf, nil)
        notify(("added %s %s %s"):format(ftype, label, f.location_text))
      end)
    end)
  end)
end

local editors = {}

--- Open the editor for a feature (fid) or a new one (fid nil, optional selection).
function M.open_editor(buf, fid, sel)
  local st = buffer().get(buf)
  local rec = buffer().record(st)
  local f = fid and rec:feature_by_id(fid)
  local lines
  if f then
    lines = M.editor_text(rec, f, ("Edit %s %s"):format(f.type, f:label()))
  else
    lines = M.editor_text(rec, nil, "New feature")
    if sel then lines[6] = "misc_feature " .. M.location_for(rec, sel, 1) end
  end
  vim.cmd("botright 14split")
  local ebuf = vim.api.nvim_create_buf(false, true)
  vim.api.nvim_win_set_buf(0, ebuf)
  vim.api.nvim_buf_set_name(ebuf, ("plasmid://feature/%d/%s/%d"):format(buf, fid or "new", ebuf))
  vim.api.nvim_buf_set_lines(ebuf, 0, -1, false, lines)
  vim.bo[ebuf].buftype = "acwrite"
  vim.bo[ebuf].bufhidden = "wipe"
  vim.bo[ebuf].filetype = "plasmidfeature"
  vim.bo[ebuf].modified = false
  editors[ebuf] = { buf = buf, fid = fid, rec = rec }
  vim.api.nvim_create_autocmd("BufWriteCmd", {
    buffer = ebuf,
    callback = function() M.apply_editor(ebuf) end,
  })
  vim.api.nvim_win_set_cursor(0, { 6, 0 })
  return ebuf
end

--- Apply the editor's content. Returns true on success.
function M.apply_editor(ebuf)
  local e = editors[ebuf]
  if not e then return false end
  local st = buffer().get(e.buf)
  if not st or buffer().record(st) ~= e.rec then
    notify("the record this editor belongs to is no longer open", vim.log.levels.ERROR)
    return false
  end
  local f, err = M.parse_editor(e.rec, vim.api.nvim_buf_get_lines(ebuf, 0, -1, false))
  if not f then
    notify(err, vim.log.levels.ERROR)
    return false
  end
  if e.fid then
    local old = e.rec:feature_by_id(e.fid)
    if not old then
      notify("feature was deleted meanwhile", vim.log.levels.ERROR)
      return false
    end
    f.color = f.color or old.color
    f.no_arrow = old.no_arrow
    f.pinned = old.pinned
    M.commit(e.buf, { op = "replace", id = e.fid, feature = f }, "edit " .. f:label())
  else
    f.pinned = true
    M.commit(e.buf, { op = "add", feature = f }, "add " .. f:label())
  end
  vim.bo[ebuf].modified = false
  editors[ebuf] = nil
  local win = vim.fn.bufwinid(ebuf)
  if win ~= -1 and #vim.api.nvim_list_wins() > 1 then vim.api.nvim_win_close(win, true) end
  notify(("%s %s %s"):format(e.fid and "updated" or "added", f.type, f.location_text))
  return true
end

local function feature_target(buf)
  local st = buffer().get(buf)
  local _, fid = buffer().cursor_info(vim.api.nvim_get_current_win())
  if not fid and st.selection and st.selection.fid then fid = st.selection.fid end
  if not fid then
    local f = select_.feature_under_cursor(vim.api.nvim_get_current_win())
    fid = f and f.id
  end
  return fid
end

local function not_overlay(fid)
  if fid and fid < 0 then
    notify("that is a proposal (ORF or detected feature), not a saved feature; ga accepts it", vim.log.levels.WARN)
    return false
  end
  return true
end

function M.edit_under_cursor()
  local buf = current()
  if not buf then return end
  local fid = feature_target(buf)
  if not fid then return notify("no feature under cursor", vim.log.levels.WARN) end
  if not not_overlay(fid) then return end
  M.open_editor(buf, fid)
end

function M.add_with_editor()
  local buf = current()
  if not buf then return end
  M.open_editor(buf, nil, select_.get(buf))
end

function M.delete(force)
  local buf, st, rec = current()
  if not buf then return end
  local fid = feature_target(buf)
  if not not_overlay(fid) then return end
  local f = fid and rec:feature_by_id(fid)
  if not f then return notify("no feature under cursor", vim.log.levels.WARN) end
  local function go()
    M.commit(buf, { op = "remove", id = f.id }, "delete " .. f:label())
    if st.selection and st.selection.fid == f.id then select_.set(buf, nil) end
    notify(("deleted %s %s (u to undo)"):format(f.type, f:label()))
  end
  if force then return go() end
  vim.ui.select({ "Delete " .. f.type .. " " .. f:label(), "Cancel" }, { prompt = "Delete feature?" },
    function(choice, idx) if idx == 1 then go() end end)
end

---------------------------------------------------------------------------
-- Personal library
---------------------------------------------------------------------------

function M.library_path()
  return config.options.library_path or (vim.fn.stdpath("data") .. "/plasmid/library.gb")
end

--- Library record for a feature: its residues 5'->3', annotated end to end.
function M.library_record(rec, f, source_name)
  local residues = select_.feature_residues(rec, f, "a")
  local label = f:label()
  local quals = {}
  for _, q in ipairs(f.qualifiers) do
    if q.key ~= "translation" and q.key ~= "codon_start" then quals[#quals + 1] = q end
  end
  local spliced = #f.parts > 1
  local name = label:gsub("[^%w%-%._]", "_"):sub(1, 16)
  local lib = record.Record.new({
    name = name, id = name,
    description = ("%s (%s) from %s"):format(label, f.type, source_name or rec.name),
    molecule_type = rec.molecule_type, alphabet = rec.alphabet, topology = "linear",
    seq = spliced and rec:extract(f) or residues, format = "genbank",
  })
  lib:add_feature(record.Feature.new({
    type = f.type, location_text = ("1..%d"):format(lib:len()), qualifiers = quals, color = f:file_color(),
  }))
  return lib
end

--- Append records to the library file (created on first use).
function M.library_append(path, recs)
  local pio = require("plasmid.io")
  local existing = {}
  local res = vim.fn.filereadable(path) == 1 and pio.read(path) or nil
  local seen = {}
  if res then
    existing = res.records
    for _, r in ipairs(existing) do seen[r.name .. "\0" .. r.seq:tostring():upper()] = true end
  end
  local added = 0
  for _, r in ipairs(recs) do
    local key = r.name .. "\0" .. r.seq:tostring():upper()
    if not seen[key] then
      existing[#existing + 1] = r
      seen[key] = true
      added = added + 1
    end
  end
  if added == 0 then return 0 end
  vim.fn.mkdir(vim.fn.fnamemodify(path, ":h"), "p")
  local ok, err = pio.write(path, existing, "genbank", res and res.file or nil, "genbank")
  if not ok then error(err) end
  return added
end

function M.save_to_library()
  local buf, st, rec = current()
  if not buf then return end
  local fid = feature_target(buf)
  local f = fid and rec:feature_by_id(fid)
  if f and f.overlay then
    -- saving an ORF to the library stores it as a CDS
    f = record.Feature.new({ type = "CDS", location_text = f.location_text, qualifiers = f.qualifiers })
  end
  local function save(feat)
    local n = M.library_append(M.library_path(), { M.library_record(rec, feat) })
    if n == 0 then
      notify(("%s is already in your library"):format(feat:label()))
    else
      notify(("saved %s to %s"):format(feat:label(), vim.fn.fnamemodify(M.library_path(), ":~")))
    end
  end
  if f then return save(f) end
  local sel = select_.get(buf)
  if not sel then return notify("no feature or selection to save", vim.log.levels.WARN) end
  vim.ui.input({ prompt = "Library label: " }, function(label)
    if not label or label == "" then return end
    vim.ui.input({ prompt = "Type: ", default = "misc_feature",
      completion = "customlist,v:lua.require'plasmid.author'.complete_types" }, function(t)
      if not t or t == "" then return end
      save(M.make_feature(rec, { type = t, label = label, sel = sel, strand = 1 }))
    end)
  end)
end

return M
