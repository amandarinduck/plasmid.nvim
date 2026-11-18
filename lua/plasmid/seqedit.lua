-- Sequence editing on a plasmid buffer: validation, the edit runner (apply,
-- impact report, cut policy, one undo step, view refresh) and the edit
-- primitives used by keymaps, commands and :DnaEdit.
local tx = require("plasmid.model.transaction")
local impact = require("plasmid.bio.impact")
local alphabet = require("plasmid.model.alphabet")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function config() return require("plasmid.config").options end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end

---------------------------------------------------------------------------
-- Residue validation
---------------------------------------------------------------------------

local NUC = "ACGTURYSWKMBDHVN"
local PROT = "ACDEFGHIKLMNPQRSTVWYBZXJUO*"

--- Clean typed or pasted residues for `rec`: whitespace, digits and "/" are
--- dropped (so GenBank ORIGIN text can be pasted), case follows the record
--- (lower-case records stay lower case; FASTA keeps what was typed), and T/U
--- follow the record. Returns residues or nil, error.
function M.clean(rec, text)
  local s = (text or ""):gsub("[%s%d/]", "")
  if s == "" then return nil, "no residues given" end
  local allowed = rec:is_protein() and PROT or NUC
  local bad = s:upper():gsub("[" .. allowed:gsub("%*", "%%*") .. "]", "")
  if bad ~= "" then
    local c = bad:sub(1, 1)
    return nil, ("'%s' is not a valid %s letter"):format(c, rec:is_protein() and "amino-acid" or "IUPAC nucleotide")
  end
  if not rec:is_protein() then
    local seq = rec.seq:tostring()
    local has_u, has_t = seq:find("[Uu]") ~= nil, seq:find("[Tt]") ~= nil
    if has_u and not has_t then s = s:gsub("T", "U"):gsub("t", "u") else s = s:gsub("U", "T"):gsub("u", "t") end
  end
  if rec.format ~= "fasta" then
    local seq = rec.seq:tostring()
    if seq ~= "" and not seq:find("%u") then s = s:lower() elseif seq:find("%u") and not seq:find("%l") then s = s:upper() end
  end
  return s
end

---------------------------------------------------------------------------
-- Runner
---------------------------------------------------------------------------

--- Protein for a CDS as GenBank writes /translation (no trailing stop), or nil
--- when the coding sequence now has an internal stop codon.
function M.translation_of(rec, f)
  local translate = require("plasmid.bio.translate")
  local res = rec:extract(f)
  local cs = tonumber(f:get("codon_start")) or 1
  local t = tonumber(f:get("transl_table")) or (rec:is_circular() and 11 or 1)
  local prot = translate.translate(res:sub(cs), { table = t, cds_start = cs == 1 })
  prot = prot:gsub("%*$", "")
  if prot:find("*", 1, true) then return nil end
  return prot
end

function M.is_seq_op(op)
  if not op then return false end
  if tx.SEQ_OPS[op.op] then return true end
  if op.op == "batch" then
    for _, o in ipairs(op.ops) do if M.is_seq_op(o) then return true end end
  end
  return false
end

-- Changes of an applied op (batch: merged per feature, first old -> last new).
local function collect(op)
  if op.op ~= "batch" then return op.changes or {}, op.moved or 0 end
  local by_id, order, moved = {}, {}, 0
  for _, o in ipairs(op.ops) do
    local ch, mv = collect(o)
    moved = moved + mv
    for _, c in ipairs(ch) do
      local id = (c.old or c.new).id
      local cur = by_id[id]
      if not cur then
        by_id[id] = { old = c.old, new = c.new, old_res = c.old_res, new_res = c.new_res, inverted = c.inverted }
        order[#order + 1] = id
      else
        cur.new, cur.new_res = c.new, c.new_res
        cur.inverted = cur.inverted or c.inverted
      end
    end
  end
  local out = {}
  for _, id in ipairs(order) do
    local c = by_id[id]
    if c.new == nil or c.new_res ~= c.old_res then out[#out + 1] = c end
  end
  return out, moved
end

--- After any change of the residues (edit, undo, redo): drop stale overlays
--- and a selection that no longer fits, refresh ORFs.
function M.after_change(buf)
  local st = buffer().get(buf)
  if not st then return end
  local rec = buffer().record(st)
  if st.detect then
    st.detect = nil -- proposals refer to the old coordinates; :DnaDetect again
    st.detect_cleared_by_edit = true
  end
  local sel = st.selection
  if sel and (sel.start > rec:len() or sel.stop > rec:len()) then st.selection = nil end
  require("plasmid.searchui").refresh(buf)
  require("plasmid.enzymes").refresh(buf)
  require("plasmid.orfs").refresh(buf)
end

--- Apply a sequence op as one undoable change.
--- opts: { label, force = bool (override the refuse policy), cursor = bp,
---         select = { start, stop } | false, quiet = bool }
--- Returns the impact report, or nil, reason.
function M.run(buf, op, opts)
  opts = opts or {}
  local st = buffer().get(buf)
  if not st then return nil, "not a plasmid buffer" end
  local rec = buffer().record(st)
  if rec.no_sequence and rec:len() == 0 then
    local msg = "this record has no residues (CON/WGS); :DnaFetch " .. (rec.id or rec.name) .. " opens it with its sequence"
    require("plasmid.util").notify("plasmid: " .. msg, vim.log.levels.WARN)
    return nil, msg
  end
  local seq_rev0, rev0 = rec.seq_revision, rec.revision
  local inverse = tx.apply(rec, op)
  local changes, moved = collect(op)
  local report = impact.build(rec, { changes = changes, moved = moved })
  local policy = (config().edit or {}).cut_policy or "truncate"
  if report.cut and policy == "refuse" and not opts.force then
    tx.apply(rec, inverse)
    -- nothing happened: a :DnaEdit buffer opened before must not look stale
    rec.seq_revision, rec.revision = seq_rev0, rev0
    local names, seen = {}, {}
    for _, e in ipairs(report.entries) do
      if (e.kind == "truncated" or e.kind == "destroyed") and not seen[e.label] then
        seen[e.label] = true
        names[#names + 1] = e.label
      end
    end
    local msg = ("edit refused: it would %s %s (cut_policy = \"refuse\"; add ! to force)"):format(
      #names == 1 and "cut" or "cut features", table.concat(names, ", "))
    notify(msg, vim.log.levels.WARN)
    return nil, msg
  end
  -- follow-up changes to affected features, in the same undo step:
  -- a /note on truncated features (cut_policy = "flag") and fresh /translation
  -- values for coding features whose bases changed
  local forward = op
  local update_tr = (config().edit or {}).update_translations ~= false
  local extra_ops, extra_inv = {}, {}
  for _, e in ipairs(report.entries) do
    local f = e.kind ~= "destroyed" and rec:feature_by_id(e.id) or nil
    if f then
      local quals, changed = {}, false
      for i, q in ipairs(f.qualifiers) do quals[i] = { key = q.key, value = q.value, quoted = q.quoted } end
      if e.kind == "truncated" and report.cut and policy == "flag" then
        quals[#quals + 1] = { key = "note", quoted = true,
          value = ("truncated by a sequence edit: %d of %d bases removed"):format(e.old_len - e.new_len, e.old_len) }
        changed = true
      end
      if update_tr and f.type == "CDS" and f:get("translation") and not f:get("transl_except") and not f:get("pseudo") then
        local prot = M.translation_of(rec, f)
        local keep = {}
        for _, q in ipairs(quals) do
          if q.key ~= "translation" then keep[#keep + 1] = q
          elseif prot then
            if q.value ~= prot then q.value = prot; changed = true; e.translation = "updated" end
            keep[#keep + 1] = q
          else
            changed = true
            e.translation = "removed"
          end
        end
        quals = keep
      end
      if changed then
        local g = require("plasmid.model.record").Feature.new({
          type = f.type, location_text = f.location_text, qualifiers = quals,
          color = f.color, no_arrow = f.no_arrow, origin = f.origin, recent_id = f.recent_id,
          primer_seq = f.primer_seq,
        })
        g.pinned = f.pinned
        local rop = { op = "replace", id = f.id, feature = g }
        extra_ops[#extra_ops + 1] = rop
        table.insert(extra_inv, 1, tx.apply(rec, rop))
      end
    end
  end
  if #extra_ops > 0 then
    forward = { op = "batch", ops = { op, unpack(extra_ops) } }
    extra_inv[#extra_inv + 1] = inverse
    inverse = { op = "batch", ops = extra_inv }
  end
  tx.record(rec, forward, inverse, opts.label or op.label or op.op)
  st.last_report = report
  st.last_report_label = opts.label or op.label or op.op
  if opts.select ~= nil then
    local s = require("plasmid.select")
    s.set(buf, opts.select and s.make(rec, opts.select[1], opts.select[2]) or nil)
  end
  M.after_change(buf)
  buffer().refresh_modified(buf)
  buffer().render(buf, { keep_bp = opts.cursor })
  if not opts.quiet then
    local unit = rec:unit()
    local text = ("%s: %s"):format(st.last_report_label, impact.summary(report, unit))
    notify(text, impact.alarming(report) and vim.log.levels.WARN or vim.log.levels.INFO)
  end
  return report
end

---------------------------------------------------------------------------
-- Primitives (positions 1-based; s > e on a circular record wraps the origin)
---------------------------------------------------------------------------

local function rec_of(buf)
  local st = buffer().get(buf)
  return st and buffer().record(st)
end

--- Insert residues before base `pos` (1..L+1).
function M.insert(buf, pos, text, opts)
  local rec = rec_of(buf)
  local s, err = M.clean(rec, text)
  if not s then notify(err, vim.log.levels.ERROR) return nil, err end
  opts = opts or {}
  return M.run(buf, { op = "splice", start = pos, stop = pos - 1, ins = s }, {
    label = opts.label or ("insert %d %s at %d"):format(#s, rec:unit(), pos),
    force = opts.force, cursor = opts.cursor or (pos + #s - 1), select = opts.select,
  })
end

local function wrap_ops(rec, s, e, ins)
  local L = rec:len()
  if s <= e then return { op = "splice", start = s, stop = e, ins = ins } end
  if not rec:is_circular() then s, e = e, s; return { op = "splice", start = s, stop = e, ins = ins } end
  -- s..L then 1..e: the replacement goes where s..L was
  return { op = "batch", ops = {
    { op = "splice", start = s, stop = L, ins = ins },
    { op = "splice", start = 1, stop = e, ins = "" },
  } }
end

local function span_len(rec, s, e)
  if s <= e then return e - s + 1 end
  return rec:len() - s + 1 + e
end

--- Delete s..e.
function M.delete(buf, s, e, opts)
  local rec = rec_of(buf)
  opts = opts or {}
  local n = span_len(rec, s, e)
  if n >= rec:len() then
    local msg = "refusing to delete the whole sequence"
    notify(msg, vim.log.levels.WARN)
    return nil, msg
  end
  local after = (s <= e) and s or 1
  return M.run(buf, wrap_ops(rec, s, e, ""), {
    label = opts.label or ("delete %d..%d (%d %s)"):format(s, e, n, rec:unit()),
    force = opts.force, cursor = opts.cursor or math.min(after, rec:len() - n), select = false,
  })
end

--- Replace s..e by residues.
function M.replace(buf, s, e, text, opts)
  local rec = rec_of(buf)
  local ins, err = M.clean(rec, text)
  if not ins then notify(err, vim.log.levels.ERROR) return nil, err end
  opts = opts or {}
  local n = span_len(rec, s, e)
  return M.run(buf, wrap_ops(rec, s, e, ins), {
    label = opts.label or ("replace %d..%d (%d → %d %s)"):format(s, e, n, #ins, rec:unit()),
    force = opts.force, cursor = opts.cursor or s,
    select = opts.select == nil and { s, s + #ins - 1 } or opts.select,
  })
end

--- Reverse-complement s..e (s <= e) or, with no range, the whole record.
function M.revcomp(buf, s, e, opts)
  local rec = rec_of(buf)
  opts = opts or {}
  if rec:is_protein() then return nil, notify("protein records have no reverse complement", vim.log.levels.WARN) end
  local whole = not s
  s, e = s or 1, e or rec:len()
  if s > e then
    local msg = "the region crosses the origin; rotate the origin first (:DnaRotate)"
    notify(msg, vim.log.levels.WARN)
    return nil, msg
  end
  return M.run(buf, { op = "revcomp", start = s, stop = e }, {
    label = opts.label or (whole and "reverse-complement the record" or ("reverse-complement %d..%d"):format(s, e)),
    force = opts.force, cursor = opts.cursor or (whole and 1 or s),
    select = (not whole) and { s, e } or nil,
  })
end

--- Make base `o` the new base 1 (circular records).
function M.rotate(buf, o, opts)
  local rec = rec_of(buf)
  opts = opts or {}
  if not rec:is_circular() then
    local msg = "only circular records have an origin to move (:DnaTopology circular)"
    notify(msg, vim.log.levels.WARN)
    return nil, msg
  end
  if o == 1 then return nil, notify("base 1 is already the origin") end
  return M.run(buf, { op = "rotate", origin = o }, {
    label = opts.label or ("move origin to %d"):format(o), cursor = 1, select = false,
  })
end

---------------------------------------------------------------------------
-- Keys (linear or circular view; the cursor's base pair is the target)
---------------------------------------------------------------------------

local function here()
  local win = vim.api.nvim_get_current_win()
  local buf = vim.api.nvim_win_get_buf(win)
  local st = buffer().get(buf)
  if not st then return nil end
  local bp = buffer().cursor_info(win)
  if not bp then notify("put the cursor on the sequence", vim.log.levels.WARN) return nil end
  return buf, st, buffer().record(st), bp
end

-- Deleted residues go to the unnamed register (and "1), as with Vim's x/d.
local function stash(rec, s, e)
  local text = s <= e and rec.seq:sub(s, e) or (rec.seq:sub(s, rec:len()) .. rec.seq:sub(1, e))
  text = require("plasmid.select").display_case(rec, text)
  vim.fn.setreg('"', text, "c")
  vim.fn.setreg("1", text, "c")
end

local function range_from(rec, bp, count)
  local L = rec:len()
  local e = bp + count - 1
  if e > L then e = rec:is_circular() and (e - L) or L end
  return bp, e
end

--- i / a: insert residues before / after the base under the cursor.
function M.key_insert(after)
  local buf, st, rec, bp = here()
  if not buf then return end
  local pos = after and bp + 1 or bp
  vim.ui.input({ prompt = ("Insert %s %d: "):format(after and "after" or "before", bp) }, function(text)
    if text and text ~= "" then M.insert(buf, pos, text) end
  end)
end

--- x / [count]x: delete bases from the cursor.
function M.key_delete(count)
  local buf, st, rec, bp = here()
  if not buf then return end
  local s, e = range_from(rec, bp, math.max(1, count or 1))
  stash(rec, s, e)
  M.delete(buf, s, e)
end

--- r / [count]r{char}: replace bases under the cursor with a typed letter.
function M.key_replace(count)
  local buf, st, rec, bp = here()
  if not buf then return end
  local ok, ch = pcall(vim.fn.getcharstr)
  if not ok or ch == "\27" or ch == "" then return end
  local n = math.max(1, count or 1)
  local s, e = range_from(rec, bp, n)
  M.replace(buf, s, e, ch:rep(n), { select = false, cursor = bp })
end

--- p / P: insert register contents after / before the cursor base.
function M.key_paste(after)
  local reg = vim.v.register
  local buf, st, rec, bp = here()
  if not buf then return end
  local text = vim.fn.getreg(reg)
  if text == "" then return notify(("register %s is empty"):format(reg), vim.log.levels.WARN) end
  local pos = after and bp + 1 or bp
  local clean = M.clean(rec, text)
  M.insert(buf, pos, text, clean and { select = { pos, pos + #clean - 1 }, cursor = pos + #clean - 1 } or nil)
end

local function visual()
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer().get(buf)
  if not st then return nil end
  local a, b = require("plasmid.select").visual_range(buf)
  vim.cmd("normal! \27")
  if not a then notify("selection is not on the sequence", vim.log.levels.WARN) return nil end
  return buf, st, buffer().record(st), a, b
end

--- visual d / x: delete the selected bases.
function M.visual_delete()
  local buf, st, rec, a, b = visual()
  if not buf then return end
  stash(rec, a, b)
  M.delete(buf, a, b)
end

--- visual c: replace the selected bases with typed residues.
function M.visual_change()
  local buf, st, rec, a, b = visual()
  if not buf then return end
  vim.ui.input({ prompt = ("Replace %d..%d with: "):format(a, b) }, function(text)
    if text and text ~= "" then M.replace(buf, a, b, text) end
  end)
end

--- visual p: replace the selected bases with register contents.
function M.visual_paste()
  local reg = vim.v.register
  local buf, st, rec, a, b = visual()
  if not buf then return end
  local text = vim.fn.getreg(reg)
  if text == "" then return notify(("register %s is empty"):format(reg), vim.log.levels.WARN) end
  M.replace(buf, a, b, text)
end

--- visual R: reverse-complement the selected bases in place.
function M.visual_revcomp()
  local buf, st, rec, a, b = visual()
  if not buf then return end
  M.revcomp(buf, a, b)
end

--- daf / dif / caf from the feature text objects (scheduled: called under textlock).
function M.feature_op(buf, f, kind, op)
  local rec = buffer().record(buffer().get(buf))
  local sel = require("plasmid.select").for_feature(rec, f)
  if op == "d" then
    if kind == "i" and #f.parts > 1 then
      -- each part, right to left, as one change
      local parts = {}
      for _, p in ipairs(f.parts) do if not p.ref and p.stop >= p.start then parts[#parts + 1] = p end end
      table.sort(parts, function(x, y) return x.start > y.start end)
      local ops = {}
      for _, p in ipairs(parts) do ops[#ops + 1] = { op = "splice", start = p.start, stop = p.stop, ins = "" } end
      local text = require("plasmid.select").feature_residues(rec, f, "i")
      vim.fn.setreg('"', require("plasmid.select").display_case(rec, text), "c")
      return M.run(buf, { op = "batch", ops = ops }, {
        label = ("delete %s (%d parts)"):format(f:label(), #parts), cursor = parts[#parts].start, select = false,
      })
    end
    stash(rec, sel.start, sel.stop)
    return M.delete(buf, sel.start, sel.stop, { label = ("delete %s"):format(f:label()) })
  elseif op == "c" then
    if kind == "i" and #f.parts > 1 then
      return notify("change the parts one at a time, or use caf for the whole span", vim.log.levels.WARN)
    end
    vim.ui.input({ prompt = ("Replace %s (%d..%d) with: "):format(f:label(), sel.start, sel.stop) }, function(text)
      if text and text ~= "" then M.replace(buf, sel.start, sel.stop, text, { label = ("change %s"):format(f:label()) }) end
    end)
  end
end

---------------------------------------------------------------------------
-- Report
---------------------------------------------------------------------------

--- Lines describing the last edit's impact.
function M.report_lines(buf)
  local st = buffer().get(buf)
  local rep = st and st.last_report
  if not rep then return nil end
  local rec = buffer().record(st)
  local lines = { st.last_report_label or "last edit", "" }
  if #rep.entries == 0 then lines[#lines + 1] = "  " .. impact.summary(rep, rec:unit()) end
  for _, e in ipairs(rep.entries) do lines[#lines + 1] = "  " .. impact.describe_entry(e, rec:unit()) end
  if rep.moved > 0 then
    lines[#lines + 1] = ""
    lines[#lines + 1] = ("  %d feature%s moved"):format(rep.moved, rep.moved == 1 and "" or "s")
  end
  return lines
end

function M.show_report(buf)
  local lines = M.report_lines(buf)
  if not lines then return notify("no sequence edit yet in this buffer") end
  local fbuf = vim.api.nvim_create_buf(false, true)
  vim.api.nvim_buf_set_lines(fbuf, 0, -1, false, lines)
  vim.bo[fbuf].modifiable = false
  vim.bo[fbuf].bufhidden = "wipe"
  local width = 20
  for _, l in ipairs(lines) do width = math.max(width, vim.fn.strdisplaywidth(l) + 2) end
  width = math.min(width, vim.o.columns - 4)
  local win = vim.api.nvim_open_win(fbuf, true, {
    relative = "editor", row = 2, col = math.max(0, math.floor((vim.o.columns - width) / 2)),
    width = width, height = math.min(#lines, vim.o.lines - 6), style = "minimal", border = "rounded",
    title = " edit impact ", title_pos = "center",
  })
  vim.keymap.set("n", "q", "<cmd>close<cr>", { buffer = fbuf, nowait = true })
  vim.keymap.set("n", "<Esc>", "<cmd>close<cr>", { buffer = fbuf, nowait = true })
  return fbuf, win
end

return M
