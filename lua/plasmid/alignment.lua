-- Sequence alignment: :DnaAlign, :DnaAlignReads, :DnaCompare, the alignment
-- view and the difference list.
--
-- The biology lives in plasmid.bio.align (the pairwise engine) and
-- plasmid.bio.variant (what a difference does to the feature it lands in).
-- This module only decides what is being aligned against what, draws the
-- alignment in blocks, lists the differences, and writes them back onto the
-- record as features. The difference list is the point of the exercise: a
-- sequencing read is read to find out what changed and whether it matters, so
-- every difference is reported with the feature it hits and, in a CDS, the
-- codon and the amino acid it becomes.
local A = require("plasmid.bio.align")
local V = require("plasmid.bio.variant")
local alphabet = require("plasmid.model.alphabet")

local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end
local function cfg() return require("plasmid.config").options.align or {} end

--- view buffer -> { src, rec, aln, diffs, ctx, targets, list_line, marks }
M.state = {}
--- read-table buffer -> { src, rec, rows, unplaced, first }
M.read_state = {}
M.last = nil         -- last alignment view
M.last_reads = nil   -- last read table

local ns = vim.api.nvim_create_namespace("plasmid_align")

M.MODES = { global = true, ["local"] = true, semi = true }
M.STRANDS = { auto = true, plus = true, minus = true }

---------------------------------------------------------------------------
-- Arguments
---------------------------------------------------------------------------

--- key=value options understood by the alignment commands.
local OPTION_KEYS = {
  mode = true, strand = true, name = true, min_identity = true,
  min_coverage = true, width = true,
}

--- Split `fargs` into key=value options and the rest. Returns opts, rest or
--- nil plus the message for an option that is not understood.
function M.parse_args(fargs)
  local o, rest = {}, {}
  for _, a in ipairs(fargs or {}) do
    local k, v = a:match("^(%a[%w_]*)=(.*)$")
    if k and OPTION_KEYS[k] then
      o[k] = v
    elseif k then
      return nil, nil, ("no such option: %s= (mode, strand, name, min_identity, min_coverage, width)"):format(k)
    else
      rest[#rest + 1] = a
    end
  end
  local c = cfg()
  o.mode = o.mode or c.mode or "semi"
  if not M.MODES[o.mode] then
    return nil, nil, ("mode must be semi, global or local, not %q"):format(o.mode)
  end
  o.strand = o.strand or "auto"
  if not M.STRANDS[o.strand] then
    return nil, nil, ("strand must be auto, plus or minus, not %q"):format(o.strand)
  end
  for _, k in ipairs({ "min_identity", "min_coverage", "width" }) do
    if o[k] then
      local n = tonumber(o[k])
      if not n then return nil, nil, ("%s= wants a number, not %q"):format(k, o[k]) end
      o[k] = n
    end
  end
  return o, rest
end

local NUC = "^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$"
local AA = "^[A-Za-z%*]+$"

--- A query argument: `selection`, a staged cloning part, an open buffer, a
--- file, a digest fragment, or residues typed inline. Everything but the
--- selection and the inline residues is resolved by plasmid.clone, which
--- already knows how to find a part.
--- Returns { name, seq, rec }, or nil and why not.
function M.resolve_query(buf, arg, rec)
  if arg == nil or arg == "" then return nil, "no query given" end
  if arg:lower() == "selection" then
    local st = buffer().get(buf)
    local sel = st and require("plasmid.select").get(buf)
    if not sel then
      return nil, "nothing is selected in this buffer (vaf, visual mode, or :DnaSelect)"
    end
    local srec = buffer().record(st)
    local S = require("plasmid.select")
    local seq = S.top_strand(srec, sel)
    if seq == "" then return nil, "the selection is empty" end
    return { name = ("%s_%d-%d"):format(srec.name, sel.start, sel.stop), seq = seq }
  end
  local p, err = require("plasmid.clone").resolve(buf, arg)
  if p then return { name = p.name, seq = p.rec.seq:tostring(), rec = p.rec } end
  -- residues typed on the command line
  local s = arg:gsub("%s+", "")
  local protein = rec and rec:is_protein()
  if #s >= 8 and ((protein and s:match(AA)) or (not protein and s:match(NUC))) then
    return { name = "query", seq = s, inline = true }
  end
  return nil, err
end

---------------------------------------------------------------------------
-- Running an alignment
---------------------------------------------------------------------------

local function score_opts(extra)
  local c = cfg()
  return vim.tbl_extend("force", {
    match = c.match or 2,
    mismatch = c.mismatch or -3,
    gap_open = c.gap_open or -5,
    gap_extend = c.gap_extend or -2,
  }, extra or {})
end

--- Place a query on one strand only. plasmid.bio.align.place always tries
--- both strands, so strand=plus / strand=minus repeat its doubling of a
--- circular reference here (see the report: an engine `strand` option would
--- make this unnecessary).
local function place_one(ref, query, opts, circular)
  if not circular then
    local aln, why = A.align(ref, query, opts)
    if not aln then return nil, why end
    aln.wraps = false
    return aln
  end
  local L = #ref
  local extra = math.min(L, #query + (opts.band or A.defaults.band))
  local aln, why = A.align(ref .. ref:sub(1, extra), query, opts)
  if not aln then return nil, why end
  if aln.ref_span > L then return nil, "the query is longer than the circular reference" end
  aln.wraps = aln.ref_stop > L
  aln.ref_start_abs, aln.ref_stop_abs = aln.ref_start, aln.ref_stop
  aln.ref_start, aln.ref_stop = (aln.ref_start - 1) % L + 1, (aln.ref_stop - 1) % L + 1
  aln.ref_length = L
  return aln
end

--- Align `seq` against the record. `o` is the parsed options table.
--- Returns the alignment (with query_seq, the residues as they were aligned)
--- or nil and why not.
function M.align_query(rec, seq, o)
  local ref = rec.seq:tostring():upper()
  local query = seq:upper():gsub("[^%a%*]", "")
  if query == "" then return nil, "the query has no residues" end
  local protein = rec:is_protein()
  local opts = score_opts({ mode = o.mode, protein = protein or nil })
  local circular = rec:is_circular() and not protein and o.mode ~= "global"
  local aln, why
  if protein or o.strand == "plus" then
    aln, why = place_one(ref, query, opts, circular)
    if aln then aln.strand = 1 end
  elseif o.strand == "minus" then
    local rc = alphabet.revcomp(query)
    aln, why = place_one(ref, rc, opts, circular)
    if aln then aln.strand, aln.query_rc = -1, rc end
  elseif circular then
    aln, why = A.place(ref, query, vim.tbl_extend("force", opts, { circular = true }))
  else
    aln, why = A.align_strands(ref, query, opts)
  end
  if not aln then return nil, why or "no alignment" end
  aln.query_seq = aln.strand == -1 and (aln.query_rc or alphabet.revcomp(query)) or query
  aln.query_input = query
  aln.circular = circular
  return aln
end

--- Differences of an alignment produced by M.align_query, in the record's own
--- coordinates (an alignment that crosses the origin is walked on the doubled
--- reference and wrapped back).
function M.differences(rec, aln)
  local ref = rec.seq:tostring():upper()
  if aln.wraps and aln.ref_start_abs then
    local L = #ref
    local ds = A.differences(ref .. ref, aln.query_seq,
      { ref_start = aln.ref_start_abs, query_start = aln.query_start, ops = aln.ops })
    for _, d in ipairs(ds) do d.pos = (d.pos - 1) % L + 1 end
    return ds
  end
  return A.differences(ref, aln.query_seq, aln)
end

---------------------------------------------------------------------------
-- The view
---------------------------------------------------------------------------

local function pad(s, w) return s .. (" "):rep(math.max(0, w - vim.fn.strdisplaywidth(s))) end
local function rpad(s, w) return (" "):rep(math.max(0, w - vim.fn.strdisplaywidth(s))) .. s end
local function pct(x) return ("%.1f%%"):format(100 * (x or 0)) end
local function clip(s, n)
  if vim.fn.strdisplaywidth(s) <= n then return s end
  return vim.fn.strcharpart(s, 0, n - 1) .. "…"
end

--- A difference's key in the column walk: substitutions and deletions are
--- keyed by every reference base they cover, insertions by the base they
--- follow.
local function diff_keys(diffs)
  local key = {}
  for idx, d in ipairs(diffs) do
    if d.kind == "ins" then
      key["i" .. d.pos] = idx
    else
      for k = 0, math.max(d.len, 1) - 1 do
        key[(d.kind == "del" and "d" or "s") .. (d.pos + k)] = idx
      end
    end
  end
  return key
end

--- One table per alignment column: the two residues, the match character, the
--- coordinates they carry and the difference they belong to.
function M.columns(ctx)
  local aln = ctx.aln
  local r, q, mid = aln.ref_aln, aln.query_aln, aln.mid
  local key = diff_keys(ctx.diffs)
  local cols = {}
  local i, j = ctx.ref_from, ctx.query_from
  local prev_ref = ctx.ref_map(i - 1)
  for c = 1, #r do
    local rc, qc = r:sub(c, c), q:sub(c, c)
    local e = { ref = rc, query = qc, mid = mid:sub(c, c) }
    if rc ~= "-" then
      e.ref_pos = ctx.ref_map(i)
      i = i + 1
    end
    if qc ~= "-" then
      e.query_pos = ctx.query_map(j)
      j = j + 1
    end
    if rc == "-" then
      e.diff, e.gap = key["i" .. prev_ref], true
    elseif qc == "-" then
      e.diff, e.gap = key["d" .. e.ref_pos], true
    elseif e.mid ~= "|" then
      e.diff = key["s" .. e.ref_pos]
    end
    if e.ref_pos then prev_ref = e.ref_pos end
    cols[c] = e
  end
  return cols
end

--- Column ranges the view draws. The whole alignment when it is short enough,
--- otherwise a window of `align.context` residues around each difference.
local function regions(cols, diffs_present)
  local c = cfg()
  local limit = c.full_limit or 20000
  if #cols <= limit then return { { 1, #cols } }, false end
  local ctx_n = c.context or 30
  local out = {}
  for i, e in ipairs(cols) do
    if e.diff then
      local a, b = math.max(1, i - ctx_n), math.min(#cols, i + ctx_n)
      local last = out[#out]
      if last and a <= last[2] + 1 then last[2] = math.max(last[2], b)
      else out[#out + 1] = { a, b } end
    end
  end
  if #out == 0 then out[1] = { 1, math.min(#cols, 2 * ctx_n) } end
  return out, true
end

--- The alignment drawn in blocks: reference row with coordinates, match row,
--- query row. Returns lines, highlights and the navigation marks
--- ({ line, col, len, diff }).
function M.block_lines(ctx, width, first_line)
  local cols = M.columns(ctx)
  local rl, ql = clip(ctx.ref_label, 14), clip(ctx.query_label, 14)
  local lw = math.max(vim.fn.strdisplaywidth(rl), vim.fn.strdisplaywidth(ql))
  local nw = math.max(#tostring(ctx.ref_len), #tostring(ctx.query_len))
  local prefix = 1 + lw + 1 + nw + 1
  local fixed = cfg().width
  local per = (fixed and fixed > 0) and fixed or math.max(20, (width or 80) - prefix - nw - 1)
  local lines, hl, marks = {}, {}, {}
  local n = first_line or 0
  local function put(s) lines[#lines + 1] = s; n = n + 1; return n end

  local list, trimmed = regions(cols, #ctx.diffs > 0)
  for ri, reg in ipairs(list) do
    if ri > 1 then put("  ⋮") end
    for a = reg[1], reg[2], per do
      local b = math.min(a + per - 1, reg[2])
      local rseq, qseq, midseq = {}, {}, {}
      local rfrom, rto, qfrom, qto
      for c = a, b do
        local e = cols[c]
        rseq[#rseq + 1], qseq[#qseq + 1], midseq[#midseq + 1] = e.ref, e.query, e.mid
        if e.ref_pos then rfrom = rfrom or e.ref_pos; rto = e.ref_pos end
        if e.query_pos then qfrom = qfrom or e.query_pos; qto = e.query_pos end
      end
      rseq, qseq, midseq = table.concat(rseq), table.concat(qseq), table.concat(midseq)
      local rline = (" %s %s %s %s"):format(pad(rl, lw), rpad(tostring(rfrom or "-"), nw), rseq, tostring(rto or "-"))
      local mline = (" %s %s %s"):format((" "):rep(lw), (" "):rep(nw), midseq)
      local qline = (" %s %s %s %s"):format(pad(ql, lw), rpad(tostring(qfrom or "-"), nw), qseq, tostring(qto or "-"))
      local rrow, mrow, qrow = put(rline), put(mline), put(qline)
      hl[#hl + 1] = { rrow - 1, 1 + lw, prefix - 1, "PlasmidAlignCoord" }
      hl[#hl + 1] = { rrow - 1, prefix + #rseq, #rline, "PlasmidAlignCoord" }
      hl[#hl + 1] = { qrow - 1, 1 + lw, prefix - 1, "PlasmidAlignCoord" }
      hl[#hl + 1] = { qrow - 1, prefix + #qseq, #qline, "PlasmidAlignCoord" }
      hl[#hl + 1] = { mrow - 1, prefix, prefix + #midseq, "PlasmidAlignMatch" }
      -- runs of differing columns, highlighted on both residue rows
      local k = a
      while k <= b do
        local e = cols[k]
        if e.diff then
          local j = k
          while j + 1 <= b and cols[j + 1].diff == e.diff and cols[j + 1].gap == e.gap do j = j + 1 end
          local c0, c1 = prefix + (k - a), prefix + (j - a) + 1
          local group = e.gap and "PlasmidAlignGap" or "PlasmidAlignMismatch"
          hl[#hl + 1] = { rrow - 1, c0, c1, group }
          hl[#hl + 1] = { qrow - 1, c0, c1, group }
          marks[#marks + 1] = { line = rrow, col = c0, len = c1 - c0, diff = e.diff }
          k = j + 1
        else
          k = k + 1
        end
      end
      put("")
    end
  end
  if trimmed then
    put(("  (%d columns; only the %d region%s around differences are shown — raise align.full_limit for the whole alignment)")
      :format(#cols, #list, #list == 1 and " is" or "s"))
  end
  return lines, hl, marks
end

--- What each difference does to the feature it lands in, as a small table:
--- feature, type, how many differences, and the verdict on the whole feature.
function M.feature_lines(entries, unit, first_line)
  local lines, hl = {}, {}
  local n = first_line or 0
  if #entries == 0 then
    lines[1] = "  (no feature is touched)"
    return lines, hl
  end
  local header = { "feature", "type", "diffs", "effect" }
  local body, alarming = {}, {}
  for i, e in ipairs(entries) do
    local parts = { e.kind }
    if e.premature_stop then parts[#parts + 1] = ("stop codon at aa %d"):format(e.premature_stop) end
    if e.lost_start then parts[#parts + 1] = "start codon lost" end
    if e.lost_stop then parts[#parts + 1] = "stop codon lost" end
    if e.aa_changes and e.aa_changes > 0 and not e.frameshift then
      parts[#parts + 1] = ("%d aa changed"):format(e.aa_changes)
    end
    if e.old_res and e.new_res and #e.old_res ~= #e.new_res then
      parts[#parts + 1] = ("%d → %d %s"):format(#e.old_res, #e.new_res, unit)
    end
    body[i] = { clip(e.label, 28), e.type, tostring(#e.diffs), table.concat(parts, ", ") }
    alarming[i] = e.frameshift or e.premature_stop or e.lost_start or e.lost_stop
  end
  local w = {}
  for c = 1, #header do
    w[c] = vim.fn.strdisplaywidth(header[c])
    for _, r in ipairs(body) do w[c] = math.max(w[c], vim.fn.strdisplaywidth(r[c])) end
  end
  local function fmt(r)
    local out = {}
    for c = 1, #header do
      out[c] = (c == 3) and rpad(r[c], w[c]) or (c < #header and pad(r[c], w[c]) or r[c])
    end
    return "  " .. table.concat(out, "  ")
  end
  lines[1] = fmt(header)
  hl[#hl + 1] = { n, 0, #lines[1], "PlasmidHeader" }
  for i, r in ipairs(body) do
    lines[#lines + 1] = fmt(r)
    if alarming[i] then hl[#hl + 1] = { n + #lines - 1, 0, #lines[#lines], "PlasmidAlignNonsense" } end
  end
  return lines, hl
end

--- "8 substitutions, 1 insertion, 2 deletions · 3 silent, 4 missense ..."
function M.summary(diffs, entries)
  local sub, ins, del = 0, 0, 0
  local effects = { silent = 0, missense = 0, nonsense = 0 }
  for _, d in ipairs(diffs) do
    if d.kind == "sub" then sub = sub + d.len
    elseif d.kind == "ins" then ins = ins + 1
    else del = del + 1 end
    if d.effect then effects[d.effect] = (effects[d.effect] or 0) + 1 end
  end
  local fs, inframe = 0, 0
  for _, e in ipairs(entries or {}) do
    if e.kind == "frameshift" then fs = fs + 1
    elseif e.kind == "in-frame indel" then inframe = inframe + 1 end
  end
  local function plural(n, one) return ("%d %s%s"):format(n, one, n == 1 and "" or "s") end
  local parts = { ("%s, %s, %s"):format(plural(sub, "substitution"), plural(ins, "insertion"), plural(del, "deletion")) }
  if effects.silent + effects.missense + effects.nonsense > 0 then
    parts[#parts + 1] = ("%d silent, %d missense, %d nonsense"):format(effects.silent, effects.missense, effects.nonsense)
  end
  if fs + inframe > 0 then
    parts[#parts + 1] = ("%s, %d in-frame"):format(plural(fs, "frameshifted feature"), inframe)
  end
  return table.concat(parts, " · "), {
    sub = sub, ins = ins, del = del, silent = effects.silent, missense = effects.missense,
    nonsense = effects.nonsense, frameshift = fs, in_frame = inframe,
  }
end

--- The difference list: one line per difference (V.describe plus the feature
--- verdict), then the summary. Returns lines, highlights and the difference
--- index of each line.
function M.difference_lines(ctx, first_line)
  local lines, hl, index = {}, {}, {}
  local n = first_line or 0
  local entries = V.by_feature(ctx.rec, ctx.diffs)
  ctx.entries = entries
  local function put(s) lines[#lines + 1] = s; n = n + 1; return n end
  if #ctx.diffs == 0 then
    put("  (none: the query matches the record over the whole alignment)")
    return lines, hl, index, nil
  end
  local w = #tostring(#ctx.diffs)
  for i, d in ipairs(ctx.diffs) do
    local row = put(("  %s  %s"):format(rpad(tostring(i), w), V.describe(d, ctx.unit)))
    index[row] = i
    local group = d.effect == "nonsense" and "PlasmidAlignNonsense"
      or (d.effect == "missense" or d.kind ~= "sub") and "PlasmidAlignMissense"
      or d.effect == "silent" and "PlasmidAlignSilent" or nil
    if group then hl[#hl + 1] = { row - 1, 2, #lines[#lines], group } end
  end
  local text, counts = M.summary(ctx.diffs, entries)
  put("")
  local srow = put("  " .. text)
  hl[#hl + 1] = { srow - 1, 0, #lines[#lines], "PlasmidHeader" }
  return lines, hl, index, counts
end

--- Header of the view: what was aligned against what, and how well.
--- `width` folds the facts onto as many lines as the window needs.
function M.header_lines(ctx, width)
  local aln = ctx.aln
  local lines = { ctx.title }
  local bits = { aln.mode }
  if not ctx.rec:is_protein() then
    bits[#bits + 1] = (aln.strand == -1 and "minus" or "plus") .. " strand"
  end
  bits[#bits + 1] = ("score %d"):format(aln.score)
  bits[#bits + 1] = ("identity %s (%d/%d)"):format(pct(aln.identity), aln.matches, aln.columns)
  bits[#bits + 1] = ("%s %d..%d%s"):format(ctx.ref_label, aln.ref_start, aln.ref_stop,
    aln.wraps and " (across the origin)" or "")
  bits[#bits + 1] = ("%s %d..%d of %d (%s)"):format(ctx.query_label, aln.query_start, aln.query_stop,
    ctx.query_len, pct(aln.query_coverage))
  if aln.rotation and aln.rotation ~= 0 then
    bits[#bits + 1] = ("rotated %d %s into phase"):format(aln.rotation, ctx.unit)
  end
  if aln.banded then bits[#bits + 1] = "banded" end
  local room = math.max(40, (width or vim.o.columns) - 2)
  local row = ""
  for _, bit in ipairs(bits) do
    local joined = row == "" and bit or (row .. " · " .. bit)
    if vim.fn.strdisplaywidth(joined) > room and row ~= "" then
      lines[#lines + 1] = " " .. row
      row = bit
    else
      row = joined
    end
  end
  if row ~= "" then lines[#lines + 1] = " " .. row end
  if aln.strand == -1 then
    lines[#lines + 1] = (" the query is shown reverse-complemented; its coordinates run along that strand")
  end
  return lines
end

--- Everything the view draws, from a context built by M.open_alignment.
function M.view_lines(ctx, width)
  local lines, hl, targets, list_line = {}, {}, {}, {}
  local function head(text)
    lines[#lines + 1] = ""
    lines[#lines + 1] = " " .. text
    hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "PlasmidHeader" }
  end
  for _, l in ipairs(M.header_lines(ctx, width)) do lines[#lines + 1] = l end
  hl[#hl + 1] = { 0, 0, #lines[1], "PlasmidTitle" }
  for i = 2, #lines do hl[#hl + 1] = { i - 1, 0, #lines[i], "Comment" } end

  head(("differences (%d)"):format(#ctx.diffs))
  local dlines, dhl, dindex, counts = M.difference_lines(ctx, #lines)
  for _, l in ipairs(dlines) do lines[#lines + 1] = l end
  for _, h in ipairs(dhl) do hl[#hl + 1] = h end
  for row, idx in pairs(dindex) do
    list_line[row] = idx
    targets[#targets + 1] = { line = row, col = 0, diff = idx }
  end
  ctx.counts = counts

  if #ctx.diffs > 0 then
    head(("features hit (%d)"):format(#(ctx.entries or {})))
    local flines, fhl = M.feature_lines(ctx.entries or {}, ctx.unit, #lines)
    for _, l in ipairs(flines) do lines[#lines + 1] = l end
    for _, h in ipairs(fhl) do hl[#hl + 1] = h end
  end

  head("alignment")
  local blines, bhl, marks = M.block_lines(ctx, width, #lines)
  for _, l in ipairs(blines) do lines[#lines + 1] = l end
  for _, h in ipairs(bhl) do hl[#hl + 1] = h end
  for _, m in ipairs(marks) do targets[#targets + 1] = m end
  table.sort(targets, function(a, b)
    if a.line ~= b.line then return a.line < b.line end
    return a.col < b.col
  end)

  local k = cfg().keymaps or {}
  lines[#lines + 1] = ""
  lines[#lines + 1] = (" %s / %s next / previous difference   %s go to that base   %s annotate as %s features   %s close")
    :format(k.next_diff or "]d", k.prev_diff or "[d", k.goto_base or "<CR>", k.annotate or "ga",
      cfg().annotate_type or "variation", k.close or "q")
  hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "Comment" }
  return lines, hl, targets, list_line, marks
end

---------------------------------------------------------------------------
-- Buffers
---------------------------------------------------------------------------

local function scratch(kind, src, state_table, filetype)
  for b, s in pairs(state_table) do
    if s.src == src and vim.api.nvim_buf_is_valid(b) then return b end
  end
  local b = vim.api.nvim_create_buf(false, true)
  vim.bo[b].bufhidden = "wipe"
  vim.bo[b].filetype = filetype
  pcall(vim.api.nvim_buf_set_name, b, ("plasmid://%s/%s"):format(kind, tostring(src or 0)))
  require("plasmid.util").case_insensitive_search(b)
  vim.api.nvim_create_autocmd("BufWipeout", { buffer = b, once = true,
    callback = function() state_table[b] = nil end })
  return b
end

local function paint(b, lines, hl)
  vim.bo[b].modifiable = true
  vim.api.nvim_buf_set_lines(b, 0, -1, false, lines)
  vim.bo[b].modifiable = false
  vim.api.nvim_buf_clear_namespace(b, ns, 0, -1)
  for _, h in ipairs(hl or {}) do
    pcall(vim.api.nvim_buf_set_extmark, b, ns, h[1], h[2], { end_col = h[3], hl_group = h[4] })
  end
end

local function open_split(b, lines, src)
  local win = vim.fn.bufwinid(b)
  if win == -1 then
    local srcwin = src and vim.fn.bufwinid(src) or -1
    if srcwin ~= -1 then vim.api.nvim_set_current_win(srcwin) end
    vim.cmd("botright " .. math.min(#lines + 1, math.floor(vim.o.lines * 0.55)) .. "split")
    vim.api.nvim_win_set_buf(0, b)
    win = vim.api.nvim_get_current_win()
    vim.wo[win].wrap = false
    vim.wo[win].number = false
    vim.wo[win].relativenumber = false
  end
  return win
end

--- Columns available to the view: the window that shows `b`, or the screen.
local function view_width(b)
  local win = vim.fn.bufwinid(b)
  if win ~= -1 then return vim.api.nvim_win_get_width(win) end
  return vim.o.columns
end

--- Draw a context in the alignment view for `src` and return the buffer.
function M.show(ctx)
  local b = scratch("align", ctx.src, M.state, "plasmid_align")
  local guessed = view_width(b)
  local lines, hl, targets, list_line = M.view_lines(ctx, guessed)
  paint(b, lines, hl)
  M.state[b] = { src = ctx.src, rec = ctx.rec, aln = ctx.aln, diffs = ctx.diffs, ctx = ctx,
    targets = targets, list_line = list_line, query_name = ctx.query_label }
  M.last = b
  local win = open_split(b, lines, ctx.src)
  -- a window of a different width than we guessed at: draw it again to fit
  if vim.api.nvim_win_get_width(win) ~= guessed then
    lines, hl, targets, list_line = M.view_lines(ctx, vim.api.nvim_win_get_width(win))
    paint(b, lines, hl)
    M.state[b].targets, M.state[b].list_line = targets, list_line
  end
  vim.api.nvim_win_set_cursor(win, { 1, 0 })
  M.set_keymaps(b)
  return b
end

function M.set_keymaps(b)
  local k = cfg().keymaps or {}
  local o = { buffer = b, nowait = true, silent = true }
  vim.keymap.set("n", k.next_diff or "]d", function() M.jump(b, 1) end,
    vim.tbl_extend("force", o, { desc = "plasmid: next difference" }))
  vim.keymap.set("n", k.prev_diff or "[d", function() M.jump(b, -1) end,
    vim.tbl_extend("force", o, { desc = "plasmid: previous difference" }))
  vim.keymap.set("n", k.goto_base or "<CR>", function() M.goto_difference(b) end,
    vim.tbl_extend("force", o, { desc = "plasmid: go to this base" }))
  vim.keymap.set("n", k.annotate or "ga", function() M.annotate(b) end,
    vim.tbl_extend("force", o, { desc = "plasmid: annotate the differences" }))
  vim.keymap.set("n", k.close or "q", "<cmd>close<cr>", o)
end

local function cursor_of(b)
  local win = vim.fn.bufwinid(b)
  if win == -1 then return 1, 0 end
  local c = vim.api.nvim_win_get_cursor(win)
  return c[1], c[2]
end

--- ]d / [d: the next (dir = 1) or previous difference, in the list first and
--- then through the alignment blocks.
function M.jump(b, dir)
  local S = M.state[b]
  if not S then return end
  if #S.targets == 0 then return notify("no differences in this alignment") end
  local line, col = cursor_of(b)
  local pick
  if dir > 0 then
    for _, t in ipairs(S.targets) do
      if t.line > line or (t.line == line and t.col > col) then pick = t break end
    end
    pick = pick or S.targets[1]
  else
    for i = #S.targets, 1, -1 do
      local t = S.targets[i]
      if t.line < line or (t.line == line and t.col < col) then pick = t break end
    end
    pick = pick or S.targets[#S.targets]
  end
  local win = vim.fn.bufwinid(b)
  if win == -1 then return end
  vim.api.nvim_win_set_cursor(win, { pick.line, pick.col })
  local d = S.diffs[pick.diff]
  if d then notify(("difference %d of %d: %s"):format(pick.diff, #S.diffs, V.describe(d, S.ctx.unit))) end
  return pick
end

--- The difference under the cursor: a difference-list row, or the column the
--- cursor sits on (or nearest to) in an alignment block.
function M.difference_at_cursor(b)
  local S = M.state[b]
  if not S then return nil end
  local line, col = cursor_of(b)
  if S.list_line[line] then return S.list_line[line] end
  local best
  for _, t in ipairs(S.targets) do
    if t.line == line and t.len then
      if col >= t.col and col < t.col + t.len then return t.diff end
      if not best or math.abs(t.col - col) < math.abs(best.col - col) then best = t end
    end
  end
  return best and best.diff
end

--- <CR>: select the difference's bases in the record's buffer and go there.
function M.goto_difference(b)
  local S = M.state[b]
  if not S then return end
  local idx = M.difference_at_cursor(b)
  local d = idx and S.diffs[idx]
  if not d then
    return notify("put the cursor on a difference (]d finds the next one)", vim.log.levels.WARN)
  end
  local src = S.src
  local st = src and vim.api.nvim_buf_is_valid(src) and buffer().get(src)
  if not st then return notify("the record's buffer is gone", vim.log.levels.WARN) end
  local win = vim.fn.bufwinid(src)
  if win == -1 then return notify("the record is not on screen", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  local a = math.max(1, d.pos)
  local z = d.kind == "ins" and a or math.min(rec:len(), d.pos + d.len - 1)
  vim.api.nvim_set_current_win(win)
  local sel = require("plasmid.select")
  sel.set(src, sel.make(rec, a, z))
  buffer().goto_bp(win, a)
  notify(("%s — %s"):format(V.describe(d, S.ctx.unit), rec.name))
  return a
end

---------------------------------------------------------------------------
-- Annotating the differences
---------------------------------------------------------------------------

local function variation_location(rec, d)
  local L = rec:len()
  if d.kind == "ins" then
    local a = d.pos
    if a < 1 then a = rec:is_circular() and L or 1 end
    if not rec:is_circular() and a >= L then a = L - 1 end
    local z = rec:is_circular() and (a % L + 1) or (a + 1)
    return ("%d^%d"):format(a, z)
  end
  local a, z = d.pos, d.pos + math.max(d.len, 1) - 1
  if z > L then
    if rec:is_circular() then return ("join(%d..%d,%d..%d)"):format(a, L, 1, z - L) end
    z = L
  end
  return ("%d..%d"):format(a, z)
end

--- Short name of a difference, used as the feature's /label.
function M.difference_label(d, unit)
  if d.kind == "sub" then return ("%s%d%s"):format(d.ref, d.pos, d.query) end
  if d.kind == "del" then return ("del %d %s @%d"):format(d.len, unit or "bp", d.pos) end
  return ("ins %d %s @%d"):format(d.len, unit or "bp", d.pos)
end

--- The /note: the change, then what it does to the feature it is in.
function M.difference_note(d, unit, query_name)
  local parts = { V.describe(d, unit) }
  if query_name then parts[#parts + 1] = ("from " .. query_name) end
  return table.concat(parts, "; ")
end

--- :DnaAlignAnnotate — write the last alignment's differences onto the record
--- as features (align.annotate_type, `variation` by default) in one undo step.
--- `b` may be the view or the record's own buffer.
function M.annotate(b)
  local S = M.state[b]
  if not S then
    for vb, s in pairs(M.state) do
      if s.src == b and vim.api.nvim_buf_is_valid(vb) then S = s break end
    end
  end
  if not S then
    return notify("no alignment to annotate here (run :DnaAlign or :DnaCompare first)", vim.log.levels.WARN)
  end
  local src = S.src
  local st = src and vim.api.nvim_buf_is_valid(src) and buffer().get(src)
  if not st then return notify("the record's buffer is gone", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  if rec ~= S.rec then
    return notify("the record changed since this alignment; run it again", vim.log.levels.WARN)
  end
  if #S.diffs == 0 then return notify("no differences to annotate") end
  local Feature = require("plasmid.model.record").Feature
  local ftype = cfg().annotate_type or "variation"
  local unit = S.ctx.unit
  local ops = {}
  for _, d in ipairs(S.diffs) do
    local quals = {
      { key = "label", value = M.difference_label(d, unit), quoted = true },
      { key = "note", value = M.difference_note(d, unit, S.query_name), quoted = true },
    }
    if d.kind ~= "ins" then
      quals[#quals + 1] = { key = "replace", value = d.query:lower(), quoted = true }
    end
    local f = Feature.new({ type = ftype, location_text = variation_location(rec, d), qualifiers = quals })
    f.pinned = true
    ops[#ops + 1] = { op = "add", feature = f }
  end
  require("plasmid.author").commit(src, { op = "batch", ops = ops },
    ("annotate %d difference%s"):format(#ops, #ops == 1 and "" or "s"))
  notify(("added %d %s feature%s (u undoes all of them)"):format(#ops, ftype, #ops == 1 and "" or "s"))
  return #ops
end

---------------------------------------------------------------------------
-- :DnaAlign
---------------------------------------------------------------------------

--- Build the view context for an alignment of `aln` (from M.align_query)
--- against `rec`.
local function align_context(src, rec, aln, diffs, query_name, query_len)
  local L = rec:len()
  local wrapped = aln.wraps and aln.ref_start_abs
  return {
    src = src, rec = rec, aln = aln, diffs = diffs,
    unit = rec:unit(),
    ref_label = rec.name, query_label = query_name,
    ref_len = L, query_len = query_len,
    ref_from = wrapped and aln.ref_start_abs or aln.ref_start,
    query_from = aln.query_start,
    ref_map = function(p) return (p - 1) % L + 1 end,
    query_map = function(p) return p end,
    title = (" Alignment: %s (%d %s %s)  ×  %s (%d %s)"):format(rec.name, L, rec:unit(), rec.topology,
      query_name, query_len, rec:unit()),
    kind = "align",
  }
end

--- :DnaAlign {query} [mode=] [strand=] [name=]
function M.align(buf, fargs)
  local st = buffer().get(buf)
  if not st then return notify("current buffer is not a sequence file", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  local o, rest, err = M.parse_args(fargs)
  if not o then return notify(err, vim.log.levels.ERROR) end
  if #rest == 0 then
    return notify(":DnaAlign {sequence|file|buffer|part|selection} [mode=semi|global|local] [strand=auto|plus|minus] [name=]",
      vim.log.levels.WARN)
  end
  local q, qerr = M.resolve_query(buf, table.concat(rest, " "), rec)
  if not q and #rest > 1 then q, qerr = M.resolve_query(buf, rest[1], rec) end
  if not q then return notify(qerr, vim.log.levels.ERROR) end
  local name = o.name or q.name
  local aln, why = M.align_query(rec, q.seq, o)
  if not aln then
    return notify(("%s does not align to %s: %s"):format(name, rec.name, tostring(why)), vim.log.levels.WARN)
  end
  local diffs = M.differences(rec, aln)
  V.annotate(rec, diffs)
  local ctx = align_context(buf, rec, aln, diffs, name, #aln.query_seq)
  local b = M.show(ctx)
  notify(("%s placed at %d..%d (%s identity, %d difference%s)"):format(name, aln.ref_start, aln.ref_stop,
    pct(aln.identity), #diffs, #diffs == 1 and "" or "s"))
  return b, aln, diffs
end

---------------------------------------------------------------------------
-- :DnaCompare
---------------------------------------------------------------------------

--- :DnaCompare {record} — whole-molecule comparison, origin and strand
--- independent on circular records.
function M.compare(buf, fargs)
  local st = buffer().get(buf)
  if not st then return notify("current buffer is not a sequence file", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  local o, rest, err = M.parse_args(fargs)
  if not o then return notify(err, vim.log.levels.ERROR) end
  if #rest == 0 then return notify(":DnaCompare {record|file|buffer|part}", vim.log.levels.WARN) end
  local q, qerr = M.resolve_query(buf, table.concat(rest, " "), rec)
  if not q and #rest > 1 then q, qerr = M.resolve_query(buf, rest[1], rec) end
  if not q then return notify(qerr, vim.log.levels.ERROR) end
  local name = o.name or q.name
  local a = rec.seq:tostring():upper()
  local b_seq = q.seq:upper():gsub("[^%a%*]", "")
  local protein = rec:is_protein()
  local opts = score_opts({ mode = o.mode == "semi" and "global" or o.mode, protein = protein or nil,
    circular = rec:is_circular() and not protein })
  local aln, diffs = A.compare_differences(a, b_seq, opts)
  if not aln then
    return notify(("%s and %s do not align: %s"):format(rec.name, name, tostring(diffs)), vim.log.levels.WARN)
  end
  aln.query_seq = aln.query_rotated or b_seq
  V.annotate(rec, diffs)
  local L, rot = aln.ref_length or #a, aln.ref_rotation or 0
  local Lq = #(aln.query_rotated or b_seq)
  local qrot = aln.rotation or 0
  local ctx = {
    src = buf, rec = rec, aln = aln, diffs = diffs, unit = rec:unit(),
    ref_label = rec.name, query_label = name,
    ref_len = L, query_len = Lq,
    ref_from = aln.ref_start, query_from = aln.query_start,
    ref_map = function(p) return (p - 1 + rot) % L + 1 end,
    query_map = function(p) return (p - 1 + qrot) % Lq + 1 end,
    title = (" Comparison: %s (%d %s %s)  ↔  %s (%d %s)"):format(rec.name, L, rec:unit(), rec.topology,
      name, Lq, rec:unit()),
    kind = "compare",
  }
  local vb = M.show(ctx)
  local how = {}
  if (aln.rotation or 0) ~= 0 then how[#how + 1] = ("rotated %d %s"):format(aln.rotation, rec:unit()) end
  if aln.strand == -1 then how[#how + 1] = "reverse-complemented" end
  notify(("%s vs %s: %s identity, %d difference%s%s"):format(rec.name, name, pct(aln.identity), #diffs,
    #diffs == 1 and "" or "s", #how > 0 and (" (the query was " .. table.concat(how, " and ") .. ")") or ""))
  return vb, aln, diffs
end

---------------------------------------------------------------------------
-- :DnaAlignReads
---------------------------------------------------------------------------

--- Every record of every file an argument names. Returns reads
--- ({ name, seq, file }) and the arguments that matched nothing.
function M.read_files(buf, args, rec)
  local reads, missing = {}, {}
  for _, a in ipairs(args) do
    -- glob expands the wildcards itself; expand() first would collapse a
    -- multi-file pattern into one space-joined string
    local pat = (a:sub(1, 1) == "~" or a:sub(1, 1) == "$") and vim.fn.expand(a) or a
    local paths = vim.fn.glob(pat, false, true)
    if #paths == 0 and vim.fn.filereadable(pat) == 1 then paths = { pat } end
    if #paths == 0 then
      local q = M.resolve_query(buf, a, rec)
      if q then reads[#reads + 1] = { name = q.name, seq = q.seq, file = a }
      else missing[#missing + 1] = a end
    end
    for _, p in ipairs(paths) do
      local res, rerr = require("plasmid.io").read(p)
      if not res or not res.records or #res.records == 0 then
        missing[#missing + 1] = ("%s (%s)"):format(a, tostring(rerr))
      else
        local base = vim.fn.fnamemodify(p, ":t")
        for _, r in ipairs(res.records) do
          reads[#reads + 1] = { name = r.name, seq = r.seq:tostring(), file = base, rec = r }
        end
      end
    end
  end
  return reads, missing
end

--- Place one read. Returns a row, or nil and the reason it is unplaced.
function M.place_read(rec, read, o)
  local c = cfg()
  local min_id = o.min_identity or c.min_identity or 0.9
  local min_cov = o.min_coverage or c.min_coverage or 0.5
  local aln, why = M.align_query(rec, read.seq, o)
  if not aln then return nil, tostring(why) end
  if aln.identity < min_id then
    return nil, ("identity %s is below align.min_identity (%s)"):format(pct(aln.identity), pct(min_id))
  end
  if aln.query_coverage < min_cov then
    return nil, ("only %s of the read aligns, below align.min_coverage (%s)")
      :format(pct(aln.query_coverage), pct(min_cov))
  end
  local diffs = M.differences(rec, aln)
  V.annotate(rec, diffs)
  return { name = read.name, file = read.file, aln = aln, diffs = diffs, seq = read.seq,
    length = #aln.query_seq }
end

local function read_table_lines(rec, rows, unplaced)
  local header = { "#", "read", "strand", "from", "to", "length", "identity", "diffs", "coverage" }
  local body = {}
  for i, r in ipairs(rows) do
    body[i] = { tostring(i), r.name, r.aln.strand == -1 and "−" or "+",
      tostring(r.aln.ref_start), tostring(r.aln.ref_stop), ("%d %s"):format(r.length, rec:unit()),
      pct(r.aln.identity), tostring(#r.diffs), pct(r.aln.query_coverage) }
  end
  local right = { true, false, false, true, true, true, true, true, true }
  local w = {}
  for c = 1, #header do
    w[c] = vim.fn.strdisplaywidth(header[c])
    for _, b in ipairs(body) do w[c] = math.max(w[c], vim.fn.strdisplaywidth(b[c])) end
  end
  local function fmt(r)
    local out = {}
    for c = 1, #header do
      out[c] = right[c] and rpad(r[c], w[c]) or (c < #header and pad(r[c], w[c]) or r[c])
    end
    return "  " .. table.concat(out, "  ")
  end
  local lines = { fmt(header) }
  for _, r in ipairs(body) do lines[#lines + 1] = fmt(r) end
  return lines
end

--- :DnaAlignReads {file|glob}...
function M.align_reads(buf, fargs)
  local st = buffer().get(buf)
  if not st then return notify("current buffer is not a sequence file", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  local o, rest, err = M.parse_args(fargs)
  if not o then return notify(err, vim.log.levels.ERROR) end
  if #rest == 0 then
    return notify(":DnaAlignReads {file|glob}... [mode=] [strand=] [min_identity=]", vim.log.levels.WARN)
  end
  local reads, missing = M.read_files(buf, rest, rec)
  if #reads == 0 then
    return notify(("no reads in %s"):format(table.concat(missing, ", ")), vim.log.levels.ERROR)
  end
  local rows, unplaced = {}, {}
  for _, read in ipairs(reads) do
    local row, why = M.place_read(rec, read, o)
    if row then rows[#rows + 1] = row else unplaced[#unplaced + 1] = { name = read.name, why = why } end
  end
  table.sort(rows, function(a, b)
    if a.aln.ref_start ~= b.aln.ref_start then return a.aln.ref_start < b.aln.ref_start end
    return a.name < b.name
  end)

  local lines = { (" Reads on %s (%d %s %s) — %d of %d placed"):format(rec.name, rec:len(), rec:unit(),
    rec.topology, #rows, #reads) }
  local hl = { { 0, 0, #lines[1], "PlasmidTitle" } }
  lines[#lines + 1] = ""
  local first = #lines + 1
  local tlines = read_table_lines(rec, rows, unplaced)
  hl[#hl + 1] = { #lines, 0, #tlines[1], "PlasmidHeader" }
  for _, l in ipairs(tlines) do lines[#lines + 1] = l end
  if #rows > 0 then
    local total = 0
    for _, r in ipairs(rows) do total = total + #r.diffs end
    lines[#lines + 1] = ""
    lines[#lines + 1] = ("  %d difference%s over %d placed read%s"):format(total, total == 1 and "" or "s",
      #rows, #rows == 1 and "" or "s")
    hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "PlasmidHeader" }
  end
  if #unplaced > 0 then
    lines[#lines + 1] = ""
    lines[#lines + 1] = " not placed"
    hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "PlasmidHeader" }
    for _, u in ipairs(unplaced) do
      lines[#lines + 1] = ("  %s — %s"):format(u.name, u.why)
      hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "ErrorMsg" }
    end
  end
  for _, m in ipairs(missing) do
    lines[#lines + 1] = ("  %s — no such file"):format(m)
    hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "ErrorMsg" }
  end
  lines[#lines + 1] = ""
  lines[#lines + 1] = " <CR> open this read's alignment   q close"
  hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "Comment" }

  local b = scratch("alignreads", buf, M.read_state, "plasmid_alignreads")
  paint(b, lines, hl)
  M.read_state[b] = { src = buf, rec = rec, rows = rows, unplaced = unplaced, first = first + 1 }
  M.last_reads = b
  local win = open_split(b, lines, buf)
  vim.api.nvim_win_set_cursor(win, { math.min(first + 1, #lines), 0 })
  local ko = { buffer = b, nowait = true, silent = true }
  vim.keymap.set("n", "<CR>", function() M.open_read(b) end,
    vim.tbl_extend("force", ko, { desc = "plasmid: open this read's alignment" }))
  vim.keymap.set("n", "q", "<cmd>close<cr>", ko)
  notify(("%d of %d read%s placed on %s"):format(#rows, #reads, #reads == 1 and "" or "s", rec.name))
  return b, rows, unplaced
end

--- <CR> in the read table: the alignment of the read under the cursor.
function M.open_read(b, n)
  local S = M.read_state[b]
  if not S then return end
  if not n then
    local win = vim.fn.bufwinid(b)
    if win == -1 then return end
    n = vim.api.nvim_win_get_cursor(win)[1] - S.first + 1
  end
  local row = S.rows[n]
  if not row then
    return notify(("put the cursor on a read row (1..%d)"):format(#S.rows), vim.log.levels.WARN)
  end
  local ctx = align_context(S.src, S.rec, row.aln, row.diffs, row.name, row.length)
  local vb = M.show(ctx)
  return vb, row
end

return M
