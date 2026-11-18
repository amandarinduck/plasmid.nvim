local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local scratch = require("plasmid.scratch")
local tx = require("plasmid.model.transaction")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local n_open = 0
local function open(rel)
  n_open = n_open + 1
  vim.cmd("silent! only")
  local dst = tmp .. "/" .. n_open .. "_" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.get(buf), buffer.record(buffer.get(buf))
end
local function cds(rec, gene)
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:get("gene") == gene then return f end end
end
local function rnd(n)
  local t = {}
  for i = 1, n do local r = math.random(4); t[i] = ("ACGT"):sub(r, r) end
  return table.concat(t)
end
-- scratch text helpers: residues as one string, and setting it back in 60-char lines
local function get_text(sbuf) return table.concat(vim.api.nvim_buf_get_lines(sbuf, 0, -1, false), "") end
local function set_text(sbuf, s)
  local lines = {}
  for i = 1, #s, 60 do lines[#lines + 1] = s:sub(i, i + 59) end
  vim.api.nvim_buf_set_lines(sbuf, 0, -1, false, #lines > 0 and lines or { "" })
end
local function last_msg() return H.messages[#H.messages] and H.messages[#H.messages].msg or "" end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("diff_reconstructs_and_is_minimal", function()
  math.randomseed(7)
  for trial = 1, 300 do
    local old = rnd(math.random(0, 400))
    local new, nedits = old, math.random(0, 6)
    for _ = 1, nedits do
      local p = math.random(0, #new)
      local r = math.random(3)
      if r == 1 then new = new:sub(1, p) .. rnd(math.random(1, 8)) .. new:sub(p + 1)
      elseif r == 2 then new = new:sub(1, p) .. new:sub(p + math.random(1, 8) + 1)
      else new = new:sub(1, p) .. rnd(math.random(1, 5)) .. new:sub(p + math.random(1, 5) + 1) end
    end
    local hunks = scratch.diff(old, new)
    H.eq(new, scratch.patch(old, hunks), "trial " .. trial)
    H.ok(#hunks <= math.max(nedits, 0) * 3 + (nedits > 0 and 1 or 0), ("trial %d: %d hunks for %d edits"):format(trial, #hunks, nedits))
  end
  -- unrelated texts beyond maxd fall back to one replacement, still correct
  local a, b = rnd(500), rnd(500)
  local hunks = scratch.diff(a, b, 10)
  H.eq(b, scratch.patch(a, hunks))
  H.eq({}, scratch.diff("ACGT", "ACGT"))
end)

case("edit_region_round_trip", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  vim.cmd("DnaEdit 450..650")
  local sbuf = vim.api.nvim_get_current_buf()
  H.ok(sbuf ~= buf, "scratch buffer opened")
  H.eq(orig:sub(450, 650):upper(), get_text(sbuf), "region shown as displayed")
  H.eq(60, #vim.api.nvim_buf_get_lines(sbuf, 0, 1, false)[1], "60 residues per line")
  H.ok(scratch.winbar(sbuf):find("editing SYNPBR322 450..650", 1, true), scratch.winbar(sbuf))
  -- insert GCT before bp 500 (offset 50 in the region) and write
  local s = get_text(sbuf)
  set_text(sbuf, s:sub(1, 50) .. "GCT" .. s:sub(51))
  H.messages = {}
  vim.cmd("write")
  H.eq(orig:sub(1, 499) .. "gct" .. orig:sub(500), rec.seq:tostring(), "applied as one insertion before 500")
  H.eq("86..1279", cds(rec, "tet").location_text, "tet extended, not truncated")
  H.ok(last_msg():find("tet  CDS  extended +3 bp", 1, true), last_msg())
  H.eq(false, vim.bo[sbuf].modified)
  H.eq(653, vim.b[sbuf].plasmid_edit.stop, "region now ends 3 later")
  -- writing again without changes does nothing
  local depth = #rec.history.undo
  vim.cmd("write")
  H.eq(depth, #rec.history.undo)
  H.ok(last_msg():find("no changes", 1, true))
  -- one undo in the sequence view
  vim.api.nvim_set_current_win(vim.fn.bufwinid(buf))
  require("plasmid.author").undo()
  H.eq(orig, rec.seq:tostring())
  -- the scratch buffer is now stale: refused
  vim.api.nvim_set_current_win(vim.fn.bufwinid(sbuf))
  set_text(sbuf, "ACGT")
  vim.cmd("write")
  H.ok(last_msg():find("sequence changed since :DnaEdit", 1, true), last_msg())
  H.eq(orig, rec.seq:tostring())
  vim.cmd("bwipeout!")
end)

case("features_inside_the_region_keep_their_bases", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  -- a small feature well inside a region
  local inner
  for _, f in ipairs(rec.features) do
    local lo, hi = f:bounds()
    if f.type ~= "source" and lo and lo > 1500 and hi < 2600 and hi - lo >= 10 and #f.parts == 1 then inner = f break end
  end
  H.ok(inner, "a feature inside 1400..2700")
  local lo, hi = inner:bounds()
  local inner_res = rec:extract(inner)
  vim.cmd("DnaEdit 1400..2700")
  local sbuf = vim.api.nvim_get_current_buf()
  local s = get_text(sbuf)
  local a, b = lo - 1400, hi - 1400 + 1 -- offsets of the feature in the region
  -- edits on both sides of the feature, none inside it
  local new = s:sub(1, 20) .. "TTTT" .. s:sub(21, a - 5) .. s:sub(a - 2, b) .. "GG" .. s:sub(b + 1, #s - 30) .. s:sub(#s - 20)
  set_text(sbuf, new)
  vim.cmd("write")
  H.eq(inner_res, rec:extract(inner.id and rec:feature_by_id(inner.id) or inner), "feature bases unchanged")
  H.eq(1400 + #new - 1, vim.b[sbuf].plasmid_edit.stop)
  H.eq(new:lower(), rec.seq:sub(1400, 1400 + #new - 1), "region now holds the new text")
  vim.cmd("bwipeout!")
end)

case("rewrite_invalid_delete_and_feature_changes", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  vim.cmd("DnaEdit 3000..3100")
  local sbuf = vim.api.nvim_get_current_buf()
  -- invalid letters: refused, buffer stays modified
  set_text(sbuf, "ACGTXX")
  H.messages = {}
  vim.cmd("write")
  H.ok(last_msg():find("not a valid", 1, true), last_msg())
  H.eq(orig, rec.seq:tostring())
  H.eq(true, vim.bo[sbuf].modified)
  -- a feature-only change in the record does not make the scratch stale
  local f = require("plasmid.model.record").Feature.new({ type = "misc_feature", location_text = "10..20",
    qualifiers = { { key = "label", value = "mine", quoted = true } } })
  tx.commit(rec, { op = "add", feature = f }, "add mine")
  -- complete rewrite with unrelated residues
  math.randomseed(11)
  local new = rnd(150)
  set_text(sbuf, new)
  vim.cmd("write")
  H.eq(orig:sub(1, 2999) .. new:lower() .. orig:sub(3101), rec.seq:tostring())
  -- delete the whole region
  set_text(sbuf, "")
  vim.cmd("write")
  H.eq(orig:sub(1, 2999) .. orig:sub(3101), rec.seq:tostring())
  H.eq(2999, vim.b[sbuf].plasmid_edit.stop, "empty region")
  -- and type into the empty region again
  set_text(sbuf, "ACGT")
  vim.cmd("write")
  H.eq(orig:sub(1, 2999) .. "acgt" .. orig:sub(3101), rec.seq:tostring())
  vim.cmd("bwipeout!")
end)

case("whole_record_and_errors", function()
  local p = tmp .. "/small.fa"
  local fh = assert(io.open(p, "w")); fh:write(">small\nACGTacgtACGT\n"); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(p))
  local buf = vim.api.nvim_get_current_buf()
  local rec = buffer.record(buffer.get(buf))
  vim.cmd("DnaEdit all")
  local sbuf = vim.api.nvim_get_current_buf()
  H.eq("ACGTacgtACGT", get_text(sbuf), "FASTA soft-masking kept in the scratch buffer")
  set_text(sbuf, "ACGTacgtACGTtt")
  vim.cmd("write")
  H.eq("ACGTacgtACGTtt", rec.seq:tostring(), "typed case kept for FASTA")
  set_text(sbuf, "")
  H.messages = {}
  vim.cmd("write")
  H.ok(last_msg():find("whole sequence", 1, true), last_msg())
  H.eq("ACGTacgtACGTtt", rec.seq:tostring())
  vim.cmd("bwipeout!")
  -- no selection, no range (the write above selected the edited region)
  vim.api.nvim_set_current_buf(buf)
  H.ok(buffer.get(buf).selection, "the written region is selected")
  vim.cmd("DnaSelectClear")
  H.messages = {}
  vim.cmd("DnaEdit")
  H.ok(last_msg():find(":DnaEdit", 1, true))
end)


case("refuse_policy_and_force_write", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local orig = rec.seq:tostring()
  require("plasmid.config").options.edit.cut_policy = "refuse"
  vim.cmd("DnaEdit 480..520")
  local sbuf = vim.api.nvim_get_current_buf()
  set_text(sbuf, get_text(sbuf):sub(1, 10) .. get_text(sbuf):sub(20)) -- cuts tet
  H.messages = {}
  vim.cmd("write")
  H.ok(last_msg():find("refused", 1, true), last_msg())
  H.eq(orig, rec.seq:tostring())
  H.eq(true, vim.bo[sbuf].modified, "still unsaved")
  vim.cmd("write!")
  H.eq(#orig - 9, rec:len(), ":w! forces it")
  H.eq(false, vim.bo[sbuf].modified)
  require("plasmid.config").options.edit.cut_policy = "truncate"
  vim.cmd("bwipeout!")
end)

return T
