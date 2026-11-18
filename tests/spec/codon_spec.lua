local H = require("helpers")
local alphabet = require("plasmid.model.alphabet")
local translate = require("plasmid.bio.translate")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local codon = require("plasmid.codon")
local sel = require("plasmid.select")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.get(buf), buffer.record(buffer.get(buf))
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end

-- independent complement for the oracle side of the checks
local COMP = { A = "T", C = "G", G = "C", T = "A", N = "N" }
local function comp(s) return (s:upper():gsub(".", function(c) return COMP[c] or c end)) end

vim.o.columns, vim.o.lines = 120, 50

-- Put the cursor on every amino acid drawn in the buffer, let the highlight
-- update, and check: (1) the three highlighted top-strand bases, read on the
-- feature's strand, are the codon reported; (2) that codon translates to the
-- residue under the cursor; (3) assembling residues by codon number reproduces
-- the file's /translation.
local function sweep(rel)
  local buf, st, rec = open(rel)
  vim.cmd("DnaTracks all")
  local win = vim.api.nvim_get_current_win()
  local lay = st.layout
  local by_f, split, checked = {}, 0, 0
  for row0 = 0, #lay.lines - 1 do
    local info = lay.rows[row0]
    if info and info.kind == "aa" then
      local line = lay.lines[row0 + 1]
      for vcol = 0, #line - 1 do
        local ch = line:sub(vcol + 1, vcol + 1)
        if ch:match("%S") then
          vim.api.nvim_win_set_cursor(win, { row0 + 1, vcol })
          codon.update(buf, win)
          local c = st.codon
          H.ok(c, ("%s: no codon at row %d col %d"):format(rel, row0, vcol))
          H.eq(ch, c.aa, "residue under cursor")
          -- (1) highlighted bases
          local marks = vim.api.nvim_buf_get_extmarks(buf, codon.ns, 0, -1, { details = true })
          local top, ncomp, naa = {}, 0, 0
          for _, m in ipairs(marks) do
            local k = lay.rows[m[2]] and lay.rows[m[2]].kind
            if k == "seq" then
              local bp = lay:bp_at(m[2], m[3])
              top[bp] = lay.lines[m[2] + 1]:sub(m[3] + 1, m[3] + 1)
            elseif k == "comp" then ncomp = ncomp + 1
            elseif k == "aa" then naa = naa + 1 end
          end
          H.eq(1, naa, "residue highlighted")
          H.eq(3, ncomp, "complement bases highlighted")
          local p = c.positions
          local read = (top[p[1]] or "?") .. (top[p[2]] or "?") .. (top[p[3]] or "?")
          if c.f.strand == -1 then read = comp(read) end
          H.eq(c.codon, read:upper(), "highlighted bases spell the codon")
          -- (2) translation
          local t = translate.get_table(c.f:get("transl_table"))
          local aa = translate.codon(t, c.codon)
          if not (c.n == 1 and c.aa == "M" and translate.is_start(t, c.codon)) then
            H.eq(aa, c.aa, ("codon %s"):format(c.codon))
          end
          if math.abs(p[3] - p[1]) ~= 2 then
            split = split + 1
            H.ok(codon.describe(st):find("bp " .. p[1] .. "," .. p[2] .. "," .. p[3], 1, true), codon.describe(st))
          end
          by_f[c.f] = by_f[c.f] or {}
          by_f[c.f][c.n] = c.aa
          checked = checked + 1
        end
      end
    end
  end
  -- (3) /translation
  local ncds = 0
  for f, res in pairs(by_f) do
    local want = f:get("translation")
    if want and f.type == "CDS" and not f:get("transl_except") then
      local got = {}
      for n = 1, #want do got[n] = res[n] or "?" end
      H.eq(want, table.concat(got), f:label() .. " residues by codon number")
      ncds = ncds + 1
    end
  end
  return checked, split, ncds
end

case("every_residue_highlights_its_codon", function()
  local total, report = 0, {}
  for _, rel in ipairs({ "ncbi/J01749.gb", "ncbi/NC_001422.gb", "ncbi/U01317.gb" }) do
    local n, split, ncds = sweep(rel)
    total = total + n
    report[#report + 1] = ("%s: %d residues, %d split codons, %d CDS match /translation"):format(
      vim.fn.fnamemodify(rel, ":t:r"), n, split, ncds)
    if rel:find("NC_001422") then H.ok(split > 0, "origin-crossing codon expected in phiX174") end
    if rel:find("U01317") then H.ok(split > 0, "intron-split codons expected in HBB") end
    H.ok(ncds > 0, rel)
  end
  -- one short line each: a headless message wider than 'columns' can crash a
  -- later redraw in Neovim 0.12.2
  for _, line in ipairs(report) do H.log("    [codon] " .. line) end
end)

case("winbar_and_clearing", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local win = vim.api.nvim_get_current_win()
  local tet
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:label() == "tet" then tet = f end end
  -- first residue of tet: the codon at 86..88
  local lay = st.layout
  local r, v = lay:pos_of(87)
  local row
  for k = r + 1, r + 12 do
    if lay.rows[k] and lay.rows[k].kind == "aa" and lay.lines[k + 1]:sub(v + 1, v + 1) == "M" then row = k break end
  end
  H.ok(row, "tet's first residue drawn")
  vim.api.nvim_win_set_cursor(win, { row + 1, v })
  codon.update(buf, win)
  local wb = buffer.winbar()
  H.ok(wb:find("codon 1/" .. math.floor((1276 - 86 + 1) / 3) .. " ATG→M 86..88", 1, true), wb)
  -- moving off the translation row clears it
  vim.api.nvim_win_set_cursor(win, { r + 1, v })
  codon.update(buf, win)
  H.eq(0, #vim.api.nvim_buf_get_extmarks(buf, codon.ns, 0, -1, {}))
  H.ok(not buffer.winbar():find("codon", 1, true))
  -- switching view clears it too
  vim.api.nvim_win_set_cursor(win, { row + 1, v })
  codon.update(buf, win)
  buffer.set_view(buf, "circular")
  H.eq(0, #vim.api.nvim_buf_get_extmarks(buf, codon.ns, 0, -1, {}), "no codon marks on the map")
  buffer.set_view(buf, "linear")
end)

case("cursor_moved_autocmd_drives_it", function()
  local buf, st = open("ncbi/J01749.gb")
  local lay = st.layout
  local row, col
  for k = 0, #lay.lines - 1 do
    if lay.rows[k] and lay.rows[k].kind == "aa" then
      local c = lay.lines[k + 1]:find("%u")
      if c then row, col = k, c - 1 break end
    end
  end
  vim.api.nvim_win_set_cursor(0, { row + 1, col })
  vim.api.nvim_exec_autocmds("CursorMoved", { buffer = buf })
  H.ok(st.codon, "autocmd set the codon")
  H.eq(7, #vim.api.nvim_buf_get_extmarks(buf, codon.ns, 0, -1, {}), "residue + 3 top + 3 bottom")
end)

case("strand_orientations", function()
  H.eq("AACG", sel.orient("AACG", "top"))
  H.eq("GCAA", sel.orient("AACG", "top35"))
  H.eq("CGTT", sel.orient("AACG", "bottom"))
  H.eq("TTGC", sel.orient("AACG", "bottom35"))
  for alias, key in pairs({ revcomp = "bottom", rc = "bottom", complement = "bottom35", reverse = "top35",
    top53 = "top", bottom53 = "bottom", TOP = "top" }) do
    H.eq(key, sel.strand_key(alias), alias)
  end
  H.eq(nil, sel.strand_key("sideways"))
end)

case("dnayank_command_all_orientations", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local plus = rec.seq:sub(10, 40):upper()
  vim.cmd("DnaSelect 10..40")
  local want = {
    top = plus, top35 = plus:reverse(),
    bottom = comp(plus):reverse(), bottom35 = comp(plus),
  }
  for key, w in pairs(want) do
    vim.cmd("DnaYank " .. key)
    H.eq(w, vim.fn.getreg('"'):upper(), key)
  end
  vim.cmd("DnaYank complement a")
  H.eq(want.bottom35, vim.fn.getreg("a"):upper(), "named register")
  -- origin-wrapping selection on a circular record
  vim.cmd("DnaSelect 4350..10")
  vim.cmd("DnaYank bottom")
  local wrap = (rec.seq:sub(4350) .. rec.seq:sub(1, 10)):upper()
  H.eq(comp(wrap):reverse(), vim.fn.getreg('"'):upper(), "wrapping bottom 5'->3'")
  H.eq(22, #wrap, "4350..4361 + 1..10")
  -- bad argument
  H.messages = {}
  vim.cmd("DnaYank sideways")
  local msgs = {}
  for _, m in ipairs(H.messages) do msgs[#msgs + 1] = m.msg end
  H.ok(table.concat(msgs, "\n"):find(":DnaYank", 1, true), "usage message")
end)

case("picker_marks_feature_strand", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local bla
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f.strand == -1 then bla = f break end end
  require("plasmid.ui").goto_feature(buf, bla.id)
  keys("vif")
  local s = sel.get(buf)
  H.ok(s and s.inner, "inner selection")
  local orig = vim.ui.select
  local shown
  vim.ui.select = function(items, opts, cb)
    shown = {}
    for i, it in ipairs(items) do shown[i] = opts.format_item(it) end
    cb(items[3]) -- bottom 5'->3'
  end
  vim.cmd("DnaYank")
  vim.ui.select = orig
  H.eq(4, #shown)
  H.ok(shown[3]:find("feature's strand", 1, true), shown[3])
  H.ok(not shown[1]:find("feature's strand", 1, true), shown[1])
  -- bottom 5'->3' of a reverse CDS's spliced residues is its coding sequence
  H.eq(rec:extract(bla):upper(), vim.fn.getreg('"'):upper())
  H.eq("ATG", vim.fn.getreg('"'):sub(1, 3):upper())
end)

case("visual_gy_picker", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  buffer.goto_bp(0, 101)
  local orig = vim.ui.select
  vim.ui.select = function(items, _, cb) cb(items[4]) end -- bottom 3'->5'
  keys("v9lgy")
  vim.ui.select = orig
  H.eq(comp(rec.seq:sub(101, 110)), vim.fn.getreg('"'):upper())
  H.eq(101, sel.get(buf).start)
end)

return T
