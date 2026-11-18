local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local tm = require("plasmid.bio.tm")
local info = require("plasmid.primerinfo")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent! only")
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if st.view ~= "linear" then buffer.set_view(buf, "linear") end
  return buf, st, buffer.record(st)
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end
local function msgs_of(fn)
  local msgs = {}
  local on = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  fn()
  vim.notify = on
  return msgs
end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("tm_command_and_options", function()
  local buf = open("ncbi/L09137.gb")
  local seq = "ATGTCCCTGCTCTTCTCTCGATGCAA"
  -- the values polymerase-tm's own tests assert for this primer
  H.eq(72, tm.tm(seq, { polymerase = "q5" }))
  H.eq(63, tm.tm(seq, { polymerase = "taq" }))
  H.eq(68, tm.tm(seq, { polymerase = "phusion_hf" }))
  local m = msgs_of(function() vim.cmd("DnaTm " .. seq) end)
  H.ok(m[1]:find("Tm 72 °C", 1, true), m[1])
  H.ok(m[1]:find("Q5 High-Fidelity", 1, true))
  -- options: polymerase, DMSO, explicit salt and concentration
  m = msgs_of(function() vim.cmd("DnaTm " .. seq .. " polymerase=taq") end)
  H.ok(m[1]:find("Tm 63 °C", 1, true), m[1])
  local plain = tm.tm_detail(seq, { polymerase = "q5" })
  local dmso = tm.tm_detail(seq, { polymerase = "q5", dmso = 5 })
  H.ok(math.abs((plain.tm - 3) - dmso.tm) < 1e-9, "5 % DMSO lowers Tm by 3 °C")
  local salted = tm.tm_detail(seq, { polymerase = "q5", salt_mM = 50, conc_nM = 200 })
  H.eq(50, salted.salt_mM); H.eq(200, salted.conc_nM)
  H.ok(salted.tm < plain.tm, "less salt lowers Tm")
  -- two sequences add Ta; errors are reported
  local pair = info.command(buf, { "GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC" })
  H.eq(3, #pair)
  H.ok(pair[3]:find("Ta 57 °C", 1, true), pair[3])
  m = msgs_of(function() vim.cmd("DnaTm ACGTRYKM") end)
  H.ok(m[1]:lower():find("sequence 1", 1, true), m[1])
  m = msgs_of(function() vim.cmd("DnaTm ACGTACGTACGTACGT polymerase=bogus") end)
  H.ok(m[1]:find("unknown polymerase", 1, true), m[1])
  -- no argument: the selection
  vim.cmd("DnaSelect 379..395")
  local out = info.command(buf, {})
  H.ok(out[1]:find("17 nt  GTAAAACGACGGCCAGT", 1, true), out[1])
end)

case("ta_rules_follow_the_polymerase", function()
  local a, b = "GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC"
  local ta_q5, t1, t2 = tm.ta(a, b, { polymerase = "q5" })
  local ta_taq = tm.ta(a, b, { polymerase = "taq" })
  local ta_phu = tm.ta(a, b, { polymerase = "phusion_hf" })
  H.ok(ta_taq < ta_q5, ("taq %d should be below q5 %d"):format(ta_taq, ta_q5))
  H.ok(ta_phu ~= ta_q5, "Phusion uses its own rule")
  -- the cap applies
  local hot = ("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")
  local capped, _, _, cinfo = tm.ta(hot, hot, { polymerase = "q5" })
  H.eq(72, capped); H.ok(cinfo.capped, "cap reported")
  local taq_cap = tm.ta(hot, hot, { polymerase = "taq" })
  H.eq(68, taq_cap)
  -- DMSO lowers Ta by 0.6 °C per %
  local base = tm.ta(a, b, { polymerase = "q5" })
  local with = tm.ta(a, b, { polymerase = "q5", dmso = 5 })
  H.ok(with < base, ("DMSO: %d < %d"):format(with, base))
end)

case("primer_float_and_polymerase_choice", function()
  local buf, st, rec = open("ncbi/L09137.gb")
  vim.cmd("DnaPrimerAdd GTAAAACGACGGCCAGT M13fwd")
  vim.cmd("DnaPrimerAdd ccggaattcCAGGAAACAGCTATGAC M13rev_Eco")
  local p = require("plasmid.primers").list(rec)[1]
  local lines = table.concat(info.primer_lines(rec, p), "\n")
  H.ok(lines:find("M13fwd  forward primer, 17 nt", 1, true), lines)
  H.ok(lines:find("Tm        62 °C", 1, true), lines)
  H.ok(lines:find("binds     1 site", 1, true), lines)
  H.ok(lines:find("pair      with M13rev_Eco: product 112 bp, Ta 57 °C", 1, true), lines)
  -- the tailed primer reports its tail and both Tms
  local q
  for _, x in ipairs(require("plasmid.primers").list(rec)) do if x.name == "M13rev_Eco" then q = x end end
  local qlines = table.concat(info.primer_lines(rec, q), "\n")
  H.ok(qlines:find("5' tail 9 nt", 1, true), qlines)
  H.ok(qlines:find("Tm full", 1, true), qlines)
  -- K on a primer base opens the primer float, not the feature float
  local row, col
  for r, i in pairs(st.layout.rows) do
    if i.kind == "primer" and i.strand == 1 and st.layout.lines[r + 1]:find("M13fwd", 1, true) then
      row = r; col = select(2, st.layout:pos_of(p.stop))
    end
  end
  H.ok(row, "primer row drawn")
  vim.api.nvim_win_set_cursor(0, { row + 1, 0 })
  vim.cmd("normal! 0" .. (col + 1) .. "|")
  local before = #vim.api.nvim_list_wins()
  keys("K")
  H.eq(before + 1, #vim.api.nvim_list_wins(), "float opened")
  local fbuf = vim.api.nvim_win_get_buf(vim.api.nvim_list_wins()[#vim.api.nvim_list_wins()])
  H.ok(table.concat(vim.api.nvim_buf_get_lines(fbuf, 0, -1, false), "\n"):find("M13fwd  forward primer", 1, true))
  vim.cmd("fclose!")
  -- :DnaPolymerase changes the model everywhere; the winbar shows it
  H.ok(buffer.winbar():find("q5", 1, true), buffer.winbar())
  local m = msgs_of(function() vim.cmd("DnaPolymerase taq") end)
  H.ok(m[1]:find("Taq DNA Polymerase", 1, true), m[1])
  H.eq("taq", info.conditions().polymerase)
  H.ok(buffer.winbar():find("taq", 1, true), buffer.winbar())
  H.ok(table.concat(info.primer_lines(rec, p), "\n"):find("Taq DNA Polymerase", 1, true))
  m = msgs_of(function() vim.cmd("DnaPolymerase nope") end)
  H.ok(m[1]:find("unknown polymerase", 1, true), m[1])
  vim.cmd("DnaPolymerase q5")
  H.ok(vim.tbl_contains(info.complete_polymerase("phu"), "phusion_hf"))
end)

return T
