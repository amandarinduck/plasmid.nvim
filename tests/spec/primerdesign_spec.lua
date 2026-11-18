local H = require("helpers")
local P = require("plasmid.bio.primer")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local design = require("plasmid.primerdesign")
local pcr = require("plasmid.pcr")

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

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("candidates_obey_the_constraints", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local tmpl = rec.seq:tostring()
  local o = design.defaults
  local checked = 0
  for _, region in ipairs({ { 86, 1276 }, { 1915, 2106 }, { 3293, 4153 } }) do
    local sel = { start = region[1], stop = region[2] }
    local list = assert(design.candidates(rec, sel))
    H.ok(#list > 0 and #list <= o.pairs, "candidate count " .. #list)
    for i, p in ipairs(list) do
      for _, x in ipairs({ p.fwd, p.rev }) do
        H.ok(x.len >= o.len[1] and x.len <= o.len[2], "length " .. x.len)
        H.ok(x.tm >= o.tm[1] - 3 and x.tm <= o.tm[2] + 3, ("Tm %.1f"):format(x.tm))
        H.eq(1, x.sites, x.seq .. " binds once")
        H.ok(x.shift >= 0 and x.shift <= o.slack, "shift " .. x.shift)
      end
      H.ok(p.dtm <= o.max_dtm, ("ΔTm %.1f"):format(p.dtm))
      if i > 1 then H.ok(list[i - 1].score <= p.score, "sorted by score") end
      -- the pair really amplifies the selected region
      local prods = pcr.products(tmpl, rec:is_circular(), p.fwd.seq, p.rev.seq)
      H.eq(1, #prods, "one product for " .. p.fwd.seq .. " / " .. p.rev.seq)
      H.eq(p.product, prods[1].length, "product length")
      H.eq(p.fwd.start, prods[1].start); H.eq(p.rev.stop, prods[1].stop)
      -- the selection lies inside the product
      local L = rec:len()
      H.ok((sel.start - prods[1].start) % L <= o.slack, "5' end covers the selection")
      H.ok((prods[1].stop - sel.stop) % L <= o.slack, "3' end covers the selection")
      checked = checked + 1
    end
  end
  H.log(("    [design] %d candidate pairs: constraints met and each amplifies its region"):format(checked))
end)

case("refusals_and_options", function()
  local buf, st, rec = open("ncbi/L09137.gb")
  local list, err = design.candidates(rec, { start = 100, stop = 120 })
  H.eq(nil, list); H.ok(err:find("only 21 bp", 1, true), err)
  -- an impossible Tm window is reported with the range that side can reach
  list, err = design.candidates(rec, { start = 100, stop = 900 }, { tm = { 20, 25 } })
  H.eq(nil, list)
  H.ok(err:find("no usable forward primer", 1, true) and err:find("that side reaches", 1, true), err)
  -- pBR322 4300..900: the two sides cannot meet within 3 °C, and the message says so
  local pbr = require("plasmid.io").read(H.fixture("ncbi/J01749.gb")).records[1]
  list, err = design.candidates(pbr, { start = 4300, stop = 900 })
  H.eq(nil, list); H.ok(err:find("closest", 1, true), err)
  -- more slack solves it
  list = assert(design.candidates(pbr, { start = 4300, stop = 900 }, { slack = 40 }))
  H.ok(#list > 0)
  -- a tighter length range is honoured
  list = assert(design.candidates(pbr, { start = 86, stop = 1276 }, { len = { 20, 24 } }))
  for _, p in ipairs(list) do
    H.ok(p.fwd.len >= 20 and p.fwd.len <= 24 and p.rev.len >= 20 and p.rev.len <= 24, "len range")
  end
  -- target Tm shifts the ranking
  local hot = assert(design.candidates(pbr, { start = 86, stop = 1276 }, { tm = { 55, 70 }, target_tm = 68 }))
  local cool = assert(design.candidates(pbr, { start = 86, stop = 1276 }, { tm = { 55, 70 }, target_tm = 57 }))
  H.ok(hot[1].fwd.tm > cool[1].fwd.tm, ("%.1f > %.1f"):format(hot[1].fwd.tm, cool[1].fwd.tm))
end)

case("accept_adds_a_pair_with_tails", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  vim.cmd("DnaSelect 86..1276")
  local picked
  local orig = vim.ui.select
  vim.ui.select = function(items, _, cb) picked = items[1]; cb(items[1]) end
  vim.cmd("DnaPrimerDesign name=tet fwd=EcoRI rev=BamHI")
  vim.ui.select = orig
  H.ok(picked, "a pair was offered")
  local list = require("plasmid.primers").list(rec)
  H.eq(2, #list)
  local f, r = list[1].strand == 1 and list[1] or list[2], list[1].strand == -1 and list[1] or list[2]
  H.eq("tet_F", f.name); H.eq("tet_R", r.name)
  H.eq("GAATTC", f.tail:upper()); H.eq("GGATCC", r.tail:upper())
  H.eq(picked.fwd.seq, f.seq:sub(7))
  H.eq(picked.rev.seq, r.seq:sub(7))
  H.eq(picked.fwd.start, f.start); H.eq(picked.rev.stop, r.stop)
  -- the tailed pair amplifies a product carrying both sites
  local prod = pcr.products(rec.seq:tostring(), true, f.seq, r.seq)[1]
  H.eq(picked.product + 12, prod.length)
  local seq = pcr.product_sequence(rec.seq:tostring(), prod):upper()
  H.eq("GAATTC", seq:sub(1, 6))
  H.eq(P.revcomp("GGATCC"), seq:sub(-6))
  -- the whole tet CDS is inside the product
  H.ok(seq:find(rec.seq:sub(86, 1276):upper(), 1, true), "tet amplified in full")
  -- one undo removes both primers
  vim.cmd("DnaUndo")
  H.eq(0, #require("plasmid.primers").list(rec))
end)

return T
