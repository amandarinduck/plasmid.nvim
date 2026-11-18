local H = require("helpers")
local pio = require("plasmid.io")
local P = require("plasmid.bio.primer")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local primers = require("plasmid.primers")

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
  return buf, st, buffer.record(st), dst
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("snapgene_annealing_calls", function()
  -- SnapGene stores where each primer anneals; the binding rule reproduces it
  local agree, total, differ = 0, 0, {}
  for _, f in ipairs(vim.fn.glob(H.fixture("snapgene") .. "/*.dna", false, true)) do
    local rec = pio.read(f).records[1]
    for _, ft in ipairs(rec.features) do
      local sg = ft.snapgene_site
      if sg and ft.primer_seq then
        total = total + 1
        local hit
        for _, s in ipairs(P.bind(rec.seq:tostring(), rec:is_circular(), ft.primer_seq)) do
          if s.strand == sg.strand and s.start == sg.a and s.stop == sg.b and s.anneal == #sg.annealed then hit = s end
        end
        if hit then agree = agree + 1 else differ[#differ + 1] = ft:get("label") end
      end
    end
  end
  H.log(("    [primer] %d/%d SnapGene primer sites reproduced (differs: %s)"):format(agree, total, table.concat(differ, ", ")))
  H.eq(20, total)
  -- o7407's stored annealed bases do not match the template (edited after the primer was placed)
  H.eq({ "o7407" }, differ)
end)

case("constructed_sites_property", function()
  math.randomseed(11)
  local B = { "A", "C", "G", "T" }
  local function rnd(n) local t = {} for i = 1, n do t[i] = B[math.random(4)] end return table.concat(t) end
  local function other(x) local y repeat y = B[math.random(4)] until y ~= x return y end
  local checked = 0
  for trial = 1, 600 do
    local circular = trial % 2 == 0
    local L = math.random(300, 1500)
    local tmpl = rnd(L)
    local strand = (trial % 4 < 2) and 1 or -1
    local anneal = math.random(12, 30)
    local tail_len = math.random(0, 8)
    -- annealed region start (may cross the origin on circular templates)
    local s = circular and math.random(1, L) or math.random(tail_len + 2, L - anneal - tail_len - 1)
    local function tpos(i) return (i - 1) % L + 1 end
    local function tb(i) local q = tpos(i) return tmpl:sub(q, q) end
    -- primer bases in primer orientation (5'->3'), k = 1 at the 3' end
    local pr = {}
    local mism = {}
    -- internal mismatches: at least 10 bases from the 3' end and 5 matches after (towards 5')
    local nm = math.random(0, 2)
    local cand = {}
    for k = 11, anneal - 5 do cand[#cand + 1] = k end
    local used = {}
    for _ = 1, nm do
      if #cand == 0 then break end
      local k = table.remove(cand, math.random(#cand))
      local okk = true
      for u in pairs(used) do if math.abs(u - k) < 6 then okk = false end end
      if okk then used[k] = true end
    end
    local n = anneal + tail_len
    for k = 1, n do
      local pos = strand == 1 and (s + anneal - 1 - (k - 1)) or (s + (k - 1))
      local x = tb(pos)
      local want = strand == 1 and x or P.revcomp(x)
      local base
      if k > anneal then base = other(want)          -- tail never pairs
      elseif used[k] then base = other(want)
      else base = want end
      pr[n - k + 1] = base
    end
    for k in pairs(used) do mism[#mism + 1] = n - k + 1 end
    table.sort(mism)
    local primer = table.concat(pr)
    local stop = tpos(s + anneal - 1)
    local found
    for _, site in ipairs(P.bind(tmpl, circular, primer)) do
      if site.strand == strand and site.start == tpos(s) and site.stop == stop then found = site end
    end
    if not found then
      error(("trial %d (%s, strand %d, s %d, anneal %d, tail %d): site not found"):format(trial,
        circular and "circular" or "linear", strand, s, anneal, tail_len))
    end
    H.eq(anneal, found.anneal, "anneal, trial " .. trial)
    H.eq(primer:sub(1, tail_len), found.tail, "tail, trial " .. trial)
    H.eq(mism, found.mismatches, "mismatches, trial " .. trial)
    checked = checked + 1
  end
  H.log(("    [primer] %d constructed sites (tails, internal mismatches, origin) found exactly"):format(checked))
end)

case("binding_rules", function()
  local tmpl = ("ACGT"):rep(5) .. "GATTACAGGCTTCCAGGTTACGATCGGAAC" .. ("TGCA"):rep(5)
  local site = tmpl:sub(21, 50) -- 30 nt
  -- a mismatch in the 3' anchor: no binding
  local bad3 = site:sub(1, 29) .. (site:sub(30, 30) == "A" and "C" or "A")
  H.eq(0, #P.bind(tmpl, false, bad3))
  -- a codon change (3 mismatches) with >= 5 matches beyond it is annealed through
  local mut = site:sub(1, 8) .. "CCC" .. site:sub(12)
  local s = assert(P.bind(tmpl, false, mut)[1])
  H.eq(30, s.anneal)
  -- a mismatch 3 bases from the 5' end is not: those bases become a tail
  local near5 = site:sub(1, 2) .. (site:sub(3, 3) == "A" and "C" or "A") .. site:sub(4)
  s = assert(P.bind(tmpl, false, near5)[1])
  H.eq(27, s.anneal); H.eq(near5:sub(1, 3), s.tail)
  -- a primer longer than a linear template end: the overhang is a tail
  local over = "GGGG" .. tmpl:sub(1, 20)
  s = assert(P.bind(tmpl, false, over)[1])
  H.eq(20, s.anneal); H.eq("GGGG", s.tail)
end)

case("add_display_save_and_reload", function()
  local buf, st, rec, path = open("ncbi/L09137.gb")
  local n0 = #rec.features
  vim.cmd("DnaPrimerAdd GTAAAACGACGGCCAGT M13fwd")
  vim.cmd("DnaPrimerAdd ccggaattcCAGGAAACAGCTATGAC M13rev_Eco")
  H.eq(n0 + 2, #rec.features)
  local list = primers.list(rec)
  H.eq(2, #list)
  local fwd, rev = list[1].strand == 1 and list[1] or list[2], list[1].strand == -1 and list[1] or list[2]
  H.eq("M13fwd", fwd.name); H.eq("M13rev_Eco", rev.name)
  H.eq("CCGGAATTC", rev.tail)
  H.ok(rev.f.location_text:find("^complement"), rev.f.location_text)
  -- the forward primer's bases are drawn over its site, label at the 3' end
  local lay = st.layout
  local found_fwd, found_rev = false, false
  for row0, info in pairs(lay.rows) do
    if info.kind == "primer" then
      local line = lay.lines[row0 + 1]
      if info.strand == 1 and line:find("›M13fwd", 1, true) then
        found_fwd = true
        -- every primer base sits in the column of the template base it copies
        for i = 1, #fwd.seq do
          local pos = fwd.stop - (#fwd.seq - i)
          if pos >= info.bs and pos <= info.be then
            local _, col = lay:pos_of(pos)
            H.eq(fwd.seq:sub(i, i), vim.fn.strcharpart(line, col, 1), "fwd base " .. i)
            H.eq(fwd.f.id, lay:feature_at(row0, col), "feature under primer base")
          end
        end
      end
      if info.strand == -1 and line:find("M13rev_Eco‹", 1, true) then found_rev = true end
    end
  end
  H.ok(found_fwd, "forward primer row"); H.ok(found_rev, "reverse primer row")
  -- undo/redo
  vim.cmd("DnaUndo")
  H.eq(1, #primers.list(rec))
  vim.cmd("DnaRedo")
  H.eq(2, #primers.list(rec))
  -- written as primer_bind with the full sequence in a note; read back the same
  vim.cmd("silent write")
  local back = pio.read(path).records[1]
  local got = {}
  for _, p in ipairs(primers.list(back)) do got[p.name] = p.seq end
  H.eq("GTAAAACGACGGCCAGT", got.M13fwd)
  H.eq("CCGGAATTCCAGGAAACAGCTATGAC", got.M13rev_Eco)
end)

case("visual_gp_and_errors", function()
  local buf, st, rec = open("ncbi/L09137.gb")
  local orig_sel, orig_in = vim.ui.select, vim.ui.input
  vim.ui.select = function(items, _, cb) cb(items[2]) end   -- "rev"
  vim.ui.input = function(_, cb) cb("myRev") end
  buffer.goto_bp(0, 500)
  keys("v")
  buffer.goto_bp(0, 519)
  keys("gp")
  vim.ui.select, vim.ui.input = orig_sel, orig_in
  local p = primers.list(rec)[1]
  H.eq("myRev", p.name); H.eq(-1, p.strand)
  H.eq(P.revcomp(rec.seq:sub(500, 519)), p.seq)
  local msgs = {}
  local on = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  vim.cmd("DnaPrimerAdd AAAAAAAAAAAAAAAAAAAACCCCCCCCCCGT nope")
  vim.cmd("DnaPrimerAdd ACGTXYZ bad")
  vim.notify = on
  H.ok(msgs[1]:find("does not bind", 1, true), msgs[1])
  H.ok(msgs[2]:find("not a nucleotide", 1, true), msgs[2])
end)

case("pcr_products_match_pydna", function()
  local pcr = require("plasmid.pcr")
  local data = vim.json.decode(H.read("tests/oracle/expected/pcr.json"))
  local recs = {}
  local same, none_ok, extra = 0, 0, 0
  for i, c in ipairs(data.cases) do
    recs[c.file] = recs[c.file] or pio.read(H.fixture("ncbi/" .. c.file .. ".gb")).records[1]
    local tmpl = recs[c.file].seq:tostring():upper()
    local prods = pcr.products(tmpl, c.topology == "circular", c.fwd, c.rev)
    if c.product == vim.NIL or c.product == nil then
      if #prods ~= 0 then error(("case %d: pydna finds no product, got %d"):format(i, #prods)) end
      none_ok = none_ok + 1
    else
      local seqs = {}
      for _, p in ipairs(prods) do seqs[pcr.product_sequence(tmpl, p)] = true end
      if not seqs[c.product] then
        error(("case %d (%s %s): pydna product %d bp not among %d products"):format(i, c.file, c.topology, #c.product, #prods))
      end
      if #prods > 1 then extra = extra + 1 end
      same = same + 1
    end
  end
  H.log(("    [pcr] %d products identical to pydna, %d no-product cases agree, %d cases with extra weak products"):format(
    same, none_ok, extra))
  H.eq(#data.cases, same + none_ok)
end)

case("pcr_command_builds_record", function()
  local buf, st, rec = open("ncbi/L09137.gb")
  vim.cmd("DnaPCR GTAAAACGACGGCCAGT ccggaattcCAGGAAACAGCTATGAC")
  local nb = vim.api.nvim_get_current_buf()
  H.ok(nb ~= buf, "product opened")
  local prod = buffer.record(buffer.get(nb))
  H.eq(112, prod:len()); H.eq("linear", prod.topology)
  local s = prod.seq:tostring():upper()
  H.eq("GTAAAACGACGGCCAGT", s:sub(1, 17))
  H.eq(P.revcomp("CCGGAATTCCAGGAAACAGCTATGAC"), s:sub(-26))
  local names = {}
  for _, p in ipairs(primers.list(prod)) do names[#names + 1] = p.name .. ":" .. p.start .. ".." .. p.stop end
  table.sort(names)
  H.eq({ "GTAAAACGACGG:1..17", "ccggaattcCAG:87..112" }, names)
  H.ok(buffer.get(nb).modified)
  vim.cmd("bwipeout! " .. nb)
  -- named primers from the record; no product when they face away
  vim.cmd("DnaPrimerAdd GTAAAACGACGGCCAGT M13fwd")
  vim.cmd("DnaPrimerAdd CAGGAAACAGCTATGAC M13rev")
  vim.cmd("DnaPCR M13fwd M13rev")
  nb = vim.api.nvim_get_current_buf()
  H.eq(103, buffer.record(buffer.get(nb)):len())
  vim.cmd("bwipeout! " .. nb)
end)

return T
