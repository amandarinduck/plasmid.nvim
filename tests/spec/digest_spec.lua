local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local R = require("plasmid.bio.restriction")
local D = require("plasmid.digest")
local gel = require("plasmid.gel")

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
local function load(stem) return assert(pio.read(H.fixture("ncbi/" .. stem .. ".gb"))).records[1] end
local function with_topology(rec, topo) return setmetatable({ topology = topo }, { __index = rec }) end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("fragment_lengths_match_oracle_cuts", function()
  -- fragments computed from the Biopython cut lists, independently of the engine
  local expected = vim.json.decode(H.read("tests/oracle/expected/enzymes.json"))
  math.randomseed(7)
  local names = {}
  for _, e in ipairs(require("plasmid.data.enzymes")) do names[#names + 1] = e.name end
  local combos = 0
  for stem, by_topo in pairs(expected) do
    local rec = load(stem)
    local L = rec:len()
    for topo, by_enz in pairs(by_topo) do
      local r = with_topology(rec, topo)
      for _ = 1, 60 do
        local pick = {}
        for _ = 1, math.random(1, 3) do pick[#pick + 1] = names[math.random(#names)] end
        local set = {}
        for _, n in ipairs(pick) do for _, c in ipairs(by_enz[n]) do set[c] = true end end
        local cuts = vim.tbl_keys(set)
        table.sort(cuts)
        local want = {}
        if #cuts == 0 then want = { L }
        elseif topo == "circular" then
          for i = 1, #cuts do
            local nxt = cuts[i % #cuts + 1]
            want[#want + 1] = (#cuts == 1) and L or ((nxt - cuts[i]) % L)
          end
        else
          local prev = 1
          for _, c in ipairs(cuts) do want[#want + 1] = c - prev; prev = c end
          want[#want + 1] = L - prev + 1
        end
        local frags = D.compute(r, pick)
        local got = vim.tbl_map(function(f) return f.length end, frags)
        if not vim.deep_equal(want, got) then
          error(("%s %s %s: want %s got %s"):format(stem, topo, table.concat(pick, "+"), vim.inspect(want), vim.inspect(got)))
        end
        local sum = 0
        for _, g in ipairs(got) do sum = sum + g end
        H.eq(L, sum)
        combos = combos + 1
      end
    end
  end
  H.log(("%d random digests (1-3 enzymes) agree with fragments from the oracle cut lists"):format(combos))
end)

case("overhang_text", function()
  local rec = load("L09137")
  local P = R.prepare(rec)
  local function first(name) return R.sites(P, R.get(name))[1] end
  H.eq("5′ AATT", D.end_text(rec, first("EcoRI")))
  H.eq("3′ TGCA", D.end_text(rec, first("PstI")))
  H.eq("5′ AGCT", D.end_text(rec, first("HindIII")))
  local s = first("BsaI")
  local seq = rec.seq:tostring():upper()
  H.eq("5′ " .. seq:sub(s.cut, s.cut + 3), D.end_text(rec, s), "BsaI overhang = the 4 top-strand bases after the cut")
end)

case("digest_table_select_and_extract", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local db, frags = D.run(buf, { "EcoRI", "HindIII" })
  H.eq(2, #frags)
  H.eq(db, vim.api.nvim_get_current_buf())
  local lines = vim.api.nvim_buf_get_lines(db, 0, -1, false)
  H.ok(lines[1]:find("2 cuts, 2 fragments", 1, true), lines[1])
  local row = D.state[db].first
  -- the 31 bp EcoRI-HindIII fragment across the origin (4360..29)
  local small
  for i, f in ipairs(frags) do if f.length == 31 then small = i end end
  H.ok(small, "31 bp fragment")
  vim.api.nvim_win_set_cursor(0, { row + small - 1, 0 })
  H.ok(vim.api.nvim_get_current_line():find("EcoRI 5′ AATT", 1, true))
  keys("\r")
  H.eq(buf, vim.api.nvim_get_current_buf(), "<CR> returns to the plasmid")
  H.eq(4360, st.selection.start); H.eq(29, st.selection.stop)
  -- extract the large fragment
  vim.api.nvim_set_current_win(vim.fn.bufwinid(db))
  local big = small == 1 and 2 or 1
  vim.api.nvim_win_set_cursor(0, { row + big - 1, 0 })
  local nb, frag = D.extract_under_cursor(db)
  H.ok(nb and nb ~= buf, "new buffer")
  local f = frags[big]
  H.eq(f.length, frag:len()); H.eq("linear", frag.topology)
  H.eq(rec.seq:sub(f.start, f.stop), frag.seq:tostring(), "top strand of the fragment")
  H.ok(frag.description:find("HindIII–EcoRI fragment", 1, true), frag.description)
  -- whole features keep their residues; none extends past the fragment
  local carried, src_whole = 0, 0
  for _, sf in ipairs(rec.features) do
    local lo, hi = sf:bounds()
    if sf.type ~= "source" and lo and #rec:extents(sf) == 1 and lo >= f.start and hi <= f.stop then
      src_whole = src_whole + 1
      local want = rec:extract(sf)
      for _, nf in ipairs(frag.features) do
        if nf.type == sf.type and nf:label() == sf:label() and frag:extract(nf) == want then carried = carried + 1; break end
      end
    end
  end
  H.ok(src_whole > 10, "whole features in the fragment: " .. src_whole)
  H.eq(src_whole, carried, "every whole feature carried with identical residues")
  for _, nf in ipairs(frag.features) do
    local lo, hi = nf:bounds()
    if lo then H.ok(lo >= 1 and hi <= frag:len(), nf:label()) end
  end
  H.ok(buffer.get(nb).modified, "new record is unsaved")
  -- :w writes a GenBank file that reads back identically
  vim.api.nvim_set_current_buf(nb)
  vim.cmd("silent write")
  local back = assert(pio.read(vim.api.nvim_buf_get_name(nb))).records[1]
  H.eq(frag.seq:tostring(), back.seq:tostring())
  H.eq("linear", back.topology)
  H.eq(#frag.features, #back.features)
  H.ok(not buffer.get(nb).modified)
  vim.cmd("bwipeout! " .. nb)
end)

case("stale_digest_and_default_names", function()
  local buf, st = open("ncbi/L09137.gb")
  local msgs = {}
  local orig = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  vim.cmd("DnaDigest")
  H.ok(msgs[#msgs]:find("name the enzymes", 1, true))
  vim.cmd("DnaEnzymes EcoRI HindIII")
  vim.cmd("DnaDigest")
  local db = vim.api.nvim_get_current_buf()
  H.eq({ "EcoRI", "HindIII" }, D.state[db].names)
  -- edit the plasmid: the digest is stale
  vim.api.nvim_set_current_win(vim.fn.bufwinid(buf))
  buffer.goto_bp(0, 100)
  keys("x")
  vim.api.nvim_set_current_win(vim.fn.bufwinid(db))
  vim.api.nvim_win_set_cursor(0, { D.state[db].first, 0 })
  local before = st.selection
  keys("\r")
  vim.notify = orig
  H.ok(msgs[#msgs]:find("changed since this digest", 1, true), msgs[#msgs])
  H.eq(before, st.selection)
  vim.cmd("close")
end)

case("gel_geometry", function()
  local lines, hl, info = gel.render({
    { label = "A", sizes = { 4000, 1000, 1000, 250, 20000 } },
    { label = "B", sizes = { 3000 } },
  }, { ladder = "1kb", height = 20 })
  local rows = {}
  for _, p in ipairs(info.placed[1]) do rows[p[1]] = p[2] end
  H.ok(rows[4000] < rows[1000], "bigger fragments migrate less")
  H.eq({ 250 }, vim.tbl_filter(function(s) return s < 1000 end, info.off[1]))
  H.ok(vim.tbl_contains(info.off[1], 20000))
  local text = table.concat(lines, "\n")
  H.ok(text:find("1000 ×2", 1, true), "co-migrating bands")
  H.ok(text:find("ran off: 250 bp", 1, true)); H.ok(text:find("above the ladder: 20000 bp", 1, true))
  -- log scale: equal ratios give equal row spacing (±1 for rounding)
  local r = info.row_of
  H.ok(math.abs((r(1000) - r(2000)) - (r(2000) - r(4000))) <= 1)
  H.eq("100bp", gel.auto_ladder({ 120, 800 }))
  H.eq("1kb", gel.auto_ladder({ 4000, 300 }))
  -- every highlight range is inside its line
  for _, h in ipairs(hl) do H.ok(h[3] <= #lines[h[1] + 1], "highlight in range") end
end)

case("gel_command_lanes", function()
  local buf = open("ncbi/J01749.gb")
  local b, lanes = require("plasmid.digest").gel(buf, "EcoRI | EcoRI HindIII | uncut ladder=1kb")
  H.eq(3, #lanes)
  H.eq({ 4361 }, lanes[1].sizes)
  H.eq(2, #lanes[2].sizes)
  H.eq("uncut", lanes[3].label)
  local text = table.concat(vim.api.nvim_buf_get_lines(b, 0, -1, false), "\n")
  H.ok(text:find("EcoRI+HindIII", 1, true))
  vim.cmd("close")
end)

return T
