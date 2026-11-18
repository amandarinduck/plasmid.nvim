-- Integration tests inside a real (headless) Neovim: autocmds, buffer state,
-- cursor <-> bp mapping in both views, writing, conversion.
local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")

local T = { _order = {} }
local function case(name, fn)
  T[name] = fn
  table.insert(T._order, name)
end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")

local function copy_fixture(rel, name)
  local dst = tmp .. "/" .. (name or vim.fn.fnamemodify(rel, ":t"))
  local fh = assert(io.open(dst, "wb"))
  fh:write(H.read(H.fixture(rel)))
  fh:close()
  return dst
end

local function open(path)
  vim.cmd("silent edit " .. vim.fn.fnameescape(path))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.get(buf)
end

vim.o.columns, vim.o.lines = 120, 50

case("bufreadcmd_opens_projection", function()
  local buf, st = open(copy_fixture("ncbi/J01749.gb"))
  H.ok(st, "no state")
  H.eq("plasmid", vim.bo[buf].filetype)
  H.eq(false, vim.bo[buf].modifiable)
  H.eq("acwrite", vim.bo[buf].buftype)
  H.eq(false, vim.bo[buf].modified)
  H.eq("linear", st.view)
  H.ok(vim.api.nvim_buf_line_count(buf) > 100)
  H.ok(vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1]:find("SYNPBR322", 1, true))
end)

case("nomodifiable_is_enforced", function()
  local ok = pcall(vim.cmd, "normal! x")
  H.ok(not ok or not vim.bo.modified, "projection buffer was edited")
end)

case("linear_cursor_roundtrip_every_bp", function()
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  local rec = buffer.record(st)
  local win = vim.api.nvim_get_current_win()
  local seq = rec.seq:tostring():upper()
  for bp = 1, rec:len() do
    buffer.goto_bp(win, bp)
    local got = buffer.cursor_info(win)
    if got ~= bp then error(("goto %d landed on %s"):format(bp, tostring(got))) end
    local pos = vim.api.nvim_win_get_cursor(win)
    local ch = vim.api.nvim_buf_get_lines(buf, pos[1] - 1, pos[1], false)[1]:sub(pos[2] + 1, pos[2] + 1)
    if ch ~= seq:sub(bp, bp) then error(("bp %d shows %q, sequence has %q"):format(bp, ch, seq:sub(bp, bp))) end
  end
end)

case("linear_feature_under_cursor_covers_bp", function()
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  local rec = buffer.record(st)
  local checked = 0
  for row0, info in pairs(st.layout.rows) do
    if info.kind == "track" then
      local line = st.layout.lines[row0 + 1]
      for _, seg in ipairs(info.segs) do
        local vcol = seg[1]
        local bp = st.layout:bp_at(row0, vcol)
        local fid = st.layout:feature_at(row0, vcol)
        H.eq(seg[3], fid)
        local f = rec:feature_by_id(fid)
        local inside = false
        for _, e in ipairs(rec:extents(f)) do
          if bp >= e[1] - 1 and bp <= e[2] + 1 then inside = true end
        end
        H.ok(inside, ("feature %s drawn at bp %d outside its extents"):format(f.location_text, bp))
        H.ok(vim.fn.byteidx(line, vcol) >= 0)
        checked = checked + 1
      end
    end
  end
  H.ok(checked > 40, "few segments checked: " .. checked)
end)

case("highlight_ranges_are_valid", function()
  local st = buffer.get(0)
  for _, view in ipairs({ "linear", "circular" }) do
    buffer.set_view(vim.api.nvim_get_current_buf(), view)
    for row0 = 0, #st.layout.lines - 1 do
      local len = #st.layout.lines[row0 + 1]
      for _, h in ipairs(st.layout:row_highlights(row0)) do
        H.ok(h[1] >= 0 and (h[2] == -1 or h[2] >= h[1]), ("bad range row %d"):format(row0))
        H.ok(h[1] <= len + 3 * 200, "range start beyond line")
        H.ok(vim.fn.hlexists(h[3]) == 1, "undefined group " .. h[3])
      end
    end
  end
  buffer.set_view(vim.api.nvim_get_current_buf(), "linear")
end)

case("circular_view_and_angle_mapping", function()
  local buf = vim.api.nvim_get_current_buf()
  local win = vim.api.nvim_get_current_win()
  buffer.goto_bp(win, 1000)
  buffer.toggle_view(buf)
  local st = buffer.get(buf)
  H.eq("circular", st.view)
  local L = buffer.record(st):len()
  -- one braille cell spans 2*pi*R / (cells around the ring); allow 2 cells of error
  local tol = math.ceil(L / (2 * math.pi * st.layout.R / 2) * 2) + 2
  for _, bp in ipairs({ 1, 250, 1000, 2180, 3333, L }) do
    buffer.goto_bp(win, bp)
    local got = buffer.cursor_info(win)
    local d = math.abs(got - bp)
    d = math.min(d, L - d)
    H.ok(d <= tol, ("circular bp %d -> %d (tol %d)"):format(bp, got, tol))
  end
  -- toggling back preserves position within tolerance
  buffer.goto_bp(win, 1000)
  buffer.toggle_view(buf)
  local back = buffer.cursor_info(win)
  H.ok(math.abs(back - 1000) <= tol, "toggle back landed on " .. back)
end)

case("circular_feature_lookup_from_label", function()
  local buf = vim.api.nvim_get_current_buf()
  buffer.set_view(buf, "circular")
  local st = buffer.get(buf)
  local rec = buffer.record(st)
  local found = false
  for r, line in ipairs(st.layout.lines) do
    local s = line:find("ROP protein", 1, true)
    if s then
      local vcol = vim.fn.charidx(line, s - 1)
      local fid = st.layout:feature_at(r - 1, vcol)
      H.ok(fid, "no feature at label")
      H.ok(rec:feature_by_id(fid):label():find("ROP"), "label maps to wrong feature")
      found = true
    end
  end
  H.ok(found, "ROP label not drawn")
  buffer.set_view(buf, "linear")
end)

case("write_unmodified_is_byte_identical", function()
  local buf = vim.api.nvim_get_current_buf()
  local path = buffer.get(buf).path
  local before = H.read(path)
  vim.cmd("silent write!")
  H.eq(before, H.read(path))
  H.eq(false, vim.bo[buf].modified)
end)

case("write_as_fasta_converts", function()
  local out = tmp .. "/converted.fasta"
  vim.cmd("silent write " .. vim.fn.fnameescape(out))
  local text = H.read(out)
  H.ok(text:match("^>J01749%.1 Cloning vector pBR322, complete sequence\n"), text:sub(1, 80))
  local res = assert(require("plasmid.io").read(out))
  H.eq(4361, res.records[1]:len())
  -- the original buffer still points at the GenBank file
  H.ok(buffer.get(0).path:match("J01749%.gb$"))
end)

case("multi_record_fasta_and_topology", function()
  local buf, st = open(copy_fixture("ncbi/multi.fasta"))
  H.eq(3, #st.records)
  vim.cmd("DnaRecord 2")
  H.eq(2, st.current)
  H.ok(vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1]:find("J01749.1", 1, true))
  H.eq("linear", buffer.record(st).topology)
  vim.cmd("DnaTopology circular")
  H.eq(true, vim.bo[buf].modified)
  vim.cmd("DnaView circular")
  H.eq("circular", st.view)
  vim.cmd("silent write")
  local res = assert(require("plasmid.io").read(st.path))
  H.eq("circular", res.records[2].topology)
  H.eq("linear", res.records[1].topology)
  H.eq(false, vim.bo[buf].modified)
end)

case("linear_record_refuses_circular_without_bang", function()
  local buf, st = open(copy_fixture("ncbi/AB000263.gb"))
  vim.cmd("DnaView circular")
  H.eq("linear", st.view)
  vim.cmd("DnaView! circular")
  H.eq("circular", st.view)
end)

case("wrap_feature_drawn_on_both_sides_of_origin", function()
  local buf, st = open(copy_fixture("synthetic/grammar.gb"))
  local rec = buffer.record(st)
  local wrap
  for _, f in ipairs(rec.features) do if f:label() == "wrapCDS" then wrap = f end end
  local seen = {}
  for row0, info in pairs(st.layout.rows) do
    for _, s in ipairs(info.segs or {}) do
      if s[3] == wrap.id and info.kind == "track" then
        seen[st.layout:bp_at(row0, s[1]) <= 30 and "head" or "tail"] = true
      end
    end
  end
  H.ok(seen.head and seen.tail, vim.inspect(seen))
  -- the translation across the origin is continuous: 60 nt -> 20 codons
  local map = require("plasmid.view.features").codon_map(rec, wrap)
  H.eq(20, vim.tbl_count(map))
end)

case("new_file_opens_empty", function()
  local buf, st = open(tmp .. "/brand_new.gb")
  H.ok(st, "no state for new file")
  H.eq(0, buffer.record(st):len())
end)

case("snapgene_dna_opens_writes_and_converts", function()
  local path = copy_fixture("snapgene/pFA-KanMX4.dna")
  local before = H.read(path)
  local buf, st = open(path)
  H.eq("snapgene", st.format)
  local rec = buffer.record(st)
  H.eq("circular", rec.topology)
  H.ok(#rec.features > 5)
  vim.cmd("DnaView circular")
  H.eq("circular", st.view)
  vim.cmd("silent write!")
  H.eq(before, H.read(path), "unmodified .dna must be written back byte-identical")
  local out = tmp .. "/pFA.gb"
  vim.cmd("silent write " .. vim.fn.fnameescape(out))
  local res = assert(require("plasmid.io").read(out))
  H.eq(rec:len(), res.records[1]:len())
  H.eq(#rec.features, #res.records[1].features)
  -- :DnaRaw shows a packet table, not binary
  vim.cmd("DnaRaw")
  local raw = vim.api.nvim_buf_get_lines(0, 0, 3, false)
  H.ok(raw[1]:match("SnapGene file"), raw[1])
  vim.cmd("close")
end)

case("snapgene_topology_change_persists", function()
  local path = copy_fixture("snapgene/sample-d.dna", "topo.dna")
  local buf, st = open(path)
  vim.cmd("DnaTopology circular")
  vim.cmd("silent write")
  H.eq("circular", assert(require("plasmid.io").read(path)).records[1].topology)
end)

case("protein_files_linear_only", function()
  for _, rel in ipairs({ "snapgene/sgffp_test.prot", "protein/proteins.faa" }) do
    local buf, st = open(copy_fixture(rel))
    local rec = buffer.record(st)
    H.eq("protein", rec.alphabet, rel)
    H.eq("linear", st.view)
    vim.cmd("DnaView! circular")
    H.eq("linear", st.view, "protein must not switch to circular")
    buffer.toggle_view(buf)
    H.eq("linear", st.view)
    H.ok(vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1]:find(" aa ", 1, true))
    local win = vim.api.nvim_get_current_win()
    buffer.goto_bp(win, 7)
    H.eq(7, (buffer.cursor_info(win)))
    H.ok(require("plasmid.buffer").winbar():find("aa 7", 1, true), require("plasmid.buffer").winbar())
  end
end)

case("track_limit_toggle_keeps_position", function()
  local buf, st = open(copy_fixture("ncbi/J01749.gb", "tracks.gb"))
  local win = vim.api.nvim_get_current_win()
  local function overflow_rows()
    local n = 0
    for _, info in pairs(st.layout.rows) do if info.kind == "overflow" then n = n + 1 end end
    return n
  end
  vim.cmd("DnaTracks 1")
  H.ok(overflow_rows() > 0, "cap 1 should overflow on pBR322")
  buffer.goto_bp(win, 1500)
  vim.cmd("normal +")
  H.eq(0, overflow_rows())
  H.eq(1500, (buffer.cursor_info(win)))
  vim.cmd("normal -")
  H.eq(nil, st.track_limit)
  H.eq(1500, (buffer.cursor_info(win)))
  -- <CR> on an overflow line expands
  vim.cmd("DnaTracks 1")
  for row0, info in pairs(st.layout.rows) do
    if info.kind == "overflow" then
      vim.api.nvim_win_set_cursor(win, { row0 + 1, 0 })
      break
    end
  end
  require("plasmid.ui").activate()
  H.eq("all", st.track_limit)
end)

case("performance_large_records", function()
  local t0 = vim.uv.hrtime()
  local _, st = open(copy_fixture("ncbi/NC_001133.gb"))
  local ms = (vim.uv.hrtime() - t0) / 1e6
  H.ok(st and #st.layout.lines > 3000)
  io.stdout:write(("    [perf] NC_001133 230 kb, %d features: open+render %.0f ms\n"):format(
    #buffer.record(st).features, ms))
  H.ok(ms < 5000, "too slow: " .. ms)

  -- 5 Mb synthetic genome, read-only viewing
  local big = tmp .. "/big.fasta"
  local fh = assert(io.open(big, "wb"))
  fh:write(">big synthetic 5 Mb\n")
  local unit = ("ACGTTGCAAGGCTTAACCGGATATCGCGAT"):rep(2) .. "\n"
  fh:write(unit:rep(math.floor(5e6 / 60)))
  fh:close()
  t0 = vim.uv.hrtime()
  local _, st2 = open(big)
  ms = (vim.uv.hrtime() - t0) / 1e6
  io.stdout:write(("    [perf] 5 Mb FASTA: open+render %.0f ms, %d lines\n"):format(ms, #st2.layout.lines))
  H.ok(ms < 15000, "too slow: " .. ms)
end)

return T
