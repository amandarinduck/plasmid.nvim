local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local config = require("plasmid.config")
local codon = require("plasmid.codon")

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

vim.o.columns, vim.o.lines = 120, 50

case("configured_group_is_used", function()
  local _, st = open("ncbi/J01749.gb")
  H.eq(tonumber(config.options.linear.group) or 10, st.layout.G)
end)

-- For several group sizes: every base pair maps to a column showing that
-- residue, the column maps back to the same base pair, lines fit the window,
-- and line widths stay multiples of 10.
case("every_bp_round_trips_for_each_group_size", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local seq = rec.seq:tostring():upper()
  local width = vim.api.nvim_win_get_width(0)
  for _, G in ipairs({ 10, 0, 3, 5, 7 }) do
    buffer.set_group(buf, G)
    local lay = st.layout
    H.eq(G, lay.G)
    H.eq(0, lay.W % 10, ("G=%d: line width %d"):format(G, lay.W))
    for bp = 1, rec:len() do
      local r, v = lay:pos_of(bp)
      H.eq(seq:sub(bp, bp), lay.lines[r + 1]:sub(v + 1, v + 1):upper(), ("G=%d bp %d drawn"):format(G, bp))
      H.eq(bp, lay:bp_at(r, v), ("G=%d bp %d round trip"):format(G, bp))
    end
    for i, line in ipairs(lay.lines) do
      H.ok(vim.fn.strdisplaywidth(line) <= width, ("G=%d line %d is %d wide"):format(G, i, vim.fn.strdisplaywidth(line)))
    end
    -- the cursor path the plugin itself uses
    for _, bp in ipairs({ 1, 9, 10, 11, 99, 100, 101, 2048, rec:len() }) do
      buffer.goto_bp(0, bp)
      H.eq(bp, (buffer.cursor_info(0)), ("G=%d cursor at %d"):format(G, bp))
    end
  end
  buffer.set_group(buf, nil)
end)

case("continuous_layout_has_no_gaps_and_a_ruler", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  buffer.set_group(buf, 0)
  local lay = st.layout
  local seq = rec.seq:tostring():upper()
  local nrulers = 0
  for row0, info in pairs(lay.rows) do
    if info.kind == "seq" then
      local r, v = lay:pos_of(info.bs)
      H.eq(row0, r)
      H.eq(seq:sub(info.bs, info.be), lay.lines[row0 + 1]:sub(v + 1, v + info.be - info.bs + 1):upper(),
        "residues are contiguous")
      local above = lay.rows[row0 - 1]
      H.eq("ruler", above and above.kind, "ruler above each sequence row")
    elseif info.kind == "ruler" then
      nrulers = nrulers + 1
      local line = lay.lines[row0 + 1]
      -- every multiple of 10 in the line is printed with its last digit over that residue
      for p = math.ceil(info.bs / 10) * 10, info.be, 10 do
        local _, v = lay:pos_of(p)
        local txt = tostring(p)
        H.eq(txt, line:sub(v - #txt + 2, v + 1), ("ruler label %d"):format(p))
      end
    end
  end
  H.eq(math.ceil(rec:len() / lay.W), nrulers)
  buffer.set_group(buf, nil)
end)

-- The point of continuous mode: within a line a codon's bases sit in three
-- adjacent columns, and successive residues of a CDS are exactly 3 apart.
case("continuous_codons_are_contiguous", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  buffer.set_group(buf, 0)
  vim.cmd("DnaTracks all")
  local lay = st.layout
  local win = vim.api.nvim_get_current_win()
  local same_line, gapped, spacing_ok, pairs_seen = 0, 0, 0, 0
  for row0 = 0, #lay.lines - 1 do
    local info = lay.rows[row0]
    if info and info.kind == "aa" then
      local line = lay.lines[row0 + 1]
      local prev = {}
      for vcol = 0, #line - 1 do
        if line:sub(vcol + 1, vcol + 1):match("%S") then
          vim.api.nvim_win_set_cursor(win, { row0 + 1, vcol })
          codon.update(buf, win)
          local c = st.codon
          H.ok(c, "codon at " .. row0 .. ":" .. vcol)
          local rows, cols = {}, {}
          for i, p in ipairs(c.positions) do rows[i], cols[i] = lay:pos_of(p) end
          if rows[1] == rows[2] and rows[2] == rows[3] then
            same_line = same_line + 1
            local lo = math.min(cols[1], cols[3])
            if not (cols[2] == lo + 1 and math.max(cols[1], cols[3]) == lo + 2) then gapped = gapped + 1 end
          end
          local p = prev[c.f.id]
          if p and math.abs(p.n - c.n) == 1 then
            pairs_seen = pairs_seen + 1
            if vcol - p.vcol == 3 then spacing_ok = spacing_ok + 1 end
          end
          prev[c.f.id] = { n = c.n, vcol = vcol }
        end
      end
    end
  end
  H.log(("    [group 0] %d codons within a line, %d gapped; %d/%d neighbouring residues 3 columns apart")
    :format(same_line, gapped, spacing_ok, pairs_seen))
  H.ok(same_line > 500, "codons checked: " .. same_line)
  H.eq(0, gapped, "codons with a gap inside")
  H.eq(pairs_seen, spacing_ok, "neighbouring residues 3 columns apart")
  buffer.set_group(buf, nil)
end)

case("dnagroup_command", function()
  local buf, st = open("ncbi/J01749.gb")
  local conf = tonumber(config.options.linear.group) or 10
  buffer.goto_bp(0, 1234)
  vim.cmd("DnaGroup")
  H.eq(conf == 0 and 10 or 0, st.layout.G, "toggle")
  H.eq(1234, (buffer.cursor_info(0)), "cursor keeps its base pair")
  vim.cmd("DnaGroup")
  H.eq(conf, st.layout.G, "toggle back")
  vim.cmd("DnaGroup 3")
  H.eq(3, st.layout.G)
  H.messages = {}
  vim.cmd("DnaGroup -1")
  local msgs = {}
  for _, m in ipairs(H.messages) do msgs[#msgs + 1] = m.msg end
  H.ok(table.concat(msgs, "\n"):find(":DnaGroup", 1, true), "usage message")
  H.eq(3, st.layout.G, "unchanged after a bad argument")
end)

return T
