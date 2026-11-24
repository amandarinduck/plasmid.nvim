-- A windowed linear layout must be the eager one, built a piece at a time:
-- same number of rows, same text, same row kinds, same answers to position
-- queries. These run over the whole fixture corpus with windowing forced on,
-- so an ordinary plasmid exercises the same code a genome does.
local H = require("helpers")

require("plasmid").setup({})
local io_ = require("plasmid.io")
local config = require("plasmid.config")
local linear = require("plasmid.view.linear")
local buffer = require("plasmid.buffer")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function opts(G, windowed)
  local o = vim.deepcopy(config.options)
  o.linear.group = G
  o.linear.window = { enabled = windowed, margin = 0, cache_blocks = 1e9 }
  return o
end

--- Build both layouts for a record and compare them row by row.
local function compare(rec, G, width, extra, tag)
  local eager = linear.build(rec, opts(G, false), width, extra)
  local win = linear.build(rec, opts(G, true), width, extra)
  H.ok(win.windowed, tag .. ": window mode was not used")
  H.eq(#eager.lines, #win.lines, tag .. ": line count")
  -- nothing is laid out until it is asked for (a record with no residues
  -- still gets its "no sequence in this file" note after the last block)
  local last = win.offsets[win.nblocks] + win:plan_of(win.nblocks - 1).n
  for i = win.header_rows + 1, last do
    if win.lines[i] ~= "" then error(tag .. ": row " .. i .. " was built before it was needed") end
  end
  for b = 0, win.nblocks - 1 do win:build_block(b) end
  for i = 1, #eager.lines do
    if eager.lines[i] ~= win.lines[i] then
      error(("%s: line %d\n  eager: %q\n  window: %q"):format(tag, i, eager.lines[i], win.lines[i]))
    end
  end
  for r = 0, #eager.lines - 1 do
    local a, c = eager.rows[r], win.rows[r]
    H.eq(a and a.kind, c and c.kind, ("%s: row %d kind"):format(tag, r))
    if a and a.bs then
      H.eq({ a.bs, a.be }, { c.bs, c.be }, ("%s: row %d span"):format(tag, r))
    end
  end
  return eager, win
end

local function fixtures()
  local out = {}
  for _, pat in ipairs({ "**/*.gb", "**/*.gbk", "**/*.fa", "**/*.dna" }) do
    for _, p in ipairs(H.fixtures(pat)) do out[#out + 1] = p end
  end
  table.sort(out)
  return out
end

case("windowed_layout_matches_the_eager_one", function()
  local files, records, rows = 0, 0, 0
  for _, path in ipairs(fixtures()) do
    local res = io_.read(path)
    if res then
      files = files + 1
      for ri, rec in ipairs(res.records) do
        records = records + 1
        for _, G in ipairs({ 10, 0 }) do
          local tag = ("%s r%d G=%d"):format(vim.fn.fnamemodify(path, ":t"), ri, G)
          local eager = compare(rec, G, 100, {}, tag)
          rows = rows + #eager.lines
        end
      end
    end
  end
  H.ok(files > 40, "expected the fixture corpus, got " .. files .. " files")
  H.log(("    [window] %d files, %d records, %d rows identical to the eager layout"):format(files, records, rows))
end)

-- The overlays are the part most likely to drift: enzyme label lanes and
-- primer lanes are assigned per line, and the plan has to count their rows
-- before any of them is drawn.
case("overlays_match_too", function()
  local checked = 0
  for _, path in ipairs(fixtures()) do
    local res = io_.read(path)
    for _, rec in ipairs(res and res.records or {}) do
      if not rec:is_protein() and rec:len() > 0 and rec:len() < 60000 then
        local extra = {}
        local entries = require("plasmid.enzymes").compute(rec, { mode = "unique" })
        if #entries > 0 then extra.enzymes = entries end
        local primers = require("plasmid.primers").list(rec)
        if #primers > 0 then
          extra.primers = primers
          extra.primer_cells = function(p) return require("plasmid.primers").cells(rec, p) end
        end
        if extra.enzymes or extra.primers then
          checked = checked + 1
          compare(rec, 10, 100, extra, "overlays " .. vim.fn.fnamemodify(path, ":t"))
          compare(rec, 0, 86, extra, "overlays " .. vim.fn.fnamemodify(path, ":t") .. " G=0")
        end
      end
    end
  end
  H.ok(checked >= 10, "expected records with enzymes or primers, got " .. checked)
  H.log(("    [window] %d records with enzyme and primer overlays identical"):format(checked))
end)

case("positions_answer_without_building", function()
  local res = assert(io_.read(H.fixture("ncbi/NC_005816.gb")))
  local rec = res.records[1]
  local eager = linear.build(rec, opts(10, false), 100, {})
  local win = linear.build(rec, opts(10, true), 100, {})
  for bp = 1, rec:len() do
    local r1, c1 = eager:pos_of(bp)
    local r2, c2 = win:pos_of(bp)
    H.eq({ r1, c1 }, { r2, c2 }, "pos_of " .. bp)
    H.eq(bp, win:bp_at(r2, c2), "bp_at " .. bp)
  end
  -- still nothing built: positions come from the row plan, not the text
  local built = 0
  for _ in pairs(win._built) do built = built + 1 end
  H.eq(0, built, "pos_of built blocks it did not need")
  -- every row maps to the block that covers it
  for row = win.header_rows, #win.lines - 1 do
    local b = win:block_at_row(row)
    H.ok(b ~= nil, "no block for row " .. row)
    local off = win.offsets[b + 1]
    H.ok(row >= off, ("row %d before its block's first row %d"):format(row, off))
    H.ok(row < off + win:plan_of(b).n, ("row %d past its block's last row"):format(row))
  end
end)

case("a_feature_resolves_on_an_unbuilt_row", function()
  local res = assert(io_.read(H.fixture("ncbi/J01749.gb")))
  local rec = res.records[1]
  local eager = linear.build(rec, opts(10, false), 100, {})
  local win = linear.build(rec, opts(10, true), 100, {})
  local found = 0
  for row = 0, #eager.lines - 1 do
    local info = eager.rows[row]
    if info and info.kind == "track" and info.segs and #info.segs > 0 then
      local col = info.segs[1][1]
      H.eq(eager:feature_at(row, col), win:feature_at(row, col), "feature_at row " .. row)
      found = found + 1
    end
  end
  H.ok(found > 10, "expected feature rows to test, got " .. found)
end)

case("scrolling_keeps_a_bounded_number_of_rows", function()
  local res = assert(io_.read(H.fixture("ncbi/NC_001133.gb")))
  local rec = res.records[1]
  local o = opts(10, true)
  o.linear.window = { enabled = true, margin = 1, cache_blocks = 64 }
  local win = linear.build(rec, o, 100, {})
  local page = 40
  local row = win.header_rows
  while row + page < #win.lines do
    win:ensure(row, row + page)
    row = row + page
  end
  local built = 0
  for _ in pairs(win._built) do built = built + 1 end
  H.ok(built <= 256, "scrolled the record and kept " .. built .. " blocks")
  local nonblank = 0
  for i = win.header_rows + 1, #win.lines do
    if win.lines[i] ~= "" then nonblank = nonblank + 1 end
  end
  H.ok(nonblank < (#win.lines - win.header_rows) / 4,
    ("%d of %d rows still laid out after scrolling through"):format(nonblank, #win.lines))
  H.log(("    [window] 230 kb scrolled end to end: %d blocks and %d rows retained of %d")
    :format(built, nonblank, #win.lines))
end)

-- The threshold is what decides the mode, and ordinary plasmids must keep the
-- eager path.
case("threshold_decides_the_mode", function()
  local small = assert(io_.read(H.fixture("ncbi/J01749.gb"))).records[1]
  local o = vim.deepcopy(config.options)
  H.ok(not linear.build(small, o, 100, {}).windowed, "a 4.4 kb plasmid should not be windowed")
  o.linear.window = { enabled = "auto", threshold = 1000 }
  H.ok(linear.build(small, o, 100, {}).windowed, "threshold 1000 should window a 4.4 kb plasmid")
  o.linear.window = { enabled = false, threshold = 1000 }
  H.ok(not linear.build(small, o, 100, {}).windowed, "enabled = false should force the eager path")
end)

-- End to end: a windowed buffer behaves like any other one.
case("a_windowed_buffer_scrolls_and_jumps", function()
  local tmp = vim.fn.tempname()
  vim.fn.mkdir(tmp, "p")
  local dst = tmp .. "/NC_001133.gb"
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture("ncbi/NC_001133.gb"))); fh:close()
  local prev = config.options.linear.window
  config.options.linear.window = { enabled = true, margin = 1, cache_blocks = 128 }
  vim.cmd("silent edit " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local win = vim.api.nvim_get_current_win()
  local st = buffer.get(buf)
  H.ok(st and st.layout.windowed, "buffer did not open windowed")
  -- the first screen is there
  local row1 = st.layout:pos_of(1)
  local shown = vim.api.nvim_buf_get_lines(buf, row1, row1 + 1, false)[1]
  H.ok(shown:match("%u%u%u"), "first sequence row not laid out: " .. vim.inspect(shown))
  -- jumping fills where it lands, and the cursor reports the right base
  for _, bp in ipairs({ 100000, 230000, 5, 123456 }) do
    buffer.goto_bp(win, bp)
    buffer.fill(buf, win)
    H.eq(bp, (buffer.cursor_info(win)), "cursor base pair after jumping to " .. bp)
    local row = vim.api.nvim_win_get_cursor(win)[1]
    local line = vim.api.nvim_buf_get_lines(buf, row - 1, row, false)[1]
    H.ok(line ~= "", "row " .. row .. " is blank after jumping to " .. bp)
    H.eq(st.layout.lines[row], line, "buffer text differs from the layout at row " .. row)
  end
  -- the buffer is still the user's: unmodified and not editable by hand
  H.eq(false, vim.bo[buf].modified, "windowed fills marked the buffer modified")
  H.eq(false, vim.bo[buf].modifiable)
  config.options.linear.window = prev
  vim.cmd("silent bwipeout!")
end)

return T
