local H = require("helpers")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local su = require("plasmid.searchui")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end
local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function open(rel)
  local dst = tmp .. "/" .. vim.fn.fnamemodify(rel, ":t")
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  vim.cmd("silent edit! " .. vim.fn.fnameescape(dst))
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if st.view ~= "linear" then buffer.set_view(buf, "linear") end
  return buf, st, buffer.record(st)
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end
local function cursor_bp() return (buffer.cursor_info(vim.api.nvim_get_current_win())) end

-- Every highlighted search cell as { [bp] = {kind=seq|comp, current=bool} } (linear view).
local function painted(st)
  local out = {}
  local lay = st.layout
  for row0, info in pairs(lay.rows) do
    if info.kind == "seq" or info.kind == "comp" then
      local line = lay.lines[row0 + 1]
      for _, h in ipairs(lay:row_highlights(row0)) do
        if h[3] == "PlasmidSearch" or h[3] == "PlasmidSearchCurrent" then
          for b = h[1], h[2] - 1 do
            local ch = line:sub(b + 1, b + 1)
            if ch:match("%a") then
              local vcol = vim.fn.strdisplaywidth(line:sub(1, b))
              local bp = lay:bp_at(row0, vcol)
              out[#out + 1] = { bp = bp, kind = info.kind, current = h[3] == "PlasmidSearchCurrent" }
            end
          end
        end
      end
    end
  end
  return out
end

vim.o.columns, vim.o.lines, vim.o.hidden = 120, 50, true

case("find_jumps_highlights_and_reports", function()
  local buf, st = open("ncbi/L09137.gb")
  vim.api.nvim_win_set_cursor(0, { 1, 0 })
  vim.cmd("DnaFind GAATTC")
  H.eq(1, #st.search.hits)
  H.eq(396, cursor_bp(), "cursor on the EcoRI site")
  local cells = painted(st)
  H.eq(6, #cells, "six bases painted")
  for i, c in ipairs(cells) do
    H.ok(c.bp >= 396 and c.bp <= 401 and c.current, "cell " .. i .. " bp " .. tostring(c.bp))
    H.eq("seq", c.kind, "palindromic hit drawn on the top strand")
  end
  H.ok(buffer.winbar():find("/GAATTC 1/1", 1, true), "winbar shows hit index")
  vim.cmd("DnaFind off")
  H.eq(nil, st.search)
  H.eq(0, #painted(st))
  H.ok(not buffer.winbar():find("GAATTC", 1, true))
end)

case("n_and_N_step_through_hits_with_counts", function()
  local buf, st = open("ncbi/J01749.gb")
  vim.api.nvim_win_set_cursor(0, { 1, 0 })
  vim.cmd("DnaFind GATC")
  local hits = st.search.hits
  H.ok(#hits > 10, "pBR322 has many GATC sites, got " .. #hits)
  H.eq(1, st.search.idx)
  keys("n"); H.eq(2, st.search.idx); H.eq(hits[2].start, cursor_bp())
  keys("3n"); H.eq(5, st.search.idx); H.eq(hits[5].start, cursor_bp())
  keys("N"); H.eq(4, st.search.idx)
  keys("4N"); H.eq(#hits, st.search.idx, "N wraps past the first hit")
  keys("n"); H.eq(1, st.search.idx, "n wraps past the last hit")
  H.ok(buffer.winbar():find(("/GATC 1/%d"):format(#hits), 1, true))
  -- ? reverses the meaning of n
  vim.api.nvim_win_set_cursor(0, { 1, 0 })
  vim.cmd("DnaFind! GATC")
  local i = st.search.idx
  keys("n")
  H.eq(i == 1 and #hits or i - 1, st.search.idx, "after ?, n goes backwards")
end)

case("bottom_strand_hits_paint_the_complement_row", function()
  local buf, st = open("ncbi/J01749.gb")
  vim.cmd("DnaFind GGTCTC") -- BsaI, not palindromic
  local by = { [1] = {}, [-1] = {} }
  for _, h in ipairs(st.search.hits) do by[h.strand][#by[h.strand] + 1] = h end
  H.ok(#by[1] + #by[-1] > 0, "no BsaI site")
  local cells = painted(st)
  local where = {}
  for _, c in ipairs(cells) do where[c.bp .. c.kind] = true end
  for _, h in ipairs(st.search.hits) do
    local kind = h.strand == 1 and "seq" or "comp"
    for bp = h.start, h.stop do H.ok(where[bp .. kind], ("bp %d on %s row"):format(bp, kind)) end
  end
  local n = 0
  for _, h in ipairs(st.search.hits) do n = n + (h.stop - h.start + 1) end
  H.eq(n, #cells, "nothing else painted")
end)

case("origin_spanning_hit_is_split_across_rows", function()
  local buf, st, rec = open("ncbi/J01749.gb")
  local L = rec:len()
  local motif = (rec.seq:sub(L - 3, L) .. rec.seq:sub(1, 4)):upper()
  vim.cmd("DnaFind " .. motif)
  local wrap
  for _, h in ipairs(st.search.hits) do if h.stop < h.start then wrap = h end end
  H.ok(wrap, "wrapping hit found")
  H.eq(L - 3, wrap.start); H.eq(4, wrap.stop)
  local seen = {}
  for _, c in ipairs(painted(st)) do seen[c.bp] = true end
  for _, bp in ipairs({ L - 3, L - 2, L - 1, L, 1, 2, 3, 4 }) do H.ok(seen[bp], "bp " .. bp .. " painted") end
end)

case("gn_selects_and_edits_refresh_hits", function()
  local buf, st, rec = open("ncbi/L09137.gb")
  vim.cmd("DnaFind GAATTC")
  keys("gn")
  H.eq(396, st.selection.start); H.eq(401, st.selection.stop)
  -- deleting a base of the site removes the hit
  buffer.goto_bp(0, 398)
  keys("x")
  H.eq(0, #st.search.hits, "hit recomputed after the edit")
  H.ok(buffer.winbar():find("/GAATTC: none", 1, true))
  keys("u")
  H.eq(1, #st.search.hits, "and back after undo")
end)

case("mismatch_and_amino_acid_search", function()
  local buf, st, rec = open("ncbi/L09137.gb")
  vim.cmd("DnaFind GAATTC 1")
  H.ok(#st.search.hits > 1)
  local exact = 0
  for _, h in ipairs(st.search.hits) do if h.mismatches == 0 then exact = exact + 1 end end
  H.eq(1, exact)
  -- the bla signal peptide, translated on the minus strand of pBR322
  buf, st, rec = open("ncbi/J01749.gb")
  vim.cmd("DnaFind aa:MSIQHFRVALIPFFAAFCLPVFA")
  H.eq(1, #st.search.hits)
  local h = st.search.hits[1]
  H.eq(-1, h.strand); H.eq(4153, h.stop, "starts at the bla start codon (complement(3293..4153))")
  H.eq(4153 - 23 * 3 + 1, h.start)
end)

case("circular_view_draws_hit_marks", function()
  local buf, st = open("ncbi/L09137.gb")
  vim.cmd("DnaFind GATC")
  buffer.set_view(buf, "circular")
  local cur, other = 0, 0
  for row0 = 0, #st.layout.lines - 1 do
    for _, h in ipairs(st.layout:row_highlights(row0)) do
      if h[3] == "PlasmidSearchCurrent" then cur = cur + 1 elseif h[3] == "PlasmidSearch" then other = other + 1 end
    end
  end
  H.ok(cur > 0, "current hit marked")
  H.ok(other > 0, "other hits marked")
  buffer.set_view(buf, "linear")
end)

case("slash_prompt_accepts_mismatch_suffix", function()
  local buf, st = open("ncbi/L09137.gb")
  local orig = vim.ui.input
  vim.ui.input = function(opts, cb) H.eq("/", opts.prompt); cb("gaattc 1") end
  keys("/")
  vim.ui.input = orig
  H.eq(1, st.search.k)
  H.eq("GAATTC", st.search.motif)
  vim.ui.input = function(_, cb) cb("BOGUS!") end
  keys("?")
  vim.ui.input = orig
  H.eq(1, st.search.k, "invalid motif keeps the previous search")
end)

return T
