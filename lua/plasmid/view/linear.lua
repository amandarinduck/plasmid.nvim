-- Linear sequence view.
--
-- Layout per block of W bases:
--      1 ATGAGCTTAA GGCATCGATT ...   60     sequence (ruler numbers both ends)
--        TACTCGAATT CCGTAGCTAA ...          complement strand
--        ━━━━━━ lacZα ━━━━━━━━━━━━▶         feature track(s), packed
--          M  S  L  R  H  R ...             translation under codon middles
--
-- Every row records which block it belongs to, so any cursor position maps
-- back to a base pair. Highlights are produced per row on demand and applied
-- by a decoration provider only for rows on screen.
local alphabet = require("plasmid.model.alphabet")
local vf = require("plasmid.view.features")
local hl = require("plasmid.view.highlights")

local M = {}

local BODY, INTRON, FWD, REV, SITE = "━", "╌", "▶", "◀", "┃"

local function utf8_chars(s)
  local out = {}
  for ch in s:gmatch("[%z\1-\127\194-\244][\128-\191]*") do out[#out + 1] = ch end
  return out
end

--- Turn a cell array into text + byte-range highlight runs + feature segments.
local function finalize(prefix_w, cells, hls, fids)
  local parts = { (" "):rep(prefix_w) }
  local runs, segs = {}, {}
  local byte = prefix_w
  local run_hl, run_b0
  local seg_fid, seg_c0
  local n = #cells
  for c = 1, n + 1 do
    local ch = cells[c]
    local h = c <= n and hls[c] or nil
    local fid = c <= n and fids[c] or nil
    if h ~= run_hl then
      if run_hl then runs[#runs + 1] = { run_b0, byte, run_hl } end
      run_hl, run_b0 = h, byte
    end
    if fid ~= seg_fid then
      if seg_fid then segs[#segs + 1] = { seg_c0, prefix_w + c - 2, seg_fid } end
      seg_fid, seg_c0 = fid, prefix_w + c - 1
    end
    if c <= n then
      parts[#parts + 1] = ch
      byte = byte + #ch
    end
  end
  return (table.concat(parts):gsub("%s+$", "")), runs, segs
end

---@param rec table Record
---@param opts table plasmid config options
---@param win_width integer
function M.build(rec, opts, win_width, extra)
  extra = extra or {}
  local o = opts.linear
  local L = rec:len()
  local numw = #tostring(math.max(L, 1))
  local prefix_w = numw + 1
  -- residues per group (a space between groups); 0 = continuous with a ruler
  local G = math.max(0, math.floor(tonumber(o.group) or 10))
  -- line widths stay a multiple of 10 (and of the group size) so the ruler and
  -- the end-of-line numbers stay round
  local step = G > 0 and (G * 10 / M.gcd(G, 10)) or 10   -- lcm(G, 10)
  local W = o.line_width
  if W == "auto" or type(W) ~= "number" then
    local avail = win_width - prefix_w - (numw + 2) - 1
    local per = G > 0 and (G + 1) / G or 1   -- columns per residue incl. gaps
    W = math.max(1, math.floor((avail + 1) / (per * step))) * step
  end
  W = math.max(step, W)
  local content_w = W + (G > 0 and math.ceil(W / G) - 1 or 0)
  local function col0(o0) return o0 + (G > 0 and math.floor(o0 / G) or 0) end

  local case = o.case
  if case == "auto" then case = rec.format == "fasta" and "preserve" or "upper" end
  local seq = rec.seq:tostring()
  if case == "upper" then seq = seq:upper() end

  local nuc = alphabet.is_nucleotide(rec.alphabet)
  -- complement in the record's own letters: GenBank mRNA records are RNA
  -- molecules written with T, so their complement uses A for T as well
  local comp_alpha = rec.alphabet
  if comp_alpha == "rna" and not seq:find("[Uu]") then comp_alpha = "dna" end
  local unit = rec:unit()
  local all_items = vf.visible(rec, opts)
  local items = all_items
  if extra.primers then
    -- primers are drawn as bases in their own lanes
    local drawn = {}
    for _, p in ipairs(extra.primers) do drawn[p.f] = true end
    items = {}
    for _, it in ipairs(all_items) do if not drawn[it.f] then items[#items + 1] = it end end
  end
  local ntracks = vf.pack(items, L, 0, opts)
  -- rank = position in priority order (1 = most important)
  for rank, it in ipairs(vf.by_priority(items, opts)) do it.rank = rank end
  local cap = o.max_tracks == "all" and math.huge or (tonumber(o.max_tracks) or 4)
  local any_overflow = false
  local ext_list = {}
  for _, it in ipairs(items) do
    for _, e in ipairs(it.extents) do ext_list[#ext_list + 1] = { e[1], e[2], it } end
  end
  table.sort(ext_list, function(a, b) return a[1] < b[1] end)

  -- restriction-site entries (sorted by top-strand cut)
  local enz = extra.enzymes and #extra.enzymes > 0 and extra.enzymes or nil
  local enz_lanes = math.max(1, tonumber((opts.enzymes or {}).lanes) or 4)
  local enz_wrapped = {}
  for _, g in ipairs(enz or {}) do if g.stop < g.start then enz_wrapped[#enz_wrapped + 1] = g end end
  local function first_entry(pos)
    local lo, hi = 1, enz and #enz + 1 or 1
    while lo < hi do
      local mid = math.floor((lo + hi) / 2)
      if enz[mid].cut < pos then lo = mid + 1 else hi = mid end
    end
    return lo
  end

  -- primers drawn as bases (forward above the sequence, reverse below)
  local primer_cells = {}
  for _, p in ipairs(extra.primers or {}) do
    local cells = extra.primer_cells(p)
    table.sort(cells, function(x, y) return x.pos < y.pos end)
    primer_cells[#primer_cells + 1] = { p = p, cells = cells }
  end
  local primer_lanes_max = math.max(1, tonumber((opts.primers or {}).lanes) or 3)


  local lines, rows = {}, {}
  local function push(text, info)
    lines[#lines + 1] = text
    rows[#lines - 1] = info
  end

  -- header
  local kind = nuc and ((rec.molecule_type or "DNA") .. "  ·  " .. rec.topology) or "protein"
  push((" %s  ·  %d %s  ·  %s  ·  %s"):format(
    rec.name, L, unit, kind, vf.count_text(all_items)
  ), { kind = "title" })
  local desc = rec.description ~= "" and rec.description or (rec.id or "")
  push(" " .. desc, { kind = "header" })
  push("", { kind = "blank" })
  local header_rows = #lines

  local block_row = {}
  local nblocks = math.max(1, math.ceil(L / W))
  local show_comp, show_aa = o.show_complement and nuc, o.show_translation and nuc
  local numfmt = "%" .. numw .. "d "
  local on_codon_map

  ---------------------------------------------------------------------------
  -- Feature index: the items overlapping a block, in time proportional to the
  -- answer rather than to the feature count. Blocks address a flat array
  -- (idx_at[idx_start[b] .. idx_start[b+1] - 1]); an extent covering more than
  -- SPAN blocks goes in `long` instead, so one chromosome-length feature
  -- cannot multiply the index by the number of blocks.
  ---------------------------------------------------------------------------
  local SPAN = 16
  local levels, huge = {}, {}
  do
    --- One level: extents covering at most SPAN of its blocks are bucketed
    --- there; the rest come back for the next, coarser level.
    local function build(mult, src)
      local size = W * mult
      local nb = math.max(1, math.ceil(nblocks / mult))
      local counts, rest = {}, {}
      local function span_of(x)
        return math.floor((math.max(x[1], 1) - 1) / size), math.floor((math.min(x[2], L) - 1) / size)
      end
      for _, x in ipairs(src) do
        local b0, b1 = span_of(x)
        if b1 - b0 > SPAN then
          rest[#rest + 1] = x
        else
          for b = b0, b1 do counts[b + 1] = (counts[b + 1] or 0) + 1 end
        end
      end
      local start, n = {}, 1
      for b = 0, nb - 1 do
        start[b + 1] = n
        n = n + (counts[b + 1] or 0)
      end
      start[nb + 1] = n
      local at_s, at_e, at_i, cur = {}, {}, {}, {}
      for _, x in ipairs(src) do
        local b0, b1 = span_of(x)
        if b1 - b0 <= SPAN then
          for b = b0, b1 do
            local k = cur[b + 1] or start[b + 1]
            at_s[k], at_e[k], at_i[k] = x[1], x[2], x[3]
            cur[b + 1] = k + 1
          end
        end
      end
      return { mult = mult, nb = nb, start = start, s = at_s, e = at_e, i = at_i }, rest
    end
    local src = {}
    if L > 0 then
      for i = 1, #items do
        for _, e in ipairs(items[i].extents) do
          if e[2] >= e[1] then src[#src + 1] = { e[1], e[2], i } end
        end
      end
    end
    local rest
    levels[1], rest = build(1, src)
    levels[2], huge = build(SPAN, rest)
  end

  --- Items overlapping block `b`, most important first. The result array and
  --- the seen-set are reused between calls: a record with tens of thousands of
  --- blocks would otherwise spend most of its layout time in the collector.
  local by_rank = function(x, y) return x.rank < y.rank end
  local scratch, stamp, token = {}, {}, 0
  local span_a, span_b, EMPTY = {}, {}, {}
  local function items_in(b, bs, be)
    local n = 0
    token = token + 1
    for li = 1, #levels do
      local lv = levels[li]
      local cb = math.floor(b / lv.mult)
      local at_s, at_e, at_i = lv.s, lv.e, lv.i
      for k = lv.start[cb + 1], lv.start[cb + 2] - 1 do
        if at_s[k] <= be and at_e[k] >= bs then
          local it = items[at_i[k]]
          if stamp[it] ~= token then stamp[it] = token; n = n + 1; scratch[n] = it end
        end
      end
    end
    for _, e in ipairs(huge) do
      if e[1] <= be and e[2] >= bs then
        local it = items[e[3]]
        if stamp[it] ~= token then stamp[it] = token; n = n + 1; scratch[n] = it end
      end
    end
    for i = #scratch, n + 1, -1 do scratch[i] = nil end
    table.sort(scratch, by_rank)
    return scratch
  end

  ---------------------------------------------------------------------------
  -- Planning a block: which rows it has, without building any text. This is
  -- what lets the view know its own length, and any row's base pair, without
  -- laying out the whole record.
  ---------------------------------------------------------------------------

  --- Lane assignment for the primers on one strand in a block, or nil.
  local function primer_lanes(strand, bs, be)
    if #primer_cells == 0 or not nuc then return nil end
    local limit = content_w + numw + 2
    local pitems = {}
    for _, pc in ipairs(primer_cells) do
      local p = pc.p
      if p.strand == strand then
        local here = {}
        for _, c in ipairs(pc.cells) do
          if c.pos >= bs and c.pos <= be then here[#here + 1] = c end
        end
        if #here > 0 then
          local c0, c1 = col0(here[1].pos - bs), col0(here[#here].pos - bs)
          -- label at the 3' end when it is on this line
          local three = strand == 1 and p.stop or p.start
          local label
          if three >= bs and three <= be then
            label = strand == 1 and ("›" .. p.name) or (p.name .. "‹")
          end
          local lw = label and vim.fn.strdisplaywidth(label) or 0
          local lo = c0 - ((strand == -1 and label) and lw or 0)
          local hi = c1 + ((strand == 1 and label) and lw or 0)
          if lo < 0 then lo, label = c0, nil; hi = c1 end
          if hi >= limit then hi, label = c1, (strand == -1) and label or nil end
          pitems[#pitems + 1] = { p = p, cells = here, lo = lo, hi = hi, label = label }
        end
      end
    end
    if #pitems == 0 then return nil end
    table.sort(pitems, function(x, y) return x.lo < y.lo end)
    local lanes, hidden = {}, 0
    for _, it in ipairs(pitems) do
      local placed = false
      for li = 1, primer_lanes_max do
        local lane = lanes[li]
        if not lane then lane = { last = -2, items = {} }; lanes[li] = lane end
        if it.lo > lane.last + 1 then
          lane.items[#lane.items + 1] = it
          lane.last = it.hi
          placed = true
          break
        end
      end
      if not placed then hidden = hidden + 1 end
    end
    return { lanes = lanes, hidden = hidden }
  end

  --- Enzyme label lanes for a block, or nil.
  local function enzyme_lanes(bs, be)
    local lanes, hidden = {}, 0
    local limit = content_w + numw + 2   -- columns after the prefix
    local function free(lane, c0, c1)
      for _, iv in ipairs(lane) do
        if c0 <= iv[2] + 1 and c1 >= iv[1] - 1 then return false end
      end
      return true
    end
    local i = first_entry(bs)
    while i <= #enz and enz[i].cut <= be do
      local g = enz[i]
      local c = col0(g.cut - bs)
      local text = g.label
      -- right of the cut: "▏Name" (bar on the left edge of the base after
      -- the cut); else left of it: "Name▕" (bar on the right edge of the
      -- base before the cut, on this line)
      local opts_ = { { c, c + #text, "right" } }
      if g.cut - 1 >= bs then
        local cl = col0(g.cut - 1 - bs)
        opts_[2] = { cl - #text, cl, "left" }
      end
      local placed = false
      for li = 1, enz_lanes do
        local lane = lanes[li] or {}
        for _, op in ipairs(opts_) do
          if op[1] >= 0 and op[2] < limit and free(lane, op[1], op[2]) then
            lane[#lane + 1] = { op[1], op[2], g = g, side = op[3], text = text }
            lanes[li] = lane
            placed = true
            break
          end
        end
        if placed then break end
        if not lanes[li] then break end
      end
      if not placed then hidden = hidden + 1 end
      i = i + 1
    end
    if hidden > 0 and #lanes == 0 then lanes[1] = {} end
    if #lanes == 0 then return nil end
    return { lanes = lanes, hidden = hidden }
  end

  --- Is there a bottom-strand cut mark under this block's complement row?
  local function has_bottom_cut(bs, be)
    for k = first_entry(bs - 64), #enz do
      local g = enz[k]
      if g.cut > be + 64 then break end
      local p = g.cut_bottom
      if p >= bs and p <= be then return true end
    end
    return false
  end

  --- The row plan for a block: kinds, counts, and the assignments the text is
  --- drawn from. Costs time proportional to what the block contains.
  local function plan(b)
    local bs = b * W + 1
    local be = math.min(bs + W - 1, L)
    local p = { b = b, bs = bs, be = be }
    local n = 0
    if G == 0 then p.ruler = true; n = n + 1 end
    if enz and nuc then
      p.enz = enzyme_lanes(bs, be)
      if p.enz then n = n + #p.enz.lanes end
    end
    p.pf = primer_lanes(1, bs, be)
    if p.pf then n = n + #p.pf.lanes end
    p.seq_i = n
    n = n + 1
    if show_comp then
      n = n + 1
      if enz and nuc and has_bottom_cut(bs, be) then p.enzcut = true; n = n + 1 end
    end
    p.pr = primer_lanes(-1, bs, be)
    if p.pr then n = n + #p.pr.lanes end

    -- feature rows for this block, filled in priority order. A feature keeps
    -- its global track when that row is free here (visual continuity across
    -- lines); otherwise it takes the first free row; if all `cap` rows are
    -- taken by more important overlapping features it is summarised below.
    local block_items = items_in(b, bs, be)
    local nitems = #block_items
    local by_track, tlist, row_sets, overflow, has_cds
    for i = 1, nitems do
      local it = block_items[i]
      local ns = 0
      for _, e in ipairs(it.extents) do
        local s0, e0 = math.max(e[1], bs), math.min(e[2], be)
        if s0 <= e0 then ns = ns + 1; span_a[ns], span_b[ns] = s0, e0 end
      end
      local row
      for attempt = 0, math.min(cap, nitems) do
        local r = attempt == 0 and it.track or attempt
        if r <= cap and (attempt > 0 or it.track <= cap) then
          local set = row_sets and row_sets[r]
          local ok = true
          if set then
            for k = 1, ns do
              -- 1 bp of air between neighbours so their arrowheads stay distinct
              if vf._overlaps(set, span_a[k] - 1, span_b[k] + 1) then ok = false break end
            end
          end
          if ok then row = r break end
        end
      end
      if row then
        row_sets = row_sets or {}
        row_sets[row] = row_sets[row] or {}
        for k = 1, ns do vf._insert(row_sets[row], span_a[k], span_b[k]) end
        by_track = by_track or {}
        if not by_track[row] then
          by_track[row] = {}
          tlist = tlist or {}
          tlist[#tlist + 1] = row
        end
        table.insert(by_track[row], it)
        if it.f.type == "CDS" or it.f.type == "ORF" then
          has_cds = has_cds or {}
          has_cds[row] = true
        end
      else
        overflow = overflow or {}
        overflow[#overflow + 1] = it
      end
    end
    if tlist then
      table.sort(tlist)
      for _, t in ipairs(tlist) do
        n = n + 1
        if show_aa and has_cds and has_cds[t] then n = n + 1 end
      end
    end
    p.tlist, p.by_track, p.has_cds = tlist or EMPTY, by_track or EMPTY, has_cds or EMPTY
    if overflow then p.overflow = overflow; n = n + 1 end
    n = n + 1   -- trailing blank
    p.n = n
    return p
  end

  ---------------------------------------------------------------------------
  -- Rendering a planned block into text and row metadata.
  ---------------------------------------------------------------------------

  local function primer_render(strand, bs, be, pl)
    local lanes, hidden = pl.lanes, pl.hidden
    local out = {}
    for li = 1, #lanes do
      local cells, hls, ids = {}, {}, {}
      local n = 0
      local function put(c, ch, h, id)
        for q = n + 1, c - 1 do cells[q] = cells[q] or " " end
        cells[c], hls[c], ids[c] = ch, h, id
        n = math.max(n, c)
      end
      for _, it in ipairs(lanes[li].items) do
        local id = it.p.f.id
        if it.label and strand == -1 then
          local k = 0
          for ch in it.label:gmatch("[%z\1-\127\194-\244][\128-\191]*") do
            put(it.lo + 1 + k, ch, ch == "‹" and "PlasmidPrimerArrow" or "PlasmidPrimerLabel", id)
            k = k + 1
          end
        end
        for _, c in ipairs(it.cells) do
          local g = c.kind == "tail" and "PlasmidPrimerTail" or (c.kind == "mismatch" and "PlasmidPrimerMismatch" or "PlasmidPrimer")
          put(col0(c.pos - bs) + 1, c.ch, g, id)
        end
        -- fill group gaps inside the primer so it reads as one oligo
        if G > 0 then
          local ca, cb = col0(it.cells[1].pos - bs) + 1, col0(it.cells[#it.cells].pos - bs) + 1
          for q = ca, cb do if not cells[q] or cells[q] == " " then cells[q], hls[q], ids[q] = " ", nil, id end end
        end
        if it.label and strand == 1 then
          local k = 0
          local c1 = col0(it.cells[#it.cells].pos - bs) + 1
          for ch in it.label:gmatch("[%z\1-\127\194-\244][\128-\191]*") do
            k = k + 1
            put(c1 + k, ch, ch == "›" and "PlasmidPrimerArrow" or "PlasmidPrimerLabel", id)
          end
        end
      end
      if li == #lanes and hidden > 0 then
        local note = (" (+%d primers)"):format(hidden)
        for k = 1, #note do put(n + 1, note:sub(k, k), "PlasmidOverflow") end
      end
      for q = 1, n do cells[q] = cells[q] or " " end
      local text, runs, segs = finalize(prefix_w, cells, hls, ids)
      out[#out + 1] = { text = text, info = { kind = "primer", bs = bs, be = be, runs = runs, segs = segs, strand = strand } }
    end
    return out
  end

  --- Emit a planned block's rows through `sink(text, info)`.
  local function render(p, sink)
    local b, bs, be = p.b, p.bs, p.be
    -- ruler (continuous mode): numbers ending over every 10th residue
    if p.ruler then
      local cells = {}
      for c = 1, content_w do cells[c] = " " end
      local first = math.ceil(bs / 10) * 10
      for pos = first, be, 10 do
        local txt = tostring(pos)
        local last = col0(pos - bs) + 1
        for k = 1, #txt do
          local c = last - #txt + k
          if c >= 1 then cells[c] = txt:sub(k, k) end
        end
      end
      sink(((" "):rep(prefix_w) .. table.concat(cells)):gsub("%s+$", ""), { kind = "ruler", b = b, bs = bs, be = be })
    end
    -- restriction-site label lanes: "▏Name" with the bar over the first base
    -- after the top-strand cut (the bar sits on the cell's left edge, i.e. on
    -- the cut itself)
    if p.enz then
      local lanes, hidden = p.enz.lanes, p.enz.hidden
      for li = #lanes, 1, -1 do
        local lane = lanes[li]
        table.sort(lane, function(x, y) return x[1] < y[1] end)
        local cells, hls, ids = {}, {}, {}
        local n = 0
        local function put(c, ch, h, id)
          for q = n + 1, c - 1 do cells[q] = " " end
          cells[c], hls[c], ids[c] = ch, h, id
          n = math.max(n, c)
        end
        for _, iv in ipairs(lane) do
          local grp = iv.g.count == 1 and "PlasmidEnzymeUnique" or "PlasmidEnzyme"
          if iv.side == "right" then
            put(iv[1] + 1, "▏", "PlasmidEnzymeCut", iv.g)
            for k = 1, #iv.text do put(iv[1] + 1 + k, iv.text:sub(k, k), grp, iv.g) end
          else
            for k = 1, #iv.text do put(iv[1] + k, iv.text:sub(k, k), grp, iv.g) end
            put(iv[2] + 1, "▕", "PlasmidEnzymeCut", iv.g)
          end
        end
        if li == #lanes and hidden > 0 then
          local note = (" (+%d)"):format(hidden)
          for k = 1, #note do put(n + 1, note:sub(k, k), "PlasmidOverflow") end
        end
        local text, runs, esegs = finalize(prefix_w, cells, hls, ids)
        sink(text, { kind = "enz", b = b, bs = bs, be = be, runs = runs, esegs = esegs })
      end
    end
    if p.pf then
      local fw = primer_render(1, bs, be, p.pf)
      for i = #fw, 1, -1 do fw[i].info.b = b; sink(fw[i].text, fw[i].info) end
    end
    -- sequence row
    local chunk = seq:sub(bs, be)
    local grouped = G > 0 and chunk:gsub("(" .. ("."):rep(G) .. ")", "%1 "):gsub(" $", "") or chunk
    local line = numfmt:format(bs) .. grouped
    line = line .. (" "):rep(prefix_w + content_w - #line) .. "  " .. be
    local info = { kind = "seq", b = b, bs = bs, be = be }
    if enz and nuc then
      -- underline recognition sites on the sequence row
      local runs = {}
      local function add(s, e)
        s, e = math.max(s, bs), math.min(e, be)
        if s <= e then runs[#runs + 1] = { prefix_w + col0(s - bs), prefix_w + col0(e - bs) + 1, "PlasmidEnzymeSite", 120 } end
      end
      for k = first_entry(bs - 64), #enz do
        local g = enz[k]
        if g.cut > be + 64 then break end
        if g.stop >= g.start then add(g.start, g.stop) end
      end
      for _, g in ipairs(enz_wrapped) do add(g.start, L); add(1, g.stop) end
      info.site_runs = runs
    end
    sink(line, info)
    if show_comp then
      sink((" "):rep(prefix_w) .. alphabet.complement_for(grouped, comp_alpha), { kind = "comp", b = b, bs = bs, be = be })
      -- bottom-strand cut marks under the complement
      if p.enzcut then
        local cells, hls = {}, {}
        local n = 0
        for k = first_entry(bs - 64), #enz do
          local g = enz[k]
          if g.cut > be + 64 then break end
          local pos = g.cut_bottom
          if pos >= bs and pos <= be then
            local c = col0(pos - bs)
            for q = n + 1, c do cells[q] = cells[q] or " " end
            cells[c + 1], hls[c + 1] = "▏", "PlasmidEnzymeCut"
            n = math.max(n, c + 1)
          end
        end
        for q = 1, n do cells[q] = cells[q] or " " end
        local text, runs = finalize(prefix_w, cells, hls, {})
        sink(text, { kind = "enzcut", b = b, bs = bs, be = be, runs = runs })
      end
    end
    if p.pr then
      for _, r in ipairs(primer_render(-1, bs, be, p.pr)) do r.info.b = b; sink(r.text, r.info) end
    end

    local function cc(pos) return col0(pos - bs) + 1 end

    for _, t in ipairs(p.tlist) do
      local cells, hls, fids = {}, {}, {}
      for c = 1, content_w do cells[c] = " " end
      for _, it in ipairs(p.by_track[t]) do
        local f, h = it.f, it.hl
        local function set(c, ch)
          cells[c], hls[c], fids[c] = ch, h, f.id
        end
        for _, g in ipairs(vf.parts_by_extent(it)) do
          local prev
          for _, part in ipairs(g) do
            if part.stop < part.start then
              if part.stop >= bs and part.stop <= be then set(cc(part.stop), SITE) end
            else
              local s0, e0 = math.max(part.start, bs), math.min(part.stop, be)
              if s0 <= e0 then
                for c = cc(s0), cc(e0) do set(c, BODY) end
              end
            end
            if prev then
              local s0, e0 = math.max(prev.stop + 1, bs), math.min(part.start - 1, be)
              if s0 <= e0 then
                for c = cc(s0), cc(e0) do
                  if fids[c] == nil then set(c, INTRON) end
                end
              end
            end
            if not prev or part.stop > prev.stop then prev = part end
          end
        end
        local apos, strand = vf.arrow(f)
        if apos and apos >= bs and apos <= be then set(cc(apos), strand == 1 and FWD or REV) end

        -- label in the longest bar run of this feature on this row
        local best_a, best_len, run_a = nil, 0, nil
        for c = 1, content_w + 1 do
          local inbar = c <= content_w and fids[c] == f.id and cells[c] == BODY
          if inbar and not run_a then run_a = c end
          if not inbar and run_a then
            if c - run_a > best_len then best_a, best_len = run_a, c - run_a end
            run_a = nil
          end
        end
        local lc = best_a and utf8_chars(it.label) or {}
        -- a label is shown whole, or truncated to >= 2 chars + "…", or not at all
        if best_a and (#lc + 2 <= best_len or best_len >= 5) then
          if #lc + 2 > best_len then
            local keep = best_len - 3
            local t2 = {}
            for k = 1, keep do t2[k] = lc[k] end
            t2[#t2 + 1] = "…"
            lc = t2
          end
          local start = best_a + math.floor((best_len - #lc) / 2)
          set(start - 1, " ")
          for k, ch in ipairs(lc) do set(start + k - 1, ch) end
          if start + #lc <= best_a + best_len - 1 then set(start + #lc, " ") end
        end
      end
      local text, runs, segs = finalize(prefix_w, cells, hls, fids)
      sink(text, { kind = "track", b = b, bs = bs, be = be, runs = runs, segs = segs, items = p.by_track[t] })

      if show_aa and p.has_cds[t] then
        local acells, ahls, afids = {}, {}, {}
        for c = 1, content_w do acells[c] = " " end
        for _, it in ipairs(p.by_track[t]) do
          if it.f.type == "CDS" or it.f.type == "ORF" then
            local map = vf.codon_map(rec, it.f)
            if on_codon_map then on_codon_map(it.f) end
            for pos = bs, be do
              local aa = map[pos]
              if aa then
                local c = cc(pos)
                acells[c], ahls[c], afids[c] = aa, it.hl, it.f.id
              end
            end
          end
        end
        local atext, aruns, asegs = finalize(prefix_w, acells, ahls, afids)
        sink(atext, { kind = "aa", b = b, bs = bs, be = be, runs = aruns, segs = asegs })
      end
    end
    if p.overflow then
      local counts, order = {}, {}
      for _, it in ipairs(p.overflow) do
        if not counts[it.f.type] then order[#order + 1] = it.f.type end
        counts[it.f.type] = (counts[it.f.type] or 0) + 1
      end
      local parts = {}
      for _, t in ipairs(order) do parts[#parts + 1] = counts[t] > 1 and (t .. " ×" .. counts[t]) or t end
      local hidden = {}
      for _, it in ipairs(p.overflow) do hidden[#hidden + 1] = it.f.id end
      sink((" "):rep(prefix_w) .. ("+%d more: %s"):format(#p.overflow, table.concat(parts, ", ")),
        { kind = "overflow", b = b, bs = bs, be = be, hidden = hidden, hidden_items = p.overflow })
    end
    sink("", { kind = "blank", b = b, bs = bs, be = be })
  end

  ---------------------------------------------------------------------------
  -- The whole record, or a window of it
  ---------------------------------------------------------------------------
  local wcfg = o.window or {}
  local threshold = tonumber(wcfg.threshold) or 300000
  local windowed = wcfg.enabled == true or (wcfg.enabled ~= false and L >= threshold)
  local any_overflow = false
  local offsets, total

  if windowed then
    offsets = {}
    total = header_rows
    for b = 0, nblocks - 1 do
      local p = plan(b)
      offsets[b + 1] = total
      total = total + p.n
      if p.overflow then any_overflow = true end
    end
    for i = header_rows + 1, total do lines[i] = "" end
  else
    for b = 0, nblocks - 1 do
      local p = plan(b)
      block_row[b] = #lines + p.seq_i
      render(p, push)
      if p.overflow then any_overflow = true end
    end
  end

  if L == 0 then
    if rec.no_sequence then
      push(("  (no sequence in this file: %s)"):format(rec.contig and ("CON record, CONTIG " .. rec.contig:sub(1, 60))
        or ("WGS master, " .. (rec.wgs or "contig records listed elsewhere"))), { kind = "header" })
      push("  :DnaFetch " .. (rec.id or rec.name) .. " opens it with its residues", { kind = "header" })
    else
      push("  (empty sequence)", { kind = "header" })
    end
  end

  if any_overflow then
    lines[1] = lines[1] .. ("   (≤%d rows per line, %s shows all)"):format(cap,
      require("plasmid.config").options.keymaps.expand_tracks or ":DnaTracks all")
  end

  local layout = {
    view = "linear", lines = lines, rows = rows, header_rows = header_rows,
    W = W, G = G, prefix_w = prefix_w, numw = numw, content_w = content_w, has_comp = show_comp,
    block_row = block_row, items = items, ntracks = ntracks, L = L,
    color_bases = o.color_bases, width = win_width,
    residue_groups = nuc and hl.base_groups or hl.residue_groups,
  }
  if windowed then
    layout.windowed = true
    layout.nblocks = nblocks
    layout.offsets = offsets
    layout.margin = math.max(0, tonumber(wcfg.margin) or 2)
    layout.cache_blocks = math.max(64, tonumber(wcfg.cache_blocks) or 0)
    layout._plan, layout._render = plan, render
    layout._built, layout._order, layout._plans, layout._maps = {}, {}, {}, {}
    on_codon_map = function(f) layout._maps[#layout._maps + 1] = f end
  end
  return setmetatable(layout, { __index = M })
end


---------------------------------------------------------------------------
-- Layout queries (methods on the returned layout)
---------------------------------------------------------------------------

---------------------------------------------------------------------------
-- Windowed layouts build their rows on demand. These answer the same
-- questions for a row whose text has not been built.
---------------------------------------------------------------------------

--- Block containing a buffer row, or nil for the header and past the end.
function M.block_at_row(self, row0)
  if not self.windowed then
    local info = self.rows[row0]
    return info and info.b
  end
  if row0 < self.header_rows or row0 >= #self.lines then return nil end
  local off = self.offsets
  local lo, hi = 1, self.nblocks
  while lo < hi do
    local mid = math.floor((lo + hi + 1) / 2)
    if off[mid] <= row0 then lo = mid else hi = mid - 1 end
  end
  return lo - 1
end

--- The row plan of a block, from the cache or freshly computed.
function M.plan_of(self, b)
  local p = self._plans[b]
  if not p then
    p = self._plan(b)
    self._plans[b] = p
  end
  return p
end

--- Build a block's text and row metadata (not the buffer). Returns its first
--- row and the lines, or nil when it was already built.
function M.build_block(self, b)
  if self._built[b] then return nil end
  local p = self:plan_of(b)
  local base = self.offsets[b + 1]
  local out, n = {}, 0
  self._render(p, function(text, info)
    n = n + 1
    out[n] = text
    self.lines[base + n] = text
    self.rows[base + n - 1] = info
  end)
  self._built[b] = true
  self._order[#self._order + 1] = b
  return base, out
end

--- Forget a block: its rows go back to blank, so a long scroll through a
--- genome does not accumulate the whole record's layout.
function M.drop_block(self, b)
  if not self._built[b] then return nil end
  local p = self._plans[b]
  local base = self.offsets[b + 1]
  local blanks = {}
  for i = 1, p.n do
    blanks[i] = ""
    self.lines[base + i] = ""
    self.rows[base + i - 1] = nil
  end
  self._built[b], self._plans[b] = nil, nil
  return base, blanks
end

--- Build every block covering rows [row0, row1] plus a margin of screens.
--- Returns the edits the buffer needs: { { row, lines }, ... }.
function M.ensure(self, row0, row1)
  if not self.windowed then return {} end
  local span = math.max(1, row1 - row0 + 1)
  local lo = math.max(0, self:block_at_row(math.max(row0 - span * self.margin, self.header_rows)) or 0)
  local hi = self:block_at_row(math.min(row1 + span * self.margin, #self.lines - 1)) or (self.nblocks - 1)
  -- a caller asking for an implausible span (a jump across the record, a
  -- window the size of a genome) gets the middle of it, not all of it
  local MAX = 4096
  if hi - lo > MAX then
    local mid = math.floor((lo + hi) / 2)
    lo, hi = math.max(0, mid - MAX / 2), math.min(self.nblocks - 1, mid + MAX / 2)
  end
  local edits = {}
  for b = lo, hi do
    local base, built = self:build_block(b)
    if base then edits[#edits + 1] = { base, built } end
  end
  -- evict the blocks built longest ago, keeping the live window
  local keep = math.max(self.cache_blocks, (hi - lo + 1) * 2)
  if #self._order > keep then
    local drop, kept = {}, {}
    for _, b in ipairs(self._order) do
      if #drop < #self._order - keep and (b < lo or b > hi) then drop[#drop + 1] = b
      else kept[#kept + 1] = b end
    end
    for _, b in ipairs(drop) do
      local base, blanks = self:drop_block(b)
      if base then edits[#edits + 1] = { base, blanks } end
    end
    self._order = kept
    -- codon maps built for rows that are gone cost more memory than they save
    if #self._maps > 512 then
      for i = 1, #self._maps - 256 do
        local f = self._maps[i]
        f._codon_map, f._codon_rev, f._codon_idx, f._codon_pos = nil, nil, nil, nil
      end
      local rest = {}
      for i = #self._maps - 255, #self._maps do rest[#rest + 1] = self._maps[i] end
      self._maps = rest
    end
  end
  return edits
end

--- Row metadata, building the row's block first if it has not been built.
function M.info(self, row0)
  local info = self.rows[row0]
  if info or not self.windowed then return info end
  local b = self:block_at_row(row0)
  if not b then return nil end
  self:build_block(b)
  return self.rows[row0]
end

--- Base pair at a display column (0-based) on a row (0-based).
function M.bp_at(self, row0, vcol0)
  local info = self.rows[row0]
  if not info and self.windowed then
    -- every row of a block covers the same bases, so an unbuilt row still
    -- knows which ones they are
    local b = self:block_at_row(row0)
    if b then info = { bs = b * self.W + 1, be = math.min((b + 1) * self.W, self.L) } end
  end
  if not info or not info.bs then return nil end
  local c = math.max(0, vcol0 - self.prefix_w)
  local G = self.G
  if G == 0 then return math.max(info.bs, math.min(info.bs + c, info.be)) end
  local g, w = math.floor(c / (G + 1)), c % (G + 1)
  if w == G then w = G - 1 end   -- the gap belongs to the residue before it
  return math.max(info.bs, math.min(info.bs + g * G + w, info.be))
end

--- Row (0-based) and display column (0-based) of a base pair on its sequence row.
function M.pos_of(self, bp)
  bp = math.max(1, math.min(bp, math.max(self.L, 1)))
  local b = math.floor((bp - 1) / self.W)
  local o0 = (bp - 1) % self.W
  local row = self.block_row[b]
  if row == nil and self.windowed then row = self.offsets[b + 1] + self:plan_of(b).seq_i end
  return row, self.prefix_w + M.col0(self, o0)
end

--- Display column (0-based, after the prefix) of residue offset o0 in a line.
function M.col0(self, o0)
  return o0 + (self.G > 0 and math.floor(o0 / self.G) or 0)
end

function M.gcd(a, b)
  while b ~= 0 do a, b = b, a % b end
  return a
end

--- Feature id drawn at a display column on a row, if any.
function M.feature_at(self, row0, vcol0)
  local info = M.info(self, row0)
  if not info or not info.segs then return nil end
  for _, s in ipairs(info.segs) do
    if vcol0 >= s[1] and vcol0 <= s[2] then return s[3] end
  end
end

-- Selection runs on a sequence/complement row (ASCII rows: byte = column).
local function selection_runs(self, info)
  local sel = self.selection
  if not sel then return nil end
  local ranges = sel.start > sel.stop and { { sel.start, self.L }, { 1, sel.stop } } or { { sel.start, sel.stop } }
  local out = {}
  for _, r in ipairs(ranges) do
    local s, e = math.max(r[1], info.bs), math.min(r[2], info.be)
    if s <= e then
      local o0, o1 = s - info.bs, e - info.bs
      out[#out + 1] = { self.prefix_w + M.col0(self, o0), self.prefix_w + M.col0(self, o1) + 1,
        "PlasmidSelection", 150 }
    end
  end
  return out
end

-- Search-hit runs on a sequence/complement row. Hits on the bottom strand are
-- shown on the complement row when there is one.
local function search_runs(self, info, kind)
  local S = self.search
  if not S or #S.hits == 0 then return nil end
  local on_comp = kind == "comp"
  local out = {}
  local function add(s, e, i)
    s, e = math.max(s, info.bs), math.min(e, info.be)
    if s > e then return end
    local o0, o1 = s - info.bs, e - info.bs
    out[#out + 1] = { self.prefix_w + M.col0(self, o0), self.prefix_w + M.col0(self, o1) + 1,
      i == S.idx and "PlasmidSearchCurrent" or "PlasmidSearch", 160 }
  end
  local function wanted(h)
    local bottom = h.strand == -1 and self.has_comp
    return bottom == on_comp
  end
  -- first hit that could reach this row (hits sorted by start; maxlen bounds the look-back)
  local lo, hi = 1, #S.hits
  local target = info.bs - S.maxlen
  while lo < hi do
    local mid = math.floor((lo + hi) / 2)
    if S.hits[mid].start < target then lo = mid + 1 else hi = mid end
  end
  for i = lo, #S.hits do
    local h = S.hits[i]
    if h.start > info.be then break end
    if h.stop >= h.start and wanted(h) then add(h.start, h.stop, i) end
  end
  for _, i in ipairs(S.wrapped) do
    local h = S.hits[i]
    if wanted(h) then add(h.start, self.L, i); add(1, h.stop, i) end
  end
  return out
end

--- Byte-range highlights for a row: { {b0, b1, group, priority?}, ... }.
function M.row_highlights(self, row0)
  local info = self.rows[row0]
  if not info then return {} end
  local k = info.kind
  if (k == "seq" or k == "comp") and (self.selection or self.search) then
    local base = M._row_highlights(self, row0, info)
    for _, r in ipairs(search_runs(self, info, k) or {}) do base[#base + 1] = r end
    for _, r in ipairs(selection_runs(self, info) or {}) do base[#base + 1] = r end
    return base
  end
  return M._row_highlights(self, row0, info)
end

function M._row_highlights(self, row0, info)
  local k = info.kind
  if k == "title" then return { { 0, -1, "PlasmidTitle" } } end
  if k == "header" then return { { 0, -1, "PlasmidHeader" } } end
  if k == "track" or k == "aa" or k == "enz" or k == "enzcut" or k == "primer" then return info.runs end
  if k == "comp" then return { { 0, -1, "PlasmidComplement" } } end
  if k == "ruler" then return { { 0, -1, "PlasmidTick" } } end
  if k == "overflow" then return { { 0, -1, "PlasmidOverflow" } } end
  if k == "seq" then
    local out = { { 0, self.numw, "PlasmidRuler" } }
    local tail = self.prefix_w + self.content_w
    if self.color_bases then
      local line = self.lines[row0 + 1]
      for c = self.prefix_w + 1, math.min(#line, tail) do
        local g = self.residue_groups[line:sub(c, c)]
        if g then out[#out + 1] = { c - 1, c, g } end
      end
    end
    out[#out + 1] = { tail, -1, "PlasmidRuler" }
    for _, r in ipairs(info.site_runs or {}) do out[#out + 1] = r end
    return out
  end
  return {}
end

return M
