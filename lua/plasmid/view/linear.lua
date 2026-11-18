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
  local ptr, active = 1, {}
  local show_comp, show_aa = o.show_complement and nuc, o.show_translation and nuc
  local numfmt = "%" .. numw .. "d "

  for b = 0, nblocks - 1 do
    local bs = b * W + 1
    local be = math.min(bs + W - 1, L)
    while ptr <= #ext_list and ext_list[ptr][1] <= be do
      active[#active + 1] = ext_list[ptr]
      ptr = ptr + 1
    end
    local j = 0
    for i = 1, #active do
      if active[i][2] >= bs then
        j = j + 1
        active[j] = active[i]
      end
    end
    for i = #active, j + 1, -1 do active[i] = nil end

    -- ruler (continuous mode): numbers ending over every 10th residue
    if G == 0 then
      local cells = {}
      for c = 1, content_w do cells[c] = " " end
      local first = math.ceil(bs / 10) * 10
      for p = first, be, 10 do
        local txt = tostring(p)
        local last = col0(p - bs) + 1
        for k = 1, #txt do
          local c = last - #txt + k
          if c >= 1 then cells[c] = txt:sub(k, k) end
        end
      end
      push(((" "):rep(prefix_w) .. table.concat(cells)):gsub("%s+$", ""), { kind = "ruler", b = b, bs = bs, be = be })
    end
    -- restriction-site label lanes: "▏Name" with the bar over the first base
    -- after the top-strand cut (the bar sits on the cell's left edge, i.e. on
    -- the cut itself)
    if enz and nuc then
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
        push(text, { kind = "enz", b = b, bs = bs, be = be, runs = runs, esegs = esegs })
      end
    end
    -- primer lanes for this line: returns rows to push (outermost first)
    local function primer_rows(strand)
      if #primer_cells == 0 or not nuc then return {} end
      local limit = content_w + numw + 2
      local items = {}
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
            items[#items + 1] = { p = p, cells = here, lo = lo, hi = hi, label = label }
          end
        end
      end
      if #items == 0 then return {} end
      table.sort(items, function(x, y) return x.lo < y.lo end)
      local lanes, hidden = {}, 0
      for _, it in ipairs(items) do
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
            local a, b = col0(it.cells[1].pos - bs) + 1, col0(it.cells[#it.cells].pos - bs) + 1
            for q = a, b do if not cells[q] or cells[q] == " " then cells[q], hls[q], ids[q] = " ", nil, id end end
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
        out[#out + 1] = { text = text, info = { kind = "primer", b = b, bs = bs, be = be, runs = runs, segs = segs, strand = strand } }
      end
      return out
    end
    do
      local fw = primer_rows(1)
      for i = #fw, 1, -1 do push(fw[i].text, fw[i].info) end
    end
    -- sequence row
    local chunk = seq:sub(bs, be)
    local grouped = G > 0 and chunk:gsub("(" .. ("."):rep(G) .. ")", "%1 "):gsub(" $", "") or chunk
    local line = numfmt:format(bs) .. grouped
    line = line .. (" "):rep(prefix_w + content_w - #line) .. "  " .. be
    block_row[b] = #lines
    local info = { kind = "seq", b = b, bs = bs, be = be }
    push(line, info)
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
    if show_comp then
      push((" "):rep(prefix_w) .. alphabet.complement_for(grouped, comp_alpha), { kind = "comp", b = b, bs = bs, be = be })
      -- bottom-strand cut marks under the complement
      if enz and nuc then
        local cells, hls = {}, {}
        local n = 0
        for k = first_entry(bs - 64), #enz do
          local g = enz[k]
          if g.cut > be + 64 then break end
          local p = g.cut_bottom
          if p >= bs and p <= be then
            local c = col0(p - bs)
            for q = n + 1, c do cells[q] = cells[q] or " " end
            cells[c + 1], hls[c + 1] = "▏", "PlasmidEnzymeCut"
            n = math.max(n, c + 1)
          end
        end
        if n > 0 then
          for q = 1, n do cells[q] = cells[q] or " " end
          local text, runs = finalize(prefix_w, cells, hls, {})
          push(text, { kind = "enzcut", b = b, bs = bs, be = be, runs = runs })
        end
      end
    end
    for _, r in ipairs(primer_rows(-1)) do push(r.text, r.info) end

    -- feature rows for this block, filled in priority order. A feature keeps
    -- its global track when that row is free here (visual continuity across
    -- lines); otherwise it takes the first free row; if all `cap` rows are
    -- taken by more important overlapping features it is summarised below.
    local block_items, seen = {}, {}
    for _, a in ipairs(active) do
      local it = a[3]
      if not seen[it] then
        seen[it] = true
        block_items[#block_items + 1] = it
      end
    end
    table.sort(block_items, function(x, y) return x.rank < y.rank end)
    local by_track, tlist, row_sets, overflow = {}, {}, {}, {}
    for _, it in ipairs(block_items) do
      local spans = {}
      for _, e in ipairs(it.extents) do
        local s0, e0 = math.max(e[1], bs), math.min(e[2], be)
        if s0 <= e0 then spans[#spans + 1] = { s0, e0 } end
      end
      local function free(r)
        local set = row_sets[r]
        if not set then return true end
        for _, sp in ipairs(spans) do
          -- 1 bp of air between neighbours so their arrowheads stay distinct
          if vf._overlaps(set, sp[1] - 1, sp[2] + 1) then return false end
        end
        return true
      end
      local row
      if it.track <= cap and free(it.track) then
        row = it.track
      else
        for r = 1, math.min(cap, #block_items) do
          if free(r) then row = r break end
        end
      end
      if row then
        row_sets[row] = row_sets[row] or {}
        for _, sp in ipairs(spans) do vf._insert(row_sets[row], sp[1], sp[2]) end
        if not by_track[row] then
          by_track[row] = {}
          tlist[#tlist + 1] = row
        end
        table.insert(by_track[row], it)
      else
        overflow[#overflow + 1] = it
      end
    end
    table.sort(tlist)

    local function cc(pos) return col0(pos - bs) + 1 end

    for _, t in ipairs(tlist) do
      local cells, hls, fids = {}, {}, {}
      for c = 1, content_w do cells[c] = " " end
      local has_cds = false
      for _, it in ipairs(by_track[t]) do
        local f, h = it.f, it.hl
        local function set(c, ch)
          cells[c], hls[c], fids[c] = ch, h, f.id
        end
        if f.type == "CDS" or f.type == "ORF" then has_cds = true end
        for _, g in ipairs(vf.parts_by_extent(it)) do
          local prev
          for _, p in ipairs(g) do
            if p.stop < p.start then
              if p.stop >= bs and p.stop <= be then set(cc(p.stop), SITE) end
            else
              local s0, e0 = math.max(p.start, bs), math.min(p.stop, be)
              if s0 <= e0 then
                for c = cc(s0), cc(e0) do set(c, BODY) end
              end
            end
            if prev then
              local s0, e0 = math.max(prev.stop + 1, bs), math.min(p.start - 1, be)
              if s0 <= e0 then
                for c = cc(s0), cc(e0) do
                  if fids[c] == nil then set(c, INTRON) end
                end
              end
            end
            if not prev or p.stop > prev.stop then prev = p end
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
      push(text, { kind = "track", b = b, bs = bs, be = be, runs = runs, segs = segs, items = by_track[t] })

      if show_aa and has_cds then
        local acells, ahls, afids = {}, {}, {}
        for c = 1, content_w do acells[c] = " " end
        for _, it in ipairs(by_track[t]) do
          if it.f.type == "CDS" or it.f.type == "ORF" then
            local map = vf.codon_map(rec, it.f)
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
        push(atext, { kind = "aa", b = b, bs = bs, be = be, runs = aruns, segs = asegs })
      end
    end
    if #overflow > 0 then
      any_overflow = true
      local counts, order = {}, {}
      for _, it in ipairs(overflow) do
        if not counts[it.f.type] then order[#order + 1] = it.f.type end
        counts[it.f.type] = (counts[it.f.type] or 0) + 1
      end
      local parts = {}
      for _, t in ipairs(order) do parts[#parts + 1] = counts[t] > 1 and (t .. " ×" .. counts[t]) or t end
      local hidden = {}
      for _, it in ipairs(overflow) do hidden[#hidden + 1] = it.f.id end
      push((" "):rep(prefix_w) .. ("+%d more: %s"):format(#overflow, table.concat(parts, ", ")),
        { kind = "overflow", b = b, bs = bs, be = be, hidden = hidden, hidden_items = overflow })
    end
    push("", { kind = "blank", b = b, bs = bs, be = be })
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
  return setmetatable(layout, { __index = M })
end

---------------------------------------------------------------------------
-- Layout queries (methods on the returned layout)
---------------------------------------------------------------------------

--- Base pair at a display column (0-based) on a row (0-based).
function M.bp_at(self, row0, vcol0)
  local info = self.rows[row0]
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
  return self.block_row[b], self.prefix_w + M.col0(self, o0)
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
  local info = self.rows[row0]
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
