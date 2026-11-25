-- The assembly plan: what a cloner writes down before touching a pipette.
--
-- plasmid.clone produces a product record and, when nothing assembles, a
-- report of every end pair it tried. Neither is a plan. A plan says which
-- parts go in the tube and where each came from, which reaction joins them,
-- what should come out, what *else* could come out, and how to tell those
-- apart once there are colonies on the plate.
--
-- The three sections that are more than bookkeeping:
--
--   by-products  Every molecule the parts can also give. plasmid.clone only
--                reports assemblies that use each part exactly once; the
--                background a cloner actually screens against is a part left
--                out (a self-closed vector), a part the other way round, a
--                part in twice, and -- for a cut destination -- the input
--                circle carried through uncut. Those are enumerated here by
--                running the same assembler over sub- and multisets of the
--                same fragments, so a by-product is a molecule the engine
--                really builds, never a guess.
--
--   diagnostic   An enzyme set whose band pattern tells the product from each
--   digest       by-product. "Tells apart" is decided by the gel model
--                (plasmid.gel): two lanes differ only when one has a band in
--                a row the other has not. Sizes that land in the same row are
--                one band on a gel, and so are not evidence, however
--                different the two numbers look in a table.
--
--   primers      A colony-PCR pair across each new junction, sized on the
--                product and on every by-product, and a sequencing primer
--                that starts far enough back for the junction to fall in
--                readable sequence. Designed by plasmid.primerdesign and
--                sized by plasmid.pcr; no design code lives here.
--
-- Where no enzyme set discriminates, the plan says so. A check that cannot
-- tell the two outcomes apart is worse than no check, because it will be run
-- and believed.
local A = require("plasmid.bio.assembly")

local M = {}

local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end
local function cfg()
  local c = (require("plasmid.config").options.cloning or {}).plan or {}
  return vim.tbl_extend("keep", c, M.defaults)
end
local function gel_opts() return require("plasmid.config").options.gel or {} end

M.defaults = {
  max_enzymes = 2,     -- largest enzyme set the diagnostic search will offer
  min_site = 6,        -- shortest recognition site considered (a diagnostic digest wants few bands)
  candidates = nil,    -- enzyme names to search; nil = every enzyme in the active set
  max_parts = 6,       -- above this many parts, only the cheap by-product families are enumerated
  tandem = true,       -- also enumerate a part taken twice
  flank = 200,         -- colony-PCR primers this far either side of a junction
  read_offset = 150,   -- sequencing primer this far back from the junction
  slack = 8,           -- bases a plan primer may move away from its nominal position
  slack_max = 60,      -- ... and how far it may move when nothing was usable at `slack`
  min_gap = 2,         -- gel rows a distinguishing band must clear to count as resolved
  gel_lanes = 6,       -- by-product lanes drawn on the plan's gel (the rest are in the table)
  max_sets = 2000,     -- cap on enzyme sets evaluated
  top_singles = 24,    -- single enzymes carried forward into the pair search
}

---------------------------------------------------------------------------
-- Small helpers
---------------------------------------------------------------------------

local function seqof(rec) return rec.seq:tostring():upper() end
local function revcomp(s) return require("plasmid.model.alphabet").revcomp(s) end
local function wrap(p, L) return (p - 1) % L + 1 end

--- Are two records the same molecule? A circle has no start and no preferred
--- strand, so a by-product that is a rotation or the reverse complement of
--- the product is the product and must not be listed as background.
local function same_molecule(a, b)
  if a:len() ~= b:len() then return false end
  local x, y = seqof(a), seqof(b)
  if a:is_circular() ~= b:is_circular() then return false end
  if not a:is_circular() then return x == y or x == revcomp(y) end
  local xx = x .. x
  return xx:find(y, 1, true) ~= nil or xx:find(revcomp(y), 1, true) ~= nil
end

--- Part offsets and total length of an assembly, by the same arithmetic
--- plasmid.bio.assembly.build uses, so the coordinates here are the product's.
function M.offsets(asm)
  local off, pos = {}, 0
  for i, p in ipairs(asm.parts) do
    local n = p.frag.rec:len()
    if i == 1 then
      off[i] = 0
      pos = n
    else
      local ov = asm.parts[i - 1].junction.len
      off[i] = pos - ov
      pos = pos + n - ov
    end
  end
  local total = pos
  if asm.circular then total = pos - asm.parts[#asm.parts].junction.len end
  return off, total
end

---------------------------------------------------------------------------
-- Parts
---------------------------------------------------------------------------

--- How a part was prepared, in the terms a bench protocol uses: the enzyme
--- that cut it, the primers that amplified it, or that it went in whole.
local function preparation(frag, part)
  local d = frag.digest
  if d then
    local names, seen = {}, {}
    for _, side in ipairs({ "left", "right" }) do
      for _, s in ipairs(d[side] or {}) do
        if not seen[s.enzyme.name] then seen[s.enzyme.name] = true; names[#names + 1] = s.enzyme.name end
      end
    end
    if #names > 0 then return ("%s digest"):format(table.concat(names, " + ")) end
  end
  local ok, primers = pcall(function() return require("plasmid.primers").list(frag.rec) end)
  if ok and #primers >= 2 then
    local names = {}
    for i = 1, math.min(2, #primers) do names[i] = primers[i].name end
    return ("PCR with %s"):format(table.concat(names, " / "))
  end
  if part and part.rec and part.rec:is_circular() then return "whole circle" end
  return "whole record"
end

local function part_rows(ctx)
  local by_name = {}
  for _, p in ipairs(ctx.parts or {}) do by_name[p.name] = p end
  local out = {}
  for i, p in ipairs(ctx.asm.parts) do
    local staged = by_name[p.frag.name]
    out[i] = {
      n = i,
      name = p.frag.name,
      length = p.frag.rec:len(),
      from = (staged and staged.from) or p.frag.source or "—",
      prepared = preparation(p.frag, staged),
      left = A.end_text(p.frag, "left"),
      right = A.end_text(p.frag, "right"),
      flipped = p.flipped or false,
    }
  end
  return out
end

local function junction_rows(ctx)
  local asm = ctx.asm
  local off, total = M.offsets(asm)
  local out = {}
  for i, p in ipairs(asm.parts) do
    local nxt = asm.parts[i % #asm.parts + 1]
    if p.junction and (i < #asm.parts or asm.circular) then
      local stop = off[i] + p.frag.rec:len()
      out[#out + 1] = {
        n = #out + 1,
        from = p.frag.name,
        to = nxt.frag.name,
        kind = p.junction.kind,
        len = p.junction.len,
        seq = p.junction.seq,
        at = wrap(stop, total),            -- last base of part i in the product
        text = require("plasmid.clone").junction_text(p.junction),
      }
    end
  end
  return out, total
end

---------------------------------------------------------------------------
-- By-products
---------------------------------------------------------------------------

-- Every family enumerated, and the wording the plan uses for it. The list is
-- fixed and stated in the report: a cloner has to know which backgrounds were
-- considered, not only which ones were found.
local FAMILIES = {
  { key = "uncut", text = "the input circle carried through uncut" },
  { key = "self", text = "one part closing on itself (a self-ligated vector)" },
  { key = "partial", text = "a part left out" },
  { key = "rearranged", text = "the same parts in another order or orientation" },
  { key = "tandem", text = "a part taken twice (a tandem insert)" },
}

-- An arrangement written the way the assembly picker writes one: part names
-- in order, each with the strand it went in on. Spelling "(reverse)" out for
-- every part makes a tandem label three lines long and tells a reader no more.
local function arrangement(items)
  local t = {}
  for _, it in ipairs(items) do t[#t + 1] = it.frag.name .. (it.flipped and "◂" or "▸") end
  return table.concat(t, " ")
end

--- Build one candidate molecule from a list of fragments, keeping only
--- circular products when the intended assembly closed (and linear ones when
--- it did not). Returns a list of records.
local function build_set(frags, ctx, want_circular)
  local c = require("plasmid.config").options.cloning or {}
  local aopts = {
    join = ctx.method == "gibson" and "homology" or "ends",
    min_homology = c.min_homology or 20,
    max_homology = c.max_homology or 200,
    max_solutions = 8,
    circular = want_circular and true or false,
  }
  local ok, asms = pcall(A.assemble, frags, aopts)
  if not ok or not asms then return {} end
  local out = {}
  for _, asm in ipairs(asms) do
    local built, rec = pcall(A.build, asm, { annotate_parts = false, annotate_junctions = false })
    if built then out[#out + 1] = { rec = rec, asm = asm } end
  end
  return out
end

--- Every molecule these parts can give other than the intended product.
--- Returns the by-products and the families that were searched.
function M.byproducts(ctx, opts)
  opts = vim.tbl_extend("keep", opts or {}, cfg())
  local frags = ctx.frags
  local n = #frags
  local want_circular = ctx.asm.circular
  local out, searched, skipped = {}, {}, {}

  local function add(kind, label, rec, asm)
    if same_molecule(rec, ctx.product) then return false end
    for _, b in ipairs(out) do
      if same_molecule(rec, b.rec) then return false end
    end
    local k = #out + 1
    out[k] = { n = k, tag = "#" .. k, kind = kind, label = label, rec = rec, asm = asm,
      length = rec:len(), topology = rec.topology }
    return true
  end

  -- An input that was circular is still in the tube: a destination vector that
  -- was cut is never cut to completion, and the uncut circle transforms far
  -- better than a ligation. It is not an assembly product, so the assembler
  -- never reports it, but it is the background most colonies are.
  for _, p in ipairs(ctx.parts or {}) do
    if p.rec and p.rec:is_circular() then
      add("uncut", ("%s carried through uncut"):format(p.name), A.copy(p.rec))
    end
  end

  -- One part closing on itself, and -- when the set is small enough to
  -- enumerate -- every proper subset of two or more parts.
  for i = 1, n do
    for _, b in ipairs(build_set({ frags[i] }, ctx, want_circular)) do
      add("self", ("%s alone (self-closure)"):format(frags[i].name), b.rec, b.asm)
    end
  end
  if n >= 3 then
    local full = n <= opts.max_parts
    local function subsets(size)
      local idx = {}
      local function rec_(start, left)
        if left == 0 then
          local sub, names = {}, {}
          for _, i in ipairs(idx) do sub[#sub + 1] = frags[i] end
          local present = {}
          for _, i in ipairs(idx) do present[i] = true end
          for i = 1, n do if not present[i] then names[#names + 1] = frags[i].name end end
          for _, b in ipairs(build_set(sub, ctx, want_circular)) do
            add("partial", ("without %s: %s"):format(table.concat(names, " and "), arrangement(b.asm.parts)),
              b.rec, b.asm)
          end
          return
        end
        for i = start, n - left + 1 do
          idx[#idx + 1] = i
          rec_(i + 1, left - 1)
          idx[#idx] = nil
        end
      end
      rec_(1, size)
    end
    if full then
      for size = 2, n - 1 do subsets(size) end
    else
      subsets(n - 1)
      skipped[#skipped + 1] = ("subsets smaller than %d parts (%d parts staged; cloning.plan.max_parts is %d)")
        :format(n - 1, n, opts.max_parts)
    end
  end

  -- The same parts in another arrangement. The assembler already enumerated
  -- these when the product was built; anything it found other than the
  -- product is background.
  for _, b in ipairs(build_set(frags, ctx, want_circular)) do
    -- Two parts have only one alternative arrangement, and it is the one every
    -- cloner screens for first, so it is named rather than drawn.
    local label = ("another arrangement: %s"):format(arrangement(b.asm.parts))
    if n == 2 and #b.asm.parts == 2 then
      -- A circular search fixes the first fragment, so both arrangements list
      -- the same part second and differ only in which way round it went in.
      local other, intended = b.asm.parts[2], ctx.asm.parts[2]
      if other.index == intended.index and other.flipped ~= intended.flipped then
        label = ("%s in the reverse orientation"):format(other.frag.name)
      end
    end
    add("rearranged", label, b.rec, b.asm)
  end

  -- A part taken twice. Ligation and Gibson both do this readily when the two
  -- ends of a part are compatible with the same partner.
  if opts.tandem and n <= opts.max_parts then
    for i = 1, n do
      local dup = {}
      for k = 1, n do dup[#dup + 1] = frags[k] end
      dup[#dup + 1] = frags[i]
      for _, b in ipairs(build_set(dup, ctx, want_circular)) do
        add("tandem", ("%s twice: %s"):format(frags[i].name, arrangement(b.asm.parts)), b.rec, b.asm)
      end
    end
  elseif opts.tandem then
    skipped[#skipped + 1] = ("tandem products (%d parts staged; cloning.plan.max_parts is %d)"):format(n, opts.max_parts)
  end

  for _, f in ipairs(FAMILIES) do
    if f.key ~= "tandem" or opts.tandem then searched[#searched + 1] = f.text end
  end
  return out, searched, skipped
end

---------------------------------------------------------------------------
-- The diagnostic digest
---------------------------------------------------------------------------

--- Candidate enzymes: the active set, cut down to the ones worth running on a
--- diagnostic gel. A site of fewer than `min_site` specific bases cuts a
--- plasmid into a smear, and an enzyme that needs methylated DNA cannot be
--- relied on at all.
function M.candidate_enzymes(opts)
  local R = require("plasmid.bio.restriction")
  opts = vim.tbl_extend("keep", opts or {}, cfg())
  local out = {}
  if opts.candidates then
    for _, name in ipairs(opts.candidates) do
      local e = R.get(name)
      if e then out[#out + 1] = e end
    end
    return out
  end
  for _, e in ipairs(R.enzymes()) do
    local specific = select(2, e.site:gsub("[ACGTacgt]", ""))
    if specific >= opts.min_site and not R.NEEDS_METHYLATION[e.name] then out[#out + 1] = e end
  end
  return out
end

--- Fragment lengths of one molecule under one enzyme set, from sites scanned
--- once per molecule. Returns the sizes and whether the set cut at all.
local function pattern(rec, sites, names)
  local R = require("plasmid.bio.restriction")
  local cuts = {}
  for _, n in ipairs(names) do
    for _, s in ipairs(sites[n] or {}) do cuts[#cuts + 1] = s end
  end
  if #cuts == 0 then return { rec:len() }, false end
  local sizes = {}
  for i, f in ipairs(R.digest(rec, cuts)) do sizes[i] = f.length end
  return sizes, true
end

--- Which gel rows each lane occupies, asked of the gel model rather than of a
--- rule invented here: `placed` is where a band actually lands, `off` is what
--- ran off the gel or sat above the ladder and so shows nothing.
--- Returns rows[lane] = { [row] = a size in it }, offs[lane], and the gel info.
function M.lane_rows(lanes, opts)
  local _, _, info = require("plasmid.gel").render(lanes, opts or gel_opts())
  local rows, offs = {}, {}
  for i = 1, #lanes do
    local set = {}
    for _, b in ipairs(info.placed[i] or {}) do set[b[2]] = b[1] end
    rows[i] = set
    offs[i] = info.off[i] or {}
  end
  return rows, offs, info
end

--- The bands that tell lane `b` from lane `a`, with how far each sits from the
--- nearest band in the other lane (rows). A difference of one row is a
--- difference the gel model can draw but a gel would struggle to read, so the
--- distance is what ranks one candidate set above another.
local function differences(rows, a, b)
  local out = {}
  local function scan(from, to, which)
    for r, size in pairs(rows[from]) do
      if rows[to][r] == nil then
        local gap = math.huge
        for r2 in pairs(rows[to]) do gap = math.min(gap, math.abs(r - r2)) end
        out[#out + 1] = { size = size, row = r, gap = gap, lane = which }
      end
    end
  end
  scan(a, b, "product")
  scan(b, a, "byproduct")
  table.sort(out, function(x, y)
    if x.gap ~= y.gap then return x.gap > y.gap end
    return x.size > y.size
  end)
  return out
end

--- Does an enzyme set tell the product from every by-product? Returns a
--- verdict table: per by-product the distinguishing bands (or nil), the number
--- told apart, and the weakest separation among them.
function M.evaluate(names, product, byproducts, sites_by_mol, gopts)
  local lanes = { { label = "product", sizes = (pattern(product, sites_by_mol[1], names)) } }
  local cut = { (select(2, pattern(product, sites_by_mol[1], names))) }
  for i, b in ipairs(byproducts) do
    local sizes, did = pattern(b.rec, sites_by_mol[i + 1], names)
    cut[i + 1] = did
    -- the lane is labelled with the by-product's number: a gel wide enough for
    -- a dozen full labels is a gel nobody can read
    lanes[i + 1] = { label = b.tag or ("#" .. i), sizes = sizes }
  end
  local rows = M.lane_rows(lanes, gopts)
  local per, told, weakest = {}, 0, math.huge
  for i = 1, #byproducts do
    local diff = differences(rows, 1, i + 1)
    if #diff > 0 then
      told = told + 1
      per[i] = { bands = diff, uncut = not cut[i + 1] }
      weakest = math.min(weakest, diff[1].gap == math.huge and 99 or diff[1].gap)
    else
      per[i] = nil
    end
  end
  local bands = 0
  for _ in pairs(rows[1]) do bands = bands + 1 end
  weakest = weakest == math.huge and 99 or weakest
  local min_gap = cfg().min_gap
  return { names = names, lanes = lanes, per = per, told = told, all = told == #byproducts,
    weakest = weakest, clear = told == #byproducts and weakest >= min_gap,
    bands = bands, cuts_product = cut[1],
    uncut_lanes = (function()
      local t = {}
      for i = 2, #cut do if not cut[i] then t[#t + 1] = byproducts[i - 1].label end end
      return t
    end)() }
end

--- The lanes to draw for a verdict. By-products whose digest gives the same
--- sizes give the same lane, so they are drawn once under both numbers: a gel
--- with a dozen identical lanes is wide and says nothing a reader cannot see
--- from one of them.
--- Returns the lanes and how many by-products are not drawn.
function M.gel_lanes(verdict, byproducts, limit)
  limit = limit or cfg().gel_lanes
  local out, index, dropped = { verdict.lanes[1] }, {}, 0
  for i = 2, #verdict.lanes do
    local key = table.concat(verdict.lanes[i].sizes, ",")
    local at = index[key]
    if at then
      out[at].label = out[at].label .. "," .. tostring(byproducts[i - 1].n)
    elseif #out <= limit then
      index[key] = #out + 1
      out[#out + 1] = { label = verdict.lanes[i].label, sizes = verdict.lanes[i].sizes }
    else
      dropped = dropped + 1
    end
  end
  return out, dropped
end

-- Which of two candidate sets a cloner would rather run. Answering the
-- question at all comes first, then whether every distinguishing band clears
-- `min_gap` rows -- a band one row from its neighbour is a band the gel model
-- can draw and a gel cannot read, and a second enzyme is a cheaper price to
-- pay than a check that is read wrong. Only then does fewness of enzymes
-- decide, and only among sets that pass the separation gate does the number
-- of bands matter.
local function better(a, b)
  if not b then return true end
  if a.all ~= b.all then return a.all end
  -- A set that leaves the product uncut has no band pattern to read: the one
  -- lane it gives is a supercoiled circle, which does not run at its linear
  -- size. "The product does not cut" is a real Golden Gate check, but it is
  -- not a diagnostic digest, so a set that cuts the product wins.
  if a.cuts_product ~= b.cuts_product then return a.cuts_product end
  if a.clear ~= b.clear then return a.clear end
  if #a.names ~= #b.names then return #a.names < #b.names end
  if a.told ~= b.told then return a.told > b.told end
  if #a.uncut_lanes ~= #b.uncut_lanes then return #a.uncut_lanes < #b.uncut_lanes end
  if a.clear then
    if a.bands ~= b.bands then return a.bands < b.bands end
    if a.weakest ~= b.weakest then return a.weakest > b.weakest end
  else
    if a.weakest ~= b.weakest then return a.weakest > b.weakest end
    if a.bands ~= b.bands then return a.bands < b.bands end
  end
  return table.concat(a.names, "+") < table.concat(b.names, "+")
end

--- Search enzyme sets for one that tells the product from every by-product.
--- Single enzymes first, then pairs of the ones that discriminated something,
--- then triples if `max_enzymes` allows; the first size that fully
--- discriminates wins, so a one-enzyme answer is never passed over for a
--- prettier two-enzyme one. Returns the best verdict, the best fully
--- discriminating verdict (or nil) and a summary of the search.
function M.diagnostic(product, byproducts, opts)
  opts = vim.tbl_extend("keep", opts or {}, cfg())
  local R = require("plasmid.bio.restriction")
  local enzymes = M.candidate_enzymes(opts)
  local mols = { product }
  for _, b in ipairs(byproducts) do mols[#mols + 1] = b.rec end
  local sites_by_mol = {}
  for i, rec in ipairs(mols) do
    local P = R.prepare(rec)
    local s = {}
    for _, e in ipairs(enzymes) do s[e.name] = R.sites(P, e) end
    sites_by_mol[i] = s
  end

  local tried, best, best_full = 0, nil, nil
  local singles = {}
  if #byproducts == 0 then
    return nil, nil, { tried = 0, candidates = #enzymes, reason = "no by-product to tell the product from" }
  end
  for _, e in ipairs(enzymes) do
    -- an enzyme that cuts nothing in the tube says nothing about what is in it
    local any = false
    for i = 1, #mols do if #(sites_by_mol[i][e.name] or {}) > 0 then any = true break end end
    if any then
      tried = tried + 1
      local v = M.evaluate({ e.name }, product, byproducts, sites_by_mol, gel_opts())
      singles[#singles + 1] = v
      if better(v, best) then best = v end
      if v.all and better(v, best_full) then best_full = v end
    end
  end

  -- A single enzyme that discriminates only by a band its neighbour nearly
  -- touches is not the end of the search: a pair that discriminates cleanly is
  -- worth the second enzyme.
  local function settled() return best_full and best_full.clear end
  if not settled() and opts.max_enzymes >= 2 then
    table.sort(singles, function(a, b) return better(a, b) end)
    local top = {}
    for i = 1, math.min(opts.top_singles, #singles) do top[i] = singles[i].names[1] end
    for i = 1, #top do
      for j = i + 1, #top do
        if tried < opts.max_sets then
          tried = tried + 1
          local v = M.evaluate({ top[i], top[j] }, product, byproducts, sites_by_mol, gel_opts())
          if better(v, best) then best = v end
          if v.all and better(v, best_full) then best_full = v end
        end
      end
    end
    if not settled() and opts.max_enzymes >= 3 then
      local k = math.min(10, #top)
      for i = 1, k do
        for j = i + 1, k do
          for l = j + 1, k do
            if tried < opts.max_sets then
              tried = tried + 1
              local v = M.evaluate({ top[i], top[j], top[l] }, product, byproducts, sites_by_mol, gel_opts())
              if better(v, best) then best = v end
              if v.all and better(v, best_full) then best_full = v end
            end
          end
        end
      end
    end
  end

  -- Which by-products nothing in the search could tell from the product: the
  -- honest answer when the search declines.
  local never = {}
  for i, b in ipairs(byproducts) do
    local any = false
    for _, v in ipairs(singles) do if v.per[i] then any = true break end end
    if not any and not (best and best.per[i]) then never[#never + 1] = b end
  end
  return best_full, best, { tried = tried, candidates = #enzymes, never = never,
    max_enzymes = opts.max_enzymes }
end

---------------------------------------------------------------------------
-- Colony PCR and sequencing primers
---------------------------------------------------------------------------

--- A primer pair for a window, asking primerdesign twice: once with its own
--- modest slack, and -- only if nothing was usable there -- once with room to
--- move. Where exactly a check primer sits does not matter, only that it spans
--- or reaches the junction, but the search cost is quadratic in the slack, so
--- paying for it on every junction would make the plan slow for no gain.
local function design(product, sel, opts)
  local PD = require("plasmid.primerdesign")
  local list, err = PD.candidates(product, sel, { slack = opts.slack })
  if not list then list, err = PD.candidates(product, sel, { slack = opts.slack_max }) end
  return list, err
end

--- Sizes a primer pair gives on a template, smallest first.
local function pcr_sizes(rec, fwd, rev)
  local PCR = require("plasmid.pcr")
  local binding = (require("plasmid.config").options.primers or {}).binding
  local ok, prods = pcall(PCR.products, seqof(rec), rec:is_circular(), fwd, rev, { binding = binding })
  if not ok or not prods then return {} end
  local out = {}
  for _, p in ipairs(prods) do out[#out + 1] = p.length end
  table.sort(out)
  return out
end

--- Can a gel tell these two sets of product sizes apart? The same question
--- the diagnostic digest asks, asked of two PCR lanes.
function M.resolvable(a, b, opts)
  if #a == 0 or #b == 0 then return #a ~= #b end
  local rows = M.lane_rows({ { label = "a", sizes = a }, { label = "b", sizes = b } }, opts or gel_opts())
  for r in pairs(rows[1]) do if rows[2][r] == nil then return true end end
  for r in pairs(rows[2]) do if rows[1][r] == nil then return true end end
  return false
end

local function sizes_text(list)
  if #list == 0 then return "no product" end
  local t = {}
  for i, s in ipairs(list) do t[i] = tostring(s) end
  return table.concat(t, ", ") .. " bp"
end

--- A colony-PCR pair across every junction, sized on the product and on each
--- by-product. The design is primerdesign's: the selection handed to it is
--- the window around the junction, so the pair it returns necessarily spans
--- the junction, and its own scoring decides the rest.
function M.colony_pcr(product, junctions, byproducts, opts)
  opts = vim.tbl_extend("keep", opts or {}, cfg())
  local L = product:len()
  local out = {}
  for _, j in ipairs(junctions) do
    local row = { junction = j.n, at = j.at, from = j.from, to = j.to }
    local flank = math.min(opts.flank, math.floor((L - 1) / 2))
    local sel
    if product:is_circular() then
      sel = { start = wrap(j.at - flank, L), stop = wrap(j.at + flank, L) }
    else
      sel = { start = math.max(1, j.at - flank), stop = math.min(L, j.at + flank) }
    end
    local list, err = design(product, sel, opts)
    if not list then
      row.error = err
    else
      local pair = list[1]
      row.pair = pair
      row.fwd, row.rev = pair.fwd.seq, pair.rev.seq
      row.tm = { pair.tm1, pair.tm2 }
      row.ta = pair.ta
      row.product_sizes = pcr_sizes(product, row.fwd, row.rev)
      row.on = {}
      for i, b in ipairs(byproducts) do
        local sizes = pcr_sizes(b.rec, row.fwd, row.rev)
        row.on[i] = { sizes = sizes, tells = M.resolvable(row.product_sizes, sizes) }
      end
    end
    out[#out + 1] = row
  end
  return out
end

--- A sequencing primer on each side of each junction, placed far enough back
--- that the junction is not in the first bases of the read, where the trace is
--- unreadable. The primer itself is primerdesign's pick for a window whose
--- edge is where the read should start.
function M.sequencing(product, junctions, opts)
  opts = vim.tbl_extend("keep", opts or {}, cfg())
  local L = product:len()
  local back = math.min(opts.read_offset, math.floor((L - 1) / 3))
  local out = {}
  local function window(a, b)
    if product:is_circular() then return { start = wrap(a, L), stop = wrap(b, L) } end
    return { start = math.max(1, a), stop = math.min(L, b) }
  end
  for _, j in ipairs(junctions) do
    local row = { junction = j.n, at = j.at, from = j.from, to = j.to }
    local fl = design(product, window(j.at - back, j.at + 60), opts)
    if fl then
      local p = fl[1].fwd
      row.fwd = { seq = p.seq, start = p.start, stop = p.stop, tm = p.tm_rounded or math.floor(p.tm + 0.5),
        reads_from = wrap(p.stop + 1, L), into = (j.at - p.stop) % L }
    end
    local rl = design(product, window(j.at - 60, j.at + back), opts)
    if rl then
      local p = rl[1].rev
      row.rev = { seq = p.seq, start = p.start, stop = p.stop, tm = p.tm_rounded or math.floor(p.tm + 0.5),
        reads_from = wrap(p.start - 1, L), into = (p.start - 1 - j.at) % L }
    end
    if not row.fwd and not row.rev then row.error = "no usable sequencing primer either side of this junction" end
    out[#out + 1] = row
  end
  return out
end

---------------------------------------------------------------------------
-- The plan
---------------------------------------------------------------------------

--- Everything the plan contains, from a finished assembly context (the one
--- plasmid.clone records as `last_assembly`).
function M.build(ctx, opts)
  opts = vim.tbl_extend("keep", opts or {}, cfg())
  assert(ctx and ctx.asm and ctx.product, "there is no finished assembly to plan")
  local C = require("plasmid.clone")
  local junctions, total = junction_rows(ctx)
  local plan = {
    method = ctx.method,
    enzyme = ctx.enzyme,
    label = C.METHODS[ctx.method] or ctx.method,
    product = ctx.product,
    length = ctx.product:len(),
    total = total,
    topology = ctx.product.topology,
    parts = part_rows(ctx),
    junctions = junctions,
    warnings = ctx.warnings or {},
    nsolutions = ctx.nsolutions,
    src = ctx.src,
  }
  plan.byproducts, plan.searched, plan.skipped = M.byproducts(ctx, opts)
  plan.diagnostic, plan.best, plan.search = M.diagnostic(ctx.product, plan.byproducts, opts)
  plan.pcr = M.colony_pcr(ctx.product, junctions, plan.byproducts, opts)
  plan.sequencing = M.sequencing(ctx.product, junctions, opts)
  return plan
end

---------------------------------------------------------------------------
-- Rendering
---------------------------------------------------------------------------

local function pad(s, w) return s .. (" "):rep(math.max(0, w - vim.fn.strdisplaywidth(s))) end
local function rpad(s, w) return (" "):rep(math.max(0, w - vim.fn.strdisplaywidth(s))) .. s end

--- Rows of a small table. `header` may be nil for a key/value block that needs
--- no column names: an empty header row is a line of trailing spaces, which
--- the written file would carry and the report buffer would not.
local function table_rows(header, rows, right)
  local n = header and #header or 0
  for _, r in ipairs(rows) do n = math.max(n, #r) end
  local w = {}
  for c = 1, n do
    w[c] = header and vim.fn.strdisplaywidth(header[c] or "") or 0
    for _, r in ipairs(rows) do w[c] = math.max(w[c], vim.fn.strdisplaywidth(tostring(r[c] or ""))) end
  end
  local function fmt(r)
    local out = {}
    for c = 1, n do
      local cell = tostring(r[c] or "")
      out[c] = (right and right[c]) and rpad(cell, w[c]) or (c < n and pad(cell, w[c]) or cell)
    end
    return "  " .. table.concat(out, "  ")
  end
  local lines = {}
  if header then lines[1] = fmt(header) end
  for _, r in ipairs(rows) do lines[#lines + 1] = fmt(r) end
  return lines
end

--- A by-product label clipped to the width of a column. The number in front
--- of it is what the gel and the other tables refer to, so a clipped label
--- still resolves.
local function clip(text, n)
  n = n or 54
  if vim.fn.strdisplaywidth(text) <= n then return text end
  return vim.fn.strcharpart(text, 0, n - 1) .. "…"
end

--- The plan as buffer lines plus highlights, in the style of the assembly
--- report: a title, then named sections of small tables.
function M.lines(plan)
  local lines, hl = {}, {}
  local function head(text)
    lines[#lines + 1] = ""
    lines[#lines + 1] = " " .. text
    hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "PlasmidHeader" }
  end
  local function put(text) lines[#lines + 1] = text end
  local function rows(header, body, right)
    for _, l in ipairs(table_rows(header, body, right)) do lines[#lines + 1] = l end
  end

  local label = plan.label
  if plan.method == "goldengate" and plan.enzyme then label = ("%s (%s)"):format(label, plan.enzyme) end
  lines[1] = (" Assembly plan — %s: %d part%s → %s, %d bp %s"):format(label, #plan.parts,
    #plan.parts == 1 and "" or "s", plan.product.name, plan.length, plan.topology)
  hl[1] = { 0, 0, #lines[1], "PlasmidTitle" }

  head("reaction")
  local rrows = {
    { "method", label },
    { "parts", tostring(#plan.parts) },
    { "product", ("%s, %d bp %s"):format(plan.product.name, plan.length, plan.topology) },
  }
  if plan.method == "gibson" then
    local c = require("plasmid.config").options.cloning or {}
    rrows[#rrows + 1] = { "joins", ("terminal homology, at least %d bp"):format(c.min_homology or 20) }
  elseif plan.method == "goldengate" then
    rrows[#rrows + 1] = { "joins", ("%s overhangs, one pot"):format(plan.enzyme) }
  else
    rrows[#rrows + 1] = { "joins", "ligase, by the ends the parts already have" }
  end
  rows(nil, rrows)

  head("parts")
  local body = {}
  for _, p in ipairs(plan.parts) do
    body[#body + 1] = { tostring(p.n), p.name, ("%d bp"):format(p.length), p.from, p.prepared,
      p.flipped and "reverse" or "forward", p.left, p.right }
  end
  rows({ "#", "part", "length", "from", "prepared by", "strand", "left end", "right end" }, body,
    { true, false, true })

  head("expected product")
  put(("  %s, %d bp %s"):format(plan.product.name, plan.length, plan.topology))
  if #plan.junctions == 0 then
    put("  a single part, no new junction")
  else
    body = {}
    for _, j in ipairs(plan.junctions) do
      body[#body + 1] = { tostring(j.n), ("%s → %s"):format(j.from, j.to), j.text, tostring(j.at) }
    end
    rows({ "#", "from → to", "join", "at bp" }, body, { true, false, false, true })
  end

  head(("by-products (%d)"):format(#plan.byproducts))
  put("  searched for:")
  for _, s in ipairs(plan.searched or {}) do put("    " .. s) end
  for _, s in ipairs(plan.skipped or {}) do put("  not enumerated: " .. s) end
  if #plan.byproducts == 0 then
    put("  none: these parts give nothing else")
  else
    body = {}
    for _, b in ipairs(plan.byproducts) do
      body[#body + 1] = { tostring(b.n), clip(b.label), ("%d bp"):format(b.length), b.topology }
    end
    rows({ "#", "what", "length", "topology" }, body, { true, false, true })
  end

  head("diagnostic digest")
  local d = plan.diagnostic
  local s = plan.search or {}
  if #plan.byproducts == 0 then
    put("  nothing to tell the product from")
  elseif d then
    put(("  %s — tells the product from %s"):format(table.concat(d.names, " + "),
      #plan.byproducts == 1 and "the one by-product" or ("all %d by-products"):format(#plan.byproducts)))
    if not d.clear then
      put(("  the bands that do it are within %d gel row%s of a neighbour: read this gel carefully, or run the PCR")
        :format(d.weakest, d.weakest == 1 and "" or "s"))
    end
    put(("  searched %d enzyme sets built from %d enzymes"):format(s.tried or 0, s.candidates or 0))
    body = {}
    for i, b in ipairs(plan.byproducts) do
      local per = d.per[i]
      local band = per and per.bands[1]
      body[#body + 1] = { tostring(b.n), clip(b.label),
        band and (("%d bp band, %d row%s clear"):format(band.size, band.gap == math.huge and 99 or band.gap,
          band.gap == 1 and "" or "s")) or "—",
        band and (band.lane == "product" and "only the product has it" or "only the by-product has it") or "—" }
    end
    rows({ "#", "by-product", "band that tells it apart", "which lane" }, body, { true })
    if not d.cuts_product then
      put("  this set does not cut the product at all: the product lane is an uncut circle, which does")
      put("  not run at its linear size, so the check is the absence of bands and not their position")
    end
    for _, u in ipairs(d.uncut_lanes or {}) do
      put(("  %s is not cut by this set: that lane is an uncut circle, which does not run at its linear size")
        :format(clip(u, 70)))
    end
  else
    put(("  no set of up to %d enzymes tells the product from every by-product"):format(s.max_enzymes or 2))
    put(("  searched %d enzyme sets built from %d enzymes"):format(s.tried or 0, s.candidates or 0))
    for _, b in ipairs(s.never or {}) do
      put(("  nothing distinguishes #%d %s (%d bp)"):format(b.n, clip(b.label), b.length))
    end
    if plan.best and plan.best.told > 0 then
      put(("  the best partial check is %s, which tells %d of %d apart; it is not a check for the rest")
        :format(table.concat(plan.best.names, " + "), plan.best.told, #plan.byproducts))
    end
    put("  use the colony PCR or sequencing below instead; a digest that cannot tell two outcomes")
    put("  apart will be run and believed")
  end

  local show = d or plan.best
  if show then
    head(("gel — %s"):format(table.concat(show.names, " + ")))
    put("  lanes are numbered as the by-product table above; lanes with the same bands are drawn once")
    local glanes, dropped = M.gel_lanes(show, plan.byproducts)
    if dropped > 0 then
      put(("  %d further by-product lane%s not drawn: the bands that tell each one apart are in the table above")
        :format(dropped, dropped == 1 and " is" or "s are"))
    end
    local glines = require("plasmid.gel").render(glanes, gel_opts())
    for _, l in ipairs(glines) do lines[#lines + 1] = l end
  end

  head("colony PCR")
  if #plan.pcr == 0 then
    put("  no junction to check")
  end
  for _, r in ipairs(plan.pcr) do
    put(("  junction %d (%s → %s, at %d bp)"):format(r.junction, r.from, r.to, r.at))
    if r.error then
      put("    no pair: " .. tostring(r.error))
    else
      put(("    fwd  5'-%s-3'  %d nt, Tm %d °C"):format(r.fwd, #r.fwd, r.tm[1]))
      put(("    rev  5'-%s-3'  %d nt, Tm %d °C"):format(r.rev, #r.rev, r.tm[2]))
      put(("    Ta %d °C"):format(r.ta))
      body = { { "", "intended product", sizes_text(r.product_sizes), "—" } }
      for i, b in ipairs(plan.byproducts) do
        body[#body + 1] = { tostring(b.n), clip(b.label), sizes_text(r.on[i].sizes),
          r.on[i].tells and "tells it apart" or "same size — not a check" }
      end
      rows({ "#", "template", "product", "verdict" }, body, { true })
    end
  end

  head("sequencing primers")
  for _, r in ipairs(plan.sequencing) do
    put(("  junction %d (%s → %s, at %d bp)"):format(r.junction, r.from, r.to, r.at))
    if r.error then put("    " .. r.error) end
    if r.fwd then
      put(("    forward  5'-%s-3'  Tm %d °C, reads from %d bp; the junction is %d bases into the read")
        :format(r.fwd.seq, r.fwd.tm, r.fwd.reads_from, r.fwd.into))
    end
    if r.rev then
      put(("    reverse  5'-%s-3'  Tm %d °C, reads back from %d bp; the junction is %d bases into the read")
        :format(r.rev.seq, r.rev.tm, r.rev.reads_from, r.rev.into))
    end
  end

  if #(plan.warnings or {}) > 0 then
    head("warnings from the assembly")
    for _, w in ipairs(plan.warnings) do
      put("  " .. w)
      hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "WarningMsg" }
    end
  end

  head("what this plan does not model")
  put("  the gel is the plugin's own model: bands are placed on a log axis next to a")
  put("  ladder, and two sizes in the same row are one band. Intensity is not modelled,")
  put("  so a doubled band is not read as evidence. Partial digestion, star activity,")
  put("  supercoiling and the transformation efficiency of each molecule are not modelled")
  put("  either; the by-products are what the parts can give, not what they will give.")

  put("")
  put(" q close")
  hl[#hl + 1] = { #lines - 1, 0, #lines[#lines], "Comment" }
  return lines, hl
end

--- The same content as plain text, for the file.
function M.text(plan)
  local lines = M.lines(plan)
  local out = {}
  for i, l in ipairs(lines) do
    -- the buffer's last line is its key hint, which a file has no use for
    if not (i == #lines and l == " q close") and not (i == #lines - 1 and lines[#lines] == " q close" and l == "") then
      out[#out + 1] = (l:gsub("%s+$", ""))
    end
  end
  return table.concat(out, "\n") .. "\n"
end

function M.write(plan, path)
  local fh, err = io.open(path, "wb")
  if not fh then return nil, tostring(err) end
  fh:write(M.text(plan))
  fh:close()
  return path
end

---------------------------------------------------------------------------
-- The command
---------------------------------------------------------------------------

M.state = {}    -- report buffer -> { src }
local ns = vim.api.nvim_create_namespace("plasmid_cloneplan")

local function report_buffer(src)
  for b, s in pairs(M.state) do
    if s.src == src and vim.api.nvim_buf_is_valid(b) then return b end
  end
  local b = vim.api.nvim_create_buf(false, true)
  vim.bo[b].bufhidden = "wipe"
  vim.bo[b].filetype = "plasmid_cloneplan"
  pcall(vim.api.nvim_buf_set_name, b, "plasmid://cloneplan/" .. tostring(src or 0))
  require("plasmid.util").case_insensitive_search(b)
  vim.api.nvim_create_autocmd("BufWipeout", { buffer = b, once = true, callback = function() M.state[b] = nil end })
  return b
end

--- Show the plan in a split next to its source buffer.
function M.show(src, lines, hl)
  local b = report_buffer(src)
  vim.bo[b].modifiable = true
  vim.api.nvim_buf_set_lines(b, 0, -1, false, lines)
  vim.bo[b].modifiable = false
  vim.api.nvim_buf_clear_namespace(b, ns, 0, -1)
  for _, h in ipairs(hl or {}) do
    pcall(vim.api.nvim_buf_set_extmark, b, ns, h[1], h[2], { end_col = h[3], hl_group = h[4] })
  end
  M.state[b] = { src = src }
  M.last = b
  local win = vim.fn.bufwinid(b)
  if win == -1 then
    vim.cmd("botright " .. math.min(#lines + 1, math.floor(vim.o.lines * 0.6)) .. "split")
    vim.api.nvim_win_set_buf(0, b)
    win = vim.api.nvim_get_current_win()
    vim.wo[win].wrap = false
    vim.wo[win].number = false
    vim.wo[win].relativenumber = false
  end
  vim.api.nvim_win_set_cursor(win, { 1, 0 })
  vim.keymap.set("n", "q", "<cmd>close<cr>", { buffer = b, nowait = true, silent = true })
  return b
end

--- The file the plan is written to when no path is given: beside the source
--- record, named after the product.
function M.default_path(plan)
  local dir = vim.fn.getcwd()
  local st = plan.src and require("plasmid.buffer").get(plan.src)
  if st and st.path then dir = vim.fn.fnamemodify(st.path, ":h") end
  return ("%s/%s_plan.txt"):format(dir, (plan.product.name:gsub("[^%w%._%-]", "_")))
end

--- :DnaClonePlan [path] — the plan for the last assembly, in a split and
--- (with a path, or with !) in a file.
function M.run(buf, path, bang)
  local ctx = require("plasmid.clone").last_assembly
  if not ctx or not ctx.product then
    return notify("no finished assembly to plan: run :DnaClone, :DnaLigate, :DnaGibson or :DnaGoldenGate first",
      vim.log.levels.WARN)
  end
  local ok, plan = pcall(M.build, ctx, {})
  if not ok then
    return notify("could not build the plan: " .. tostring(plan), vim.log.levels.ERROR)
  end
  local lines, hl = M.lines(plan)
  local src = plan.src or (require("plasmid.buffer").get(buf) and buf) or nil
  M.show(src, lines, hl)
  local written
  if (path and path ~= "") or bang then
    local p = (path and path ~= "") and vim.fn.expand(path) or M.default_path(plan)
    local wrote, err = M.write(plan, p)
    if not wrote then
      notify("could not write " .. p .. ": " .. tostring(err), vim.log.levels.ERROR)
    else
      written = p
      notify(("plan written to %s (%d by-products, %s)"):format(vim.fn.fnamemodify(p, ":~:."),
        #plan.byproducts, plan.diagnostic and ("diagnostic digest " .. table.concat(plan.diagnostic.names, " + "))
          or "no discriminating digest"))
    end
  end
  return plan, written
end

return M
