-- Partial digests: every fragment an incomplete digestion can leave.
--
-- The independent answer is built by simulation rather than by the same
-- arithmetic: an incomplete digestion cuts some subset of the sites, so for
-- every subset of the Biopython cut list the complete-digest code (checked
-- against Biopython in digest_spec) is run and the products are collected.
-- The union over all 2^n subsets is what R.partial has to return.
local H = require("helpers")
local pio = require("plasmid.io")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")
local R = require("plasmid.bio.restriction")
local D = require("plasmid.digest")
local A = require("plasmid.bio.assembly")
local C = require("plasmid.clone")

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
local function load(stem) return assert(pio.read(H.fixture("ncbi/" .. stem .. ".gb"))).records[1] end
local function with_topology(rec, topo) return setmetatable({ topology = topo }, { __index = rec }) end
local function key(fr) return ("%d:%d:%s"):format(fr.start, fr.length, tostring(fr.uncut == true)) end

vim.o.columns, vim.o.lines, vim.o.hidden = 140, 50, true

-- Sites of a set, from the committed Biopython cut lists, deduplicated by
-- top-strand position the way a digest counts them.
local function biopython_cuts(by_enz, names)
  local set = {}
  for _, n in ipairs(names) do
    for _, c in ipairs(by_enz[n] or {}) do set[c] = true end
  end
  local out = {}
  for p in pairs(set) do out[#out + 1] = p end
  table.sort(out)
  return out
end

-- Every fragment reachable by cutting some subset of `positions`, as a set of
-- keys: the complete digest of each subset, unioned. 2^n subsets.
local function by_simulation(rec, positions)
  local seen = {}
  for mask = 0, 2 ^ #positions - 1 do
    local cuts = {}
    for i = 1, #positions do
      if math.floor(mask / 2 ^ (i - 1)) % 2 == 1 then cuts[#cuts + 1] = { cut = positions[i] } end
    end
    for _, fr in ipairs(R.digest(rec, cuts)) do seen[key(fr)] = true end
  end
  return seen
end

local SETS = {
  { "EcoRI" }, { "BsaI" }, { "BbsI" }, { "EcoRI", "HindIII" }, { "EcoRV", "ScaI" },
  { "PstI", "SalI" }, { "BamHI", "EcoRI", "HindIII" }, { "AatII", "NarI", "PvuII" },
  { "BbsI", "PstI", "EcoRV" },
}

---------------------------------------------------------------------------
-- Enumeration
---------------------------------------------------------------------------

case("partial_fragments_are_every_subset_digest", function()
  local expected = vim.json.decode(H.read("tests/oracle/expected/enzymes.json"))
  local stems = vim.tbl_keys(expected)
  table.sort(stems)
  local cases, frags_checked, max_sites_seen = 0, 0, 0
  -- the named sets, plus seeded random ones so that the site counts spread
  local all_names = {}
  for _, e in ipairs(require("plasmid.data.enzymes")) do all_names[#all_names + 1] = e.name end
  math.randomseed(11)
  local sets = vim.deepcopy(SETS)
  for _ = 1, 60 do
    local pick = {}
    for _ = 1, math.random(1, 3) do pick[#pick + 1] = all_names[math.random(#all_names)] end
    sets[#sets + 1] = pick
  end
  for _, stem in ipairs(stems) do
    local rec = load(stem)
    for _, topo in ipairs({ "circular", "linear" }) do
      local by_enz = expected[stem][topo]
      if by_enz then
        local r = with_topology(rec, topo)
        for _, names in ipairs(sets) do
          local pos = biopython_cuts(by_enz, names)
          -- 2^n subsets: keep the simulation affordable, and the sets varied
          if #pos <= 8 then
            local got = D.compute(r, names, { partial = true, max_fragments = 100000 })
            local want = by_simulation(r, pos)
            local seen = {}
            for _, fr in ipairs(got) do
              local k = key(fr)
              H.ok(want[k], ("%s %s %s: %s is not reachable by cutting any subset")
                :format(stem, topo, table.concat(names, "+"), k))
              H.ok(not seen[k], ("%s %s %s: %s listed twice"):format(stem, topo, table.concat(names, "+"), k))
              seen[k] = true
            end
            local missing = {}
            for k in pairs(want) do if not seen[k] then missing[#missing + 1] = k end end
            table.sort(missing)
            H.eq({}, missing, ("%s %s %s: fragments a subset gives but the partial digest does not list")
              :format(stem, topo, table.concat(names, "+")))
            cases = cases + 1
            frags_checked = frags_checked + #got
            max_sites_seen = math.max(max_sites_seen, #pos)
          end
        end
      end
    end
  end
  H.ok(cases > 100, ("only %d cases ran"):format(cases))
  H.log(("%d partial digests (%d fixtures x 2 topologies, %d enzyme sets of 1-3, up to %d sites): "
    .. "%d fragments, every one reachable by cutting a subset of the Biopython sites and none missing")
    :format(cases, #stems, #sets, max_sites_seen, frags_checked))
end)

case("the_count_is_the_square_of_the_site_count", function()
  -- n sites on a circle bound n² spans (n of them the whole molecule,
  -- linearised at one site) and leave the intact circle; n sites on a linear
  -- molecule bound C(n+2, 2) spans over its n+2 boundaries.
  local rec = load("J01749")
  local rows = {}
  for _, names in ipairs({ { "BsaI" }, { "BbsI" }, { "EcoRI", "HindIII" }, { "BamHI", "EcoRI", "HindIII" } }) do
    for _, topo in ipairs({ "circular", "linear" }) do
      local r = with_topology(rec, topo)
      local frags, _, _, _, info = D.compute(r, names, { partial = true })
      local n = info.sites
      local want = topo == "circular" and (n * n + 1) or ((n + 2) * (n + 1) / 2)
      H.eq(want, #frags, ("%s %s: %d sites"):format(table.concat(names, "+"), topo, n))
      -- the complete digest is exactly the fragments that missed no site
      local complete = D.compute(r, names)
      local none_missed = 0
      for _, f in ipairs(frags) do if f.skipped == 0 and not f.uncut then none_missed = none_missed + 1 end end
      H.eq(#complete, none_missed, ("%s %s: complete-digest fragments"):format(table.concat(names, "+"), topo))
      rows[#rows + 1] = ("%s %s %d sites → %d"):format(table.concat(names, "+"), topo:sub(1, 4), n, #frags)
    end
  end
  H.log("counts: " .. table.concat(rows, ", "))
end)

case("skipped_counts_the_sites_read_through", function()
  local rec = load("J01749")
  for _, topo in ipairs({ "circular", "linear" }) do
    local r = with_topology(rec, topo)
    local names = { "BamHI", "EcoRI", "HindIII" }
    local frags, cuts = D.compute(r, names, { partial = true })
    local pos = {}
    for _, c in ipairs(cuts) do pos[c.cut] = true end
    local L = r:len()
    for _, f in ipairs(frags) do
      -- count the cuts strictly inside the span, independently of the engine
      local inside = 0
      for p in pairs(pos) do
        local off = (p - f.start) % L
        if off > 0 and off < f.length then inside = inside + 1 end
      end
      if f.uncut then
        local n = 0
        for _ in pairs(pos) do n = n + 1 end
        H.eq(n, f.skipped, "an intact molecule missed every site")
      else
        H.eq(inside, f.skipped, ("%d..%d (%d bp)"):format(f.start, f.stop, f.length))
      end
    end
  end
  H.log("skipped: every fragment's count equals the cuts lying inside it, on both topologies")
end)

---------------------------------------------------------------------------
-- Sequence and ends
---------------------------------------------------------------------------

case("a_partial_fragment_is_its_complete_fragments_ligated", function()
  -- A span bounded by two cuts is what the complete-digest fragments between
  -- them make when they are joined again, so building it that way is an
  -- independent route to the sequence and to both ends.
  local alphabet = require("plasmid.model.alphabet")
  local checked = 0
  for _, stem in ipairs({ "J01749", "L09137", "NC_005816" }) do
    local rec = load(stem)
    for _, names in ipairs({ { "BbsI" }, { "EcoRI", "HindIII" }, { "PstI", "SalI" } }) do
      local complete = D.compute(rec, names)
      local whole = (not complete[1].uncut) and A.from_digest(rec, complete, rec.name) or {}
      local by_start = {}
      for _, f in ipairs(whole) do by_start[f.digest.start] = f end
      local parts = D.compute(rec, names, { partial = true })
      local frags = A.from_digest(rec, parts, rec.name)
      for i, f in ipairs(frags) do
        local fr = parts[i]
        if not fr.uncut then
          -- walk the complete fragments from this one's start, in order
          local chain, at, guard = {}, fr.start, 0
          while by_start[at] and guard <= #whole do
            chain[#chain + 1] = by_start[at]
            at = (at + by_start[at].digest.length - 1) % rec:len() + 1
            guard = guard + 1
            if at == ((fr.start + fr.length - 1) % rec:len() + 1) then break end
          end
          if #chain > 0 then
            local seq = A.sequence(chain[1])
            for k = 2, #chain do
              local j = assert(A.compatible(chain[k - 1], chain[k]), "neighbours must rejoin")
              seq = seq .. A.sequence(chain[k]):sub(j.len + 1)
            end
            H.eq(seq:upper(), A.sequence(f):upper(),
              ("%s %s: fragment %d (%d..%d) from %d pieces"):format(stem, table.concat(names, "+"),
                i, fr.start, fr.stop, #chain))
            H.eq(A.end_text(chain[1], "left"), A.end_text(f, "left"))
            H.eq(A.end_text(chain[#chain], "right"), A.end_text(f, "right"))
            checked = checked + 1
          end
        end
      end
    end
  end
  H.ok(checked > 40, ("only %d fragments rebuilt"):format(checked))
  H.log(("%d partial fragments equal their complete-digest pieces ligated back together, ends included")
    :format(checked))
end)

case("features_follow_a_partial_fragment", function()
  local rec = load("J01749")
  local parts = D.compute(rec, { "BamHI", "EcoRI", "HindIII" }, { partial = true })
  -- the longest read-through fragment that does not cross the origin, so that
  -- "inside the fragment" is an interval the fixture's features can be in
  local biggest
  for _, f in ipairs(parts) do
    if not f.uncut and f.skipped == 1 and f.start + f.length - 1 <= rec:len()
      and (not biggest or f.length > biggest.length) then biggest = f end
  end
  H.ok(biggest, "a read-through fragment")
  local copy = D.extract(rec, biggest)
  H.eq(biggest.length, copy:len())
  H.ok(copy.description:find("partial-digest", 1, true), copy.description)
  H.ok(copy.description:find("reading through 1 uncut site", 1, true), copy.description)
  -- every whole feature of the span is carried with the same residues
  local carried, whole = 0, 0
  for _, sf in ipairs(rec.features) do
    local lo, hi = sf:bounds()
    if sf.type ~= "source" and lo and #rec:extents(sf) == 1
      and lo >= biggest.start and hi <= biggest.start + biggest.length - 1 then
      whole = whole + 1
      local want = rec:extract(sf)
      for _, nf in ipairs(copy.features) do
        if nf.type == sf.type and nf:label() == sf:label() and copy:extract(nf) == want then
          carried = carried + 1; break
        end
      end
    end
  end
  H.ok(whole > 10, ("only %d whole features in the span"):format(whole))
  H.eq(whole, carried)
  H.log(("a %d bp read-through fragment carries all %d of its whole features"):format(biggest.length, whole))
end)

---------------------------------------------------------------------------
-- The view, the gel and the part list
---------------------------------------------------------------------------

case("the_table_says_it_is_partial_and_shows_the_ladder", function()
  local buf = open("ncbi/J01749.gb")
  local db, frags = D.run(buf, { "BbsI" }, { partial = true })
  H.eq(10, #frags, "3 sites on a circle: 3² spans plus the intact plasmid")
  local lines = vim.api.nvim_buf_get_lines(db, 0, -1, false)
  local text = table.concat(lines, "\n")
  H.ok(lines[1]:find("Partial digest of", 1, true), lines[1])
  H.ok(lines[1]:find("3 sites, 10 fragments", 1, true), lines[1])
  H.ok(text:find("missed", 1, true), "the table counts the sites read through")
  H.ok(text:find("BbsI partial", 1, true), "the gel lane is labelled partial")
  H.ok(text:find("not what one reaction yields", 1, true), "the gel says what it is")
  H.ok(D.state[db].partial, "the digest state remembers")
  -- the bands are the fragment lengths
  local sizes = {}
  for _, f in ipairs(frags) do sizes[f.length] = true end
  H.ok(sizes[4361] and sizes[863] and sizes[3614], "full length, a complete fragment and a read-through")
  vim.cmd("close")
  H.log(("partial table: %d rows, lengths %d..%d"):format(#frags, 747, 4361))
end)

case("a_partial_fragment_stages_as_a_cloning_part", function()
  C.staged = {}
  local buf = open("ncbi/J01749.gb")
  local db, frags = D.run(buf, { "BbsI" }, { partial = true })
  -- the 3614 bp read-through fragment
  local n
  for i, f in ipairs(frags) do if f.length == 3614 then n = i end end
  H.ok(n, "the read-through fragment is in the table")
  vim.api.nvim_set_current_win(vim.fn.bufwinid(db))
  vim.api.nvim_win_set_cursor(0, { D.state[db].first + n - 1, 0 })
  vim.cmd("DnaCloneAdd")
  H.eq(1, #C.staged)
  local p = C.staged[1]
  H.eq(3614 + 4, p.rec:len(), "the fragment plus its two 4 base overhangs minus the shared ends")
  H.ok(p.from:find("partial digest of", 1, true), p.from)
  H.ok(p.from:find("1 site uncut", 1, true), p.from)
  C.list(buf)
  local text = table.concat(vim.api.nvim_buf_get_lines(C.last, 0, -1, false), "\n")
  H.ok(text:find("partial digest of", 1, true), "the part list says where it came from")
  C.staged = {}
  H.log(("staged a 3614 bp read-through fragment: %s"):format(p.from))
end)

case("a_partial_lane_on_the_gel_command", function()
  local buf = open("ncbi/J01749.gb")
  local b, lanes = D.gel(buf, "BbsI | partial BbsI | uncut")
  H.eq(3, #lanes)
  H.eq(3, #lanes[1].sizes)
  H.eq(10, #lanes[2].sizes)
  H.ok(lanes[2].label:find("partial", 1, true), lanes[2].label)
  H.ok(lanes[2].note:find("every one an incomplete", 1, true), lanes[2].note)
  vim.cmd("close")
  H.log(("gel lanes: complete %d bands, partial %d bands"):format(#lanes[1].sizes, #lanes[2].sizes))
end)

case("too_many_sites_is_reported_not_enumerated", function()
  local buf = open("ncbi/L09137.gb")
  local msgs = {}
  local orig = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  -- AluI cuts pUC19 many times: the square of that is not a table anyone reads
  vim.cmd("DnaDigestPartial AluI")
  vim.notify = orig
  local m = msgs[#msgs]
  H.ok(m:find("possible fragments", 1, true), m)
  H.ok(m:find("digest.partial.max_sites", 1, true), m)
  local sites, total = m:match("with (%d+) sites has (%d+) possible")
  H.eq(tonumber(sites) ^ 2 + 1, tonumber(total), "the refusal counts what it would have listed")
  -- raising the limit lists them
  local rec = load("L09137")
  local frags, _, _, _, info = D.compute(rec, { "AluI" }, { partial = true, max_sites = 100, max_fragments = 40 })
  H.eq(40, #frags, "max_fragments caps the list")
  H.eq(true, info.capped)
  H.eq(tonumber(total), info.total)
  H.log(("AluI on pUC19: %s sites, %s possible fragments — refused above max_sites=%d")
    :format(sites, total, 12))
end)

return T