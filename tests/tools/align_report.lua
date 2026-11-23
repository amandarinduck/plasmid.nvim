-- Numbers for docs/alignment_report.md. Not part of the test suite.
--
--   /opt/nvim-linux-x86_64/bin/nvim -l tests/tools/align_report.lua [section]
--
-- Sections: oracle (agreement with Biopython), band (banded vs full matrix),
-- reads (simulated-read recovery), timing (speed and memory). Default: all.
-- Everything reads the committed JSON; no network, no Python.
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({
  root .. "/lua/?.lua", root .. "/lua/?/init.lua", root .. "/tests/?.lua", package.path,
}, ";")
vim.opt.runtimepath:prepend(root)
vim.notify = function() end
require("plasmid").setup({})

local A = require("plasmid.bio.align")
local pio = require("plasmid.io")
local alphabet = require("plasmid.model.alphabet")
local O = dofile(root .. "/tests/tools/align_oracle.lua")
local R = dofile(root .. "/tests/tools/align_reads.lua")

local function out(...) io.stdout:write(string.format(...), "\n") end
local want = arg[1]
local function want_section(name) return not want or want == name or (want == "limit" and name == "timing") end

---------------------------------------------------------------------------
-- Agreement with Biopython
---------------------------------------------------------------------------

local doc = O.load(root)

if want_section("oracle") then
  out("## corpus")
  local fams, modes, schemes = {}, {}, {}
  for _, c in ipairs(doc.cases) do
    fams[c.family] = (fams[c.family] or 0) + 1
    modes[c.mode] = (modes[c.mode] or 0) + 1
    schemes[c.scheme] = (schemes[c.scheme] or 0) + 1
  end
  out("cases\t%d\tbiopython\t%s", #doc.cases, doc.biopython)
  for _, t in ipairs({ { "family", fams }, { "mode", modes }, { "scheme", schemes } }) do
    local ks = {}
    for k in pairs(t[2]) do ks[#ks + 1] = k end
    table.sort(ks)
    for _, k in ipairs(ks) do out("%s\t%s\t%d", t[1], k, t[2][k]) end
  end

  local t0 = vim.uv.hrtime()
  local counts, bad, by_family = O.run(doc)
  out("replay_ms\t%.0f", (vim.uv.hrtime() - t0) / 1e6)
  out("\n## agreement")
  local ks = {}
  for k in pairs(counts) do ks[#ks + 1] = k end
  table.sort(ks)
  for _, k in ipairs(ks) do out("%s\t%d", k, counts[k]) end
  out("\n## agreement by family")
  local fs = {}
  for k in pairs(by_family) do fs[#fs + 1] = k end
  table.sort(fs)
  for _, f in ipairs(fs) do
    local c = by_family[f]
    out("%s\texact %d\tcooptimal %d\tnil_both %d\tscore %d\tsuboptimal %d",
      f, c.exact or 0, c.cooptimal or 0, c.nil_both or 0, c.score or 0, c.suboptimal or 0)
  end
  out("\n## disagreements\t%d", #bad)
  for i = 1, math.min(#bad, 20) do
    local b = bad[i]
    out("%s\t%s\t%s/%s/%s\t%s", b.case.id, b.result.status, b.case.family, b.case.mode,
      b.case.scheme, tostring(b.result.detail))
  end
  -- how often Biopython had a unique optimal alignment, and whether we hit it
  local uniq, uniq_exact = 0, 0
  for _, c in ipairs(doc.cases) do
    if c.n_optimal == 1 then
      uniq = uniq + 1
      if O.check(doc, c).status == "exact" then uniq_exact = uniq_exact + 1 end
    end
  end
  out("\nunique_optimum_cases\t%d\tengine_matched\t%d", uniq, uniq_exact)
  -- strand selection on the reverse-complement family
  local rcn, rcok = 0, 0
  for _, c in ipairs(doc.cases) do
    if c.rc_score then
      rcn = rcn + 1
      local o = O.opts(doc, c)
      local aln = A.align_strands(c.ref, c.query, o)
      local best = math.max(c.score or -math.huge, c.rc_score or -math.huge)
      if aln and aln.score == best and aln.strand == (c.rc_score > (c.score or -math.huge) and -1 or 1) then
        rcok = rcok + 1
      end
    end
  end
  out("strand_cases\t%d\tcorrect\t%d", rcn, rcok)
end

---------------------------------------------------------------------------
-- The banded path against the full matrix
---------------------------------------------------------------------------

if want_section("band") then
  out("\n## band vs full matrix (reads and plasmids, the sizes that really band)")
  local b = R.band_check(root, 12)
  out("pairs\t%d\tsame_score\t%d\tband_worse\t%d\tdeclined\t%d\tworst_gap\t%d",
    b.n, b.same, b.worse, b.declined, b.worst_gap)
  for _, c in ipairs(b.cases) do out("band_case\t%s", c) end

  out("\n## band vs full matrix (corpus cases, max_cells forced to 1)")
  local n, same, worse, nilled = 0, 0, 0, 0
  local by_ident = {}
  for _, c in ipairs(doc.cases) do
    if not c.protein and c.mode ~= "global" and #c.ref >= 40 and #c.query >= 40 then
      local o = O.opts(doc, c)
      local fullaln = A.align(c.ref, c.query, o)
      local o2 = vim.deepcopy(o)
      o2.max_cells = 1
      local band = A.align(c.ref, c.query, o2)
      n = n + 1
      local bucket = fullaln and (fullaln.identity >= 0.85 and "high" or "low") or "none"
      by_ident[bucket] = by_ident[bucket] or { n = 0, same = 0 }
      by_ident[bucket].n = by_ident[bucket].n + 1
      if not band then
        nilled = nilled + 1
      elseif fullaln and band.score == fullaln.score then
        same = same + 1
        by_ident[bucket].same = by_ident[bucket].same + 1
      else
        worse = worse + 1
      end
    end
  end
  out("banded_cases\t%d\tsame_score\t%d\tdifferent\t%d\tdeclined\t%d", n, same, worse, nilled)
  for _, k in ipairs({ "high", "low", "none" }) do
    local b = by_ident[k]
    if b then out("identity_%s\t%d/%d", k, b.same, b.n) end
  end
end

---------------------------------------------------------------------------
-- Simulated reads
---------------------------------------------------------------------------

if want_section("reads") then
  out("\n## simulated reads")
  local rows = R.run(root)
  out("plasmid\terror_rate\treads\tplaced\tstrand\twrap\tpos\tdiffs_exact\tmed_err\ttp\tfp\tfn\trecall")
  for _, r in ipairs(rows.per_rate) do
    out("%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f", r.plasmid, r.rate, r.n,
      r.placed, r.strand_ok, r.wrap_ok, r.pos_ok, r.diffs_ok, r.median_err,
      r.tp, r.fp, r.fn, r.tp / math.max(1, r.tp + r.fn))
  end
  -- pooled over the three plasmids, by error rate
  out("")
  local pool = {}
  for _, r in ipairs(rows.per_rate) do
    local p = pool[r.rate] or { n = 0, placed = 0, pos_ok = 0, wrap_ok = 0, diffs_ok = 0,
      tp = 0, fp = 0, fn = 0 }
    for _, k in ipairs({ "n", "placed", "pos_ok", "wrap_ok", "diffs_ok", "tp", "fp", "fn" }) do
      p[k] = p[k] + r[k]
    end
    pool[r.rate] = p
  end
  local rates = {}
  for k in pairs(pool) do rates[#rates + 1] = k end
  table.sort(rates)
  out("pooled_rate\treads\tplaced\tpos_ok\twrap_ok\tdiffs_exact\ttp\tfp\tfn\trecall\tprecision")
  for _, k in ipairs(rates) do
    local p = pool[k]
    out("%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f", k, p.n, p.placed, p.pos_ok,
      p.wrap_ok, p.diffs_ok, p.tp, p.fp, p.fn,
      p.tp / math.max(1, p.tp + p.fn), p.tp / math.max(1, p.tp + p.fp))
  end
  out("totals\t-\t%d\t%d\t%d\t%d\t%d\t%d", rows.total.n, rows.total.placed,
    rows.total.strand_ok, rows.total.wrap_ok, rows.total.pos_ok, rows.total.diffs_ok)
  out("origin_spanning_reads\t%d\tdetected\t%d", rows.total.wrapping, rows.total.wrap_detected)
  for _, f in ipairs(rows.failures) do out("fail\t%s", f) end
end

---------------------------------------------------------------------------
-- Timings and memory
---------------------------------------------------------------------------

if want_section("timing") then
  out("\n## timings")
  local function load_gb(stem)
    return assert(pio.read(root .. "/tests/fixtures/ncbi/" .. stem .. ".gb")).records[1]
      .seq:tostring():upper()
  end
  local pbr = load_gb("J01749")          -- 4361 bp, circular
  local yp = load_gb("NC_005816")        -- 9609 bp, circular
  local hbb = load_gb("U01317")          -- 73308 bp
  local prots = {}
  for _, line in ipairs(vim.split(R.read(root .. "/tests/fixtures/protein/proteins.faa"), "\n")) do
    if line:sub(1, 1) == ">" then prots[#prots + 1] = {}
    elseif line ~= "" then prots[#prots][#prots[#prots] + 1] = line end
  end
  local p1, p2 = table.concat(prots[3]), table.concat(prots[4])

  local function timed(label, ref, query, o)
    collectgarbage("collect")
    local before = collectgarbage("count")
    local t0 = vim.uv.hrtime()
    local aln, why = A.align(ref, query, o)
    local ms = (vim.uv.hrtime() - t0) / 1e6
    local peak = collectgarbage("count") - before
    collectgarbage("collect")
    out("%s\tref %d\tquery %d\tcells %.3g\t%s\t%.1f ms\tlua_kb %.0f\t%s", label, #ref, #query,
      #ref * #query, aln and (aln.banded and "band" or "full") or "nil", ms, peak,
      aln and ("score " .. aln.score) or tostring(why))
    return ms, aln
  end

  math.randomseed(11)
  for _, len in ipairs({ 200, 400, 800, 1200, 2000 }) do
    local at = math.random(1, #pbr - len)
    timed("read_vs_pBR322", pbr, R.mutate(pbr:sub(at, at + len - 1), 0.02, 0.004, 0.004),
      { mode = "semi" })
  end
  timed("pBR322_vs_pBR322", pbr, R.mutate(pbr, 0.01, 0, 0), { mode = "global" })
  timed("NC_005816_vs_itself", yp, R.mutate(yp, 0.01, 0.001, 0.001), { mode = "global" })
  timed("protein_vs_protein", p1, p2, { mode = "global", protein = true, gap_open = -10, gap_extend = -1 })
  timed("U01317_vs_itself_band", hbb, R.mutate(hbb, 0.005, 0, 0), { mode = "global" })
  -- The same pair through the full matrix takes minutes and ~16 GB of Lua
  -- heap; it is run only when asked for explicitly.
  if arg[1] == "limit" then
    timed("U01317_vs_itself_full", hbb, R.mutate(hbb, 0.005, 0, 0),
      { mode = "global", max_cells = math.huge })
  end
  timed("U01317_read_1kb", hbb, R.mutate(hbb:sub(30000, 31000), 0.02, 0, 0), { mode = "semi" })

  out("\n## full-matrix memory")
  -- three traceback strings of #ref * #query bytes, plus six number rows
  for _, pair in ipairs({ { 1000, 1000 }, { 2000, 2000 }, { 4361, 4361 }, { 9609, 9609 } }) do
    local n, m = pair[1], pair[2]
    out("full_matrix\t%d x %d\ttraceback_MB\t%.1f\trows_MB\t%.2f", n, m,
      3 * n * m / 1048576, 6 * m * 8 / 1048576)
  end
  out("default_max_cells\t%g", A.defaults.max_cells)
end
