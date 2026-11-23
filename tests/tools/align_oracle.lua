-- Replay machinery for the Biopython alignment oracle.
--
-- tests/oracle/expected/alignments.json holds, per case, the two sequences,
-- the mode and scoring scheme, Biopython's optimal score and one optimal
-- alignment (as a cigar plus the starts, in the engine's convention). This
-- module runs the engine on a case and classifies the result:
--
--   exact         same score, same cigar and starts as Biopython's alignment
--   cooptimal     same score, a different path that rescores to the optimum
--   suboptimal    the engine's path does not rescore to its own reported score
--   score         the scores differ
--   nil_both      both declined (empty query, no positive local segment)
--   nil_engine    Biopython aligned, the engine declined
--   nil_oracle    the engine aligned, Biopython declined
--
-- Used by tests/spec/align_oracle_spec.lua and tests/tools/align_report.lua.
local O = {}

local A = require("plasmid.bio.align")

--- Read the committed corpus. `root` is the repository root.
function O.load(root)
  local path = root .. "/tests/oracle/expected/alignments.json"
  local fh = assert(io.open(path, "rb"), "missing " .. path)
  local text = fh:read("*a")
  fh:close()
  return vim.json.decode(text, { luanil = { object = true, array = true } })
end

--- Engine options for a case: the mode's option set plus its scoring scheme.
function O.opts(doc, case)
  local o = vim.deepcopy(doc.modes[case.mode].engine_opts)
  local scheme = (case.protein and doc.aa_schemes or doc.nuc_schemes)[case.scheme]
  for k, v in pairs(scheme) do o[k] = v end
  o.protein = case.protein and true or false
  -- every corpus case is small enough for the full matrix; the banded path is
  -- exercised separately
  o.max_cells = math.huge
  return o
end

--- Score an op list the way the mode does: every op present is charged, and
--- the end gaps the mode makes free are simply absent from the list.
function O.rescore(ref, query, ops, ref_start, query_start, o)
  local score_of = A._scorer(o)
  local total = 0
  local i, j = ref_start, query_start
  for _, op in ipairs(ops) do
    if op.op == "M" then
      for k = 0, op.n - 1 do
        total = total + score_of(ref:byte(i + k), query:byte(j + k))
      end
      i, j = i + op.n, j + op.n
    elseif op.op == "D" then
      total = total + o.gap_open + o.gap_extend * op.n
      i = i + op.n
    else
      total = total + o.gap_open + o.gap_extend * op.n
      j = j + op.n
    end
  end
  return total, i, j
end

local function cigar_of(ops)
  local parts = {}
  for _, op in ipairs(ops) do parts[#parts + 1] = ("%d%s"):format(op.n, op.op) end
  return table.concat(parts)
end

--- Run one case. Returns { status, detail, aln, score, want }.
function O.check(doc, case)
  local o = O.opts(doc, case)
  local aln, why = A.align(case.ref, case.query, o)
  local want = case.score
  if not aln then
    if want == nil or case.alignments == 0 then
      return { status = "nil_both", detail = why }
    end
    return { status = "nil_engine", detail = ("declined (%s), Biopython scored %s")
      :format(tostring(why), tostring(want)), want = want }
  end
  if want == nil or case.alignments == 0 then
    return { status = "nil_oracle", detail = ("engine scored %d, Biopython declined")
      :format(aln.score), aln = aln }
  end
  if aln.score ~= want then
    return { status = "score", detail = ("score %d, Biopython %s"):format(aln.score, tostring(want)),
      aln = aln, want = want, score = aln.score }
  end
  local cig = cigar_of(aln.ops)
  if cig == case.cigar and aln.ref_start == case.ref_start and aln.query_start == case.query_start then
    return { status = "exact", aln = aln, score = aln.score, want = want }
  end
  local re, iend, jend = O.rescore(case.ref, case.query, aln.ops, aln.ref_start, aln.query_start, o)
  if re ~= want then
    return { status = "suboptimal", aln = aln, score = aln.score, want = want,
      detail = ("reported %d, its own path rescores to %d, Biopython's optimum is %s; engine %s@%d/%d vs oracle %s@%d/%d")
        :format(aln.score, re, tostring(want), cig, aln.ref_start, aln.query_start,
          case.cigar, case.ref_start, case.query_start) }
  end
  return { status = "cooptimal", aln = aln, score = aln.score, want = want,
    detail = ("%s@%d/%d vs %s@%d/%d"):format(cig, aln.ref_start, aln.query_start,
      case.cigar, case.ref_start, case.query_start), ref_end = iend, query_end = jend }
end

--- Replay every case (optionally filtered). Returns counts by status and the
--- list of cases that are not exact or cooptimal.
function O.run(doc, filter)
  local counts, bad, by_family = {}, {}, {}
  for _, case in ipairs(doc.cases) do
    if not filter or filter(case) then
      local r = O.check(doc, case)
      counts[r.status] = (counts[r.status] or 0) + 1
      local f = by_family[case.family] or {}
      f[r.status] = (f[r.status] or 0) + 1
      by_family[case.family] = f
      if r.status == "score" or r.status == "suboptimal"
        or r.status == "nil_engine" or r.status == "nil_oracle" then
        bad[#bad + 1] = { case = case, result = r }
      end
    end
  end
  return counts, bad, by_family
end

O.cigar_of = cigar_of

return O
