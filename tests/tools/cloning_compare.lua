-- Development driver: run every case of the cloning oracle through the engine
-- and write a report of the disagreements.
--
--   nvim -l tests/tools/cloning_compare.lua [out.json] [filter]
--
-- The committed spec (tests/spec/cloning_oracle_spec.lua) asserts the same
-- comparison; this tool exists to look at the failures in bulk while the
-- engine is being fixed.
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = table.concat({
  root .. "/lua/?.lua", root .. "/lua/?/init.lua", root .. "/tests/?.lua", package.path,
}, ";")
vim.opt.runtimepath:prepend(root)
vim.notify = function() end
require("plasmid").setup({})

local D = require("cloning_driver")
local out_path = arg[1] or "/tmp/cloning_report.json"
local filter = arg[2]

local fh = assert(io.open(root .. "/tests/oracle/expected/cloning.json", "rb"))
local doc = vim.json.decode(fh:read("*a"), { luanil = { object = true, array = true } })
fh:close()

local report = { counts = doc.counts, bad = {}, ok = 0, errors = 0, timing = {}, asm_counts = {} }
local t0 = vim.uv.hrtime()
for _, case in ipairs(doc.cases) do
  if not filter or case.id:find(filter, 1, true) or (case.method or ""):find(filter, 1, true) then
    local t1 = vim.uv.hrtime()
    local ok, bad, res = pcall(D.compare, case, root)
    local ms = (vim.uv.hrtime() - t1) / 1e6
    if not ok then
      report.errors = report.errors + 1
      report.bad[#report.bad + 1] = { id = case.id, label = case.label, method = case.method,
        designed = case.designed, oracle = case.oracle, error = tostring(bad) }
    elseif #bad > 0 then
      report.bad[#report.bad + 1] = { id = case.id, label = case.label, method = case.method,
        designed = case.designed, oracle = case.oracle, complaints = bad,
        want_circular = #case.expect.circular, got_circular = #res.circular,
        problems = res.problems, frag_lens = res.frag_lens, oracle_lens = case.frag_lens }
    else
      report.ok = report.ok + 1
      if case.n_circ_raw then
        report.asm_counts[#report.asm_counts + 1] = { id = case.id, method = case.method,
          oracle = case.n_circ_raw, engine = res.ncircular_assemblies,
          products = #res.circular }
      end
    end
    report.timing[#report.timing + 1] = { id = case.id, ms = ms }
  end
end
report.total_ms = (vim.uv.hrtime() - t0) / 1e6
local w = assert(io.open(out_path, "w"))
w:write(vim.json.encode(report))
w:close()
io.stdout:write(("%d ok, %d disagree, %d errored (%.0f ms) -> %s\n")
  :format(report.ok, #report.bad - report.errors, report.errors, report.total_ms, out_path))
