-- Dev tool: dump the plugin's reading of GenBank/GenPept files as JSON, in the
-- same shape as tests/oracle/dump_ncbi_real.py, plus parse time and Lua heap.
--   nvim -l tests/tools/ncbi_real_dump.lua OUT.json FILE...
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = root .. "/lua/?.lua;" .. root .. "/lua/?/init.lua;" .. package.path
local pio = require("plasmid.io")

local function quals(f)
  local out = {}
  for _, q in ipairs(f.qualifiers) do out[#out + 1] = { q.key, q.value == true and "" or q.value } end
  table.sort(out, function(a, b)
    if a[1] ~= b[1] then return a[1] < b[1] end
    return a[2] < b[2]
  end)
  return out
end

local res_all = {}
for i = 2, #arg do
  local path = arg[i]
  local name = vim.fn.fnamemodify(path, ":t")
  local fh = assert(io.open(path, "rb")); local text = fh:read("*a"); fh:close()
  collectgarbage("collect")
  local kb0 = collectgarbage("count")
  local t0 = vim.uv.hrtime()
  local ok, res, err = pcall(pio.parse, text, path)
  local parse_ms = (vim.uv.hrtime() - t0) / 1e6
  local entry = { parse_ms = parse_ms, bytes = #text }
  if not ok or not res then
    entry.error = tostring(ok and err or res)
  else
    collectgarbage("collect")
    entry.heap_kb = collectgarbage("count") - kb0
    local t1 = vim.uv.hrtime()
    local out = pio.format(res.records, res.format, res.file)
    entry.format_ms = (vim.uv.hrtime() - t1) / 1e6
    entry.roundtrip_identical = out == text
    if out ~= text then
      local la, lb = vim.split(text, "\n", { plain = true }), vim.split(out, "\n", { plain = true })
      for k = 1, math.max(#la, #lb) do
        if la[k] ~= lb[k] then entry.roundtrip_first_diff = { k, la[k] or vim.NIL, lb[k] or vim.NIL }; break end
      end
    end
    entry.records = {}
    for _, r in ipairs(res.records) do
      local seq = r.seq:tostring()
      local feats = {}
      for _, f in ipairs(r.features) do
        local parts = {}
        for _, p in ipairs(f.parts) do parts[#parts + 1] = { p.start, p.stop, p.strand, p.ref or vim.NIL } end
        feats[#feats + 1] = { type = f.type, location = f.location_text, loc_error = f.loc_error or vim.NIL,
          parts = parts, qualifiers = quals(f) }
      end
      entry.records[#entry.records + 1] = {
        name = r.name, id = r.id, description = r.description, length = r:len(),
        declared_length = r.declared_length or vim.NIL, alphabet = r.alphabet,
        sha256 = vim.fn.sha256(seq:upper()), topology = r.topology,
        molecule_type = r.molecule_type or vim.NIL, division = r.division or vim.NIL, features = feats,
      }
    end
  end
  res_all[name] = entry
  io.stdout:write(name, "\t", entry.error or #entry.records, "\t", ("%.1f ms"):format(parse_ms), "\n")
end
local fh = assert(io.open(arg[1], "wb"))
fh:write(vim.json.encode(res_all))
fh:close()
