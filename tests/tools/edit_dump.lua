-- Apply random sequence edits to a GenBank file with the plugin, write the
-- result, and dump what the plugin believes: residues and each feature's
-- location and bases. Used by tests/tools/edit_oracle.py.
--   nvim -l tests/tools/edit_dump.lua in.gb out.gb out.json seed n_edits
local root = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":p:h:h:h")
package.path = root .. "/lua/?.lua;" .. root .. "/lua/?/init.lua;" .. package.path
local pio = require("plasmid.io")
local tx = require("plasmid.model.transaction")
local inp, outp, jsonp, seed, n = arg[1], arg[2], arg[3], tonumber(arg[4]), tonumber(arg[5])
math.randomseed(seed)
local file = assert(pio.read(inp))
local rec = file.records[1]
local function bases(k)
  local t = {}
  for i = 1, k do local r = math.random(4); t[i] = ("acgt"):sub(r, r) end
  return table.concat(t)
end
for _ = 1, n do
  local L, r = rec:len(), math.random()
  local s = math.random(1, L)
  if r < 0.3 then tx.commit(rec, { op = "splice", start = s, stop = s - 1, ins = bases(math.random(1, 9)) })
  elseif r < 0.55 then tx.commit(rec, { op = "splice", start = s, stop = math.min(L, s + math.random(0, 20)), ins = "" })
  elseif r < 0.8 then tx.commit(rec, { op = "splice", start = s, stop = math.min(L, s + math.random(0, 10)), ins = bases(math.random(1, 12)) })
  elseif r < 0.92 then tx.commit(rec, { op = "revcomp", start = s, stop = math.min(L, s + math.random(5, 300)) })
  elseif rec:is_circular() then tx.commit(rec, { op = "rotate", origin = s }) end
end
assert(pio.write(outp, file.records, "genbank", file.file, file.format))
local feats = {}
for _, f in ipairs(rec.features) do
  feats[#feats + 1] = { type = f.type, location = f.location_text,
    residues = (f.loc and not require("plasmid.model.location").has_remote(f.parts)) and rec:extract(f) or vim.NIL }
end
local fh = assert(io.open(jsonp, "w"))
fh:write(vim.json.encode({ seq = rec.seq:tostring(), features = feats }))
fh:close()
