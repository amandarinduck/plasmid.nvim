-- Dev tool: print a view to stdout without a UI.
--   nvim -l tests/tools/render.lua FILE linear|circular WIDTH [HEIGHT] [MAXLINES]
package.path = "lua/?.lua;lua/?/init.lua;" .. package.path
local config = require("plasmid.config"); config.setup({})
require("plasmid.view.highlights").setup()
local res = assert(require("plasmid.io").read(arg[1]))
local rec = res.records[1]
local lay
if arg[2] == "circular" then
  lay = require("plasmid.view.circular").build(rec, config.options, tonumber(arg[3]), tonumber(arg[4]))
else
  lay = require("plasmid.view.linear").build(rec, config.options, tonumber(arg[3]))
end
local n = tonumber(arg[5]) or #lay.lines
for i = 1, math.min(n, #lay.lines) do io.stdout:write(lay.lines[i], "\n") end
