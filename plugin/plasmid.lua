if vim.g.loaded_plasmid then return end
vim.g.loaded_plasmid = true

if vim.fn.has("nvim-0.10") == 0 then
  vim.notify("plasmid.nvim requires Neovim >= 0.10", vim.log.levels.ERROR)
  return
end

-- Zero-config: sequence files open in plasmid view out of the box.
-- require("plasmid").setup({...}) later reconfigures; or set vim.g.plasmid = {...}
-- before the plugin loads.
require("plasmid")._autoload()
