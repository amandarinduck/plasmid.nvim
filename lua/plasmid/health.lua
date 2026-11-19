local M = {}

function M.check()
  local h = vim.health
  h.start("plasmid.nvim")
  if vim.fn.has("nvim-0.10") == 1 then
    h.ok("Neovim " .. tostring(vim.version()))
  else
    h.error("Neovim >= 0.10 required")
  end
  local ok_bit = pcall(require, "bit")
  if ok_bit then h.ok("LuaJIT bit library available (braille canvas)") else h.error("`bit` module missing") end
  if vim.o.termguicolors then
    h.ok("'termguicolors' is on (feature colours exact)")
  else
    h.warn("'termguicolors' is off; feature colours fall back to the 256-colour cube")
  end
  if vim.o.ambiwidth == "single" then
    h.ok("'ambiwidth' is single (box-drawing and braille glyphs align)")
  else
    h.error("'ambiwidth' is '" .. vim.o.ambiwidth .. "': tracks and the circular map will misalign; set it to single")
  end
  local cfgc = require("plasmid.config").options.yank_clipboard
  if cfgc == false then
    h.info("yank_clipboard = false: yanks stay in Neovim registers")
  elseif vim.fn.has("clipboard") == 1 then
    local name = vim.g.clipboard and (type(vim.g.clipboard) == "table" and (vim.g.clipboard.name or "g:clipboard") or tostring(vim.g.clipboard))
      or (vim.fn.exists("*provider#clipboard#Executable") == 1 and vim.fn["provider#clipboard#Executable"]()) or "?"
    h.ok(("clipboard provider: %s (yanks reach the system clipboard)"):format(name ~= "" and name or "?"))
  else
    h.warn("no clipboard provider: yanks are sent to the terminal as OSC 52, which works in most terminals "
      .. "(kitty, WezTerm, iTerm2, Alacritty, foot, Windows Terminal; tmux needs `set -g set-clipboard on`)",
      { "Install xclip/xsel (X11) or wl-clipboard (Wayland), or see :help clipboard" })
  end
  local cfg = require("plasmid.config").options
  local mcfg = (require("plasmid.config").options or {}).mouse or {}
  if mcfg.enabled == false then
    h.info("mouse integration off (config mouse.enabled = false)")
  elseif vim.o.mouse == "" then
    h.warn("'mouse' is empty: drag-select, double-click and Ctrl-click primer pairs need :set mouse=a")
  else
    h.ok(("mouse '%s': drag selects, double click selects a feature, Ctrl-click pairs primers")
      :format(vim.o.mouse))
  end
  if vim.fn.executable("curl") == 1 then
    h.ok("curl found (NCBI fetch: :DnaFetch, :DnaSearch, :DnaGene)")
  else
    h.warn("curl not found: NCBI fetch is unavailable; everything else works")
  end
  local ncfg = (require("plasmid.config").options or {}).ncbi or {}
  if ncfg.api_key or (vim.env.NCBI_API_KEY or "") ~= "" then
    h.info("NCBI API key set: up to 10 requests/s")
  else
    h.info("no NCBI API key: requests are spaced to NCBI's 3/s limit")
  end
  h.info(("default view: %s, circular.aspect = %.2f (check with :DnaCalibrate)"):format(
    cfg.default_view, cfg.circular.aspect))
  h.info("Glyphs used: ━ ╌ ▶ ◀ ┃ and braille U+2800–U+28FF; any font with box-drawing + braille works")
  local n = vim.tbl_count(require("plasmid.buffer").state)
  h.info(("%d sequence buffer%s open"):format(n, n == 1 and "" or "s"))
end

return M
