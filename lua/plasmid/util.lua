-- Small shared helpers.
local M = {}

--- vim.notify that keeps info/warning messages within v:echospace, so they
--- never trigger a hit-enter prompt (and avoid a Neovim 0.12 headless crash
--- when redrawing after a message wider than the screen). Errors are shown
--- in full.
function M.notify(msg, level)
  level = level or vim.log.levels.INFO
  msg = tostring(msg)
  if level < vim.log.levels.ERROR then
    local room = math.max(20, (vim.v.echospace or vim.o.columns) - 1)
    if vim.fn.strdisplaywidth(msg) > room then
      msg = vim.fn.strcharpart(msg, 0, room - 1) .. "…"
    end
  end
  vim.notify(msg, level)
end

local repaint_ns = vim.api.nvim_create_namespace("plasmid_repaint")

--- Make Neovim redraw every line of `buf` (so the decoration provider runs
--- again), in whichever windows show it. Setting and deleting a highlighted
--- extmark over the buffer does this through the public API; nvim__redraw
--- {buf = ...} crashes Neovim 0.12.2 ("double free") when another window is
--- current, and :redraw! can crash it headless after a wide message.
--- Report buffers (digests, gels, assembly reports, primer and part lists)
--- are plain text, so `/` there is Vim's own search, which is case-sensitive
--- unless 'ignorecase' is set. Sequence buffers search motifs case-blind, so
--- make these match: `/` and `?` start the pattern with \c.
function M.case_insensitive_search(buf)
  if vim.o.ignorecase then return end
  for _, lhs in ipairs({ "/", "?" }) do
    vim.keymap.set("n", lhs, lhs .. "\\c", { buffer = buf, silent = false, nowait = true,
      desc = "plasmid: search, ignoring case" })
  end
end

function M.repaint(buf)
  if not vim.api.nvim_buf_is_valid(buf) then return end
  local n = vim.api.nvim_buf_line_count(buf)
  local ok, id = pcall(vim.api.nvim_buf_set_extmark, buf, repaint_ns, 0, 0,
    { end_row = math.max(n - 1, 0), end_col = 0, hl_group = "PlasmidRepaint" })
  if ok then pcall(vim.api.nvim_buf_del_extmark, buf, repaint_ns, id) end
end

return M
