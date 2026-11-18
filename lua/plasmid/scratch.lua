-- :DnaEdit — edit a region of a record as plain text.
--
-- The region opens in an ordinary buffer (all Vim motions, native undo). On :w
-- the text is validated, compared with the original residues, and the
-- differences are applied to the record as one change made of minimal
-- insertions/deletions/replacements, so features inside the region keep
-- following their bases. The scratch buffer remembers its source buffer,
-- range and the record's sequence revision; if the sequence changed in the
-- meantime the write is refused rather than applied to stale coordinates.
local M = {}

local function buffer() return require("plasmid.buffer") end
local function notify(msg, level) require("plasmid.util").notify("plasmid: " .. msg, level or vim.log.levels.INFO) end

---------------------------------------------------------------------------
-- Diff (pure)
---------------------------------------------------------------------------

-- Myers O(ND) edit script between a and b; nil if more than maxd edits.
-- Returns a list of ops "=", "-" (delete a char), "+" (insert b char) in order.
local function myers(a, b, maxd)
  local n, m = #a, #b
  local v = { [1] = 0 }
  local trace = {}
  for d = 0, math.min(n + m, maxd) do
    local snap = {}
    for k, x in pairs(v) do snap[k] = x end
    trace[d] = snap
    for k = -d, d, 2 do
      local x
      if k == -d or (k ~= d and (v[k - 1] or -1) < (v[k + 1] or -1)) then
        x = v[k + 1] or 0
      else
        x = (v[k - 1] or 0) + 1
      end
      local y = x - k
      while x < n and y < m and a:byte(x + 1) == b:byte(y + 1) do x, y = x + 1, y + 1 end
      v[k] = x
      if x >= n and y >= m then
        -- backtrack
        local script = {}
        local cx, cy = n, m
        for dd = d, 0, -1 do
          local tv = trace[dd]
          local kk = cx - cy
          local pk
          if kk == -dd or (kk ~= dd and (tv[kk - 1] or -1) < (tv[kk + 1] or -1)) then pk = kk + 1 else pk = kk - 1 end
          local px = tv[pk] or 0
          local py = px - pk
          while cx > px and cy > py do
            script[#script + 1] = "="
            cx, cy = cx - 1, cy - 1
          end
          if dd > 0 then
            if cx == px then script[#script + 1] = "+"; cy = cy - 1 else script[#script + 1] = "-"; cx = cx - 1 end
          end
        end
        -- reverse
        local out = {}
        for i = #script, 1, -1 do out[#out + 1] = script[i] end
        return out
      end
    end
  end
  return nil
end

--- Hunks turning `old` into `new`: { { at = offset0 in old, del = n, ins = "..." }, ... }
--- in increasing order. Common prefix/suffix are trimmed first; if the middle
--- needs more than `maxd` edits it becomes one replacement.
function M.diff(old, new, maxd)
  maxd = maxd or 1000
  if old == new then return {} end
  local p = 0
  local lim = math.min(#old, #new)
  while p < lim and old:byte(p + 1) == new:byte(p + 1) do p = p + 1 end
  local q = 0
  while q < lim - p and old:byte(#old - q) == new:byte(#new - q) do q = q + 1 end
  local a, b = old:sub(p + 1, #old - q), new:sub(p + 1, #new - q)
  local script = myers(a, b, maxd)
  if not script then return { { at = p, del = #a, ins = b } } end
  local hunks, i, j = {}, 0, 0 -- positions consumed in a and b
  local cur
  for _, op in ipairs(script) do
    if op == "=" then
      if cur then hunks[#hunks + 1] = cur; cur = nil end
      i, j = i + 1, j + 1
    else
      cur = cur or { at = p + i, del = 0, ins = {} }
      if op == "-" then cur.del = cur.del + 1; i = i + 1
      else cur.ins[#cur.ins + 1] = b:sub(j + 1, j + 1); j = j + 1 end
    end
  end
  if cur then hunks[#hunks + 1] = cur end
  for _, h in ipairs(hunks) do h.ins = type(h.ins) == "table" and table.concat(h.ins) or h.ins end
  return hunks
end

--- Apply hunks to a string (for tests and sanity checks).
function M.patch(old, hunks)
  local out, pos = {}, 0
  for _, h in ipairs(hunks) do
    out[#out + 1] = old:sub(pos + 1, h.at)
    out[#out + 1] = h.ins
    pos = h.at + h.del
  end
  out[#out + 1] = old:sub(pos + 1)
  return table.concat(out)
end

---------------------------------------------------------------------------
-- Scratch buffer
---------------------------------------------------------------------------

local WIDTH = 60

local function chunk(s)
  local lines = {}
  for i = 1, #s, WIDTH do lines[#lines + 1] = s:sub(i, i + WIDTH - 1) end
  if #lines == 0 then lines[1] = "" end
  return lines
end

local function winbar(sbuf)
  local p = vim.b[sbuf].plasmid_edit
  if not p then return "" end
  local st = buffer().get(p.src)
  local name = st and buffer().record(st).name or "?"
  return ("%%#PlasmidTitle# editing %s %d..%d%%*  (%d %s)  │  :w applies to the record, :q discards"):format(
    (name:gsub("%%", "%%%%")), p.start, p.stop, p.stop - p.start + 1, p.unit)
end
M.winbar = winbar

--- Open s..e of the record in `buf` for editing (s <= e).
function M.open(buf, s, e)
  local st = buffer().get(buf)
  if not st then return notify("not a plasmid buffer", vim.log.levels.WARN) end
  local rec = buffer().record(st)
  if s > e then
    return notify("the region crosses the origin; rotate the origin (:DnaRotate) or edit it in two parts",
      vim.log.levels.WARN)
  end
  local raw = rec.seq:sub(s, e)
  local shown = require("plasmid.select").display_case(rec, raw)
  vim.cmd("botright new")
  local sbuf = vim.api.nvim_get_current_buf()
  local name = ("plasmid://edit/%d/%s"):format(buf, rec.name)
  pcall(vim.api.nvim_buf_set_name, sbuf, name)
  vim.bo[sbuf].buftype = "acwrite"
  vim.bo[sbuf].bufhidden = "wipe"
  vim.bo[sbuf].swapfile = false
  vim.bo[sbuf].filetype = "plasmidedit"
  vim.api.nvim_buf_set_lines(sbuf, 0, -1, false, chunk(shown))
  vim.bo[sbuf].modified = false
  vim.b[sbuf].plasmid_edit = {
    src = buf, start = s, stop = e, revision = rec.seq_revision or 0, original = raw, unit = rec:unit(),
  }
  local win = vim.api.nvim_get_current_win()
  vim.wo[win].winbar = ("%%{%%v:lua.require'plasmid.scratch'.winbar(%d)%%}"):format(sbuf)
  vim.wo[win].wrap = false
  vim.wo[win].spell = false
  vim.api.nvim_create_autocmd("BufWriteCmd", {
    buffer = sbuf,
    callback = function() M.apply(sbuf, { force = vim.v.cmdbang == 1 }) end, -- :w! forces past cut_policy=refuse
  })
  return sbuf
end

--- Apply the scratch buffer to its record (the BufWriteCmd handler).
function M.apply(sbuf, opts)
  opts = opts or {}
  local p = vim.b[sbuf].plasmid_edit
  if not p then return notify("not a :DnaEdit buffer", vim.log.levels.ERROR) end
  local st = buffer().get(p.src)
  if not st or not vim.api.nvim_buf_is_valid(p.src) then
    return notify("the sequence buffer this region came from is gone", vim.log.levels.ERROR)
  end
  local rec = buffer().record(st)
  if (rec.seq_revision or 0) ~= p.revision then
    return notify(("the sequence changed since :DnaEdit opened %d..%d; not applied (reopen :DnaEdit)"):format(
      p.start, p.stop), vim.log.levels.ERROR)
  end
  local text = table.concat(vim.api.nvim_buf_get_lines(sbuf, 0, -1, false), "")
  local se = require("plasmid.seqedit")
  local new
  if text:gsub("[%s%d/]", "") == "" then
    new = ""
  else
    local err
    new, err = se.clean(rec, text)
    if not new then return notify(err, vim.log.levels.ERROR) end
  end
  if new == p.original then
    vim.bo[sbuf].modified = false
    return notify("no changes")
  end
  if #new == 0 and p.stop - p.start + 1 >= rec:len() then
    return notify("refusing to delete the whole sequence", vim.log.levels.ERROR)
  end
  local hunks = M.diff(p.original, new)
  local ops = {}
  for i = #hunks, 1, -1 do -- right to left: earlier positions stay valid
    local h = hunks[i]
    local s = p.start + h.at
    ops[#ops + 1] = { op = "splice", start = s, stop = s + h.del - 1, ins = h.ins }
  end
  local op = #ops == 1 and ops[1] or { op = "batch", ops = ops }
  local report = se.run(p.src, op, {
    label = (":DnaEdit %d..%d (%d change%s)"):format(p.start, p.stop, #hunks, #hunks == 1 and "" or "s"),
    force = opts.force, cursor = p.start, select = #new > 0 and { p.start, p.start + #new - 1 } or false,
  })
  if not report then return end -- refused by the cut policy; the buffer stays modified
  p.stop = p.start + #new - 1
  p.original = new
  p.revision = rec.seq_revision or 0
  vim.b[sbuf].plasmid_edit = p
  vim.bo[sbuf].modified = false
  return report
end

return M
