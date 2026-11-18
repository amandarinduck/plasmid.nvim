-- GenBank/EMBL/DDBJ feature location grammar. Pure Lua.
--
-- The AST mirrors the source text so that format(parse(s)) == s for every
-- supported construct (lossless), while flatten() gives the biological meaning.
--
-- Node kinds
--   range    { start, stop, start_mark?, stop_mark?, ref? }   100..200, <1..>50
--   site     { pos, mark?, ref? }                             467
--   between  { left, right, ref? }                            23^24
--   within   { start, stop, ref? }                            102.110
--   complement { child }
--   join / order { parts = { node, ... } }
--
-- Coordinates are 1-based inclusive throughout.
local M = {}

---------------------------------------------------------------------------
-- Tokenizer
---------------------------------------------------------------------------

local function tokenize(s)
  local toks, i, n = {}, 1, #s
  while i <= n do
    local c = s:sub(i, i)
    if c:match("%s") then
      i = i + 1
    elseif c:match("%d") then
      local j = i
      while j <= n and s:sub(j, j):match("%d") do j = j + 1 end
      toks[#toks + 1] = { t = "num", v = tonumber(s:sub(i, j - 1)) }
      i = j
    elseif c:match("[%a_]") then
      local j = i
      while j <= n and s:sub(j, j):match("[%w_%.%-]") do j = j + 1 end
      -- An identifier never ends with '.', so "J00194.1:" lexes cleanly and
      -- trailing dots are handed back.
      while s:sub(j - 1, j - 1) == "." do j = j - 1 end
      toks[#toks + 1] = { t = "id", v = s:sub(i, j - 1) }
      i = j
    elseif c == "." then
      if s:sub(i + 1, i + 1) == "." then
        toks[#toks + 1] = { t = ".." }
        i = i + 2
      else
        toks[#toks + 1] = { t = "." }
        i = i + 1
      end
    elseif c:match("[%(%),:<>%^]") then
      toks[#toks + 1] = { t = c }
      i = i + 1
    else
      error(("unexpected character %q at %d"):format(c, i), 0)
    end
  end
  toks[#toks + 1] = { t = "eof" }
  return toks
end

---------------------------------------------------------------------------
-- Recursive-descent parser
---------------------------------------------------------------------------

-- bond(a,b) is not INSDC nucleotide syntax but NCBI's GenPept form for
-- disulfide bonds and heterogen contacts (UniProt- and PDB-derived proteins);
-- it is read like order() of single sites and written back as bond().
local FUNCS = { join = true, order = true, complement = true, bond = true }

local function parser(toks)
  local p = { i = 1 }

  local function peek(k) return toks[p.i + (k or 0)] end
  local function next_()
    local t = toks[p.i]
    p.i = p.i + 1
    return t
  end
  local function expect(tt)
    local t = next_()
    if t.t ~= tt then
      error(("expected %s, got %s"):format(tt, t.v or t.t), 0)
    end
    return t
  end

  local function position()
    local mark
    local t = peek()
    if t.t == "<" or t.t == ">" then
      mark = t.t
      next_()
    end
    return expect("num").v, mark
  end

  local parse_loc

  local function simple(ref)
    local a, amark = position()
    local t = peek()
    if t.t == ".." then
      next_()
      local b, bmark = position()
      return { kind = "range", start = a, stop = b, start_mark = amark, stop_mark = bmark, ref = ref }
    elseif t.t == "^" then
      next_()
      local b = expect("num").v
      if amark then error("fuzzy marker not allowed in between-location", 0) end
      return { kind = "between", left = a, right = b, ref = ref }
    elseif t.t == "." then
      next_()
      local b = expect("num").v
      if amark then error("fuzzy marker not allowed in within-location", 0) end
      return { kind = "within", start = a, stop = b, ref = ref }
    end
    return { kind = "site", pos = a, mark = amark, ref = ref }
  end

  parse_loc = function()
    local t = peek()
    if t.t == "id" then
      if peek(1).t == "(" and FUNCS[t.v] then
        next_()
        next_()
        local node
        if t.v == "complement" then
          node = { kind = "complement", child = parse_loc() }
        else
          local parts = { parse_loc() }
          while peek().t == "," do
            next_()
            parts[#parts + 1] = parse_loc()
          end
          node = { kind = t.v, parts = parts }
        end
        expect(")")
        return node
      elseif peek(1).t == ":" then
        next_()
        next_()
        return simple(t.v)
      end
      error(("unsupported location operator %q"):format(t.v), 0)
    end
    return simple(nil)
  end

  return function()
    local node = parse_loc()
    expect("eof")
    return node
  end
end

--- Parse a location string. Returns node or nil, err.
function M.parse(s)
  local ok, res = pcall(function()
    return parser(tokenize(s))()
  end)
  if ok then return res end
  return nil, res
end

---------------------------------------------------------------------------
-- Formatting (inverse of parse)
---------------------------------------------------------------------------

function M.format(node)
  local k = node.kind
  local ref = node.ref and (node.ref .. ":") or ""
  if k == "range" then
    return ref .. (node.start_mark or "") .. node.start .. ".." .. (node.stop_mark or "") .. node.stop
  elseif k == "site" then
    return ref .. (node.mark or "") .. node.pos
  elseif k == "between" then
    return ref .. node.left .. "^" .. node.right
  elseif k == "within" then
    return ref .. node.start .. "." .. node.stop
  elseif k == "complement" then
    return "complement(" .. M.format(node.child) .. ")"
  else
    local out = {}
    for i, p in ipairs(node.parts) do out[i] = M.format(p) end
    return k .. "(" .. table.concat(out, ",") .. ")"
  end
end

---------------------------------------------------------------------------
-- Semantics
---------------------------------------------------------------------------

--- Flatten into parts in biological (transcription) order.
-- Each part: { start, stop, strand = 1|-1, ref?, start_fuzzy?, stop_fuzzy?, kind }
-- `start <= stop` always holds numerically; a "between" part is empty (stop = start - 1).
function M.flatten(node, strand, out)
  strand = strand or 1
  out = out or {}
  local k = node.kind
  if k == "range" then
    out[#out + 1] = {
      start = node.start, stop = node.stop, strand = strand, ref = node.ref, kind = k,
      start_fuzzy = node.start_mark, stop_fuzzy = node.stop_mark,
    }
  elseif k == "site" then
    out[#out + 1] = { start = node.pos, stop = node.pos, strand = strand, ref = node.ref, kind = k,
      start_fuzzy = node.mark, stop_fuzzy = node.mark }
  elseif k == "within" then
    out[#out + 1] = { start = node.start, stop = node.stop, strand = strand, ref = node.ref, kind = k }
  elseif k == "between" then
    out[#out + 1] = { start = node.left + 1, stop = node.left, strand = strand, ref = node.ref, kind = k }
  elseif k == "complement" then
    local sub = M.flatten(node.child, -strand, {})
    for i = #sub, 1, -1 do out[#out + 1] = sub[i] end
  else
    for _, p in ipairs(node.parts) do M.flatten(p, strand, out) end
  end
  return out
end

--- Overall strand: 1, -1, or 0 (mixed).
function M.strand(parts)
  local s
  for _, p in ipairs(parts) do
    if s == nil then s = p.strand elseif s ~= p.strand then return 0 end
  end
  return s or 1
end

--- Total residues covered (sum of part lengths, local parts only).
function M.length(parts)
  local n = 0
  for _, p in ipairs(parts) do
    if not p.ref then n = n + (p.stop - p.start + 1) end
  end
  return n
end

--- True if any part lives on another record (remote accession).
function M.has_remote(parts)
  for _, p in ipairs(parts) do if p.ref then return true end end
  return false
end

--- Covered extents on the forward-strand axis, for layout.
-- A feature that crosses the origin of a circular record yields two extents
-- ({first.start, seq_len} and {1, last.stop}); otherwise one {min, max}.
-- Origin crossing is read from biological order: on the + strand, a part that
-- starts before its predecessor ends means the feature wrapped.
---@return { [1]: integer, [2]: integer }[]
function M.extents(parts, seq_len, circular)
  local local_parts = {}
  for _, p in ipairs(parts) do
    if not p.ref and p.stop >= p.start then local_parts[#local_parts + 1] = p end
  end
  if #local_parts == 0 then
    -- empty (between-only) features still need a position to draw a tick
    for _, p in ipairs(parts) do
      if not p.ref then
        local pos = math.max(1, math.min(p.stop, seq_len))
        return { { pos, pos } }
      end
    end
    return {}
  end
  local lo, hi = math.huge, -math.huge
  for _, p in ipairs(local_parts) do
    lo = math.min(lo, p.start)
    hi = math.max(hi, p.stop)
  end
  local strand = M.strand(local_parts)
  if not circular or strand == 0 or #local_parts < 2 then return { { lo, hi } } end
  -- left-to-right order on the forward axis
  local ordered = {}
  if strand == 1 then
    ordered = local_parts
  else
    for i = #local_parts, 1, -1 do ordered[#ordered + 1] = local_parts[i] end
  end
  for i = 2, #ordered do
    if ordered[i].start < ordered[i - 1].stop then
      -- wrap between i-1 and i
      local first_start = math.huge
      for j = 1, i - 1 do first_start = math.min(first_start, ordered[j].start) end
      local last_stop = -math.huge
      for j = i, #ordered do last_stop = math.max(last_stop, ordered[j].stop) end
      return { { first_start, seq_len }, { 1, last_stop } }
    end
  end
  return { { lo, hi } }
end

--- Extract the feature's residues (biological orientation) from a ResidueStore.
--- Remote parts are skipped.
function M.extract(parts, store)
  local alphabet = require("plasmid.model.alphabet")
  local out = {}
  for _, p in ipairs(parts) do
    if not p.ref and p.stop >= p.start then
      local s = store:sub(p.start, p.stop)
      if p.strand == -1 then s = alphabet.revcomp(s) end
      out[#out + 1] = s
    end
  end
  return table.concat(out)
end

--- Positions in transcription order (for codon mapping). Remote parts skipped.
function M.positions(parts)
  local pos = {}
  for _, p in ipairs(parts) do
    if not p.ref then
      if p.strand == -1 then
        for x = p.stop, p.start, -1 do pos[#pos + 1] = x end
      else
        for x = p.start, p.stop do pos[#pos + 1] = x end
      end
    end
  end
  return pos
end

return M
