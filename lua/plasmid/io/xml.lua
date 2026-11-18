-- Minimal XML parser for the small, well-formed documents embedded in
-- SnapGene files (features, primers, notes). Pure Lua.
--
-- parse(s) -> root node { tag, attrs = {..}, children = {node|string, ...} }
-- Handles: declarations, comments, CDATA, self-closing tags, attributes in
-- either quote style, predefined and numeric character entities.
local M = {}

local ENT = { amp = "&", lt = "<", gt = ">", quot = '"', apos = "'" }

local function utf8_char(cp)
  if cp < 0x80 then return string.char(cp) end
  if cp < 0x800 then
    return string.char(0xC0 + math.floor(cp / 64), 0x80 + cp % 64)
  end
  if cp < 0x10000 then
    return string.char(0xE0 + math.floor(cp / 4096), 0x80 + math.floor(cp / 64) % 64, 0x80 + cp % 64)
  end
  return string.char(0xF0 + math.floor(cp / 262144), 0x80 + math.floor(cp / 4096) % 64,
    0x80 + math.floor(cp / 64) % 64, 0x80 + cp % 64)
end

M.utf8_char = utf8_char

function M.escape(s)
  return (tostring(s):gsub("[&<>\"\n\r]", {
    ["&"] = "&amp;", ["<"] = "&lt;", [">"] = "&gt;", ['"'] = "&quot;", ["\n"] = "&#10;", ["\r"] = "&#13;",
  }))
end

function M.unescape(s)
  return (s:gsub("&(#?x?)(%w+);", function(kind, body)
    if kind == "#" then return utf8_char(tonumber(body, 10) or 63) end
    if kind == "#x" then return utf8_char(tonumber(body, 16) or 63) end
    return ENT[body] or ("&" .. body .. ";")
  end))
end

local function parse_attrs(s)
  local attrs = {}
  for k, q, v in s:gmatch("([%w_:%-%.]+)%s*=%s*(['\"])(.-)%2") do
    attrs[k] = M.unescape(v)
    local _ = q
  end
  return attrs
end

function M.parse(s)
  local root = { tag = "#document", attrs = {}, children = {} }
  local stack = { root }
  local i, n = 1, #s
  while i <= n do
    local lt = s:find("<", i, true)
    local top = stack[#stack]
    if not lt then
      local text = s:sub(i)
      if text:match("%S") then top.children[#top.children + 1] = M.unescape(text) end
      break
    end
    if lt > i then
      local text = s:sub(i, lt - 1)
      if text:match("%S") or #stack > 1 then top.children[#top.children + 1] = M.unescape(text) end
    end
    if s:sub(lt, lt + 3) == "<!--" then
      local e = assert(s:find("-->", lt + 4, true), "unterminated comment")
      i = e + 3
    elseif s:sub(lt, lt + 8) == "<![CDATA[" then
      local e = assert(s:find("]]>", lt + 9, true), "unterminated CDATA")
      top.children[#top.children + 1] = s:sub(lt + 9, e - 1)
      i = e + 3
    elseif s:sub(lt, lt + 1) == "<?" or s:sub(lt, lt + 1) == "<!" then
      local e = assert(s:find(">", lt, true), "unterminated declaration")
      i = e + 1
    elseif s:sub(lt, lt + 1) == "</" then
      local e = assert(s:find(">", lt, true), "unterminated end tag")
      local tag = s:sub(lt + 2, e - 1):match("^%s*(.-)%s*$")
      if top.tag ~= tag then
        error(("mismatched </%s>, expected </%s>"):format(tag, top.tag), 0)
      end
      top.e = e -- byte span s..e covers the element including its end tag
      stack[#stack] = nil
      i = e + 1
    else
      -- start tag: find the closing '>' outside quoted attribute values
      local j, quote = lt + 1, nil
      while j <= n do
        local c = s:sub(j, j)
        if quote then
          if c == quote then quote = nil end
        elseif c == '"' or c == "'" then
          quote = c
        elseif c == ">" then
          break
        end
        j = j + 1
      end
      assert(j <= n, "unterminated start tag")
      local body = s:sub(lt + 1, j - 1)
      local self_closing = body:sub(-1) == "/"
      if self_closing then body = body:sub(1, -2) end
      local tag, rest = body:match("^([^%s/>]+)(.*)$")
      local node = { tag = tag, attrs = parse_attrs(rest or ""), children = {}, s = lt }
      if self_closing then node.e = j end
      top.children[#top.children + 1] = node
      if not self_closing then stack[#stack + 1] = node end
      i = j + 1
    end
  end
  if #stack ~= 1 then error("unclosed <" .. stack[#stack].tag .. ">", 0) end
  return root
end

--- All descendant elements named `tag` (document order).
function M.find_all(node, tag, out)
  out = out or {}
  for _, c in ipairs(node.children) do
    if type(c) == "table" then
      if c.tag == tag then out[#out + 1] = c end
      M.find_all(c, tag, out)
    end
  end
  return out
end

--- First descendant named `tag`.
function M.find(node, tag)
  return M.find_all(node, tag)[1]
end

--- Concatenated direct text content of an element.
function M.text(node)
  if not node then return nil end
  local out = {}
  for _, c in ipairs(node.children) do
    if type(c) == "string" then out[#out + 1] = c end
  end
  return table.concat(out)
end

return M
