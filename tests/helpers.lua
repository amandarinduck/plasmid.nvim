local H = {}

--- Report a line from a spec. Writes to stdout directly: print() in `nvim -l`
--- goes through the message area, where lines accumulate, and a redraw after
--- they outgrow 'columns' crashes Neovim 0.12.2.
function H.log(...)
  local parts = {}
  for i = 1, select("#", ...) do parts[i] = tostring((select(i, ...))) end
  io.stdout:write(table.concat(parts, "\t"), "\n")
end

function H.eq(expected, actual, msg)
  if not vim.deep_equal(expected, actual) then
    error(("%sexpected: %s\n     got: %s"):format(
      msg and (msg .. "\n") or "", vim.inspect(expected), vim.inspect(actual)), 2)
  end
end

function H.ok(cond, msg)
  if not cond then error(msg or "assertion failed", 2) end
end

function H.read(path)
  local fh = assert(io.open(path, "rb"))
  local s = fh:read("*a")
  fh:close()
  return s
end

function H.fixture(rel)
  return H.root .. "/tests/fixtures/" .. rel
end

function H.fixtures(pattern)
  local out = vim.fn.globpath(H.root .. "/tests/fixtures", pattern, false, true)
  table.sort(out)
  return out
end

function H.oracle(path)
  local rel = path:sub(#(H.root .. "/tests/fixtures/") + 1):gsub("/", "__")
  local p = H.root .. "/tests/oracle/expected/" .. rel .. ".json"
  if vim.fn.filereadable(p) == 0 then return nil end
  return vim.json.decode(H.read(p), { luanil = { object = true, array = true } })
end

--- First line that differs between two texts (for readable diffs).
function H.first_diff(a, b)
  local la = vim.split(a, "\n", { plain = true })
  local lb = vim.split(b, "\n", { plain = true })
  for i = 1, math.max(#la, #lb) do
    if la[i] ~= lb[i] then
      return ("line %d\n  want: %s\n  got:  %s"):format(i, tostring(la[i]), tostring(lb[i]))
    end
  end
  return nil
end

--- Independent count of feature key lines (5-space indent + key) in every
--- GenBank fixture, used to check that no feature is silently dropped.
function H.raw_feature_count()
  local n = 0
  for _, path in ipairs(H.fixtures("**/*.gb")) do
    for line in H.read(path):gmatch("[^\n]+") do
      if line:match("^     [%a_%-'][^%s]*") then n = n + 1 end
    end
  end
  return n
end

return H
