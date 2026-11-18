local H = require("helpers")
local rebase = require("plasmid.rebase")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local function fx(name) return H.fixture("rebase/" .. name) end
local function by_name(list)
  local m = {}
  for _, e in ipairs(list) do m[e.name] = e end
  return m
end

-- every case that writes uses a private storage path
local function with_store(fn)
  local saved = rebase.path
  local dir = vim.fn.tempname()
  rebase.path = dir .. "/plasmid/rebase_enzymes.json"
  rebase.reset()
  local ok, err = pcall(fn, rebase.path)
  rebase.clear()
  vim.fn.delete(dir, "rf")
  rebase.path = saved
  rebase.reset()
  if not ok then error(err, 0) end
end

case("emboss_e_converts_cuts_like_biopython", function()
  local list, skipped = rebase.parse_emboss_e(H.read(fx("emboss_e.test")))
  local m = by_name(list)
  -- name = { site, cut5, cut3 } in the convention of plasmid.data.enzymes
  local want = {
    EcoRI = { "GAATTC", 1, 5 },              -- G^AATTC
    KpnI = { "GGTACC", 5, 1 },               -- GGTAC^C, 3' overhang
    BsaI = { "GGTCTC", 7, 11 },              -- GGTCTC(1/5), downstream Type IIS
    TspRI = { "CASTG", 7, -2 },              -- bottom strand cut upstream (EMBOSS -3)
    MboI = { "GATC", 0, 4 },                 -- ^GATC: EMBOSS -1 -> before first base
    AvaI = { "CYCGRG", 1, 5 },               -- degenerate
    BglI = { "GCCNNNNNGGC", 7, 4 },          -- interrupted palindrome
    EcoRV = { "GATATC", 3, 3 },              -- lower-case input is upper-cased
    ["CviKI-1"] = { "RGCY", 2, 2 },          -- name kept verbatim
    DupI = { "GACCGA", 17, 15 },             -- first of two records kept
  }
  H.eq(vim.tbl_count(want), #list)
  for name, w in pairs(want) do
    local e = m[name]
    H.ok(e, "missing " .. name)
    H.eq(w, { e.site, e.cut5, e.cut3 }, name)
    local keys = vim.tbl_keys(e)
    table.sort(keys)
    H.eq({ "cut3", "cut5", "name", "site" }, keys, name .. " fields")
  end
  H.eq({ two_cuts = 1, no_cut = 1, bad_site = 1, duplicate = 1, malformed = 1 }, skipped)
  H.ok(not m.TwoCutI and not m.NoCutI and not m.BadSiteI, "skipped records leaked")
end)

case("emboss_e_matches_bundled_table", function()
  local bundled = {}
  for _, e in ipairs(require("plasmid.data.enzymes")) do bundled[e.name] = e end
  local n = 0
  for _, e in ipairs((rebase.parse_emboss_e(H.read(fx("emboss_e.test"))))) do
    local b = bundled[e.name]
    if b then
      H.eq({ b.site, b.cut5, b.cut3 }, { e.site, e.cut5, e.cut3 }, e.name)
      n = n + 1
    end
  end
  H.eq(8, n) -- EcoRI KpnI BsaI TspRI MboI AvaI BglI EcoRV
end)

case("emboss_e_handles_crlf_and_empty", function()
  local text = H.read(fx("emboss_e.test")):gsub("\n", "\r\n")
  local list = rebase.parse_emboss_e(text)
  H.eq(10, #list)
  H.eq("GAATTC", by_name(list).EcoRI.site)
  local none, sk = rebase.parse_emboss_e("")
  H.eq({}, none)
  H.eq(0, sk.malformed)
end)

case("emboss_r_suppliers", function()
  local sup = rebase.parse_emboss_r(H.read(fx("emboss_r.test")))
  H.eq("NRF", sup.EcoRI)
  H.eq("NF", sup.KpnI)
  H.eq("N", sup.BsaI)       -- record with a reference line
  H.eq("N", sup["CviKI-1"])
  H.eq("", sup.DupI)        -- not commercially available
  H.eq(nil, sup.EcoRV)      -- not in the file
  H.eq(9, vim.tbl_count(sup))
end)

case("emboss_s_supplier_names", function()
  local s = rebase.parse_emboss_s(H.read(fx("emboss_s.test")))
  H.eq({ F = "Thermo Fisher Scientific", N = "New England Biolabs", R = "Promega Corporation" }, s)
end)

case("import_writes_sorted_list_and_summary", function()
  with_store(function(path)
    local list, summary = rebase.import(fx("emboss_e.test"), fx("emboss_r.test"), fx("emboss_s.test"))
    H.ok(list, tostring(summary))
    H.eq(1, vim.fn.filereadable(path))
    local names = {}
    for i, e in ipairs(list) do
      names[i] = e.name
      if i > 1 then H.ok(list[i - 1].name < e.name, "sorted at " .. e.name) end
    end
    H.eq({ "AvaI", "BglI", "BsaI", "CviKI-1", "DupI", "EcoRI", "EcoRV", "KpnI", "MboI", "TspRI" }, names)
    H.eq(10, summary.count)
    H.eq({ two_cuts = 1, no_cut = 1, bad_site = 1, duplicate = 1, malformed = 1 }, summary.skipped)
    H.eq(8, summary.commercial)
    local m = by_name(list)
    H.eq("NRF", m.EcoRI.suppliers)
    H.eq("", m.EcoRV.suppliers)
    H.eq("", m.DupI.suppliers)
    H.eq("New England Biolabs", rebase.supplier_names().N)
  end)
end)

case("import_without_r_file_has_no_suppliers", function()
  with_store(function()
    local list, summary = rebase.import(fx("emboss_e.test"))
    H.ok(list, tostring(summary))
    H.eq(0, summary.commercial)
    for _, e in ipairs(list) do H.eq(nil, e.suppliers, e.name) end
    H.eq(nil, rebase.supplier_names())
  end)
end)

case("load_cache_reset_clear", function()
  with_store(function(path)
    H.eq(nil, rebase.load())
    local list = assert(rebase.import(fx("emboss_e.test"), fx("emboss_r.test")))
    H.ok(rebase.load() == list, "import primes the cache")
    rebase.reset()
    local again = rebase.load()
    H.ok(again ~= list, "reset drops the cache")
    H.eq(#list, #again)
    for i, e in ipairs(list) do
      H.eq({ e.name, e.site, e.cut5, e.cut3, e.suppliers },
        { again[i].name, again[i].site, again[i].cut5, again[i].cut3, again[i].suppliers })
    end
    H.ok(rebase.load() == again, "cached")
    H.eq(true, rebase.clear())
    H.eq(0, vim.fn.filereadable(path))
    H.eq(nil, rebase.load())
    H.eq(false, rebase.clear())
  end)
end)

case("import_errors", function()
  with_store(function(path)
    local list, err = rebase.import(fx("does_not_exist"))
    H.eq(nil, list); H.ok(err:match("cannot read"), err)
    list, err = rebase.import(fx("emboss_s.test")) -- no e-format records
    H.eq(nil, list); H.ok(err:match("no usable enzymes"), err)
    list, err = rebase.import(fx("emboss_e.test"), fx("missing_r"))
    H.eq(nil, list); H.ok(err:match("cannot read"), err)
    H.eq(0, vim.fn.filereadable(path), "nothing written on error")
  end)
end)

case("restriction_uses_imported_set", function()
  local R = require("plasmid.bio.restriction")
  with_store(function()
    assert(rebase.import(fx("emboss_e.test"), fx("emboss_r.test")))
    R.reset()
    local list = R.enzymes()
    H.eq(10, #list)
    local e = R.get("cviki-1")
    H.ok(e and e.suppliers == "N", "imported enzyme visible via restriction")
    H.eq(-2, R.get("TspRI").cut3)
  end)
  R.reset()
  H.ok(#R.enzymes() >= 200, "bundled table back after clear")
end)

return T
