local H = require("helpers")
local tm = require("plasmid.bio.tm")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local TOL = 0.01

local function fixture()
  return vim.json.decode(H.read(H.root .. "/tests/oracle/expected/tm.json"), { luanil = { object = true, array = true } })
end

--- NEB Ta rounding: to 0.1 °C, then to an integer, ties up
local function neb_round(x)
  return tm._round(tm._round(x * 10) / 10)
end

local function near(a, b, msg)
  H.ok(type(a) == "number" and math.abs(a - b) < TOL, ("%s: got %s, want %.4f"):format(msg, tostring(a), b))
end

local function label(c)
  return ("[%s %s] %s%s %s"):format(c.source, c.kind, c.seq, c.seq2 and ("/" .. c.seq2) or "",
    vim.inspect(c.opts or { conc_nM = c.conc_nM }, { newline = "", indent = "" }))
end

--- Check every oracle case from one source; returns the number checked.
local function check_source(source)
  local n = 0
  for _, c in ipairs(fixture().cases) do
    if c.source == source then
      n = n + 1
      local e, msg = c.expect, label(c)
      if c.kind == "raw" then
        local _, _, t = tm.raw_nn(c.seq, c.conc_nM)
        near(t, e.raw, msg)
      elseif c.kind == "tm" then
        local d, err = tm.tm_detail(c.seq, c.opts)
        H.ok(d, msg .. ": " .. tostring(err))
        if e.tm then near(d.tm, e.tm, msg) end
        if e.tm_rounded then H.eq(e.tm_rounded, d.tm_rounded, msg) end
        if e.tm_rounded then H.eq(e.tm_rounded, tm.tm(c.seq, c.opts), msg) end
      elseif c.kind == "ta" then
        local ta, t1, t2, info = tm.ta(c.seq, c.seq2, c.opts)
        H.ok(ta, msg .. ": " .. tostring(t1))
        if e.tm1 then H.eq(e.tm1, t1, msg .. " tm1") end
        if e.tm2 then H.eq(e.tm2, t2, msg .. " tm2") end
        if e.ta then H.eq(e.ta, ta, msg .. " ta") end
        if e.ta_decimal then near(info.ta_decimal, e.ta_decimal, msg .. " ta_decimal") end
        if e.tm1_unrounded then near(tm.tm_detail(c.seq, c.opts).tm, e.tm1_unrounded, msg .. " tm1 unrounded") end
        if e.tm2_unrounded then near(tm.tm_detail(c.seq2, c.opts).tm, e.tm2_unrounded, msg .. " tm2 unrounded") end
      end
    end
  end
  return n
end

case("oracle_polymerase_tm", function()
  local n = check_source("polymerase-tm")
  H.ok(n >= 150, "too few polymerase-tm cases: " .. n)
end)

case("oracle_biopython_unrounded", function()
  local n = check_source("biopython")
  H.ok(n >= 100, "too few biopython cases: " .. n)
end)

case("oracle_neb_api", function()
  local n = check_source("neb")
  H.ok(n >= 100, "too few NEB cases: " .. n)
end)

case("oracle_neb_api_long_format_unrounded", function()
  local n = check_source("neb-long")
  H.ok(n >= 20, "too few NEB long-format cases: " .. n)
end)

case("oracle_neb_procedure", function()
  local n = check_source("neb-procedure")
  H.ok(n >= 22, "too few NEB-procedure cases: " .. n)
end)

case("oracle_rounding_ties", function()
  local n = check_source("rounding")
  H.ok(n >= 1, "no rounding cases")
end)

case("fixture_covers_every_polymerase_and_buffer", function()
  local seen_p, seen_b = {}, {}
  for _, c in ipairs(fixture().cases) do
    if c.source == "polymerase-tm" then
      seen_p[c.opts.polymerase] = true
      if c.opts.buffer then seen_b[c.opts.buffer] = true end
    end
  end
  for _, p in ipairs(tm.POLYMERASES) do H.ok(seen_p[p.key], "no polymerase-tm case for " .. p.key) end
  for b in pairs(tm.BUFFERS) do H.ok(seen_b[b], "no polymerase-tm case for buffer " .. b) end
end)

-- NEB API documentation example (tmapi.neb.com/docs/reference): Q5, 0.5 uM
case("neb_docs_example", function()
  local ta, t1, t2 = tm.ta("AGCGGATAACAATTTCACACAGGA", "GTAAAACGACGGCCAGT", { polymerase = "q5" })
  H.eq({ 63, 66, 62 }, { ta, t1, t2 })
end)

case("polymerase_table", function()
  for i = 2, #tm.POLYMERASES do
    H.ok(tm.POLYMERASES[i - 1].key < tm.POLYMERASES[i].key, "POLYMERASES not sorted by key")
  end
  for _, p in ipairs(tm.POLYMERASES) do
    H.ok(tm.BUFFERS[p.buffer], p.key .. ": unknown buffer " .. p.buffer)
    H.ok(p.ta_rule == "offset" or p.ta_rule == "phusion", p.key)
    H.ok(p.method == "owczarzy" or p.method == "schildkraut", p.key)
    H.ok(type(p.neb_code) == "string", p.key .. ": no neb_code")
    H.ok(type(p.name) == "string" and p.conc_nM > 0 and p.ta_cap > 0 and p.min_len >= 0, p.key)
  end
  H.eq("q5_hot_start", tm.polymerase("Q5 Hot-Start").key)
  H.eq("phusion_hf", tm.polymerase("PHUSION-HF").key)
  H.eq(nil, tm.polymerase("pfu"))
  H.eq(nil, tm.polymerase(nil))
end)

case("valid_and_invalid_letters", function()
  H.eq(true, (tm.valid("ACGTacgt")))
  for _, bad in ipairs({ "", "ACGN", "ACGU", "ACG T", "ACGR", "ACGT-" }) do
    local ok, err = tm.valid(bad)
    H.ok(ok == false and type(err) == "string", "accepted " .. vim.inspect(bad))
    local v, e2 = tm.tm(bad)
    H.ok(v == nil and type(e2) == "string", "tm accepted " .. vim.inspect(bad))
    local ta, e3 = tm.ta("ACGTACGTACGTACGTAC", bad)
    H.ok(ta == nil and e3:find("primer 2", 1, true), "ta accepted " .. vim.inspect(bad))
  end
  H.eq(false, (tm.valid(nil)))
  H.eq(false, (tm.valid(42)))
end)

case("lowercase_equals_uppercase", function()
  local s = "AGCGGATAACAATTTCACACAGGA"
  H.eq(tm.tm_detail(s), tm.tm_detail(s:lower()))
  H.eq({ tm.raw_nn(s, 500) }, { tm.raw_nn(s:lower(), 500) })
  H.eq({ tm.ta(s, "gtaaaacgacggccagt") }, { tm.ta(s:lower(), "GTAAAACGACGGCCAGT") })
end)

case("raw_nn_hand_computed", function()
  -- AAAAAAAA: 7 x AA (-7.9, -22.2) + 2 x terminal A (2.3, 4.1)
  local dH, dS, t = tm.raw_nn("AAAAAAAA", 500)
  near(dH, -50.7, "dH")
  near(dS, -147.2, "dS")
  near(t, -50700 / (-147.2 + 1.987 * math.log(500e-9)) - 273.15, "Tm")
end)

case("ta_min_len_and_short_primers", function()
  -- Q5: min_len 8, so two 8-mers still get +1
  local a, b = "GCGCGCGC", "GGCCGGCC"
  local low = math.min(tm.tm_detail(a).tm, tm.tm_detail(b).tm)
  local ta, t1, t2, info = tm.ta(a, b, { polymerase = "q5" })
  H.eq(neb_round(low + 1), ta)
  H.eq("offset", info.rule)
  -- polymerase-tm compat: rule on the integer Tms
  local tc, c1, c2 = tm.ta(a, b, { polymerase = "q5", compat = "polymerase-tm" })
  H.eq(math.min(c1, c2) + 1, tc)
  H.eq({ t1, t2 }, { c1, c2 }) -- Q5 salt is the same in both tables
  -- Vent: min_len 21; a 20-mer gets no -2 offset, a 21-mer does
  local p20 = "ATGACCATGATTACGCCAAG"
  local p21 = "ATGACCATGATTACGCCAAGC"
  local v = { polymerase = "vent" }
  local ta20, _, _, i20 = tm.ta(p20, p21, v)
  H.eq(neb_round(math.min(tm.tm_detail(p20, v).tm, tm.tm_detail(p21, v).tm)), ta20)
  H.eq("short", i20.rule)
  local ta21, _, _, i21 = tm.ta(p21, p21, v)
  H.eq(neb_round(tm.tm_detail(p21, v).tm - 2), ta21)
  H.eq("offset", i21.rule)
  -- Phusion ignores length; rule on the unrounded lower Tm
  local ph = { polymerase = "phusion_hf" }
  local tp, _, _, ip = tm.ta(a, b, ph)
  H.eq(neb_round(0.93 * math.min(tm.tm_detail(a, ph).tm, tm.tm_detail(b, ph).tm) + 7.5), tp)
  H.eq("phusion", ip.rule)
end)

case("ta_caps", function()
  local g1 = "GCGGCGGCCGCGGCGGCCGCGGCGGCGGCGG"
  local g2 = "CCGCCGCCGCGGCCGCCGCGGCCGCCGCCGC"
  for key, cap in pairs({ q5 = 72, taq = 68, longamp = 65, phusion_hf = 72, vent = 72 }) do
    local ta, _, _, info = tm.ta(g1, g2, { polymerase = key })
    H.eq(cap, ta, key)
    H.eq(true, info.capped, key)
    H.ok(info.note:find("capped", 1, true), key)
  end
  local _, _, _, info = tm.ta("ATGACCATGATTACGCCAAGC", "GTAAAACGACGGCCAGT", { polymerase = "q5" })
  H.eq(false, info.capped)
end)

case("dmso", function()
  local s = "ATGACCATGATTACGCCAAGC"
  local d0, d5 = tm.tm_detail(s), tm.tm_detail(s, { dmso = 5 })
  near(d5.tm, d0.tm - 3.0, "tm with 5% DMSO")
  H.eq(d0.raw, d5.raw)
  -- Ta: DMSO after the cap; tm1/tm2 are the DMSO-free Tms
  local g = "GCGGCGGCCGCGGCGGCCGCGGCGGCGGCGG"
  local ta, t1, t2, info = tm.ta(g, g, { polymerase = "q5", dmso = 5 })
  H.eq(69, ta)
  H.eq(true, info.capped)
  H.eq(tm.tm(g), t1)
  H.eq(t1, t2)
  -- exact .5 tie (2.5 % DMSO -> 1.5 °C) rounds half up (NEB/JavaScript), not to even
  local ta2, a1, a2 = tm.ta(g, g, { polymerase = "q5", dmso = 2.5 })
  H.eq(71, ta2) -- 72 - 1.5 = 70.5 -> 71 (Python round() would give 70)
  H.eq(a1, a2)
  H.eq(71, tm._round(70.5))
  H.eq(-2, tm._round(-2.5))
  H.eq(-3, tm._round(-2.51))
  H.eq(nil, (tm.tm(s, { dmso = -1 })))
end)

case("custom_salt_conc_buffer_method", function()
  local s = "AGCGGATAACAATTTCACACAGGA"
  local d = tm.tm_detail(s, { salt_mM = 50, conc_nM = 250 })
  H.eq(50, d.salt_mM)
  H.eq(250, d.conc_nM)
  H.eq("owczarzy", d.method)
  local _, _, raw250 = tm.raw_nn(s, 250)
  near(d.raw, raw250, "raw at 250 nM")
  -- salt override beats buffer, buffer beats polymerase default
  H.eq(40, tm.tm_detail(s, { buffer = "thermopol" }).salt_mM)
  H.eq(40, tm.tm_detail(s, { buffer = "ThermoPol" }).salt_mM)
  H.eq(77, tm.tm_detail(s, { buffer = "thermopol", salt_mM = 77 }).salt_mM)
  H.eq(nil, (tm.tm_detail(s, { buffer = "nope" })))
  -- Phusion: Schildkraut with conc/4
  local p = tm.tm_detail(s, { polymerase = "phusion_hf" })
  H.eq("schildkraut", p.method)
  H.eq(125, p.conc_eff_nM)
  near(p.tm, p.raw + 16.6 * math.log10(0.222), "Schildkraut")
  local p1 = tm.tm_detail(s, { polymerase = "phusion_hf", conc_nM = 1000 })
  H.eq(250, p1.conc_eff_nM)
  -- method override
  H.eq("owczarzy", tm.tm_detail(s, { polymerase = "phusion_hf", method = "owczarzy" }).method)
  H.eq(nil, (tm.tm_detail(s, { method = "wallace" })))
  -- owczarzy2008: Mg-free equals 2004; default 2 mM Mg
  near(tm.tm_detail(s, { method = "owczarzy2008", mg_mM = 0 }).tm, tm.tm_detail(s).tm, "2008 Mg=0")
  H.eq(2, tm.tm_detail(s, { method = "owczarzy2008" }).mg_mM)
  for _, bad in ipairs({ 0, -5, 0 / 0, "50" }) do
    H.eq(nil, (tm.tm_detail(s, { salt_mM = bad })), "salt " .. tostring(bad))
    H.eq(nil, (tm.tm_detail(s, { conc_nM = bad })), "conc " .. tostring(bad))
  end
  H.eq(nil, (tm.tm(s, { polymerase = "pfu" })))
  H.eq(nil, (tm.ta(s, s, { polymerase = "pfu" })))
end)

case("neb_vs_polymerase_tm_salts", function()
  local s = "ATGACCATGATTACGCCAAGC"
  for key, pair in pairs({ taq = { 55, 51 }, onetaq = { 54, 51 }, onetaq_gc = { 80, 51 }, q5u_hot_start = { 170, 180 },
    crimson_taq = { 55, 51 }, q5 = { 150, 150 }, phusion_hf = { 222, 222 } }) do
    H.eq(pair[1], tm.tm_detail(s, { polymerase = key }).salt_mM, key)
    H.eq(pair[2], tm.tm_detail(s, { polymerase = key, compat = "polymerase-tm" }).salt_mM, key)
  end
  H.eq("polymerase-tm", tm.tm_detail(s, { compat = "polymerase-tm" }).compat)
  H.eq("neb", tm.tm_detail(s).compat)
  H.eq(nil, (tm.tm(s, { compat = "primer3" })))
  H.eq(nil, (tm.ta(s, s, { compat = "primer3" })))
end)

case("neb_two_step_ta_rounding", function()
  -- NEB batch API (v0.19.0): Q5, tm1 = 63, tm2 = 75, ta = 65 (not 64): the lower
  -- unrounded Tm 63.465 -> +1 = 64.465 -> 64.5 -> 65
  local s1, s2
  for _, c in ipairs(fixture().cases) do
    if c.source == "neb" and c.neb_prodcode == "q5-0" and c.expect.ta and c.expect.ta == math.min(c.expect.tm1, c.expect.tm2) + 2 then
      s1, s2 = c.seq, c.seq2
      break
    end
  end
  H.ok(s1, "fixture has no NEB Q5 case with Ta = lower Tm + 2")
  local ta, t1, t2, info = tm.ta(s1, s2, { polymerase = "q5" })
  H.eq(math.min(t1, t2) + 2, ta)
  local low = math.min(tm.tm_detail(s1).tm, tm.tm_detail(s2).tm)
  H.ok(low % 1 >= 0.45 and low % 1 < 0.5, "lower Tm fraction " .. low % 1)
  near(info.ta_decimal, tm._round((low + 1) * 10) / 10, "ta_decimal")
  local tc = tm.ta(s1, s2, { polymerase = "q5", compat = "polymerase-tm" })
  H.eq(math.min(t1, t2) + 1, tc)
end)

return T
