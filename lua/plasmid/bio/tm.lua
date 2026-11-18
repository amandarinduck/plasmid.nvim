-- Primer melting (Tm) and annealing (Ta) temperatures, NEB Tm Calculator model.
-- Pure Lua, no dependencies.
--
-- Port of polymerase-tm 2.3.2 (https://github.com/Dioskurides/polymerase-tm),
-- src/polymerase_tm/core.py and constants.py, which is distributed under the
-- MIT Licence. Its LICENSE file, including the copyright notice and the
-- permission notice, is reproduced verbatim in docs/polymerase-tm-LICENSE as
-- the licence requires; it is not paraphrased here.
-- polymerase-tm reproduces the NEB Tm Calculator; the NN model, primer
-- concentrations, Ta rules and caps below come from it. Four buffer salts and
-- the Ta rounding procedure were corrected to NEB's own Tm API
-- (tmapi.neb.com, v0.19.0, queried 2026-11-16; see tm_report.md). The exact
-- polymerase-tm behaviour is available with opts.compat = "polymerase-tm".
--
-- Model
--   * Nearest-neighbour duplex enthalpy/entropy, unified parameters and
--     terminal A·T / G·C initiation terms: SantaLucia J. (1998) PNAS 95:1460.
--     No self-complementarity (symmetry) term; total primer concentration C
--     used as-is in Tm = dH / (dS + R ln C) (NEB / polymerase-tm convention).
--   * Monovalent salt correction: Owczarzy R. et al. (2004) Biochemistry
--     43:3537 (eq. 22), 1/Tm = 1/Tm(1 M) + (4.29e-5 fGC - 3.95e-5) ln[Na+]
--     + 9.40e-6 ln^2[Na+].
--   * Phusion family: Schildkraut C. & Lifson S. (1965) Biopolymers 3:195,
--     Tm = Tm(1 M) + 16.6 log10[Na+], primer concentration divided by 4
--     (NEB's Phusion convention; the NEB API echoes conc = 0.125 uM for 0.5).
--   * Optional Mg2+-aware correction ("owczarzy2008"): Owczarzy R. et al.
--     (2008) Biochemistry 47:5336, eq. 16 with the [Na+]-dependent a, d, g of
--     eqs. 18-20 for 0.22 <= sqrt[Mg2+]/[Na+] < 6. This follows the published
--     equation (as Biopython MeltingTemp saltcorr=7 does); polymerase-tm's
--     owczarzy_bivariate() omits the [Na+]-dependent a, d, g in that regime,
--     so the two differ there (see tm_report.md). No free-Mg/dNTP chelation.
--   * DMSO: -0.6 °C per % (v/v), as in polymerase-tm (NEB's API has no DMSO).
--
-- Rounding (NEB, verified against tmapi.neb.com; see tm_report.md)
--   * Tm: the unrounded Tm is rounded once to the nearest integer.
--   * Ta: computed from the UNROUNDED lower Tm, rule applied, capped, then
--     rounded to 0.1 °C (NEB's long-format value) and that is rounded to an
--     integer. Ties go up (JavaScript Math.round): e.g. lower Tm 63.47 with
--     the Q5 rule -> 64.47 -> 64.5 -> 65.
--   polymerase-tm instead applies the Ta rule to the integer Tms and uses
--   Python round() (half to even); compat = "polymerase-tm" reproduces that.
--
-- Ta rules
--   offset : Ta = lower Tm + offset when both primers are >= min_len nt,
--            else Ta = lower Tm; then capped.
--   phusion: Ta = 0.93 * lower Tm + 7.5; then capped.
--   DMSO is subtracted after the cap.

local M = {}

-- SantaLucia 1998 unified NN parameters: dH kcal/mol, dS cal/(mol K)
local NN = {
  AA = { -7.9, -22.2 }, TT = { -7.9, -22.2 },
  AT = { -7.2, -20.4 }, TA = { -7.2, -21.3 },
  CA = { -8.5, -22.7 }, TG = { -8.5, -22.7 },
  GT = { -8.4, -22.4 }, AC = { -8.4, -22.4 },
  CT = { -7.8, -21.0 }, AG = { -7.8, -21.0 },
  GA = { -8.2, -22.2 }, TC = { -8.2, -22.2 },
  CG = { -10.6, -27.2 }, GC = { -9.8, -24.4 },
  GG = { -8.0, -19.9 }, CC = { -8.0, -19.9 },
}
-- terminal initiation
local TERM = { A = { 2.3, 4.1 }, T = { 2.3, 4.1 }, G = { 0.1, -2.8 }, C = { 0.1, -2.8 } }
local R = 1.987 -- cal/(mol K)
local K0 = 273.15

--- NEB buffer effective monovalent salt (mM), as used by NEB's Tm API v0.19.0
--- (fitted exactly from its integer Tms; polymerase-tm values in BUFFERS_PTM).
M.BUFFERS = {
  standard_taq = 55,
  thermopol = 40,
  hemo_klentaq = 70,
  crimson_taq = 55,
  longamp = 100,
  multiplex = 90,
  phusion_hf = 222,
  phusion_gc = 222,
  phusionflex_hf = 222,
  phusionflex_gc = 222,
  onetaq_std = 54,
  onetaq_gc = 80,
  q5 = 150,
  q5mm = 150,
  q5u = 170,
  q5bd = 150,
  longamphs = 100,
}

--- polymerase-tm 2.3.2 BUFFERS (from the NEB Tm Calculator JavaScript).
M.BUFFERS_PTM = {
  standard_taq = 51,
  thermopol = 40,
  hemo_klentaq = 70,
  crimson_taq = 51,
  longamp = 100,
  multiplex = 90,
  phusion_hf = 222,
  phusion_gc = 222,
  phusionflex_hf = 222,
  phusionflex_gc = 222,
  onetaq_std = 51,
  onetaq_gc = 51,
  q5 = 150,
  q5mm = 150,
  q5u = 180,
  q5bd = 150,
  longamphs = 100,
}

local function poly(key, name, buffer, conc, rule, offset, cap, min_len, method, neb)
  return {
    key = key, name = name, buffer = buffer, conc_nM = conc, ta_rule = rule,
    ta_offset = offset, ta_cap = cap, min_len = min_len, method = method, neb_code = neb,
  }
end

--- Polymerase products (sorted by key). neb_code is the tmapi.neb.com prodcode.
--- Ta rules/caps are polymerase-tm's (NEB web calculator) and match NEB's API
--- for every product EXCEPT crimson_taq, multiplex and q5u_hot_start: for
--- ctaq-0, multiplex-0 and q5u-0 the API (v0.19.0) returns Ta = lower Tm with
--- no offset and no cap (observed up to 106 °C), i.e. it has no Ta rule for
--- them. That is not adopted; these three keep polymerase-tm's rules
--- (-5/cap 68, -5/cap 68, +2/cap 72), which are NOT confirmed by the API.
--- Their Tm (buffer salt) does match the API. See tm_report.md, finding 5.
M.POLYMERASES = {
  poly("crimson_taq", "Crimson Taq DNA Polymerase", "crimson_taq", 200, "offset", -5, 68, 8, "owczarzy", "ctaq-0"),
  poly("deep_vent", "Deep Vent DNA Polymerase", "thermopol", 200, "offset", -2, 72, 21, "owczarzy", "deepvent-0"),
  poly("epimark", "EpiMark Hot Start Taq DNA Polymerase", "standard_taq", 200, "offset", -5, 68, 8, "owczarzy", "epimarkhs-0"),
  poly("hemo_klentaq", "Hemo KlenTaq", "hemo_klentaq", 200, "offset", -5, 68, 8, "owczarzy", "hktaq-0"),
  poly("hot_start_taq", "Hot Start Taq DNA Polymerase", "standard_taq", 200, "offset", -5, 68, 8, "owczarzy", "hstaq-0"),
  poly("longamp", "LongAmp Taq DNA Polymerase", "longamp", 400, "offset", -5, 65, 8, "owczarzy", "lataq-0"),
  poly("longamp_hot_start", "LongAmp Hot Start Taq DNA Polymerase", "longamphs", 400, "offset", -5, 65, 8, "owczarzy", "lahstaq-0"),
  poly("multiplex", "Multiplex PCR 5X Master Mix", "multiplex", 200, "offset", -5, 68, 8, "owczarzy", "multiplex-0"),
  poly("onetaq", "OneTaq DNA Polymerase (Standard Buffer)", "onetaq_std", 200, "offset", -5, 68, 8, "owczarzy", "onetaq-0"),
  poly("onetaq_gc", "OneTaq DNA Polymerase (GC Buffer)", "onetaq_gc", 200, "offset", -5, 68, 8, "owczarzy", "onetaq-1"),
  poly("onetaq_hot_start", "OneTaq Hot Start DNA Polymerase", "onetaq_std", 200, "offset", -5, 68, 8, "owczarzy", "onetaqhs-0"),
  poly("phusion_flex_gc", "Phusion Hot Start Flex DNA Polymerase (GC Buffer)", "phusionflex_gc", 500, "phusion", 0, 72, 0, "schildkraut", "phusionflex-1"),
  poly("phusion_flex_hf", "Phusion Hot Start Flex DNA Polymerase (HF Buffer)", "phusionflex_hf", 500, "phusion", 0, 72, 0, "schildkraut", "phusionflex-0"),
  poly("phusion_gc", "Phusion High-Fidelity DNA Polymerase (GC Buffer)", "phusion_gc", 500, "phusion", 0, 72, 0, "schildkraut", "phusion-1"),
  poly("phusion_hf", "Phusion High-Fidelity DNA Polymerase (HF Buffer)", "phusion_hf", 500, "phusion", 0, 72, 0, "schildkraut", "phusion-0"),
  poly("q5", "Q5 High-Fidelity DNA Polymerase", "q5", 500, "offset", 1, 72, 8, "owczarzy", "q5-0"),
  poly("q5_blood_direct", "Q5 Blood Direct 2X Master Mix", "q5bd", 500, "offset", 1, 72, 8, "owczarzy", "q5bd-0"),
  poly("q5_hot_start", "Q5 Hot Start High-Fidelity DNA Polymerase", "q5", 500, "offset", 1, 72, 8, "owczarzy", "q5hs-0"),
  poly("q5_master_mix", "Q5 High-Fidelity 2X Master Mix", "q5mm", 500, "offset", 1, 72, 8, "owczarzy", "q5-1"),
  poly("q5u_hot_start", "Q5U Hot Start High-Fidelity DNA Polymerase", "q5u", 500, "offset", 2, 72, 8, "owczarzy", "q5u-0"),
  poly("taq", "Taq DNA Polymerase", "standard_taq", 200, "offset", -5, 68, 8, "owczarzy", "taq-0"),
  poly("vent", "Vent DNA Polymerase", "thermopol", 200, "offset", -2, 72, 21, "owczarzy", "vent-0"),
}

local BY_KEY = {}
for _, p in ipairs(M.POLYMERASES) do BY_KEY[p.key] = p end

local function norm(key)
  return (tostring(key):lower():gsub("^%s+", ""):gsub("%s+$", ""):gsub("[%s%-]", "_"))
end

--- Polymerase entry by key (case, space and hyphen insensitive), or nil.
function M.polymerase(key)
  if key == nil then return nil end
  return BY_KEY[norm(key)]
end

--- true, or false and a message. A/C/G/T only (any case); empty is invalid.
function M.valid(seq)
  if type(seq) ~= "string" then return false, "sequence must be a string" end
  if #seq == 0 then return false, "empty sequence" end
  local pos = seq:find("[^ACGTacgt]")
  if pos then
    return false, ("invalid base %q at position %d (only A, C, G, T allowed)"):format(seq:sub(pos, pos), pos)
  end
  return true
end

--- Round to nearest, ties towards +inf (JavaScript Math.round; NEB).
local function round_half_up(x)
  local f = math.floor(x)
  if x - f >= 0.5 then return f + 1 end
  return f
end
--- Round to nearest, ties to even (Python round(); polymerase-tm).
local function round_half_even(x)
  local f = math.floor(x)
  local d = x - f
  if d > 0.5 then return f + 1 end
  if d < 0.5 then return f end
  return (f % 2 == 0) and f or f + 1
end
M._round = round_half_up
M._round_even = round_half_even

--- Nearest-neighbour dH (kcal/mol), dS (cal/mol/K) and Tm (°C) at 1 M NaCl
--- for total primer concentration conc_nM (default 500). seq must be valid.
function M.raw_nn(seq, conc_nM)
  seq = seq:upper()
  conc_nM = conc_nM or 500
  local dH, dS = 0, 0
  for i = 1, #seq - 1 do
    local p = NN[seq:sub(i, i + 1)]
    dH = dH + p[1]
    dS = dS + p[2]
  end
  local a, b = TERM[seq:sub(1, 1)], TERM[seq:sub(-1)]
  dH = dH + a[1] + b[1]
  dS = dS + a[2] + b[2]
  local tm = (dH * 1000) / (dS + R * math.log(conc_nM * 1e-9)) - K0
  return dH, dS, tm
end

local function gc_fraction(seq)
  local _, n = seq:gsub("[GCgc]", "")
  return n / #seq
end

local function owczarzy2004(raw, fgc, salt_mM)
  local ln = math.log(salt_mM / 1000)
  local corr = (4.29e-5 * fgc - 3.95e-5) * ln + 9.40e-6 * ln * ln
  return 1 / (1 / (raw + K0) + corr) - K0
end

local function owczarzy2008(raw, fgc, n, mono_mM, mg_mM)
  if mg_mM <= 0 then return owczarzy2004(raw, fgc, mono_mM) end
  local mg = mg_mM / 1000
  local a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5
  if mono_mM > 0 then
    local mono = mono_mM / 1000
    local ratio = math.sqrt(mg) / mono
    if ratio < 0.22 then return owczarzy2004(raw, fgc, mono_mM) end
    if ratio < 6.0 then
      local lnm = math.log(mono)
      a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mono) * lnm)
      d = 1.42e-5 * (1.279 - 4.03e-3 * lnm - 8.03e-3 * lnm * lnm)
      g = 8.31e-5 * (0.486 - 0.258 * lnm + 5.25e-3 * lnm * lnm * lnm)
    end
  end
  local lg = math.log(mg)
  local corr = a + b * lg + fgc * (c + d * lg) + (1 / (2 * (n - 1))) * (e + f * lg + g * lg * lg)
  return 1 / (1 / (raw + K0) + corr) - K0
end

local METHODS = { owczarzy = true, schildkraut = true, owczarzy2008 = true }

local function positive(v, name)
  if v == nil then return true end
  if type(v) ~= "number" or v ~= v or v <= 0 or v == math.huge then
    return false, ("%s must be a positive number"):format(name)
  end
  return true
end

local function compat_mode(opts)
  local c = opts.compat
  if c == nil or c == "neb" then return false end
  if c == "polymerase-tm" then return true end
  return nil, ("unknown compat %q (use \"neb\" or \"polymerase-tm\")"):format(tostring(c))
end

--- Full Tm calculation. Returns a table (see module docs) or nil, err.
function M.tm_detail(seq, opts)
  opts = opts or {}
  local ok, err = M.valid(seq)
  if not ok then return nil, err end
  seq = seq:upper()
  local compat, cerr = compat_mode(opts)
  if compat == nil then return nil, cerr end
  local buffers = compat and M.BUFFERS_PTM or M.BUFFERS
  local p = M.polymerase(opts.polymerase or "q5")
  if not p then return nil, ("unknown polymerase %q"):format(tostring(opts.polymerase)) end

  local salt
  if opts.salt_mM ~= nil then
    ok, err = positive(opts.salt_mM, "salt_mM")
    if not ok then return nil, err end
    salt = opts.salt_mM
  elseif opts.buffer ~= nil then
    salt = buffers[norm(opts.buffer)]
    if not salt then return nil, ("unknown buffer %q"):format(tostring(opts.buffer)) end
  else
    salt = buffers[p.buffer]
  end
  ok, err = positive(opts.conc_nM, "conc_nM")
  if not ok then return nil, err end
  local conc = opts.conc_nM or p.conc_nM
  local method = opts.method or p.method
  if not METHODS[method] then return nil, ("unknown method %q"):format(tostring(method)) end
  local dmso = opts.dmso or 0
  if type(dmso) ~= "number" or dmso ~= dmso or dmso < 0 then return nil, "dmso must be a number >= 0" end
  local mg = opts.mg_mM
  if mg ~= nil and (type(mg) ~= "number" or mg ~= mg or mg < 0) then return nil, "mg_mM must be a number >= 0" end

  local raw, tm, conc_eff
  if method == "schildkraut" then
    conc_eff = conc / 4 -- NEB's Phusion convention
    local _, _, r = M.raw_nn(seq, conc_eff)
    raw = r
    tm = raw + 16.6 * math.log10(salt / 1000)
  else
    conc_eff = conc
    local _, _, r = M.raw_nn(seq, conc)
    raw = r
    local fgc = gc_fraction(seq)
    if method == "owczarzy2008" then
      tm = owczarzy2008(raw, fgc, #seq, salt, mg or 2)
    else
      tm = owczarzy2004(raw, fgc, salt)
    end
  end
  tm = tm - 0.6 * dmso
  return {
    tm = tm,
    tm_rounded = compat and round_half_even(tm) or round_half_up(tm),
    raw = raw,
    salt_mM = salt,
    conc_nM = conc,
    conc_eff_nM = conc_eff,
    method = method,
    mg_mM = method == "owczarzy2008" and (mg or 2) or nil,
    dmso = dmso,
    polymerase = p.key,
    compat = compat and "polymerase-tm" or "neb",
  }
end

--- Integer Tm (°C, NEB rounding) or nil, err.
function M.tm(seq, opts)
  local d, err = M.tm_detail(seq, opts)
  if not d then return nil, err end
  return d.tm_rounded
end

--- Recommended annealing temperature for a primer pair.
--- Returns ta, tm1, tm2, info  (or nil, err). tm1/tm2 are the integer buffer
--- Tms without DMSO; DMSO lowers Ta after the cap. info = { capped, rule
--- ("offset" | "phusion" | "short"), note, dmso, ta_decimal (Ta to 0.1 °C as
--- NEB's long format reports it; unrounded in polymerase-tm compat mode) }.
function M.ta(seq1, seq2, opts)
  opts = opts or {}
  local compat, cerr = compat_mode(opts)
  if compat == nil then return nil, cerr end
  local p = M.polymerase(opts.polymerase or "q5")
  if not p then return nil, ("unknown polymerase %q"):format(tostring(opts.polymerase)) end
  local dmso = opts.dmso or 0
  if type(dmso) ~= "number" or dmso ~= dmso or dmso < 0 then return nil, "dmso must be a number >= 0" end
  local o = {}
  for k, v in pairs(opts) do o[k] = v end
  o.dmso = nil
  local d1, e1 = M.tm_detail(seq1, o)
  if not d1 then return nil, "primer 1: " .. e1 end
  local d2, e2 = M.tm_detail(seq2, o)
  if not d2 then return nil, "primer 2: " .. e2 end
  local t1, t2 = d1.tm_rounded, d2.tm_rounded

  -- NEB: rule on the unrounded lower Tm; polymerase-tm: on the integer Tms
  local low = compat and math.min(t1, t2) or math.min(d1.tm, d2.tm)
  local shown = compat and tostring(low) or ("%.1f"):format(low)
  local minlen = math.min(#seq1, #seq2)
  local x, rule, note
  if p.ta_rule == "phusion" then
    x = 0.93 * low + 7.5
    rule = "phusion"
    note = ("Ta = 0.93 x %s + 7.5 (Phusion rule)"):format(shown)
  elseif minlen >= (p.min_len or 8) then
    local off = p.ta_offset or 0
    x = low + off
    rule = "offset"
    note = ("Ta = lower Tm %s %s %d °C"):format(shown, off < 0 and "-" or "+", math.abs(off))
  else
    x = low
    rule = "short"
    note = ("Ta = lower Tm %s (primer shorter than %d nt: no offset)"):format(shown, p.min_len)
  end
  local capped = x > p.ta_cap
  if capped then
    x = p.ta_cap
    note = note .. (", capped at %d °C"):format(p.ta_cap)
  end
  if dmso > 0 then
    x = x - 0.6 * dmso
    note = note .. (", -%.1f °C for %g%% DMSO"):format(0.6 * dmso, dmso)
  end
  local ta, dec
  if compat then
    ta, dec = round_half_even(x), x
  else
    dec = round_half_up(x * 10) / 10
    ta = round_half_up(dec)
  end
  return ta, t1, t2, { capped = capped, rule = rule, note = note, dmso = dmso, ta_decimal = dec }
end

return M
