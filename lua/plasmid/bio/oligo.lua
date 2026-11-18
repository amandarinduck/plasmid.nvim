--- Oligo (primer) analysis: secondary-structure and sequence-quality estimates.
---
--- Pure Lua, no dependencies. Everything here is an *estimate* from a
--- two-state nearest-neighbour (NN) model; see oligo_report.md for how the
--- numbers compare with primer3 (thal) on several thousand oligos.
---
--- Model (implemented from the published tables, not from any program's code):
---   * Watson–Crick stacks: SantaLucia & Hicks (2004) Annu. Rev. Biophys.
---     Biomol. Struct. 33:415, Table 1 (duplex initiation +0.2 / -5.7,
---     terminal A·T penalty +2.2 / +6.9).
---   * single internal mismatches: Allawi & SantaLucia (1997, 1998 a-c),
---     Peyret et al. (1999); terminal mismatches: SantaLucia & Peyret (2001);
---     dangling ends: Bommarito, Peyret & SantaLucia (2000). Values transcribed from the
---     Biopython Bio.SeqUtils.MeltingTemp tables DNA_IMM1, DNA_TMM1, DNA_DE1,
---     which cite those papers.
---   * loop initiation (internal, bulge, hairpin; purely entropic, dH = 0):
---     SantaLucia & Hicks (2004) Table 4, with the Jacobson–Stockmayer
---     extrapolation dG(n) = dG(x) + 2.44·R·310.15·ln(n/x) between/after
---     tabulated lengths; internal-loop asymmetry +0.3 kcal/mol per nt.
---   * salt: Na-equivalent [Na+] + 120·sqrt([Mg2+] - [dNTP]) (von Ahsen et
---     al. 2001), entropy correction 0.368·ln[Na+] per nearest-neighbour step
---     (SantaLucia 1998).
---
--- Search: exact dynamic programming over single-helix structures — stacks,
--- 1x1 mismatches, bulges and internal loops up to opts.max_loop (default 30)
--- unpaired bases per loop, plus dangling ends / terminal mismatches at helix
--- ends. Dimers are one antiparallel duplex between the two strands
--- (bimolecular Tm with Ct/4, or Ct for self-dimers); hairpin folds are one or
--- more consecutive stem-loops, each stem possibly interrupted, each closing a
--- hairpin loop of 3..max_loop nt (unimolecular Tm = dH/dS).
--- The structure minimising dG at opts.temp_c is reported (primer3's thal
--- chooses structures differently, so it sometimes reports a less stable
--- structure than this search finds; see the report).
--- Units: dg/dh in kcal/mol, ds in cal/(K mol), tm in °C.
--- Approximations (vs. a full thal/mfold-style model): no tri/tetraloop
--- bonus table, no special hairpin-mismatch table (duplex terminal-mismatch
--- values are used), internal loops scored as initiation + asymmetry +
--- terminal mismatches at both closing pairs, tandem mismatches scored as
--- 2x2 internal loops, no multi-branch (closed multiloop) structures, no
--- G·T pairs as helix pairs (only as mismatches), no symmetry correction,
--- non-ACGT letters never pair.
local M = {}

local R = 1.98720425864083 -- cal / (K mol)
local KELVIN = 273.15

M.defaults = { mv_mM = 50, dv_mM = 1.5, dntp_mM = 0.6, dna_nM = 50, temp_c = 37, max_loop = 30 }

--- Default warning thresholds for M.assess (all overridable via opts.thresholds).
--- These are conventional cut-offs in common use for primer QC, not measured
--- limits, and no single authority defines them; adjust them to your assay.
--- Only struct_tm and end_stability correspond to a named program's default
--- (primer3).
---   dimer_dg = -9      a commonly used cut-off for self-/hetero-dimers (dG
---                      weaker than about -9 kcal/mol); a convention, not a
---                      measured limit — set it to suit your assay
---   dimer_3p_dg = -5   3'-end dimers are judged more harshly than internal ones;
---   hairpin_dg = -2    hairpins below about -2 are usually flagged (3' end in the
---                      stem) / -3 (internal); we apply the stricter -2 when the
---                      hairpin involves the 3'-terminal 5 nt, else -3
---   hairpin_int_dg = -3
---   struct_tm = 47     primer3 defaults PRIMER_MAX_HAIRPIN_TH / _SELF_ANY_TH /
---                      _SELF_END_TH = 47 °C
---   gc_min, gc_max     0.40 / 0.60, the conventional primer-design window
---                      (0.35-0.80 is also commonly quoted)
---   len_min, len_max   18 / 30 nt, conventional PCR primer lengths
---   max_gc_clamp = 3   no more than 3 G/C in the last 5 bases;
---   min_gc_clamp = 1   at least one G/C in the last 5 bases
---   max_run = 4        homopolymer runs longer than 4 (
---                      primer3's PRIMER_MAX_POLY_X default is 5)
---   max_repeat = 4     dinucleotide repeats longer than 4 units
---   end_stability = -9 3'-pentamer dG more negative than -9 kcal/mol (primer3's
---                      PRIMER_MAX_END_STABILITY; its default 100 disables the
---                      check, 9 is the commonly used value)
M.thresholds = {
  dimer_dg = -9, dimer_3p_dg = -5, hairpin_dg = -2, hairpin_int_dg = -3, struct_tm = 47,
  gc_min = 0.40, gc_max = 0.60, len_min = 18, len_max = 30,
  min_gc_clamp = 1, max_gc_clamp = 3, max_run = 4, max_repeat = 4, end_stability = -9,
}

---------------------------------------------------------------------------
-- parameter tables: key "XY/AB" = 5'-XY-3' over 3'-AB-5'; {dH kcal/mol, dS cal/(K mol)}
---------------------------------------------------------------------------
-- SantaLucia & Hicks 2004 Table 1
local NN = {
  ["AA/TT"] = { -7.6, -21.3 }, ["AT/TA"] = { -7.2, -20.4 }, ["TA/AT"] = { -7.2, -21.3 },
  ["CA/GT"] = { -8.5, -22.7 }, ["GT/CA"] = { -8.4, -22.4 }, ["CT/GA"] = { -7.8, -21.0 },
  ["GA/CT"] = { -8.2, -22.2 }, ["CG/GC"] = { -10.6, -27.2 }, ["GC/CG"] = { -9.8, -24.4 },
  ["GG/CC"] = { -8.0, -19.9 },
}
local INIT = { 0.2, -5.7 }
local AT_PEN = { 2.2, 6.9 }
-- single internal mismatches (Allawi & SantaLucia 1997-98; Peyret et al. 1999)
local IMM = {
  ["AG/TT"] = { 1, 0.9 }, ["AT/TG"] = { -2.5, -8.3 }, ["CG/GT"] = { -4.1, -11.7 }, ["CT/GG"] = { -2.8, -8 },
  ["GG/CT"] = { 3.3, 10.4 }, ["GG/TT"] = { 5.8, 16.3 }, ["GT/CG"] = { -4.4, -12.3 }, ["GT/TG"] = { 4.1, 9.5 },
  ["TG/AT"] = { -0.1, -1.7 }, ["TG/GT"] = { -1.4, -6.2 }, ["TT/AG"] = { -1.3, -5.3 }, ["AA/TG"] = { -0.6, -2.3 },
  ["AG/TA"] = { -0.7, -2.3 }, ["CA/GG"] = { -0.7, -2.3 }, ["CG/GA"] = { -4, -13.2 }, ["GA/CG"] = { -0.6, -1 },
  ["GG/CA"] = { 0.5, 3.2 }, ["TA/AG"] = { 0.7, 0.7 }, ["TG/AA"] = { 3, 7.4 }, ["AC/TT"] = { 0.7, 0.2 },
  ["AT/TC"] = { -1.2, -6.2 }, ["CC/GT"] = { -0.8, -4.5 }, ["CT/GC"] = { -1.5, -6.1 }, ["GC/CT"] = { 2.3, 5.4 },
  ["GT/CC"] = { 5.2, 13.5 }, ["TC/AT"] = { 1.2, 0.7 }, ["TT/AC"] = { 1, 0.7 }, ["AA/TC"] = { 2.3, 4.6 },
  ["AC/TA"] = { 5.3, 14.6 }, ["CA/GC"] = { 1.9, 3.7 }, ["CC/GA"] = { 0.6, -0.6 }, ["GA/CC"] = { 5.2, 14.2 },
  ["GC/CA"] = { -0.7, -3.8 }, ["TA/AC"] = { 3.4, 8 }, ["TC/AA"] = { 7.6, 20.2 }, ["AA/TA"] = { 1.2, 1.7 },
  ["CA/GA"] = { -0.9, -4.2 }, ["GA/CA"] = { -2.9, -9.8 }, ["TA/AA"] = { 4.7, 12.9 }, ["AC/TC"] = { 0, -4.4 },
  ["CC/GC"] = { -1.5, -7.2 }, ["GC/CC"] = { 3.6, 8.9 }, ["TC/AC"] = { 6.1, 16.4 }, ["AG/TG"] = { -3.1, -9.5 },
  ["CG/GG"] = { -4.9, -15.3 }, ["GG/CG"] = { -6, -15.8 }, ["TG/AG"] = { 1.6, 3.6 }, ["AT/TT"] = { -2.7, -10.8 },
  ["CT/GT"] = { -5, -15.8 }, ["GT/CT"] = { -2.2, -8.4 }, ["TT/AT"] = { 0.2, -1.5 },
}
-- terminal mismatches (SantaLucia & Peyret 2001, WO 01/94611, as tabulated in
-- Biopython DNA_TMM1), mismatch on the right of each key. "GG/CG" is listed
-- there as {-0.7, -19.2} (dG37 +5.25, far outside the -1.5..-0.1 of every
-- other entry and destabilising every C·G pair flanked by G/G); we assume a
-- misplaced decimal point in dH and use -7.0 (dG37 -1.05), which is also what
-- the residuals against primer3 suggest (see oligo_report.md).
local TMM = {
  ["AA/TA"] = { -3.1, -7.8 }, ["TA/AA"] = { -2.5, -6.3 }, ["CA/GA"] = { -4.3, -10.7 }, ["GA/CA"] = { -8, -22.5 },
  ["AC/TC"] = { -0.1, 0.5 }, ["TC/AC"] = { -0.7, -1.3 }, ["CC/GC"] = { -2.1, -5.1 }, ["GC/CC"] = { -3.9, -10.6 },
  ["AG/TG"] = { -1.1, -2.1 }, ["TG/AG"] = { -1.1, -2.7 }, ["CG/GG"] = { -3.8, -9.5 }, ["GG/CG"] = { -7.0, -19.2 },
  ["AT/TT"] = { -2.4, -6.5 }, ["TT/AT"] = { -3.2, -8.9 }, ["CT/GT"] = { -6.1, -16.9 }, ["GT/CT"] = { -7.4, -21.2 },
  ["AA/TC"] = { -1.6, -4 }, ["AC/TA"] = { -1.8, -3.8 }, ["CA/GC"] = { -2.6, -5.9 }, ["CC/GA"] = { -2.7, -6 },
  ["GA/CC"] = { -5, -13.8 }, ["GC/CA"] = { -3.2, -7.1 }, ["TA/AC"] = { -2.3, -5.9 }, ["TC/AA"] = { -2.7, -7 },
  ["AC/TT"] = { -0.9, -1.7 }, ["AT/TC"] = { -2.3, -6.3 }, ["CC/GT"] = { -3.2, -8 }, ["CT/GC"] = { -3.9, -10.6 },
  ["GC/CT"] = { -4.9, -13.5 }, ["GT/CC"] = { -3, -7.8 }, ["TC/AT"] = { -2.5, -6.3 }, ["TT/AC"] = { -0.7, -1.2 },
  ["AA/TG"] = { -1.9, -4.4 }, ["AG/TA"] = { -2.5, -5.9 }, ["CA/GG"] = { -3.9, -9.6 }, ["CG/GA"] = { -6, -15.5 },
  ["GA/CG"] = { -4.3, -11.1 }, ["GG/CA"] = { -4.6, -11.4 }, ["TA/AG"] = { -2, -4.7 }, ["TG/AA"] = { -2.4, -5.8 },
  ["AG/TT"] = { -3.2, -8.7 }, ["AT/TG"] = { -3.5, -9.4 }, ["CG/GT"] = { -3.8, -9 }, ["CT/GG"] = { -6.6, -18.7 },
  ["GG/CT"] = { -5.7, -15.9 }, ["GT/CG"] = { -5.9, -16.1 }, ["TG/AT"] = { -3.9, -10.5 }, ["TT/AG"] = { -3.6, -9.8 },
}
-- dangling ends (Bommarito et al. 2000); "." marks the missing partner
local DE = {
  ["AA/.T"] = { 0.2, 2.3 }, ["AC/.G"] = { -6.3, -17.1 }, ["AG/.C"] = { -3.7, -10 }, ["AT/.A"] = { -2.9, -7.6 },
  ["CA/.T"] = { 0.6, 3.3 }, ["CC/.G"] = { -4.4, -12.6 }, ["CG/.C"] = { -4, -11.9 }, ["CT/.A"] = { -4.1, -13 },
  ["GA/.T"] = { -1.1, -1.6 }, ["GC/.G"] = { -5.1, -14 }, ["GG/.C"] = { -3.9, -10.9 }, ["GT/.A"] = { -4.2, -15 },
  ["TA/.T"] = { -6.9, -20 }, ["TC/.G"] = { -4, -10.9 }, ["TG/.C"] = { -4.9, -13.8 }, ["TT/.A"] = { -0.2, -0.5 },
  [".A/AT"] = { -0.7, -0.8 }, [".C/AG"] = { -2.1, -3.9 }, [".G/AC"] = { -5.9, -16.5 }, [".T/AA"] = { -0.5, -1.1 },
  [".A/CT"] = { 4.4, 14.9 }, [".C/CG"] = { -0.2, -0.1 }, [".G/CC"] = { -2.6, -7.4 }, [".T/CA"] = { 4.7, 14.2 },
  [".A/GT"] = { -1.6, -3.6 }, [".C/GG"] = { -3.9, -11.2 }, [".G/GC"] = { -3.2, -10.4 }, [".T/GA"] = { -4.1, -13.1 },
  [".A/TT"] = { 2.9, 10.4 }, [".C/TG"] = { -4.4, -13.1 }, [".G/TC"] = { -5.2, -15 }, [".T/TA"] = { -3.8, -12.6 },
}
-- loop initiation dG37 (kcal/mol), SantaLucia & Hicks 2004 Table 4
local LOOP_N = { 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30 }
local LOOP_TAB = {
  internal = { 3.2, 3.6, 4.0, 4.4, 4.6, 4.8, 4.9, 4.9, 5.2, 5.4, 5.6, 5.8, 5.9, 6.3, 6.6 },
  bulge = { 3.1, 3.2, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.8, 5.0, 5.2, 5.3, 5.6, 5.9 },
  hairpin = { 3.5, 3.5, 3.3, 4.0, 4.2, 4.3, 4.5, 4.6, 5.0, 5.1, 5.3, 5.5, 5.7, 6.1, 6.3 },
}
local BULGE_SMALL = { 4.0, 2.9 } -- bulges of 1 and 2 nt
local ASYM = 0.3 -- kcal/mol per nt of internal-loop asymmetry
-- end-stability pentamer: SantaLucia 1998 unified dG37 (same dG37 as 2004)
-- with initiation 1.96 + 0.05 per terminal A·T (the 2004 form of the 1998 terms)
M.params = { NN = NN, IMM = IMM, TMM = TMM, DE = DE, INIT = INIT, AT_PEN = AT_PEN, LOOP_TAB = LOOP_TAB }

---------------------------------------------------------------------------
-- table compilation: flat arrays indexed by base codes (0 = none, A1 C2 G3 T4)
---------------------------------------------------------------------------
local CODE = { A = 1, C = 2, G = 3, T = 4, U = 4, ["."] = 0 }

local function key4(x, y, a, b) return ((x * 5 + y) * 5 + a) * 5 + b end

local H_NN, S_NN, H_MM, S_MM, H_TM, S_TM, H_DE, S_DE = {}, {}, {}, {}, {}, {}, {}, {}
local function compile(src, H, S)
  for k, v in pairs(src) do
    local x, y, a, b = CODE[k:sub(1, 1)], CODE[k:sub(2, 2)], CODE[k:sub(4, 4)], CODE[k:sub(5, 5)]
    -- the same window read from the other strand: "XY/AB" == "BA/YX"
    H[key4(x, y, a, b)], S[key4(x, y, a, b)] = v[1] * 1000, v[2]
    H[key4(b, a, y, x)], S[key4(b, a, y, x)] = v[1] * 1000, v[2]
  end
end

--- Interpolated loop initiation dG37 for a loop of n nt (kind: internal/bulge/hairpin)
local function loop_dg(kind, n)
  if kind == "bulge" and n <= 2 then return BULGE_SMALL[n] end
  local tab = LOOP_TAB[kind]
  if n < 3 then n = 3 end
  local lo = 1
  for i = 1, #LOOP_N do if LOOP_N[i] <= n then lo = i end end
  local x = LOOP_N[lo]
  return tab[lo] + 2.44 * R * 310.15 / 1000 * math.log(n / x)
end

local LOOP_S = { internal = {}, bulge = {}, hairpin = {} }
local bound_cache = {} -- see loop_bound()
local function compile_all()
  bound_cache = {}
  for _, t in ipairs({ H_NN, S_NN, H_MM, S_MM, H_TM, S_TM, H_DE, S_DE }) do
    for k = 0, 624 do t[k] = false end -- dense (array part); false = no value
  end
  compile(NN, H_NN, S_NN); compile(IMM, H_MM, S_MM); compile(TMM, H_TM, S_TM); compile(DE, H_DE, S_DE)
  for kind, t in pairs(LOOP_S) do
    for n = 1, 1000 do t[n] = -loop_dg(kind, n) * 1000 / 310.15 end
  end
end
compile_all()
--- Recompile the flat tables after editing M.params (development hook).
M._recompile = compile_all

---------------------------------------------------------------------------
-- helpers
---------------------------------------------------------------------------
local function encode(seq)
  local s = seq:upper():gsub("[^%a]", "")
  local c = {}
  for i = 1, #s do c[i] = CODE[s:sub(i, i)] or 0 end -- non-ACGT: never pairs
  return c, s
end

local function options(opts)
  local o = {}
  for k, v in pairs(M.defaults) do o[k] = v end
  if opts then for k, v in pairs(opts) do o[k] = v end end
  return o
end

--- Na+-equivalent concentration (M), von Ahsen et al. 2001
function M.na_eq(o)
  o = options(o)
  local free_mg = math.max(0, o.dv_mM - o.dntp_mM)
  return (o.mv_mM + 120 * math.sqrt(free_mg)) / 1000
end

local function salt_ds(o)
  return 0.368 * math.log(M.na_eq(o))
end

-- end term for a helix end at pair (p, q): "x p / y q", 0 = absent neighbour
local function end_term(x, p, y, q)
  if x > 0 and y > 0 then
    local k = key4(x, p, y, q)
    local h = H_TM[k]
    if h then return h, S_TM[k] end
    return 0, 0 -- the flanking bases could pair: no terminal-mismatch value
  elseif x > 0 or y > 0 then
    local k = key4(x, p, y, q)
    local h = H_DE[k]
    if h then return h, S_DE[k] end
  end
  return 0, 0
end

local function at_pen(p, q)
  if p == 1 or p == 4 then return AT_PEN[1] * 1000, AT_PEN[2] end
  return 0, 0
end

---------------------------------------------------------------------------
-- dimer DP: a (top, 5'->3') against b (bottom, given 5'->3', read 3'->5')
---------------------------------------------------------------------------
-- loop between outer pair (a[p], r[q]) and inner pair (a[i], r[k]) with l1/l2
-- unpaired bases on top/bottom; returns dH (cal/mol), dS (cal/(K mol)) of the
-- loop, including the salt term of any nearest-neighbour stacks it contains
local function loop_energy(ap, ap1, rq, rq1, ai1, ai, rk1, rk, l1, l2, salt)
  if l1 == 0 or l2 == 0 then
    local L = l1 + l2
    local s = LOOP_S.bulge[L]
    if L == 1 then
      local k = key4(ap, ai, rq, rk)
      return H_NN[k], S_NN[k] + s + salt
    end
    local h1, s1 = at_pen(ap, rq)
    local h2, s2 = at_pen(ai, rk)
    return h1 + h2, s1 + s2 + s
  end
  if l1 == 1 and l2 == 1 then
    local k1, k2 = key4(ap, ap1, rq, rq1), key4(ai1, ai, rk1, rk)
    local h1, h2 = H_MM[k1], H_MM[k2]
    if h1 and h2 then return h1 + h2, S_MM[k1] + S_MM[k2] + 2 * salt end
  end
  local n = l1 + l2
  local s = LOOP_S.internal[n] - ASYM * math.abs(l1 - l2) * 1000 / 310.15
  local k1, k2 = key4(ap, ap1, rq, rq1), key4(rk, rk1, ai, ai1)
  local h1, s1 = H_TM[k1] or 0, S_TM[k1] or 0
  local h2, s2 = H_TM[k2] or 0, S_TM[k2] or 0
  return h1 + h2, s1 + s2 + s
end

local INF = math.huge
-- DP scratch arrays, reused across calls to avoid per-call allocation/GC
local W1, W2, W3, W4, W5, W6 = {}, {}, {}, {}, {}, {}

-- Structure choice: minimise dG at temp_c (exact DP over the structure class).
-- Scores are maximised, score = -dG = T*dS - dH. (A greedy "highest Tm per DP
-- cell" variant was tried and agreed less well with primer3; see report.)
local function make_score(o)
  local T = o.temp_c + KELVIN
  return function(h, s) return T * s - h end
end

-- Upper bound on the score (-dG) any single bulge/internal loop can add at
-- temperature T: lets the DP skip predecessors that cannot win (exact pruning).
local function loop_bound(T, salt)
  local key = T .. ":" .. salt
  if bound_cache[key] then return bound_cache[key] end
  local function best(H, S)
    local b = -INF
    for k = 0, 624 do if H[k] then local x = T * S[k] - H[k]; if x > b then b = x end end end
    return b
  end
  local mm, nn, tm = best(H_MM, S_MM), best(H_NN, S_NN), best(H_TM, S_TM)
  local at = math.max(0, T * AT_PEN[2] - AT_PEN[1] * 1000)
  local b = 2 * (mm + T * salt)
  b = math.max(b, nn + T * (LOOP_S.bulge[1] + salt))
  for L = 2, #LOOP_S.bulge do
    b = math.max(b, T * LOOP_S.bulge[L] + 2 * at)
    b = math.max(b, T * LOOP_S.internal[L] + 2 * math.max(0, tm))
  end
  b = b + 1e-6
  bound_cache[key] = b
  return b
end

local function rc_term(o, self)
  local ct = o.dna_nM * 1e-9
  return R * math.log(self and ct or ct / 4)
end

-- One antiparallel duplex of a (top) with r (= bottom strand reversed, i.e.
-- read 3'->5'). Returns the best structure overall and the best one whose
-- last pair is a's 3'-terminal base (anchored).
local function duplex(a, r, o, rc)
  local n, m = #a, #r
  local T = o.temp_c + KELVIN
  local salt = salt_ds(o)
  local maxloop = math.min(o.max_loop, 1000)
  local score = make_score(o)
  local LB = loop_bound(T, salt)
  local FH, FS, FX, BP = W1, W2, W3, W4 -- reused scratch arrays
  local NEG = -INF
  for c = 1, n * m do FX[c], FH[c], FS[c], BP[c] = NEG, 0, 0, 0 end
  local best, bc, bh, bs = NEG, 0, 0, 0
  local abest, ac, ah, as = NEG, 0, 0, 0
  local ih, is = INIT[1] * 1000, INIT[2]
  -- for every row p, the columns q (descending) where a[p] pairs with r[q]
  local WQ = {}
  for p = 1, n do
    local t, ap = {}, a[p]
    for q = m, 1, -1 do if ap > 0 and ap + r[q] == 5 then t[#t + 1] = q end end
    WQ[p] = t
  end
  for i = 1, n do
    local ai = a[i]
    local ai1 = a[i - 1] or 0
    local xa = a[i + 1] or 0
    for k = 1, m do
      local rk = r[k]
      if ai > 0 and ai + rk == 5 then
        local c = (i - 1) * m + k
        local rk1 = r[k - 1] or 0
        -- start a helix here (left end)
        local h0, s0 = end_term(ai1, ai, rk1, rk)
        local ph, ps = at_pen(ai, rk)
        local gh, gs = ih + h0 + ph, is + s0 + ps
        local gx, gb = score(gh, gs), 0
        if i > 1 and k > 1 then
          -- stack on (i-1, k-1)
          local pc = c - m - 1
          if FX[pc] > NEG then
            local kk = key4(ai1, ai, rk1, rk)
            local h, s = FH[pc] + H_NN[kk], FS[pc] + S_NN[kk] + salt
            local x = score(h, s)
            if x > gx then gh, gs, gx, gb = h, s, x, pc end
          end
          -- bulges and internal loops closed by (p, q) and (i, k)
          local pmin = i - 1 - maxloop
          if pmin < 1 then pmin = 1 end
          for p = i - 1, pmin, -1 do
            local l1 = i - p - 1
            local ap, ap1 = a[p], a[p + 1]
            local qmin = k - 1 - (maxloop - l1)
            if qmin < 1 then qmin = 1 end
            local base = (p - 1) * m
            local wq = WQ[p]
            for t = 1, #wq do
              local q = wq[t]
              if q < qmin then break end
              local pc2 = base + q
              if q < k and FX[pc2] + LB > gx then
                local l2 = k - q - 1
                if l1 + l2 > 0 then
                  local h, s = loop_energy(ap, ap1, r[q], r[q + 1], ai1, ai, rk1, rk, l1, l2, salt)
                  h, s = FH[pc2] + h, FS[pc2] + s
                  local x = score(h, s)
                  if x > gx then gh, gs, gx, gb = h, s, x, pc2 end
                end
              end
            end
          end
        end
        FH[c], FS[c], FX[c], BP[c] = gh, gs, gx, gb
        -- close the helix here (right end): window "ai x / rk y"
        local y = r[k + 1] or 0
        local h1, s1 = 0, 0
        local kk = key4(ai, xa, rk, y)
        if xa > 0 and y > 0 then
          h1, s1 = H_TM[kk] or 0, S_TM[kk] or 0
        elseif xa > 0 or y > 0 then
          h1, s1 = H_DE[kk] or 0, S_DE[kk] or 0
        end
        local th, ts = gh + h1 + ph, gs + s1 + ps
        local tx = score(th, ts)
        if tx > best then best, bc, bh, bs = tx, c, th, ts end
        if i == n and tx > abest then abest, ac, ah, as = tx, c, th, ts end
      end
    end
  end
  local function trace(c0, h, s)
    if c0 == 0 then return nil end
    local prs = {}
    local c = c0
    while c > 0 do
      local i = math.floor((c - 1) / m) + 1
      table.insert(prs, 1, { i, c - (i - 1) * m })
      c = BP[c]
    end
    return { dh = h, ds = s, dg = h - T * s, pairs = prs }
  end
  return trace(bc, bh, bs), trace(ac, ah, as)
end

local function tm_bimolecular(dh, ds, rc)
  local denom = ds + rc
  if dh >= 0 or denom >= 0 then return nil end
  return dh / denom - KELVIN
end

local function finish_dimer(d, na, m, rc)
  if not d then return nil end
  local out = { dg = d.dg / 1000, dh = d.dh / 1000, ds = d.ds, pairs = {} }
  local three = false
  for _, pk in ipairs(d.pairs) do
    local i, j = pk[1], m - pk[2] + 1
    out.pairs[#out.pairs + 1] = { i, j }
    if i == na or j == m then three = true end
  end
  out.three_prime = three
  out.offset = out.pairs[1][1] - out.pairs[1][2]
  out.tm = tm_bimolecular(d.dh, d.ds, rc)
  return out
end

-- returns best dimer, best dimer with a's 3' end paired
local function dimer(sa, sb, o, self)
  local a = encode(sa)
  local b = encode(sb)
  local m = #b
  local r = {}
  for k = 1, m do r[k] = b[m - k + 1] end
  local rc = rc_term(o, self)
  local d, e = duplex(a, r, o, rc)
  return finish_dimer(d, #a, m, rc), finish_dimer(e, #a, m, rc)
end

--- Best self-dimer of seq (two copies, antiparallel).
--- Returns { dg, dh, ds, tm, offset, pairs = {{i, j}, ...}, three_prime } or nil.
--- dg/dh in kcal/mol, ds in cal/(K mol), tm in °C (nil if not defined).
--- pairs: base i of one copy pairs with base j of the other (1-based, 5'->3');
--- offset = i - j of the first pair (0 = the copies are in register);
--- three_prime: a 3'-terminal base of either strand is paired, i.e. the
--- dimer is extensible by a polymerase.
function M.homodimer(seq, opts)
  return (dimer(seq, seq, options(opts), true))
end

--- Best heterodimer of a with b: pairs {i in a, j in b}; same shape as homodimer.
function M.heterodimer(a, b, opts)
  return (dimer(a, b, options(opts), a == b))
end

--- Best dimer in which the 3'-terminal base of `a` is paired with `b` (the
--- structure a polymerase could extend; cf. primer3's 3'-anchored alignment
--- used for PRIMER_MAX_SELF_END_TH / calc_end_stability). b defaults to a.
function M.end_dimer(a, b, opts)
  local _, e = dimer(a, b or a, options(opts), b == nil or b == a)
  return e
end

---------------------------------------------------------------------------
-- hairpin DP
---------------------------------------------------------------------------
--- Best intramolecular fold. Returns
---   { dg, dh, ds, tm, stem = { i, j, len }, loop, pairs, domains } or nil
--- dg/dh/ds/tm describe the whole fold, which may consist of several
--- consecutive stem-loops (`domains`, each { dg, dh, ds, stem, loop, pairs });
--- stem/loop are those of the most stable domain. stem.i/stem.j: its outermost
--- pair (1-based, i < j); stem.len: base pairs in the stem (bulges/internal
--- loops allowed); loop: unpaired bases of its hairpin loop (>= 3).
--- When no fold has dG < 0 the least unstable single hairpin is returned
--- (dg > 0); nil only when no stem-loop is possible at all.
function M.hairpin(seq, opts)
  local o = options(opts)
  local s = encode(seq)
  local n = #s
  local T = o.temp_c + KELVIN
  local salt = salt_ds(o)
  local maxloop = math.min(o.max_loop, 1000)
  local score = make_score(o)
  local LB = loop_bound(T, salt)
  local HH, HS, HX, BP, CH, CS = W1, W2, W3, W4, W5, W6 -- reused scratch arrays
  local NEG = -INF
  for c = 1, n * n do HX[c], HH[c], HS[c], BP[c], CH[c], CS[c] = NEG, 0, 0, 0, false, 0 end
  local best, bc = NEG, nil
  local WQ = {}
  for p = 1, n do
    local t, sp = {}, s[p]
    for q = n, p + 4, -1 do if sp > 0 and sp + s[q] == 5 then t[#t + 1] = q end end
    WQ[p] = t
  end
  for d = 4, n - 1 do
    for i = 1, n - d do
      local j = i + d
      local si, sj = s[i], s[j]
      if si > 0 and si + sj == 5 then
        local c = (i - 1) * n + j
        local si1, sj1 = s[i + 1], s[j - 1]
        local gh, gs, gx, gb = 0, 0, NEG, 0
        -- hairpin loop closed by (i, j)
        local L = j - i - 1
        if L <= maxloop then
          local h, sx
          if L > 3 then
            local kk = key4(si, si1, sj, sj1)
            h, sx = H_TM[kk] or 0, S_TM[kk] or 0
          else
            h, sx = at_pen(si, sj)
          end
          gh, gs = h, sx + LOOP_S.hairpin[L]
          gx = score(gh, gs)
        end
        -- stack on (i+1, j-1)
        local ic = c + n - 1
        if HX[ic] > NEG then
          local kk = key4(si, si1, sj, sj1)
          local h, sx = HH[ic] + H_NN[kk], HS[ic] + S_NN[kk] + salt
          local x = score(h, sx)
          if x > gx then gh, gs, gx, gb = h, sx, x, ic end
        end
        -- bulges / internal loops to an inner pair (p, q)
        local pmax = i + 1 + maxloop
        if pmax > j - 5 then pmax = j - 5 end
        for p = i + 1, pmax do
          local l1 = p - i - 1
          local sp, sp0 = s[p], s[p - 1]
          local qmin = j - 1 - (maxloop - l1)
          if qmin < p + 4 then qmin = p + 4 end
          local base = (p - 1) * n
          local wq = WQ[p]
          for t = 1, #wq do
            local q = wq[t]
            if q < qmin then break end
            local pc = base + q
            if q < j and HX[pc] + LB > gx then
              local l2 = j - q - 1
              if l1 + l2 > 0 then
                local h, sx = loop_energy(si, si1, sj, sj1, sp0, sp, s[q + 1], s[q], l1, l2, salt)
                h, sx = HH[pc] + h, HS[pc] + sx
                local x = score(h, sx)
                if x > gx then gh, gs, gx, gb = h, sx, x, pc end
              end
            end
          end
        end
        if gx > NEG then
          HH[c], HS[c], HX[c], BP[c] = gh, gs, gx, gb
          -- closed as an outermost pair (with its end terms)
          local h1, s1 = end_term(s[i - 1] or 0, si, s[j + 1] or 0, sj)
          local ph, ps = at_pen(si, sj)
          CH[c], CS[c] = gh + h1 + ph, gs + s1 + ps
          local tx = score(CH[c], CS[c])
          if tx > best then best, bc = tx, c end
        end
      end
    end
  end
  if not bc then return nil end
  -- exterior loop: the best set of non-overlapping stem-loop domains
  -- (a sequence can fold into several consecutive hairpins; their energies add)
  local EH, ES, EX, EB = { [0] = 0 }, { [0] = 0 }, { [0] = 0 }, { [0] = false }
  for j = 1, n do
    EH[j], ES[j], EX[j], EB[j] = EH[j - 1], ES[j - 1], EX[j - 1], false
    for i = 1, j - 4 do
      local c = (i - 1) * n + j
      local ch = CH[c]
      if ch then
        local h, sx = EH[i - 1] + ch, ES[i - 1] + CS[c]
        local x = score(h, sx)
        if x > EX[j] then EH[j], ES[j], EX[j], EB[j] = h, sx, x, c end
      end
    end
  end
  local tops = {}
  local j = n
  while j > 0 do
    local c = EB[j]
    if c then
      table.insert(tops, 1, c)
      j = math.floor((c - 1) / n) -- i - 1
    else
      j = j - 1
    end
  end
  local th, ts
  if #tops == 0 then
    tops = { bc } -- nothing is stable: report the least unstable single hairpin
    th, ts = CH[bc], CS[bc]
  else
    th, ts = EH[n], ES[n]
  end
  local domains, prs, main = {}, {}, nil
  for _, top in ipairs(tops) do
    local dp = {}
    local c = top
    while c > 0 do
      local i = math.floor((c - 1) / n) + 1
      dp[#dp + 1] = { i, c - (i - 1) * n }
      prs[#prs + 1] = dp[#dp]
      c = BP[c]
    end
    local dom = {
      dg = (CH[top] - T * CS[top]) / 1000, dh = CH[top] / 1000, ds = CS[top],
      stem = { i = dp[1][1], j = dp[1][2], len = #dp }, loop = dp[#dp][2] - dp[#dp][1] - 1, pairs = dp,
    }
    domains[#domains + 1] = dom
    if not main or dom.dg < main.dg then main = dom end
  end
  local out = {
    dg = (th - T * ts) / 1000, dh = th / 1000, ds = ts, pairs = prs,
    stem = main.stem, loop = main.loop, domains = domains,
  }
  if th < 0 and ts < 0 then out.tm = th / ts - KELVIN end
  return out
end

---------------------------------------------------------------------------
-- sequence-quality measures
---------------------------------------------------------------------------
-- SantaLucia 1998 unified NN dG37 (kcal/mol), used for the 3'-pentamer
local DG37 = {
  AA = -1.00, TT = -1.00, AT = -0.88, TA = -0.58, CA = -1.45, TG = -1.45, GT = -1.44, AC = -1.44,
  CT = -1.28, AG = -1.28, GA = -1.30, TC = -1.30, CG = -2.17, GC = -2.24, GG = -1.84, CC = -1.84,
}

--- dG37 (kcal/mol, negative = stable) of the duplex formed by the 3'-terminal
--- pentamer with its perfect complement — primer3's "3' end stability"
--- (primer3 reports the magnitude, i.e. -dG). 1 M Na+, no salt correction
--- (as in primer3). Shorter sequences use all their bases.
function M.end_stability(seq, opts)
  local s = seq:upper():gsub("[^%a]", ""):gsub("U", "T")
  local tail = s:sub(-5)
  if #tail < 2 then return nil end
  local dg = 1.96
  local f, l = tail:sub(1, 1), tail:sub(-1)
  if f == "A" or f == "T" then dg = dg + 0.05 end
  if l == "A" or l == "T" then dg = dg + 0.05 end
  for i = 1, #tail - 1 do
    local v = DG37[tail:sub(i, i + 1)]
    if not v then return nil end
    dg = dg + v
  end
  return dg
end

--- GC fraction of the A/C/G/T(/U) bases (0 for an empty sequence)
function M.gc(seq)
  local s = seq:upper()
  local gc, n = 0, 0
  for i = 1, #s do
    local ch = s:sub(i, i)
    if ch == "G" or ch == "C" then gc = gc + 1; n = n + 1
    elseif ch == "A" or ch == "T" or ch == "U" then n = n + 1 end
  end
  return n > 0 and gc / n or 0
end

--- number of G/C among the last 5 bases
function M.gc_clamp(seq)
  local s = seq:upper():gsub("[^%a]", ""):sub(-5)
  local _, n = s:gsub("[GC]", "")
  return n
end

--- longest homopolymer run (length in nt)
function M.max_run(seq)
  local s = seq:upper():gsub("[^%a]", "")
  local best, run = 0, 0
  for i = 1, #s do
    if i > 1 and s:byte(i) == s:byte(i - 1) then run = run + 1 else run = 1 end
    if run > best then best = run end
  end
  return best
end

--- longest dinucleotide repeat, in repeat units (e.g. ATATATAT -> 4).
--- Homopolymers (AA) are not dinucleotide repeats (0 for "AAAA"); any other
--- dinucleotide counts 1.
function M.max_repeat(seq)
  local s = seq:upper():gsub("[^%a]", "")
  local best = 0
  for start = 1, 2 do
    local i, run = start, 0
    local prev
    while i + 1 <= #s do
      local di = s:sub(i, i + 1)
      if di:sub(1, 1) ~= di:sub(2, 2) and di == prev then
        run = run + 1
      else
        run = (di:sub(1, 1) ~= di:sub(2, 2)) and 1 or 0
      end
      prev = di
      if run > best then best = run end
      i = i + 2
    end
  end
  return best
end

--- Full assessment with warnings. opts: thermodynamic options (M.defaults)
--- plus opts.thresholds overriding M.thresholds.
function M.assess(seq, opts)
  local o = options(opts)
  local th = {}
  for k, v in pairs(M.thresholds) do th[k] = v end
  if opts and opts.thresholds then for k, v in pairs(opts.thresholds) do th[k] = v end end
  local s = seq:upper():gsub("[^%a]", "")
  local n = #s
  local hd, ed = dimer(s, s, o, true)
  local res = {
    length = n, gc = M.gc(s), gc_clamp = M.gc_clamp(s), max_run = M.max_run(s),
    max_repeat = M.max_repeat(s), end_stability = M.end_stability(s, o),
    hairpin = M.hairpin(s, o), homodimer = hd, end_dimer = ed,
    thresholds = th,
  }
  local w = {}
  local function warn(fmt, ...) w[#w + 1] = fmt:format(...) end
  if n < th.len_min then warn("short: %d nt (< %d)", n, th.len_min) end
  if n > th.len_max then warn("long: %d nt (> %d)", n, th.len_max) end
  if res.gc < th.gc_min or res.gc > th.gc_max then
    warn("GC %.0f%% outside %.0f-%.0f%%", res.gc * 100, th.gc_min * 100, th.gc_max * 100)
  end
  if res.gc_clamp > th.max_gc_clamp then warn("%d G/C in last 5 bases (> %d)", res.gc_clamp, th.max_gc_clamp) end
  if res.gc_clamp < th.min_gc_clamp then warn("no G/C in last 5 bases") end
  if res.max_run > th.max_run then warn("homopolymer run of %d", res.max_run) end
  if res.max_repeat > th.max_repeat then warn("dinucleotide repeat x%d", res.max_repeat) end
  if res.end_stability and res.end_stability < th.end_stability then
    warn("3' end very stable (dG %.1f kcal/mol)", res.end_stability)
  end
  local hp = res.hairpin
  if hp then
    local at3 = hp.stem.j > n - 5
    local lim = at3 and th.hairpin_dg or th.hairpin_int_dg
    if hp.dg < lim then warn("hairpin dG %.1f kcal/mol%s", hp.dg, at3 and " (3' end in stem)" or "") end
    if hp.tm and hp.tm > th.struct_tm then warn("hairpin Tm %.0f °C", hp.tm) end
  end
  local hd = res.homodimer
  if hd then
    if hd.dg < th.dimer_dg then warn("self-dimer dG %.1f kcal/mol", hd.dg) end
    if hd.tm and hd.tm > th.struct_tm then warn("self-dimer Tm %.0f °C", hd.tm) end
  end
  local ed = res.end_dimer
  if ed and ed.dg < th.dimer_3p_dg then warn("3' self-dimer dG %.1f kcal/mol", ed.dg) end
  res.warnings = w
  return res
end

return M
