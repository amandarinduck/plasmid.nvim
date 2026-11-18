#!/usr/bin/env python3
"""Development-time oracles for lua/plasmid/bio/tm.lua (never used at runtime).

Sources
  polymerase-tm  pip install polymerase-tm primer3-py   (MIT; NEB Tm Calculator port)
  biopython      Bio.SeqUtils.MeltingTemp.Tm_NN, DNA_NN3, saltcorr 1/6/7
  neb            https://tmapi.neb.com (/docs/productcodes, POST /tm/batch, POST /tm)

Usage (from the plugin root; NVIM=/path/to/nvim if nvim is not on PATH)
  python tests/tools/tm_oracle.py bulk  [--n 5000] [--out DIR]
        large comparison of the Lua module against polymerase-tm and Biopython;
        writes DIR/bulk_summary.json and DIR/bulk_mismatches.json
  python tests/tools/tm_oracle.py neb   --email-file F | --no-email  [--out DIR]
        queries tmapi.neb.com (default HTTP client headers, ~1 request/s) and
        writes DIR/neb_raw.json (every request/response) and DIR/neb_summary.json
  python tests/tools/tm_oracle.py fixtures [--neb DIR/neb_raw.json]
        writes tests/oracle/expected/tm.json (compact, few hundred cases). NEB
        values are taken from neb_raw.json if given, else carried over from the
        existing tm.json (NEB values are never recomputed or invented).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NVIM = os.environ.get("NVIM", "nvim")
FIXTURE = ROOT / "tests" / "oracle" / "expected" / "tm.json"
SEED = 20261116

# ----------------------------------------------------------------------------
# Lua evaluator: runs the plugin module in headless Neovim on a JSON case list
# ----------------------------------------------------------------------------
LUA_EVAL = r"""
local root, inp, outp = arg[1], arg[2], arg[3]
package.path = root .. "/lua/?.lua;" .. package.path
local M = require("plasmid.bio.tm")
local fh = assert(io.open(inp, "rb")); local cases = vim.json.decode(fh:read("*a"), { luanil = { object = true } }); fh:close()
local out = assert(io.open(outp, "wb"))
local function f(x)
  if x == nil then return "nil" end
  if type(x) == "boolean" then return x and "true" or "false" end
  if type(x) == "number" then return string.format("%.17g", x) end
  return tostring(x)
end
for _, c in ipairs(cases) do
  local o = c.o or {}
  if c.k == "tm" then
    local d, err = M.tm_detail(c.s, o)
    if d then out:write(f(d.tm), "\t", f(d.tm_rounded), "\t", f(d.raw), "\n")
    else out:write("ERR\t", err, "\n") end
  elseif c.k == "ta" then
    local ta, t1, t2, info = M.ta(c.s, c.s2, o)
    if ta then out:write(f(ta), "\t", f(t1), "\t", f(t2), "\t", f(info.capped), "\t", f(info.ta_decimal), "\n")
    else out:write("ERR\t", t1, "\n") end
  elseif c.k == "raw" then
    local dH, dS, tm = M.raw_nn(c.s, c.c)
    out:write(f(dH), "\t", f(dS), "\t", f(tm), "\n")
  end
end
out:close()
"""


def lua_eval(cases: list[dict]) -> list[list]:
    with tempfile.TemporaryDirectory() as td:
        script, inp, outp = (Path(td) / n for n in ("eval.lua", "in.json", "out.tsv"))
        script.write_text(LUA_EVAL)
        inp.write_text(json.dumps(cases))
        subprocess.run([NVIM, "-l", str(script), str(ROOT), str(inp), str(outp)], check=True)
        rows = []
        for line in outp.read_text().splitlines():
            parts = line.split("\t")
            if parts[0] == "ERR":
                rows.append(["ERR", parts[1]])
                continue
            row = []
            for p in parts:
                row.append(None if p == "nil" else (p == "true") if p in ("true", "false") else float(p))
            rows.append(row)
    assert len(rows) == len(cases), (len(rows), len(cases))
    return rows


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------
REAL_PRIMERS = {
    "M13F(-20)": "GTAAAACGACGGCCAGT",
    "M13F(-40)": "GTTTTCCCAGTCACGAC",
    "M13R(-27)": "CAGGAAACAGCTATGAC",
    "M13R(-48)": "AGCGGATAACAATTTCACACAGGA",
    "T7": "TAATACGACTCACTATAGGG",
    "T7term": "GCTAGTTATTGCTCAGCGG",
    "SP6": "ATTTAGGTGACACTATAG",
    "T3": "ATTAACCCTCACTAAAGGG",
    "BGH-R": "TAGAAGGCACAGTCGAGG",
    "CMV-F": "CGCAAATGGGCGGTAGGCGTG",
    "pGEX-5": "GGGCTGGCAAGCCACGTTTGGTG",
    "pGEX-3": "CCGGGAGCTGCATGTGTCAGAGG",
    "EGFP-N": "CGTCGCCGTCCAGCTCGACCAG",
    "EGFP-C": "CATGGTCCTGCTGGAGTTCGTG",
    "LacZ-alpha": "ATGACCATGATTACGCCAAGC",
    "pBAD-F": "ATGCCATAGCATTTTTATCC",
    "pBAD-R": "GATTTAATCTGTATCAGG",
    "AmpR-F": "ATGAGTATTCAACATTTCCGTG",
    "AmpR-R": "TTACCAATGCTTAATCAGTGAG",
    "EcoRI-site8": "GAATTCGC",
}

POLYMERASE_TM_TEST_SEQS = [
    # sequences used in polymerase-tm's own test-suite (tests/test_*.py)
    "ATGTCCCTGCTCTTCTCTCGATGCAA", "ATCGATCGATCG", "TTAGGTGGCGGTACTTGGGTCGATATCAAAGTG",
    "AGCGGATAACAATTTCACACAGGA", "GTAAAACGACGGCCAGT",
]


def load_ptm_test_seqs() -> list[str]:
    """All A/C/G/T literals (>= 8 nt) in polymerase-tm's tests, if the repo checkout is given."""
    d = os.environ.get("POLYMERASE_TM_TESTS")
    seqs = list(POLYMERASE_TM_TEST_SEQS)
    if d and Path(d).is_dir():
        import re
        for p in sorted(Path(d).glob("test_*.py")):
            seqs += re.findall(r'"([ACGT]{8,100})"', p.read_text())
    out, seen = [], set()
    for s in seqs:
        if s not in seen:
            seen.add(s); out.append(s)
    return out


def rand_seq(rng: random.Random, n: int, gc: float) -> str:
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n))


def revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def rand_palindrome(rng: random.Random, n: int, gc: float) -> str:
    half = rand_seq(rng, n // 2, gc)
    return half + revcomp(half)


def make_sequences(n: int, rng: random.Random) -> list[str]:
    seqs = list(REAL_PRIMERS.values()) + load_ptm_test_seqs()
    seqs += ["A" * 8, "T" * 8, "G" * 8, "C" * 8, "A" * 60, "G" * 60, "GC" * 30, "AT" * 30, "ACGT" * 15]
    for i in range(n):
        L = rng.randint(8, 60)
        gc = 0.0 if i % 25 == 0 else 1.0 if i % 25 == 1 else rng.random()
        seqs.append(rand_palindrome(rng, L - L % 2, gc) if i % 50 == 7 else rand_seq(rng, L, gc))
    out, seen = [], set()
    for s in seqs:
        if s not in seen:
            seen.add(s); out.append(s)
    return out


# NEB Tm API v0.19.0 effective monovalent salt per polymerase-tm buffer key,
# fitted exactly from the API's integer Tms (see `neb` / tm_report.md).
NEB_BUFFERS = {
    "standard_taq": 55, "thermopol": 40, "hemo_klentaq": 70, "crimson_taq": 55, "longamp": 100,
    "multiplex": 90, "phusion_hf": 222, "phusion_gc": 222, "phusionflex_hf": 222, "phusionflex_gc": 222,
    "onetaq_std": 54, "onetaq_gc": 80, "q5": 150, "q5mm": 150, "q5u": 170, "q5bd": 150, "longamphs": 100,
}
# products whose Ta the NEB API returns as the lower Tm, uncapped (no rule in the API)
NEB_RULELESS = {"ctaq-0", "ctaq-1", "multiplex-0", "q5u-0", "taqmaster-0", "taqmaster-1", "taqquickload-0", "taqpcrkit-0"}


def neb_unrounded(core, constants, seq, pkey, buffer=None, salt=None, method=None):
    """Unrounded Tm with NEB API salts (polymerase-tm core functions)."""
    poly = constants.POLYMERASES[pkey]
    s = salt if salt is not None else NEB_BUFFERS[buffer or poly["buffer"]]
    m = method or poly.get("tm_method", "owczarzy")
    if m == "schildkraut":
        return core.schildkraut_correction(core.calc_nn_raw(seq, poly["conc"] / 4), s)
    return core.owczarzy_correction(core.calc_nn_raw(seq, poly["conc"]), seq, s)


def neb_ta(constants, u1, u2, pkey, l1, l2, dmso=0.0):
    """NEB Ta procedure: rule on the unrounded lower Tm, cap, DMSO, round to 0.1, round half up."""
    poly = constants.POLYMERASES[pkey]
    low = min(u1, u2)
    if poly["ta_rule"] == "phusion":
        x = 0.93 * low + 7.5
    elif min(l1, l2) >= poly.get("min_len", 8):
        x = low + poly.get("ta_offset", 0)
    else:
        x = low
    x = min(x, poly["ta_cap"]) - 0.6 * dmso
    dec = round_half_up(x * 10) / 10
    return round_half_up(dec), dec


def frac_tie(x: float) -> bool:
    return x - math.floor(x) == 0.5


def round_half_up(x: float) -> int:
    f = math.floor(x)
    return f + 1 if x - f >= 0.5 else f


# ----------------------------------------------------------------------------
# polymerase-tm and Biopython references
# ----------------------------------------------------------------------------
def ptm_modules():
    import polymerase_tm
    from polymerase_tm import core, constants
    return polymerase_tm, core, constants


def ptm_unrounded(core, constants, seq, pkey, buffer=None, salt=None, method=None):
    poly = constants.POLYMERASES[pkey]
    s = core._resolve_salt(pkey, buffer, salt)
    m = method or poly.get("tm_method", "owczarzy")
    if m == "schildkraut":
        return core.schildkraut_correction(core.calc_nn_raw(seq, poly["conc"] / 4), s)
    return core.owczarzy_correction(core.calc_nn_raw(seq, poly["conc"]), seq, s)


def ptm_ta_prerounded(constants, t1, t2, pkey, l1, l2, dmso):
    poly = constants.POLYMERASES[pkey]
    mn = min(t1, t2)
    if poly["ta_rule"] == "phusion":
        r = 0.93 * mn + 7.5
    elif min(l1, l2) >= poly.get("min_len", 8):
        r = mn + poly.get("ta_offset", 0)
    else:
        r = mn
    r = min(r, poly["ta_cap"])
    return r - 0.6 * dmso


def bio_tm(seq, conc_nM, saltcorr, Na=50.0, Mg=0.0, selfcomp=False):
    from Bio.SeqUtils import MeltingTemp as mt
    return mt.Tm_NN(seq, nn_table=mt.DNA_NN3, dnac1=conc_nM, dnac2=0, Na=Na, K=0, Tris=0, Mg=Mg, dNTPs=0,
                    saltcorr=saltcorr, selfcomp=selfcomp)


def o2008_regime(mono_mM, mg_mM):
    if mg_mM <= 0:
        return "mono-only"
    r = math.sqrt(mg_mM / 1000) / (mono_mM / 1000)
    return "R<0.22" if r < 0.22 else "0.22<=R<6" if r < 6 else "R>=6"


# ----------------------------------------------------------------------------
# bulk comparison
# ----------------------------------------------------------------------------
def cmd_bulk(args):
    import Bio
    polymerase_tm, core, constants = ptm_modules()
    rng = random.Random(SEED)
    seqs = make_sequences(args.n, rng)
    pkeys = sorted(constants.POLYMERASES)
    bkeys = sorted(constants.BUFFERS)
    print(f"{len(seqs)} sequences, {len(pkeys)} polymerases, {len(bkeys)} buffers", file=sys.stderr)

    # --- (a) polymerase-tm: Tm ---
    tm_cases, tm_ref = [], []
    for s in seqs:
        combos = [dict(polymerase=p) for p in pkeys]
        combos += [dict(polymerase=p, buffer=b) for p in ("q5", "phusion_hf", "taq") for b in bkeys]
        combos += [dict(polymerase=p, salt_mM=rng.randint(5, 1000)) for p in ("q5", "phusion_hf", "taq", "longamp")]
        combos += [dict(polymerase="q5", method="schildkraut"), dict(polymerase="phusion_hf", method="owczarzy")]
        for o in combos:
            p = o["polymerase"]
            ref_int = polymerase_tm.tm(s, p, buffer=o.get("buffer"), salt_mM=o.get("salt_mM"), method=o.get("method"))
            ref_f = ptm_unrounded(core, constants, s, p, o.get("buffer"), o.get("salt_mM"), o.get("method"))
            neb_f = neb_unrounded(core, constants, s, p, o.get("buffer"), o.get("salt_mM"), o.get("method"))
            tm_cases.append({"k": "tm", "s": s, "o": dict(o, compat="polymerase-tm")})
            tm_ref.append(("ptm", ref_int, ref_f))
            tm_cases.append({"k": "tm", "s": s, "o": o})
            tm_ref.append(("neb", round_half_up(neb_f), neb_f, ref_int))
    print(f"tm cases: {len(tm_cases)}", file=sys.stderr)
    tm_out = lua_eval(tm_cases)
    tm_mm, tm_maxdiff = [], {"ptm": 0.0, "neb": 0.0}
    neb_vs_ptm = {"same_salt": 0, "same_salt_differ": 0, "changed_salt": 0, "changed_salt_differ": 0}
    for c, ref, row in zip(tm_cases, tm_ref, tm_out):
        mode, ri, rf = ref[0], ref[1], ref[2]
        if row[0] == "ERR":
            tm_mm.append({"case": c, "err": row[1]}); continue
        tm_maxdiff[mode] = max(tm_maxdiff[mode], abs(row[0] - rf))
        if int(row[1]) != ri:
            tm_mm.append({"mode": mode, "case": c, "lua": int(row[1]), "ref": ri, "unrounded": rf, "tie": frac_tie(rf)})
        if mode == "neb":
            o = c["o"]
            changed = o.get("salt_mM") is None and (
                NEB_BUFFERS[o.get("buffer") or constants.POLYMERASES[o["polymerase"]]["buffer"]]
                != constants.BUFFERS[o.get("buffer") or constants.POLYMERASES[o["polymerase"]]["buffer"]])
            k = "changed_salt" if changed else "same_salt"
            neb_vs_ptm[k] += 1
            neb_vs_ptm[k + "_differ"] += int(row[1]) != ref[3]

    # --- (a) polymerase-tm: Ta ---
    dmsos = [0, 0, 0, 0, 1, 2.5, 3, 5, 7.5, 10]
    perm = list(range(len(seqs))); rng.shuffle(perm)
    ta_cases, ta_ref = [], []
    for i, s in enumerate(seqs):
        s2 = seqs[perm[i]]
        for j, p in enumerate(pkeys):
            o = {"polymerase": p}
            d = dmsos[(i + j) % len(dmsos)]
            if d:
                o["dmso"] = d
            if (i + j) % 7 == 0:
                o["buffer"] = bkeys[(i * 3 + j) % len(bkeys)]
            elif (i + j) % 7 == 1:
                o["salt_mM"] = rng.randint(5, 1000)
            r = polymerase_tm.ta(s, s2, p, dmso_pct=d, buffer=o.get("buffer"), salt_mM=o.get("salt_mM"))
            pre = ptm_ta_prerounded(constants, r[1], r[2], p, len(s), len(s2), d)
            ta_cases.append({"k": "ta", "s": s, "s2": s2, "o": dict(o, compat="polymerase-tm")})
            ta_ref.append(("ptm", r, pre))
            u1 = neb_unrounded(core, constants, s, p, o.get("buffer"), o.get("salt_mM"))
            u2 = neb_unrounded(core, constants, s2, p, o.get("buffer"), o.get("salt_mM"))
            nta, _ = neb_ta(constants, u1, u2, p, len(s), len(s2), d)
            ta_cases.append({"k": "ta", "s": s, "s2": s2, "o": o})
            ta_ref.append(("neb", (nta, round_half_up(u1), round_half_up(u2)), r))
    print(f"ta cases: {len(ta_cases)}", file=sys.stderr)
    ta_out = lua_eval(ta_cases)
    ta_mm, ta_ties, ta_neb_vs_ptm = [], 0, 0
    for c, (mode, r, extra), row in zip(ta_cases, ta_ref, ta_out):
        if mode == "ptm":
            ta_ties += frac_tie(extra)
        else:
            ta_neb_vs_ptm += tuple(r) != tuple(extra)
        if row[0] == "ERR" or (int(row[0]), int(row[1]), int(row[2])) != tuple(r):
            ta_mm.append({"mode": mode, "case": c, "lua": row, "ref": list(r)})

    # --- Owczarzy 2008 vs polymerase-tm calc_sdm_tm and Biopython saltcorr=7 ---
    grid = [(m, g) for m in (10, 20, 50, 150) for g in (0, 0.5, 1.5, 2, 6, 20, 100)]
    sub = seqs[: min(len(seqs), 1500)]
    o8_cases, o8_ref = [], []
    for s in sub:
        for mono, mg in grid:
            o8_cases.append({"k": "tm", "s": s, "o": {"method": "owczarzy2008", "salt_mM": mono, "mg_mM": mg}})
            o8_ref.append((core.calc_sdm_tm(s, mono, mg, 500), bio_tm(s, 500, 7, Na=mono, Mg=mg), mono, mg,
                           s == revcomp(s)))
    o8_out = lua_eval(o8_cases)
    o8 = {}
    for (ptm_v, bio_v, mono, mg, pal), row in zip(o8_ref, o8_out):
        reg = o2008_regime(mono, mg)
        e = o8.setdefault(reg, {"n": 0, "max_vs_ptm": 0.0, "max_vs_bio_nonpal": 0.0, "n_pal": 0, "int_mismatch_vs_ptm": 0})
        e["n"] += 1
        e["max_vs_ptm"] = max(e["max_vs_ptm"], abs(row[0] - ptm_v))
        e["int_mismatch_vs_ptm"] += int(row[1]) != round_half_up(ptm_v)
        if pal:
            e["n_pal"] += 1
        else:
            e["max_vs_bio_nonpal"] = max(e["max_vs_bio_nonpal"], abs(row[0] - bio_v))

    # --- (b) Biopython: raw NN and salt corrections ---
    bio_cases, bio_ref = [], []
    for s in seqs:
        pal = s == revcomp(s)
        for conc in (500, 400, 200, 125):
            bio_cases.append({"k": "raw", "s": s, "c": conc})
            bio_ref.append(("raw", conc, None, bio_tm(s, conc, 0), pal))
            if pal:
                bio_cases.append({"k": "raw", "s": s, "c": conc})
                bio_ref.append(("raw selfcomp=True", conc, None, bio_tm(s, conc, 0, selfcomp=True), pal))
        for p in ("q5", "q5u_hot_start", "taq", "hemo_klentaq", "longamp", "multiplex", "vent"):
            poly = constants.POLYMERASES[p]
            salt = constants.BUFFERS[poly["buffer"]]
            bio_cases.append({"k": "tm", "s": s, "o": {"polymerase": p, "salt_mM": salt}})
            bio_ref.append(("owczarzy2004", poly["conc"], salt, bio_tm(s, poly["conc"], 6, Na=salt), pal))
        bio_cases.append({"k": "tm", "s": s, "o": {"polymerase": "phusion_hf"}})
        bio_ref.append(("schildkraut", 125, 222, bio_tm(s, 125, 1, Na=222), pal))
    bio_out = lua_eval(bio_cases)
    bio = {}
    for (kind, conc, salt, ref, pal), row in zip(bio_ref, bio_out):
        v = row[2] if kind.startswith("raw") else row[0]
        key = kind + (" (self-complementary)" if pal and "selfcomp" not in kind else "")
        e = bio.setdefault(key, {"n": 0, "max_abs_diff": 0.0, "mean_diff": 0.0})
        e["n"] += 1
        d = v - ref
        e["max_abs_diff"] = max(e["max_abs_diff"], abs(d))
        e["mean_diff"] += d
    for e in bio.values():
        e["mean_diff"] /= e["n"]

    summary = {
        "versions": {"polymerase_tm": polymerase_tm.__version__, "biopython": Bio.__version__, "seed": SEED},
        "n_sequences": len(seqs), "n_palindromes": sum(s == revcomp(s) for s in seqs),
        "len_range": [min(map(len, seqs)), max(map(len, seqs))],
        "polymerase_tm": {
            "tm_cases_per_mode": len(tm_cases) // 2, "ta_cases_per_mode": len(ta_cases) // 2,
            "tm_int_mismatches_compat": sum(m.get("mode") == "ptm" for m in tm_mm),
            "tm_int_mismatches_neb_vs_python_neb_reference": sum(m.get("mode") == "neb" for m in tm_mm),
            "tm_max_unrounded_diff": tm_maxdiff,
            "ta_mismatches_compat": sum(m.get("mode") == "ptm" for m in ta_mm),
            "ta_mismatches_neb_vs_python_neb_reference": sum(m.get("mode") == "neb" for m in ta_mm),
            "ta_compat_cases_with_exact_half_before_rounding": ta_ties,
            "neb_mode_tm_differs_from_polymerase_tm": neb_vs_ptm,
            "neb_mode_ta_differs_from_polymerase_tm": ta_neb_vs_ptm,
        },
        "owczarzy2008": {"cases": len(o8_cases), "by_regime": o8},
        "biopython": bio,
    }
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    (out / "bulk_summary.json").write_text(json.dumps(summary, indent=1))
    (out / "bulk_mismatches.json").write_text(json.dumps({"tm": tm_mm, "ta": ta_mm}, indent=0))
    print(json.dumps(summary, indent=1))


# ----------------------------------------------------------------------------
# NEB Tm API
# ----------------------------------------------------------------------------
NEB = "https://tmapi.neb.com"


def neb_session():
    import requests
    return requests.Session()  # library-default headers; no User-Agent override


def neb_post(sess, path, body, log, pause=1.0):
    time.sleep(pause)
    r = sess.post(NEB + path, json=body, timeout=60)
    entry = {"path": path, "body": {k: v for k, v in body.items() if k != "email"}, "status": r.status_code}
    try:
        entry["response"] = r.json()
    except ValueError:
        entry["response_text"] = r.text[:500]
    log.append(entry)
    if r.status_code in (401, 403, 429):
        raise SystemExit(f"NEB refused scripted access: HTTP {r.status_code} on {path}")
    return entry


def neb_pairs(rng, n, lmin=15, lmax=35):
    out = []
    for _ in range(n):
        a = rand_seq(rng, rng.randint(lmin, lmax), rng.uniform(0.25, 0.75))
        b = rand_seq(rng, rng.randint(lmin, lmax), rng.uniform(0.25, 0.75))
        out.append([a, b])
    return out


def cmd_neb(args):
    import requests  # noqa: F401
    polymerase_tm, core, constants = ptm_modules()
    email = None
    if args.email_file:
        email = Path(args.email_file).read_text().strip() or None
    rng = random.Random(SEED + 1)
    sess = neb_session()
    log: list[dict] = []

    time.sleep(0.5)
    r = sess.get(NEB + "/docs/productcodes", timeout=60)
    log.append({"path": "/docs/productcodes", "status": r.status_code,
                "response": r.json() if r.ok else r.text[:500]})
    if r.status_code in (401, 403, 429):
        raise SystemExit(f"NEB refused: HTTP {r.status_code}")
    products = r.json()["data"]

    # Lua polymerase table -> NEB product codes
    lua = lua_polymerases()
    mapped = {p["neb_code"]: p for p in lua if p.get("neb_code")}

    # primer pairs: realistic 15-35 nt, plus near-half Tm probes and short primers
    base = neb_pairs(rng, 300)
    real = list(REAL_PRIMERS.values())
    base += [[real[i], real[(i + 3) % len(real)]] for i in range(len(real))]
    # rounding probes: sequences whose unrounded Q5 Tm is within 0.08 of x.5
    probes = []
    while len(probes) < 60:
        s = rand_seq(rng, rng.randint(14, 40), rng.uniform(0.2, 0.8))
        f = ptm_unrounded(core, constants, s, "q5") % 1
        if abs(f - 0.5) < 0.08:
            probes.append(s)
    probe_pairs = [[probes[i], probes[i + 1]] for i in range(0, len(probes), 2)]
    short = neb_pairs(rng, 30, 8, 24)  # min_len behaviour (Vent: 21 nt; 8-nt primers)

    def batch(prodcode, pairs, conc, monosalt=None, chunk=100):
        res = []
        for i in range(0, len(pairs), chunk):
            body = {"prodcode": prodcode, "conc": conc, "seqpairs": pairs[i:i + chunk]}
            if monosalt is not None:
                body["monosalt"] = monosalt
            if email:
                body["email"] = email
            res.append(neb_post(sess, "/tm/batch", body, log))
        return res

    # mapped products
    for code, p in sorted(mapped.items()):
        conc = p["conc_nM"] / 1000
        if code in ("q5-0", "phusion-0", "taq-0"):
            pairs = base + probe_pairs + short
        else:
            pairs = base[:40] + probe_pairs[:10] + short[:10]
        print(code, len(pairs), file=sys.stderr)
        batch(code, pairs, conc)
    # primer-concentration and custom-salt overrides
    batch("q5-0", base[:50], 0.25)
    batch("phusion-0", base[:50], 1.0)
    batch("custom", base[:50], 0.5, monosalt=75)
    batch("custom", base[:50], 0.2, monosalt=40)
    batch("default", base[:50], 0.5)
    # NEB products that polymerase-tm lacks: 30 pairs each to infer their buffer/rules
    unmapped = [x["prodcode"] for x in products if x["prodcode"] not in mapped and x["prodcode"] not in ("custom", "default")]
    for code in unmapped:
        fam_conc = 0.5 if code.startswith(("q5", "phusion")) else 0.4 if code.startswith(("lataq", "lahstaq")) else 0.2
        batch(code, base[:20] + short[:10], fam_conc)
    # long format for notes (one per mapped family + custom/default)
    for code in ("q5-0", "phusion-0", "taq-0", "vent-0", "lataq-0", "custom", "default"):
        body = {"prodcode": code, "conc": 0.5, "seq1": REAL_PRIMERS["M13R(-48)"], "seq2": REAL_PRIMERS["M13F(-20)"], "fmt": "long"}
        if code == "custom":
            body["monosalt"] = 75
        if email:
            body["email"] = email
        neb_post(sess, "/tm", body, log)

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    (out / "neb_raw.json").write_text(json.dumps({"products": products, "log": log}, indent=0))
    print(f"{len(log)} requests logged", file=sys.stderr)
    neb_phase2(sess, email, out / "neb_raw.json")
    neb_compare(out / "neb_raw.json", out / "neb_summary.json")


def neb_phase2(sess, email, raw_path):
    """Long-format calls (unrounded Tm, one-decimal Ta) for every mapped product,
    and high-Tm pairs to locate Ta caps. Appends to neb_raw.json."""
    raw = json.loads(Path(raw_path).read_text())
    log = raw["log"]
    rng = random.Random(SEED + 3)
    lua = lua_polymerases()
    mapped = {p["neb_code"]: p for p in lua if p.get("neb_code")}
    # long format: M13 pair + one random pair per mapped product
    pairs = [[REAL_PRIMERS["M13R(-48)"], REAL_PRIMERS["M13F(-20)"]]]
    for code, p in sorted(mapped.items()):
        for s1, s2 in pairs + neb_pairs(rng, 1):
            body = {"prodcode": code, "conc": p["conc_nM"] / 1000, "seq1": s1, "seq2": s2, "fmt": "long"}
            if email:
                body["email"] = email
            neb_post(sess, "/tm", body, log)
    # cap probes: GC-rich pairs, lower Tm spread over ~60-95 °C
    probes = []
    for i in range(20):
        L = 16 + i * 2
        probes.append([rand_seq(rng, L, 0.55 + 0.02 * i), rand_seq(rng, L + 2, 0.6 + 0.02 * i)])
    cap_codes = sorted(set(mapped) | {"taqmaster-0", "taqquickload-0", "ctaq-1", "q5xt-0", "taq-2"})
    for code in cap_codes:
        conc = mapped[code]["conc_nM"] / 1000 if code in mapped else (0.5 if code.startswith("q5") else 0.2)
        body = {"prodcode": code, "conc": conc, "seqpairs": probes}
        if email:
            body["email"] = email
        neb_post(sess, "/tm/batch", body, log)
    Path(raw_path).write_text(json.dumps(raw, indent=0))
    print(f"{len(log)} requests logged (after phase 2)", file=sys.stderr)


def lua_polymerases():
    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "p.lua"
        outp = Path(td) / "p.json"
        script.write_text(r"""
package.path = arg[1] .. "/lua/?.lua;" .. package.path
local M = require("plasmid.bio.tm")
local fh = assert(io.open(arg[2], "wb")); fh:write(vim.json.encode({ polymerases = M.POLYMERASES, buffers = M.BUFFERS })); fh:close()
""")
        subprocess.run([NVIM, "-l", str(script), str(ROOT), str(outp)], check=True)
        d = json.loads(outp.read_text())
    return d["polymerases"]


def neb_rows(raw):
    """(prodcode, request body, response row) for every successful batch row."""
    out = []
    for e in raw["log"]:
        if e["path"] == "/tm/batch" and e["status"] == 200 and e["response"].get("success"):
            for row in e["response"]["data"]:
                out.append((e["body"]["prodcode"], e["body"], row))
    return out


def neb_opts(code, body, mapped):
    if code in mapped:
        return {"polymerase": mapped[code]["key"], "conc_nM": body["conc"] * 1000}
    if code == "custom":
        return {"polymerase": "q5", "salt_mM": body["monosalt"], "conc_nM": body["conc"] * 1000}
    if code == "default":
        return {"polymerase": "q5", "salt_mM": 50, "conc_nM": body["conc"] * 1000}
    return None


def neb_compare(raw_path, summary_path):
    """Compare every NEB result with the Lua module (NEB mode)."""
    raw = json.loads(Path(raw_path).read_text())
    lua = lua_polymerases()
    mapped = {p["neb_code"]: p for p in lua if p.get("neb_code")}
    rows = [(c, b, r, neb_opts(c, b, mapped)) for c, b, r in neb_rows(raw)]
    rows = [x for x in rows if x[3] is not None]
    cases = []
    for code, b, r, o in rows:
        cases.append({"k": "ta", "s": r["seq1"].upper(), "s2": r["seq2"].upper(), "o": o})
    out = lua_eval(cases)
    stats, mism, ruleless = {}, [], {}
    for (code, b, r, o), row in zip(rows, out):
        key = f"{code}@{b['conc']}" + (f"/salt{b['monosalt']}" if "monosalt" in b else "")
        st = stats.setdefault(key, {"pairs": 0, "tm_compared": 0, "tm_mismatch": 0, "ta_compared": 0, "ta_mismatch": 0})
        st["pairs"] += 1
        st["tm_compared"] += 2
        bad_tm = (int(row[1]), int(row[2])) != (r["tm1"], r["tm2"])
        st["tm_mismatch"] += bad_tm
        bad_ta = False
        if code in NEB_RULELESS:
            low = min(r["tm1"], r["tm2"])
            e = ruleless.setdefault(code, {"pairs": 0, "neb_ta_equals_lower_tm": 0, "max_neb_ta": None})
            e["pairs"] += 1
            e["neb_ta_equals_lower_tm"] += abs(r["ta"] - low) <= 1  # lower Tm rounded via 0.1 step
            e["max_neb_ta"] = max(e["max_neb_ta"] or 0, r["ta"])
        elif code in mapped:
            st["ta_compared"] += 1
            bad_ta = int(row[0]) != r["ta"]
            st["ta_mismatch"] += bad_ta
        if bad_tm or bad_ta:
            mism.append({"product": key, "neb": r, "lua": row[:3]})
    # long format: unrounded Tm and 0.1 °C Ta
    long_cases, long_meta = [], []
    for e in raw["log"]:
        if e["path"] != "/tm" or e["status"] != 200:
            continue
        b, d = e["body"], e["response"]["data"]
        o = neb_opts(b["prodcode"], b, mapped)
        if o is None:
            continue
        # NEB's long format echoes the concentration it used; for Phusion this is
        # conc/4 and the endpoint then divides by 4 AGAIN (batch divides once).
        # Evaluate at the echoed concentration so the double division is explicit.
        o = dict(o, conc_nM=d["conc"] * 1000)
        p1 = d["p1"][0] if isinstance(d["p1"], list) else d["p1"]
        p2 = d["p2"][0] if isinstance(d["p2"], list) else d["p2"]
        long_cases += [{"k": "tm", "s": p1["seq"].upper(), "o": o}, {"k": "tm", "s": p2["seq"].upper(), "o": o},
                       {"k": "ta", "s": p1["seq"].upper(), "s2": p2["seq"].upper(), "o": o}]
        long_meta.append((b["prodcode"], p1["tm"], p2["tm"], d.get("ta"), d.get("conc")))
    lo = lua_eval(long_cases) if long_cases else []
    long_stats = {"note": "evaluated at the conc NEB echoes (Phusion: 0.125 uM for 0.5 sent, then /4 again)",
                  "calls": len(long_meta), "tm_values": 0, "tm_max_abs_diff": 0.0, "ta_compared": 0,
                  "ta_decimal_mismatch": 0, "ta_skipped_ruleless_or_custom": 0, "examples": []}
    for i, (code, n1, n2, nta, nconc) in enumerate(long_meta):
        t1, t2, ta = lo[3 * i], lo[3 * i + 1], lo[3 * i + 2]
        long_stats["tm_values"] += 2
        long_stats["tm_max_abs_diff"] = max(long_stats["tm_max_abs_diff"], abs(t1[0] - n1), abs(t2[0] - n2))
        if code in mapped and code not in NEB_RULELESS:
            long_stats["ta_compared"] += 1
            long_stats["ta_decimal_mismatch"] += abs(ta[4] - nta) > 1e-9
        else:
            long_stats["ta_skipped_ruleless_or_custom"] += 1
        if len(long_stats["examples"]) < 8:
            long_stats["examples"].append({"prodcode": code, "neb_tm": [n1, n2], "lua_tm": [t1[0], t2[0]],
                                           "neb_ta": nta, "lua_ta_decimal": ta[4], "neb_conc_uM": nconc})
    inferred = infer_unmapped(raw, mapped)
    summary = {"neb_version": next((e["response"]["meta"]["version"] for e in raw["log"]
                                    if isinstance(e.get("response"), dict) and e["response"].get("meta")), None),
               "requests": len(raw["log"]),
               "http_status_counts": _count(e["status"] for e in raw["log"]),
               "rows_compared": len(rows),
               "tm_values_compared": sum(v["tm_compared"] for v in stats.values()),
               "tm_mismatches": sum(v["tm_mismatch"] for v in stats.values()),
               "ta_values_compared": sum(v["ta_compared"] for v in stats.values()),
               "ta_mismatches": sum(v["ta_mismatch"] for v in stats.values()),
               "by_product": stats, "ruleless_products": ruleless, "long_format": long_stats,
               "mismatches": mism, "unmapped_inferred": inferred}
    Path(summary_path).write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k not in ("mismatches", "unmapped_inferred", "by_product")}, indent=1))


def _count(it):
    d = {}
    for x in it:
        d[str(x)] = d.get(str(x), 0) + 1
    return d


def infer_unmapped(raw, mapped):
    """NEB products without a plugin entry: which plugin polymerase + buffer reproduces NEB's Tm1/Tm2/Ta?"""
    lua = lua_polymerases()
    buffers = sorted(NEB_BUFFERS)
    rows = {}
    for code, b, r in neb_rows(raw):
        if code in mapped or code in ("custom", "default"):
            continue
        rows.setdefault(code, []).append((b["conc"], r))
    result = {}
    for code, rr in sorted(rows.items()):
        combos = [(p["key"], bf) for p in lua for bf in buffers]
        cases = []
        for pk, bf in combos:
            for conc, r in rr:
                cases.append({"k": "ta", "s": r["seq1"].upper(), "s2": r["seq2"].upper(),
                              "o": {"polymerase": pk, "buffer": bf, "conc_nM": conc * 1000}})
        out = lua_eval(cases)
        tm_hits, full_hits, i = set(), [], 0
        for pk, bf in combos:
            tm_ok = ta_ok = True
            for conc, r in rr:
                ta, t1, t2 = int(out[i][0]), int(out[i][1]), int(out[i][2])
                i += 1
                tm_ok &= (t1, t2) == (r["tm1"], r["tm2"])
                ta_ok &= ta == r.get("ta")
            if tm_ok:
                tm_hits.add(bf)
                if ta_ok:
                    full_hits.append(f"{pk}+{bf}")
        name = next((x["name"] for x in raw["products"] if x["prodcode"] == code), code)
        result[code] = {"name": name, "pairs": len(rr), "conc_uM_sent": rr[0][0],
                        "salt_mM": sorted({NEB_BUFFERS[bf] for bf in tm_hits}),
                        "plugin_polymerase+buffer_reproducing_tm_and_ta": full_hits,
                        "ruleless_ta": code in NEB_RULELESS}
    return result


# ----------------------------------------------------------------------------
# compact fixture tests/oracle/expected/tm.json
# ----------------------------------------------------------------------------
def cmd_fixtures(args):
    import Bio
    polymerase_tm, core, constants = ptm_modules()
    rng = random.Random(SEED + 2)
    cases = []
    seqs = make_sequences(400, rng)
    real = list(REAL_PRIMERS.values()) + load_ptm_test_seqs()
    pick = real + rng.sample(seqs[len(real):], 120)
    pkeys = sorted(constants.POLYMERASES)
    bkeys = sorted(constants.BUFFERS)
    # polymerase-tm (compat mode): every polymerase and every buffer appears
    for i, s in enumerate(pick):
        o = {"polymerase": pkeys[i % len(pkeys)], "compat": "polymerase-tm"}
        if i % 4 == 1:
            o["buffer"] = bkeys[i % len(bkeys)]
        elif i % 4 == 2:
            o["salt_mM"] = rng.choice([10, 25, 50, 75, 100, 200, 300, 500])
        rf = ptm_unrounded(core, constants, s, o["polymerase"], o.get("buffer"), o.get("salt_mM"))
        cases.append({"source": "polymerase-tm", "kind": "tm", "seq": s, "opts": o,
                      "expect": {"tm_rounded": polymerase_tm.tm(s, o["polymerase"], buffer=o.get("buffer"),
                                                                salt_mM=o.get("salt_mM")),
                                 "tm": rf}})
    n_tie = 0
    for i in range(130):
        s1, s2 = rng.sample(pick, 2)
        p = pkeys[i % len(pkeys)]
        d = [0, 0, 1, 2.5, 3, 5, 7.5, 10][i % 8]
        r = polymerase_tm.ta(s1, s2, p, dmso_pct=d)
        pre = ptm_ta_prerounded(constants, r[1], r[2], p, len(s1), len(s2), d)
        n_tie += frac_tie(pre)
        o = {"polymerase": p, "compat": "polymerase-tm"}
        if d:
            o["dmso"] = d
        cases.append({"source": "polymerase-tm", "kind": "ta", "seq": s1, "seq2": s2, "opts": o,
                      "expect": {"ta": r[0], "tm1": r[1], "tm2": r[2]}})
    # rounding: DMSO ties. NEB mode rounds half up (expected from the independent
    # Python NEB procedure neb_ta); polymerase-tm/Python round() goes to even.
    n_tie = 0
    for s1 in pick:
        for s2 in pick[:15]:
            for d in (2.5, 7.5):
                u1, u2 = (neb_unrounded(core, constants, x, "q5") for x in (s1, s2))
                pre = min(u1, u2) + 1
                x = min(pre, 72) - 0.6 * d
                if pre > 72 and frac_tie(x) and n_tie < 12:
                    n_tie += 1
                    nta, dec = neb_ta(constants, u1, u2, "q5", len(s1), len(s2), d)
                    r = polymerase_tm.ta(s1, s2, "q5", dmso_pct=d)
                    cases.append({"source": "rounding", "kind": "ta", "seq": s1, "seq2": s2,
                                  "opts": {"polymerase": "q5", "dmso": d},
                                  "expect": {"ta": nta, "tm1": round_half_up(u1), "tm2": round_half_up(u2)},
                                  "note": f"capped 72 - {0.6 * d:.1f} = {x}; half up -> {nta}; "
                                          f"polymerase-tm (Python round) gives {r[0]}"})
    # NEB-mode values from the independent Python NEB procedure, all polymerases
    for i in range(66):
        s1, s2 = rng.sample(pick, 2)
        p = pkeys[i % len(pkeys)]
        u1, u2 = (neb_unrounded(core, constants, x, p) for x in (s1, s2))
        nta, dec = neb_ta(constants, u1, u2, p, len(s1), len(s2))
        cases.append({"source": "neb-procedure", "kind": "ta", "seq": s1, "seq2": s2, "opts": {"polymerase": p},
                      "expect": {"ta": nta, "tm1": round_half_up(u1), "tm2": round_half_up(u2), "ta_decimal": dec,
                                 "tm1_unrounded": u1, "tm2_unrounded": u2}})
    # Biopython: raw (1 M) and corrected, unrounded, non-self-complementary
    bio_pick = [s for s in pick if s != revcomp(s)][:60]
    for i, s in enumerate(bio_pick):
        conc = [500, 200, 400][i % 3]
        cases.append({"source": "biopython", "kind": "raw", "seq": s, "conc_nM": conc,
                      "expect": {"raw": bio_tm(s, conc, 0)}})
        salt = [150, 51, 70, 100, 90, 40, 180][i % 7]
        cases.append({"source": "biopython", "kind": "tm", "seq": s,
                      "opts": {"polymerase": "q5", "salt_mM": salt, "conc_nM": conc},
                      "expect": {"tm": bio_tm(s, conc, 6, Na=salt)}, "saltcorr": 6})
        if i % 2 == 0:
            cases.append({"source": "biopython", "kind": "tm", "seq": s,
                          "opts": {"polymerase": "phusion_hf", "conc_nM": conc},
                          "expect": {"tm": bio_tm(s, conc / 4, 1, Na=222)}, "saltcorr": 1})
        mono, mg = [(50, 2), (150, 1.5), (20, 6), (10, 100), (150, 0.5), (50, 0)][i % 6]
        cases.append({"source": "biopython", "kind": "tm", "seq": s,
                      "opts": {"method": "owczarzy2008", "salt_mM": mono, "mg_mM": mg, "conc_nM": conc},
                      "expect": {"tm": bio_tm(s, conc, 7, Na=mono, Mg=mg)}, "saltcorr": 7,
                      "regime": o2008_regime(mono, mg)})
    # NEB values (only ever copied from real API responses)
    neb_cases = []
    if args.neb:
        raw = json.loads(Path(args.neb).read_text())
        lua = lua_polymerases()
        mapped = {p["neb_code"]: p for p in lua if p.get("neb_code")}
        per = {}
        for e in raw["log"]:
            if e["path"] != "/tm/batch" or e["status"] != 200 or not e["response"].get("success"):
                continue
            b = e["body"]; code = b["prodcode"]
            if code in mapped:
                o = {"polymerase": mapped[code]["key"], "conc_nM": b["conc"] * 1000}
            elif code in ("custom", "default"):
                o = {"polymerase": "q5", "salt_mM": b.get("monosalt", 50), "conc_nM": b["conc"] * 1000}
            else:
                continue
            tag = f"{code}@{b['conc']}/{b.get('monosalt')}"
            for row in e["response"]["data"]:
                # rows where NEB's Ta differs from the rule applied to the integer Tms
                # (two-step rounding evidence) are always kept, up to 6 per tag
                interesting = False
                if code in mapped and code not in NEB_RULELESS:
                    pk = mapped[code]["key"]
                    pre = ptm_ta_prerounded(constants, row["tm1"], row["tm2"], pk, len(row["seq1"]), len(row["seq2"]), 0)
                    interesting = round(pre) != row["ta"]
                if interesting and per.get(tag + "!", 0) < 6:
                    per[tag + "!"] = per.get(tag + "!", 0) + 1
                elif per.get(tag, 0) < (12 if code in ("q5-0", "phusion-0", "taq-0") else 6):
                    per[tag] = per.get(tag, 0) + 1
                else:
                    continue
                exp = {"tm1": row["tm1"], "tm2": row["tm2"]}
                if code in mapped and code not in NEB_RULELESS:
                    exp["ta"] = row["ta"]
                neb_cases.append({"source": "neb", "kind": "ta", "seq": row["seq1"].upper(),
                                  "seq2": row["seq2"].upper(), "opts": o, "neb_prodcode": code,
                                  "expect": exp, "neb_version": e["response"]["meta"].get("version")})
        for e in raw["log"]:
            if e["path"] != "/tm" or e["status"] != 200:
                continue
            b, d = e["body"], e["response"]["data"]
            if b["prodcode"] in mapped:
                o = {"polymerase": mapped[b["prodcode"]]["key"], "conc_nM": d["conc"] * 1000}
            elif b["prodcode"] in ("custom", "default"):
                o = {"polymerase": "q5", "salt_mM": b.get("monosalt", 50), "conc_nM": d["conc"] * 1000}
            else:
                continue
            p1 = d["p1"][0] if isinstance(d["p1"], list) else d["p1"]
            p2 = d["p2"][0] if isinstance(d["p2"], list) else d["p2"]
            exp = {"tm1_unrounded": p1["tm"], "tm2_unrounded": p2["tm"]}
            if b["prodcode"] in mapped and b["prodcode"] not in NEB_RULELESS and not b["prodcode"].startswith("phusion"):
                exp["ta_decimal"] = d["ta"]
            neb_cases.append({"source": "neb-long", "kind": "ta", "seq": p1["seq"].upper(), "seq2": p2["seq"].upper(),
                              "opts": o, "neb_prodcode": b["prodcode"], "expect": exp,
                              "note": "conc_nM = concentration echoed by NEB (Phusion: already /4; endpoint divides again)"})
    elif FIXTURE.exists():
        neb_cases = [c for c in json.loads(FIXTURE.read_text())["cases"] if c["source"] in ("neb", "neb-long")]
    cases += neb_cases
    doc = {"generated_by": "tests/tools/tm_oracle.py fixtures",
           "versions": {"polymerase_tm": polymerase_tm.__version__, "biopython": Bio.__version__},
           "tolerance_unrounded_C": 0.01,
           "counts": _count(c["source"] for c in cases),
           "cases": cases}
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE.write_text(json.dumps(doc, separators=(",", ":")) + "\n")
    print(json.dumps(doc["counts"]), file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bulk"); b.add_argument("--n", type=int, default=5000); b.add_argument("--out", default="tm_oracle_out")
    n = sub.add_parser("neb"); n.add_argument("--email-file"); n.add_argument("--no-email", action="store_true")
    n.add_argument("--out", default="tm_oracle_out")
    c = sub.add_parser("compare-neb"); c.add_argument("raw"); c.add_argument("--out", default="tm_oracle_out")
    n2 = sub.add_parser("neb-phase2"); n2.add_argument("raw"); n2.add_argument("--email-file")
    f = sub.add_parser("fixtures"); f.add_argument("--neb")
    a = ap.parse_args()
    if a.cmd == "bulk":
        cmd_bulk(a)
    elif a.cmd == "neb":
        cmd_neb(a)
    elif a.cmd == "compare-neb":
        neb_compare(a.raw, Path(a.out) / "neb_summary.json")
    elif a.cmd == "neb-phase2":
        email = Path(a.email_file).read_text().strip() if a.email_file else None
        neb_phase2(neb_session(), email or None, a.raw)
    else:
        cmd_fixtures(a)


if __name__ == "__main__":
    main()
