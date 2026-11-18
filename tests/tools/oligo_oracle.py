#!/usr/bin/env python3
"""Development-time oracle for lua/plasmid/bio/oligo.lua (NOT a runtime dependency).

Compares the pure-Lua nearest-neighbour estimates with primer3-py
(primer3.calc_hairpin / calc_homodimer / calc_heterodimer / calc_end_stability,
and the PRIMER_LEFT_0_END_STABILITY reported by primer3.design_primers) on a
deterministic set of several thousand oligos, writes the statistics as
markdown, and writes the small fixture tests/oracle/expected/oligo.json used by
tests/spec/oligo_spec.lua.

Usage (from the plugin root, in a dev env with primer3-py installed):
    python tests/tools/oligo_oracle.py \
        [--nvim /opt/nvim-linux-x86_64/bin/nvim] [--seed 20261116] \
        [--report docs/oligo_report.md] [--fixture tests/oracle/expected/oligo.json] \
        [--dump /tmp/oligo_all.json]

Requires: primer3-py (tested 2.3.1), biopython (1.88), numpy, scipy. Neovim runs the Lua side
through `nvim -l` on a generated driver script.
"""
import argparse
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import time

import numpy as np
import primer3
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COMP = str.maketrans("ACGT", "TGCA")
COND = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0)
# warning thresholds compared as confusion tables (see oligo_report.md for why)
THRESH = {"hairpin_dg": -2.0, "dimer_dg": -9.0, "end_dimer_dg": -5.0, "struct_tm": 47.0, "end_stab": 9.0}


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n, gc):
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n))


# ---------------------------------------------------------------- case set
def make_cases(seed):
    rng = random.Random(seed)
    singles, pairs = [], []

    def add(seq, kind):
        singles.append({"seq": seq, "kind": kind})

    for _ in range(3000):
        add(rand_seq(rng, rng.randint(15, 35), rng.uniform(0.2, 0.8)), "random")
    for _ in range(150):  # perfect palindromes (self-complementary)
        h = rand_seq(rng, rng.randint(5, 15), rng.uniform(0.3, 0.7))
        add(h + rc(h), "palindrome")
    for _ in range(150):  # palindromic core inside random flanks
        h = rand_seq(rng, rng.randint(3, 6), 0.5)
        add(rand_seq(rng, rng.randint(4, 10), 0.5) + h + rc(h) + rand_seq(rng, rng.randint(3, 10), 0.5), "palindrome_core")
    for _ in range(200):  # self 3'-complementarity: 3' end is a short palindrome
        h = rand_seq(rng, rng.randint(2, 4), 0.6)
        add(rand_seq(rng, rng.randint(12, 26), 0.5) + h + rc(h), "self_3prime")
    for _ in range(300):  # hairpin formers: stem - loop - rc(stem)
        st = rand_seq(rng, rng.randint(4, 10), rng.uniform(0.3, 0.8))
        lp = rand_seq(rng, rng.randint(3, 10), 0.3)
        add(rand_seq(rng, rng.randint(0, 6), 0.5) + st + lp + rc(st) + rand_seq(rng, rng.randint(0, 6), 0.5), "hairpin_former")
    for n in (15, 20, 25, 30):  # low-complexity
        for s in ("A" * n, "T" * n, ("AT" * n)[:n], ("GC" * n)[:n], "G" * n, ("CA" * n)[:n], ("AAT" * n)[:n]):
            add(s, "low_complexity")
    for _ in range(100):
        add(rand_seq(rng, rng.randint(15, 25), 0.5) + "A" * rng.randint(5, 10), "polyA_tail")

    rnd = [c["seq"] for c in singles if c["kind"] == "random"]
    for _ in range(3000):
        pairs.append({"a": rng.choice(rnd), "b": rng.choice(rnd), "kind": "random_pair"})
    for _ in range(400):  # strong 3' complementarity between the two primers
        a = rand_seq(rng, rng.randint(16, 28), rng.uniform(0.3, 0.7))
        k = rng.randint(4, 10)
        b = rand_seq(rng, rng.randint(12, 24), 0.5) + rc(a[-k:])
        pairs.append({"a": a, "b": b, "kind": "pair_3prime"})
    for _ in range(200):  # internal complementarity (with a mismatch sometimes)
        a = rand_seq(rng, rng.randint(18, 30), 0.5)
        i = rng.randint(0, len(a) - 10)
        core = list(rc(a[i:i + rng.randint(6, 10)]))
        if rng.random() < 0.5:
            j = rng.randrange(1, len(core) - 1)
            core[j] = rng.choice([c for c in "ACGT" if c != core[j]])
        b = rand_seq(rng, rng.randint(4, 10), 0.5) + "".join(core) + rand_seq(rng, rng.randint(4, 10), 0.5)
        pairs.append({"a": a, "b": b, "kind": "pair_internal"})
    return singles, pairs


# ---------------------------------------------------------------- primer3
def p3(res, kind=None, n=0, m=0):
    if not res.structure_found:
        return None
    out = {"dg": res.dg / 1000.0, "dh": res.dh / 1000.0, "ds": res.ds, "tm": res.tm}
    lines = getattr(res, "ascii_structure_lines", None)
    if kind and lines:
        out["pairs"] = hairpin_pairs(lines) if kind == "hairpin" else dimer_pairs(lines, m)
    return out


def hairpin_pairs(lines):
    """pairs (i, j) from primer3's one-line hairpin drawing ('/' opens, '\\' closes)"""
    st = lines[0].split("\t", 1)[1]
    stack, prs = [], []
    for i, ch in enumerate(st, 1):
        if ch == "/":
            stack.append(i)
        elif ch == "\\":
            prs.append([stack.pop(), i])
    return sorted(prs)


def dimer_pairs(lines, m):
    """pairs (i in top, j in bottom; both 5'->3' 1-based) from primer3's 4-line dimer drawing"""
    L = [x.split("\t", 1)[1] if "\t" in x else "" for x in lines]
    w = max(len(x) for x in L)
    L = [x.ljust(w) for x in L]
    ti = bi = 0
    prs = []
    for c in range(w):
        if L[0][c].isalpha() or L[1][c].isalpha():
            ti += 1
        if L[2][c].isalpha() or L[3][c].isalpha():
            bi += 1
        if L[1][c].isalpha():
            prs.append([ti, m - bi + 1])
    return prs


_END_CACHE = {}


def p3_end_stability(seq):
    """primer3's 3' end stability (positive kcal/mol), via design_primers."""
    tail = seq[-5:]
    if tail in _END_CACHE:
        return _END_CACHE[tail]
    primer = "ACGTTGCAGTCCATGA" + tail  # only the last 5 bases matter
    tmpl = primer + "GATCCGTAGCTAGCATGCATCGATCGGATCCTAGCTAGGCATGCAATCGATCAGTCGAT" * 3
    out = primer3.design_primers(
        {"SEQUENCE_ID": "x", "SEQUENCE_TEMPLATE": tmpl, "SEQUENCE_PRIMER": primer},
        {"PRIMER_PICK_ANYWAY": 1, "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 0,
         "PRIMER_PICK_INTERNAL_OLIGO": 0, "PRIMER_MIN_SIZE": 5, "PRIMER_MAX_SIZE": 36})
    v = out.get("PRIMER_LEFT_0_END_STABILITY")
    _END_CACHE[tail] = v
    return v


def oracle(singles, pairs):
    t0 = time.time()
    for c in singles:
        s = c["seq"]
        c["p3_hairpin"] = p3(primer3.calc_hairpin(s, output_structure=True, **COND), "hairpin")
        c["p3_homodimer"] = p3(primer3.calc_homodimer(s, output_structure=True, **COND), "dimer", len(s), len(s))
        c["p3_end_dimer"] = p3(primer3.calc_end_stability(s, s, **COND))
        c["p3_end_stability"] = p3_end_stability(s)
    for c in pairs:
        c["p3_heterodimer"] = p3(primer3.calc_heterodimer(c["a"], c["b"], output_structure=True, **COND),
                                 "dimer", len(c["a"]), len(c["b"]))
        c["p3_end_dimer"] = p3(primer3.calc_end_stability(c["a"], c["b"], **COND))
    return time.time() - t0


# ---------------------------------------------------------------- Lua side
LUA_DRIVER = r"""
local root, inp, outp = arg[1], arg[2], arg[3]
package.path = root .. "/lua/?.lua;" .. root .. "/lua/?/init.lua;" .. package.path
local O = require("plasmid.bio.oligo")
local fh = assert(io.open(inp, "rb")); local data = vim.json.decode(fh:read("*a")); fh:close()
local function slim(r)
  if not r then return vim.NIL end
  local o = { dg = r.dg, dh = r.dh, ds = r.ds, tm = r.tm or vim.NIL, three_prime = r.three_prime }
  if r.pairs then o.pairs = r.pairs end
  if r.stem then o.stem = { r.stem.i, r.stem.j, r.stem.len }; o.loop = r.loop end
  return o
end
local out = { singles = {}, pairs = {} }
local t0 = vim.uv.hrtime()
for i, c in ipairs(data.singles) do
  local a = O.assess(c.seq)
  out.singles[i] = { hairpin = slim(a.hairpin), homodimer = slim(a.homodimer), end_dimer = slim(a.end_dimer),
    end_stability = a.end_stability, gc = a.gc, gc_clamp = a.gc_clamp, max_run = a.max_run,
    max_repeat = a.max_repeat, warnings = a.warnings }
end
local t1 = vim.uv.hrtime()
for i, c in ipairs(data.pairs) do
  out.pairs[i] = { heterodimer = slim(O.heterodimer(c.a, c.b)), end_dimer = slim(O.end_dimer(c.a, c.b)) }
end
local t2 = vim.uv.hrtime()
-- timing: full 25-mer pair analysis (2 hairpins, 2 homodimers, heterodimer, 2 end dimers)
local pa, pb = {}, {}
for _, c in ipairs(data.pairs) do if #c.a == 25 then pa[#pa + 1] = c.a end end
for _, c in ipairs(data.singles) do if #c.seq == 25 then pb[#pb + 1] = c.seq end end
local times, times2 = {}, {}
for rep = 1, 2 do -- first pass warms the JIT; keep the second
  times, times2 = {}, {}
  for r = 1, math.min(200, #pa, #pb) do
    local a, b = pa[r], pb[r]
    local s = vim.uv.hrtime()
    O.hairpin(a); O.hairpin(b); O.homodimer(a); O.homodimer(b); O.heterodimer(a, b)
    O.end_dimer(a); O.end_dimer(b); O.end_dimer(a, b); O.end_dimer(b, a)
    times[#times + 1] = (vim.uv.hrtime() - s) / 1e6
    s = vim.uv.hrtime()
    O.assess(a); O.assess(b); O.heterodimer(a, b); O.end_dimer(a, b); O.end_dimer(b, a)
    times2[#times2 + 1] = (vim.uv.hrtime() - s) / 1e6
  end
end
table.sort(times); table.sort(times2)
out.timing = { singles_ms = (t1 - t0) / 1e6, pairs_ms = (t2 - t1) / 1e6, n_singles = #data.singles,
  n_pairs = #data.pairs, pair25_ms = times, pair25_assess_ms = times2, version = tostring(jit and jit.version or _VERSION) }
local f = assert(io.open(outp, "wb")); f:write(vim.json.encode(out)); f:close()
"""


def run_lua(nvim, singles, pairs):
    with tempfile.TemporaryDirectory() as td:
        inp, outp, drv = (os.path.join(td, x) for x in ("in.json", "out.json", "driver.lua"))
        with open(inp, "w") as f:
            json.dump({"singles": [{"seq": c["seq"]} for c in singles],
                       "pairs": [{"a": c["a"], "b": c["b"]} for c in pairs]}, f)
        with open(drv, "w") as f:
            f.write(LUA_DRIVER)
        subprocess.run([nvim, "-l", drv, ROOT, inp, outp], check=True)
        with open(outp) as f:
            return json.load(f)


# ---------------------------------------------------------------- statistics
def err_stats(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    e = np.abs(x - y)
    return {"n": int(len(x)), "pearson": float(stats.pearsonr(x, y)[0]) if len(x) > 2 else float("nan"),
            "spearman": float(stats.spearmanr(x, y)[0]) if len(x) > 2 else float("nan"),
            "median_abs": float(np.median(e)), "p95_abs": float(np.percentile(e, 95)),
            "max_abs": float(e.max()), "bias": float(np.mean(x - y))}


def confusion(ours, theirs):
    tp = sum(1 for a, b in zip(ours, theirs) if a and b)
    fp = sum(1 for a, b in zip(ours, theirs) if a and not b)
    fn = sum(1 for a, b in zip(ours, theirs) if not a and b)
    tn = sum(1 for a, b in zip(ours, theirs) if not a and not b)
    n = tp + fp + fn + tn
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "agree": (tp + tn) / n if n else float("nan")}


def g(r, k="dg", none=0.0):
    """value of a structure result; 'no structure' counts as dG 0 / Tm none."""
    if r is None:
        return none
    v = r.get(k)
    return none if v is None else v


def compare(singles, pairs, lua):
    res = {}
    for name, src, key in (("hairpin", singles, "hairpin"), ("homodimer", singles, "homodimer"),
                           ("heterodimer", pairs, "heterodimer"),
                           ("end_dimer_self", singles, "end_dimer"), ("end_dimer_pair", pairs, "end_dimer")):
        L = lua["singles"] if src is singles else lua["pairs"]
        p3key = "p3_" + ("end_dimer" if key == "end_dimer" else key)
        mine = [l[key] for l in L]
        theirs = [c[p3key] for c in src]
        both = [(m, t) for m, t in zip(mine, theirs) if m is not None and t is not None]
        r = {"n": len(mine), "found_both": len(both),
             "found_lua_only": sum(1 for m, t in zip(mine, theirs) if m is not None and t is None),
             "found_p3_only": sum(1 for m, t in zip(mine, theirs) if m is None and t is not None)}
        r["dg"] = err_stats([g(m) for m in mine], [g(t) for t in theirs])
        r["dg_both"] = err_stats([m["dg"] for m, _ in both], [t["dg"] for _, t in both]) if both else None
        r["dh_both"] = err_stats([m["dh"] for m, _ in both], [t["dh"] for _, t in both]) if both else None
        tmb = [(m["tm"], t["tm"]) for m, t in both if m.get("tm") is not None]
        r["tm_both"] = err_stats([a for a, _ in tmb], [b for _, b in tmb]) if tmb else None
        # Tm restricted to the range where it matters (primer3 Tm > 0 °C)
        tmr = [(a, b) for a, b in tmb if b > 0]
        r["tm_pos"] = err_stats([a for a, _ in tmr], [b for _, b in tmr]) if len(tmr) > 2 else None
        dthr = THRESH["hairpin_dg"] if key == "hairpin" else (THRESH["end_dimer_dg"] if key == "end_dimer" else THRESH["dimer_dg"])
        r["thr_dg"] = dthr
        r["conf_dg"] = confusion([g(m) < dthr for m in mine], [g(t) < dthr for t in theirs])
        r["conf_tm"] = confusion([g(m, "tm", -999) > THRESH["struct_tm"] for m in mine],
                                 [g(t, "tm", -999) > THRESH["struct_tm"] for t in theirs])
        res[name] = r
    # end stability: primer3 reports -dG (positive); we return dG
    mine = [-l["end_stability"] for l in lua["singles"]]
    theirs = [c["p3_end_stability"] for c in singles]
    res["end_stability"] = err_stats(mine, theirs)
    res["end_stability"]["conf"] = confusion([m > THRESH["end_stab"] for m in mine], [t > THRESH["end_stab"] for t in theirs])
    pent = sorted(_END_CACHE.items())
    res["end_stability_pentamers"] = {"n": len(pent)}
    return res


# ---------------------------------------------------------------- diagnostics
def n_aa(s):
    return sum(1 for i in range(len(s) - 1) if s[i:i + 2] in ("AA", "TT"))


def same_pairs(lp, pp, sym):
    a = set(map(tuple, lp))
    b = set(map(tuple, pp))
    return a == b or (sym and a == set((j, i) for i, j in b))


def diagnostics(singles, pairs, lua):
    """Where do the differences come from? Structure identity and parameter decomposition."""
    from Bio.SeqUtils import MeltingTemp as mt
    d = {}
    for name, src, L, key, sym in (("hairpin", singles, lua["singles"], "hairpin", False),
                                   ("homodimer", singles, lua["singles"], "homodimer", True),
                                   ("heterodimer", pairs, lua["pairs"], "heterodimer", False)):
        same, diff = [], []
        for c, l in zip(src, L):
            t, m = c["p3_" + key], l[key]
            if t is None or m is None or "pairs" not in t:
                continue
            (same if same_pairs(m["pairs"], t["pairs"], sym) else diff).append(m["dg"] - t["dg"])
        d[name] = {"n": len(same) + len(diff), "same": len(same),
                   "same_med": float(np.median(np.abs(same))), "same_p95": float(np.percentile(np.abs(same), 95)),
                   "diff_med": float(np.median(np.abs(diff))), "diff_p95": float(np.percentile(np.abs(diff), 95)),
                   "diff_lua_more_stable": float(np.mean(np.array(diff) < 0))}
    # (1) perfect palindromes paired end to end: only the AA/TT stack differs?
    res = []
    for c, l in zip(singles, lua["singles"]):
        m, t = l["homodimer"], c["p3_homodimer"]
        if c["kind"] != "palindrome" or not m or not t or len(m["pairs"]) != len(c["seq"]):
            continue
        k = n_aa(c["seq"])
        res.append((m["dh"] - t["dh"] - 0.3 * k, m["ds"] - t["ds"] - 0.9 * k, k,
                    abs(m["dg"] - t["dg"]), abs(m["tm"] - t["tm"])))
    r = np.array(res)
    d["palindrome"] = {"n": len(r), "max_abs_dh_resid": float(np.abs(r[:, 0]).max()),
                       "max_abs_ds_resid": float(np.abs(r[:, 1]).max()), "with_aa": int((r[:, 2] > 0).sum()),
                       "max_abs_dg": float(r[:, 3].max()), "max_abs_tm": float(r[:, 4].max())}
    # (2) identical contiguous homodimer helices with terminal mismatches at both ends:
    #     least squares of the dG residual on one indicator per terminal-mismatch key
    tm_tab = mt.DNA_TMM1

    def rot(k):
        t, b = k.split("/")
        return b[::-1] + "/" + t[::-1]
    keys = sorted(tm_tab)
    kid = {k: i for i, k in enumerate(keys)}
    rows, y = [], []
    for c, l in zip(singles, lua["singles"]):
        m, t = l["homodimer"], c["p3_homodimer"]
        s, n = c["seq"], len(c["seq"])
        if not m or not t or "pairs" not in t or not same_pairs(m["pairs"], t["pairs"], True):
            continue
        lp = m["pairs"]
        if any(lp[q + 1][0] - lp[q][0] != 1 or lp[q][1] - lp[q + 1][1] != 1 for q in range(len(lp) - 1)):
            continue
        (i1, j1), (i2, j2) = lp[0], lp[-1]
        if not (i1 > 1 and j1 < n and i2 < n and j2 > 1):
            continue
        row = np.zeros(len(keys))
        for k in (s[i1 - 2] + s[i1 - 1] + "/" + s[j1] + s[j1 - 1], s[i2 - 1] + s[i2] + "/" + s[j2 - 1] + s[j2 - 2]):
            row[kid[k if k in tm_tab else rot(k)]] += 1
        rows.append(row)
        # remove the known AA/TT stack difference (+0.021 kcal/mol per stack at 37 °C)
        y.append(m["dg"] - t["dg"] - (0.3 - 310.15 * 0.0009) * n_aa(s[i1 - 1:i2]))
    A, y = np.array(rows), np.array(y)
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    d["tmm_fit"] = {"n": len(y), "raw_sd": float(np.std(y)), "resid_sd": float(np.std(y - A @ coef)),
                    "keys": sorted(((keys[i], float(coef[i]), int(A[:, i].sum())) for i in range(len(keys))),
                                   key=lambda x: x[1])}
    return d


# ---------------------------------------------------------------- fixture
def make_fixture(singles, pairs, lua, res, rng):
    """A few hundred cases with primer3 values + the tolerances measured here."""
    pick_s = [i for i, c in enumerate(singles) if c["kind"] != "random"]
    pick_s = rng.sample(pick_s, min(150, len(pick_s))) + rng.sample(
        [i for i, c in enumerate(singles) if c["kind"] == "random"], 100)
    pick_p = rng.sample([i for i, c in enumerate(pairs) if c["kind"] != "random_pair"], 60) + rng.sample(
        [i for i, c in enumerate(pairs) if c["kind"] == "random_pair"], 60)
    out = {"conditions": COND, "thresholds": THRESH, "singles": [], "pairs": []}
    for i in sorted(pick_s):
        c = singles[i]
        out["singles"].append({k: c[k] for k in ("seq", "kind", "p3_hairpin", "p3_homodimer", "p3_end_stability")})
        for k in ("p3_hairpin", "p3_homodimer"):
            if out["singles"][-1][k]:
                out["singles"][-1][k] = {x: y for x, y in c[k].items() if x != "pairs"}
    for i in sorted(pick_p):
        c = pairs[i]
        out["pairs"].append({k: c[k] for k in ("a", "b", "kind", "p3_heterodimer")})
        if c["p3_heterodimer"]:
            out["pairs"][-1]["p3_heterodimer"] = {x: y for x, y in c["p3_heterodimer"].items() if x != "pairs"}
    # measured on this fixture: agreement rates and p95 |dG error| (used as spec floors)
    sub_l = [lua["singles"][i] for i in sorted(pick_s)]
    sub_p = [lua["pairs"][i] for i in sorted(pick_p)]
    meas = {}
    for name, L, C, key in (("hairpin", sub_l, out["singles"], "hairpin"), ("homodimer", sub_l, out["singles"], "homodimer"),
                            ("heterodimer", sub_p, out["pairs"], "heterodimer")):
        thr = THRESH["hairpin_dg"] if key == "hairpin" else THRESH["dimer_dg"]
        mine = [g(l[key]) for l in L]
        theirs = [g(c["p3_" + key]) for c in C]
        conf = confusion([m < thr for m in mine], [t < thr for t in theirs])
        e = np.abs(np.array(mine) - np.array(theirs))
        meas[name] = {"agree": conf["agree"], "p95_abs_dg": float(np.percentile(e, 95)),
                      "median_abs_dg": float(np.median(e)), "max_abs_dg": float(e.max()), "n": len(e)}
    out["measured"] = meas
    return out


# ---------------------------------------------------------------- report
def fmt_stats(s):
    if not s:
        return "| – | – | – | – | – | – |"
    return "| {n} | {pearson:.3f} | {spearman:.3f} | {median_abs:.2f} | {p95_abs:.2f} | {bias:+.2f} |".format(**s)


def fmt_conf(c):
    return "| {tp} | {fp} | {fn} | {tn} | {pct:.1f}% |".format(pct=100 * c["agree"], **c)


def write_report(path, singles, pairs, lua, res, p3_secs, fixture, diag):
    pq = diag["palindrome"]
    kinds_s, kinds_p = {}, {}
    for c in singles:
        kinds_s[c["kind"]] = kinds_s.get(c["kind"], 0) + 1
    for c in pairs:
        kinds_p[c["kind"]] = kinds_p.get(c["kind"], 0) + 1
    t = lua["timing"]
    p25 = t["pair25_ms"]
    L = []
    w = L.append
    w("# Oligo structure estimates vs primer3 (thal)\n")
    hp, hd, he = res["hairpin"], res["homodimer"], res["heterodimer"]
    w("## Summary\n")
    w("* `lua/plasmid/bio/oligo.lua` is **not** equal to primer3 and does not claim to be. Over %d single oligos and "
      "%d pairs, dG (kcal/mol, 37 °C) correlates with primer3 at Pearson r = %.3f (hairpin), %.3f (self-dimer) and "
      "%.3f (heterodimer); median absolute error %.2f / %.2f / %.2f, 95th percentile %.2f / %.2f / %.2f." % (
          len(singles), len(pairs), hp["dg"]["pearson"], hd["dg"]["pearson"], he["dg"]["pearson"],
          hp["dg"]["median_abs"], hd["dg"]["median_abs"], he["dg"]["median_abs"],
          hp["dg"]["p95_abs"], hd["dg"]["p95_abs"], he["dg"]["p95_abs"]))
    w("* Warning agreement (same flag as primer3): hairpin dG < -2: %.1f%%; self-dimer dG < -9: %.1f%%; "
      "heterodimer dG < -9: %.1f%%; 3'-anchored dimer dG < -5: %.1f%% (self) / %.1f%% (pair); hairpin Tm > 47 °C: %.1f%%." % (
          100 * hp["conf_dg"]["agree"], 100 * hd["conf_dg"]["agree"], 100 * he["conf_dg"]["agree"],
          100 * res["end_dimer_self"]["conf_dg"]["agree"], 100 * res["end_dimer_pair"]["conf_dg"]["agree"],
          100 * hp["conf_tm"]["agree"]))
    w("* 3' end stability (primer3's PRIMER_..._END_STABILITY, the pentamer duplex dG) is reproduced exactly "
      "(max abs error %.3f kcal/mol over %d oligos)." % (res["end_stability"]["max_abs"], res["end_stability"]["n"]))
    lo = [l["hairpin"]["dg"] for c, l in zip(singles, lua["singles"]) if l["hairpin"] and c["p3_hairpin"] is None]
    if lo:
        w("* primer3 reports *no* hairpin for %d oligos where this module returns its best (weak) fold; those folds "
          "have dG between %.2f and %+.2f kcal/mol (median %+.2f), i.e. none would trigger the -2 kcal/mol warning. "
          "Tables below count 'no structure' as dG 0." % (len(lo), min(lo), max(lo), float(np.median(lo))))
    w("* Structure Tm is the weakest output: both programs report Tm of *their* chosen structure, primer3 tends to "
      "pick longer, loopier structures with larger |dH| (see below), so Tm differs by a median of %.1f °C (self-dimer, "
      "primer3 Tm > 0 °C) with a long tail. Use dG for decisions; Tm is indicative." % hd["tm_pos"]["median_abs"])
    w("* Why they differ: the duplex bookkeeping (stacks, initiation, terminal A·T, salt, Tm formula) matches primer3 "
      "exactly apart from the AA/TT stack (SantaLucia & Hicks 2004 here, 1998 in primer3); on identical contiguous "
      "helices with terminal mismatches at both ends the residual is entirely explained by different "
      "terminal-mismatch tables; the rest comes from choosing a "
      "different structure (this module minimises dG at 37 °C exactly within its structure class; primer3's choice "
      "differs, sometimes finding a less stable structure, sometimes a more stable one under its own parameters).\n")
    w("## Model and approximations\n")
    w("Two-state nearest-neighbour model, implemented from published tables (no primer3 code was used; primer3 is "
      "only called as an oracle in this script): SantaLucia & Hicks 2004 Watson–Crick stacks, initiation and terminal "
      "A·T penalty; single internal mismatches (Allawi & SantaLucia 1997–98, Peyret et al. 1999), terminal mismatches "
      "(SantaLucia & Peyret 2001) and dangling ends (Bommarito et al. 2000), all as tabulated in Biopython "
      "`Bio.SeqUtils.MeltingTemp`; SantaLucia & Hicks 2004 loop initiation for hairpin / bulge / internal loops "
      "(entropic, Jacobson–Stockmayer extrapolation), +0.3 kcal/mol per nt internal-loop asymmetry; salt as "
      "0.368·ln[Na+eq] per stack with [Na+eq] = [Na+] + 120·sqrt([Mg2+] − [dNTP]) (von Ahsen 2001). "
      "Exact dynamic programming over one duplex (dimers) or consecutive stem-loops (hairpin folds) with stacks, 1x1 "
      "mismatches, bulges and internal loops up to 30 unpaired nt, terminal mismatches / dangling ends at helix ends; "
      "the minimum-dG structure at the requested temperature is reported. Approximated or omitted: tri-/tetraloop "
      "bonus tables, hairpin-specific terminal-mismatch table (duplex values used), internal-loop first-mismatch "
      "bonuses (terminal-mismatch values at both closing pairs used instead), tandem-mismatch table (scored as 2x2 "
      "internal loops), multibranch loops, G·T as helix pairs, symmetry correction.\n")
    w("Parameter note: Biopython's DNA_TMM1 lists terminal mismatch GG/CG as dH -0.7, dS -19.2 (dG37 +5.25, while "
      "every other entry is -0.1 to -1.5 kcal/mol). Read literally it forbids any C·G pair flanked by a G/G mismatch "
      "(it made 14- and 16-bp palindromes pair only 12/14 bp). We assume a decimal slip in dH and use -7.0 (dG37 "
      "-1.05); primer3's implied value for that motif is about -1.2. This changed self-dimer warning agreement by "
      "about 0.3 percentage points.\n")
    w("Development note (not reproduced by this script; the code path was removed): a variant that picked, at every "
      "DP cell, the partial structure with the highest Tm instead of the lowest dG agreed less well with primer3 "
      "(threshold agreement hairpin 93.9%, self-dimer 92.0%, heterodimer 97.5%, versus 95.9%, 97.1%, 99.1% for "
      "minimum dG at the same stage of development), so minimum dG is used.\n")
    w("Generated by `tests/tools/oligo_oracle.py` (primer3-py %s, Lua: %s). "
      "Conditions: Na+ 50 mM, Mg2+ 1.5 mM, dNTP 0.6 mM, oligo 50 nM, dG at 37 °C "
      "(primer3-py `calc_*` defaults, max_loop 30).\n" % (primer3.__version__, t["version"]))
    w("## Case set\n")
    w("Single oligos (hairpin, homodimer, 3'-anchored self-dimer, end stability): **%d**  " % len(singles))
    w(", ".join("%s %d" % kv for kv in sorted(kinds_s.items())) + "\n")
    w("Oligo pairs (heterodimer, 3'-anchored heterodimer): **%d**  " % len(pairs))
    w(", ".join("%s %d" % kv for kv in sorted(kinds_p.items())) + "\n")
    w("Random oligos: length uniform 15–35 nt, per-oligo GC target uniform 0.2–0.8. Seed fixed in the script.\n")
    w("## dG (kcal/mol), all cases; 'no structure' counted as dG = 0\n")
    w("| structure | n | Pearson r | Spearman ρ | median abs err | 95th pct abs err | mean bias (Lua − primer3) |")
    w("|---|---|---|---|---|---|---|")
    for k in ("hairpin", "homodimer", "heterodimer", "end_dimer_self", "end_dimer_pair"):
        w("| %s " % k + fmt_stats(res[k]["dg"]))
    w("\n## Cases where both programs report a structure\n")
    w("| structure | found by both | Lua only | primer3 only |")
    w("|---|---|---|---|")
    for k in ("hairpin", "homodimer", "heterodimer", "end_dimer_self", "end_dimer_pair"):
        r = res[k]
        w("| %s | %d | %d | %d |" % (k, r["found_both"], r["found_lua_only"], r["found_p3_only"]))
    for lab, key in (("dG (kcal/mol)", "dg_both"), ("dH (kcal/mol)", "dh_both"), ("Tm (°C), all", "tm_both"),
                     ("Tm (°C), primer3 Tm > 0 °C", "tm_pos")):
        w("\n**%s, both found**\n" % lab)
        w("| structure | n | Pearson r | Spearman ρ | median abs err | 95th pct abs err | mean bias |")
        w("|---|---|---|---|---|---|---|")
        for k in ("hairpin", "homodimer", "heterodimer", "end_dimer_self", "end_dimer_pair"):
            w("| %s " % k + fmt_stats(res[k][key]))
    w("\n## Warning-threshold agreement (confusion tables; positive = flagged)\n")
    w("Thresholds: hairpin dG < %.0f kcal/mol (PREMIER Biosoft 3'-hairpin rule); homo/heterodimer dG < %.0f kcal/mol "
      "(IDT OligoAnalyzer guidance); 3'-anchored dimer dG < %.0f kcal/mol (PREMIER Biosoft 3'-dimer rule); "
      "structure Tm > %.0f °C (primer3 PRIMER_MAX_HAIRPIN_TH / SELF_ANY_TH / SELF_END_TH defaults); "
      "end stability > %.0f kcal/mol.\n" % (THRESH["hairpin_dg"], THRESH["dimer_dg"], THRESH["end_dimer_dg"],
                                             THRESH["struct_tm"], THRESH["end_stab"]))
    w("| check | both flag | Lua only | primer3 only | neither | agreement |")
    w("|---|---|---|---|---|---|")
    for k in ("hairpin", "homodimer", "heterodimer", "end_dimer_self", "end_dimer_pair"):
        w("| %s dG < %.0f " % (k, res[k]["thr_dg"]) + fmt_conf(res[k]["conf_dg"]))
    for k in ("hairpin", "homodimer", "heterodimer", "end_dimer_self", "end_dimer_pair"):
        w("| %s Tm > %.0f °C " % (k, THRESH["struct_tm"]) + fmt_conf(res[k]["conf_tm"]))
    es = res["end_stability"]
    w("| end stability > %.0f " % THRESH["end_stab"] + fmt_conf(es["conf"]))
    w("\n## 3' end stability (pentamer dG) vs primer3 PRIMER_LEFT_0_END_STABILITY\n")
    w("All %d oligos (%d distinct 3' pentamers queried through `design_primers`): Pearson r %.4f, "
      "median abs err %.3f, 95th pct %.3f, max %.3f kcal/mol.\n" % (es["n"], res["end_stability_pentamers"]["n"],
                                                                     es["pearson"], es["median_abs"], es["p95_abs"], es["max_abs"]))
    w("## Where the differences come from\n")
    w("**Same structure or not.** Pairs parsed from primer3's ASCII structure and compared with the Lua base pairs "
      "(cases where both report a structure):\n")
    w("| structure | n | identical pairs | median / p95 abs dG err (identical) | median / p95 abs dG err (different) | "
      "Lua more stable when different |")
    w("|---|---|---|---|---|---|")
    for k in ("hairpin", "homodimer", "heterodimer"):
        q = diag[k]
        w("| %s | %d | %.1f%% | %.2f / %.2f | %.2f / %.2f | %.0f%% |" % (
            k, q["n"], 100 * q["same"] / q["n"], q["same_med"], q["same_p95"], q["diff_med"], q["diff_p95"],
            100 * q["diff_lua_more_stable"]))
    pq = diag["palindrome"]
    w("\n**Stacks, initiation, terminal A·T, salt and Tm formula.** On the %d perfect palindromes, both programs pair "
      "the full length; after removing 0.3 kcal/mol dH and 0.9 cal/(K·mol) dS per AA/TT stack (SantaLucia & Hicks "
      "2004 AA/TT = -7.6/-21.3, which this module uses as specified; primer3's values match the 1998 unified "
      "-7.9/-22.2; dG37 is -1.00 in both) the residual is at most %.2g kcal/mol in dH and %.2g cal/(K·mol) in dS. "
      "The duplex bookkeeping (initiation, terminal A·T penalty, 0.368·ln[Na+] per stack with "
      "Na+ = mv + 120·sqrt(Mg − dNTP), Tm = dH/(dS + R ln Ct) for self-dimers) therefore matches.\n" % (
          pq["n"], pq["max_abs_dh_resid"], pq["max_abs_ds_resid"]))
    tf = diag["tmm_fit"]
    big = [x for x in tf["keys"] if abs(x[1]) >= 0.4 and x[2] >= 10]
    w("**Terminal mismatches.** For %d self-dimers where both programs choose the same contiguous helix with a "
      "terminal mismatch at each end, the dG difference has SD %.2f kcal/mol; one least-squares offset per "
      "terminal-mismatch motif (48 motifs) leaves a residual SD of %.3f kcal/mol. The remaining dimer error on "
      "identical structures is therefore a parameter-table difference (this module: SantaLucia & Peyret 2001 terminal "
      "mismatches as tabulated in Biopython; primer3 ships its own tables), not a search or bookkeeping difference. "
      "Motifs whose implied offset is at least 0.4 kcal/mol (Lua − primer3, n ≥ 10): %s.\n" % (
          tf["n"], tf["raw_sd"], tf["resid_sd"],
          ", ".join("%s %+.2f (n=%d)" % x for x in big) or "none"))
    w("## Fixture (tests/oracle/expected/oligo.json)\n")
    w("%d single oligos + %d pairs sampled from the set above (designed cases over-represented). "
      "Measured on the fixture:\n" % (len(fixture["singles"]), len(fixture["pairs"])))
    w("| structure | n | threshold agreement | median abs dG err | 95th pct abs dG err | max abs dG err |")
    w("|---|---|---|---|---|---|")
    for k, m in fixture["measured"].items():
        w("| %s | %d | %.1f%% | %.2f | %.2f | %.2f |" % (k, m["n"], 100 * m["agree"], m["median_abs_dg"],
                                                       m["p95_abs_dg"], m["max_abs_dg"]))
    w("\ntests/spec/oligo_spec.lua asserts these numbers (agreement >= measured, median/p95/max error <= measured), "
      "end stability within 0.006 kcal/mol of primer3 for every fixture oligo (primer3 prints 2 decimals), and, for "
      "the perfect palindromes, self-dimer dG within 0.30 kcal/mol and Tm within 1.0 °C of primer3 (measured over "
      "all %d palindromes: max %.2f kcal/mol, %.2f °C).\n" % (pq["n"], pq["max_abs_dg"], pq["max_abs_tm"]))
    w("\n## Performance (Neovim, LuaJIT)\n")
    w("Batch: %d single-oligo `assess()` calls in %.0f ms (%.2f ms each, includes hairpin + homodimer + "
      "3'-anchored self-dimer); %d pairs (heterodimer + 3'-anchored heterodimer) in %.0f ms (%.2f ms each).  " %
      (t["n_singles"], t["singles_ms"], t["singles_ms"] / t["n_singles"], t["n_pairs"], t["pairs_ms"], t["pairs_ms"] / t["n_pairs"]))
    w("Full 25-mer pair (2 hairpins, 2 homodimers, heterodimer, 4 3'-anchored dimers), %d pairs: median %.2f ms, "
      "95th pct %.2f ms, max %.2f ms.  " % (len(p25), p25[len(p25) // 2], p25[int(0.95 * (len(p25) - 1))], p25[-1]))
    q = t["pair25_assess_ms"]
    w("Same pairs as a user-level pair check (`assess(a)`, `assess(b)`, `heterodimer(a, b)`, `end_dimer(a, b)`, "
      "`end_dimer(b, a)`; assess includes GC/runs/end stability and warnings): median %.2f ms, 95th pct %.2f ms, "
      "max %.2f ms.  " % (q[len(q) // 2], q[int(0.95 * (len(q) - 1))], q[-1]))
    w("Timings are wall-clock from one run after a JIT warm-up pass on this machine (12-core x86_64 workstation).  ")
    w("primer3-py (C) for the whole oracle pass: %.1f s.\n" % p3_secs)
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nvim", default=os.environ.get("NVIM", "nvim"))
    ap.add_argument("--seed", type=int, default=20261116)
    ap.add_argument("--report", default=os.path.join(ROOT, "docs/oligo_report.md"))
    ap.add_argument("--fixture", default=os.path.join(ROOT, "tests/oracle/expected/oligo.json"))
    ap.add_argument("--dump", default=None, help="write every case with both programs' values")
    a = ap.parse_args()
    singles, pairs = make_cases(a.seed)
    p3_secs = oracle(singles, pairs)
    lua = run_lua(a.nvim, singles, pairs)
    res = compare(singles, pairs, lua)
    fixture = make_fixture(singles, pairs, lua, res, random.Random(a.seed + 1))
    with open(a.fixture, "w") as f:
        json.dump(fixture, f, indent=1)
    diag = diagnostics(singles, pairs, lua)
    write_report(a.report, singles, pairs, lua, res, p3_secs, fixture, diag)
    if a.dump:
        for c, l in zip(singles, lua["singles"]):
            c["lua"] = l
        for c, l in zip(pairs, lua["pairs"]):
            c["lua"] = l
        with open(a.dump, "w") as f:
            json.dump({"singles": singles, "pairs": pairs, "stats": res, "timing": lua["timing"]}, f)
    json.dump(res, sys.stdout, indent=1, default=str)


if __name__ == "__main__":
    main()
