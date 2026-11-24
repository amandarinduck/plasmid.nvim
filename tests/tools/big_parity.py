#!/usr/bin/env python3
"""Dev-only check: the windowed linear layout against the eager one, on records
too large to commit.

tests/spec/window_spec.lua already does this over the fixture corpus in one
process. The records this script is for are bigger than that approach can
afford -- the eager layout of E. coli K-12 is about half a gigabyte of Lua heap
-- so each mode runs in its own `nvim -l tests/tools/big_dump.lua` and writes a
JSON dump of the same windows of the record. This compares the two dumps row by
row: text, row kind, base span, highlight runs (which is where selection and
search hits live), feature segments, enzyme segments, overflow lists, and the
answers to pos_of / bp_at / feature_at for the window's bases.

A difference is printed as the first row that differs, with the two texts, so a
divergence can be cut down to a single block.

usage: python tests/tools/big_parity.py FILE [key=value ...] [--nvim PATH]
       python tests/tools/big_parity.py --matrix          every case, every record
"""
import argparse
import json
import os
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]

# The axes the layout has: both sequence layouts, two widths, the track limit
# including "all", and each overlay on its own so a failure names one of them.
MATRIX = [
    ("G10", "G=10 width=120"),
    ("G0", "G=0 width=120"),
    ("narrow", "G=10 width=80"),
    ("wide", "G=10 width=200"),
    ("tracks1", "G=10 width=120 tracks=1"),
    ("tracks_all", "G=10 width=120 tracks=all"),
    ("enz_unique", "G=10 width=120 enz=unique"),
    ("enz_all", "G=0 width=120 enz=all"),
    ("primers", "G=10 width=120 primers=1"),
    ("orfs", "G=10 width=120 orfs=300"),
    ("search", "G=10 width=120 search=GGATCC"),
    ("sel", "G=10 width=120 sel=1000..4000"),
    ("churn", "G=10 width=120 churn=1"),
    ("wrap", "G=10 width=120 wrap=1 windows=origin"),
]


def dump(nvim, path, mode, args, out):
    cmd = [nvim, "-l", "tests/tools/big_dump.lua", str(path), str(out), "mode=" + mode] + args
    p = subprocess.run(cmd, cwd=ROOT, stdin=subprocess.DEVNULL,
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit("%s %s failed:\n%s\n%s" % (path, mode, p.stdout, p.stderr))
    with open(out) as fh:
        d = json.load(fh)
    d["stdout"] = p.stdout.strip()
    return d


def diff_window(we, ww):
    """First difference between the eager and the windowed dump of one window."""
    if (we["lo_row"], we["hi_row"]) != (ww["lo_row"], ww["hi_row"]):
        return "window %s spans rows %s..%s eager, %s..%s windowed" % (
            we["name"], we["lo_row"], we["hi_row"], ww["lo_row"], ww["hi_row"])
    for a, b in zip(we["rows"], ww["rows"]):
        if a != b:
            keys = [k for k in set(a) | set(b) if a.get(k) != b.get(k)]
            lines = ["window %s row %d differs in %s" % (we["name"], a["row"], ", ".join(sorted(keys)))]
            for k in sorted(keys):
                lines.append("    eager  %-10s %s" % (k, json.dumps(a.get(k))[:400]))
                lines.append("    window %-10s %s" % (k, json.dumps(b.get(k))[:400]))
            return "\n".join(lines)
    for key in ("pos", "feats"):
        for a, b in zip(we[key], ww[key]):
            if a != b:
                return "window %s %s: eager %s, windowed %s" % (we["name"], key, a, b)
        if len(we[key]) != len(ww[key]):
            return "window %s %s: %d eager, %d windowed" % (we["name"], key, len(we[key]), len(ww[key]))
    return None


def compare(nvim, path, args, keep=None):
    with tempfile.TemporaryDirectory() as td:
        e = dump(nvim, path, "eager", args, os.path.join(td, "eager.json"))
        w = dump(nvim, path, "window", args, os.path.join(td, "window.json"))
    problems = []
    for k in ("nlines", "header_rows", "W", "G", "prefix_w", "content_w", "ntracks", "title", "len"):
        if e[k] != w[k]:
            problems.append("%s: eager %r, windowed %r" % (k, e[k], w[k]))
    # position queries over the whole record, answered before a row exists
    if w["built_after_positions"] != 0:
        problems.append("%d blocks were built by position queries alone" % w["built_after_positions"])
    for a, b in zip(e["global"], w["global"]):
        if a != b:
            problems.append("bp %d: eager row/col/bp_at/block %s, windowed %s" % (a[0], a[1:], b[1:]))
            break
    for we, ww in zip(e["windows"], w["windows"]):
        d = diff_window(we, ww)
        if d:
            problems.append(d)
    rows = sum(len(x["rows"]) for x in w["windows"])
    # a case that compared no rows is a hole in the harness, not a pass
    if rows == 0:
        problems.append("no rows were compared: windows %r" % [x["name"] for x in w["windows"]])
    return problems, rows, e, w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", nargs="?")
    ap.add_argument("args", nargs="*")
    ap.add_argument("--nvim", default=os.environ.get("NVIM", "nvim"))
    ap.add_argument("--matrix", action="store_true")
    ap.add_argument("--big", default=str(ROOT / "big"))
    a = ap.parse_args()

    cases = []
    if a.matrix:
        records = [f for pat in ("*.gb", "*.gbk", "*.fasta", "*.fa", "*.dna")
                   for f in pathlib.Path(a.big).glob(pat)]
        for f in sorted(records):
            # a multi-record file is checked record by record
            recs = [1, 2, 3, 4] if f.name == "phages.gbk" else [1]
            for name, args in MATRIX:
                # every site of every enzyme over a megabase is an overlay no
                # user asks for and minutes of restriction search per case; the
                # dense-overlay axis is covered by the smaller records
                if name == "enz_all" and f.stat().st_size > 2_000_000:
                    continue
                for r in recs:
                    cases.append((f, name, args.split() + ["rec=%d" % r]))
    else:
        cases.append((pathlib.Path(a.file), "cli", a.args))

    bad = 0
    total_rows = 0
    for f, name, args in cases:
        problems, rows, e, w = compare(a.nvim, f, args)
        total_rows += rows
        tag = "%s %s %s" % (f.name, name, " ".join(args))
        if problems:
            bad += 1
            print("DIFFER %s" % tag)
            for p in problems[:4]:
                print("  " + p.replace("\n", "\n  "))
        else:
            print("ok     %-44s %6d rows  eager %5.0f ms/%4.0f MB  window %5.0f ms/%4.0f MB" % (
                tag, rows, e["build_ms"], e["heap_mb"], w["build_ms"], w["heap_mb"]))
        sys.stdout.flush()
    print("%d cases, %d differ, %d rows compared" % (len(cases), bad, total_rows))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
