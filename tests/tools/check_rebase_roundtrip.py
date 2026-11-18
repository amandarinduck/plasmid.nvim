#!/usr/bin/env python3
"""Dev-only check for lua/plasmid/rebase.lua against Biopython.

Takes every enzyme in Biopython's Bio.Restriction.Restriction_Dictionary and
writes it back out as EMBOSS emboss_e / emboss_r records. This is the inverse
of Biopython's Scripts/Restriction/ranacompiler.py:

    fst5 = c1 (+1 if c1 < 0)          ->  c1 = cut5 if cut5 > 0 else cut5 - 1
    fst3 = c2 - size (+1 if c2 < 0)   ->  c2 = cut3 if cut3 > 0 else cut3 - 1
                                          with cut3 = size + fst3
    TwoCuts -> ncuts 4, c3/c4 from scd5/scd3;  NoCut -> ncuts 0, all cuts 0

The records are then imported with rebase.import() under `nvim -l`. The check
confirms that:
  * every single-cut enzyme with an IUPAC site comes back with
    site == Enzyme.site, cut5 == fst5, cut3 == size + fst3, and
    suppliers == "".join(suppl);
  * every two-cut, no-cut and non-IUPAC enzyme is skipped and counted under
    the right reason;
  * the commercial count equals the number of enzymes that have suppliers.
The generated files are temporary and are not kept. No REBASE file is read.

usage: python tests/tools/check_rebase_roundtrip.py [--root PLUGIN_ROOT] [--nvim PATH]
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

import Bio
from Bio.Restriction.Restriction_Dictionary import rest_dict

IUPAC = set("ACGTRYSWKMBDHVN")

LUA = r"""
local root, e_path, r_path, store, out = _G.arg[1], _G.arg[2], _G.arg[3], _G.arg[4], _G.arg[5]
package.path = root .. "/lua/?.lua;" .. package.path
local rebase = require("plasmid.rebase")
rebase.path = store
local list, summary = rebase.import(e_path, r_path)
assert(list, summary)
rebase.reset()
local reloaded = rebase.load()
local fh = assert(io.open(out, "w"))
fh:write(vim.json.encode({ list = list, summary = summary, reloaded = #reloaded }))
fh:close()
rebase.clear()
"""


def emboss_pos(cut):
    """Plugin/Biopython gap index -> EMBOSS position (which has no 0)."""
    return cut if cut > 0 else cut - 1


def to_emboss(name, d):
    site, size = d["site"], d["size"]
    if d["fst5"] is None:
        return f"{name}\t{site}\t{size}\t0\t0\t0\t0\t0\t0", "no_cut"
    c1 = emboss_pos(d["fst5"])
    c2 = emboss_pos(size + d["fst3"])
    blunt = int(d["fst5"] == size + d["fst3"])
    if d["scd5"] is not None:
        c3 = emboss_pos(d["scd5"])
        c4 = emboss_pos(size + d["scd3"])
        return f"{name}\t{site}\t{size}\t4\t{blunt}\t{c1}\t{c2}\t{c3}\t{c4}", "two_cuts"
    kind = "one" if set(site) <= IUPAC else "bad_site"
    return f"{name}\t{site}\t{size}\t2\t{blunt}\t{c1}\t{c2}\t0\t0", kind


def r_record(name, d):
    return "\n".join([name, "organism", "", "", "source", "".join(d["suppl"]), "0", "//"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nvim", default=os.environ.get("NVIM_BIN", "/opt/nvim-linux-x86_64/bin/nvim"))
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    root = os.path.abspath(a.root)

    e_lines, r_blocks, kinds = [], [], {}
    for name in sorted(rest_dict):
        d = rest_dict[name]
        line, kind = to_emboss(name, d)
        e_lines.append(line)
        r_blocks.append(r_record(name, d))
        kinds[name] = kind

    with tempfile.TemporaryDirectory() as td:
        e_path, r_path = os.path.join(td, "emboss_e.gen"), os.path.join(td, "emboss_r.gen")
        lua, out = os.path.join(td, "run.lua"), os.path.join(td, "out.json")
        with open(e_path, "w") as fh:
            fh.write("# generated from Biopython Restriction_Dictionary (dev check)\n")
            fh.write("\n".join(e_lines) + "\n")
        with open(r_path, "w") as fh:
            fh.write("# generated from Biopython Restriction_Dictionary (dev check)\n")
            fh.write("\n".join(r_blocks) + "\n")
        with open(lua, "w") as fh:
            fh.write(LUA)
        subprocess.run([a.nvim, "--clean", "-l", lua, root, e_path, r_path,
                        os.path.join(td, "store", "rebase_enzymes.json"), out],
                       check=True, capture_output=True, timeout=120)
        res = json.load(open(out))

    got = {e["name"]: e for e in res["list"]}
    summary = res["summary"]
    errors = []
    expect_one = [n for n, k in kinds.items() if k == "one"]
    compared = 0
    for n in expect_one:
        d, e = rest_dict[n], got.get(n)
        if e is None:
            errors.append(f"{n}: missing after import")
            continue
        compared += 1
        want = (d["site"], d["fst5"], d["size"] + d["fst3"], "".join(d["suppl"]))
        have = (e["site"], e["cut5"], e["cut3"], e.get("suppliers"))
        if want != have:
            errors.append(f"{n}: want {want} got {have}")
    extra = set(got) - set(expect_one)
    if extra:
        errors.append(f"imported enzymes that should have been skipped: {sorted(extra)}")
    want_skipped = {k: sum(1 for v in kinds.values() if v == k) for k in ("two_cuts", "no_cut", "bad_site")}
    for k, v in want_skipped.items():
        if summary["skipped"].get(k) != v:
            errors.append(f"skipped.{k}: want {v} got {summary['skipped'].get(k)}")
    for k in ("duplicate", "malformed"):
        if summary["skipped"].get(k) != 0:
            errors.append(f"skipped.{k}: want 0 got {summary['skipped'].get(k)}")
    want_comm = sum(1 for n in expect_one if rest_dict[n]["suppl"])
    if summary["commercial"] != want_comm:
        errors.append(f"commercial: want {want_comm} got {summary['commercial']}")
    names = [e["name"] for e in res["list"]]
    if names != sorted(names):  # Python str order == byte order == Lua string <
        errors.append("imported list not sorted by name")
    if res["reloaded"] != len(names):
        errors.append(f"load() after reset returned {res['reloaded']} entries, want {len(names)}")

    bad = [n for n, k in kinds.items() if k == "bad_site"]
    print(f"Biopython {Bio.__version__}: {len(rest_dict)} enzymes in Restriction_Dictionary")
    print(f"converted to EMBOSS records: one-cut {len(expect_one)}, two-cut {want_skipped['two_cuts']}, "
          f"no-cut {want_skipped['no_cut']}, non-IUPAC site {want_skipped['bad_site']} "
          f"{[(n, rest_dict[n]['site']) for n in bad]}")
    print(f"rebase.import summary: count {summary['count']}, skipped {summary['skipped']}, "
          f"commercial {summary['commercial']}")
    print(f"single-cut enzymes compared (site, cut5 == fst5, cut3 == size + fst3, suppliers): {compared}")
    for e in errors:
        print("FAIL", e)
    print("PASS" if not errors else f"FAILED ({len(errors)})")
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
