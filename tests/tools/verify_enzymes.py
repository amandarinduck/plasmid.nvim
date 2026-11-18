#!/usr/bin/env python3
"""Dev-only checker for lua/plasmid/data/enzymes.lua (plasmid.nvim).

Loads the Lua table through Neovim (`nvim -l`, falls back to a regex parser),
then compares every entry with Biopython's Bio.Restriction:

    site == Enzyme.site
    cut5 == Enzyme.fst5
    cut3 == Enzyme.size + Enzyme.fst3

It also rejects enzymes Biopython does not know, enzymes Biopython classes as
TwoCuts / NoCut / Unknown cleavage, duplicate names, unsorted order and
unknown fields.  Exit status 0 = all entries pass.

Optional:  --oracle FIXTURE_DIR OUT.json  writes the Biopython cut-position
oracle (Enzyme.search, linear and circular) for every enzyme in the table on
every *.gb file in FIXTURE_DIR, keyed by file stem.

usage: python verify_enzymes.py [enzymes.lua] [--nvim PATH] [--oracle DIR OUT]
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

from Bio import Restriction as R
from Bio.Restriction import AllEnzymes

ALLOWED = {"name", "site", "cut5", "cut3", "dam", "dcm", "cpg"}
DEFAULT_NVIM = "/opt/nvim-linux-x86_64/bin/nvim"

DUMP_LUA = r"""
local path, out = _G.arg[1], _G.arg[2]
local t = dofile(path)
local rows = {}
for i, e in ipairs(t) do
  local keys = {}
  for k, _ in pairs(e) do keys[#keys + 1] = k end
  e = vim.deepcopy(e)
  e._keys = keys
  rows[i] = e
end
local f = assert(io.open(out, "w"))
f:write(vim.json.encode(rows))
f:close()
"""


def load_with_nvim(lua_path, nvim):
    with tempfile.TemporaryDirectory() as td:
        script = os.path.join(td, "dump.lua")
        out = os.path.join(td, "out.json")
        open(script, "w").write(DUMP_LUA)
        subprocess.run([nvim, "--clean", "-l", script, os.path.abspath(lua_path), out],
                       check=True, capture_output=True, timeout=60)
        rows = json.load(open(out))
    for r in rows:
        r["_keys"] = sorted(r["_keys"])
    return rows


ENTRY = re.compile(r"^\s*\{(.*)\},\s*$")
FIELD = re.compile(r'(\w+)\s*=\s*("([^"]*)"|-?\d+|true|false)')


def load_with_regex(lua_path):
    rows = []
    for line in open(lua_path):
        if line.lstrip().startswith("--"):
            continue
        m = ENTRY.match(line)
        if not m:
            continue
        e = {}
        for k, raw, s in FIELD.findall(m.group(1)):
            e[k] = s if raw.startswith('"') else (raw == "true" if raw in ("true", "false") else int(raw))
        e["_keys"] = sorted(e)
        rows.append(e)
    return rows


def biopython_enzyme(name):
    if name not in AllEnzymes.as_string():
        return None
    return getattr(R, name)


def check(rows):
    errors, unknown = [], []
    names = [r.get("name") for r in rows]
    if names != sorted(names):
        errors.append("table is not sorted by name")
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        errors.append(f"duplicate names: {sorted(dup)}")
    for r in rows:
        n = r.get("name")
        extra = set(r["_keys"]) - ALLOWED
        if extra:
            errors.append(f"{n}: unknown fields {sorted(extra)}")
        for flag in ("dam", "dcm", "cpg"):
            if flag in r and r[flag] is not True:
                errors.append(f"{n}: {flag} must be true or absent")
        e = biopython_enzyme(n)
        if e is None:
            unknown.append(n)
            continue
        cls = {c.__name__ for c in type(e).__mro__}
        if cls & {"TwoCuts", "NoCut", "Unknown", "NotDefined"}:
            errors.append(f"{n}: Biopython class {sorted(cls & {'TwoCuts','NoCut','Unknown','NotDefined'})} not allowed")
            continue
        exp = dict(site=e.site, cut5=e.fst5, cut3=e.size + e.fst3)
        for k, v in exp.items():
            if r.get(k) != v:
                errors.append(f"{n}: {k}={r.get(k)!r} but Biopython gives {v!r}")
    return errors, unknown


def oracle(rows, fixture_dir):
    from Bio import SeqIO
    out = {}
    for fn in sorted(os.listdir(fixture_dir)):
        if not fn.endswith((".gb", ".gbk", ".genbank")):
            continue
        rec = SeqIO.read(os.path.join(fixture_dir, fn), "genbank")
        seq = rec.seq.upper()
        stem = os.path.splitext(fn)[0]
        out[stem] = {"linear": {}, "circular": {}}
        for r in rows:
            e = getattr(R, r["name"])
            out[stem]["linear"][r["name"]] = sorted(e.search(seq, linear=True))
            out[stem]["circular"][r["name"]] = sorted(e.search(seq, linear=False))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lua", nargs="?", default="lua/plasmid/data/enzymes.lua")
    ap.add_argument("--nvim", default=os.environ.get("NVIM_BIN", DEFAULT_NVIM))
    ap.add_argument("--oracle", nargs=2, metavar=("FIXTURE_DIR", "OUT_JSON"))
    a = ap.parse_args()

    import Bio
    if a.nvim and (os.path.exists(a.nvim) or shutil.which(a.nvim)):
        rows, how = load_with_nvim(a.lua, a.nvim), "nvim -l"
        regex_rows = load_with_regex(a.lua)
        strip = lambda rs: [{k: v for k, v in r.items() if k != "_keys"} for r in rs]
        if strip(rows) != strip(regex_rows):
            print("WARNING: regex parser disagrees with nvim loader", file=sys.stderr)
    else:
        rows, how = load_with_regex(a.lua), "regex"
    print(f"loaded {len(rows)} entries from {a.lua} via {how}; Biopython {Bio.__version__}")

    errors, unknown = check(rows)
    for u in unknown:
        print(f"UNKNOWN to Biopython: {u}")
    for e in errors:
        print(f"FAIL {e}")
    nflags = {f: sum(1 for r in rows if r.get(f)) for f in ("dam", "dcm", "cpg")}
    print(f"methylation flags: {nflags}")

    if a.oracle:
        data = oracle(rows, a.oracle[0])
        with open(a.oracle[1], "w") as fh:
            json.dump(data, fh, indent=1, sort_keys=True)
        print(f"oracle written: {a.oracle[1]} ({', '.join(data)})")

    ok = not errors and not unknown
    print("PASS" if ok else f"FAILED: {len(errors)} mismatches, {len(unknown)} unknown")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
