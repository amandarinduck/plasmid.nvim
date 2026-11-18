#!/usr/bin/env python3
"""Biopython check of GenBank files written after sequence edits (development tool).

For each fixture and seed, the plugin applies random edits and writes GenBank
(tests/tools/edit_dump.lua). Biopython then parses the written file; its
sequence and every feature's extracted bases must equal what the plugin holds.

    python tests/tools/edit_oracle.py [seeds]      # needs biopython
"""
import json, pathlib, subprocess, sys, tempfile
from Bio import SeqIO

ROOT = pathlib.Path(__file__).resolve().parents[2]
FIXTURES = ["ncbi/J01749.gb", "ncbi/NC_001422.gb", "ncbi/U01317.gb", "ncbi/NC_005816.gb"]

def main():
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    files = feats = 0
    with tempfile.TemporaryDirectory() as tmp:
        for rel in FIXTURES:
            for seed in range(1, seeds + 1):
                out, js = f"{tmp}/out.gb", f"{tmp}/out.json"
                subprocess.run(["nvim", "-l", str(ROOT / "tests/tools/edit_dump.lua"), str(ROOT / "tests/fixtures" / rel),
                                out, js, str(seed), "25"], check=True, capture_output=True)
                plugin = json.load(open(js))
                rec = SeqIO.read(out, "genbank")
                assert str(rec.seq).upper() == plugin["seq"].upper(), f"{rel} seed {seed}: sequence differs"
                assert len(rec.features) == len(plugin["features"]), f"{rel} seed {seed}: feature count"
                for i, (bf, pf) in enumerate(zip(rec.features, plugin["features"])):
                    if pf["residues"] is None:
                        continue
                    got = str(bf.extract(rec.seq)).upper()
                    if got != pf["residues"].upper():
                        raise AssertionError(f"{rel} seed {seed}: feature {i} {pf['type']} {pf['location']}: "
                                             f"Biopython reads {bf.location}")
                    feats += 1
                files += 1
    print(f"{files} edited GenBank files, {feats} features: Biopython agrees on sequence and every feature's bases")

if __name__ == "__main__":
    main()
