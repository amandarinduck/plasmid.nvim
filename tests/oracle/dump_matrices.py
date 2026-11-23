"""Write lua/plasmid/data/blosum62.lua from Biopython's copy of the NCBI matrix.

Development-only: the plugin never depends on Python. Run from the repo root:

    python tests/oracle/dump_matrices.py
"""
import pathlib

from Bio.Align import substitution_matrices

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "lua" / "plasmid" / "data" / "blosum62.lua"


def main():
    m = substitution_matrices.load("BLOSUM62")
    alphabet = "".join(m.alphabet)
    rows = [
        '  ["{}"] = {{ {} }},'.format(a, ",".join(str(int(m[a, b])) for b in m.alphabet))
        for a in m.alphabet
    ]
    OUT.write_text(
        "-- BLOSUM62 substitution matrix, generated from Biopython's copy of the\n"
        "-- NCBI matrix by tests/oracle/dump_matrices.py. Values are integers; the\n"
        "-- column order is M.alphabet.\n"
        "local M = {}\n\n"
        'M.alphabet = "%s"\n\n' % alphabet
        + "M.scores = {\n"
        + "\n".join(rows)
        + "\n}\n\nreturn M\n"
    )
    print("wrote", OUT.relative_to(ROOT), len(m.alphabet), "symbols")


if __name__ == "__main__":
    main()
