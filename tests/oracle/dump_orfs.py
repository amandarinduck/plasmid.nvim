"""Independent ORF oracle (development only).

Brute force, deliberately different from the plugin's frame scanner: from
every start codon on both strands, walk codon by codon (wrapping on circular
records) to the first stop; keep the most upstream start per stop. Codon
classes come from Biopython's NCBI tables.

    python tests/oracle/dump_orfs.py
"""
import json
import pathlib

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

ROOT = pathlib.Path(__file__).resolve().parents[1]
CASES = [
    ("ncbi/J01749.gb", "genbank"), ("ncbi/L09137.gb", "genbank"), ("ncbi/NC_001422.gb", "genbank"),
    ("ncbi/NC_005816.gb", "genbank"), ("ncbi/U01317.gb", "genbank"), ("snapgene/pFA-KanMX4.dna", "snapgene"),
]
SETTINGS = [(75, "atg", 11), (30, "atg", 11), (30, "alternative", 11), (50, "alternative", 4)]


def orfs(seq, circular, min_aa, starts, table):
    t = CodonTable.unambiguous_dna_by_id[table]
    stops = set(t.stop_codons)
    start_set = {"ATG"} if starts == "atg" else set(t.start_codons)
    L = len(seq)
    out = {}
    for strand, s in ((1, seq), (-1, str(Seq(seq).reverse_complement()))):
        def codon(p):  # 0-based position on this strand
            if circular:
                return "".join(s[(p + k) % L] for k in range(3))
            return s[p:p + 3] if p + 3 <= L else None
        best = {}  # stop position (0-based, strand coords, mod L) -> (length, start)
        for p in range(L if circular else L - 2):
            if codon(p) not in start_set:
                continue
            q, n = p, 0
            while True:
                c = codon(q)
                if c is None or n > L:  # ran off a linear end / no stop round the circle
                    q = None
                    break
                if c in stops:
                    break
                q, n = q + 3, n + 3
            if q is None:
                continue
            length = n + 3
            if length > L or n // 3 < min_aa:
                continue
            key = q % L if circular else q
            if key not in best or length > best[key][0]:
                best[key] = (length, p)
        for key, (length, p) in best.items():
            s1, e1 = p + 1, p + length  # 1-based on this strand (may exceed L when circular)
            if strand == 1:
                a, b = (s1 - 1) % L + 1, (e1 - 1) % L + 1
            else:
                a, b = (L - e1) % L + 1, (L - s1) % L + 1
            out[(strand, a, b)] = length // 3 - 1
    return sorted([[a, b, st, n] for (st, a, b), n in out.items()], key=lambda x: (x[0], -x[2]))


def main():
    res = {}
    for rel, fmt in CASES:
        rec = next(SeqIO.parse(ROOT / "fixtures" / rel, fmt))
        circ = rec.annotations.get("topology") == "circular"
        seq = str(rec.seq).upper()
        res[rel] = {f"{m}:{st}:{tb}": orfs(seq, circ, m, st, tb) for m, st, tb in SETTINGS}
        print(rel, len(seq), circ, {k: len(v) for k, v in res[rel].items()})
    (ROOT / "oracle" / "expected" / "orfs.json").write_text(json.dumps(res))


if __name__ == "__main__":
    main()
