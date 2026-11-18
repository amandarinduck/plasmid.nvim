"""Generate oracle JSON from Biopython for every GenBank/FASTA fixture.

Development-only: the plugin never depends on Python. The JSON written here is
committed to tests/oracle/expected/ and the Lua test-suite compares against it.

    python tests/oracle/dump_biopython.py
"""
import hashlib
import warnings
import json
import pathlib

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIX = ROOT / "fixtures"
OUT = ROOT / "oracle" / "expected"
OUT.mkdir(parents=True, exist_ok=True)


def parts(loc):
    if loc is None:  # Biopython could not parse it (e.g. legacy 102.110)
        return None
    ps = loc.parts if isinstance(loc, CompoundLocation) else [loc]
    out = []
    for p in ps:
        # 1-based inclusive; a between-site (23^24) has start1 = stop1 + 1
        out.append([int(p.start) + 1, int(p.end), p.strand if p.strand is not None else 1, p.ref])
    return out


def quals(f):
    out = []
    for k, vs in f.qualifiers.items():
        for v in vs:
            out.append([k, str(v)])
    return sorted(out)


def dump(path, fmt):
    recs = []
    for r in SeqIO.parse(path, fmt):
        seq = str(r.seq)
        d = {
            "name": r.name,
            "id": r.id,
            "length": len(seq),
            "sha256": hashlib.sha256(seq.encode()).hexdigest(),
        }
        if fmt in ("genbank", "snapgene"):
            d["topology"] = r.annotations.get("topology", "linear")
            d["description"] = r.description
            d["features"] = [
                {"type": f.type, "parts": parts(f.location), "qualifiers": quals(f)}
                for f in r.features
            ]
        else:
            d["description"] = r.description
        recs.append(d)
    return recs


def main():
    warnings.simplefilter("ignore")
    n = 0
    for p in sorted(FIX.rglob("*")):
        if p.suffix in (".gb", ".gbk"):
            fmt = "genbank"
        elif p.suffix in (".fa", ".fasta", ".faa"):
            fmt = "fasta"
        elif p.suffix == ".dna":
            fmt = "snapgene"
        else:
            continue
        try:
            recs = dump(p, fmt)
        except Exception as e:  # fixture Biopython can't read: no oracle
            print(f"skip {p.name}: {e}")
            continue
        rel = p.relative_to(FIX).as_posix().replace("/", "__")
        (OUT / f"{rel}.json").write_text(json.dumps({"format": fmt, "records": recs}, indent=1))
        n += 1
    print(f"wrote {n} oracle files")


if __name__ == "__main__":
    main()
