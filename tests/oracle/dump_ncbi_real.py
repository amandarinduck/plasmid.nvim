"""Biopython oracle for the real NCBI corpus in tests/fixtures/ncbi_real/.

Development-only. Writes tests/oracle/expected/ncbi_real.json, keyed by fixture
file name, which tests/spec/ncbi_records_spec.lua compares against.

    python tests/oracle/dump_ncbi_real.py                      # committed fixtures
    python tests/oracle/dump_ncbi_real.py --out X.json F1 F2   # any files (e.g. big records)

Records without residues (CON as rettype=gb, WGS masters) have sha256 = null and
"defined" = false: Biopython knows their length but not their sequence.
"""
import argparse
import hashlib
import json
import pathlib
import sys
import time
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import UndefinedSequenceError

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from dump_biopython import parts, quals  # noqa: E402  (same location/qualifier model)

# Codons with IUPAC ambiguity codes that still translate to one amino acid
# (or stop), as NCBI does; X otherwise. Found via AH003512's GCN (Ala).
CODONS = ["GCN", "GGN", "CTN", "TCN", "CCN", "ACN", "GTN", "CGN", "AAR", "GAY", "TAR", "TRA",
          "YTR", "MGR", "ATH", "NNN", "GAN"]

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIX = ROOT / "fixtures" / "ncbi_real"
OUT = ROOT / "oracle" / "expected" / "ncbi_real.json"


def dump_file(path):
    recs = []
    t0 = time.perf_counter()
    parsed = list(SeqIO.parse(path, "genbank"))
    secs = time.perf_counter() - t0
    for r in parsed:
        try:
            seq = str(r.seq)
            sha, defined = hashlib.sha256(seq.upper().encode()).hexdigest(), True
        except UndefinedSequenceError:
            sha, defined = None, False
        a = r.annotations
        recs.append({
            "name": r.name,
            "id": r.id,
            "accessions": a.get("accessions", []),
            "description": r.description,
            "length": len(r.seq),
            "defined": defined,
            "sha256": sha,
            "topology": a.get("topology", "linear"),
            "molecule_type": a.get("molecule_type"),
            "data_file_division": a.get("data_file_division"),
            "contig": a.get("contig"),
            "features": [
                {"type": f.type, "location": str(f.location) if f.location is not None else None,
                 "parts": parts(f.location), "qualifiers": quals(f)}
                for f in r.features
            ],
        })
    return {"records": recs, "biopython_parse_ms": round(secs * 1000, 1)}


def key(p):
    return p.relative_to(FIX).as_posix() if FIX in p.resolve().parents else p.name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("files", nargs="*")
    args = ap.parse_args()
    warnings.simplefilter("ignore")
    files = [pathlib.Path(f) for f in args.files] or sorted(
        p for p in FIX.rglob("*") if p.suffix in (".gbk", ".gp"))
    res = {}
    for p in files:
        try:
            res[key(p)] = dump_file(p)
        except Exception as e:  # recorded, never silently skipped
            res[key(p)] = {"error": f"{type(e).__name__}: {e}"}
        print(key(p), res[key(p)].get("error") or len(res[key(p)]["records"]), flush=True)
    if not args.files:
        res["_codons"] = {c: str(Seq(c).translate()) for c in CODONS}
    pathlib.Path(args.out).write_text(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
