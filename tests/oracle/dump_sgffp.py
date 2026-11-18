"""Oracle for SnapGene protein/RNA files, which Biopython does not read.

Uses sgffp (MIT, https://github.com/merv1n34k/sgffp) in a development venv:
    pip install sgffp && python tests/oracle/dump_sgffp.py
"""
import hashlib
import json
import pathlib

from sgffp import SgffReader

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "oracle" / "expected"
ALPHA = {1: "dna", 2: "protein", 7: "rna"}

for p in sorted((ROOT / "fixtures").rglob("*")):
    if p.suffix not in (".prot", ".rna"):
        continue
    o = SgffReader.from_file(str(p))
    seq = o.sequence.value
    rec = {
        "length": len(seq),
        "sha256": hashlib.sha256(seq.encode()).hexdigest(),
        "topology": o.sequence.topology,
        "alphabet": ALPHA[o.cookie.type_of_sequence],
        "nfeatures": len(o.features or []),
    }
    rel = p.relative_to(ROOT / "fixtures").as_posix().replace("/", "__")
    (OUT / f"{rel}.json").write_text(json.dumps({"format": "snapgene_seq", "records": [rec]}, indent=1))
    print("wrote", rel)
