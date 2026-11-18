"""Biopython oracle for the flat files recorded in tests/fixtures/ncbi_http/.

Development-only. For every recorded efetch response that is a GenBank/GenPept
flat file (body starting with LOCUS), write
tests/oracle/expected/ncbi_http__<case>.json in the same layout as
dump_biopython.py; tests/spec/ncbi_spec.lua compares the plugin parser with it.

    python tests/oracle/dump_ncbi_http.py
"""
import json
import pathlib
import sys
import warnings

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from dump_biopython import dump  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
HTTP = ROOT / "fixtures" / "ncbi_http"
OUT = ROOT / "oracle" / "expected"


def main():
    warnings.simplefilter("ignore")
    manifest = json.loads((HTTP / "manifest.json").read_text())
    n = 0
    for case in manifest["cases"]:
        for r in case["requests"]:
            if r["status"] != 200 or "efetch.fcgi" not in r["url"]:
                continue
            body = HTTP / f"{r['key']}.body"
            if not body.read_bytes().lstrip().startswith(b"LOCUS"):
                continue
            try:
                recs = dump(body, "genbank")
            except Exception as e:  # e.g. rettype=gb CON record: no sequence
                print(f"skip {case['name']}: {type(e).__name__}: {e}")
                continue
            (OUT / f"ncbi_http__{case['name']}.json").write_text(
                json.dumps({"format": "genbank", "key": r["key"], "records": recs}, indent=1))
            n += 1
    print(f"wrote {n} oracle files")


if __name__ == "__main__":
    main()
