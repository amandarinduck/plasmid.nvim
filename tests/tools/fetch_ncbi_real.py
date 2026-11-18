"""Fetch the real-world NCBI record corpus used by tests/spec/ncbi_records_spec.lua.

Development-only (Biopython Entrez). The plugin itself never runs Python.

    python tests/tools/fetch_ncbi_real.py OUTDIR

Writes OUTDIR/<file> for every entry in RECORDS plus OUTDIR/manifest.json.
Policy: tool=plasmid.nvim, no API key, no contact email, library-default
User-Agent, at most ~2.5 requests/second, one efetch per record (multi-record
files use one request with comma-separated ids).
"""
import json
import pathlib
import sys
import time
import warnings

from Bio import Entrez

Entrez.tool = "plasmid.nvim"
warnings.simplefilter("ignore")

N, P = "nuccore", "protein"
# (file, db, ids, rettype, extra efetch params, reason, commit)
RECORDS = [
    # circular plasmids / small circular genomes
    ("U55762.gbk", N, "U55762.1", "gbwithparts", {}, "pEGFP-N1: commercial circular vector (sparse 2003 annotation: 4 features)", True),
    ("X52328.gbk", N, "X52328.1", "gbwithparts", {}, "pBluescript II: classic circular cloning vector, EMBL-origin 1990 record", True),
    ("U13850.gbk", N, "U13850.1", "gbwithparts", {}, "pGEX-2T: circular GST expression vector, 1994 SYN record with 13 features", True),
    ("X06403.gbk", N, "X06403.1", "gbwithparts", {}, "pACYC184: pre-2000 circular plasmid", True),
    ("NC_001669.gbk", N, "NC_001669.1", "gbwithparts", {}, "SV40 RefSeq: circular, origin-spanning joins (join(...,1..N))", True),
    ("NC_012920.gbk", N, "NC_012920.1", "gbwithparts", {}, "human mtDNA rCRS: circular, D-loop spans origin, tRNA/rRNA, transl_except, overlapping genes", True),
    ("NC_000932.gbk", N, "NC_000932.1", "gbwithparts", {}, "Arabidopsis chloroplast: circular 154 kb, trans-spliced rps12, many complement(join())", True),
    ("NC_002483.gbk", N, "NC_002483.1", "gbwithparts", {}, "E. coli F plasmid: bacterial plasmid with many features", True),
    # phages / viruses
    ("NC_001416.gbk", N, "NC_001416.1", "gbwithparts", {}, "phage lambda: 48.5 kb phage genome", True),
    ("NC_045512.gbk", N, "NC_045512.2", "gbwithparts", {}, "SARS-CoV-2: ss-RNA, ribosomal-slippage join with overlapping parts, mat_peptide, long notes", True),
    ("NC_001802.gbk", N, "NC_001802.1", "gbwithparts", {}, "HIV-1: RNA virus, LTRs, join() CDSs", True),
    ("influenza_PR8.gbk", N, "NC_002016.1,NC_002017.1,NC_002018.1,NC_002019.1,NC_002020.1,NC_002021.1,NC_002022.1,NC_002023.1",
     "gbwithparts", {}, "influenza A/PR/8/34: 8-segment genome as one multi-record file (/segment qualifiers)", True),
    # mRNA / ncRNA
    ("NM_000546.gbk", N, "NM_000546.6", "gbwithparts", {}, "TP53 mRNA (RefSeq NM_): many misc_feature, variation, regulatory", True),
    ("XM_062871673.gbk", N, "XM_062871673.1", "gbwithparts", {}, "predicted model mRNA (XM_, Gnomon): mouse Srsx", True),
    ("NM_001267550.gbk", N, "NM_001267550.2", "gbwithparts", {}, "TTN mRNA: ~105 kb, /translation of ~35,000 aa (very long qualifier)", True),
    ("NM_000581.gbk", N, "NM_000581.4", "gbwithparts", {}, "GPX1 mRNA: selenocysteine /transl_except", True),
    ("NR_003286.gbk", N, "NR_003286.4", "gbwithparts", {}, "18S rRNA (NR_): non-coding RefSeq", True),
    # RefSeq genomic region extracts (seq_start/seq_stop, both strands)
    ("NC_000017_TP53_rev.gbk", N, "NC_000017.11", "gbwithparts", {"seq_start": 7668402, "seq_stop": 7687550, "strand": 2},
     "TP53 locus extract, strand=2 (gene on minus strand shown forward), spliced mRNA/CDS joins", True),
    ("NC_000011_HBB_fwd.gbk", N, "NC_000011.10", "gbwithparts", {"seq_start": 5225464, "seq_stop": 5229395, "strand": 1},
     "HBB locus extract, plus strand (gene on minus: complement(join()))", True),
    ("NC_000011_HBB_rev.gbk", N, "NC_000011.10", "gbwithparts", {"seq_start": 5225464, "seq_stop": 5229395, "strand": 2},
     "same HBB extract, strand=2", True),
    ("NC_000012_GAPDH_fwd.gbk", N, "NC_000012.12", "gbwithparts", {"seq_start": 6534000, "seq_stop": 6539000},
     "GAPDH locus extract, features truncated by the window get < / >", True),
    # CON / WGS / scaffold / segmented
    ("NZ_JCEXQL010000002_con.gbk", N, "NZ_JCEXQL010000002.1", "gb", {}, "RefSeq CON record as rettype=gb: CONTIG join(), no ORIGIN", True),
    ("NZ_JCEXQL010000002_parts.gbk", N, "NZ_JCEXQL010000002.1", "gbwithparts", {}, "same CON record expanded with gbwithparts", True),
    ("NW_028328278_con.gbk", N, "NW_028328278.1", "gb", {}, "RefSeq scaffold (NW_) as rettype=gb", True),
    ("CGRMEV000000000.gbk", N, "CGRMEV000000000.2", "gb", {}, "WGS master record: no sequence, WGS / WGS_SCAFLD lines", True),
    ("WGS_CONTIG", N, None, "gbwithparts", {}, "WGS contig of the master above", True),
    ("AH005336.gbk", N, "AH005336.3", "gbwithparts", {}, "segmented set (AH_): several source features, gap features, many secondary accessions", True),
    ("AH003512.gbk", N, "AH003512.2", "gbwithparts", {}, "segmented set with between-base misc_recomb sites (769^770)", True),
    # unusual location syntax / old records
    ("S75221.gbk", N, "S75221.1", "gbwithparts", {}, "between-base location 212^213 (misc_recomb), MMTV int-6-LTR fusion, 1990s record", True),
    ("AF084526.gbk", N, "AF084526.1", "gbwithparts", {}, "between-base misc_feature 2638^2639 (RNA-editing insertion sites), Physarum mitochondrial cox3", True),
    ("X61312.gbk", N, "X61312.1", "gbwithparts", {}, "order() gene location with > partial", True),
    ("J01636.gbk", N, "J01636.1", "gbwithparts", {}, "E. coli lac operon: 1980s record", True),
    ("U49845.gbk", N, "U49845.1", "gbwithparts", {}, "NCBI's own sample record: <1..>206 partial CDS, complement()", True),
    # proteins
    ("NP_000537.gp", P, "NP_000537.3", "gp", {}, "RefSeq protein p53: Region/Site features, order() sites", True),
    ("NP_005219.gp", P, "NP_005219.2", "gp", {}, "RefSeq protein EGFR: many Region/Site, long /note", True),
    ("AAB02572.gp", P, "AAB02572.1", "gp", {}, "GenPept (INSDC, non-RefSeq) protein: EGFP from cloning vector pEGFP-1, 1996 record", True),
    ("P01308.gp", P, "P01308.1", "gp", {}, "UniProtKB-derived GenPept (insulin): Bond features bond(a,b), long COMMENT", True),
    ("P04637.gp", P, "P04637.2", "gp", {}, "UniProtKB-derived GenPept (P53_HUMAN): hundreds of Site/Region features", True),
    ("4HHB_A.gp", P, "4HHB_A", "gp", {}, "PDB-derived GenPept: SecStr and Het features (bond())", True),
    # big ones: measured, not committed
    ("NC_000908.gbk", N, "NC_000908.2", "gbwithparts", {}, "Mycoplasma genitalium chromosome 580 kb (scaling point)", False),
    ("NC_000913.gbk", N, "NC_000913.3", "gbwithparts", {}, "E. coli K-12 MG1655 chromosome 4.6 Mb (big record)", False),
]


def fetch(db, ids, rettype, extra):
    h = Entrez.efetch(db=db, id=ids, rettype=rettype, retmode="text", **extra)
    text = h.read()
    h.close()
    time.sleep(0.4)
    return text


def main():
    out = pathlib.Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for fname, db, ids, rettype, extra, reason, commit in RECORDS:
        if fname == "WGS_CONTIG":
            master = (out / "CGRMEV000000000.gbk").read_text()
            first = next(l.split()[1].split("-")[0] for l in master.splitlines() if l.startswith("WGS "))
            ids, fname = first, first.split(".")[0] + ".gbk"
        target = out / fname
        if not target.exists() or target.stat().st_size == 0:
            t0 = time.time()
            text = fetch(db, ids, rettype, extra)
            target.write_text(text)
            print(f"{fname}: {len(text)} bytes in {time.time() - t0:.1f}s", flush=True)
        manifest.append({"file": fname, "db": db, "ids": ids, "rettype": rettype, "params": extra,
                         "reason": reason, "commit": commit, "bytes": target.stat().st_size})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
