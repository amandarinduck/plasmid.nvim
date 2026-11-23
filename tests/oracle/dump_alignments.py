"""Generate the pairwise-alignment oracle from Biopython's PairwiseAligner.

Development-only: the plugin never depends on Python. The JSON written here is
committed to tests/oracle/expected/alignments.json and tests/spec/
align_oracle_spec.lua replays it offline.

    conda activate primerdev && python tests/oracle/dump_alignments.py

Parameter mapping (verified in docs/alignment_report.md, "Mapping"):

    engine                          Biopython PairwiseAligner
    ------------------------------  ---------------------------------------
    match / mismatch                match_score / mismatch_score
    gap_open + gap_extend           open_gap_score
    gap_extend                      extend_gap_score
    mode = "global"                 mode = "global"
    mode = "local"                  mode = "local"
    mode = "semi" (free ref ends)   mode = "global", end_deletion_score = 0
    mode = "semi", free_query_ends  mode = "global", end_insertion_score = 0
    protein = true                  substitution_matrix = BLOSUM62

Biopython calls a gap in the *query* a deletion and a gap in the *target* an
insertion, so freeing end deletions is what lets the reference hang off the
ends -- which is the engine's `semi` with free_ref_ends.

Alignment representation.  Biopython's global mode always spells the free end
gaps out as columns; the engine omits them and reports ref_start/ref_stop
instead.  The cigar recorded here is therefore trimmed of the end gaps that
the mode makes free, and ref_start/query_start are moved accordingly, so the
committed form is directly comparable with A.cigar(aln).
"""
import hashlib
import json
import pathlib
import random

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIX = ROOT / "fixtures"
OUT = ROOT / "oracle" / "expected" / "alignments.json"

BLOSUM62 = substitution_matrices.load("BLOSUM62")

# --- scoring schemes ------------------------------------------------------
# Nucleotide.  gap_open is charged once per gap *in addition to* gap_extend
# per base, so a one-base gap costs gap_open + gap_extend.
NUC_SCHEMES = {
    "blastn": dict(match=2, mismatch=-3, gap_open=-5, gap_extend=-2),
    "unit": dict(match=1, mismatch=-1, gap_open=-2, gap_extend=-1),
    "longgap": dict(match=5, mismatch=-4, gap_open=-10, gap_extend=-1),
    "linear": dict(match=2, mismatch=-1, gap_open=0, gap_extend=-2),
    "steep": dict(match=3, mismatch=-2, gap_open=-7, gap_extend=-3),
}
# Protein, BLOSUM62.  "blastp" is NCBI's 11/1 written in the engine's
# convention (a one-residue gap costs 11, each further residue 1).
AA_SCHEMES = {
    "blastp": dict(gap_open=-10, gap_extend=-1),
    "cheap": dict(gap_open=-5, gap_extend=-2),
    "near_linear": dict(gap_open=-1, gap_extend=-1),
}

# mode name -> (biopython mode, extra PairwiseAligner attributes, engine opts)
MODES = {
    "global": ("global", {}, {"mode": "global"}),
    "local": ("local", {}, {"mode": "local"}),
    "semi": ("global", {"end_deletion_score": 0.0},
             {"mode": "semi", "free_ref_ends": True, "free_query_ends": False}),
    "semi_q": ("global", {"end_insertion_score": 0.0},
               {"mode": "semi", "free_ref_ends": False, "free_query_ends": True}),
    "overlap": ("global", {"end_gap_score": 0.0},
                {"mode": "semi", "free_ref_ends": True, "free_query_ends": True}),
}

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def aligner(mode, scheme, protein):
    bio_mode, extra, _ = MODES[mode]
    a = Align.PairwiseAligner()
    a.mode = bio_mode
    if protein:
        a.substitution_matrix = BLOSUM62
    else:
        a.match_score = scheme["match"]
        a.mismatch_score = scheme["mismatch"]
    a.open_gap_score = scheme["gap_open"] + scheme["gap_extend"]
    a.extend_gap_score = scheme["gap_extend"]
    for k, v in extra.items():
        setattr(a, k, v)
    return a


def to_ops(alignment):
    """[(op, n), ...] from a Biopython alignment's coordinates."""
    co = alignment.coordinates
    ops = []
    for c in range(co.shape[1] - 1):
        dt = int(co[0][c + 1] - co[0][c])
        dq = int(co[1][c + 1] - co[1][c])
        if dt > 0 and dq > 0:
            assert dt == dq, (dt, dq)
            op, n = "M", dt
        elif dt > 0:
            op, n = "D", dt
        elif dq > 0:
            op, n = "I", dq
        else:
            continue
        if ops and ops[-1][0] == op:
            ops[-1][1] += n
        else:
            ops.append([op, n])
    return ops, int(co[0][0]), int(co[1][0])


def trim(ops, ref_start0, query_start0, mode):
    """Drop the end gaps this mode makes free, moving the starts with them.

    Returns (ops, ref_start, query_start) 1-based, in the engine's convention.
    """
    # "local" never emits end gaps; listed for completeness only.
    free = {"local": "", "semi": "D", "semi_q": "I", "overlap": "DI",
            "global": ""}[mode]
    ops = [list(o) for o in ops]
    while ops and ops[0][0] in free:
        op, n = ops.pop(0)
        if op == "D":
            ref_start0 += n
        else:
            query_start0 += n
    while ops and ops[-1][0] in free:
        ops.pop()
    return ops, ref_start0 + 1, query_start0 + 1


def cigar(ops):
    return "".join("%d%s" % (n, op) for op, n in ops)


def count_optimal(a, ref, query, limit=10000):
    try:
        alns = a.align(ref, query)
        n = len(alns)
        return int(n) if n <= limit else None
    except (OverflowError, MemoryError, ValueError):
        return None


def oracle(ref, query, mode, scheme, protein):
    """Biopython's answer for one pair, in the engine's convention."""
    if not ref or not query:
        return {"score": None, "why": "empty sequence"}
    a = aligner(mode, scheme, protein)
    score = a.score(ref, query)
    try:
        alns = a.align(ref, query)
        best = alns[0]
    except (IndexError, ValueError):
        # local mode with no positive-scoring segment
        return {"score": float(score), "alignments": 0}
    ops, r0, q0 = to_ops(best)
    ops, r1, q1 = trim(ops, r0, q0, mode)
    return {
        "score": float(score),
        "cigar": cigar(ops),
        "ref_start": r1,
        "query_start": q1,
        "n_optimal": count_optimal(a, ref, query),
    }


# --- sequence material ----------------------------------------------------

def load_fixtures():
    seqs = {}
    for name in ("AB000263", "J01749", "L09137", "M81224", "NC_001422",
                 "NC_005816", "X55053"):
        rec = next(SeqIO.parse(FIX / "ncbi" / f"{name}.gb", "genbank"))
        seqs[name] = str(rec.seq).upper()
    prots = {}
    for rec in SeqIO.parse(FIX / "protein" / "proteins.faa", "fasta"):
        prots[rec.id.split("|")[-1] if "|" in rec.id else rec.id] = str(rec.seq)
    return seqs, prots


def mutate(rng, s, sub_rate=0.0, ins_rate=0.0, del_rate=0.0, alphabet="ACGT"):
    out = []
    for ch in s:
        r = rng.random()
        if r < del_rate:
            continue
        if r < del_rate + ins_rate:
            out.append(rng.choice(alphabet))
            out.append(ch)
            continue
        if r < del_rate + ins_rate + sub_rate:
            out.append(rng.choice([c for c in alphabet if c != ch]))
            continue
        out.append(ch)
    return "".join(out)


def rand_dna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def build_cases():
    rng = random.Random(20240817)
    seqs, prots = load_fixtures()
    fixture_names = sorted(seqs)
    prot_names = sorted(prots)
    cases = []

    def add(family, mode, scheme_name, ref, query, protein=False, strands=False,
            note=None):
        scheme = (AA_SCHEMES if protein else NUC_SCHEMES)[scheme_name]
        c = {
            "id": "c%04d" % (len(cases) + 1),
            "family": family,
            "mode": mode,
            "scheme": scheme_name,
            "protein": protein,
            "ref": ref,
            "query": query,
        }
        if note:
            c["note"] = note
        c.update(oracle(ref, query, mode, scheme, protein))
        if strands:
            rc = revcomp(query)
            o = oracle(ref, rc, mode, scheme, protein)
            c["rc_score"] = o["score"]
        cases.append(c)

    nuc_modes = list(MODES)
    nuc_schemes = list(NUC_SCHEMES)

    # 1. random pairs: unrelated, and one derived from the other
    for i in range(4):
        for mode in nuc_modes:
            for sname in nuc_schemes:
                n = rng.randint(40, 160)
                ref = rand_dna(rng, n)
                add("random_unrelated", mode, sname, ref, rand_dna(rng, rng.randint(40, 160)))
                q = mutate(rng, ref, sub_rate=0.08, ins_rate=0.02, del_rate=0.02)
                add("random_related", mode, sname, ref, q)

    # 2. fixture-derived pairs with substitutions / indels
    for name in fixture_names:
        s = seqs[name]
        for mode in nuc_modes:
            for sname in nuc_schemes:
                n = rng.randint(120, 320)
                start = rng.randrange(0, len(s) - n)
                ref = s[start:start + n]
                rates = rng.choice([(0.01, 0.0, 0.0), (0.05, 0.0, 0.0),
                                    (0.03, 0.01, 0.01), (0.10, 0.02, 0.02)])
                q = mutate(rng, ref, *rates)
                add("fixture_mutated", mode, sname, ref, q,
                    note="%s:%d-%d sub/ins/del %.2f/%.2f/%.2f" % (
                        name, start + 1, start + n, *rates))

    # 3. reverse-complemented queries: the forward score and the rc score, so
    #    the spec can check align_strands picks the better one
    for name in fixture_names:
        s = seqs[name]
        for mode in ("local", "semi"):
            for sname in ("blastn", "unit", "longgap"):
                n = rng.randint(120, 260)
                start = rng.randrange(0, len(s) - n)
                ref = s[start:start + n]
                inner = ref[rng.randint(5, 20):n - rng.randint(5, 20)]
                q = revcomp(mutate(rng, inner, sub_rate=0.04, del_rate=0.01, ins_rate=0.01))
                add("revcomp", mode, sname, ref, q, strands=True, note=name)

    # 4. queries that hang off an end of the reference
    for name in fixture_names:
        s = seqs[name]
        for side in ("left", "right", "both"):
            for sname in ("blastn", "unit", "steep", "linear"):
                for mode in ("semi", "semi_q", "global", "local"):
                    n = rng.randint(120, 240)
                    start = rng.randrange(20, len(s) - n - 20)
                    ref = s[start:start + n]
                    core = ref[:rng.randint(60, 100)] if side in ("right", "both") \
                        else ref[-rng.randint(60, 100):]
                    tail = rand_dna(rng, rng.randint(8, 30))
                    if side == "left":
                        q = tail + core
                    elif side == "right":
                        q = core + tail
                    else:
                        q = tail + ref[20:80] + rand_dna(rng, rng.randint(8, 30))
                    add("hangoff_" + side, mode, sname, ref,
                        mutate(rng, q, sub_rate=0.02))

    # 5. queries much shorter than the reference
    for name in fixture_names:
        s = seqs[name]
        for sname in ("blastn", "unit", "longgap", "steep"):
            for mode in ("semi", "local"):
                w = min(400, len(s))
                ref = s[rng.randrange(0, len(s) - w + 1):][:w]
                k = rng.randint(14, 40)
                at = rng.randrange(0, w - k)
                add("short_query", mode, sname, ref,
                    mutate(rng, ref[at:at + k], sub_rate=0.03), note="%s %d bp" % (name, k))

    # 6. protein pairs, BLOSUM62
    pairs = [(a, b) for i, a in enumerate(prot_names) for b in prot_names[i + 1:]]
    for a, b in pairs:
        for mode in nuc_modes:
            for sname in AA_SCHEMES:
                add("protein_pair", mode, sname, prots[a], prots[b], protein=True,
                    note="%s vs %s" % (a, b))
    for name in prot_names:
        p = prots[name]
        for mode in nuc_modes:
            for sname in AA_SCHEMES:
                q = mutate(rng, p, sub_rate=0.08, ins_rate=0.01, del_rate=0.01,
                           alphabet="ACDEFGHIKLMNPQRSTVWY")
                add("protein_mutated", mode, sname, p, q, protein=True, note=name)
                add("protein_fragment", mode, sname, p,
                    p[rng.randint(10, 30):len(p) - rng.randint(10, 30)],
                    protein=True, note=name + " fragment")

    # 7. pairs that must not align
    for sname in nuc_schemes:
        for mode in nuc_modes:
            add("disjoint", mode, sname, "A" * 60, "C" * 60)
            add("disjoint", mode, sname, "AC" * 40, "GT" * 40)
            add("empty_query", mode, sname, seqs["X55053"][:200], "")
    for sname in AA_SCHEMES:
        add("disjoint", "local", sname, "W" * 40, "P" * 40, protein=True)

    return cases


def main():
    cases = build_cases()
    doc = {
        "biopython": __import__("Bio").__version__,
        "generator": "tests/oracle/dump_alignments.py",
        "nuc_schemes": NUC_SCHEMES,
        "aa_schemes": AA_SCHEMES,
        "modes": {k: {"biopython_mode": v[0], "biopython_attrs": v[1],
                      "engine_opts": v[2]} for k, v in MODES.items()},
        "cases": cases,
    }
    text = json.dumps(doc, separators=(",", ":"), sort_keys=False)
    OUT.write_text(text)
    fams = {}
    for c in cases:
        fams[c["family"]] = fams.get(c["family"], 0) + 1
    print("wrote %s: %d cases, %.0f kB" % (OUT.name, len(cases), len(text) / 1024))
    for k in sorted(fams):
        print("  %-18s %d" % (k, fams[k]))
    print("sha256", hashlib.sha256(text.encode()).hexdigest()[:16])


if __name__ == "__main__":
    main()
