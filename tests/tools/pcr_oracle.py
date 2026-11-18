"""Generate PCR expectations with pydna (development oracle only).

Random primer pairs are cut from the NCBI fixtures (forward = top strand,
reverse = reverse complement), with 15-25 annealing bases and optional random
5' tails, on linear and circular templates (products across the origin
included), plus pairs facing away from each other (no product). pydna's
pcr() gives the expected product sequence.

    python tests/tools/pcr_oracle.py > tests/oracle/expected/pcr.json
"""
import json
import random
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from pydna.amplify import pcr
from pydna.dseqrecord import Dseqrecord

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "ncbi"
rng = random.Random(5)


def rc(s):
    return str(Seq(s).reverse_complement())


def tail(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


cases = []
for stem in ["L09137", "J01749", "NC_001422"]:
    seq = str(SeqIO.read(FIX / f"{stem}.gb", "genbank").seq).upper()
    L = len(seq)
    for topo in ["linear", "circular"]:
        circ = topo == "circular"
        n = 0
        while n < 40:
            fl, rl = rng.randint(15, 25), rng.randint(15, 25)
            size = rng.randint(60, 2000)
            if circ:
                a = rng.randint(0, L - 1)
            else:
                a = rng.randint(0, L - size - 1)
            ext = seq + seq if circ else seq
            region = ext[a:a + size]
            fwd = region[:fl]
            rev = rc(region[-rl:])
            if rng.random() < 0.5:
                fwd = tail(rng.randint(3, 10)) + fwd
            if rng.random() < 0.5:
                rev = tail(rng.randint(3, 10)) + rev
            away = rng.random() < 0.15
            if away:
                # swap roles so the primers face away from each other
                fwd, rev = rc(region[:fl]), region[-rl:]
            tmpl = Dseqrecord(seq, circular=circ)
            try:
                prod = pcr(fwd, rev, tmpl, limit=13)
                product = str(prod.seq).upper()
            except Exception as e:  # no product or several
                msg = str(e)
                if "No" in msg or "no" in msg:
                    product = None
                else:
                    continue  # ambiguous (several products): skip the case
            cases.append({"file": stem, "topology": topo, "fwd": fwd, "rev": rev, "product": product,
                          "away": away})
            n += 1
json.dump({"source": "pydna.amplify.pcr (limit=13)", "cases": cases}, sys.stdout, indent=0)
