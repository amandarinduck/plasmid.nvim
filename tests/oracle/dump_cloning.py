"""Generate the cloning oracle: pydna's answer for several hundred assemblies.

Development-only. The plugin never depends on Python: the JSON written here is
committed to tests/oracle/expected/cloning.json and tests/spec/
cloning_oracle_spec.lua compares the Lua engine against it offline.

    conda activate primerdev        # pydna 5.5.16, biopython 1.88
    python tests/oracle/dump_cloning.py

Three families of case:

  ligation    fragments of the test plasmids (whole records, or linear slices
              of them) cut with one or more enzymes and religated, plus
              vector+insert cases that take named fragments from two digests.
              Oracle: pydna.assembly2.ligation_assembly.
  gibson      2-5 linear parts carrying terminal homology arms.
              Oracle: pydna.assembly2.gibson_assembly.
  goldengate  2-6 parts flanked by inward-facing Type IIS sites.
              Oracle: pydna.assembly2.golden_gate_assembly.

Products are compared through a canonical form so that rotation and strand do
not matter: for a circular product the lexicographically least rotation of the
sequence and of its reverse complement, whichever is smaller; for a linear one
min(seq, revcomp(seq)). Only the sha256 of that string and its length are
committed, which keeps the file small.

A handful of cases cannot be expressed with pydna (a Golden Gate part with no
site at all, or with mixed overhang lengths: pydna simply finds no join and
says nothing, where the plugin is expected to explain why). Those carry
"oracle": "hand" and a reason, and their expectation is derived from the
model, not from pydna.
"""

import hashlib
import json
import pathlib
import random
import warnings

from Bio import SeqIO
from Bio import Restriction
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import ligation_assembly, gibson_assembly, golden_gate_assembly

warnings.simplefilter("ignore")

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIX = ROOT / "fixtures"
OUT = ROOT / "oracle" / "expected" / "cloning.json"

SEED = 20260101
RNG = random.Random(SEED)

COMP = str.maketrans("ACGTacgtNnRrYyKkMmSsWwBbVvDdHh", "TGCAtgcaNnYyRrMmKkSsWwVvBbHhDd")


def rc(s):
    return s.translate(COMP)[::-1]


def least_rotation(s):
    """Booth's algorithm: the lexicographically least rotation of s."""
    n = len(s)
    t = s + s
    f = [-1] * (2 * n)
    k = 0
    for j in range(1, 2 * n):
        sj = t[j]
        i = f[j - k - 1]
        while i != -1 and sj != t[k + i + 1]:
            if sj < t[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != t[k + i + 1]:
            if sj < t[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return t[k:k + n]


def canonical(seq, circular):
    s = seq.upper()
    if circular:
        a, b = least_rotation(s), least_rotation(rc(s))
    else:
        a, b = s, rc(s)
    return a if a <= b else b


def product(rec):
    s = str(rec.seq).upper()
    return {"len": len(s), "circular": bool(rec.circular),
            "canon": hashlib.sha256(canonical(s, rec.circular).encode()).hexdigest()}


def expectation(records):
    """pydna's products, split by topology and deduplicated by canonical form."""
    circ, lin, seen = [], [], set()
    for r in records:
        p = product(r)
        key = (p["circular"], p["canon"])
        if key in seen:
            continue
        seen.add(key)
        (circ if p["circular"] else lin).append({"len": p["len"], "canon": p["canon"]})
    circ.sort(key=lambda p: (p["len"], p["canon"]))
    lin.sort(key=lambda p: (p["len"], p["canon"]))
    return {"circular": circ, "linear": lin}


# ---------------------------------------------------------------------------
# Random sequence that is clean of the enzymes we are about to use
# ---------------------------------------------------------------------------

def rnd(n, gc=0.5):
    at, gcb = "AT", "GC"
    return "".join(RNG.choice(gcb) if RNG.random() < gc else RNG.choice(at)
                   for _ in range(n))


def clean(n, forbidden, gc=0.5):
    """Random sequence containing none of `forbidden` (or their rc)."""
    bad = set()
    for f in forbidden:
        bad.add(f.upper())
        bad.add(rc(f.upper()))
    for _ in range(400):
        s = rnd(n, gc)
        if not any(b in s for b in bad):
            return s
    raise RuntimeError("could not draw a clean sequence")


# ---------------------------------------------------------------------------
# Templates
# ---------------------------------------------------------------------------

TEMPLATES = {
    "pUC19": "ncbi/L09137.gb",
    "pBR322": "ncbi/J01749.gb",
    "phiX174": "ncbi/NC_001422.gb",
    "pGRAMMAR": "synthetic/grammar.gb",
}

_seqs = {}


def tseq(key):
    if key not in _seqs:
        r = SeqIO.read(FIX / TEMPLATES[key], "genbank")
        _seqs[key] = str(r.seq).upper()
    return _seqs[key]


def enz(names):
    return [getattr(Restriction, n) for n in names]


cases = []
notes = []
skipped = []


def safe(fn, *a, **kw):
    """pydna refuses to enumerate above its own path/assembly caps. When it
    does, the oracle has no answer for that slice of the case: record None and
    the comparison skips it."""
    try:
        return fn(*a, **kw)
    except ValueError:
        return None


def add(case):
    cases.append(case)
    return case


# ---------------------------------------------------------------------------
# 1. Restriction - ligation
# ---------------------------------------------------------------------------

# enzyme sets, grouped by the kind of end they leave
FIVE = ["EcoRI", "HindIII", "BamHI", "XhoI", "SalI", "XbaI", "NcoI", "NdeI",
        "SpeI", "AvrII", "BglII", "NheI", "MfeI", "AflII"]
THREE = ["PstI", "KpnI", "SacI", "SphI", "ApaI", "BsiWI", "NsiI", "AatII", "BsrGI"]
BLUNT = ["EcoRV", "SmaI", "ScaI", "PvuII", "HincII", "DraI", "StuI", "NruI", "SspI"]

SETS = []
for n in FIVE:
    SETS.append([n])
for n in THREE:
    SETS.append([n])
for n in BLUNT:
    SETS.append([n])
# pairs and triples, mixed end chemistry
PAIRS = [
    ["EcoRI", "HindIII"], ["EcoRI", "BamHI"], ["BamHI", "SalI"], ["XbaI", "SpeI"],
    ["PstI", "SacI"], ["KpnI", "SacI"], ["PstI", "EcoRI"], ["KpnI", "BamHI"],
    ["EcoRV", "SmaI"], ["ScaI", "PvuII"], ["EcoRV", "EcoRI"], ["ScaI", "PstI"],
    ["SphI", "SalI"], ["NdeI", "XhoI"], ["AatII", "NcoI"], ["DraI", "HindIII"],
    ["NruI", "ApaI"], ["BsiWI", "SpeI"], ["SspI", "BglII"], ["StuI", "NsiI"],
]
TRIPLES = [
    ["EcoRI", "HindIII", "PstI"], ["EcoRI", "BamHI", "SalI"],
    ["PstI", "KpnI", "SacI"], ["EcoRV", "SmaI", "ScaI"],
    ["EcoRI", "EcoRV", "PstI"], ["XbaI", "SpeI", "NheI"],
    ["HindIII", "ScaI", "KpnI"], ["BamHI", "PvuII", "SphI"],
    ["NdeI", "NcoI", "XhoI"], ["AflII", "DraI", "ApaI"],
]
SETS += PAIRS + TRIPLES

MAX_FRAGS = 6
lig_id = 0


def ligation_case(label, frag_specs, dseq_frags, extra=None):
    """frag_specs: what the Lua side must build. dseq_frags: pydna fragments."""
    global lig_id
    circ = safe(ligation_assembly, dseq_frags, allow_blunt=True, circular_only=True)
    if circ is None:
        skipped.append((label, "pydna path limit"))
        return None
    lin = safe(ligation_assembly, dseq_frags, allow_blunt=True, circular_only=False)
    lin = None if lin is None else [r for r in lin if not r.circular]
    lig_id += 1
    c = {
        "id": "lig%04d" % lig_id,
        "method": "ligation",
        "oracle": "pydna",
        "label": label,
        "sources": frag_specs,
        "frag_lens": [len(f) for f in dseq_frags],
        "n_circ_raw": len(circ),
        "expect": expectation(circ),
    }
    c["expect"]["linear"] = None if lin is None else expectation(lin)["linear"]
    if extra:
        c.update(extra)
    return add(c)


def build_ligation_cases():
    for key in TEMPLATES:
        s = tseq(key)
        for names in SETS:
            d = Dseqrecord(s, circular=True)
            try:
                frs = list(d.cut(*enz(names)))
            except Exception:
                continue
            if len(frs) < 1 or len(frs) > MAX_FRAGS:
                continue
            if len(frs) == 1 and frs[0].seq.ovhg == 0 and len(frs[0]) == len(s):
                continue  # no cut at all
            ligation_case(
                "%s circular / %s" % (key, "+".join(names)),
                [{"file": TEMPLATES[key], "topology": "circular",
                  "enzymes": names, "take": "all"}],
                frs)

    # linear slices of the same plasmids
    slices = [("pUC19", 1, 1200), ("pUC19", 400, 2100), ("pBR322", 1, 1500),
              ("pBR322", 900, 3200), ("phiX174", 1, 2000), ("phiX174", 1500, 4800),
              ("pGRAMMAR", 1, 240)]
    for key, a, b in slices:
        sub = tseq(key)[a - 1:b]
        for names in SETS:
            d = Dseqrecord(sub, circular=False)
            try:
                frs = list(d.cut(*enz(names)))
            except Exception:
                continue
            if len(frs) < 2 or len(frs) > MAX_FRAGS:
                continue
            ligation_case(
                "%s[%d..%d] linear / %s" % (key, a, b, "+".join(names)),
                [{"file": TEMPLATES[key], "topology": "linear", "slice": [a, b],
                  "enzymes": names, "take": "all"}],
                frs)

    # vector + insert: one fragment from each of two digests
    combos = [
        ("pUC19", ["EcoRI", "HindIII"], "pBR322", ["EcoRI", "HindIII"]),
        ("pUC19", ["BamHI", "SalI"], "pBR322", ["BamHI", "SalI"]),
        ("pUC19", ["EcoRI", "PstI"], "pBR322", ["EcoRI", "PstI"]),
        ("pUC19", ["KpnI", "SacI"], "pGRAMMAR", ["KpnI", "SacI"]),
        ("pBR322", ["EcoRV", "ScaI"], "phiX174", ["EcoRV", "ScaI"]),
        ("pUC19", ["EcoRI"], "pBR322", ["EcoRI"]),
        ("pBR322", ["HindIII"], "phiX174", ["HindIII"]),
        ("pUC19", ["SmaI"], "pBR322", ["EcoRV"]),
        ("pUC19", ["PstI"], "pBR322", ["PstI"]),
        ("pUC19", ["XbaI", "SphI"], "pBR322", ["XbaI", "SphI"]),
    ]
    # ends that cannot anneal: a 3' overhang against a 5' one, and two
    # different 5' overhangs
    mismatched = [
        ("pUC19", ["EcoRI"], "pBR322", ["PstI"]),
        ("pUC19", ["HindIII"], "pBR322", ["PstI"]),
        ("pUC19", ["EcoRI"], "pBR322", ["SalI"]),
        ("pUC19", ["KpnI"], "pBR322", ["EcoRI"]),
        ("pUC19", ["SacI"], "pBR322", ["HindIII"]),
        ("pUC19", ["EcoRI"], "pBR322", ["AatII"]),
        ("pUC19", ["XbaI"], "pBR322", ["SphI"]),
        ("pUC19", ["SphI"], "pBR322", ["BamHI"]),
    ]
    for ka, ea, kb, eb in mismatched:
        da = Dseqrecord(tseq(ka), circular=True)
        db = Dseqrecord(tseq(kb), circular=True)
        fa, fb = list(da.cut(*enz(ea))), list(db.cut(*enz(eb)))
        if len(fa) != 1 or len(fb) != 1:
            continue
        ligation_case("%s[%s] + %s[%s], ends do not anneal"
                      % (ka, "+".join(ea), kb, "+".join(eb)),
                      [{"file": TEMPLATES[ka], "topology": "circular", "enzymes": ea, "take": [1]},
                       {"file": TEMPLATES[kb], "topology": "circular", "enzymes": eb, "take": [1]}],
                      [fa[0], fb[0]], extra={"designed": "incompatible_ends"})

    for ka, ea, kb, eb in combos:
        da = Dseqrecord(tseq(ka), circular=True)
        db = Dseqrecord(tseq(kb), circular=True)
        try:
            fa = list(da.cut(*enz(ea)))
            fb = list(db.cut(*enz(eb)))
        except Exception:
            continue
        if not fa or not fb or len(fa) > 4 or len(fb) > 4:
            continue
        for i in range(len(fa)):
            for j in range(len(fb)):
                frs = [fa[i], fb[j]]
                if sum(len(f) for f in frs) > 12000:
                    continue
                ligation_case(
                    "%s[%s]#%d + %s[%s]#%d" % (ka, "+".join(ea), i + 1, kb, "+".join(eb), j + 1),
                    [{"file": TEMPLATES[ka], "topology": "circular", "enzymes": ea,
                      "take": [i + 1]},
                     {"file": TEMPLATES[kb], "topology": "circular", "enzymes": eb,
                      "take": [j + 1]}],
                    frs)


# ---------------------------------------------------------------------------
# 2. Gibson / HiFi
# ---------------------------------------------------------------------------

gib_id = 0


def gibson_case(label, parts, limit, designed=None):
    global gib_id
    gib_id += 1
    ds = [Dseqrecord(p["seq"], name=p["name"]) for p in parts]
    circ = safe(gibson_assembly, ds, limit=limit, circular_only=True)
    if circ is None:
        skipped.append((label, "pydna path limit"))
        gib_id -= 1
        return None
    allp = safe(gibson_assembly, ds, limit=limit, circular_only=False)
    lin = None if allp is None else [r for r in allp if not r.circular]
    c = {
        "id": "gib%04d" % gib_id,
        "method": "gibson",
        "oracle": "pydna",
        "label": label,
        "min_homology": limit,
        "parts": parts,
        "n_circ_raw": len(circ),
        "expect": expectation(circ),
    }
    c["expect"]["linear"] = None if lin is None else expectation(lin)["linear"]
    if designed:
        c["designed"] = designed
    elif not c["expect"]["circular"]:
        raise AssertionError("%s: a well-formed Gibson set gave no product" % label)
    return add(c)


def gibson_parts(n, hlen, body, gc=0.5, circular=True, arms=None):
    """n parts in a cycle (or a chain) sharing terminal homology arms."""
    if arms is None:
        arms = [clean(hlen, [], gc) for _ in range(n if circular else n - 1)]
    bodies = [clean(body, []) for _ in range(n)]
    parts = []
    for i in range(n):
        left = arms[i - 1] if (circular or i > 0) else ""
        right = arms[i] if (circular or i < n - 1) else ""
        parts.append({"name": "p%d" % (i + 1), "seq": left + bodies[i] + right})
    return parts


def build_gibson_cases():
    for n in (2, 3, 4, 5):
        for hlen in (15, 18, 20, 25, 30, 40, 60):
            for gc in (0.25, 0.5, 0.75):
                parts = gibson_parts(n, hlen, RNG.choice([120, 180, 260]), gc)
                limit = min(hlen, 20)
                gibson_case("circular %d parts, %d bp arms, GC %.0f%%" % (n, hlen, gc * 100),
                            parts, limit)
    # linear outcomes: a chain, no closing arm
    for n in (2, 3, 4):
        for hlen in (20, 25, 35):
            parts = gibson_parts(n, hlen, 150, circular=False)
            gibson_case("linear chain, %d parts, %d bp arms" % (n, hlen), parts, 20,
                        designed="linear_chain")
    # homology below the floor: must not assemble
    for hlen in (8, 10, 12, 14, 16, 19):
        for n in (2, 3):
            parts = gibson_parts(n, hlen, 200)
            gibson_case("arms of %d bp against a floor of 20" % hlen, parts, 20,
                        designed="homology_too_short")
    # the same parts with the floor lowered: must assemble
    for hlen in (10, 12, 14, 16):
        parts = gibson_parts(3, hlen, 200)
        gibson_case("arms of %d bp against a floor of %d" % (hlen, hlen), parts, hlen,
                    designed="lowered_floor")
    # a repeated arm: the order is ambiguous
    for n in (3, 4):
        for hlen in (25, 30):
            arms = [clean(hlen, []) for _ in range(n)]
            arms[1] = arms[0]  # two joins share the same overlap
            parts = gibson_parts(n, hlen, 180, arms=arms)
            gibson_case("repeated %d bp arm among %d parts" % (hlen, n), parts, 20,
                        designed="repeated_arm")
    # parts with nothing in common: nothing can join
    for n in (2, 3, 4):
        for _ in range(3):
            parts = [{"name": "p%d" % (i + 1), "seq": clean(200, [])} for i in range(n)]
            gibson_case("%d unrelated parts" % n, parts, 20, designed="no_homology")
    # one part in the wrong orientation: assembly must flip it
    for n in (2, 3, 4):
        parts = gibson_parts(n, 25, 160)
        k = n - 1
        parts[k] = {"name": parts[k]["name"], "seq": rc(parts[k]["seq"])}
        gibson_case("%d parts, part %d supplied reverse" % (n, k + 1), parts, 20,
                    designed="flipped_part")


# ---------------------------------------------------------------------------
# 3. Golden Gate
# ---------------------------------------------------------------------------

# site and the number of spacer bases between it and the 4-base overhang:
# BsaI GGTCTC(1/5) and BsmBI CGTCTC(1/5) leave one, BbsI GAAGAC(2/6) two.
IIS = {"BsaI": ("GGTCTC", 1), "BsmBI": ("CGTCTC", 1), "BbsI": ("GAAGAC", 2),
       "Esp3I": ("CGTCTC", 1), "BpiI": ("GAAGAC", 2)}
SITES = sorted(set(v[0] for v in IIS.values()))

# 4-base overhangs for a one-pot reaction: none palindromic, and no two the
# reverse complement of each other -- two parts whose overhangs are reverse
# complements can be joined in either order once one of them is flipped, which
# makes the design ambiguous rather than defined.
OVERHANGS = []
for _a in "ACGT":
    for _b in "ACGT":
        for _c in "ACGT":
            for _d in "ACGT":
                _o = _a + _b + _c + _d
                if _o != rc(_o) and rc(_o) not in OVERHANGS:
                    OVERHANGS.append(_o)
RNG.shuffle(OVERHANGS)
OVERHANGS = OVERHANGS[:24]
assert all(o != rc(o) for o in OVERHANGS)
assert len({min(o, rc(o)) for o in OVERHANGS}) == len(OVERHANGS)

gg_id = 0


def sites_in(s):
    """(forward, reverse) counts of every Type IIS site we use, anywhere."""
    f = sum(s.count(x) for x in SITES)
    r = sum(s.count(rc(x)) for x in SITES)
    return f, r


def gg_part(enzyme, ov_left, ov_right, body_len, pad=8, extra_site=False):
    """A part flanked by inward-facing Type IIS sites. Returns (part, body).

    Redrawn until the finished molecule carries exactly the sites it is meant
    to carry: a site can otherwise be created across the junction between an
    overhang and the body, which would silently change the part.
    """
    site, spacer = IIS[enzyme]
    gap = "A" * spacer
    want_f, want_r = (2, 1) if extra_site else (1, 1)
    for _ in range(4000):
        body = clean(body_len, SITES)
        if extra_site:
            mid = len(body) // 2
            body = body[:mid] + site + gap + "CCGT" + body[mid:]
        s = (clean(pad, SITES) + site + gap + ov_left + body
             + ov_right + gap + rc(site) + clean(pad, SITES))
        if sites_in(s) != (want_f, want_r):
            continue
        # the middle piece must really be ov_left + body + ov_right with a
        # 4-base 5' overhang at each end
        pieces = safe(lambda: list(Dseqrecord(s).cut(enz([enzyme])[0])))
        if pieces is None or len(pieces) != (3 if not extra_site else 4):
            continue
        mid = pieces[1]
        if mid.seq.ovhg != -4 or str(mid.seq)[:4].upper() != ov_left:
            continue
        if not extra_site and str(mid.seq)[-4:].upper() != ov_right:
            continue
        return s, body
    raise RuntimeError("could not draw a clean part")


def goldengate_case(label, parts, enzyme, designed=None, oracle="pydna",
                    expect=None, reason=None):
    global gg_id
    gg_id += 1
    c = {
        "id": "gg%04d" % gg_id,
        "method": "goldengate",
        "oracle": oracle,
        "label": label,
        "enzyme": enzyme,
        "parts": parts,
    }
    if oracle == "pydna":
        ds = [Dseqrecord(p["seq"], name=p["name"]) for p in parts]
        circ = safe(golden_gate_assembly, ds, enz([enzyme]), allow_blunt=False,
                    circular_only=True)
        if circ is None:
            skipped.append((label, "pydna path limit"))
            gg_id -= 1
            return None
        c["n_circ_raw"] = len(circ)
        c["expect"] = expectation(circ)
        c["expect"]["linear"] = None    # see the report: linear intermediates
    else:                               # are not Golden Gate products
        c["expect"] = expect
        c["reason"] = reason
    if designed:
        c["designed"] = designed
    elif oracle == "pydna" and len(c["expect"]["circular"]) != 1:
        raise AssertionError("%s: a well-formed Golden Gate set gave %d products"
                             % (label, len(c["expect"]["circular"])))
    return add(c)


def build_goldengate_cases():
    for enzyme in ("BsaI", "BsmBI", "BbsI"):
        for n in (2, 3, 4, 5, 6):
            for body in (60, 120, 200):
                ovs = OVERHANGS[:n]
                parts = [{"name": "p%d" % (i + 1),
                          "seq": gg_part(enzyme, ovs[i], ovs[(i + 1) % n], body)[0]}
                         for i in range(n)]
                goldengate_case("%s, %d parts, %d bp bodies" % (enzyme, n, body),
                                parts, enzyme)
        # one part supplied on the other strand
        for n in (3, 4):
            ovs = OVERHANGS[:n]
            parts = [{"name": "p%d" % (i + 1),
                      "seq": gg_part(enzyme, ovs[i], ovs[(i + 1) % n], 100)[0]}
                     for i in range(n)]
            parts[1] = {"name": parts[1]["name"], "seq": rc(parts[1]["seq"])}
            goldengate_case("%s, %d parts, part 2 reverse" % (enzyme, n), parts, enzyme,
                            designed="flipped_part")
        # duplicate overhang: two joins compete
        for n in (3, 4, 5):
            bad = list(OVERHANGS[:n])
            bad[1] = bad[0]
            parts = [{"name": "p%d" % (i + 1),
                      "seq": gg_part(enzyme, bad[i - 1], bad[i], 100)[0]}
                     for i in range(n)]
            goldengate_case("%s, %d parts, duplicate overhang %s" % (enzyme, n, bad[0]),
                            parts, enzyme, designed="duplicate_overhang")
        # palindromic overhang
        for pal in ("AATT", "GGCC", "ATAT", "CATG"):
            ovs = [pal, OVERHANGS[1], OVERHANGS[2]]
            parts = [{"name": "p%d" % (i + 1),
                      "seq": gg_part(enzyme, ovs[i - 1], ovs[i], 90)[0]} for i in range(3)]
            goldengate_case("%s, palindromic overhang %s" % (enzyme, pal), parts, enzyme,
                            designed="palindromic_overhang")
        # a part carrying an extra internal site
        for n in (3, 4):
            ovs = OVERHANGS[:n]
            parts = []
            for i in range(n):
                s = gg_part(enzyme, ovs[i], ovs[(i + 1) % n], 120, extra_site=(i == 1))[0]
                parts.append({"name": "p%d" % (i + 1), "seq": s})
            goldengate_case("%s, %d parts, part 2 has an internal site" % (enzyme, n),
                            parts, enzyme, designed="internal_site")

    # ---- hand-built: cases pydna has no way to express ----
    enzyme = "BsaI"
    site, gap = IIS[enzyme][0], "A" * IIS[enzyme][1]
    ovs = OVERHANGS[:3]
    # (a) a part with no BsaI site at all
    parts = [{"name": "p1", "seq": gg_part(enzyme, ovs[0], ovs[1], 100)[0]},
             {"name": "p2", "seq": gg_part(enzyme, ovs[1], ovs[2], 100)[0]},
             {"name": "p3", "seq": clean(300, SITES)}]
    goldengate_case("BsaI, part 3 has no site", parts, "BsaI",
                    designed="no_site", oracle="hand",
                    expect={"circular": [], "linear": None},
                    reason=("pydna's golden_gate_assembly has no way to say 'this part "
                            "cannot take part': it finds no join and silently returns "
                            "nothing, the same answer it gives for a set that is merely "
                            "incomplete. From the model: a part with no site is never cut, "
                            "so it contributes no ligatable end and no circle can contain "
                            "all three parts. Expected: no circular product, and the "
                            "engine must name the part."))
    # (b) a forward site so close to the 3' end that its overhang runs off it
    body = clean(100, SITES)
    p3 = clean(8, SITES) + site + gap + ovs[2] + body + site + gap + "GC"
    parts = [{"name": "p1", "seq": gg_part(enzyme, ovs[0], ovs[1], 100)[0]},
             {"name": "p2", "seq": gg_part(enzyme, ovs[1], ovs[2], 100)[0]},
             {"name": "p3", "seq": p3}]
    goldengate_case("BsaI, p3 ends two bases into its right overhang",
                    parts, "BsaI", designed="truncated_overhang", oracle="hand",
                    expect={"circular": [], "linear": None},
                    reason=("p3's second BsaI site points off the 3' end of the molecule: "
                            "the top-strand cut lands two bases before the end and the "
                            "bottom-strand cut falls past it, so only two of the four "
                            "overhang bases exist. pydna/Biopython and the plugin need not "
                            "even agree on whether that counts as a cut, and pydna has no "
                            "way to report the situation. Either way the model gives the "
                            "same answer: a 2-base end (or an uncut part) cannot anneal to "
                            "the 4-base ends of p1 and p2, so no circle can contain all "
                            "three parts. Expected: no circular product, and the engine "
                            "must report ends of different lengths or an unusable part."))
    # (c) one part whose two overhangs are the same: it closes on itself
    ov = OVERHANGS[0]
    one, body = gg_part(enzyme, ov, ov, 150)
    goldengate_case("BsaI, one part, both overhangs %s" % ov,
                    [{"name": "p1", "seq": one}], "BsaI",
                    designed="self_circularisation", oracle="hand",
                    expect={"circular": [{"len": len(ov) + len(body),
                                          "canon": hashlib.sha256(
                                              canonical(ov + body, True).encode()).hexdigest()}],
                            "linear": None},
                    reason=("pydna routes a one-fragment reaction through "
                            "SingleFragmentAssembly, which also reports insertion and "
                            "inversion products, so its product list is not comparable. "
                            "The answer follows from the construction instead: cutting "
                            "the part releases a middle fragment whose top strand is "
                            "ov + body + ov with a 4-base 5' overhang at each end, and "
                            "closing it writes the shared overhang once, giving the "
                            "circle ov + body."))


# ---------------------------------------------------------------------------

def main():
    build_ligation_cases()
    build_gibson_cases()
    build_goldengate_cases()
    counts = {}
    for c in cases:
        counts[c["method"]] = counts.get(c["method"], 0) + 1
    hand = sum(1 for c in cases if c["oracle"] == "hand")
    for label, why in skipped:
        print("skip %s: %s" % (label, why))
    doc = {
        "generator": "tests/oracle/dump_cloning.py",
        "seed": SEED,
        "pydna": __import__("pydna").__version__,
        "biopython": __import__("Bio").__version__,
        "counts": counts,
        "hand_built": hand,
        "skipped": [{"label": a, "reason": b} for a, b in skipped],
        "cases": cases,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, separators=(",", ":"), sort_keys=False))
    print("%d cases (%s), %d hand-built -> %s (%.1f kB)"
          % (len(cases), counts, hand, OUT, OUT.stat().st_size / 1024))


if __name__ == "__main__":
    main()
