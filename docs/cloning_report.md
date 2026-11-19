# The cloning engine against pydna

`lua/plasmid/bio/assembly.lua` was checked against an independent
implementation: **pydna 5.5.16** on **Biopython 1.88**. pydna models
double-stranded DNA as a `Dseq` (watson strand, crick strand, overhang) and
has first-class routines for all three methods this milestone implements —
`assembly2.ligation_assembly`, `assembly2.gibson_assembly` and
`assembly2.golden_gate_assembly` — so it is a genuine second opinion and not a
re-statement of the same model.

Python is a development tool only. `tests/oracle/dump_cloning.py` writes
`tests/oracle/expected/cloning.json` (468 cases, 329 kB); the suite reads only
that file. `tests/spec/cloning_oracle_spec.lua` runs every case through the
engine offline and compares. Nothing in the loop touches the network.

**Result: 468 of 468 cases agree** — on the product sequences, on the number
of assemblies, on the topology, on the fragment lengths and order coming out
of each digest, and on the cases where nothing should assemble. Two engine
defects were found on the way and fixed (below); neither changed a product
sequence.

## How a product is compared

A circular molecule has no start and no preferred strand, so comparing
sequences literally would report differences that are not differences. Each
product is reduced to a canonical form:

* circular — the lexicographically least rotation (Booth's algorithm, O(n)) of
  the sequence and of its reverse complement, whichever is smaller;
* linear — `min(seq, revcomp(seq))`.

Only the length and the sha256 of that string are committed, which is what
keeps the corpus small. The Lua and Python implementations of Booth's
algorithm are independent; the Lua one is checked against a naive search over
every rotation on 200 random strings in
`canonical_form_ignores_rotation_and_strand`, and the Python one is checked
transitively by the 420 cases whose hashes match.

Per case the spec compares:

1. the fragment lengths of each digest, element by element (ligation cases) —
   this also pins the **order** pydna and the plugin number fragments in;
2. the number of circular assemblies found;
3. the set of distinct circular products;
4. the set of distinct linear products, for the 40 cases where neither side
   produces a circle and pydna's own enumeration did not hit its caps;
5. for Golden Gate, whether the engine reported a problem, against whether the
   case was designed to fail.

## The corpus

468 cases, written by a seeded generator (`SEED = 20260101`) so that it
reproduces byte for byte.

| method | cases | of which designed to fail |
|---|---|---|
| restriction–ligation | 262 | 10 |
| Gibson / HiFi | 125 | 30 |
| Golden Gate | 81 | 8 |

("designed to fail" counts the cases that expect no circular product at all;
a further 39 are designed to be awkward without failing — a part supplied
reverse-complemented, a repeated homology arm, a floor lowered to match short
arms, a duplicate or palindromic overhang, a single part closing on itself.)

420 cases expect at least one circular product (50 expect two, one expects
four, 18 expect eight distinct products); 48 expect none. Three cases are
hand-built rather than taken from pydna.

**Restriction–ligation** — 237 single-digest cases and 25 vector+insert cases
that take one named fragment from each of two digests. Templates: pUC19
(`L09137`, 2686 bp), pBR322 (`J01749`, 4361 bp), phiX174 (`NC_001422`,
5386 bp) and the synthetic `pGRAMMAR` (240 bp), each as the whole circle and
as linear slices of it (126 cases on whole circles, 136 on linear slices).
30 enzymes across 56 sets (24 single, 22 pairs, 10 triples):
single and multiple cutters, 5′ overhangs (EcoRI, HindIII, BamHI, XhoI, SalI,
XbaI, NcoI, NdeI, SpeI, AvrII, BglII, NheI, MfeI, AflII), 3′ overhangs (PstI,
KpnI, SacI, SphI, ApaI, BsiWI, NsiI, AatII, BsrGI), blunt cutters (EcoRV,
SmaI, ScaI, PvuII, HincII, DraI, StuI, NruI, SspI) and mixed sets of two and
three. Digests of 1 to 5 fragments (54 / 151 / 46 / 10 / 1). Eight cases pair
fragments whose ends cannot anneal (a 3′ overhang against a 5′ one, or two
different 5′ overhangs) and must give nothing.

**Gibson / HiFi** — 2 to 5 parts (34 / 40 / 30 / 21), homology arms of 15, 18,
20, 25, 30, 40 and 60 bp at 25 %, 50 % and 75 % GC, bodies of 120–260 bp.
Circular designs (84), linear chains with no closing arm (9), arms below the
floor that must not assemble (12), the same parts with the floor lowered so
that they must (4), a repeated arm (4), a part supplied reverse-complemented
(3), and parts with nothing in common at all (9).

**Golden Gate** — BsaI (29), BsmBI (26) and BbsI (26); 1 to 6 parts
(1 / 9 / 32 / 18 / 12 / 9) with bodies of 60–200 bp. Well-formed sets (45)
plus designed failures: a duplicate overhang (9), a palindromic overhang (12),
a part carrying an extra internal site (6), a part with no site (1), a site
whose overhang runs off the end of the molecule (1), a part supplied reversed
(6) and a single part that closes on itself (1).

### What is *not* in the corpus

* Feature carry-over and junction annotation. The oracle compares molecules;
  pydna's feature handling is its own convention and is not a useful reference
  for ours. Those are covered by `assembly_spec.lua` instead.
* Partial sticky-end overlaps (pydna's `allow_partial_overlap`). The engine
  requires ends of equal length, which is what a ligase does without
  fill-in; the option is not modelled and no case uses it.
* One case was dropped: `phiX174[1..2000]` cut with HincII gives seven blunt
  fragments and pydna refuses to enumerate (`Too many possible paths (more
  than 10000)`). It is recorded in `skipped` in the JSON.

## Disagreements

### 1. Overhang accounting ignored that a part can be supplied either way round — **plugin bug, fixed**

`A.overhang_problems` bucketed each end by the text it reads along the top
strand of the fragment as handed in. A part handed in on the other strand has
both of its ends written as their reverse complements and swapped left for
right, so none of them matched their partners and every end of the reversed
part — and of its two neighbours — was reported as having no partner.

Minimal reproduction: three BsaI parts with overhangs A→B, B→C, C→A, with
part 2 reverse-complemented (`gg0016`, `gg0017`, `gg0042`, `gg0043`,
`gg0068`, `gg0069`).

    overhang GCTT (p1 right) has no partner
    overhang CGAA (p2 left) has no partner
    overhang AAGC (p2 right) has no partner
    overhang TTCG (p3 left) has no partner

The search itself was right throughout: `A.assemble` flips the part and finds
exactly the product pydna finds. Only the report was wrong, and it was wrong
in the worst way — four false alarms on a design that works.

The fix folds an overhang and its reverse complement into one overhang, and
keeps the one thing flipping cannot change: which half of a join an end can
be. Written in the canonical direction, an end is either the left half or the
right half of a join, and two ends that would both be the left half do not
anneal however they are turned. That last check did not exist before and is
new behaviour: two parts that both end in the same overhang used to pass
silently as "used twice".

This also removed a false negative in the corpus itself. The first draft drew
Golden Gate overhangs from a list containing both `TTCG` and `CGAA`, which are
reverse complements, so the designs were ambiguous and pydna returned two
products where one was intended. The engine's new message named it. The
generator now draws overhangs that are neither palindromic nor reverse
complements of one another, which is the rule a real Golden Gate design
follows.

### 2. `max_solutions` capped the answer but not the work — **plugin bug, fixed**

`A.assemble` walks every ordering and orientation of the fragments.
`max_solutions` stopped it *recording* solutions past the cap but not
*searching*, so the cost stayed factorial whatever the cap. Whenever many ends
are interchangeable — a single-enzyme digest, or blunt fragments, where every
end joins every other — this is the difference between milliseconds and never
finishing. The first benchmark run of this track was still inside that series
after 2.7 hours and was killed.

The measured exhaustive cost is in the two-row table below (blunt fragments
and a single-enzyme digest, which differ only in the cost of one
compatibility test); it is also the pre-fix cost, because before the fix every
search was exhaustive. The fix is a flag
set when the cap is reached and checked in the search; the solutions returned
are the same ones the cap would have kept, and all 468 oracle cases and the
295 pre-existing specs are unchanged by it.

### 3. Linear intermediates of a Golden Gate reaction — **convention, recorded**

For Golden Gate the oracle records circular products only. `restriction_
ligation_assembly` filters products that contain a partially digested
sequence, and with `only_adjacent_edges` it does not report chains that merely
run out of compatible ends; the engine reports every linear chain that uses
all the parts. On the six `internal_site` cases the engine finds a linear
product using all three parts (the two ends that cannot meet are left free)
where pydna reports nothing. Neither is wrong: Golden Gate is a
circularisation, a linear chain is not a product of it, and the engine already
says why the circle did not close. The comparison therefore asserts circular
products for this method, and the engine's problem list instead of pydna's
silence.

### 4. pydna's own enumeration caps — **recorded**

pydna raises rather than enumerate past its limits: `Too many possible paths
(more than 10000)` while filtering cycles, and `Too many assemblies (96
pre-validation)` while listing linear assemblies. Where the circular
enumeration raised, the case was dropped (one case). Where only the linear
enumeration raised, the circular expectation is kept and the linear one is
recorded as "no answer" and not compared.

### Non-disagreements worth recording

* **Fragment order and length.** For all 262 ligation cases the engine's
  fragments match pydna's one for one, including the length. Both number a
  circular digest from the first cut and both count the overhang into the
  fragment: pUC19 cut once with EcoRI gives one fragment of 2690 nt, not 2686.
* **Assembly counts.** On all 465 pydna-backed cases the engine finds exactly
  as many circular assemblies as pydna does. Both deduplicate by fixing the
  first fragment of a cycle in one orientation, which removes rotations and
  the mirror image.
* **Terminal homology length.** The engine takes the longest exact terminal
  match; pydna's `terminal_overlap` can return several lengths for one join.
  They cannot differ for a non-repetitive arm: an exact terminal match of
  length *n* extends to *n+1* only when the overlap is a homopolymer, so the
  two always chose the same join in this corpus.
* **Blunt ligation.** pydna's `ligation_assembly` defaults to
  `allow_blunt=False`; the oracle passes `allow_blunt=True` throughout,
  because the engine treats two blunt ends as compatible (which is what T4
  ligase does). This is a default, not a disagreement.

## The hand-built cases

Three Golden Gate cases have no pydna answer. Each carries its reasoning in
the JSON (`reason`), and the expectation is derived from the model, never from
what the engine does.

1. **`gg0079` — a part with no BsaI site.** pydna has no way to say "this part
   cannot take part": it finds no join and returns nothing, the same answer it
   gives for a set that is merely incomplete. From the model, a part with no
   site is never cut, contributes no ligatable end, and no circle can contain
   all three parts. Expected: no circular product, and the engine must name
   the part. It does — `p3 has no BsaI site`.

2. **`gg0080` — a site whose overhang runs off the end of the molecule.** p3's
   second BsaI site sits two bases from the 3′ end, so the top-strand cut
   lands inside the molecule and the bottom-strand cut falls past it: only two
   of the four overhang bases exist. Biopython and the plugin need not agree
   on whether that counts as a cut at all, and pydna cannot report the
   situation either way. Both readings give the same answer: a 2-base end — or
   an uncut part — cannot anneal to the 4-base ends of p1 and p2, so there is
   no circular product.

3. **`gg0081` — one part whose two overhangs are the same.** pydna routes a
   single-fragment reaction through `SingleFragmentAssembly`, which also
   reports insertion and inversion products, so its list is not comparable.
   The answer follows from the construction: cutting the part releases a
   middle fragment whose top strand is `ov + body + ov` with a 4-base 5′
   overhang at each end, and closing it writes the shared overhang once,
   giving the circle `ov + body`. The expected hash is computed from that
   string directly.

## Performance

Neovim 0.12.2, headless, one core; every figure below is from one run of
`tests/tools/cloning_bench.lua`. "joins" counts calls to `A.compatible` /
`A.homology`, which is what the search spends its time on.

### The reactions a user actually runs

| | time | joins | result |
|---|---|---|---|
| Golden Gate, 5 parts of ~700 bp — cut, search and build | **6.9 ms** | 21 | one product, 3520 bp |
| of which the search alone | 0.9 ms | 21 | |
| religate phiX174 / HaeIII (11 fragments) | 12.7 ms | 150 | capped at 50 |
| religate phiX174 / HinfI (21 fragments) | 17.4 ms | 532 | capped at 50 |
| religate phiX174 / AluI (24 fragments) | 16.0 ms | 176 | capped at 50 |
| religate pBR322 / HinfI (10 fragments) | 9.7 ms | 514 | capped at 50 |
| religate NC_001133 230 kb / SbfI (5 fragments) | 78.1 ms | 380 | capped at 50 |
| religate NC_001133 230 kb / EcoRI (80 fragments) | 114.0 ms | 532 | capped at 50 |

The search is not the cost on a large record; cutting it is.

| digest | time | fragments |
|---|---|---|
| phiX174 / HaeIII | 27.3 ms | 11 |
| phiX174 / AluI | 54.6 ms | 24 |
| pBR322 / HinfI | 31.4 ms | 10 |
| NC_001133 230 kb / NotI | 132.9 ms | 2 |
| NC_001133 230 kb / SbfI | 314.9 ms | 5 |
| NC_001133 230 kb / EcoRI | **4892 ms** | 80 |

Cutting a 230 kb record into 80 fragments takes nearly five seconds — about
61 ms a fragment — because every fragment is extracted as a record of its own,
and `A.copy` and `A.flip` round-trip a record through the GenBank writer and
parser. The cost is per fragment and per feature, not per base: the same
record cut into 5 pieces takes 315 ms.

### How the search grows

A defined design — every overhang used exactly twice — prunes to a single
path, and the work is quadratic in the number of parts, not factorial:

| parts | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|
| joins | 3 | 7 | 13 | 21 | 31 | 57 | 91 | 133 |
| time (ms) | 0.2 | 0.2 | 0.3 | 0.5 | 0.4 | 0.5 | 0.5 | 0.6 |
| solutions | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |

The worst case is fragments whose ends are all interchangeable. Two series
were measured, because the cost per join is not the same for both: blunt
fragments, where a compatibility test is two integer comparisons, and a
single-enzyme (EcoRI) digest religated, where it also compares the overhang
text and so runs about twice as long for the same number of tests. The search
itself is identical — same joins, same orderings. This is the behaviour
before the `max_solutions` fix, since every search was exhaustive then.

| fragments | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| joins | 4 | 20 | 126 | 1 016 | 10 170 | 122 052 | 1 708 742 |
| orderings found | 2 | 8 | 48 | 384 | 3 840 | 46 080 | 645 120 |
| all-blunt, time (ms) | 0.7 | 0.5 | 0.4 | 1.6 | 11.5 | 208.4 | 3384.2 |
| one-enzyme digest, time (ms) | 0.1 | 0.2 | 0.5 | 2.8 | 23.3 | 303.7 | 4037.7 |

Each added fragment multiplies the join count by 12–14 and the time by 13–18.
At the sticky-end rate of ~2.4 µs a join, nine fragments would be about a
minute and ten about a quarter of an hour; twelve is out of reach. With the
cap in force the same series is flat, because the search stops at 50
orderings:

| fragments | 2 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|
| all-blunt, time (ms) | 0.1 | 0.9 | 0.9 | 0.9 | 0.8 | 1.0 | 2.2 | 8.2 |
| one-enzyme digest, time (ms) | 0.3 | 0.7 | 0.7 | 0.9 | 1.2 | 1.4 | 1.7 | 1.9 |

(the joins are 4 to 168 across that whole row: the cap is reached after a few
dozen tests however many fragments are on the bench.)

## Limits a user should know about

* **The engine enumerates orderings, and `max_solutions` (default 50) is what
  keeps that bounded.** When a set has many interchangeable ends the 50
  orderings reported are the first 50 found, not a sample and not all of them.
  A design whose answer should be unique — distinct overhangs, or distinct
  homology arms — is unaffected: it prunes to one path.
* **Up to about 8 fragments with interchangeable ends can be enumerated
  exhaustively** (3.4 s blunt, 4.0 s for a single-enzyme digest); past that,
  raising `max_solutions` is not viable. A
  one-enzyme digest of more than a handful of fragments has a combinatorial
  number of religation products and listing them is not a useful question to
  ask.
* **Defined designs scale well.** A 12-part Golden Gate with distinct
  overhangs searches in 0.6 ms with 133 compatibility tests, so the practical
  cap on part count is the biology, not the engine.
* **Cutting dominates on large records.** Budget roughly 60 ms per fragment
  when a record of a few hundred kb is cut into many pieces. A plasmid-sized
  digest is immediate.
* **Golden Gate parts are selected one per record.** `A.type_iis_parts` keeps
  the piece that no longer carries a site for the enzyme, and when a record
  yields several it takes the longest and says so. A part with an internal
  site is therefore reported, not silently resolved — pydna would instead
  consider every piece as a candidate.
* **Gibson joins are exact terminal matches**, longest first, with a floor of
  `min_homology` (default 20) and a ceiling of `max_homology` (default 200).
  Mismatched or partially annealing arms are not modelled.
* **Ends must be the same length to anneal.** Partial sticky-end overlaps and
  fill-in/chew-back are not modelled.

## Reproducing

    conda activate primerdev          # pydna 5.5.16, biopython 1.88
    python tests/oracle/dump_cloning.py          # rewrites the committed JSON
    nvim -l tests/run.lua cloning_oracle         # the comparison, offline
    nvim -l tests/tools/cloning_compare.lua out.json   # every disagreement, in bulk
    nvim -l tests/tools/cloning_bench.lua bench.json   # the numbers above
