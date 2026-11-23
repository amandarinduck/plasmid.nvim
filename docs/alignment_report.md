# The alignment engine against Biopython

`lua/plasmid/bio/align.lua` is checked against **Biopython 1.88**'s
`PairwiseAligner`, which implements the same affine-gap model from an
independent codebase. Python is a development tool only:
`tests/oracle/dump_alignments.py` writes
`tests/oracle/expected/alignments.json` (1,097 cases, 586 kB) and
`tests/spec/align_oracle_spec.lua` replays it offline. Nothing in the loop
touches the network.

**Result: 1,097 of 1,097 cases agree on score — 0 disagreements.** 870 of them
reproduce Biopython's alignment column for column; 189 produce a different
alignment of the same (optimal) score; 38 are cases where both sides decline
to align at all. Of the 414 cases where Biopython reports that the optimum is
unique, the engine reproduces all 414 exactly.

## Mapping the two parameter sets

| engine | Biopython `PairwiseAligner` |
|---|---|
| `match` / `mismatch` | `match_score` / `mismatch_score` |
| `gap_open + gap_extend` | `open_gap_score` |
| `gap_extend` | `extend_gap_score` |
| `mode = "global"` | `mode = "global"` |
| `mode = "local"` | `mode = "local"` |
| `mode = "semi"` (free reference ends) | `mode = "global"`, `end_deletion_score = 0` |
| `mode = "semi", free_query_ends` | `mode = "global"`, `end_insertion_score = 0` |
| `protein = true` | `substitution_matrix = BLOSUM62` |

Biopython charges a gap of length *n* as `open_gap_score + (n-1) *
extend_gap_score`, the engine as `gap_open + n * gap_extend`; the first row of
the table is what makes the two identical. Biopython calls a gap in the query
a *deletion* and a gap in the target an *insertion*, so freeing end deletions
is what lets the reference hang off the ends — the engine's `semi`.

BLOSUM62 itself is not transcribed by hand: `tests/oracle/dump_matrices.py`
writes `lua/plasmid/data/blosum62.lua` from Biopython's copy of the NCBI
matrix.

## The corpus

1,097 cases from a seeded generator, so it reproduces byte for byte.

| family | cases | | mode | cases | | scoring scheme | cases |
|---|---|---|---|---|---|---|---|
| fixture_mutated | 175 | | global | 216 | | blastn (2/-3/-5/-2) | 202 |
| hangoff_both | 112 | | local | 268 | | unit (1/-1/-1/-1) | 202 |
| hangoff_left | 112 | | semi | 265 | | steep | 188 |
| hangoff_right | 112 | | semi_q | 216 | | linear | 174 |
| random_related | 100 | | overlap | 132 | | longgap | 118 |
| random_unrelated | 100 | | | | | blastp | 71 |
| protein_pair | 90 | | | | | cheap | 71 |
| protein_fragment / protein_mutated | 60 / 60 | | | | | near_linear | 71 |
| short_query | 56 | | | | | | |
| disjoint | 53 | | | | | | |
| revcomp | 42 | | | | | | |
| empty_query | 25 | | | | | | |

`fixture_mutated` pairs come from the test plasmids (`tests/fixtures/ncbi/`)
with substitutions and indels; `protein_*` from
`tests/fixtures/protein/proteins.faa` scored through BLOSUM62. Per family and
per mode the agreement is exact: no slice has a score mismatch or a
sub-optimal alignment. Replaying the whole corpus takes 3.7 s.

## What the comparison found

Three defects, all fixed, none of which survived into a release:

1. **A gap could not follow a gap of the other kind.** The recurrence let the
   insertion and deletion states be entered from a match, but not from each
   other. With cheap gap opens or an expensive mismatch, "insert then delete"
   can beat aligning two runs against each other, and Biopython took that
   path where the engine could not.
2. **Semi-global alignments missed end states.** With free reference ends the
   alignment may end in a gap state, and the whole-query-inserted case (the
   answer when two sequences share nothing) was unreachable; both are now end
   candidates.
3. **Local alignments were not trimmed.** Smith–Waterman's traceback can leave
   a leading or trailing run that contributes nothing — a mismatch cancelled
   by the match beside it under a scheme where they cost the same. Such an
   alignment is co-optimal, but every other implementation prints the trimmed
   form, and it was the only thing separating the engine from Biopython on six
   unique-optimum cases.

A fourth change came out of the read simulation rather than the corpus: the
seed step now groups well-voted diagonals into clusters and takes the heaviest
cluster, instead of widening the band to cover every good diagonal. A circular
reference is doubled before a read is placed, so a read crossing the origin
matches the start of the second copy as well as the junction; the two peaks
tie, and merging them would have made the band span the whole molecule.

### Conventions, not defects

* Where several alignments share the optimal score the engine may print a
  different one (189 cases). The spec rescores the engine's alignment under
  the same scheme and requires the optimum, so a co-optimal answer passes and
  a worse one does not.
* Biopython spells free end gaps out as columns; the engine omits them and
  reports `ref_start`/`ref_stop` instead. The committed CIGARs are trimmed of
  the end gaps the mode makes free so the two forms are comparable.

## The banded path

Above `max_cells` (4,000,000 by default) the engine aligns inside a band
around the diagonals that shared 11-mers vote for.

* On 39 read-against-plasmid and plasmid-against-plasmid pairs — the sizes
  that actually take this path — the band reproduces the full matrix's score
  exactly, with nothing declined.
* Forced onto the whole corpus (`max_cells = 1`), the band reproduces the full
  matrix on 312 of 335 high-identity cases and declines 225 more rather than
  guess. It is not meant for unrelated sequences: at low identity it matches
  on 76 of 292, and on the 10 cases with no shared k-mer at all it declines
  every time. That is the intended behaviour — the band is for sequences that
  are versions of each other.

## Simulated reads

288 reads cut from pBR322, pUC19 and phiX174 (200–1,200 bp, half
reverse-complemented, 119 of them crossing the origin), with substitutions and
indels at four error rates, placed with `A.place` against the circular record:

| substitution rate | reads | placed | right strand | wrap flag right | differences exactly right | recall | precision |
|---|---|---|---|---|---|---|---|
| 0 % | 72 | 72 | 72 | 72 | 72 | — | — |
| 0.5 % | 72 | 72 | 72 | 72 | 70 | 99.2 % | 99.2 % |
| 2 % | 72 | 72 | 72 | 72 | 52 | 97.0 % | 97.2 % |
| 5 % | 72 | 72 | 72 | 72 | 15 | 91.7 % | 92.5 % |

Every read is placed within 2 bp of where it came from, on the right strand,
with the right wrap flag, at every error rate — including all 119 that cross
the origin. "Differences exactly right" counts reads where the recovered list
matches the introduced mutations exactly; recall and precision count
individual differences. The losses at 5 % are what the model implies rather
than a bug: when mutations fall within a few bases of each other an aligner
cannot tell three substitutions from a short indel, and both descriptions have
the same score.

## Speed and memory

Neovim 0.12.2, headless, one core. "cells" is `#ref * #query`.

| pair | reference | query | cells | path | time | Lua heap |
|---|---|---|---|---|---|---|
| read vs pBR322 | 4,361 | 200 | 8.7e5 | full | 107 ms | 3.1 MB |
| read vs pBR322 | 4,361 | 400 | 1.7e6 | full | 211 ms | 5.8 MB |
| read vs pBR322 | 4,361 | 799 | 3.5e6 | full | 412 ms | 11 MB |
| read vs pBR322 | 4,361 | 1,203 | 5.2e6 | band | 13 ms | 0.9 MB |
| read vs pBR322 | 4,361 | 2,006 | 8.8e6 | band | 17 ms | 1.2 MB |
| pBR322 vs pBR322 | 4,361 | 4,361 | 1.9e7 | band | 25 ms | 2.2 MB |
| pPCP1 vs itself | 9,609 | 9,605 | 9.2e7 | band | 59 ms | 3.3 MB |
| protein vs protein | 396 | 286 | 1.1e5 | full | 12 ms | 0.6 MB |
| 73 kb region vs itself | 73,308 | 73,308 | 5.4e9 | band | 315 ms | 23 MB |
| 1 kb read vs 73 kb region | 73,308 | 1,001 | 7.3e7 | band | 110 ms | 3.8 MB |

The full matrix keeps one traceback byte per cell per state:

| matrix | traceback | working rows |
|---|---|---|
| 1,000 × 1,000 | 2.9 MB | 0.05 MB |
| 2,000 × 2,000 | 11.4 MB | 0.09 MB |
| 4,361 × 4,361 | 54.4 MB | 0.20 MB |
| 9,609 × 9,609 | 264 MB | 0.44 MB |

which is why `max_cells` defaults to 4e6: above that the band is both faster
and far smaller.

## Limits a user should know about

* **The band needs shared 11-mers.** Two sequences that are versions of each
  other band well; unrelated sequences above `max_cells` are declined with a
  message rather than aligned badly. Lower `max_cells` only if you are willing
  to pay the memory in the table above.
* **Reads are placed, not assembled.** Each read is aligned independently;
  there is no consensus, no quality handling and no trace reading (`.ab1`).
* **A difference list is one description among several.** Where mutations
  cluster, an equally scoring alternative exists; at 5 % error about one
  difference in twelve is described differently from the mutation that caused
  it.
* **No multiple alignment.** Everything here is pairwise.
* The whole oracle spec (1,097 replayed cases plus 288 simulated reads) runs
  in about 73 s, which is most of the suite's running time.

## Reproducing

    conda activate primerdev            # biopython 1.88
    python tests/oracle/dump_matrices.py           # BLOSUM62 -> lua/plasmid/data
    python tests/oracle/dump_alignments.py         # rewrites the committed JSON
    nvim -l tests/run.lua align_oracle             # the comparison, offline
    nvim -l tests/tools/align_report.lua           # every number in this report
