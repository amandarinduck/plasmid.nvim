# plasmid.nvim: GenBank/GenPept parser checked against real NCBI records

> Status note (v0.6 release): all six bugs pinned by this report have since
> been fixed in the plugin; their specs in `tests/spec/ncbi_records_spec.lua`
> are now named `fixed_*` and pass. The performance figures are unchanged.

Track: Record checks (sequence-fetch milestone). Plugin runtime code was not modified. Oracle: Biopython 1.88 (`SeqIO.parse(..., "genbank")`), in a dedicated dev environment. Records fetched on 17 Nov 2026 with Biopython Entrez (`tool=plasmid.nvim`, no API key, no email, default User-Agent, at most about 2.5 requests per second).

## Deviations

1. **one-of() and remote references (`ACC.1:1..20`).** I scanned about 420 candidate records in four discovery batches (pre-1999 exon, partial-cds, immunoglobulin and mitochondrial records, plus vectors) and found neither form in any record NCBI serves today. NCBI has rewritten them: the classic `one-of(1888,1901)` in U18266 is now `misc_feature 1901` in segmented set AH005336. Both forms are therefore tested with **synthetic** records (`synthetic/one_of.gbk`, `synthetic/remote_refs.gbk`). Their expected values still come from Biopython.
2. **Old-style headers (GI, NID, pre-2003 LOCUS columns, BASE COUNT).** efetch returns every record in the current layout. I fetched pre-2000 records (J01636, X06403, U13850, X52328, S75221, AAB02572) in that current layout. The GI, NID and old-LOCUS variants are **synthetic** edits of U49845 (`synthetic/old_header_gi.gbk`, `synthetic/old_header_nid.gbk`).
3. **"Multi-segment" record.** No record I fetched has a SEGMENT header line. I used two substitutes: segmented sets AH005336 and AH003512 (several `source` features, `gap` features, secondary accessions) and the 8-segment influenza PR8 genome as one multi-record file.
4. **CON records.** Each CON/scaffold record I used has a single-component `CONTIG join(...)`. A CONTIG line containing `gap()` was not tested.
5. **Big records.** E. coli NC_000913 (4.6 Mb, 11.4 MB file) and M. genitalium NC_000908 (580 kb) were timed but not committed, as the task asked.
6. **Timing conditions.** Open and render were timed in a headless Neovim (`--headless --clean`, 120×40) with no UI attached, so the decoration provider (highlighting) did not run. Memory is the Lua heap delta (`collectgarbage("count")`) and the process RSS after the open, not peak RSS. All timings come from one run on one machine (12 cores) and should be treated as indicative.
7. **Scope additions.**
   - The translation check found a bug in `plasmid.bio.translate`, which is not the parser. It is pinned as a spec case because real NCBI records trigger it.
   - The spec includes a timing-ratio case for quadratic parsing. It has a wide margin (it fails at ×8, the measured value is ×14.5), but it remains timing-based.
8. **Fixture extension.** Fixtures use `.gbk`/`.gp` rather than `.gb`. The existing `**/*.gb` checks in `io_spec`/`location_spec` would otherwise pick them up and, because the oracle here is a single JSON file, report a raw feature count mismatch.
9. **Discovery traffic.** About 430 extra records were downloaded in 6 batched requests (about 1.5 MB) to search for rare location syntax. None of them were committed.

## Summary

- **Corpus:** 41 fetched files containing 48 records, 13,454 features in total. The committed spec corpus is 50 files and 57 records: 39 real files, 4 minimal reproductions and 7 synthetic grammar cases. Fixtures total 2.5 MB and the oracle JSON 1.9 MB, so 4.4 MB is committed, under the 8 MB budget.
- **What already agrees with Biopython:**
  - every residue of every record that has residues;
  - name, id, description, topology and molecule type (including `cRNA`, `ss-RNA`, `rRNA` and protein);
  - every feature count and feature type;
  - 3,063 of 3,071 locations (the other 8 are the `bond()` and `one-of()` locations, counted in both the real and the minimal fixtures);
  - all 12,467 qualifier values, including TTN's 35,991-character `/translation`, quoted values that start with `/`, inner `""` quotes, and multi-line unquoted values.
- **Other existing invariants that hold:**
  - 44 of 50 files round-trip byte-for-byte;
  - 395 CDS `/translation` values are reproduced from the sequence. This covers origin-spanning joins (SV40, mtDNA D-loop), trans-spliced rps12, SARS-CoV-2's ribosomal-slippage join, `strand=2` extracts and `<`/`>` partials;
  - the canonical writer (every feature regenerated) re-reads identically in Biopython for every record that has a sequence.
- **Six distinct bugs**, each pinned by one failing case in `tests/spec/ncbi_records_spec.lua`:
  - four come from records NCBI serves (GenPept `bond()`; CON records without a sequence: writer and circular view; WGS master `rc` units);
  - one is in translation of IUPAC-degenerate codons;
  - one is an O(n²) parse of multi-line qualifiers.
- **Performance:** files up to about 300 kB open and render in under 0.3 s. E. coli K-12 (11.4 MB) parses in 0.85 s but takes **6.1 s** to open and render and uses **636 MB of Lua heap (945 MB RSS)**. Almost all of that is the linear layout of 9,285 features; the same sequence as FASTA takes 1.5 s and 7 MB.

Spec status: `nvim -l tests/run.lua ncbi_records` gives 7 passed and 6 failed (the 6 bug cases). The full suite gives 241 passed and 6 failed, so no existing spec changed status. `PLASMID_TEST_GROUP=0` gives the same result.

## Records and why they were chosen

| # | File | Accession(s), efetch parameters | Why chosen |
|---|---|---|---|
| 1 | U55762.gbk | U55762.1 (nuccore, rettype=gbwithparts) | pEGFP-N1: commercial circular vector (sparse 2003 annotation: 4 features) |
| 2 | X52328.gbk | X52328.1 (nuccore, rettype=gbwithparts) | pBluescript II: classic circular cloning vector, EMBL-origin 1990 record |
| 3 | U13850.gbk | U13850.1 (nuccore, rettype=gbwithparts) | pGEX-2T: circular GST expression vector, 1994 SYN record with 13 features |
| 4 | X06403.gbk | X06403.1 (nuccore, rettype=gbwithparts) | pACYC184: pre-2000 circular plasmid |
| 5 | NC_001669.gbk | NC_001669.1 (nuccore, rettype=gbwithparts) | SV40 RefSeq: circular, origin-spanning joins (join(...,1..N)) |
| 6 | NC_012920.gbk | NC_012920.1 (nuccore, rettype=gbwithparts) | human mtDNA rCRS: circular, D-loop spans origin, tRNA/rRNA, transl_except, overlapping genes |
| 7 | NC_000932.gbk | NC_000932.1 (nuccore, rettype=gbwithparts) | Arabidopsis chloroplast: circular 154 kb, trans-spliced rps12, many complement(join()) |
| 8 | NC_002483.gbk | NC_002483.1 (nuccore, rettype=gbwithparts) | E. coli F plasmid: bacterial plasmid with many features |
| 9 | NC_001416.gbk | NC_001416.1 (nuccore, rettype=gbwithparts) | phage lambda: 48.5 kb phage genome |
| 10 | NC_045512.gbk | NC_045512.2 (nuccore, rettype=gbwithparts) | SARS-CoV-2: ss-RNA, ribosomal-slippage join with overlapping parts, mat_peptide, long notes |
| 11 | NC_001802.gbk | NC_001802.1 (nuccore, rettype=gbwithparts) | HIV-1: RNA virus, LTRs, join() CDSs |
| 12 | influenza_PR8.gbk | NC_002016.1–NC_002023.1 (nuccore, rettype=gbwithparts) | influenza A/PR/8/34: 8-segment genome as one multi-record file (/segment qualifiers) |
| 13 | NM_000546.gbk | NM_000546.6 (nuccore, rettype=gbwithparts) | TP53 mRNA (RefSeq NM_): many misc_feature, variation, regulatory |
| 14 | XM_062871673.gbk | XM_062871673.1 (nuccore, rettype=gbwithparts) | predicted model mRNA (XM_, Gnomon): mouse Srsx |
| 15 | NM_001267550.gbk | NM_001267550.2 (nuccore, rettype=gbwithparts) | TTN mRNA: ~105 kb, /translation of ~35,000 aa (very long qualifier) |
| 16 | NM_000581.gbk | NM_000581.4 (nuccore, rettype=gbwithparts) | GPX1 mRNA: selenocysteine /transl_except |
| 17 | NR_003286.gbk | NR_003286.4 (nuccore, rettype=gbwithparts) | 18S rRNA (NR_): non-coding RefSeq |
| 18 | NC_000017_TP53_rev.gbk | NC_000017.11 (nuccore, rettype=gbwithparts, seq_start=7668402, seq_stop=7687550, strand=2) | TP53 locus extract, strand=2 (gene on minus strand shown forward), spliced mRNA/CDS joins |
| 19 | NC_000011_HBB_fwd.gbk | NC_000011.10 (nuccore, rettype=gbwithparts, seq_start=5225464, seq_stop=5229395, strand=1) | HBB locus extract, plus strand (gene on minus: complement(join())) |
| 20 | NC_000011_HBB_rev.gbk | NC_000011.10 (nuccore, rettype=gbwithparts, seq_start=5225464, seq_stop=5229395, strand=2) | same HBB extract, strand=2 |
| 21 | NC_000012_GAPDH_fwd.gbk | NC_000012.12 (nuccore, rettype=gbwithparts, seq_start=6534000, seq_stop=6539000) | GAPDH locus extract, features truncated by the window get < / > |
| 22 | NZ_JCEXQL010000002_con.gbk | NZ_JCEXQL010000002.1 (nuccore, rettype=gb) | RefSeq CON record as rettype=gb: CONTIG join(), no ORIGIN |
| 23 | NZ_JCEXQL010000002_parts.gbk | NZ_JCEXQL010000002.1 (nuccore, rettype=gbwithparts) | same CON record expanded with gbwithparts |
| 24 | NW_028328278_con.gbk | NW_028328278.1 (nuccore, rettype=gb) | RefSeq scaffold (NW_) as rettype=gb |
| 25 | CGRMEV000000000.gbk | CGRMEV000000000.2 (nuccore, rettype=gb) | WGS master record: no sequence, WGS / WGS_SCAFLD lines |
| 26 | CGRMEV020000001.gbk | CGRMEV020000001 (nuccore, rettype=gbwithparts) | WGS contig of the master above |
| 27 | AH005336.gbk | AH005336.3 (nuccore, rettype=gbwithparts) | segmented set (AH_): several source features, gap features, many secondary accessions |
| 28 | AH003512.gbk | AH003512.2 (nuccore, rettype=gbwithparts) | segmented set with between-base misc_recomb sites (769^770) |
| 29 | S75221.gbk | S75221.1 (nuccore, rettype=gbwithparts) | between-base location 212^213 (misc_recomb), MMTV int-6-LTR fusion, 1990s record |
| 30 | AF084526.gbk | AF084526.1 (nuccore, rettype=gbwithparts) | between-base misc_feature 2638^2639 (RNA-editing insertion sites), Physarum mitochondrial cox3 |
| 31 | X61312.gbk | X61312.1 (nuccore, rettype=gbwithparts) | order() gene location with > partial |
| 32 | J01636.gbk | J01636.1 (nuccore, rettype=gbwithparts) | E. coli lac operon: 1980s record |
| 33 | U49845.gbk | U49845.1 (nuccore, rettype=gbwithparts) | NCBI's own sample record: <1..>206 partial CDS, complement() |
| 34 | NP_000537.gp | NP_000537.3 (protein, rettype=gp) | RefSeq protein p53: Region/Site features, order() sites |
| 35 | NP_005219.gp | NP_005219.2 (protein, rettype=gp) | RefSeq protein EGFR: many Region/Site, long /note |
| 36 | AAB02572.gp | AAB02572.1 (protein, rettype=gp) | GenPept (INSDC, non-RefSeq) protein: EGFP from cloning vector pEGFP-1, 1996 record |
| 37 | P01308.gp | P01308.1 (protein, rettype=gp) | UniProtKB-derived GenPept (insulin): Bond features bond(a,b), long COMMENT |
| 38 | P04637.gp | P04637.2 (protein, rettype=gp) | UniProtKB-derived GenPept (P53_HUMAN): hundreds of Site/Region features |
| 39 | 4HHB_A.gp | 4HHB_A (protein, rettype=gp) | PDB-derived GenPept: SecStr and Het features (bond()) |
| 40 | NC_000908.gbk | NC_000908.2 (nuccore, rettype=gbwithparts) | Mycoplasma genitalium chromosome 580 kb (scaling point) |
| 41 | NC_000913.gbk | NC_000913.3 (nuccore, rettype=gbwithparts) | E. coli K-12 MG1655 chromosome 4.6 Mb (big record) |

Spec-only fixtures (not fetched as-is):

- `minimal/bond_genpept.gp`: P01308 header and sequence, carrying one `Bond bond(31,96)` from P01308 and one `Het bond(87)` from 4HHB_A.
- `minimal/con_contig.gbk` and `minimal/con_features.gbk`: NZ_JCEXQL010000002 cut down to its LOCUS, DEFINITION, ACCESSION and VERSION lines, a source feature (plus one gene in the second file) and the CONTIG line.
- `minimal/wgs_master.gbk`: CGRMEV000000000 cut down to LOCUS `807 rc`, one source feature and the WGS line.
- `synthetic/*`: `one_of`, `remote_refs`, `between_circular` (`300^1` on a circular record), `mixed_strand_order`, `quote_edge` (`""` at the end of a line, empty `/note=""`, a repeated `/pseudo`), `old_header_gi` and `old_header_nid`.

## Results per record

Columns:
- "Features equal" compares feature type plus the location parts as Biopython interprets them, using the existing `{start, stop, strand, ref}` model from `dump_biopython.py`.
- "Qualifier diffs" counts features whose sorted (key, value) lists differ.
- "open + first render" is the time for `:edit` (read, parse, initial linear layout, set lines). "Other view" is toggling to the circular view.

| # | File | Recs | Length | Features | Parse ok | Header fields | Sequence equal | Features equal (type + location) | Qualifier diffs | Byte round-trip | Plugin parse ms | Biopython parse ms | open + first render ms | other view ms | Lua heap after open MB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | U55762.gbk | 1 | 4,733 bp | 4 | yes | equal | yes | yes | 0 | yes | 0.5 | 0.3 | 10 | 5 | 0.7 |
| 2 | X52328.gbk | 1 | 2,961 bp | 2 | yes | equal | yes | yes | 0 | yes | 0.3 | 0.2 | 8 | 3 | 0.4 |
| 3 | U13850.gbk | 1 | 4,948 bp | 13 | yes | equal | yes | yes | 0 | yes | 0.7 | 0.4 | 15 | 5 | 0.9 |
| 4 | X06403.gbk | 1 | 4,245 bp | 8 | yes | equal | yes | yes | 0 | yes | 0.5 | 0.3 | 11 | 4 | 0.4 |
| 5 | NC_001669.gbk | 1 | 5,243 bp | 62 | yes | equal | yes | yes | 0 | yes | 2.9 | 2.6 | 27 | 7 | 1.1 |
| 6 | NC_012920.gbk | 1 | 16,569 bp | 77 | yes | equal | yes | yes | 0 | yes | 4.6 | 2.3 | 30 | 9 | 1.3 |
| 7 | NC_000932.gbk | 1 | 154,478 bp | 259 | yes | equal | yes | yes | 0 | yes | 18.3 | 9.0 | 264 | 11 | 4.5 |
| 8 | NC_002483.gbk | 1 | 99,159 bp | 213 | yes | equal | yes | yes | 0 | yes | 16.1 | 6.9 | 126 | 9 | 13.2 |
| 9 | NC_001416.gbk | 1 | 48,502 bp | 285 | yes | equal | yes | yes | 0 | yes | 11.2 | 6.9 | 98 | 8 | 8.9 |
| 10 | NC_045512.gbk | 1 | 29,903 bp | 57 | yes | equal | yes | yes | 0 | yes | 9.0 | 2.0 | 75 | 8 | 2.1 |
| 11 | NC_001802.gbk | 1 | 9,181 bp | 43 | yes | equal | yes | yes | 0 | yes | 2.8 | 1.3 | 33 | 5 | 2.6 |
| 12 | influenza_PR8.gbk | 8 | 13,588 bp | 36 | yes | equal | yes | yes | 0 | yes | 3.2 | 3.3 | 15 | 5 | 0.8 |
| 13 | NM_000546.gbk | 1 | 2,512 bp | 70 | yes | equal | yes | yes | 0 | yes | 2.6 | 1.8 | 18 | 5 | 0.8 |
| 14 | XM_062871673.gbk | 1 | 1,491 bp | 3 | yes | equal | yes | yes | 0 | yes | 0.3 | 0.3 | 10 | 4 | 0.5 |
| 15 | NM_001267550.gbk | 1 | 109,224 bp | 369 | yes | equal | yes | yes | 0 | yes | 140.2 | 7.9 | 303 | 8 | 19.6 |
| 16 | NM_000581.gbk | 1 | 899 bp | 16 | yes | equal | yes | yes | 0 | yes | 0.8 | 0.6 | 12 | 5 | 0.5 |
| 17 | NR_003286.gbk | 1 | 1,869 bp | 4 | yes | equal | yes | yes | 0 | yes | 0.5 | 0.3 | 9 | 3 | 0.4 |
| 18 | NC_000017_TP53_rev.gbk | 1 | 19,149 bp | 432 | yes | equal | yes | yes | 0 | yes | 18.7 | 13.7 | 88 | 11 | 4.7 |
| 19 | NC_000011_HBB_fwd.gbk | 1 | 3,932 bp | 73 | yes | equal | yes | yes | 0 | yes | 4.3 | 2.1 | 30 | 6 | 0.9 |
| 20 | NC_000011_HBB_rev.gbk | 1 | 3,932 bp | 73 | yes | equal | yes | yes | 0 | yes | 3.2 | 1.9 | 26 | 6 | 0.8 |
| 21 | NC_000012_GAPDH_fwd.gbk | 1 | 5,001 bp | 143 | yes | equal | yes | yes | 0 | yes | 6.3 | 4.0 | 36 | 8 | 2.1 |
| 22 | NZ_JCEXQL010000002_con.gbk | 1 | 39,953 bp (declared) | 91 | yes | equal except length (plugin 0) | n/a (no residues) | yes | 0 | NO | 4.6 | 2.5 | 13 | 1560 | 0.8 |
| 23 | NZ_JCEXQL010000002_parts.gbk | 1 | 39,953 bp | 91 | yes | equal | yes | yes | 0 | yes | 6.4 | 2.8 | 53 | 9 | 5.5 |
| 24 | NW_028328278_con.gbk | 1 | 281,759,365 bp (declared) | 1 | yes | equal except length (plugin 0) | n/a (no residues) | yes | 0 | NO | 0.2 | 0.2 | 5 | 2 | 0.3 |
| 25 | CGRMEV000000000.gbk | 1 | 807 bp (declared) | 1 | yes | equal except length (plugin 0) | n/a (no residues) | yes | 0 | NO | 0.3 | 0.4 | 5 | 1 | 0.2 |
| 26 | CGRMEV020000001.gbk | 1 | 1,012 bp | 1 | yes | equal | yes | yes | 0 | yes | 0.2 | 0.3 | 6 | 1 | 0.2 |
| 27 | AH005336.gbk | 1 | 13,855 bp | 28 | yes | equal | yes | yes | 0 | yes | 3.1 | 1.5 | 29 | 4 | 0.7 |
| 28 | AH003512.gbk | 1 | 2,398 bp | 22 | yes | equal | yes | yes | 0 | yes | 1.6 | 1.1 | 14 | 5 | 0.6 |
| 29 | S75221.gbk | 1 | 687 bp | 6 | yes | equal | yes | yes | 0 | yes | 0.3 | 0.3 | 8 | 4 | 0.4 |
| 30 | AF084526.gbk | 1 | 3,505 bp | 36 | yes | equal | yes | yes | 0 | yes | 2.9 | 1.1 | 13 | 6 | 0.7 |
| 31 | X61312.gbk | 1 | 679 bp | 7 | yes | equal | yes | yes | 0 | yes | 0.3 | 0.3 | 8 | 5 | 0.4 |
| 32 | J01636.gbk | 1 | 7,477 bp | 20 | yes | equal | yes | yes | 0 | yes | 1.2 | 1.0 | 16 | 6 | 1.4 |
| 33 | U49845.gbk | 1 | 5,028 bp | 9 | yes | equal | yes | yes | 0 | yes | 0.7 | 0.4 | 14 | 4 | 1.0 |
| 34 | NP_000537.gp | 1 | 393 aa | 63 | yes | equal | yes | yes | 0 | yes | 2.1 | 1.6 | 12 | 0 | 0.5 |
| 35 | NP_005219.gp | 1 | 1,210 aa | 60 | yes | equal | yes | yes | 0 | yes | 2.4 | 1.7 | 15 | 0 | 0.6 |
| 36 | AAB02572.gp | 1 | 239 aa | 4 | yes | equal | yes | yes | 0 | yes | 0.2 | 0.2 | 6 | 0 | 0.2 |
| 37 | P01308.gp | 1 | 110 aa | 53 | yes | equal | yes | 50/53 (3 bond) | 0 | yes | 2.2 | 1.8 | 11 | 0 | 0.5 |
| 38 | P04637.gp | 1 | 393 aa | 290 | yes | equal | yes | yes | 0 | yes | 8.0 | 7.3 | 21 | 0 | 1.2 |
| 39 | 4HHB_A.gp | 1 | 141 aa | 11 | yes | equal | yes | 10/11 (1 bond) | 0 | yes | 0.4 | 0.5 | 8 | 0 | 0.3 |
| 40 | NC_000908.gbk (measured only) | 1 | 580,076 bp | 1133 | yes | equal | yes | yes | 0 | yes | 101.7 | 42.5 | 754 | 23 | 85.1 |
| 41 | NC_000913.gbk (measured only) | 1 | 4,641,652 bp | 9285 | yes | equal | yes | yes | 0 | yes | 844.7 | 365.1 | 6124 | 103 | 635.7 |

Synthetic and minimal fixtures, compared in the same spec:
- All 7 synthetic files agree with Biopython on header fields, sequence, locations and qualifiers. The exceptions are the two one-of locations (by design, see below) and the repeated value-less `/pseudo`: Biopython drops the duplicate, while the plugin keeps the file as written. I consider the plugin right there, and the spec normalises Biopython's behaviour.
- The 4 minimal files reproduce bugs 1–4.

## Bugs

Each bug has one failing case in `tests/spec/ncbi_records_spec.lua`. The `corpus_*` cases skip only the items a bug case pins, and they count the skips.

### 1. GenPept `bond()` locations are rejected (case `bug_genpept_bond_location`)

- **Where it shows up:** P01308 (UniProtKB-derived; 3 `Bond` features) and 4HHB_A (PDB-derived; `Het bond(87)`). NCBI uses this form for every UniProt or PDB protein that has disulfides or heterogen contacts.
- **Plugin:** `loc_error = 'unsupported location operator "bond"'`, `parts = {}`. The feature is kept and written back verbatim, but it is not drawn and cannot be navigated or searched.
- **Biopython:** `bond{[30:31], [95:96]}` → parts `{31,31}`, `{96,96}`; `bond(87)` → `{87,87}`.
- **Minimal reproduction** (`minimal/bond_genpept.gp`):
  ```
       Bond            bond(31,96)
                       /bond_type="disulfide"
       Het             bond(87)
                       /heterogen="(HEM,4385)"
  ```
- **Who is right:** Biopython, in practice. `bond` is not part of the INSDC nucleotide feature table, whose operators are complement, join and order. It is NCBI's protein flat-file form for a Seq-loc of type bond, and the plugin accepts `.gp` GenPept files.
- **Suggested fix:**
  - Add `bond` to the location parser's function list as a node kind that lists one or two sites.
  - Flatten each argument to a site part.
  - Format it back as `bond(a,b)` so writing stays lossless.
  - In the protein view, draw the two residues as linked ticks (or as two sites).

### 2. Writing a record that has no sequence adds an `ORIGIN` block (case `bug_con_record_write_adds_origin`)

- **Where it shows up:** efetch `rettype=gb` on CON or scaffold records (NZ_JCEXQL010000002, NW_028328278) and on WGS masters (CGRMEV000000000). Each ends `CONTIG …` or `WGS …` and then `//`, with no ORIGIN line.
- **Plugin:** writing the unmodified record inserts `ORIGIN` before `//`. This breaks the plugin's own byte-for-byte policy. The saved file then has a 39,953 bp LOCUS line above an empty sequence block. Biopython still treats that as an undefined sequence, but stricter readers can take it as 0 bp. With `locus_dirty` (any header edit), the regenerated LOCUS line says `0 bp`, and Biopython then reads a defined empty sequence: the record has silently become a different one.
- **Minimal reproduction** (`minimal/con_contig.gbk`, 10 lines): LOCUS `39953 bp … CON`, one source feature `1..39953`, then `CONTIG      join(JCEXQL010000002.1:1..39953)` and `//`.
- **Who is right:** the input file. Biopython keeps the record as undefined with length 39,953, plus `annotations["contig"]`.
- **Suggested fix:**
  - In `format_record`, emit `ORIGIN` only when `raw.origin_line` was present or the record has residues.
  - Mark sequence-less records at parse time (for example `rec.raw.no_origin`, or `rec.contig` taken from `pre_origin`), and use that mark both to refuse sequence edits and to regenerate LOCUS from `declared_length` rather than `rec:len()`.
  - For the fetch milestone, request `rettype=gbwithparts` for nucleotide records. With it, the same CON record arrives with its 39,953 residues and parses identically to Biopython (row 23).

### 3. The circular view draws features past the end of the sequence (case `bug_circular_view_draws_features_past_sequence_end`)

- **Where it shows up:** NZ_JCEXQL010000002 fetched as `rettype=gb` (0 residues, 45 features up to position 39,953, circular).
- **Plugin:** the linear view shows `(empty sequence)` and draws no features. The circular view uses `L = max(len, 1)` and draws every feature as a multi-turn spiral with labels. Toggling to that view takes **1.5–1.6 s**; every other record toggles in under 25 ms except E. coli (about 100 ms).
- **Minimal reproduction** (`minimal/con_features.gbk`): the CON record of bug 2 plus `gene 1000..2000 /gene="repA"`. `circular.build` then draws `repA` on a 0 bp record.
- **Who is right:** there is no Biopython counterpart for rendering. The expected behaviour is defined by the plugin's own linear view: a feature beyond the sequence cannot be placed on it.
- **Suggested fix:**
  - In `circular.build`, use the same rule as the linear view: when `rec:len() == 0`, draw no features and print a "no sequence (CON/WGS record)" note.
  - More generally, clip parts to `[1, L]` (or skip parts outside it) before computing angles, so that layout cost depends on L and not on the feature coordinates.
  - A load-time warning for features that extend past the sequence end would also catch truncated files.

### 4. A WGS master record opens without any warning (case `bug_wgs_master_opens_without_warning`)

- **Where it shows up:** CGRMEV000000000, whose LOCUS line is `CGRMEV020000000          807 rc    DNA     linear   ENV`. Here `rc` is a record count, not base pairs: it gives the number of contig records in the project, listed on the `WGS` line.
- **Plugin:** `parse_locus` recognises only `bp`/`aa`, so `declared_length` is nil. The record opens as an empty 0 bp sequence, without the "LOCUS declares … but sequence has 0" warning that CON records get.
- **Biopython:** reports an undefined sequence of length 807. Taking a record count as a length is itself questionable.
- **Minimal reproduction:** `minimal/wgs_master.gbk` (11 lines).
- **Who is right:** NCBI's convention (`rc` = record count). Neither tool should present 807 as a length. The plugin's defect is that it gives no explanation.
- **Suggested fix:**
  - Accept `rc` in `parse_locus` (for example `units = "rc"`, `record_count = 807`).
  - On open, warn along the lines of "WGS master: 807 contig records, no sequence; see WGS line CGRMEV020000001-CGRMEV020000807".
  - In the fetch feature, follow the WGS range or refuse a master record with that message.
  - Keep the LOCUS line verbatim on write. Bug 2's fix covers the missing ORIGIN.

### 5. Codons with IUPAC codes become X even when only one amino acid is possible (case `bug_translate_degenerate_codons`)

This is not a parser bug. The translation cross-check on real records found it.

- **Where it shows up:**
  - AH003512, CDS `join(<456..833,2196..>2398)`: codon 117 is `GCN`, annotated as A; the plugin gives X.
  - AH005336, TMPO CDS `join(2201..2479,…,10475..>10521)`: the 3'-partial CDS ends in an incomplete codon `TC`. NCBI translates it as S, because `TCN` is serine whatever the third base; the plugin drops it.
- **Plugin:** returns X for all 15 such codons in the oracle table (`GCN`, `GGN`, `CTN`, `TCN`, `CCN`, `ACN`, `GTN`, `CGN`, `AAR`, `GAY`, `TAR`, `TRA`, `YTR`, `MGR`, `ATH`). `NNN` and `GAN` correctly give X.
- **Biopython and NCBI:** A, G, L, S, P, T, V, R, K, D, \*, \*, L, R, I.
- **Who is right:** NCBI and Biopython. The codon table is fully determined in these cases.
- **Suggested fix:**
  - For a codon containing IUPAC codes, expand it to the set of concrete codons it can stand for (at most 64) and return the amino acid if all expansions agree (or `*` if all are stops), X otherwise.
  - For the final partial codon of a 3'-partial CDS, pad with N and apply the same rule.
  - Return X rather than B/Z/J for mixed sets, to stay compatible with NCBI.
- **Effect on other modules:** `plasmid.bio.translate` is also used by `bio/impact.lua`, which recomputes `/translation` after an edit, and by `seqedit`, `orfs`, `search` and `view/features`. An edit-impact report on a CDS containing such a codon would therefore show X where the record has a residue.

### 6. Multi-line qualifier parsing is O(lines²) (case `bug_multiline_qualifier_parse_is_quadratic`)

- **Where it shows up:** every long multi-line qualifier. Real example: TTN NM_001267550, `/translation` of 35,991 aa over 621 lines. Its 322 kB file takes **140 ms** to parse, while the 302 kB chloroplast takes 18 ms.
- **Measurement:** `/translation` of 10k, 20k, 40k and 80k aa parses in 10.6, 31, 120 and 452 ms. Each doubling costs about ×3.8, and 20k→80k costs ×14.5.
- **Cause:** for each continuation line, the reader runs `table.concat(cur_qual, " ")` and rescans the whole value for balanced quotes (`genbank.lua`, the `open = cur_qual and not quotes_balanced(...)` line).
- **Who is right:** not a correctness question. The output agrees with Biopython; only the cost is wrong.
- **Suggested fix:** keep a running quote parity per qualifier. Count `"` only in the new line, and treat the value as open while the count since the opening quote is odd. A doubled `""` adds 2 and leaves parity unchanged. On well-formed values this gives the same result as the current `quotes_balanced`. It differs only on malformed input such as text after the closing quote, which the current check keeps open.
- **Impact:** TTN is about the longest protein that exists, so real records cost at most about 0.15 s. The case is pinned because the growth is algorithmic, and a GenPept file containing many TTN isoforms would multiply the cost.

### Disagreements that are not plugin bugs

- **`one-of()`:** `one-of(1888,1901)..2200` is rejected with a location error. The feature, its qualifiers and its text are kept, and writing is lossless. Biopython accepts it (start 1888). The current INSDC location operators are complement, join and order, and NCBI no longer serves one-of. Rejecting it is acceptable; the case `legacy_one_of_degrades_gracefully` passes and pins this.
- **Repeated value-less qualifiers:** `/pseudo /pseudo` is kept twice by the plugin and once by Biopython. The plugin is lossless, so it is right.
- **Length of CON / WGS records:** Biopython reports the declared length (39,953; 281,759,365; 807) for an undefined sequence, and the plugin reports 0. This is the root of bugs 2–4. I did not count it as a separate bug.
- **Gap-spanning CDSs in segmented sets:** CDSs with a fuzzy end in the middle of a join (a gap of unknown size) cannot be checked against `/translation`. The corpus check skips them by rule; in this corpus, the only such CDS fails anyway for the reason in bug 5.

## Performance

| Record | File | Parse (plugin / Biopython) | Open + linear render | Circular view | Lua heap after open | RSS |
|---|---|---|---|---|---|---|
| NC_000932 chloroplast, 154 kb, 259 features | 302 kB | 18 / 9 ms | 264 ms | 11 ms | 4.5 MB | 52 MB |
| NM_001267550 TTN, 109 kb, 35,991-aa `/translation` | 322 kB | **140** / 8 ms | 303 ms | 8 ms | 19.6 MB | 48 MB |
| NC_000908 M. genitalium, 580 kb, 1,133 features | 1.5 MB | 102 / 43 ms | 754 ms | 23 ms | 85 MB | 152 MB |
| NC_000913 E. coli K-12, 4.64 Mb, 9,285 features | 11.4 MB | 845 / 365 ms | **6,124 ms** | 103 ms | **636 MB** | **945 MB** |
| same E. coli sequence as FASTA (no features) | 4.7 MB | – | 1,481 ms | 29 ms | 6.7 MB | 150 MB |

Observations:

- **Parsing scales linearly in file size**, at roughly 70 ms per MB, apart from bug 6. It is about 2× slower than Biopython, which is fine.
- **The cost of big records is in the linear layout, not the parser.**
  - For E. coli, `linear.build` alone takes 4.2 s and allocates 511 MB, producing 256,745 lines and a `rows` table of 256,744 entries. The parsed record itself is 46 MB.
  - Without features (the FASTA), the same layout takes 1.0 s and 39 MB.
  - So each annotated row costs about 2 kB of Lua tables.
  - The existing `performance_large_records` budget (15 s for 5 Mb) still passes, but 945 MB RSS for one bacterial genome will be a problem for users fetching whole chromosomes.
  - Suggested direction (not patched): build rows lazily per visible window (or per block of N rows) and keep only line strings for the rest; alternatively, limit feature tracks per block for records over about 1 Mb.
- **The circular view is cheap** for every record (at most 103 ms), except bug 3.
- **Small and medium records** (up to about 230 kB, under 300 features) all open in under 100 ms. The slowest committed fixtures are the chloroplast (264 ms, many overlapping features) and TTN (303 ms: 140 ms of parsing plus a long layout).

## Files

The tarball is laid out relative to the plugin root:

- `tests/spec/ncbi_records_spec.lua`: 7 passing corpus cases and 6 failing bug cases.
- `tests/fixtures/ncbi_real/`:
  - 39 real records (`.gbk` nucleotide, `.gp` protein);
  - `manifest.json`, with accession, efetch parameters, reason, size and whether the file is committed, including the two big records that are not;
  - `minimal/` (4 files) and `synthetic/` (7 files).
- `tests/oracle/dump_ncbi_real.py`: the Biopython oracle. It reuses `parts()` and `quals()` from `dump_biopython.py`, adds a degenerate-codon table, and writes `tests/oracle/expected/ncbi_real.json`.
- `tests/tools/fetch_ncbi_real.py`: refetches the corpus (`python tests/tools/fetch_ncbi_real.py OUTDIR`).
- `tests/tools/ncbi_real_dump.lua`: the plugin side of the comparison, as JSON in the oracle's shape, plus parse time, heap and round-trip result.
- `tests/tools/ncbi_open_time.lua`: headless open, render and toggle timing plus memory, for one file.
