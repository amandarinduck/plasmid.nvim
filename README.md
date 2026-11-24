# plasmid.nvim

DNA sequence viewing, editing and annotation inside Neovim. Open a GenBank,
FASTA or SnapGene file and get an annotated sequence view, or a circular plasmid
map. Toggle between them with `<Tab>`. Edit residues with Vim-style keys or as
plain text (`:DnaEdit`) while features follow their bases, select features as
text objects, add your own, find open reading frames, detect known elements
against an imported feature database, search IUPAC motifs, show restriction
sites, simulate digests with a text gel, and work with primers: annotate them,
calculate NEB-style Tm and Ta, design pairs and simulate PCR, and assemble
constructs by Gibson, Golden Gate or restriction–ligation, and align
sequencing reads or whole plasmids against a record. Fetch records from
NCBI by accession, search query or gene symbol. Protein sequences open too.
Pure Lua; `curl` is needed only for NCBI fetch.

```
nvim pBR322.gb
```

## Install

lazy.nvim:

```lua
{ "amandarinduck/plasmid.nvim", lazy = false, opts = {} }
```

Keep `lazy = false`: the plugin reads and writes sequence files itself
(`BufReadCmd`), so it has to be loaded before one is opened, including
`nvim pBR322.gb` from the shell. Filetype-based lazy loading does not work,
because these buffers get no filetype until the plugin has read them. Loading
takes about 15 ms. Add `version = "*"` to follow tagged releases only.

vim-plug: `Plug 'amandarinduck/plasmid.nvim'`. Built-in packages:
`git clone https://github.com/amandarinduck/plasmid.nvim
~/.local/share/nvim/site/pack/plasmid/start/plasmid.nvim`, then
`:helptags ALL`.

Without `setup()`, the plugin auto-configures with defaults (or from
`vim.g.plasmid = { ... }` set before load). It opens `.fa`/`.fasta` and other
FASTA files too; set `patterns` (see `lua/plasmid/config.lua`) to limit which
extensions it takes over.

Requirements: Neovim ≥ 0.10, `ambiwidth=single`, a font with box-drawing and
braille glyphs (nearly all coding fonts), and `termguicolors` recommended.
Run `:checkhealth plasmid`. Feature detection needs an imported database (see
below); importing pLannotate's BLAST database additionally needs NCBI BLAST+
(`blastdbcmd`) on PATH once, at import time. NCBI fetch needs `curl` on PATH
and network access to `eutils.ncbi.nlm.nih.gov`.

## Screenshots
<details> 
<summary> Click to expand </summary>
  
| | |
|---|---|
| ![Detection proposals on pUC19](docs/screenshots/detect_proposals.png) | ![Track limit in a crowded region](docs/screenshots/priority_track_limit.png) |
| `:DnaDetect` on a bare pUC19 FASTA: 11 proposals from the pLannotate set, in italics until accepted | `:DnaTracks 2` in the MCS region: the CDS and primer sites keep their rows; the MCS (`misc_feature`) is summarised as "+1 more" |
| ![Feature selection and ORFs](docs/screenshots/feature_selection_orfs.png) | ![Circular map with overlays](docs/screenshots/puc19_circular_overlays.png) |
| `vaf` on AmpR (selection 885..1745 in the winbar) with the ORF overlay (`go`) showing an 88 aa reverse-frame ORF | The same record on the circular map: accepted features, remaining fragment proposals and ORFs |
| ![Custom feature editor](docs/screenshots/custom_feature_editor.png) | ![Codon highlight](docs/screenshots/codon_highlight.png) |
| A promoter added with `ga` over 3110..3175 of pBR322, open in the `gE` editor | Cursor on M55 of TetA: its codon ATG (248..250) is highlighted on both strands; the winbar reads `codon 55/397 ATG→M 248..250` |
| ![Codon highlight on a reverse-strand CDS](docs/screenshots/codon_highlight_reverse.png) | ![Strand picker](docs/screenshots/yank_strand_picker.png) |
| On reverse-strand β-lactamase, K252 highlights CTT/GAA at 3398..3400, read as AAG on the bottom strand | Visual `gy`: yank the selection as either strand, in either direction, with previews |
| ![Continuous layout](docs/screenshots/continuous_layout.png) | ![Insert inside a CDS](docs/screenshots/edit_insert_codon.png) |
| `:DnaGroup 0`: continuous sequence with a ruler; codons are never split by a gap and residues sit three columns apart | `i` + `GCT` at 500: tet is extended in frame and its new codon 139 (GCT→A) is highlighted across the line break |
| ![Impact report](docs/screenshots/edit_impact_report.png) | ![DnaEdit round trip](docs/screenshots/dnaedit_round_trip.png) |
| Deleting 4140..4158 cuts bla's start codon: `:DnaEditReport` lists the truncated features, the frameshift, the new stop and the removed `/translation` | `:DnaEdit 450..650`: the region as plain text; `:w` applied one 6-bp insertion and extended tet in frame |
| ![Motif search](docs/screenshots/search_hits.png) | ![Restriction sites, linear view](docs/screenshots/enzymes_linear.png) |
| `/GGTCTC 1`: 25 BsaI-like sites with ≤ 1 mismatch on pUC19; bottom-strand hits on the complement row, one across the origin (2685..2) | `gz`: unique 6+ cutters. The bar sits on the cut (top strand above, bottom strand below the complement); isoschizomers share a label |
| ![Restriction sites, circular map](docs/screenshots/enzymes_circular.png) | ![Digest table](docs/screenshots/digest_table.png) |
| The same set on the circular map: cut ticks across the backbone, polylinker labels fanned out | `:DnaDigest EcoRI PstI HindIII AvaI` on pBR322: fragments with their ends and overhangs, features carried, and the lane on a gel |
| ![Gel lanes](docs/screenshots/gel_lanes.png) | ![Primers on the sequence](docs/screenshots/primers_linear.png) |
| `:DnaGel EcoRI \| EcoRI HindIII \| PvuII PstI \| uncut`: one lane per digest next to a 1 kb ladder; the 31 bp fragment runs off | M13 primers on pUC19: forward above the top strand, reverse below the complement, the EcoRI tail in lowercase |
| ![Primer details](docs/screenshots/primer_info.png) | ![Primer design](docs/screenshots/primer_design.png) |
| `K` on a primer: Tm for the current polymerase, GC and 3' clamp, end stability, dimers, and the Ta for the pair it makes | `:DnaPrimerDesign name=tet fwd=EcoRI rev=BamHI` over pBR322's tet CDS: pairs ranked by Tm match, with what each one costs |
| ![PCR product](docs/screenshots/pcr_product.png) | ![Fetched record](docs/screenshots/fetch_record.png) |
| The chosen pair amplified: a new 1212 bp record starting with the EcoRI tail, tet carried across and both primers re-annotated | `:DnaFetch NM_000518.5`: the human β-globin mRNA from NCBI, 628 nt and 47 features, opened as an unsaved buffer |
| ![Links from a feature](docs/screenshots/fetch_links.png) | ![Gene region](docs/screenshots/gene_region.png) |
| `K` on the HBB CDS: the `links` line lists the protein (NP_000509.1) and NCBI Gene 3043; `gF` opens them | `:DnaGene thrA Escherichia coli str. K-12 substr. MG1655 flank=60`: only the gene region (337..2799 + 60 bp each side) is downloaded |
| ![Search picker](docs/screenshots/search_picker.png) | ![Golden Gate report](docs/screenshots/clone_report.png) |
| `:DnaSearch pBR322[title] AND complete sequence[title]`: hits with size, topology and title; pick one to open it | `:DnaGoldenGate BsmBI` on two staged parts: the product above, and the report with parts, ends, junctions and carried features below |
| ![Assembled plasmid](docs/screenshots/clone_product.png) | ![Nothing assembles](docs/screenshots/clone_failure.png) |
| The 3,361 bp product on the circular map: both parts and both junction overhangs annotated beside the carried tet and bla | Parts whose ends cannot anneal: every ordered pair of ends with the reason, and what to change |

</details>

## Status (v0.9)
<details> 
<summary> Click to expand </summary>

| | |
|---|---|
| Read/write GenBank (`.gb .gbk .gbff .ape`) | yes. Unmodified files round-trip byte-for-byte |
| Read/write FASTA (`.fa .fasta .fna .ffn .fas .fsa`) | yes. Multi-record, soft-masking, CRLF preserved |
| Read SnapGene `.dna` / `.rna` / `.prot` | yes. Features, primers, notes, colours; RNA shown with U |
| Write SnapGene | files opened from SnapGene: unmodified saves are byte-identical; added/edited features regenerate the Features packet with untouched features' XML kept verbatim. Creating `.dna` from other formats is not supported yet |
| Protein FASTA (`.faa .pep`, or detected by content) and GenPept (`.gp .gpff`) | yes. Residue numbering in aa, no complement/translation/arrows, linear view only |
| Format conversion | `:w copy.fasta`, `:w copy.gb` from any buffer (including `.dna` and `.prot`) |
| Full location grammar | `join`, `order`, `complement`, nested forms, `<`/`>` fuzzy ends, `a^b`, `a.b`, remote `ACC.1:x..y`, origin-spanning joins |
| Linear view | residues in groups of 10 or continuous with a ruler (`:DnaGroup`), complement strand, packed feature tracks with strand arrows and introns, per-CDS translation (own frame, `/transl_table`, `/codon_start`) |
| Circular view | braille-canvas map, feature rings with arrowheads, tick scale, labels with leader lines |
| Crowded regions | features placed by a configurable priority score; per-line track limit with a "+N more" summary |
| Selection | `af`/`if` text objects (`if` = spliced, strand-aware residues); yank the top or bottom strand in either direction |
| Codon highlight | cursor on an amino acid highlights its codon on both strands (also across line breaks, introns and the origin) |
| Feature authoring | add (`ga`), edit (`gE`, GenBank syntax), delete, undo/redo, personal feature library |
| ORFs | six-frame ORF overlay (`go`), any NCBI genetic code, across the origin; promote an ORF to a CDS |
| Feature detection | `:DnaDetect` against imported databases (e.g. pLannotate's SnapGene-derived set) and your library; review and accept proposals |
| Sequence editing | insert, delete, replace, reverse-complement, rotate the origin; Vim-style keys, commands, `:DnaEdit` plain-text editing; features follow their bases; impact report with frameshift/stop detection; `/translation` kept current; one undo step per edit |
| Search | IUPAC motifs on both strands and across the origin, with mismatches (`/`, `?`, `n`, `N`, `gn`, `:DnaFind`); amino-acid motifs in six-frame translations or protein records |
| Restriction enzymes | 222 common commercial enzymes bundled (Type IIS included); unique / 1–2 / all / named sets in both views, cuts on both strands, overhangs, Dam/Dcm warnings where known; import the full REBASE set with `:DnaEnzymeImport` |
| Digests | `:DnaDigest`: fragment table with ends, overhangs and features; select or extract a fragment as a new record; text gel with ladders and several lanes (`:DnaGel`) |
| Primers | annotate where a primer binds (3'-anchored, mismatches, 5' tails), drawn base-by-base over the template; primer details with `K` |
| Tm / Ta | the NEB Tm Calculator model for 22 polymerase products, verified against NEB's own API; `:DnaTm`, `:DnaPolymerase`, DMSO, custom salt and primer concentration |
| Primer design | `:DnaPrimerDesign` ranks pairs for the selection (Tm window, ΔTm, GC, 3' clamp, unique binding), optional restriction-site tails |
| PCR | `:DnaPCR` builds the product as a new record, with tails, mismatches and the template's features |
| Cloning | `:DnaClone` (method inferred), `:DnaLigate`, `:DnaGoldenGate`, `:DnaGibson`: parts staged from records, selections, digest fragments or PCR products; sticky, blunt and homology joins in both orientations; features carried and re-joined, parts and junctions annotated; a report that says why nothing assembled |
| Assembly-aware primer design | `:DnaPrimerDesign into={part}` adds the homology or Type IIS tails the assembly needs, then amplifies and assembles the chosen pair to check it closes |
| Large records | a bacterial genome opens in about a second: above `linear.window.threshold` the view lays out only the rows a window is showing, and `:DnaOverview` draws the whole record one row per span, with a jump into the sequence |
| Alignment | `:DnaAlign` places a query on the record (either strand, across the origin), `:DnaAlignReads` places a whole run of reads, `:DnaCompare` compares two molecules and ignores a different origin; differences are listed with the feature, codon and effect, and can be written on as `variation` features |
| Mouse | drag to select a region, double click to select a feature, Ctrl-click two primers to select the product they make |
| NCBI fetch | `:DnaFetch` by accession (nucleotide or protein, optional sub-range and strand), `:DnaSearch` with a picker, `:DnaGene` for a gene region by symbol, `gF` follows a feature's protein/gene links; local cache, size guard, NCBI rate limits honoured |

### Formats in detail

- **SnapGene**: the cookie identifies DNA (type 1), protein (2) and RNA (7).
  Features (`Segment` ranges, wrapped origin ranges, `gap` segments), primer
  binding sites (0-based, filtered by the file's hybridisation parameters) and
  notes are read with the same semantics as Biopython's `SnapGeneIO`. Everything
  else (history, enzymes, display settings, traces) is preserved byte for byte.
  One deliberate difference: a reverse-strand feature crossing the origin is
  stored as `complement(join(a..L,1..b))`, so its parts are in biological
  order. SnapGene text fields are HTML: the plugin decodes entities (`&amp;` →
  `&`) for display and encodes them on write, which Biopython does not.
  `:DnaRaw` on a SnapGene file shows a packet table and the XML packets.
  Edited primers cannot be written to SnapGene yet; such a write is refused
  rather than dropping them.
- **`.prot`** means a SnapGene protein file. MacVector also uses `.prot`; those
  files are reported as unsupported rather than misread.
- **Alphabet detection** (FASTA): letters outside the IUPAC nucleotide codes
  (E, F, I, L, P, Q, …) or < 60 % A/C/G/T/U/N mark a protein. Degenerate DNA
  (`NNKNNK…`, primers with R/Y/W) stays DNA, and `-` and `*` are ignored. `.faa`,
  `.pep` and `.gp` force protein.
- **File colours**: SnapGene segment colours and ApE/Benchling
  `/ApEinfo_fwdcolor` / `/ApEinfo_revcolor` are used when legible on light and
  dark backgrounds (`use_file_colors = false` to colour by type only).
</details>

## Use
<details> 
<summary> Click to expand </summary>


| Key (buffer-local) | Action |
|---|---|
| `<Tab>` | toggle linear ↔ circular (cursor position is kept as a bp coordinate) |
| `<CR>` | linear: feature details · circular: jump to that position in the linear view |
| `K` | details of the feature under the cursor (qualifiers, translation) |
| `]f` / `[f` | next / previous feature |
| `gO` | feature picker (`vim.ui.select`, so telescope/fzf-lua pick it up) |
| `+` / `-` | show every feature track / back to the track limit |
| `vaf` `yaf` | select / yank the whole span of the feature under the cursor |
| `vif` `yif` | the feature's own residues, spliced and 5'→3' on its strand |
| visual `y` / `Y` | yank the bases under the selection / their reverse complement |
| visual `gy` | yank in a chosen orientation: top 5'→3', top 3'→5', bottom 5'→3', bottom 3'→5' |
| `gy` / `gY` | yank the current selection / its reverse complement |
| `<Esc>` | clear the selection |
| `go` | toggle the ORF overlay |
| `ga` | add a feature over the selection; on a proposal, accept it; on an ORF, make it a CDS |
| `gE` | edit the feature under the cursor (`:w` applies, `:q` discards) |
| `u` / `<C-r>` | undo / redo (features and residues) |
| `i` / `a` | insert residues before / after the cursor base |
| `[count]x` · `[count]r{letter}` | delete bases · replace bases |
| `p` / `P` | insert register contents after / before the cursor base |
| visual `d` `x` · `c` · `p` · `R` | delete · replace (prompt) · replace with a register · reverse-complement the selection |
| `daf` / `dif` · `caf` | delete a feature (`dif`: only its own parts, e.g. exons) · replace its span |

| Command | |
|---|---|
| `:DnaView [linear\|circular]` | switch view; `!` draws a linear molecule as a circle anyway |
| `:DnaToggle` · `:DnaGoto {bp}` | toggle view · go to a base pair |
| `:DnaFeatures` / `:DnaFeatureInfo` | picker / details |
| `:DnaRecord [n]` | choose a record in a multi-record file |
| `:DnaTopology linear\|circular` | set topology (written to LOCUS, or as `[topology=circular]` in FASTA) |
| `:DnaInfo` · `:DnaRaw` · `:DnaCalibrate` | summary · raw file text · circle aspect test |
| `:DnaTracks [n\|all]` | track limit for this buffer (no argument: back to config) |
| `:DnaGroup [n]` | residues per group for this buffer; `0` = continuous with a ruler; no argument toggles |
| `:DnaSelect {start}..{stop}` / `:DnaSelect {label}` · `:DnaSelectClear` | select a range or a feature |
| `:DnaYank [top\|top35\|bottom\|bottom35] [register]` | yank the selection in that orientation (aliases `revcomp`, `complement`, `reverse`); no argument opens a picker |
| `:DnaFeatureAdd[!] [{s}..{e}]` · `:DnaFeatureEdit` · `:DnaFeatureDelete[!]` | author features (`!` on add opens the full editor) |
| `:DnaFeatureSave` · `:DnaLibrary` | save to / open your personal feature library |
| `:DnaUndo` · `:DnaRedo` | undo / redo feature changes |
| `:DnaInsert {pos} {residues}` · `:DnaDelete [{s}..{e}]` · `:DnaReplace [{s}..{e}] {residues}` | sequence edits (ranges default to the selection; `!` overrides `cut_policy = "refuse"`) |
| `:DnaRevComp [{s}..{e}\|all]` · `:DnaRotate [{bp}]` | reverse-complement a range or the record · move the origin |
| `:DnaEdit [{s}..{e}\|all]` · `:DnaEditReport` | edit residues as plain text · features affected by the last edit |
| `:DnaOrfs [min_aa\|off]` · `:DnaOrfAdd` | ORF overlay · ORF under the cursor → CDS |
| `:DnaFeatureDbImport {path} [name] [descriptions]` · `:DnaFeatureDbList` | import / list feature databases |
| `:DnaDetect [db…\|off]` · `:DnaDetectAccept[!] [fragments]` · `:DnaDetectReview` | propose, accept, review detected features |
| `:DnaFetch[!] {acc}[:{from}-{to}][:-]` · `:DnaFetch cancel` | open an NCBI record (`!` re-downloads and lifts the size limit) |
| `:DnaSearch[!] {query}` · `:DnaGene {symbol} [organism] [flank=N] [plus]` | search nuccore (`!`: protein) · fetch a gene region |
| `:DnaFetchCache [clear]` | list or clear cached NCBI records |
| `:DnaCloneAdd [name]` · `:DnaCloneList` · `:DnaCloneClear` | stage the record, the selection or a digest fragment as a part (`gA`) |
| `:DnaClone [parts…]` · `:DnaLigate [parts…]` | assemble, method inferred · join by ends |
| `:DnaGoldenGate [enzyme] [parts…]` · `:DnaGibson [parts…]` | one-pot Type IIS assembly · join by terminal homology |
| `:DnaAlign {query} [mode=] [strand=]` · `:DnaAlignReads {file\|glob}…` | align a query · place a run of reads |
| `:DnaOverview [bases-per-row]` · `:DnaRegion {start}..{stop}` | the whole record on one screen · go to a region and select it |
| `:DnaCompare {record}` · `:DnaAlignAnnotate` | compare two molecules · write the differences on as features |

The winbar shows the record, view, base pair under the cursor, the selection,
the codon under the cursor, and the feature under the cursor. Full reference:
`:help plasmid`.

</details>


## Features 

<details> 
<summary> Click to expand </summary>

### Strands and codons

With the cursor on an amino acid of a CDS or ORF translation, the three bases
of its codon are highlighted on both strands, and the winbar shows the codon
number, sense codon, residue and position (`codon 252/287 AAG→K
complement(3398..3400)`). Codons split by a line break, an intron or the origin
are highlighted in every piece.

Grouping by ten splits codons across gaps for most reading frames, and no fixed
group size avoids that for every gene. `:DnaGroup 0` (or `linear = { group = 0 }`)
shows the sequence continuously with a ruler numbering every 10th residue, as
SnapGene and Benchling do: no codon is split by a gap, and each residue sits
under the middle of its codon three columns from the next. Codons can still
break at the end of a line.

To copy a region, select it (`vaf`, `vif`, visual mode, `:DnaSelect`) and
choose the orientation: `top` (5'→3', as shown on the top row), `top35`
(reversed), `bottom` (5'→3', the reverse complement) or `bottom35` (3'→5', the
complement as shown on the bottom row). Visual `y` and `Y` are shortcuts for
`top` and `bottom`.

Yanks go to the system clipboard as well as Neovim's unnamed register, so they
can be pasted into other programs (BLAST, primer tools, a lab notebook). They
use the case the view shows (upper case for GenBank; FASTA soft-masking kept).
A named register (`"ay`) stays inside Neovim. The copy goes through Neovim's
clipboard provider (xclip, wl-copy, pbcopy, …) or, if there is none, through
the terminal as OSC 52, which works in most terminals and over SSH (tmux needs
`set -g set-clipboard on`). `:checkhealth plasmid` shows which is used;
`yank_clipboard = "auto"` follows your `'clipboard'` setting instead, and
`false` keeps yanks inside Neovim.

### Editing sequences

The view itself stays read-only; every edit goes through the record, so
annotations move with the bases they describe, and each edit is one undo step.
Use the keys above, the commands, or `:DnaEdit` for anything larger: it opens a
range as ordinary text where every Vim motion, search and macro works, and `:w`
applies only the differences, so features inside the region keep their bases.

How features follow an edit: deleted bases leave every feature (a feature that
loses all of them is removed); an insertion strictly inside a feature extends
it, while one at its first base or after its last shifts it; replacement bases
belong to the features that covered the bases they replace; reverse-complementing
a region flips features wholly inside it; moving the origin splits features that
cross it.

After each edit a message says what happened to the features, and
`:DnaEditReport` lists them: truncated, extended, destroyed, changed bases, and
for coding sequences FRAMESHIFT, new stop codons, lost start or stop codons and
the number of changed amino acids ("silent" if none). `/translation` is
rewritten for changed CDS features, or removed when the frame now hits a stop.
`edit.cut_policy` decides what happens when a deletion cuts into a feature:
`"truncate"` (report it), `"flag"` (also add a `/note` to the feature) or
`"refuse"` (cancel unless forced with `!`, or `:w!` in a `:DnaEdit` buffer).
Deleted residues go to the unnamed register as with `x`/`d`, so `xp` swaps two
bases. Letters are checked against the IUPAC codes (amino acids for proteins);
pasted GenBank ORIGIN text works.

Saving: GenBank files get a new ORIGIN and LOCUS length, and features that only
moved keep their original qualifier lines. A file opened from SnapGene is saved
as GenBank once its residues have changed (`:w name.gb`).

### Crowded regions

Each feature gets a score: a per-type value (`CDS` 100, `promoter` 80,
`primer_bind` 55, `misc_feature` 30, …), +10 if it has a real name, and +1000 if
you added or accepted it in this session. The highest scores are placed first,
nearest the sequence; ties go to the longer feature. Each line draws at most
`linear.max_tracks` (4) feature rows and summarises the rest. Placement is
decided per line, so a feature squeezed out of one crowded line still appears
wherever there is room.

### Your own features

Select a range (`vaf`, visual mode, or `:DnaSelect 3110..3175`), press `ga`, and
give a label, type and strand. `gE` opens any feature as text: the first line is
`TYPE LOCATION` in GenBank syntax (`complement(…)`, `join(…)`), then one
`/qualifier=value` per line. Every change is one undo step. `:DnaFeatureSave`
copies a feature into your library (`stdpath("data")/plasmid/library.gb`),
which `:DnaDetect` always searches, so an element you annotate once is found in
your other plasmids.

### ORFs

`go` shows open reading frames in all six frames: the longest ORF for each stop
codon, ATG starts, ≥ 75 aa by default (SnapGene's default), using the record's
`/transl_table` or else code 11 for circular and code 1 for linear records.
ORFs that end at an annotated CDS's stop codon are hidden. Each ORF shows its
translation; `ga` on one turns it into a CDS.

### Feature detection

No database ships with the plugin. To use pLannotate's SnapGene-derived feature
set (GPL-3, [mmcguffi/pLannotate](https://github.com/mmcguffi/pLannotate)):

1. Download `plannotate-databases-v2.tar.gz` from the pLannotate releases and unpack it.
2. In Neovim: `:DnaFeatureDbImport /path/to/BLAST_dbs`.
   This reads `snapgene.*` with `blastdbcmd` and the descriptions from the
   bundled SQLite file (via `sqlite3`, or Python's `sqlite3` module); about 1,370
   references are imported in well under a second.
3. Open a plasmid and run `:DnaDetect`.

Without BLAST+, export the sequences on another machine
(`blastdbcmd -db BLAST_dbs/snapgene -entry all > snapgene.fasta`) and run
`:DnaFeatureDbImport snapgene.fasta snapgene descriptions.csv`. Any FASTA file
(optionally with a CSV of id → name, type, description) or any GenBank/`.dna`
file or folder (every annotated feature ≥ 8 bp becomes a reference) can be
imported the same way.

Search: references ≥ 20 bp are found by shared 12-mers and banded
Smith-Waterman alignment on both strands, across the origin of circular
records; shorter ones (primer sites, tags) by exact search. Full-length hits
need ≥ 95 % identity over ≥ 95 % of the reference; partial matches (≥ 50 % of
the reference and ≥ 50 bp) are proposed as "(fragment)". Proposals are drawn in
italics and change nothing until accepted (`ga`, `:DnaDetectAccept`, `!` for all
full-length, `! fragments` for everything). Accepted features carry the
reference's description and a note such as
`detected by plasmid.nvim: snapgene/AmpR, 100.0% identity, 100% of the reference`.

### Search

`/` asks for a motif and `?` searches backwards. IUPAC codes match the bases
they stand for (`R` matches A or G). A sequence letter matches only if every
base it could be is allowed, so an `N` in the sequence matches only `N`. Both
strands are searched, palindromic sites are reported once, and circular records
are searched across the origin. Add a number for mismatches: `GGTCTC 1`.
Prefix `aa:` to search amino acids in the six-frame translation (`aa:MSIQHF`)
or just type the motif in a protein record. `n`/`N` step through the hits
(counts work), `gn` selects the current hit for `gy`, `c`, `:DnaEdit`, ….
Hits are highlighted on the strand they are on and marked on the circular map.
They are recomputed after every edit. `:DnaFind {motif} [mismatches]` does the
same from the command line; `:DnaFind off` clears.

### Restriction enzymes and digests

`gz` shows unique 6+ cutters. `:DnaEnzymes dual` shows enzymes cutting once or
twice, `:DnaEnzymes all` shows every cutting enzyme, and
`:DnaEnzymes EcoRI BamHI` shows only those. `gZ` lists every enzyme with its
number of sites and toggles enzymes in and out of the named set. In the linear
view a label's bar `▏` sits on the left edge of the base after the top-strand
cut, i.e. on the cut itself (`▕` when the label had to go left). A second row
under the complement marks the bottom-strand cut, so sticky ends show as a
stagger. Recognition sites are underlined. Enzymes that cut at the same place
share a label (`BsmBI, Esp3I`). `K` on a label shows the site, both cuts, the
overhang, the number of sites, and a warning when a Dam (GATC) or Dcm (CCWGG)
site overlaps a site of an enzyme known to be sensitive.

The bundled table holds 222 commercially available Type II enzymes, from
palindromic cutters through degenerate sites to Type IIS enzymes that cut
outside their site. Every site and cut position was checked against
Biopython's `Bio.Restriction`. For the full set, download REBASE's EMBOSS files
(`emboss_e.NNN`, and optionally `emboss_r.NNN` for suppliers and `emboss_s.NNN`
for supplier names) and run
`:DnaEnzymeImport emboss_e.NNN emboss_r.NNN emboss_s.NNN`. The import is stored
under `stdpath("data")/plasmid/` and replaces the bundled table until
`:DnaEnzymeImport clear`. With supplier data the `unique`/`dual`/`all` sets keep
to commercially available enzymes. REBASE data is not distributed with the
plugin.

`:DnaDigest EcoRI HindIII` (or just `:DnaDigest` after choosing enzymes with
`gZ`) opens a fragment table with each fragment's span, length, ends with
their overhang sequence (`EcoRI 5′ AATT`, `PstI 3′ TGCA`, `blunt`), and the
features it contains (`…` = cut). A gel lane with a ladder sits below it. `<CR>`
selects a fragment in the map. `e` opens it as a new linear GenBank record:
the top strand from one cut to the next, with whole features kept and cut ones
truncated by the same rules as editing. `:DnaGel EcoRI | EcoRI HindIII | uncut`
draws one lane per digest. The ladder is chosen automatically (1 kb or
100 bp), or set with `ladder=100bp`. Migration is drawn on a log scale; uncut
circular DNA is drawn at its linear size.

### Primers

`:DnaPrimerAdd GTAAAACGACGGCCAGT M13fwd` finds where a primer binds and
annotates it as a `primer_bind` feature, keeping the full oligo (in SnapGene's
primer record, or a `/note` in GenBank). A primer binds where its 3'-terminal
10 bases match exactly; the annealed region then extends 5' through mismatches
(a changed codon, say) as long as at least 5 matching bases follow, and
anything beyond that is a 5' tail. The rule reproduces SnapGene's own
annealing calls on 19 of the 20 primer sites in the test files. In visual mode
`gp` makes a primer from the selection in either direction.

Primers are drawn base by base over the template: forward primers above the
sequence, reverse below the complement, tails in lowercase, mismatched bases
highlighted, and the name at the 3' end (`›M13fwd`, `M13rev‹`). `K` on one
shows its sequence, where it anneals, Tm, binding-site count and, with a
facing primer, the product size and annealing temperature.

### Tm and annealing temperature

`:DnaTm GTAAAACGACGGCCAGT` gives the melting temperature for the current
polymerase; with two sequences it adds the recommended annealing temperature.
Options: `polymerase=taq`, `dmso=5`, `conc=200` (nM), `salt=50` (mM). With no
sequence it uses the selection. `:DnaPolymerase` lists the 22 NEB products
with their buffer salt and primer concentration.

The model is the one behind the NEB Tm Calculator: SantaLucia 1998
nearest-neighbour thermodynamics with the Owczarzy 2004 salt correction (or
Schildkraut–Lifson for Phusion, with NEB's conc/4 convention), per-product
buffer salts, and each product's Ta rule and cap. It is a Lua port of
[polymerase-tm](https://github.com/Dioskurides/polymerase-tm) (MIT; its licence
is included as `docs/polymerase-tm-LICENSE`), checked
against NEB's own Tm API: 5,940 Tm values across 22 products and 2,580 Ta
values across the 19 that have a Ta rule, no mismatches. Three corrections to polymerase-tm came out of that
comparison — four buffer salts, the Ta rounding procedure, and the published
form of the Owczarzy 2008 Mg²⁺ correction — so the default differs slightly
from polymerase-tm; `compat = "polymerase-tm"` reproduces it exactly.

### Designing primers and simulating PCR

Select a region and run `:DnaPrimerDesign`. Candidates vary in length and may
start up to `slack` bases outside the selection (the product grows; the
selection stays inside), which is what makes a Tm match possible at a GC-rich
boundary. Each primer is scored on its distance from the target Tm, GC
content, a 3'-terminal G/C, homopolymer runs and whether it binds the template
only once; pairs are scored on the Tm difference between them. Options:
`name=`, `tm=`, `len=18-30`, `slack=`, `dtm=`, and `fwd=EcoRI` / `rev=BamHI`
(or a sequence) to add 5' tails. Choosing a pair adds both primers in one undo
step. When nothing fits, the message says what to change and what Tm that side
can actually reach.

`:DnaPCR M13fwd M13rev` (names or sequences; the record's primers by default)
finds the product on a linear or circular template, including across the
origin, and opens it as a new linear record: 5' tails and primer mismatches
are written into the product, the template's features come along (truncated
where the product cuts them) and both primers are re-annotated full length.
Products were checked against [pydna](https://github.com/pydna-group/pydna):
232 products identical, and the 8 primer pairs pydna says cannot amplify give
no product here either.

### Cloning

Stage the parts, then assemble them. A part can be a whole record, the current
selection, a fragment of a digest table (`:DnaCloneAdd 2`), a PCR product, an
open buffer or a file:

```vim
:DnaCloneAdd insert      " in a sequence buffer: the selection, or the record
:DnaCloneAdd             " in a digest table: the fragment under the cursor
:DnaGoldenGate BsmBI     " or :DnaGibson, :DnaLigate, :DnaClone
```

`:DnaClone` infers the method: Golden Gate when every part carries a site for
one Type IIS enzyme, Gibson when parts share terminal homology of at least
`cloning.min_homology` (20 bp), otherwise ligation. Sticky and blunt ends join
in either orientation, so a blunt pair reports both products and offers a
picker; a set with distinct overhangs or distinct homology arms has one answer.

The product opens as a new unsaved record named after its parts. Features come
across from every part — a feature a cut had split is re-joined when the
assembly puts it back together — and each part and junction is annotated. A
report split lists the parts in order with their ends, every junction with its
overhang or overlap, the features carried, and warnings: Golden Gate parts
without a site, overhangs that are palindromic, unpaired or used more than
twice, a product that stayed linear, parts left over. **When nothing assembles
the same report lists every ordered pair of ends and why each was refused** —
a 5′ overhang meeting a 3′ overhang, overhangs that are not complementary, a
blunt end meeting a 4-base end — with what to change.

To design primers for an assembly, select the region and name the destination:

```vim
:DnaPrimerDesign name=tet into=pBR_dest enzyme=BsmBI     " Golden Gate tails
:DnaPrimerDesign name=tet into=pVec_linear homology=25   " Gibson tails
```

The tails are taken from the destination: its terminal sequence for Gibson, or
the enzyme site, spacer and the fusion overhangs its own cut leaves for Golden
Gate (`overhangs=AATG/GCTT` to set them yourself). A site for that enzyme
inside the region to amplify is flagged, and the pair you choose is amplified
and assembled to confirm it closes.

### Large records

A bacterial genome is a quarter of a million rows of sequence. Above
`linear.window.threshold` (300 kb) the linear view works out how long it is
and then lays out only the rows a window is showing, building them as you
scroll and dropping them again behind you. Nothing else changes: the
coordinates, the feature tracks, selection, search, enzymes, primers,
alignment and editing all behave as they do on a plasmid, and
`window_spec` checks row for row that the windowed layout is the ordinary
one built a piece at a time.

```vim
:DnaOverview        " the whole record on one screen, one row per span
:DnaRegion 150000..151000
```

`:DnaOverview` (`gm`) trades the residues for the shape of the record: one
row per span — 10 kb a row for a 230 kb chromosome, 200 kb a row for *E.
coli* — with a strip of cells showing which strands carry features where,
coloured by the most important feature in each cell, the per-strand counts,
and as many feature names as fit. The span is a round number picked so the
record fits the window, and the header says what it is. `<CR>` or a click
jumps to the region under the cursor and `q` puts the window back.

![The whole genome at 200 kb per row](docs/screenshots/overview_genome.png)

![Yeast chromosome I at 10 kb per row](docs/screenshots/overview_chromosome.png)

### Alignment

`:DnaAlign {query}` places a query on the open record and opens an alignment
view. The query can be residues typed inline, a file, an open buffer, a staged
cloning part, a digest fragment or `selection`. Both strands are tried, and on
a circular record the alignment may run through the origin; `mode=` picks
between `semi` (the query placed inside the record, the default), `global` and
`local`, and `strand=plus|minus` pins one strand.

```vim
:DnaAlign sanger/well_A.fa          " one read
:DnaAlignReads sanger/*.fa          " a whole run, as a table
:DnaCompare designed.gb             " two molecules, end to end
```

The view shows the alignment in blocks — reference row, match row, query row,
with mismatches and gaps highlighted — above a difference list that names the
feature, the codon and the effect (`1255 G>A  tet codon 390 TGG>TGA W390*
nonsense`), a table of the features hit with their verdicts (silent, missense,
nonsense, frameshift with the premature stop, in-frame indel), and a summary
line. `]d` and `[d` walk the differences, `<CR>` jumps to one in the record,
and `:DnaAlignAnnotate` writes them onto the record as `variation` features in
one undo step.

`:DnaAlignReads` places every record of every file (globs welcome) and lists
them by position with strand, span, identity, difference count and coverage;
reads that do not place are listed with the reason. `:DnaCompare` rotates a
circular molecule into phase before comparing, so a different origin is not
reported as a difference — the header says what rotation and strand it had to
apply.

### Fetching from NCBI

`:DnaFetch L09137.2` downloads a record from NCBI and opens it in a new,
unsaved buffer named after the accession (`L09137.2.gb`, `NP_000509.1.gp`) in
`fetch.dir` or the current directory; `:w` keeps it. Nucleotide accessions are
fetched as GenBank with all parts (so contig records arrive with their
sequence), protein accessions (`NP_`, `XP_`, UniProt such as `P69905`, PDB
chains such as `4HHB_A`) as GenPept. A range fetches only that region,
numbered from 1: `:DnaFetch NC_000913.3:190-2800`, and a trailing `:-` gives
the reverse complement. GI numbers and accessions without a version are
resolved to the current version, which is shown in the message and the file
name.

`:DnaSearch pBR322[title] AND complete sequence[title]` runs an Entrez query
(full Entrez syntax) and lists hits with length, topology and title;
`:DnaSearch!` searches protein. `:DnaGene BRCA1` looks the symbol up in NCBI
Gene (organism defaults to `fetch.organism`, otherwise human) and fetches the
gene's genomic region on its own strand, so minus-strand genes come
reverse-complemented and read 5'→3' like the mRNA; `plus` keeps the reference
orientation and `flank=500` adds bases on both sides.

On a feature with cross-references (`/protein_id`, `/db_xref="GeneID:…"`,
UniProt), the `K` window shows a `links` line and `gF` opens the linked record,
or for a gene its region.

Records are cached under `stdpath("cache")/plasmid/ncbi` for `cache_days`
(30); `:DnaFetchCache clear` empties it. Before downloading, the record's
size is checked against `ncbi.max_length` (5 Mb) so a stray chromosome
accession does not pull hundreds of megabytes; fetch a sub-range, or use
`:DnaFetch!`, which re-downloads and lifts the limit. Requests go through `curl`
asynchronously (Neovim stays usable; `:DnaFetch cancel` stops them) and are
spaced to NCBI's limit of 3 per second, or 10 with an API key. Nothing
identifying is sent unless you configure `ncbi.email` / `ncbi.api_key`;
NCBI asks heavy users to do both (see its E-utilities usage guidelines). The
API key is passed to curl on standard input, not on the command line.

</details>

## Configuration

<details> 
<summary> Click to expand </summary>

Defaults, abridged (see `lua/plasmid/config.lua` and `:help plasmid-config`):

```lua
require("plasmid").setup({
  default_view = "linear",        -- or "circular" / "auto" (circular for circular records)
  hidden_types = { source = true },
  label_qualifiers = { "label", "gene", "product", "locus_tag", "standard_name", "allele", "note" },
  priority = {
    types = { CDS = 100, gene = 95, promoter = 80, rep_origin = 80, primer_bind = 55,
              misc_feature = 30, --[[ ... ]] },
    default = 40, labelled_bonus = 10, pinned_bonus = 1000,
  },
  linear = { line_width = "auto", show_complement = true, show_translation = true,
             color_bases = false, max_tracks = 4, group = 10 },   -- group = 0: continuous + ruler
  circular = { aspect = 1.0, max_tracks = 6, label_width = 22 },
  orfs = { min_aa = 75, starts = "atg", table = nil, hide_annotated = true },
  detect = { min_identity = 0.95, min_coverage = 0.95, fragments = true,
             min_fragment_coverage = 0.5, min_fragment_length = 50, nested_same_type = false },
  yank_clipboard = true,          -- also copy plain yanks to the system clipboard; "auto" | false
  edit = { cut_policy = "truncate", update_translations = true },  -- "flag" | "refuse"
  enzymes = { default_set = "unique", min_site = 6, lanes = 4, commercial_only = true },
  gel = { ladder = "auto", height = 24 },      -- "1kb" | "100bp"
  primers = {
    polymerase = "q5", dmso = 0, conc_nM = nil, salt_mM = nil,
    show = true, lanes = 3,
    binding = { anchor = 10, max_mismatches = 3, resume = 5, min_anneal = 10 },
    design = { len = { 18, 30 }, tm = { 57, 68 }, max_dtm = 3, gc = { 0.35, 0.65 },
               max_run = 4, slack = 8, pairs = 8 },
  },
  cloning = { min_homology = 20, max_homology = 200, enzyme = "BsaI",
              annotate_parts = true, annotate_junctions = true, max_solutions = 50 },
  -- above `threshold` residues the linear view lays out only what a window shows
  linear = { window = { enabled = "auto", threshold = 300000, margin = 2, cache_blocks = 512 } },
  overview = { span = "auto", min_span = 10, labels = true, counts = true },
  align = { match = 2, mismatch = -3, gap_open = -5, gap_extend = -2, mode = "semi",
            min_identity = 0.9, min_coverage = 0.5, width = 0, annotate_type = "variation" },
  mouse = { enabled = true, drag_select = true, double_click_feature = true, primer_pair = true },
  ncbi = { api_key = nil, email = nil, max_length = 5000000, cache = true, cache_days = 30 },
  fetch = { dir = nil, search_results = 20, gene_flank = 0, organism = nil },  -- dir nil: cwd
  featuredb_dir = nil,            -- stdpath("data")/plasmid/featuredb
  library_path = nil,             -- stdpath("data")/plasmid/library.gb
  keymaps = { toggle_view = "<Tab>", feature_outer = "af", feature_inner = "if",
              toggle_orfs = "go", add_feature = "ga", edit_feature = "gE",
              search_forward = "/", search_backward = "?", search_next = "n", search_prev = "N",
              select_hit = "gn", toggle_enzymes = "gz", pick_enzymes = "gZ",
              add_primer = "gp", follow_link = "gF", clone_add = "gA", --[[ ... ]] },
  colors = { CDS = "#f2b134", promoter = "#4fbf75", --[[ ... ]] },
})
```

All highlight groups are `default` links and can be overridden: `PlasmidTitle`,
`PlasmidRuler`, `PlasmidComplement`, `PlasmidBackbone`, `PlasmidTick`,
`PlasmidLeader`, `PlasmidBase{A,C,G,T}`, `PlasmidSelection`, `PlasmidProposed`,
`PlasmidSearch`, `PlasmidSearchCurrent`, `PlasmidEnzyme`, `PlasmidEnzymeUnique`,
`PlasmidEnzymeCut`, `PlasmidEnzymeSite`, `PlasmidGel{Ladder,Band,BandStrong,Label,Well}`,
`PlasmidPrimer`, `PlasmidPrimerTail`, `PlasmidPrimerMismatch`,
`PlasmidPrimerArrow`, `PlasmidPrimerLabel`,
and `PlasmidFeature<type>` (e.g. `PlasmidFeatureCDS`).

## Design

The `Record` model is the source of truth. The buffer is a read-only
projection (`nomodifiable`, enforced) rendered from it, so annotations cannot
drift from the residues they describe. Every change goes through a transaction
layer (`model/transaction.lua`): sequence edits remap feature locations on the
location syntax tree (`model/edit.lua`), so a feature keeps its form (`join`,
`complement`, fuzzy ends), and each change records an exact inverse for undo.
ORFs and detection proposals are overlays: drawn like features but never part
of the record until accepted.

```
lua/plasmid/
  model/   record, location (grammar AST), store (ResidueStore), alphabet, transaction, edit (remapping)
  io/      genbank, fasta, snapgene, xml, format detection
  bio/     translation, NCBI codon tables, ORF finder, feature detection (seeds + banded SW), edit impact
  view/    linear, circular, braille canvas, track packing and priority, highlights
  buffer.lua  ui.lua  select.lua  author.lua  seqedit.lua  scratch.lua  codon.lua
  orfs.lua  detectui.lua  featuredb.lua  search.lua  enzymes.lua  digest.lua
  primers.lua  primerdesign.lua  pcr.lua  ncbi.lua (E-utilities client)  fetch.lua
  bio/assembly.lua (fragment ends, ligation, Golden Gate, Gibson)  clone.lua  mouse.lua
  bio/align.lua (affine-gap pairwise, seeded band)  bio/variant.lua  alignment.lua
  view/overview.lua (whole-record map)  overview.lua
  init.lua  config.lua  health.lua
```

Residues sit behind a `ResidueStore` interface. The current backend is a plain
Lua string: comfortable for editing plasmid-sized records (≤ ~100 kb) and for
viewing records up to a few Mb (a 5 Mb genome opens in < 1 s). Highlights come
from a decoration provider that touches only on-screen rows.

</details>

## Tests
<details> 
<summary> Click to expand </summary>

```
nvim -l tests/run.lua                        # all 403 specs, headless, no dependencies
nvim -l tests/run.lua detect                 # filter by spec name
PLASMID_TEST_GROUP=0 nvim -l tests/run.lua   # the whole suite in the continuous layout
```

- **Oracle:** every record, feature location and qualifier in the corpus is
  compared with Biopython's interpretation (`tests/oracle/expected/*.json`,
  regenerated by `tests/oracle/dump_biopython.py`; Python is a development tool
  only, never a runtime dependency).
- **Round-trip:** each fixture is written back byte-identical. A forced
  canonical rewrite re-parses to the same model. Features written to SnapGene
  files are re-read by Biopython.
- **Translation:** every annotated `/translation` in the corpus is reproduced
  from the sequence.
- **Location grammar:** exact format/parse round-trip over the corpus and 2,000
  random ASTs.
- **ORFs:** identical to an independent brute-force six-frame scan
  (`tests/oracle/dump_orfs.py`) for 6 records × 4 settings (> 3,000 ORFs).
  Of 20 annotated CDS that start with ATG and meet the length limit, 17 are
  found exactly and 3 lie inside a longer ORF with an upstream in-frame ATG.
  Six frames of the 230 kb yeast chromosome take ~30 ms.
- **Codons:** the cursor is put on every residue drawn for pBR322, phiX174
  and the human β-globin region (4,010 residues, including 2 codons split by
  the phiX174 origin and 7 split by introns). Each time, the highlighted bases
  must spell the reported codon and translate to the residue under the cursor,
  and the residues assembled by codon number must reproduce each CDS's
  `/translation` (20 CDS).
- **Sequence editing:** a property test tags every base of six records
  (circular, origin-crossing, spliced, SnapGene) and applies random insertions,
  deletions, replacements, region inversions and origin moves; after each edit
  every feature must cover exactly the tagged bases an independent model of the
  rules predicts, and undoing everything must give back the original objects
  and a byte-identical GenBank file. The suite runs one seed (347 edits);
  `PLASMID_EDIT_SEEDS=200` ran 69,866 edits with 2.9 million feature checks.
  Writing this test found four remapping bugs before release. Unit tests cover
  every rule, the impact report on pBR322's genes (in-frame insertion,
  frameshift, new stop, silent change, lost start codon), the cut policies,
  every key and command, `:DnaEdit` (diff, stale refusal, forced writes) and
  saving. `tests/tools/edit_oracle.py` writes randomly edited GenBank files and
  has Biopython re-read them: 40 files, 2,422 features, all agreeing.
- **Clipboard:** unit tests cover every `yank_clipboard` setting, named and
  clipboard registers, case handling and the no-clipboard warning.
  `tests/tools/osc52_check.py` runs Neovim in a real pseudo-terminal with no
  clipboard tool installed and decodes the OSC 52 sequences it writes: all four
  strand orientations of pBR322 10..40 match Biopython's.
- **Layouts:** CI runs the whole suite twice, with groups of 10 and with the
  continuous layout. For group sizes 10, 0, 3, 5 and 7 every base pair of
  pBR322 maps to a column showing that residue and back. In the continuous
  layout every multiple of 10 on the ruler sits over its residue, no codon
  within a line has a gap (713 checked), and neighbouring residues are always
  three columns apart (706 pairs).
- **Detection:** a synthetic plasmid with planted exact, reverse,
  origin-crossing, mismatched, gapped, partial and 15 bp elements is recovered
  exactly, and an unrelated 20 kb sequence gives no hits. 6 kb against 1,510
  references takes ~0.1 s. During development, results on pLannotate's 10
  control plasmids were compared with pLannotate 1.2.5: 83 of 84 full-length
  SnapGene-set hits were reproduced (the miss, "FRT (minimal)", lies inside a
  longer FRT that is reported).
- **Search:** 120 random IUPAC motifs (with N, degenerate codes and 0–2
  mismatches), each searched on a linear and a circular record (240 searches),
  give the same 24,052 hits as a brute-force scan with its own IUPAC table. Exact search of the 230 kb yeast
  chromosome takes ~7 ms, one mismatch ~25 ms.
- **Restriction enzymes:** `tests/tools/verify_enzymes.py` checks every entry of
  the bundled table against Biopython's `Bio.Restriction` (site, both cuts).
  Introduced errors are caught. Cut positions for all 222 enzymes on pUC19,
  pBR322 and phiX174, linear and circular, are identical to Biopython's
  `search()`: 1,332 cases, 5,527 cuts, including sites across the origin and
  Type IIS cuts outside the site. Scanning all enzymes over 230 kb takes
  ~0.25 s, using a 4-mer index that agrees with the plain scan on every enzyme.
- **REBASE import:** fixture files cover every skip rule (two-cut, no cut,
  bad site, duplicate, malformed). `tests/tools/check_rebase_roundtrip.py`
  writes all 1,088 enzymes in Biopython's dictionary back out as EMBOSS files
  and imports them: the 729 single-cut enzymes come back with identical sites,
  cuts and suppliers.
- **Digests and gel:** 360 random digests (1–3 enzymes, both topologies)
  give the fragments predicted from Biopython's cut lists. An extracted pBR322
  fragment keeps every whole feature with identical residues and saves as
  GenBank that reads back the same. Gel placement is checked for order, log
  spacing, co-migration and off-range sizes.
- **Tm / Ta:** every value is checked against three independent sources
  (`tests/tools/tm_oracle.py`, fixture `tests/oracle/expected/tm.json`):
  NEB's public Tm API (5,940 Tm values over 22 products and 2,580 Ta values
  over the 19 with a Ta rule, no mismatches), polymerase-tm in compatibility mode (400,846 Tm and 111,628 Ta
  values, no mismatches), and Biopython's `MeltingTemp` for the raw and
  salt-corrected values (agreement to 3e-13 °C). The comparison with NEB
  pinned down the buffer salts, the two-step Ta rounding and NEB's round-half-up
  behaviour; `tm_report.md` records what differs from polymerase-tm and why.
- **Primer binding:** the annealing rule reproduces SnapGene's own calls on 19
  of 20 primer sites in the test files (the exception stores annealing that no
  longer matches its template). 600 constructed sites — random tails, internal
  mismatches, both strands, across the origin — are found with exactly the
  right span, tail and mismatch list.
- **PCR:** 240 random primer pairs on pUC19, pBR322 and phiX174, linear and
  circular, with random 5' tails and some pairs facing away from each other,
  are compared with pydna's `pcr()`: 232 identical products, 8 agreed
  no-product cases. (On 7 of them the plugin also lists longer products from
  additional weak binding sites, which pydna's `limit=13` rule excludes.)
- **Primer design:** every ranked pair must satisfy the stated constraints and
  is fed back through the PCR engine to confirm it amplifies the selected
  region exactly once.
- **UI (real headless Neovim):** goto/cursor round-trip on every base pair,
  feature lookups, circular angle mapping, text objects, yanks, authoring with
  undo, ORF and detection overlays, writes, conversion, multi-record files, and
  timing on a 230 kb chromosome and a 5 Mb genome.
- **Large records:** the windowed layout is compared with the ordinary one
  row by row — text, row kinds, spans and position queries — across the whole
  fixture corpus (71 files, 78 records, 172,227 rows), and on a fetched corpus
  of real genomes (*E. coli* K-12, yeast chromosome IV, human cytomegalovirus,
  a multi-record phage file, an 8.7 Mb *Streptomyces* FASTA) at five window
  positions each, in both sequence layouts, at three widths, with the enzyme,
  primer, ORF, search and selection overlays: **109 cases, 71,858 rows, no
  differences**. `budget_spec` generates a 400 kb record and fails if opening
  it starts costing what laying out every row used to.
  `docs/large_records_report.md` has the corpus and the measurements.
- **Alignment:** the engine is compared with Biopython 1.88's
  `PairwiseAligner` over 1,097 generated cases — five mode variants, eight
  scoring schemes, nucleotide and protein — and **all 1,097 agree on score**,
  with every case that has a unique optimum reproduced column for column. 288
  simulated reads (119 of them across the origin) are placed within 2 bp on the
  right strand at every error rate tested. The comparison found three engine
  defects, all fixed. `docs/alignment_report.md` has the corpus, the
  read-recovery table and the timings.
- **Cloning:** the assembly engine is compared with
  [pydna](https://github.com/pydna-group/pydna) 5.5.16 over a generated corpus
  of 468 cases — 262 restriction–ligation, 125 Gibson, 81 Golden Gate,
  including 48 designed to fail. Products are compared in a canonical form
  (least rotation of the sequence or its reverse complement), along with
  fragment lengths and order, the number of assemblies and the no-product
  cases: **468 of 468 agree**. The comparison found two engine bugs, both
  fixed: an overhang census that miscounted when a part was supplied on the
  other strand, and a solution cap that bounded what was reported but not what
  was searched. `docs/cloning_report.md` has the corpus, the conventions where
  pydna differs, and the performance measurements.
- **NCBI client:** 30 cases (accession parsing, GI numbers, sub-ranges on
  both strands, UniProt ids, CON records, the size guard, searches, gene
  coordinates, NCBI's error formats, retries with backoff, rate-limit
  spacing, cancellation, the cache, and that an API key never reaches a
  command line, file or history) run offline on recorded NCBI responses and a
  local test server; with `PLASMID_LIVE=1` the recorded requests are repeated
  against NCBI and must give the same results. 44 accession forms are accepted
  and 17 malformed ones rejected. See `docs/ncbi_client_report.md`.
- **Real NCBI records:** 48 records in 41 files (plasmids, viral and organelle
  genomes, mRNAs, ncRNA, contig and WGS master records, segmented sets,
  GenPept, a PDB chain, a protein with `bond()` sites; 13,454 features) are
  compared with Biopython: sequences, headers, feature types,
  3,063 / 3,071 locations (the other 8 are `bond()` and `one-of()` forms,
  which Biopython represents differently) and 12,467 qualifier values agree. The check pinned six
  bugs, all fixed: `bond()` locations, contig/WGS records without sequence,
  out-of-range features on the circular map, edits on sequence-less records,
  multi-line qualifiers with embedded quotes, and translation of ambiguous
  codons (every codon containing an IUPAC ambiguity code was translated to `X`
  even where it resolves, e.g. `CTN` → L). See `docs/ncbi_records_report.md`.
- `tests/tools/tui_capture.py` drives Neovim in a pseudo-terminal and saves
  screenshots (`--script` replays a key sequence), for checking rendering by eye.

SnapGene fixtures come from Biopython's test suite (Biopython licence) and from
sgffp (MIT, github.com/merv1n34k/sgffp). `.dna` files are checked against
Biopython; `.prot` and `.rna` against sgffp (`tests/oracle/dump_sgffp.py`),
since Biopython does not read those. Protein FASTA fixtures are NCBI/UniProt
entries for haemoglobin α and β and pBR322's TetA and β-lactamase. The TetA entry
is checked against the plugin's own translation of pBR322. No pLannotate data is
included in this repository.

NCBI fixtures (besides the fetch corpus above): pUC19 (L09137), pBR322 (J01749), phiX174
(NC_001422, origin-spanning genes), pPCP1 (NC_005816), human β-globin region
(U01317), yeast chr I (NC_001133), three mRNA records, a multi-FASTA, and a
synthetic record covering the rest of the grammar.

</details>

## Known issues

- Neovim 0.12.2 headless (`nvim --headless`, `nvim -l`) can crash in
  `:redraw!` after a message wider than `'columns'`, independent of this plugin.
  plasmid.nvim truncates its own messages to the available width.
- Methylation warnings are partial: only enzymes whose Dam/Dcm/CpG sensitivity
  could be confirmed from NEB's published information carry a flag (41 of 222).
  No warning does not mean an enzyme is insensitive.
- `:DnaEnzymeImport` follows REBASE's documented EMBOSS layout (as Biopython's
  converter reads it). It is tested on fixtures and on files generated from
  Biopython's enzyme data, not yet on a downloaded REBASE release.
- The gel is a schematic: sizes are placed on a log scale with ~15 % size steps
  per row, and band intensity and topology effects are not modelled.
- The windowed layout knows how long the buffer is, so it plans every block
  when the record opens: that is 264 ms of the 1.1 s *E. coli* takes, and it
  grows with the record. Resizing the window re-plans (615 ms on that genome),
  which is the one operation that still feels like work. The circular view
  builds in full whatever the record's size.
- The overview's picture saturates on a dense genome: at 200 kb a row nearly
  every cell has features on both strands, so what you read is where the
  coding stops. The per-strand counts and the names carry the detail.
- Alignment is pairwise: there is no multiple alignment, no consensus from a
  set of reads, and no trace (`.ab1`) reading. Above `align.max_cells` (4e6)
  the engine aligns inside a seeded band, which needs the two sequences to
  share 11-mers; unrelated long sequences are declined rather than aligned
  badly. Where differences cluster, the list is one description among several
  of equal score.
- Assembly enumerates orderings. With many interchangeable ends (a one-enzyme
  digest of several fragments) only the first `cloning.max_solutions` (50)
  orderings are reported, and exhaustive enumeration is practical to about 8
  such fragments; a design with distinct overhangs or homology arms prunes to
  one answer. Cutting dominates on large records: roughly 60 ms per fragment
  when a few hundred kb is cut into many pieces.
- Gibson joins are exact terminal matches; mismatched or partially annealing
  arms, chew-back and fill-in are not modelled, and ends of different lengths
  never anneal.
- NCBI occasionally answers with HTTP 429 or a transient "Search Backend
  failed"; requests are retried with backoff, but a busy period can still
  surface an error. Try again, or set an API key.
- Neovim 0.12.2 can crash ("double free or corruption") when the experimental
  `nvim__redraw({ buf = … })` targets a buffer while another window is current;
  plasmid.nvim repaints through extmarks instead.

## Roadmap

- **Alignment, further:** reading Sanger traces, and a consensus across a set
  of reads.
- **Large records, further:** a circular view that scales to a genome, and
  laying out the blocks lazily enough that opening does not have to plan them
  all.
- **Cloning, further:** partial digests, chew-back/fill-in, and writing an
  assembly plan (which parts, which reaction, which checks) to a file.
- **Later:** writing edited primers back to SnapGene files, and a piece-table
  store for editing records larger than ~100 kb.

## License

MIT, see [LICENSE](LICENSE). Bundled third-party material keeps its own
licence: the Tm model is ported from polymerase-tm (MIT,
`docs/polymerase-tm-LICENSE`); test fixtures from Biopython's test suite
(Biopython licence) and sgffp (MIT); NCBI records are public data.
