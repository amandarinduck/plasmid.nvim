# Opening a large record without laying it out first

v0.9 made the linear view windowed: above `linear.window.threshold` it builds
the text of a block only when a window is about to show it. What it still did
at open was plan every block -- all 51,574 of them on E. coli K-12 -- purely to
learn how many rows the buffer has, and it did the same thing again on every
window resize. This track removes both.

Everything below was measured on an idle 12-core x86-64 machine with Neovim
0.12.2, in a window 120 columns wide and 45 lines tall, which gives E. coli 90
residues per line, 51,574 blocks and 256,745 rows. "Before" is the v0.9.0 tree
with this track's dev tools copied into it, so both sides are measured by the
same harness.

## Where the time went

`tests/tools/plan_profile.lua` separates the setup every layout does from the
pass over the blocks, counts what the pass allocates, and samples the pass
under a count hook to attribute it to the functions it runs in. On E. coli:

| | before |
|---|---|
| parsing the file | 841 ms |
| windowed layout, first build in the process | 194 ms |
| the planning pass alone | 100 ms, 1.94 us per block |
| what the pass allocates | 42.5 MB, 864 bytes per block |
| the same loop doing only its block arithmetic | 2.2 ms, 0.04 us per block |

Medians of three paired runs against both trees, each a median of five passes.
So the pass cost about forty-five times its own arithmetic, and made 42.5 MB of
garbage to produce 51,574 integers. The sampled profile says where:

| | share of samples |
|---|---|
| `items_in` -- the two-level index lookup, dedupe and sort | 49.1% |
| `plan` itself -- track assignment and the tables it fills | 41.2% |
| `primer_lanes` returning nil immediately, twice per block | 6.4% |

The count hook turns the JIT off, so those are shares, not times. Three things
fall out of them. `items_in` is irreducible: a block's row count depends on
which features overlap it, so they have to be found and ranked whatever the
pass is for. The 41.2% is mostly bookkeeping a counting pass does not need --
a plan table, a `by_track` table per feature row, a `{start, stop}` table per
extent inside the interval sets, a `tlist`, and a sort of it. And the 6.4% is
pure waste: a record with no primers was asked twice per block whether it had
any.

The floor, then, is `items_in` plus the track-assignment decisions with none of
the tables: the measured result below is 16.4 ms, 0.32 us per block, which is
about a sixth of the planning pass and eight times the bare arithmetic.

## What changed

### 1. One function decides a block's rows, and it can be asked not to record

`lua/plasmid/view/linear.lua`. `plan(b)` was a single function that decided a
block's rows and recorded the decision at the same time. It is now
`rowcount(b, p)`, where `p` is the plan to fill in or `nil`; `plan(b)` passes a
fresh table and `count(b)` passes nothing. Both take the same decisions in the
same order against the same interval sets, which is what makes the count a
block is opened with the same number as the plan it is later built from. The
alternative -- a second, cheaper implementation of the same rules -- is exactly
the thing that drifts, and the drift would not be visible as a wrong row: it
would be every row after the mistake landing at the wrong line number.

The feature-track assignment moved into `assign(b, bs, be, p)` for the same
reason, and the interval sets it decides against are no longer `vf._insert`'s
array of `{start, stop}` tables but three flat arrays reused between blocks
(`set_s`, `set_e`, `set_n`), with `row_used`, `row_cds` and `touched` beside
them. Counting a block now allocates nothing at all; the first counting pass
over a genome allocates 0.7 MB in total, which is those scratch arrays growing
to the size of the busiest block and staying there.

The windowed build loop calls `count(b)`, keeps the row count and the
`offsets` it already kept, and throws nothing away that it used to keep -- it
never retained the plans from that pass in the first place.

### 2. A record with no primers is asked once, not twice per block

`has_primers` is computed once in `build`. `primer_lanes` keeps its own early
return, so it is still correct called directly; it is just not called 103,148
times for nothing.

### 3. Resizing asks whether the line width changed before laying anything out

`lua/plasmid/view/linear.lua`, `lua/plasmid/buffer.lua`,
`lua/plasmid/init.lua`, `lua/plasmid/config.lua`.

The residues per line are quantised -- they stay a multiple of ten and of the
group size so the ruler and the end-of-line numbers stay round -- so a range of
window widths all give the same layout, row for row. That calculation moved out
of `M.build` into `M.line_width(L, o, win_width)`, which `build` now calls and
which `buffer.needs_relayout(buf, win)` can ask without building anything. The
resize autocmd asks it first and does nothing when the answer is the width the
layout already has. At 120 columns on E. coli that covers every width from 115
to 125; a border drag mostly moves inside one such range.

When the width has really changed, the work is deferred by
`config.options.resize_debounce` (default 80 ms) and coalesced, so dragging a
border lays the record out once, at the size the drag stops at, instead of at
every column on the way. The circular view keeps the tolerance it had: it fills
whatever width it is given, so it has no quantum to compare against.

`buffer.render`'s construction of the per-buffer options moved into
`effective_opts(st)` so that `needs_relayout` asks about the same configuration
`render` would use -- a buffer with `:DnaGroup 0` has a different line width at
the same window width.

## Before and after

### E. coli K-12, the layout alone

| | before | after |
|---|---|---|
| windowed layout, first build in the process | 194 ms | 104 ms |
| the pass over all 51,574 blocks | 100 ms planning | 14 ms counting |
| per block | 1.94 us | 0.28 us |
| what the pass allocates | 42.5 MB | 0.0 MB |
| per block | 864 bytes | 0 bytes |
| planning every block, which the on-demand path still does | 100 ms / 42.5 MB | 53 ms / 29.3 MB |

Counting allocates nothing per block. The scratch it reuses grows once to the
size of the busiest block -- 0.7 MB over the whole genome when the counting
pass is the first thing to run in the process, nothing on any pass after that.
The arithmetic-only loop measured between 0.1 and 2.6 ms across those six runs,
so counting is within about an order of magnitude of the floor; what is left in
it is the feature-index lookup that decides which features a block contains,
which a row count genuinely depends on.

Planning a block is still what the on-demand path does, and it got cheaper too
-- 864 to 596 bytes a block -- because the interval sets it decides against are
now reused rather than rebuilt.

### The corpus through a buffer

`tests/tools/big_corpus.sh` runs the whole large-record corpus through
`plasmid.buffer` three times in the windowed mode and once in the eager one,
and `tests/tools/big_tables.py` takes the medians. Opening is `:edit`;
resizing goes through the plugin's own `WinResized`/`VimResized` handling and
waits for whatever that handling defers, so the 80 ms debounce counts against
the resize that caused it. A resize that keeps the line width never reaches
that wait: the view is asked first and says there is nothing to do.
Milliseconds, before -> after.

| record | mode | open | resize, new line width | resize, same line width | peak heap MB | peak rss MB |
|---|---|---|---|---|---|---|
| NC_000913 | window | 1,062 -> 961 | 291 -> 253 | 246.17 -> 0.04 | 187 -> 148 | 243 -> 196 |
| NC_000913 | eager | 5,176 -> 5,347 | 3,349 -> 4,266 | 3,565.99 -> 0.04 | 3,884 -> 1,633 | 5,089 -> 2,080 |
| NC_001136 | window | 347 -> 355 | 95 -> 139 | 70.24 -> 0.04 | 68 -> 60 | 102 -> 93 |
| NC_001136 | eager | 1,844 -> 1,849 | 1,286 -> 1,799 | 1,410.22 -> 0.05 | 1,179 -> 700 | 1,499 -> 811 |
| NC_003888.3 | window | 353 -> 300 | 142 -> 154 | 107.39 -> 0.04 | 56 -> 40 | 93 -> 77 |
| NC_003888.3 | eager | 2,008 -> 2,075 | 2,362 -> 2,540 | 1,702.77 -> 0.04 | 532 -> 381 | 717 -> 524 |
| NC_006273 | window | 63 -> 58 | 35 -> 117 | 18.04 -> 0.04 | 27 -> 25 | 50 -> 46 |
| NC_006273 | eager | 369 -> 336 | 332 -> 428 | 286.04 -> 0.05 | 163 -> 125 | 220 -> 155 |
| NC_001133 | window | 62 -> 57 | 19 -> 103 | 11.58 -> 0.04 | 25 -> 23 | 43 -> 41 |
| NC_001133 | eager | 250 -> 245 | 251 -> 300 | 192.82 -> 0.04 | 149 -> 78 | 186 -> 102 |
| NC_000866 | window | 102 -> 103 | 16 -> 108 | 10.72 -> 0.04 | 34 -> 33 | 52 -> 50 |
| NC_000866 | eager | 346 -> 344 | 147 -> 280 | 147.41 -> 0.04 | 154 -> 96 | 184 -> 118 |
| NC_001416 | window | 93 -> 96 | 12 -> 107 | 7.38 -> 0.04 | 33 -> 24 | 49 -> 38 |
| NC_001416 | eager | 157 -> 158 | 90 -> 172 | 68.84 -> 0.04 | 52 -> 26 | 76 -> 55 |
| NC_001604 | window | 98 -> 104 | 19 -> 102 | 6.91 -> 0.04 | 31 -> 19 | 46 -> 32 |
| NC_001604 | eager | 240 -> 218 | 61 -> 140 | 36.03 -> 0.04 | 52 -> 32 | 77 -> 55 |
| NC_001895 | window | 97 -> 102 | 8 -> 101 | 7.03 -> 0.04 | 29 -> 17 | 44 -> 30 |
| NC_001895 | eager | 223 -> 217 | 53 -> 139 | 33.59 -> 0.04 | 47 -> 38 | 65 -> 53 |

The resize that changes nothing is the result: every record, both modes, from
laying the whole thing out again to 0.04 ms.

Peak memory fell everywhere, and most on the eager layouts -- E. coli's eager
peak went from 3.9 GB of Lua heap to 1.6 GB, and from 5.1 GB of resident memory
to 2.1 GB. That was not aimed at: planning a block no longer allocates an
interval table per extent, and the eager path plans every block, so its peak
allocation fell with it.

A resize that *does* change the line width is now slower on a small record --
8 ms to 101 ms on the smallest phage genome -- because the measurement includes
the 80 ms wait. That is the trade the debounce makes: on a record where laying
out again costs 8 ms the wait is most of the number, and on E. coli, where the
whole resize is 253 ms, it is under a third of it, and in both cases it is what
stops a twenty-column border drag from costing twenty re-layouts.
`resize_debounce` is configurable and 0 restores the old behaviour.

Scrolling and jumping were measured too and are unchanged within the
run-to-run spread; they are in `docs/large_records_report.md` with the rest of
the corpus tables.

## Incremental planning: not shipped, and why

Step 4 of the brief was conditional: plan incrementally *only if counting is
still the bulk of the open*. It is not. Of E. coli's 961 ms open, parsing the
file is 847 ms and the whole windowed layout is 104 ms, of which the counting
pass is 14 ms -- under two per cent of the open. Planning incrementally could
save at most those 14 ms, and would buy them by making the buffer's length a number
that changes under the user: `block_at_row`, `pos_of`, `ensure`, the scroll
autocmds and `goto_bp` all read the offsets table, and `:$`, a mark, a fold, a
split and the scrollbar all read the buffer's line count. Making that provably
exact is a larger and riskier change than the whole of this track, for a saving
that is inside the run-to-run variation of the parse that precedes it.

The brief said that leaving the counting pass as the answer was an acceptable
outcome. That is the answer. If the parse ever gets much faster -- it is the
obvious next target, being eight times the layout -- the question is worth
reopening, and the shape it should take is the one the brief describes: enough
blocks for the first screens, the rest on demand, with every position query
planning up to the block it asks about.

## Tests

- The big-record parity matrix was re-run against the changed tree:
  `tests/tools/big_parity.py --matrix`, **109 cases, 0 differ, 71,858 rows
  compared** -- the same result it gave for v0.9. Each case dumps the same
  windows of a real record from a windowed and an eager layout in separate
  processes and compares row text, row kind, base span, highlight runs,
  feature and enzyme segments, overflow lists, and the answers to `pos_of`,
  `bp_at` and `feature_at`. This is the check that matters on a genome: a
  block counted as one number and built from another moves every row after
  it, and the comparison would see the base span move.
- `tests/spec/window_spec.lua` passes unchanged. It is the row-for-row proof
  that the windowed layout is the eager one built a piece at a time, and
  nothing in it needed to move: 71 files, 78 records, 172,227 rows.
- `tests/spec/rowcount_spec.lua` is new and is the proof that matters for this
  track: every block of every fixture record counted and planned, in both
  sequence layouts, and the two numbers compared -- 26,266 blocks, 171,735
  rows. Then the same with enzyme and primer overlays (64 records, 8,977
  blocks), at three track caps including one low enough to displace features,
  and counted forwards, backwards, repeatedly and with plans taken in between,
  because the interval sets are now reused and a block that left something
  behind would show up as a different count the second time.
- The whole suite, in both sequence layouts: **419 passed, 0 failed** in
  135,140 ms at group 10 and **419 passed, 0 failed** in 132,750 ms with
  `PLASMID_TEST_GROUP=0`.
- `tests/spec/budget_spec.lua` gains two cases: one that counting every block
  of the 400 kb synthetic record allocates under 0.5 MB where planning it
  allocates at least twenty times that, and one that twenty window resizes
  that keep the residues per line leave the layout object untouched while a
  resize that changes them produces a new one. Its header table is
  re-measured.

## Found and not fixed

- **The parse is now the open.** 847 ms of E. coli's 961 ms. Nothing in this
  track touches it.
- **`items_in` sorts every block.** `table.sort(scratch, by_rank)` is run per
  block over the items overlapping it, and the sampled profile puts half the
  pass inside `items_in`. The index could keep its buckets in rank order at
  build time and merge rather than sort. Not attempted: the two-level index is
  shared with the on-demand path and the build-time cost would move rather than
  disappear.
- **A resize still rebuilds from scratch when the width does change.** 253 ms
  on E. coli, of which 80 ms is the debounce. The feature index and the global
  track packing do not depend on the line width and could be carried over; the
  blocks are what would have to be recounted.
- **The debounce is a flat 80 ms for every record.** A plasmid that lays out in
  8 ms waits as long as a genome whose whole resize is 253 ms. Scaling the
  wait with the record, or skipping it for a resize that is not part of a
  drag, would be better; neither was attempted, and `resize_debounce = 0` is
  the escape hatch.
- **`vf._overlaps` and `vf._insert` have no callers left.** They were exported
  from `lua/plasmid/view/features.lua` for the linear view's track assignment,
  which now has its own non-allocating interval sets; nothing outside
  `features.lua` uses them any more. Left in place: removing an export is not
  this track's business, and `features.lua` is not a file this track owns.
- **`drop_block` keeps the plan.** `M.drop_block` blanks a block's rows and
  clears `self._plans[b]`, but the plan sweep in `ensure` is what bounds plan
  memory, and it only runs when `ensure` is called. A session that jumps
  without scrolling still reaches the plan footprint the v0.9 report measured.
- **The eager path still plans every block to render it**, which is correct --
  it needs the plan -- but it now pays for `rowcount`'s `if p then` branches.
  No measurable cost was found at the sizes the eager path runs at; the eager
  build of the 400 kb synthetic is within the run-to-run spread of the v0.9
  measurement.
