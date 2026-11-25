# Large records

The linear view lays out a record a block at a time above
`linear.window.threshold` (300,000 residues): opening a file counts each
block's rows, which is what tells the buffer how many rows it has, and builds
the text of a block only when a window is about to show it. Resizing the window
asks whether the residues per line have changed before it lays anything out
again. `tests/spec/window_spec.lua` proves the windowed layout is the eager one
built a piece at a time, over the whole fixture corpus, and
`tests/spec/rowcount_spec.lua` proves the count a block is opened with is the
count it is later built with. This report is the part the fixtures cannot
reach: real records of a megabase and more, the overlays on them, and what they
cost.

Everything below was measured on an idle 12-core x86-64 machine with Neovim
0.12.2, in a window 120 columns wide and 45 lines tall, which gives the view 90
residues per line. The measurements were taken again for v0.10, against the
same corpus and the same harness, and where a number is given as "before ->
after" the before side is the v0.9.0 tree measured by this tree's tools.

## The corpus

Large records are not committed -- the E. coli GenBank file alone is 11.5 MB --
so `tests/tools/fetch_big.sh` rebuilds them from NCBI E-utilities into `big/`,
which git ignores. The sizes and feature counts are what the parser reports,
not what the records' headers claim.

| record | file | bp | features | topology | line width | rows | blocks |
|---|---|---|---|---|---|---|---|
| NC_000913 | NC_000913.gbk | 4,641,652 | 9,285 | circular | 90 | 256,745 | 51,574 |
| NC_001136 | NC_001136.gbk | 1,531,933 | 2,464 | linear | 90 | 92,835 | 17,022 |
| NC_003888.3 | NC_003888.fasta | 8,667,507 | 0 | linear | 90 | 288,921 | 96,306 |
| NC_006273 | NC_006273.gbk | 235,646 | 441 | linear | 90 | 16,872 | 2,619 |
| NC_001133 | NC_001133.gb | 230,218 | 325 | linear | 90 | 12,955 | 2,558 |
| NC_000866 | phages.gbk | 168,903 | 591 | linear | 90 | 9,935 | 1,877 |
| NC_001416 | phages.gbk | 48,502 | 285 | linear | 90 | 3,702 | 539 |
| NC_001604 | phages.gbk | 39,937 | 162 | linear | 90 | 2,454 | 444 |
| NC_001895 | phages.gbk | 33,593 | 217 | linear | 90 | 2,190 | 374 |

The four phage genomes are one file: a multi-record GenBank, so that switching
records inside a buffer is covered. `NC_003888.3` is a FASTA with no features
at all, which is the case where every row of a block is sequence and complement
and nothing else. `NC_001133` (yeast chromosome I) is already a committed
fixture and is included for continuity with the measurements the previous track
reported.

File sizes: NC_000913.gbk 11.5 MB, NC_003888.fasta 8.8 MB, NC_001136.gbk
4.2 MB, phages.gbk 1.0 MB, NC_006273.gbk 0.5 MB; 25.9 MB in total.

## Method

Eight dev tools, none of them part of the suite:

- `tests/tools/fetch_big.sh` fetches the corpus and skips what is already
  there.
- `tests/tools/big_dump.lua` builds one layout of one record -- windowed or
  eager, one mode per process -- and writes selected windows as JSON: each
  row's text, kind, base span, highlight runs, feature segments, enzyme
  segments and overflow list, plus the answers to `pos_of`, `bp_at`,
  `block_at_row` and `feature_at`.
- `tests/tools/big_parity.py` runs both modes and compares the two dumps row by
  row. Two processes rather than one because the eager layout of E. coli K-12
  is over half a gigabyte of Lua heap, and the windowed layout of the same
  record has to be measured without it in memory.
- `tests/tools/big_budget.lua` opens a record through `plasmid.buffer` in a
  headless Neovim and measures opening, scrolling, jumping and resizing.
  Scrolling is a real viewport move (`winrestview`, then the fill the autocmds
  would do), not a direct call into the layout.
- `tests/tools/big_cursor.lua` jumps to 309 base pairs of a record through the
  buffer and checks each one lands where it was asked to.
- `tests/tools/big_corpus.sh` runs `big_budget.lua` over the whole corpus, one
  process per record per mode, three runs of the windowed mode and one of the
  eager one, appending a line of JSON each time.
- `tests/tools/big_tables.py` reduces repeated runs of the same record to
  their median and prints the tables below, so none of them is a single noisy
  run and none of the numbers in them was typed by hand.
- `tests/tools/plan_profile.lua` separates the setup every layout does from the
  pass over the blocks, times the planning and counting passes against each
  other in turn, counts what each allocates, and samples the pass under a count
  hook to attribute it to the functions it runs in.

Resizing is measured through the plugin's own `WinResized`/`VimResized`
handling rather than by calling the renderer, and the measurement waits for
whatever that handling defers, so a layout that arrives 80 ms after the resize
counts against the resize that caused it.

### What is compared, and where

Each record is compared at five window positions: the start, the middle, the
end, the region where the most feature extents begin -- where the track packer
has the most to do -- and, for a circular record, both sides of the origin. A
window is 24 blocks, which came to between 120 and 1,308 rows per case
depending on how much is drawn there.

Before any of that, and before a single row has been built, the windowed layout
is asked for `pos_of`, `bp_at` and `block_at_row` at about 4,000 positions
spread over the whole record, and the answers are compared with the eager
layout's. A plan that counted a block's rows wrongly shows up there rather than
only inside a window; the comparison also asserts that those queries built
nothing.

### Axes

Every record is compared across both sequence layouts (group 10 and the
continuous ruler layout), two window widths either side of the default, the
track limit at 1 and at `all`, and each overlay on its own so a difference
names one of them: restriction enzymes (unique cutters, and every site for
records under 2 MB), primers, ORFs, a search, a selection. Two more cases are
about the windowed path itself. `churn` builds the windows, scrolls past them
with a cache too small to hold them, and compares the rows that come back:
a block that is laid out, evicted and laid out again has to come back the same.
`wrap` adds a feature across the origin of a circular record, because real
genomes rarely annotate one and the seam is otherwise never drawn -- with it,
the feature is drawn on 15 rows at the end of E. coli and 29 at the start. A
linear record has no seam, so there it adds the feature to the last 4,000 bases
instead and compares the final block, which is the short one (24 rows of the
added feature on NC_006273).

That last part was a correction. The first run of this matrix asked for the
origin window on every record, and on the seven linear ones it got no window
back and compared no rows at all -- seven cases reporting a pass having
examined nothing. `tests/tools/big_parity.py` now fails a case that compares
zero rows, which is what turns that class of hole into a result rather than a
silence.

## Parity result

**109 cases, 0 differences, 71,858 rows compared**, every case comparing
between 120 and 1,308 rows. Per file: NC_000913 13 cases, NC_001136 13,
NC_003888 13, NC_006273 14 (it is small enough for the every-site enzyme
overlay), phages 56 (14 cases across each of its 4 records).

No divergence between the windowed and the eager layout was found, on any
record, at any window position, under any overlay, in either sequence layout.

The whole matrix was run again for v0.10, after the open stopped planning
every block and started counting them, and gave the same result: 109 cases,
0 differ, 71,858 rows compared. That is the check that matters for the
counting pass on real genomes -- a block counted as one number and built from
another would move every row after it, and the matrix compares the base span
of every row it looks at as well as the text.

The comparison was checked against injected faults before being trusted, since
a harness that never fails proves nothing. Three windowed-only faults were
introduced into `lua/plasmid/view/linear.lua` one at a time. Only the first was
caught on the first attempt; the other two went unnoticed, and each miss turned
out to be the fault itself being vacuous rather than the harness being blind.
The eviction fault blanked a block's trailing row, which is empty anyway, and
was then not reached because the scroll that was meant to evict the compared
blocks never visited them; the `pos_of` fault removed an offset that is zero in
group 10 with no overlays, so there was nothing to get wrong. Reformulated --
blanking a block's *first* row, and running in the continuous layout where the
ruler makes that offset non-zero -- and with the eviction pass made to evict
what it compares, all three are caught:

- `plan` counting one row too few for one block in a thousand -- caught on the
  buffer's total length, on the row span of every window after it, and on the
  whole-record position sample (`bp 755: eager row 42, windowed row 41`);
- `pos_of` ignoring the rows a block has above its sequence row -- caught by
  the position sample in the continuous layout, where the ruler is such a row;
- a block coming back with a blanked row after eviction -- caught by the
  `churn` case.

Two further injections were made into code that *both* layouts run (a label
centred one column off, a row count changed in `plan` under a condition that
never fired) and, correctly, produced no difference: the two layouts share the
planner and the renderer, so only the windowed path's own arithmetic can
diverge.

Jumping through the buffer was checked separately on five records with
`tests/tools/big_cursor.lua` for v0.9 and not run again for v0.10, where the
parity matrix above covers the same arithmetic: 309 jumps each, 1,545 in
total. Every jump landed on the base it asked for, with the row present in the
buffer and identical to the layout's, and every feature the track rows
reported under the cursor's column really covered that base (298 to 718 such
resolutions per record). No problems.

## Measurements

Opening is `:edit` through the plugin: read, parse, lay out, and fill the first
screen. Heap is Lua heap; RSS is the whole process. "Peak" is the highest value
seen across the whole run, which includes the resize, where two layouts of the
record are briefly alive at once.

| record | mode | open ms | heap after open MB | rss after open MB | peak heap MB | peak rss MB |
|---|---|---|---|---|---|---|
| NC_000913 | window | 961 | 112 | 130 | 148 | 196 |
| NC_000913 | eager | 5,347 | 710 | 862 | 1,633 | 2,080 |
| NC_001136 | window | 355 | 21 | 49 | 60 | 93 |
| NC_001136 | eager | 1,849 | 255 | 321 | 700 | 811 |
| NC_003888.3 | window | 300 | 24 | 47 | 40 | 77 |
| NC_003888.3 | eager | 2,075 | 118 | 207 | 381 | 524 |
| NC_006273 | window | 58 | 7 | 19 | 25 | 46 |
| NC_006273 | eager | 336 | 42 | 64 | 125 | 155 |
| NC_001133 | window | 57 | 7 | 18 | 23 | 41 |
| NC_001133 | eager | 245 | 30 | 48 | 78 | 102 |
| NC_000866 | window | 103 | 11 | 24 | 33 | 50 |
| NC_000866 | eager | 344 | 40 | 60 | 96 | 118 |
| NC_001416 | window | 96 | 9 | 23 | 24 | 38 |
| NC_001416 | eager | 158 | 17 | 31 | 26 | 55 |
| NC_001604 | window | 104 | 10 | 23 | 19 | 32 |
| NC_001604 | eager | 218 | 20 | 39 | 32 | 55 |
| NC_001895 | window | 102 | 10 | 23 | 17 | 30 |
| NC_001895 | eager | 217 | 19 | 39 | 38 | 53 |

Latency, in milliseconds, as median / maximum. The two resize columns are a
resize that lands on a different number of residues per line, which needs the
record laid out again, and one that does not, which needs nothing. The first
includes the 80 ms the plugin waits for a border drag to stop before it lays
anything out; the second never reaches that wait, because the view is asked
first and says there is nothing to do.

| record | mode | one page | long scroll, per page | jump to a base pair | resize, new line width | resize, same line width |
|---|---|---|---|---|---|---|
| NC_000913 | window | 0.9 / 1.5 | 0.9 / 16.4 | 6.8 / 7.4 | 253 | 0.04 / 0.11 |
| NC_000913 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 4,266 | 0.04 / 0.12 |
| NC_001136 | window | 0.6 / 0.9 | 0.8 / 4.4 | 6.7 / 19.3 | 139 | 0.04 / 0.13 |
| NC_001136 | eager | 0.0 / 0.0 | 0.0 / 0.1 | 0.0 / 0.0 | 1,799 | 0.05 / 0.12 |
| NC_003888.3 | window | 0.3 / 0.5 | 0.4 / 6.2 | 3.4 / 5.9 | 154 | 0.04 / 0.11 |
| NC_003888.3 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 2,540 | 0.04 / 0.10 |
| NC_006273 | window | 0.6 / 1.0 | 0.7 / 3.9 | 6.8 / 16.7 | 117 | 0.04 / 0.12 |
| NC_006273 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 428 | 0.05 / 0.08 |
| NC_001133 | window | 0.7 / 1.5 | 0.5 / 4.9 | 5.6 / 15.6 | 103 | 0.04 / 0.10 |
| NC_001133 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 300 | 0.04 / 0.07 |
| NC_000866 | window | 0.8 / 1.1 | 0.7 / 3.6 | 4.6 / 19.9 | 108 | 0.04 / 0.12 |
| NC_000866 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 280 | 0.04 / 0.15 |
| NC_001416 | window | 0.7 / 1.0 | 0.0 / 1.0 | 0.0 / 5.4 | 107 | 0.04 / 0.11 |
| NC_001416 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 172 | 0.04 / 0.12 |
| NC_001604 | window | 1.0 / 3.4 | 0.0 / 0.9 | 0.0 / 2.4 | 102 | 0.04 / 0.14 |
| NC_001604 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 140 | 0.04 / 0.11 |
| NC_001895 | window | 0.8 / 1.4 | 0.0 / 2.5 | 0.0 / 2.4 | 101 | 0.04 / 0.10 |
| NC_001895 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 139 | 0.04 / 0.10 |

Scrolling was sampled as 300 consecutive pages of 42 rows from a third of the
way into the record (100 pages for the phage records), and jumping as 20 jumps
to pseudo-random base pairs. "One page" is a single page-down from a cold
position, ten times. The resize columns are one resize four columns narrower
and back, twenty times, for the same-width case, and one resize twenty columns
narrower for the other.

### What changed in v0.10

| record | mode | open ms | resize, new line width | resize, same line width | peak heap MB | peak rss MB |
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

Four things to read out of that table. The resize that changes nothing is the
headline: it went from laying the whole record out again -- 246 ms on E. coli,
3.6 s for its eager layout -- to 0.04 ms, on every record in both modes. Peak
memory fell everywhere, most of all for the eager layouts, because planning a
block no longer allocates an interval table per extent: E. coli's eager layout
peaked at 3.9 GB of Lua heap and now peaks at 1.6 GB. E. coli's windowed open
is about a tenth faster, which is the counting pass replacing the planning one.

And a resize that *does* change the line width is now slower on a small record
-- 8 ms to 101 ms on the smallest phage -- because the measurement includes the
80 ms the plugin now waits for a drag to finish. On a record where laying out
again costs 8 ms that wait is most of the number; on E. coli, where the whole
resize is 253 ms, it is under a third of it, and it is what stops a
twenty-column drag from costing twenty re-layouts. `resize_debounce` is
configurable, and 0 restores the old behaviour of laying out on every resize
event that changes the width.

### Where the open goes

E. coli again, this time the layout on its own rather than through a buffer:
the median of three paired runs of `tests/tools/plan_profile.lua`, each a
median of five passes, against both trees.

| | v0.9.0 | v0.10 |
|---|---|---|
| parsing the file | 841 ms | 847 ms |
| the windowed layout, first build in the process | 194 ms | 104 ms |
| the pass over all 51,574 blocks | 100 ms planning | 14 ms counting |
| per block | 1.94 us | 0.28 us |
| what the pass allocates | 42.5 MB | 0.0 MB |
| per block | 864 bytes | 0 bytes |
| planning every block, where the on-demand path still does | 100 ms / 42.5 MB | 53 ms / 29.3 MB |

The counting pass allocates nothing per block. The scratch it reuses -- the
interval sets the track assignment decides against -- grows once to the size of
the busiest block: 0.7 MB over the whole genome when counting is the first
thing to run, and nothing on any pass after that. The same loop doing only its
block arithmetic and no work at all measured between 0.1 and 2.6 ms across
those six runs, so the counting pass is within about an order of magnitude of
the floor, and what remains in it is the feature-index lookup that decides
which features a block contains -- which a row count genuinely depends on.

Planning is still what the on-demand path does for a block a window is about to
show, and it got cheaper too, from 864 to 596 bytes a block, because the
interval sets it decides against are now reused rather than rebuilt.

The sampled profile of the v0.9.0 planning pass, under a count hook (which
turns the JIT off, so these are shares of samples and not of time): 49.1% in
`items_in`, 41.2% in `plan` itself, 6.4% in the two `primer_lanes` calls that
returned nothing because the record has no primers. The first is irreducible,
the second is the bookkeeping the counting pass does not do, and the third is
now asked once per record instead of twice per block.

The eager layout's scroll and jump columns are zero because there is nothing
left to do: it paid for every row when the file was opened. That is the whole
trade. What the windowed view buys with those sub-millisecond-to-20 ms
latencies is the first two columns of the table above: E. coli opens in 1.0 s
instead of 5.3 s, holding 112 MB of heap instead of 710 MB, and never goes
above 196 MB of resident memory against 2.1 GB.

The worst single page in the whole sweep was 16.4 ms and the worst single jump
19.9 ms, both on records where a jump lands somewhere with nothing laid out and
has to build a screenful plus its margin before the cursor can be placed. Both
are inside a 30 Hz frame. Median page latency is under 1 ms on every record,
and median jump latency between 0 and 6.8 ms.

Resizing is the one place the windowed view is still slow in absolute terms
when it has work to do: 253 ms for E. coli, of which 80 ms is the debounce and
the rest is laying the record out at the new line width. It is 17x faster than
the eager re-layout of the same record (4.3 s), and a resize that keeps the
line width -- which most resizes are -- costs 0.04 ms.

### What stays laid out

After opening, 300 pages of scrolling, 20 jumps and twenty-two resizes, this is
what the view is still holding:

| record | blocks built | rows laid out | of total rows |
|---|---|---|---|
| NC_000913 | 131 | 514 | 256,745 |
| NC_001136 | 117 | 535 | 92,835 |
| NC_003888.3 | 213 | 426 | 288,921 |
| NC_006273 | 94 | 540 | 16,872 |
| NC_001133 | 112 | 536 | 12,955 |
| NC_000866 | 113 | 523 | 9,935 |
| NC_001416 | 97 | 402 | 3,702 |
| NC_001604 | 82 | 402 | 2,454 |
| NC_001895 | 82 | 367 | 2,190 |

Between 367 and 540 rows of text, out of up to 288,921, and between 82 and 213
blocks. These are higher than the v0.9 report's 142 to 182 rows because the
sweep now ends with twenty-two resizes rather than two, and each one fills the
window it lands in; the cache is still doing what it was built to do.

## The engine changes

### v0.9: the plan cache is swept

`lua/plasmid/view/linear.lua`, in `M.ensure`, plus a counter in `M.plan_of` and
`M.drop_block`.

A position query -- jumping to a base pair, asking what is under the cursor --
plans the block it lands in and keeps the plan, so that the next query about
the same block is free. The eviction in `ensure` only forgot blocks that had
been *built*, so plans made by position queries alone were never dropped.
Measured on E. coli: a position query in every one of its 51,574 blocks left
51,574 plans and grew the heap by 27.9 MB, from 14.2 MB to 42.1 MB, and
scrolling did not take any of it back. With the sweep, one screen of scrolling
brings it to 26 plans and 12.1 MB. Re-measured on the v0.10 tree, holding a
plan for every one of the 51,574 blocks still costs 29.8 MB; what plans contain
did not change.

### v0.10: counting instead of planning, and a resize that asks first

`plan(b)` decided a block's rows and recorded the decision at the same time. It
is now `rowcount(b, p)`, where `p` is the plan to fill in or `nil`; `plan(b)`
passes a fresh table and `count(b)` passes nothing, so both take the same
decisions in the same order and the count a block is opened with is the count
it is later built with by construction rather than by coincidence. The
feature-track assignment moved into `assign` for the same reason, and the
interval sets it decides against are flat arrays reused between blocks instead
of a `{start, stop}` table per extent. The windowed build calls `count`.

`M.line_width(L, o, win_width)` is the residues-per-line calculation lifted out
of `M.build`, so `buffer.needs_relayout(buf, win)` can ask what a window width
would give without building anything. The resize autocmd in
`lua/plasmid/init.lua` asks it first and does nothing when the answer is the
width the layout already has, and defers the rest by
`config.options.resize_debounce` (default 80 ms) so a border drag lays the
record out once, at the width it stops at.

The parity matrix above was run against the engine with both changes in it and
reported no differences. `tests/spec/window_spec.lua` passes unchanged, and
`tests/spec/rowcount_spec.lua` is the counting pass's own spec: it counts and
plans every block of every fixture record and insists the two numbers agree.

## Budgets in the suite

`tests/spec/budget_spec.lua` keeps this from regressing without the network or
a download: it generates a 400 kb record with 4,000 overlapping features --
joins, complements, enough overlap to fill the feature tracks -- and measures
the same things at a tenth of E. coli's scale. Its limits, and what they were
measured at, are in the comment at the top of that file. The short version:

- the windowed layout of that record takes 400 ms and 15 MB of heap at most
  (measured 32 ms and 2.6 MB), and must lay out exactly zero rows;
- counting every block allocates under 0.5 MB (measured 0.00 MB), and planning
  them must cost at least 20x that allocation and 1.5x the time (the lowest
  seen over six runs was 341x and 3.0x);
- the eager layout of the same record must cost at least 3x the time and 5x the
  heap (the lowest seen over six runs was 12.6x and 17.0x) -- a ratio, so a slow
  CI runner cannot flake it, and the assertion that would fail if the view went
  back to laying everything out up front;
- one screen under 150 ms, 200 pages under 3 s, 400 position queries under
  300 ms and building 200 cold rows under 1.5 s;
- a position query in every block, then one screen of scrolling, must leave
  under an eighth of the plans and under half the heap they cost;
- twenty window resizes that keep the residues per line must leave the layout
  object untouched and cost under 200 ms in total (measured 0.9 ms), and one
  that changes them must produce a new layout with a different line width;
- through a real buffer: opening under 3 s and 60 MB, 200 pages under 3 s, 40
  jumps under 2 s with the cursor reporting the right base every time, two
  resizes under 4 s.

The spec takes about 3 s and runs in both sequence layouts.
`tests/spec/rowcount_spec.lua` is the correctness half of the same change: it
counts and plans every block of every fixture record, in both sequence layouts
and with the enzyme and primer overlays, and insists the two numbers agree.

## What this does not cover

- **The circular view.** Only the linear view is windowed. `:DnaView circular`
  on a megabase record still builds the whole thing, and it is not measured
  here: it is measured in `docs/ring_report.md`, which covers what that view
  draws above `circular.atlas.threshold` and what it costs (322 ms for the
  first build of E. coli K-12 at 120x40, 51 ms for the next).
- **Editing a large record.** Every measurement is read-only viewing. An edit
  rebuilds the layout, which costs what a resize that changes the line width
  costs -- 253 ms on E. coli, debounce included -- and nothing here measures a
  sequence of edits.
- **The parse.** At 847 ms of E. coli's 961 ms open it is now the whole of the
  open; nothing in this report touches it.
- **Highlighting.** The decoration provider only runs for rows on screen, and
  no measurement here forces a redraw: `:redraw` is not safe in a headless
  Neovim 0.12.2, so the cost of drawing highlights on a genome is untested.
  What is tested is that the highlight *runs* the provider would apply are
  identical between the two layouts.
- **Peak heap during a long jump sequence.** The plan sweep runs on the next
  `ensure`, so a user who jumps many times without scrolling can still hold a
  plan for every block before anything is swept -- 29.8 MB on E. coli,
  re-measured for v0.10.
- **Records above 10 Mb.** The largest measured here is 8.7 Mb. Nothing on the
  scale of a human chromosome was fetched or tried, so where the current
  approach stops being comfortable is not known.
- **Real terminals.** Every measurement is headless. Scrolling is a viewport
  move plus the fill the autocmds would do, which is what the plugin does, but
  not what a terminal's own redraw costs.
