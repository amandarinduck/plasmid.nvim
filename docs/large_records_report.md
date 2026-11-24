# Large records

The linear view lays out a record a block at a time above
`linear.window.threshold` (300,000 residues): opening a file plans every block,
which is what tells the buffer how many rows it has, and builds the text of a
block only when a window is about to show it. `tests/spec/window_spec.lua`
proves that layout is the eager one built a piece at a time, over the whole
fixture corpus. This report is the part the fixtures cannot reach: real records
of a megabase and more, the overlays on them, and what they cost.

Everything below was measured on an idle 12-core x86-64 machine with Neovim
0.12.2, in a window 120 columns wide and 45 lines tall, which gives the view 90
residues per line.

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

Five dev tools, none of them part of the suite:

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
`tests/tools/big_cursor.lua`: 309 jumps each, 1,545 in total. Every jump landed
on the base it asked for, with the row present in the buffer and identical to
the layout's, and every feature the track rows reported under the cursor's
column really covered that base (298 to 718 such resolutions per record). No
problems.

## Measurements

Opening is `:edit` through the plugin: read, parse, lay out, and fill the first
screen. Heap is Lua heap; RSS is the whole process. "Peak" is the highest value
seen across the whole run, which includes the resize, where two layouts of the
record are briefly alive at once.

| record | mode | open ms | heap after open MB | rss after open MB | peak heap MB | peak rss MB |
|---|---|---|---|---|---|---|
| NC_000913 | window | 1078 | 128 | 165 | 154 | 233 |
| NC_000913 | eager | 5103 | 724 | 877 | 2634 | 3113 |
| NC_001136 | window | 369 | 37 | 52 | 60 | 105 |
| NC_001136 | eager | 1766 | 261 | 327 | 777 | 1112 |
| NC_003888.3 | window | 337 | 27 | 50 | 45 | 94 |
| NC_003888.3 | eager | 1947 | 118 | 206 | 438 | 597 |
| NC_006273 | window | 65 | 7 | 19 | 24 | 47 |
| NC_006273 | eager | 310 | 44 | 66 | 101 | 206 |
| NC_001133 | window | 63 | 5 | 19 | 21 | 41 |
| NC_001133 | eager | 245 | 31 | 49 | 132 | 161 |
| NC_000866 | window | 104 | 13 | 26 | 32 | 51 |
| NC_000866 | eager | 329 | 41 | 60 | 118 | 142 |
| NC_001416 | window | 91 | 10 | 23 | 20 | 38 |
| NC_001416 | eager | 150 | 17 | 31 | 27 | 67 |
| NC_001604 | window | 102 | 12 | 24 | 23 | 37 |
| NC_001604 | eager | 203 | 20 | 39 | 46 | 62 |
| NC_001895 | window | 98 | 12 | 24 | 21 | 34 |
| NC_001895 | eager | 199 | 19 | 39 | 40 | 57 |

Latency, in milliseconds, as median / maximum:

| record | mode | one page | long scroll, per page | jump to a base pair | resize re-render |
|---|---|---|---|---|---|
| NC_000913 | window | 0.9 / 16.1 | 0.7 / 1.1 | 7.6 / 36.7 | 615 |
| NC_000913 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 8817 |
| NC_001136 | window | 0.7 / 3.4 | 0.7 / 5.3 | 6.4 / 27.9 | 177 |
| NC_001136 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 3287 |
| NC_003888.3 | window | 0.3 / 0.5 | 0.4 / 2.9 | 3.5 / 9.3 | 255 |
| NC_003888.3 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 3855 |
| NC_006273 | window | 0.7 / 1.6 | 0.7 / 4.5 | 6.8 / 16.9 | 56 |
| NC_006273 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 629 |
| NC_001133 | window | 0.9 / 2.1 | 0.6 / 2.8 | 5.9 / 16.8 | 48 |
| NC_001133 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 408 |
| NC_000866 | window | 0.8 / 1.2 | 0.7 / 5.1 | 5.5 / 15.3 | 46 |
| NC_000866 | eager | 0.0 / 0.1 | 0.0 / 0.0 | 0.0 / 0.0 | 323 |
| NC_001416 | window | 0.7 / 2.3 | 0.0 / 0.9 | 0.0 / 5.2 | 40 |
| NC_001416 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 151 |
| NC_001604 | window | 1.0 / 1.7 | 0.0 / 0.8 | 0.0 / 2.3 | 19 |
| NC_001604 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 97 |
| NC_001895 | window | 0.8 / 2.2 | 0.0 / 1.0 | 0.0 / 2.2 | 17 |
| NC_001895 | eager | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.1 | 76 |

Scrolling was sampled as 300 consecutive pages of 42 rows from a third of the
way into the record (100 pages for the phage records), and jumping as 20 jumps
to pseudo-random base pairs. "One page" is a single page-down from a cold
position, ten times.

The eager layout's scroll and jump columns are zero because there is nothing
left to do: it paid for every row when the file was opened. That is the whole
trade. What the windowed view buys with those sub-millisecond-to-37 ms
latencies is the first two columns of the table above: E. coli opens in 1.1 s
instead of 5.1 s, holding 128 MB of heap instead of 724 MB, and never goes
above 233 MB of resident memory against 3.1 GB.

The worst single page in the whole sweep was 16.1 ms and the worst single jump
36.7 ms, both on E. coli: a jump lands somewhere with nothing laid out and has
to build a screenful plus its margin before the cursor can be placed. The page
is inside a 30 Hz frame; that worst jump is just outside one, and it is the
slowest thing the windowed view does. Median page latency is under 1 ms on
every record, and median jump latency between 0 and 7.6 ms.

Resizing is the one place the windowed view is slower in absolute terms than it
looks: 615 ms for E. coli, because a resize re-plans every block at the new line
width. It is still 14x faster than the eager re-layout of the same record
(8.8 s), and it is the operation a user performs least often.

### What stays laid out

After opening, 300 pages of scrolling, 20 jumps and two resizes, this is what
the view is still holding:

| record | blocks built | rows laid out | of total rows |
|---|---|---|---|
| NC_000913 | 45 | 170 | 256,745 |
| NC_001136 | 37 | 179 | 92,835 |
| NC_003888.3 | 71 | 142 | 288,921 |
| NC_006273 | 28 | 178 | 16,872 |
| NC_001133 | 37 | 182 | 12,955 |
| NC_000866 | 41 | 172 | 9,935 |
| NC_001416 | 42 | 174 | 3,702 |
| NC_001604 | 39 | 177 | 2,454 |
| NC_001895 | 38 | 178 | 2,190 |

Between 142 and 182 rows of text, out of up to 288,921, and between 28 and 71
blocks. The cache is doing what it was built to do.

## The one engine change

`lua/plasmid/view/linear.lua`, in `M.ensure`, plus a counter in `M.plan_of` and
`M.drop_block`: the plan cache is now swept.

A position query -- jumping to a base pair, asking what is under the cursor --
plans the block it lands in and keeps the plan, so that the next query about
the same block is free. The eviction in `ensure` only forgot blocks that had
been *built*, so plans made by position queries alone were never dropped.
Measured on E. coli: a position query in every one of its 51,574 blocks left
51,574 plans and grew the heap by 27.9 MB, from 14.2 MB to 42.1 MB, and
scrolling did not take any of it back. With the sweep, one screen of scrolling
brings it to 26 plans and 12.1 MB.

Plans are cheap to make again (2 microseconds a block, measured as 103 ms for
all 51,574), so the ones outside the live window are dropped when there are
more than four times as many as the window could account for. Nothing else
changed: the planner and the renderer are untouched, and the parity matrix
above was run against the engine with this change in it.

## Budgets in the suite

`tests/spec/budget_spec.lua` keeps this from regressing without the network or
a download: it generates a 400 kb record with 4,000 overlapping features --
joins, complements, enough overlap to fill the feature tracks -- and measures
the same things at a tenth of E. coli's scale. Its limits, and what they were
measured at, are in the comment at the top of that file. The short version:

- the windowed layout of that record takes 400 ms and 15 MB of heap at most
  (measured 47 ms and 2.6 MB), and must lay out exactly zero rows;
- the eager layout of the same record must cost at least 3x the time and 5x the
  heap (the lowest seen over six runs was 9.5x and 17.2x) -- a ratio, so a slow
  CI runner cannot flake it, and the assertion that would fail if the view went
  back to laying everything out up front;
- one screen under 150 ms, 200 pages under 3 s, 400 position queries under
  300 ms and building 200 cold rows under 1.5 s;
- a position query in every block, then one screen of scrolling, must leave
  under an eighth of the plans and under half the heap they cost;
- through a real buffer: opening under 3 s and 60 MB, 200 pages under 3 s, 40
  jumps under 2 s with the cursor reporting the right base every time, two
  resizes under 4 s.

The spec takes about 2.5 s and runs in both sequence layouts.

## What this does not cover

- **The circular view.** Only the linear view is windowed. `:DnaView circular`
  on a megabase record still builds the whole thing; it was not measured here.
- **Editing a large record.** Every measurement is read-only viewing. An edit
  rebuilds the layout, which costs what opening costs -- a full re-plan, 1.1 s
  on E. coli -- and nothing here measures a sequence of edits.
- **Highlighting.** The decoration provider only runs for rows on screen, and
  no measurement here forces a redraw: `:redraw` is not safe in a headless
  Neovim 0.12.2, so the cost of drawing highlights on a genome is untested.
  What is tested is that the highlight *runs* the provider would apply are
  identical between the two layouts.
- **Peak heap during a long jump sequence.** The plan sweep runs on the next
  `ensure`, so a user who jumps many times without scrolling can still reach
  the 42 MB seen above on E. coli before anything is swept.
- **Records above 10 Mb.** The largest measured here is 8.7 Mb. Nothing on the
  scale of a human chromosome was fetched or tried, so where the current
  approach stops being comfortable is not known.
- **Real terminals.** Every measurement is headless. Scrolling is a viewport
  move plus the fill the autocmds would do, which is what the plugin does, but
  not what a terminal's own redraw costs.
