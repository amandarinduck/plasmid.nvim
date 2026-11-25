#!/bin/sh
# Run tests/tools/big_budget.lua over the whole large-record corpus, one
# process per record so no two layouts are ever live at once.
#
#   tests/tools/big_corpus.sh OUT.jsonl [DIR]     default DIR is ./big
#
# The corpus comes from tests/tools/fetch_big.sh, plus the yeast chromosome
# that is already a committed fixture. Every line of OUT.jsonl is one record in
# one mode in one run.
#
# The windowed mode is run REPS times (default 3) so the report can quote a
# median; the eager mode once, because the eager layout of E. coli K-12 is
# three gigabytes of resident memory and it is there for contrast, not for a
# latency figure.
set -eu

out="${1:?usage: big_corpus.sh OUT.jsonl [DIR]}"
dir="${2:-big}"
reps="${REPS:-3}"
nvim="${NVIM:-/opt/nvim-linux-x86_64/bin/nvim}"
: >"$out"

run() {
  PLASMID_BUDGET="$1 $2 $out" "$nvim" --headless --clean -u NONE \
    -c 'luafile tests/tools/big_budget.lua' </dev/null
}

# record file, extra keys; NC_001133 is a committed fixture, the rest come
# from fetch_big.sh
corpus() {
  mode="$1"; pages="$2"
  run "$dir/NC_000913.gbk"            "mode=$mode width=120 height=45 pages=$pages jumps=20"
  run "$dir/NC_001136.gbk"            "mode=$mode width=120 height=45 pages=$pages jumps=20"
  run "$dir/NC_003888.fasta"          "mode=$mode width=120 height=45 pages=$pages jumps=20"
  run "$dir/NC_006273.gbk"            "mode=$mode width=120 height=45 pages=$pages jumps=20"
  run tests/fixtures/ncbi/NC_001133.gb "mode=$mode width=120 height=45 pages=$pages jumps=20"
  for r in 1 2 3 4; do
    run "$dir/phages.gbk"             "mode=$mode width=120 height=45 pages=100 jumps=20 rec=$r"
  done
}

i=1
while [ "$i" -le "$reps" ]; do
  corpus window 300
  i=$((i + 1))
done
corpus eager 300
