#!/bin/sh
# Fetch the large-record corpus the windowed linear view is measured against.
#
#   tests/tools/fetch_big.sh [DIR]      default DIR is ./big
#
# These records are far too large to commit (the E. coli GenBank file alone is
# 11 MB), so the tree keeps this script instead and the corpus is rebuilt on
# demand. tests/spec/budget_spec.lua deliberately does not need any of them: it
# generates its own large record, so the suite still runs offline.
#
# Downloads go through NCBI E-utilities without an API key, so requests are
# spaced out. The five records came to 25.9 MB when this was last run; each
# line below prints what it got. $NCBI_API_KEY is used when it is set, and
# never echoed.
set -eu

dir="${1:-big}"
mkdir -p "$dir"
base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
key=""
[ -n "${NCBI_API_KEY:-}" ] && key="&api_key=$NCBI_API_KEY"

# accession(s) | rettype | file | what it is here for
fetch() {
  acc="$1"; typ="$2"; out="$dir/$3"; why="$4"
  if [ -s "$out" ]; then
    printf '%-22s %10s bytes  (have it)  %s\n' "$3" "$(wc -c <"$out")" "$why"
    return 0
  fi
  curl -sS --fail --retry 3 --retry-delay 2 \
    "$base?db=nuccore&id=$acc&rettype=$typ&retmode=text$key" -o "$out.part"
  mv "$out.part" "$out"
  printf '%-22s %10s bytes  %s\n' "$3" "$(wc -c <"$out")" "$why"
  sleep 1
}

fetch NC_000913.3 gbwithparts NC_000913.gbk "4.6 Mb bacterial genome, feature dense"
fetch NC_006273.2 gbwithparts NC_006273.gbk "235 kb virus, features on top of each other"
fetch NC_001136.10 gbwithparts NC_001136.gbk "1.5 Mb eukaryotic chromosome, spliced CDS"
fetch NC_001416.1,NC_001604.1,NC_000866.4,NC_001895.1 gbwithparts phages.gbk "four phage genomes in one file"
fetch NC_003888.3 fasta NC_003888.fasta "8.7 Mb of residues and no features at all"
