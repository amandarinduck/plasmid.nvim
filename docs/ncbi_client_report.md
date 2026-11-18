# NCBI client: recordings, parser check, NCBI quirks

`lua/plasmid/ncbi.lua` is an asynchronous NCBI E-utilities client (efetch / esummary /
esearch over `curl` via `vim.system`). This report covers the fixtures recorded for
its offline tests, request counts and timings, a parser check on the fetched records,
and the E-utilities behaviour the client has to handle. Recorded on 2026-11-17 from
`eutils.ncbi.nlm.nih.gov`, without an API key or email, with `tool=plasmid.nvim` and
curl's default User-Agent (curl 8.5.0).

## Tests

| command | result |
|---|---|
| `nvim -l tests/run.lua ncbi` (replay, offline) | 30 cases passed, ~9 s (most of that is real rate-limit and backoff waiting) |
| `PLASMID_LIVE=1 nvim -l tests/run.lua ncbi` | 30 passed; 17 fetch/search/gene cases re-run against NCBI and compared with the recordings |
| `nvim -l tests/run.lua` (whole suite) | 264 passed, 0 failed; `PLASMID_TEST_GROUP=0`: 264 passed |

The spec covers:

- `parse_accession`: 44 accepted and 17 rejected inputs.
- Fetched records compared with the committed fixtures and with Biopython.
- The mRNA CDS translation compared with the fetched protein.
- Range fetches on both strands, plus the reverse-complement relation between them.
- A CON record.
- A GI number.
- An unknown accession.
- The size guard, with and without force.
- Clipping of ranges that run past the end of a record.
- The cache: hit, alias, refresh, expiry, list and clear.
- POSTed summaries and their ordering.
- Searches, including one with no hits.
- Gene coordinates checked against NCBI's own gene records, plus the 0-based → 1-based conversion.
- 12 error-message cases.
- curl missing.
- Rate-limit spacing, measured on a stub: 350/350/350 ms without a key, 110/110/110 ms with one.
- Retries: Retry-After, exponential backoff, the old 200-with-rate-limit-JSON response, and transient esearch failures.
- Cancel: in flight, while queued, during a backoff, and between the result and the scheduled callback.
- Callbacks on the main loop: never in a fast event and never synchronous; they can write to a buffer.
- That the API key never reaches fixture keys, files or the history.
- `reset()`.
- The real curl transport against a local `vim.uv` HTTP server: GET, POST body and content type, 429 with Retry-After, a 1 s timeout, connection refused, and a kill on cancel.

## What was recorded

Recorder: `nvim -l tests/tools/ncbi_record.lua [name-filter]`. The fixed case list is in
`tests/ncbi_cases.lua`, which is shared with the spec and with live mode.

Each response is stored under `tests/fixtures/ncbi_http/<key>` as two files:

- `<key>.json`: the metadata, with api_key and email scrubbed from the URL and body.
- `<key>.body`: the raw body.

The key is the first 24 hex digits of sha256 over `METHOD url-without-api_key/email \n body`. `manifest.json` lists every case with the requests it made.

Size: 33 responses, 291 kB on disk in total (limit 3 MB). The recordings do not use the `.gb` extension, so the corpus-wide oracle and round-trip specs do not pick them up.

| case | request | requests | body bytes | ms | result |
|---|---|---|---|---|---|
| L09137 | fetch `L09137` | esummary, efetch | 9,007 | 862 | L09137.2, 2,686 bp |
| J01749 | fetch `J01749` | 2 | 25,496 | 764 | J01749.1, 4,361 bp |
| NC_001422 | fetch `NC_001422` | 2 | 27,023 | 766 | NC_001422.1, 5,386 bp |
| NM_000518.5 | fetch (HBB mRNA) | 2 | 23,407 | 714 | 628 bp |
| thr_plus | fetch `NC_000913.3:190-2800` | 2 | 14,789 | 1,346 | 2,611 bp |
| thr_minus | fetch `NC_000913.3:2800-190` | 2 | 14,855 | 1,059 | 2,611 bp, strand=2 |
| NP_000509.1 | fetch (HBB protein) | 2 | 18,613 | 1,042 | 147 aa |
| P69905 | fetch (UniProt accession) | 2 | 114,734 | 1,101 | P69905.2, 142 aa |
| con_LN625376 | fetch `LN625376.1` (CON) | 2 | 5,146 | 770 | 1,765 bp with sequence |
| unknown | fetch `XX999999` | esummary only | 114 | 245 | error |
| big_refused | fetch `NC_000913.3`, max_length 1e6 | esummary only | 884 | 353 | refused (4,641,652 bp) |
| gi_number | fetch `20141090` | esummary (+ shared efetch) | 9,007 | 811 | L09137.2 |
| summary_post | esummary POST, 3 ids | 1 POST | 1,981 | 277 | ordered |
| search_puc19 | `pUC19[title] AND complete sequence[title]` | esearch, esummary | 1,906 | 756 | count 2 |
| search_protein | `hemoglobin subunit beta[title] AND Homo sapiens[orgn] AND refseq[filter]` | 3 (one 429, retried) | 1,338 | 2,827 | NP_000509.1 |
| search_nohit | `zzqxqzzq[title]` | esearch | 347 | 351 | count 0 + warnings |
| gene_BRCA1 | `BRCA1`, Homo sapiens | esearch, esummary | 4,921 | 689 | 672, NC_000017.11:43,044,295-43,170,327, strand 2, 17q21.31 |
| gene_thrA | `thrA`, E. coli K-12 MG1655 | esearch, esummary | 1,666 | 662 | 945803, NC_000913.3:337-2799, strand 1 |
| raw_con_gb | efetch LN625376.1 `rettype=gb` | 1 | 2,004 | 730 | CONTIG line, no ORIGIN |
| raw_efetch_unknown | efetch XX999999 | 1 | 72 | 426 | efetch error format |
| raw_gene_672 / raw_gene_945803 | efetch db=gene text | 1 + 1 | 938 | 791 | NCBI's own Annotation lines |

The last recording run made 37 requests: 18 esummary, 14 efetch and 5 esearch. Bodies totalled 278 kB and the run took 17.3 s, most of it in the rate-limit queue. The typical request took 220–250 ms for an esummary and 360–1,000 ms for an efetch (the E. coli range extract took the longest).

A re-recording produced byte-identical `.body` files. Only the timing and date fields in the `.json` metadata changed.

The live re-check took 0.6–1.5 s per case (two requests each) and 26 s for all 17 cases. The live flat files were byte-identical to the recordings.

Total NCBI traffic during development was roughly 250 requests over about an hour:

- 2 recording runs plus 1 aborted run.
- 5 live runs.
- About 15 hand-run curl probes, spaced at least 0.5 s apart.

## Parser check on real NCBI records

Every fetched flat file was parsed with `require("plasmid.io").parse(text, path, "genbank")`. Each was compared with Biopython 1.88, using `tests/oracle/dump_ncbi_http.py`, which reuses `dump_biopython.dump`. The comparison covers:

- Sequence length and sha256.
- id, name, topology and description.
- Every feature's type and parts (start, stop, strand, ref).
- All qualifiers.

Each file is also written back byte-for-byte. The test is `fetched_records_match_biopython`.

| record | type | length | features | vs Biopython | round-trip |
|---|---|---|---|---|---|
| L09137.2 | circular DNA | 2,686 | 1 | identical | byte-identical |
| J01749.1 | circular DNA | 4,361 | 50 | identical | byte-identical |
| NC_001422.1 | circular ss-DNA | 5,386 | 32 | identical | byte-identical |
| NM_000518.5 | mRNA | 628 | 48 | identical | byte-identical |
| NC_000913.3 190..2800 | range extract | 2,611 | 5 | identical | byte-identical |
| NC_000913.3 complement(190..2800) | range extract | 2,611 | 5 | identical | byte-identical |
| NP_000509.1 | GenPept | 147 | 45 | identical | byte-identical |
| P69905.2 | GenPept (Swiss-Prot) | 142 | 217 | identical | byte-identical |
| LN625376.1 | CON via gbwithparts | 1,765 | 1 | identical | byte-identical |
| L09137.2 via GI 20141090 | circular DNA | 2,686 | 1 | identical | byte-identical |

That is 405 features compared, with no differences. These are the first GenPept records in the test corpus; the parser handled them unchanged.

**Overlap with the existing fixtures:** fetched L09137, J01749 and NC_001422 are byte-identical to `tests/fixtures/ncbi/*.gb`. The sequences and feature counts match (1, 50 and 32 features). NCBI has not released new versions since the fixtures were committed: L09137.2 dates from 2002, J01749.1 from 2008 and NC_001422.1 from 2023.

**Cross-checks between records:**

- The NM_000518.5 CDS has `/protein_id="NP_000509.1"`, and its `/translation` equals the fetched NP_000509.1 sequence.
- The minus-strand extract is exactly the reverse complement of the plus-strand extract.
- thrA (gene 945803, 337..2799) appears as CDS 148..2610 in the 190..2800 extract, and as complement(2..2464) in the reverse extract.

## NCBI quirks found

1. **Unknown accession in esummary:** HTTP 200 with
   `{"error":"Invalid uid XX999999 at position= 0","result":{"uids":[]}}`. The client
   reports `NCBI does not know nuccore accession 'XX999999' (NCBI: Invalid uid …)`.
   With a mixed list, the error is set but the known ids are still returned. The
   client returns those and passes the unknown ones as a third callback argument
   (`missing`).
2. **Unknown accession in efetch:** HTTP 200 and `text/plain`, with the letters spaced
   out: `Error: F a i l e d  t o  u n d e r s t a n d  i d :  X X 9 9 9 9 9 9`. Words are
   separated by two spaces. The client undoes this to `Failed to understand id:
   XX999999`. fetch never gets this far for unknown ids, because the esummary fails
   first.
3. **Rate limit:** HTTP 429 with a JSON body:
   `{"error":"API rate limit exceeded","api-key":"<client IP address>","count":"4","limit":"3"}`.
   The `api-key` field echoes the caller's IP address when no key is used. The
   recorder redacts it from `manifest.json`, and 429 bodies are never stored as
   fixtures.
   - With starts spaced exactly 334 ms apart, NCBI still sent 429 for 2 of 38
     requests. At 350 ms it sent 429 for 1 of 37 requests in one run and none in
     another. Arrival jitter seems to put four requests inside NCBI's one-second
     window. The client therefore spaces requests 350 ms apart (110 ms with a key)
     and retries.
   - When the test suite reset the client between live calls without keeping the
     rate-limit clock, NCBI sent four 429s in one run. `reset()` now keeps the time of
     the last request.
   - Whether NCBI sent a `Retry-After` header on these 429s was not captured. The
     client honours a numeric `Retry-After`, as the stub and local-server tests show.
4. **Transient esearch failure:** HTTP 200 with
   `{"esearchresult":{"ERROR":"Search Backend failed: \t\t\n"}}`. For the pUC19 title
   query it persisted for several minutes, across 3 retries and 3 hand probes, then
   cleared. The client retries it like a 5xx and reports `NCBI esearch: Search Backend
   failed:`. Live mode logs and skips a case when NCBI is in this state.
5. **gbwithparts vs gb for CON records:** with `rettype=gb`, LN625376.1 (division CON,
   1,765 bp) comes back as a `CONTIG join(CBZL010000100.1:1..1765)` line with no
   ORIGIN. Biopython raises `UndefinedSequenceError` on it. With `rettype=gbwithparts`
   the same record has its sequence. The client always uses gbwithparts for nuccore,
   and rejects a CON reply that has no sequence, as a guard.
6. **Range and strand extracts:** with `seq_start/seq_stop/strand`:
   - LOCUS shows the extract length.
   - VERSION keeps `NC_000913.3`.
   - ACCESSION gains `REGION: 190..2800` or `REGION: complement(190..2800)`.
   - Features are clipped and renumbered to the extract.

   `strand=2` without a range also works and returns the whole record reverse-complemented (`REGION: complement(1..2686)`, seen during exploration).
7. **genomicinfo coordinates:** `chrstart` and `chrstop` are 0-based, and
   `chrstop < chrstart` means the minus strand.
   - BRCA1: `chrstart 43170326`, `chrstop 43044294`. Converted, that is
     NC_000017.11 43,044,295..43,170,327, strand 2. This matches NCBI's own gene
     record: `Annotation: Chromosome 17 NC_000017.11 (43044295..43170327,
     complement)`. The current annotation spans 126 kb. That is longer than the often
     quoted 43,044,295..43,125,483, which comes from an earlier annotation. The client
     reports `genomicinfo`, the current annotation, and ignores `locationhist` (16
     entries for BRCA1).
   - thrA: `336 / 2798`. Converted, that is 337..2799 on the plus strand, which matches
     `Annotation:  NC_000913.3 (337..2799)`. Bacterial genes have no "Chromosome"
     prefix and two spaces after `Annotation:`.
   - The gene esummary `status` field is an empty string, not `"0"`, for these live
     genes. So "current only" is enforced in the search instead, with `AND alive[prop]`.
     Genes without any `genomicinfo` are skipped.
8. **esummary and accessions:** ids may be accessions, but results are keyed by GI uid.
   The client maps them back to the requested order by uid, `accessionversion` or the
   unversioned `caption`. `slen` is a JSON number and uids are strings. P69905
   (Swiss-Prot) is in db protein as `P69905.2`. Its GenPept record is 114 kB, mostly
   from 217 features and a long comment.
9. **esearch with no hits:** HTTP 200, count 0, with
   `errorlist.phrasesnotfound = ["zzqxqzzq[title]"]` and
   `warninglist.outputmessages = ["No items found."]`. Both are passed on as
   `res.warnings`.
10. **Proxies:** behind an HTTPS proxy, `curl -D` writes the proxy's
    `HTTP/1.1 200 Connection established` block before the real headers. The transport
    keeps the last `HTTP/` block, which also covers redirects.
11. **The pUC19 title search returns only 2 hits** (`pUC19[title] AND complete
    sequence[title]`). L09137 ("pUC19c") is not among them, because title search
    matches whole words.

## Design notes

- **Cancel:** `handle:cancel()` kills curl, or removes the request from the queue, or
  stops a backoff timer. The callback is never called afterwards, including when the
  result has been computed but its scheduled callback has not run yet. cancel returns
  true when it prevented a callback.
- **Retries:** 429, 5xx, network errors (curl exit codes other than 28),
  200-with-rate-limit JSON and "Search Backend failed" are retried. The backoff is
  `retry_base_ms * 2^n`, at least the Retry-After value and at most `retry_max_ms`.
  Timeouts are reported without a retry, so the worst case stays at one `timeout_ms`.
- **The API key:** the URL and POST body are passed to curl on stdin as a `-K -` config
  file, so the key is not visible in the process list. Fixture keys, manifests,
  `M.history` and cache file names are computed from scrubbed URLs.
- **Cache layout:**
  - Files are named `<ACC.VER>[_FROM-TO][_rc].gb|.gp`, each with a `.json` sidecar
    holding the summary fields.
  - `aliases.json` maps unversioned or GI input to the versioned file. A second fetch
    of `L09137` or `20141090` therefore needs no request at all.
  - Files older than `cache_days` are re-fetched.
- **POST:** a request is POSTed when it has more than `post_ids` (200) ids, or when its
  URL would be longer than 2,000 characters.
- **Range clipping:** a range that starts past the end of a record is an error. A range
  that ends past the end is clipped, and `res.to` shows the clipped end.
