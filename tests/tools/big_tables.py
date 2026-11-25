#!/usr/bin/env python3
"""Turn tests/tools/big_corpus.sh output into the tables of
docs/large_records_report.md.

    tests/tools/big_tables.py AFTER.jsonl [BEFORE.jsonl]

One line of the JSONL is one record in one mode in one run; repeated runs of
the same (file, record, mode) are reduced to their median, so a table never
quotes a single noisy run. With a second file, the before/after columns are
printed as well, which is how a change to the open or the resize path is
reported.
"""
import json
import statistics
import sys
from collections import defaultdict

ORDER = ["NC_000913", "NC_001136", "NC_003888.3", "NC_006273", "NC_001133",
         "NC_000866", "NC_001416", "NC_001604", "NC_001895"]


def load(path):
    runs = defaultdict(list)
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        runs[(r["name"], r["mode"])].append(r)
    out = {}
    for key, rs in runs.items():
        merged = dict(rs[0])
        for field, value in rs[0].items():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                merged[field] = statistics.median(r[field] for r in rs)
        merged["runs"] = len(rs)
        out[key] = merged
    return out


def rank(name):
    return ORDER.index(name) if name in ORDER else len(ORDER)


def rows(data):
    return sorted(data, key=lambda k: (rank(k[0]), k[1]))


def fmt(x, nd=0):
    return f"{x:,.{nd}f}"


def main():
    after = load(sys.argv[1])
    before = load(sys.argv[2]) if len(sys.argv) > 2 else {}

    print("### corpus\n")
    print("| record | file | bp | features | topology | line width | rows | blocks |")
    print("|---|---|---|---|---|---|---|---|")
    for key in rows(after):
        if key[1] != "window":
            continue
        r = after[key]
        print(f"| {r['name']} | {r['file']} | {fmt(r['len'])} | {fmt(r['features'])} | "
              f"{r['topology']} | {fmt(r['W'])} | {fmt(r['rows'])} | {fmt(r.get('nblocks') or 0)} |")

    print("\n### opening and memory\n")
    print("| record | mode | open ms | heap after open MB | rss after open MB | peak heap MB | peak rss MB |")
    print("|---|---|---|---|---|---|---|")
    for key in rows(after):
        r = after[key]
        print(f"| {r['name']} | {r['mode']} | {fmt(r['open_ms'])} | {fmt(r['heap_after_open_mb'])} | "
              f"{fmt(r['rss_after_open_mb'])} | {fmt(r['peak_heap_mb'])} | {fmt(r['peak_rss_mb'])} |")

    print("\n### latency, median / maximum ms\n")
    print("| record | mode | one page | long scroll, per page | jump to a base pair | "
          "resize that changes the line width | resize that does not |")
    print("|---|---|---|---|---|---|---|")
    for key in rows(after):
        r = after[key]
        print(f"| {r['name']} | {r['mode']} | {r['page_ms']:.1f} / {r['page_max_ms']:.1f} | "
              f"{r['long_page_ms']:.1f} / {r['long_page_max_ms']:.1f} | "
              f"{r['jump_ms']:.1f} / {r['jump_max_ms']:.1f} | {fmt(r['resize_ms'])} | "
              f"{r['resize_same_ms']:.2f} / {r['resize_same_max_ms']:.2f} |")

    print("\n### what is still laid out afterwards\n")
    print("| record | blocks built | rows laid out | of total rows |")
    print("|---|---|---|---|")
    for key in rows(after):
        r = after[key]
        if key[1] != "window":
            continue
        print(f"| {r['name']} | {fmt(r['blocks_built'])} | {fmt(r['rows_live'])} | {fmt(r['rows'])} |")

    if not before:
        return
    print("\n### before and after\n")
    print("| record | mode | open ms | resize ms | resize that changes nothing ms | peak heap MB | peak rss MB |")
    print("|---|---|---|---|---|---|---|")
    for key in rows(after):
        if key not in before:
            continue
        a, b = after[key], before[key]
        def pair(field, nd=0):
            return f"{b[field]:,.{nd}f} -> {a[field]:,.{nd}f}"
        print(f"| {a['name']} | {a['mode']} | {pair('open_ms')} | {pair('resize_ms')} | "
              f"{pair('resize_same_ms', 2)} | {pair('peak_heap_mb')} | {pair('peak_rss_mb')} |")


if __name__ == "__main__":
    main()
