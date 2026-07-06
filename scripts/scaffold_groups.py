#!/usr/bin/env python3
"""Group compounds by EXACT Murcko scaffold identity (not fuzzy Butina
clustering) and produce:
  1. A frequency table of scaffold groups, sorted descending, to stderr/stdout
  2. An output CSV = input CSV + a `scaffold_group` column, where the top-N
     most frequent scaffolds get group ids 1..N (largest = 1) and everything
     else is labeled "other".

This sidesteps Butina similarity-threshold tuning entirely: two compounds are
in the same group iff their canonical Murcko SMILES string is identical.
Threshold-free, deterministic, and won't split/merge scaffolds based on
peripheral substituent noise.

Input: the CSV produced by classify_seeds.py (must have `murcko` and `cmpd`
columns; any additional columns, e.g. score/mw/clogp, are passed through).

Usage:
  python scaffold_groups.py classified_0.40.csv -o scaffold_labeled.csv --top 10
  python scaffold_groups.py classified_0.40.csv --top 10   # table only, no -o
"""
import sys, csv, argparse
from collections import Counter, defaultdict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("-o", "--out", default=None,
                     help="write input CSV + scaffold_group column here")
    ap.add_argument("--top", type=int, default=10,
                     help="number of top scaffolds to keep as named groups; rest -> 'other'")
    ap.add_argument("--scaffold-col", default="murcko",
                     help="column holding the canonical scaffold SMILES")
    ap.add_argument("--min-size", type=int, default=1,
                     help="drop the frequency-table print for groups smaller than this "
                          "(does not affect scaffold_group labeling)")
    args = ap.parse_args()

    with open(args.infile, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        if args.scaffold_col not in fieldnames:
            sys.exit(f"error: column {args.scaffold_col!r} not found in header: {fieldnames}")
        rows = list(reader)

    # Count exact scaffold occurrences (skip blank/unparseable scaffolds)
    counts = Counter()
    members = defaultdict(list)
    n_blank = 0
    for r in rows:
        scaf = (r.get(args.scaffold_col) or "").strip()
        if not scaf:
            n_blank += 1
            continue
        counts[scaf] += 1
        members[scaf].append(r.get("cmpd", "?"))

    ranked = counts.most_common()
    top_scaffolds = [s for s, _ in ranked[:args.top]]
    group_id_of = {s: i + 1 for i, s in enumerate(top_scaffolds)}

    # --- frequency table ---
    total = sum(counts.values())
    sys.stderr.write(f"\n{len(rows)} compounds total"
                      + (f" ({n_blank} with blank/unparseable scaffold)" if n_blank else "")
                      + f"\n{len(counts)} distinct exact Murcko scaffolds\n\n")
    sys.stderr.write(f"{'grp':>4}  {'n':>4}  {'%':>5}  scaffold\n")
    covered = 0
    for i, (scaf, n) in enumerate(ranked, 1):
        if n < args.min_size:
            break
        covered += n
        marker = f"{i:>4}" if i <= args.top else "   -"
        sys.stderr.write(f"{marker}  {n:>4}  {100*n/total:5.1f}  {scaf}\n")
        if i == args.top and len(ranked) > args.top:
            n_other = total - covered
            sys.stderr.write(f"{'other':>4}  {n_other:>4}  {100*n_other/total:5.1f}  "
                              f"[{len(ranked)-args.top} distinct scaffolds, singletons/small groups]\n")
            break
    sys.stderr.write(f"\ntop-{args.top} groups cover {covered}/{total} compounds "
                      f"({100*covered/total:.1f}%)\n")

    # --- labeled output CSV ---
    if args.out:
        out_fields = fieldnames + ["scaffold_group"]
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=out_fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                scaf = (r.get(args.scaffold_col) or "").strip()
                if scaf in group_id_of:
                    r["scaffold_group"] = f"scaffold_{group_id_of[scaf]}"
                elif scaf:
                    r["scaffold_group"] = "other"
                else:
                    r["scaffold_group"] = "unparseable"
                w.writerow(r)
        sys.stderr.write(f"\nwrote {args.out}\n")

if __name__ == "__main__":
    main()
