#!/usr/bin/env python3
"""Merge two or more classify_seeds.py output CSVs into one pooled file with
a `source` column, for correlation analysis that needs wider docking-score
variance than a single filtered pool provides (e.g. seed baseline vs.
top-scoring generated hits).

Each input CSV must already have the classify_seeds.py column set
(cmpd, score, cluster, murcko, generic_scaffold, mw, clogp, n_rings,
n_arom_rings, has_CF3, n_halogen, has_amide, biaryl, smiles).

Usage:
  python merge_with_source.py --label seed seeds_classified.csv \\
                               --label generated tc_merged_classified.csv \\
                               -o pooled.csv

Then feed pooled.csv straight into plot_correlations.py with --group-col
source -- this reuses its existing leave-one-group-out sensitivity check to
test whether any global correlation survives dropping either population
(the standard guard against a pooled-population batch-effect confound).
"""
import sys, csv, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", action="append", nargs=2, metavar=("NAME", "PATH"),
                     required=True,
                     help="repeatable: --label seed seeds.csv --label generated hits.csv")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    all_rows = []
    fieldnames = None
    for name, path in args.label:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                missing = set(fieldnames) - set(reader.fieldnames or [])
                extra = set(reader.fieldnames or []) - set(fieldnames)
                sys.stderr.write(
                    f"warning: {path} header differs from first file "
                    f"(missing: {missing or 'none'}, extra: {extra or 'none'}) "
                    f"-- proceeding, mismatched columns will be blank\n")
            n = 0
            for r in reader:
                r["source"] = name
                all_rows.append(r)
                n += 1
            sys.stderr.write(f"{name}: {n} compounds from {path}\n")

    out_fields = fieldnames + ["source"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    sys.stderr.write(f"\nwrote {len(all_rows)} total compounds to {args.out}\n")

if __name__ == "__main__":
    main()
