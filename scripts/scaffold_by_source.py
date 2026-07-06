#!/usr/bin/env python3
"""Test whether scaffold diversity differs by generation source (job_id),
before deciding whether to split downstream analysis by pool.

For each distinct value of --source-col (default: job_id), reports:
  n compounds, n distinct Murcko scaffolds, singleton rate (%),
  and the single most common scaffold's share of that source's compounds.

A source pool with a high singleton rate (most scaffolds appear once) is
behaving as a diversity generator. A source pool with a low singleton rate
(scaffolds repeat a lot) is behaving as a conservative/local elaborator.
If these numbers look similar across sources, splitting the analysis by
source won't buy you much. If they're very different, that's your answer.

Input: a CSV with an id column and a smiles column (defaults ID / SMILES --
override with --id-col/--smiles-col to match your header).

Source is determined by the first letter of the id column:
  starts with 'i' -> inpaint (DiffSBDD)
  starts with 'm' -> mol2mol (REINVENT4)
  anything else   -> other
Override with --source-col if you'd rather group by an explicit column
(e.g. job_id) instead of the id-prefix convention.

Usage:
  python scaffold_by_source.py tc_merged.csv
  python scaffold_by_source.py tc_merged.csv --smiles-col SMILES
  python scaffold_by_source.py tc_merged.csv --source-col job_id   # explicit column instead of id prefix
"""
import sys, csv, argparse
from collections import Counter, defaultdict
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

PREFIX_MAP = {"i": "inpaint", "m": "mol2mol"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--id-col", default="ID")
    ap.add_argument("--smiles-col", default="SMILES")
    ap.add_argument("--source-col", default=None,
                     help="if set, group by this column's exact value instead of id-prefix")
    args = ap.parse_args()

    with open(args.infile, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []

        def find_col(name):
            for c in fieldnames:
                if c.strip().lower() == name.strip().lower():
                    return c
            sys.exit(f"error: column {name!r} not found in header: {fieldnames}")

        id_key = find_col(args.id_col)
        smi_key = find_col(args.smiles_col)
        src_key = find_col(args.source_col) if args.source_col else None

        rows = list(reader)

    by_source_scaffolds = defaultdict(list)
    n_unparseable = 0
    n_unrecognized_prefix = 0
    for r in rows:
        smi = (r.get(smi_key) or "").strip()
        if src_key:
            src = (r.get(src_key) or "").strip() or "(blank)"
        else:
            mol_id = (r.get(id_key) or "").strip()
            first = mol_id[:1].lower()
            src = PREFIX_MAP.get(first)
            if src is None:
                n_unrecognized_prefix += 1
                src = f"other(prefix={first or 'blank'})"
        if not smi:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None:
            n_unparseable += 1
            continue
        murcko = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
        by_source_scaffolds[src].append(murcko)

    if n_unparseable:
        sys.stderr.write(f"warning: {n_unparseable} unparseable SMILES skipped\n")
    if n_unrecognized_prefix:
        sys.stderr.write(f"warning: {n_unrecognized_prefix} ids didn't start with 'i' or 'm' "
                          f"-- check the other(prefix=...) rows below\n")
    if n_unparseable or n_unrecognized_prefix:
        sys.stderr.write("\n")

    sys.stderr.write(f"{'source':<40}{'n':>6}{'distinct':>10}{'singleton%':>12}{'top_scaf%':>12}\n")
    total_n, total_distinct = 0, 0
    for src, scaffolds in sorted(by_source_scaffolds.items(), key=lambda kv: -len(kv[1])):
        n = len(scaffolds)
        counts = Counter(scaffolds)
        n_distinct = len(counts)
        n_singleton = sum(1 for c in counts.values() if c == 1)
        singleton_pct = 100 * n_singleton / n_distinct if n_distinct else 0
        top_scaf_pct = 100 * counts.most_common(1)[0][1] / n if n else 0
        sys.stderr.write(f"{src:<40}{n:>6}{n_distinct:>10}{singleton_pct:>11.1f}%{top_scaf_pct:>11.1f}%\n")
        total_n += n
        total_distinct += n_distinct

    sys.stderr.write(f"\n{'TOTAL':<40}{total_n:>6}{total_distinct:>10}\n")
    sys.stderr.write("\nsingleton% = fraction of that source's DISTINCT scaffolds that appear exactly once\n"
                      "top_scaf%  = share of that source's COMPOUNDS covered by its single most common scaffold\n"
                      "Low singleton% / high top_scaf% => conservative/local elaboration.\n"
                      "High singleton% / low top_scaf% => diversity generation.\n")

if __name__ == "__main__":
    main()
