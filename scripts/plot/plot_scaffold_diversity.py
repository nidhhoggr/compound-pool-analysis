#!/usr/bin/env python3
"""Compare scaffold diversity between two (or more) generation sources --
e.g. DiffSBDD inpaint vs REINVENT4 mol2mol -- using exact Murcko scaffold
identity (see scaffold_groups.py / scaffold_by_source.py for why exact
identity, not Butina clustering, is the right tool here).

Produces:
  1. A rank-abundance plot: scaffold rank (x, log) vs. compound count for
     that scaffold (y, log), one curve per source.
  2. A summary bar chart: singleton-scaffold rate and top-scaffold coverage
     per source.
  3. A rarefaction check: distinct_per_compound (distinct scaffolds /
     compound count) is biased upward for smaller pools purely from sample
     size -- with fewer draws, you have fewer chances to redraw the same
     scaffold, independent of any real diversity difference. To control for
     this, every source larger than the smallest is repeatedly subsampled
     (without replacement) down to the smallest source's size, and the
     resulting distinct_per_compound distribution is plotted alongside the
     smallest source's actual (unrarefied) value. If the smallest source's
     value falls comfortably inside the rarefied distribution of a larger
     source, the apparent diversity gap is sample-size noise, not a real
     difference. If it sits clearly outside, that's a real effect.

Input: a classify_seeds.py output CSV (has `cmpd` and `murcko` columns).
Source is derived from the first letter of `cmpd` by default (i=inpaint,
m=mol2mol) -- override with --source-col to use an explicit column instead.

Usage:
  python plot_scaffold_diversity.py tc_merged_classified.csv --out-prefix results/diversity
  python plot_scaffold_diversity.py tc_merged_classified.csv --source-col source --out-prefix results/diversity
"""
import sys, argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PREFIX_MAP = {"i": "inpaint", "m": "mol2mol"}

def derive_source(cmpd):
    first = str(cmpd)[:1].lower()
    return PREFIX_MAP.get(first, f"other(prefix={first or 'blank'})")

def rarefy(scaffolds, target_n, n_iter, rng):
    """Repeatedly subsample `scaffolds` (a list/array of scaffold labels)
    down to target_n without replacement; return array of distinct/target_n
    ratios, one per iteration."""
    scaffolds = np.asarray(scaffolds)
    ratios = np.empty(n_iter)
    for i in range(n_iter):
        sample = rng.choice(scaffolds, size=target_n, replace=False)
        ratios[i] = len(set(sample)) / target_n
    return ratios

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--id-col", default="cmpd")
    ap.add_argument("--scaffold-col", default="murcko")
    ap.add_argument("--source-col", default=None,
                     help="if set, group by this column instead of id-prefix")
    ap.add_argument("--out-prefix", default="scaffold_diversity")
    ap.add_argument("--rarefy-iters", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    df = pd.read_csv(args.infile)
    for col in [args.id_col, args.scaffold_col]:
        if col not in df.columns:
            sys.exit(f"error: column {col!r} not found. Columns present: {list(df.columns)}")

    if args.source_col:
        if args.source_col not in df.columns:
            sys.exit(f"error: --source-col {args.source_col!r} not found.")
        df["source"] = df[args.source_col]
    else:
        df["source"] = df[args.id_col].apply(derive_source)

    df = df.dropna(subset=[args.scaffold_col])
    df = df[df[args.scaffold_col].astype(str).str.strip() != ""]

    sources = sorted(df["source"].unique())
    cmap = plt.get_cmap("tab10")
    color_of = {s: cmap(i % 10) for i, s in enumerate(sources)}

    # --- 1. rank-abundance curve ---
    fig, ax = plt.subplots(figsize=(6, 4.5))
    summary_rows = []
    scaffolds_by_source = {}
    for s in sources:
        sub = df[df["source"] == s]
        scaffolds_by_source[s] = sub[args.scaffold_col].tolist()
        counts = sub[args.scaffold_col].value_counts().sort_values(ascending=False)
        n = len(sub)
        n_distinct = len(counts)
        n_singleton = int((counts == 1).sum())
        singleton_pct = 100 * n_singleton / n_distinct if n_distinct else 0
        top_pct = 100 * counts.iloc[0] / n if n else 0
        summary_rows.append(dict(source=s, n=n, n_distinct=n_distinct,
                                  singleton_pct=round(singleton_pct, 1),
                                  top_scaffold_pct=round(top_pct, 1),
                                  distinct_per_compound=round(n_distinct / n, 3) if n else np.nan))
        ax.plot(range(1, len(counts) + 1), counts.values, marker="o", markersize=3,
                linewidth=1, color=color_of[s], label=f"{s} (n={n})")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("scaffold rank (log)")
    ax.set_ylabel("compounds sharing that scaffold (log)")
    ax.set_title("Scaffold rank-abundance by generation source")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_rank_abundance.png", dpi=200)
    plt.close(fig)

    # --- 2. summary bar chart ---
    summary = pd.DataFrame(summary_rows).set_index("source")
    sys.stderr.write("\n" + summary.to_string() + "\n")
    sys.stderr.write("\nsingleton_pct: % of that source's distinct scaffolds appearing exactly once\n"
                      "top_scaffold_pct: % of that source's compounds covered by its single most common scaffold\n"
                      "distinct_per_compound: distinct scaffolds / compounds (closer to 1.0 = more diverse)\n")

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    x = np.arange(len(summary))
    axes[0].bar(x, summary["singleton_pct"], color=[color_of[s] for s in summary.index])
    axes[0].set_xticks(x); axes[0].set_xticklabels(summary.index)
    axes[0].set_ylabel("% distinct scaffolds that are singletons")
    axes[0].set_ylim(0, 100)
    for i, v in enumerate(summary["singleton_pct"]):
        axes[0].text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=8)

    axes[1].bar(x, summary["top_scaffold_pct"], color=[color_of[s] for s in summary.index])
    axes[1].set_xticks(x); axes[1].set_xticklabels(summary.index)
    axes[1].set_ylabel("% compounds in top scaffold")
    axes[1].set_ylim(0, max(20, summary["top_scaffold_pct"].max() * 1.3))
    for i, v in enumerate(summary["top_scaffold_pct"]):
        axes[1].text(i, v + 0.3, f"{v:.1f}%", ha="center", fontsize=8)

    fig.suptitle("Scaffold diversity summary by generation source")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_summary_bars.png", dpi=200)
    plt.close(fig)

    summary.to_csv(f"{args.out_prefix}_summary.csv")

    # --- 3. rarefaction check ---
    sizes = {s: len(scaffolds_by_source[s]) for s in sources}
    min_source = min(sizes, key=sizes.get)
    target_n = sizes[min_source]
    rng = np.random.default_rng(args.seed)

    rarefied = {}
    for s in sources:
        if s == min_source:
            continue
        rarefied[s] = rarefy(scaffolds_by_source[s], target_n, args.rarefy_iters, rng)

    actual_small = summary.loc[min_source, "distinct_per_compound"]

    sys.stderr.write(f"\nRarefaction: subsampling every larger source down to "
                      f"{min_source}'s size (n={target_n}), {args.rarefy_iters} iterations\n")
    fig, ax = plt.subplots(figsize=(6, 4))
    positions, labels, data = [], [], []
    for i, s in enumerate(sources):
        if s == min_source:
            continue
        r = rarefied[s]
        data.append(r); positions.append(i); labels.append(s)
        pct_below = 100 * (r < actual_small).mean()
        sys.stderr.write(
            f"  {s}: rarefied distinct_per_compound mean={r.mean():.3f} "
            f"(95% range {np.percentile(r,2.5):.3f}-{np.percentile(r,97.5):.3f}) | "
            f"{min_source}'s actual value ({actual_small:.3f}) exceeds {100-pct_below:.1f}% "
            f"of rarefied {s} samples\n")

    vp = ax.violinplot(data, positions=positions, showmeans=True)
    ax.axhline(actual_small, color="black", linestyle="--", linewidth=1.2,
               label=f"{min_source} actual (n={target_n}): {actual_small:.3f}")
    ax.set_xticks(positions); ax.set_xticklabels(labels)
    ax.set_ylabel("distinct scaffolds / compound (rarefied)")
    ax.set_title(f"Rarefaction to n={target_n} -- is the diversity gap real or sample-size noise?")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_rarefaction.png", dpi=200)
    plt.close(fig)

    sys.stderr.write(f"\nwrote {args.out_prefix}_rank_abundance.png, "
                      f"{args.out_prefix}_summary_bars.png, {args.out_prefix}_summary.csv, "
                      f"{args.out_prefix}_rarefaction.png\n")

if __name__ == "__main__":
    main()
