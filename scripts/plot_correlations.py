#!/usr/bin/env python3
"""Correlation analysis of docking score vs. molecular descriptors, with
scaffold-family highlighting for the Results-section figure.

Input: the CSV produced by scaffold_groups.py (classify_seeds.py output +
`scaffold_group` column). Requires: score, mw, clogp, n_rings, n_arom_rings,
has_CF3, n_halogen, has_amide, biaryl, smiles, scaffold_group.

Computes:
  - Ligand efficiency (score / n_heavy_atoms) via RDKit on the `smiles` column,
    since raw docking score is confounded with molecular size (more atoms ->
    more possible favorable contacts, independent of true per-atom binding
    quality). LE is included alongside raw score in all outputs so it's easy
    to check whether a correlation survives size-normalization.
  - Spearman (not Pearson) rank correlation between score/LE and each
    descriptor, since docking-score and descriptor distributions are rarely
    normal and are easily dominated by outliers.
  - A global correlation matrix, PLUS the same matrix recomputed with each
    of the top named scaffold groups excluded one at a time, so you can see
    whether any single scaffold family is driving a correlation
    (Simpson's-paradox-style check).

Outputs (written next to --out-prefix):
  {prefix}_correlations.csv   - Spearman rho + p-value, global and per-exclusion
  {prefix}_heatmap.png        - Spearman correlation heatmap (score & LE rows)
  {prefix}_scatter_<desc>.png - scatter of score vs. each descriptor, colored
                                 by scaffold_group (named groups colored,
                                 'other' in muted grey, drawn first/behind)

Usage:
  python plot_correlations.py scaffold_labeled.csv --out-prefix results/fig
"""
import sys, argparse
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, linregress
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from rdkit import Chem
except ImportError:
    Chem = None

DESCRIPTORS = ["mw", "clogp", "n_rings", "n_arom_rings",
               "has_CF3", "n_halogen", "has_amide", "biaryl"]

def heavy_atom_count(smi):
    if Chem is None:
        return np.nan
    m = Chem.MolFromSmiles(smi)
    return m.GetNumHeavyAtoms() if m else np.nan

def load(path):
    df = pd.read_csv(path)
    for col in ["score"] + DESCRIPTORS:
        if col not in df.columns:
            sys.exit(f"error: required column {col!r} not found. Columns present: {list(df.columns)}")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    for col in DESCRIPTORS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "smiles" in df.columns and Chem is not None:
        df["n_heavy_atoms"] = df["smiles"].apply(heavy_atom_count)
        # Standard (Hopkins) ligand efficiency convention: LE = -score / N_heavy,
        # a POSITIVE number where higher = more efficient per heavy atom.
        # (docking score itself is negative for favorable binding, so this
        # flips sign relative to raw score/N_heavy.)
        df["ligand_efficiency"] = -df["score"] / df["n_heavy_atoms"]
    else:
        df["ligand_efficiency"] = np.nan
    return df.dropna(subset=["score"])

def spearman_table(df, target_cols, desc_cols, label):
    rows = []
    for tgt in target_cols:
        for desc in desc_cols:
            sub = df[[tgt, desc]].dropna()
            if len(sub) < 5:
                rows.append(dict(subset=label, target=tgt, descriptor=desc,
                                  n=len(sub), rho=np.nan, p=np.nan))
                continue
            rho, p = spearmanr(sub[tgt], sub[desc])
            rows.append(dict(subset=label, target=tgt, descriptor=desc,
                              n=len(sub), rho=round(rho, 3), p=round(p, 4)))
    return rows

BINARY_DESCRIPTORS = {"has_CF3", "has_amide", "biaryl"}  # true 0/1 flags.
# n_halogen is a COUNT (0,1,2,3+ halogens), not binary -- its slope is
# kcal/mol per +1 halogen atom, same interpretation as mw/clogp/n_rings.

def regression_table(df, target_cols, desc_cols, label):
    """OLS slope (target units per 1 unit of descriptor) + intercept + Pearson r/p.
    For a strictly binary 0/1 descriptor (has_CF3, has_amide, biaryl), the slope
    is exactly the difference in group means (present - absent), which is the
    natural effect-size unit for a binary predictor. For continuous/count
    descriptors (mw, clogp, n_rings, n_arom_rings, n_halogen), the slope is
    "change in target per +1 unit of descriptor" (e.g. kcal/mol per Da for mw).
    Note: OLS is sensitive to outliers -- docking-score distributions often
    have a heavy tail, so treat the slope as indicative, not exact, and check
    the scatter plot before quoting it as a precise effect size.
    """
    rows = []
    for tgt in target_cols:
        for desc in desc_cols:
            sub = df[[tgt, desc]].dropna()
            if len(sub) < 5 or sub[desc].nunique() < 2:
                rows.append(dict(subset=label, target=tgt, descriptor=desc, n=len(sub),
                                  slope=np.nan, intercept=np.nan, r=np.nan, p=np.nan,
                                  kind="binary" if desc in BINARY_DESCRIPTORS else "continuous"))
                continue
            res = linregress(sub[desc], sub[tgt])
            rows.append(dict(subset=label, target=tgt, descriptor=desc, n=len(sub),
                              slope=round(res.slope, 4), intercept=round(res.intercept, 3),
                              r=round(res.rvalue, 3), p=round(res.pvalue, 4),
                              kind="binary" if desc in BINARY_DESCRIPTORS else "continuous"))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--out-prefix", default="corr_out")
    ap.add_argument("--group-col", default="scaffold_group")
    args = ap.parse_args()

    df = load(args.infile)
    targets = ["score", "ligand_efficiency"]
    targets = [t for t in targets if df[t].notna().sum() >= 5]

    all_rows = spearman_table(df, targets, DESCRIPTORS, "ALL")
    reg_rows = regression_table(df, targets, DESCRIPTORS, "ALL")

    # leave-one-scaffold-out check on the named groups (Simpson's-paradox guard)
    if args.group_col in df.columns:
        named_groups = sorted(g for g in df[args.group_col].unique()
                               if g not in ("other", "unparseable"))
        for g in named_groups:
            sub = df[df[args.group_col] != g]
            all_rows += spearman_table(sub, targets, DESCRIPTORS, f"EXCL_{g}")

    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(f"{args.out_prefix}_correlations.csv", index=False)

    reg_df = pd.DataFrame(reg_rows)
    reg_df.to_csv(f"{args.out_prefix}_regression.csv", index=False)

    # flag descriptors where sign/magnitude of rho shifts a lot when a group is dropped
    global_rho = out_df[out_df.subset == "ALL"].set_index(["target", "descriptor"])["rho"]
    sys.stderr.write("\nGlobal Spearman correlations (score / ligand_efficiency vs descriptors):\n")
    sys.stderr.write(out_df[out_df.subset == "ALL"].to_string(index=False) + "\n")

    global_reg = reg_df[reg_df.subset == "ALL"]
    sys.stderr.write("\nOLS regression (real units -- kcal/mol per unit descriptor for continuous, "
                      "kcal/mol group difference for binary):\n")
    sys.stderr.write(global_reg.to_string(index=False) + "\n")
    sys.stderr.write("\nFor 'continuous' rows: slope = change in target per +1 unit of descriptor "
                      "(e.g. kcal/mol per Da of MW).\n"
                      "For 'binary' rows: slope = mean(target | descriptor=1) - mean(target | descriptor=0).\n"
                      "OLS is sensitive to outliers -- cross-check against the scatter plots before "
                      "quoting a slope as a precise effect size.\n")

    if args.group_col in df.columns and named_groups:
        sys.stderr.write("\nSensitivity check (|change in rho| > 0.15 when scaffold group excluded):\n")
        flagged = False
        for g in named_groups:
            sub = out_df[out_df.subset == f"EXCL_{g}"].set_index(["target", "descriptor"])["rho"]
            for key in global_rho.index:
                if key in sub.index and pd.notna(global_rho[key]) and pd.notna(sub[key]):
                    diff = sub[key] - global_rho[key]
                    if abs(diff) > 0.15:
                        flagged = True
                        sys.stderr.write(f"  excluding {g}: {key[0]} vs {key[1]} "
                                          f"rho {global_rho[key]:+.3f} -> {sub[key]:+.3f}\n")
        if not flagged:
            sys.stderr.write("  none -- correlations look stable across scaffold subsets.\n")

    # --- heatmap ---
    pivot = out_df[out_df.subset == "ALL"].pivot(index="target", columns="descriptor", values="rho")
    pivot = pivot[DESCRIPTORS]
    fig, ax = plt.subplots(figsize=(1.1 * len(DESCRIPTORS) + 2, 2.2))
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(pivot.columns))); ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index))); ax.set_yticklabels(pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if abs(v) > 0.5 else "black", fontsize=9)
    fig.colorbar(im, ax=ax, label="Spearman rho", shrink=0.8)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_heatmap.png", dpi=200)
    plt.close(fig)

    # --- scatter plots, scaffold-colored ---
    if args.group_col in df.columns:
        groups = sorted(g for g in df[args.group_col].unique() if g != "other")
        cmap = plt.get_cmap("tab10")
        color_of = {g: cmap(i % 10) for i, g in enumerate(groups)}
        color_of["other"] = (0.75, 0.75, 0.75, 0.6)

        for desc in DESCRIPTORS:
            if df[desc].nunique() <= 1:
                continue
            fig, ax = plt.subplots(figsize=(5.5, 4.5))
            # draw 'other' first so named groups sit on top
            for g in ["other"] + groups:
                sub = df[df[args.group_col] == g]
                if sub.empty:
                    continue
                ax.scatter(sub[desc], sub["score"], s=22 if g != "other" else 14,
                           color=color_of[g], label=g, alpha=0.9 if g != "other" else 0.5,
                           edgecolors="none")
            # overlay OLS fit line (score vs this descriptor) if available
            fit_row = global_reg[(global_reg.target == "score") & (global_reg.descriptor == desc)]
            if not fit_row.empty and pd.notna(fit_row.iloc[0]["slope"]):
                slope, intercept = fit_row.iloc[0]["slope"], fit_row.iloc[0]["intercept"]
                xs = np.linspace(df[desc].min(), df[desc].max(), 50)
                ax.plot(xs, slope * xs + intercept, "k--", linewidth=1, alpha=0.7,
                        label=f"OLS slope={slope:+.3f}")
            ax.set_xlabel(desc)
            ax.set_ylabel("docking score (kcal/mol)")
            ax.legend(fontsize=7, loc="best", framealpha=0.9)
            fig.tight_layout()
            fig.savefig(f"{args.out_prefix}_scatter_{desc}.png", dpi=200)
            plt.close(fig)

    sys.stderr.write(f"\nwrote {args.out_prefix}_correlations.csv, {args.out_prefix}_regression.csv, "
                      f"{args.out_prefix}_heatmap.png, and per-descriptor scatter PNGs\n")

if __name__ == "__main__":
    main()
