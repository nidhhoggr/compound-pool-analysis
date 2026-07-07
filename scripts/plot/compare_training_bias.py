#!/usr/bin/env python3
"""Test whether ChEMBL augmentation (scripts/tl/03_augment_from_chembl.py)
introduced a structural bias into the mol2mol training corpus, by comparing
descriptor distributions between the original hand-picked seeds and the
ChEMBL-pulled analogs that got merged in alongside them.

No docking scores involved -- this is a pure structure/descriptor comparison,
since training compounds were never docked (see README: ChEMBL hits are
structural neighbors, not validated actives).

Since the shipped compounds_train.smi/compounds_val.smi have no ChEMBL-ID
tag column (just plain SMILES), the seed/analog split is reconstructed by
canonical-SMILES exact match against references.smi (the original seed
pool): anything in train+val matching a reference is 'seed', everything
else is 'chembl_analog'.

Optionally also compares against a third group: your already-classified
generated top-hits (classify_seeds.py output), filtered to the mol2mol
source only (id prefix 'm'), to see whether the generated output's
descriptor profile tracks the analogs' shift (if any) rather than the
original seeds'.

Usage:
  python compare_training_bias.py --references references.smi \\
      --train compounds_train.smi --val compounds_val.smi \\
      --out-prefix results/training_bias

  # with the generated/mol2mol comparison group:
  python compare_training_bias.py --references references.smi \\
      --train compounds_train.smi --val compounds_val.smi \\
      --generated tc_merged_classified.csv \\
      --out-prefix results/training_bias
"""
import sys, argparse
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ks_2samp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

CF3 = Chem.MolFromSmarts("[CX4](F)(F)F")
HALO = Chem.MolFromSmarts("[F,Cl,Br,I]")

def read_plain_smi(path):
    smis = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smis.append(line.split()[0])  # first whitespace-separated token
    return smis

def canonical(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else None

def descriptors(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return dict(
        mw=round(Descriptors.MolWt(m), 1),
        clogp=round(Crippen.MolLogP(m), 2),
        n_rings=rdMolDescriptors.CalcNumRings(m),
        n_arom_rings=rdMolDescriptors.CalcNumAromaticRings(m),
        has_CF3=int(m.HasSubstructMatch(CF3)),
        n_halogen=len(m.GetSubstructMatches(HALO)),
        n_heavy_atoms=m.GetNumHeavyAtoms(),
    )

DESC_COLS = ["mw", "clogp", "n_rings", "n_arom_rings", "n_halogen", "n_heavy_atoms"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--references", required=True, help="original seed pool, plain SMILES")
    ap.add_argument("--train", required=True, help="compounds_train.smi")
    ap.add_argument("--val", default=None, help="compounds_val.smi (optional, pooled with train)")
    ap.add_argument("--generated", default=None,
                     help="classify_seeds.py output CSV of generated hits (optional third group)")
    ap.add_argument("--generated-source-prefix", default="m",
                     help="only keep generated rows whose cmpd id starts with this letter (default 'm' = mol2mol)")
    ap.add_argument("--out-prefix", default="training_bias")
    args = ap.parse_args()

    ref_smis = read_plain_smi(args.references)
    ref_canon = set(filter(None, (canonical(s) for s in ref_smis)))
    sys.stderr.write(f"references: {len(ref_smis)} seeds, {len(ref_canon)} valid/canonicalized\n")

    corpus_smis = read_plain_smi(args.train)
    if args.val:
        corpus_smis += read_plain_smi(args.val)
    sys.stderr.write(f"training corpus (train+val): {len(corpus_smis)} compounds\n")

    rows = []
    n_bad = 0
    for smi in corpus_smis:
        c = canonical(smi)
        d = descriptors(smi)
        if c is None or d is None:
            n_bad += 1
            continue
        d["smiles"] = c
        d["group"] = "seed" if c in ref_canon else "chembl_analog"
        rows.append(d)
    if n_bad:
        sys.stderr.write(f"warning: {n_bad} unparseable SMILES skipped\n")

    df = pd.DataFrame(rows)
    n_seed = (df["group"] == "seed").sum()
    n_analog = (df["group"] == "chembl_analog").sum()
    sys.stderr.write(f"split: {n_seed} matched original seeds, {n_analog} chembl_analog "
                      f"(expected seed count <= {len(ref_canon)}; if much lower, some seeds "
                      f"may have failed to survive standardization into the training file)\n\n")

    groups = ["seed", "chembl_analog"]

    if args.generated:
        gdf = pd.read_csv(args.generated)
        for needed in ("cmpd", "smiles"):
            if needed not in gdf.columns:
                sys.exit(f"error: {args.generated} has no {needed!r} column")
        gdf = gdf[gdf["cmpd"].astype(str).str[:1].str.lower() == args.generated_source_prefix.lower()].copy()
        sys.stderr.write(f"generated (source prefix '{args.generated_source_prefix}'): {len(gdf)} compounds\n")
        # classify_seeds.py output has no heavy-atom count -- compute it here for a fair comparison
        gdf["n_heavy_atoms"] = gdf["smiles"].apply(
            lambda s: (Chem.MolFromSmiles(s).GetNumHeavyAtoms()
                       if Chem.MolFromSmiles(s) is not None else np.nan))
        gdf["group"] = "generated_mol2mol"
        keep = [c for c in DESC_COLS if c in gdf.columns] + ["smiles", "group"]
        df = pd.concat([df, gdf[keep]], ignore_index=True)
        groups.append("generated_mol2mol")

    # --- statistical comparisons ---
    sys.stderr.write("\nPairwise comparisons (Mann-Whitney U, robust to outliers; KS as a second check):\n")
    stat_rows = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1, g2 = groups[i], groups[j]
            for col in DESC_COLS:
                if col not in df.columns:
                    continue
                a = df.loc[df.group == g1, col].dropna()
                b = df.loc[df.group == g2, col].dropna()
                if len(a) < 3 or len(b) < 3:
                    continue
                u, p_mw = mannwhitneyu(a, b, alternative="two-sided")
                d, p_ks = ks_2samp(a, b)
                stat_rows.append(dict(group_a=g1, group_b=g2, descriptor=col,
                                       n_a=len(a), n_b=len(b),
                                       median_a=round(a.median(), 2), median_b=round(b.median(), 2),
                                       mannwhitney_p=round(p_mw, 4), ks_p=round(p_ks, 4)))
    stat_df = pd.DataFrame(stat_rows)
    sys.stderr.write(stat_df.to_string(index=False) + "\n")
    stat_df.to_csv(f"{args.out_prefix}_stats.csv", index=False)

    # --- plots: one boxplot per descriptor, groups side by side ---
    present_cols = [c for c in DESC_COLS if c in df.columns]
    fig, axes = plt.subplots(1, len(present_cols), figsize=(3.2 * len(present_cols), 4))
    if len(present_cols) == 1:
        axes = [axes]
    rng = np.random.default_rng(0)
    for ax, col in zip(axes, present_cols):
        data = [df.loc[df.group == g, col].dropna().values for g in groups]
        ax.boxplot(data, showmeans=True, showfliers=False)  # fliers replaced by the strip below
        for i, vals in enumerate(data, start=1):
            jitter = rng.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, s=10, alpha=0.5,
                       color="steelblue", zorder=3, edgecolors="none")
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups, rotation=30, fontsize=7, ha="right")
        ax.set_title(col, fontsize=9)
    fig.suptitle("Descriptor distributions: original seeds vs ChEMBL analogs"
                 + (" vs generated" if args.generated else ""))
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_boxplots.png", dpi=200)
    plt.close(fig)

    df.to_csv(f"{args.out_prefix}_labeled.csv", index=False)
    sys.stderr.write(f"\nwrote {args.out_prefix}_stats.csv, {args.out_prefix}_boxplots.png, "
                      f"{args.out_prefix}_labeled.csv\n")

if __name__ == "__main__":
    main()
