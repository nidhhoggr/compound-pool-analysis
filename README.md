# compound-pool-analysis — Pipeline Guide

This README walks through generating every file currently under `output/`,
starting from the raw inputs in `input/`. Commands are grouped by pipeline
stage; run them in order from the repo root.

```
input/
├── compounds_train.smi     # mol2mol TL training corpus (post-ChEMBL-augmentation)
├── compounds_val.smi       # mol2mol TL validation split
├── receptor.pdbqt          # 7KEW receptor, prepped for AutoDock Vina
├── references.smi          # the 17 original hand-picked seed compounds (plain SMILES)
├── seeds.smi                # the 16 mol2mol seeds, mol_id/score/smiles format
└── tc_merged.csv            # top-hits shortlist (docking score <= -11.5 kcal/mol),
                              # pooled across both DiffSBDD-inpaint and REINVENT4-mol2mol jobs

scripts/
├── classify_seeds.py        # per-compound descriptors + Murcko scaffold + Butina cluster
├── merge_with_source.py     # concatenate 2+ classify_seeds.py outputs with a source tag
├── multithreshold-classifier.sh   # sweeps classify_seeds.py --sim across thresholds
├── scaffold_by_source.py    # exact-scaffold diversity diagnostic, split by ID-prefix source
├── scaffold_groups.py       # groups compounds by exact Murcko scaffold identity
├── plot/
│   ├── plot_correlations.py       # Spearman + OLS regression, score/LE vs descriptors
│   ├── plot_scaffold_diversity.py # rank-abundance + rarefaction, diversity by source
│   └── compare_training_bias.py   # descriptor drift: seeds vs training corpus vs generated
└── redock/                  # AutoDock Vina redocking + PyMOL candidate review (upstream of
                              # everything above -- not covered in this README, see its own docs)
```

---

## Stage 0 — What you need before starting

Everything below assumes `tc_merged.csv` already exists (i.e. redocking + score
filtering already happened via `scripts/redock/`) and that `seeds.smi` is in
the corrected `mol_id  score  smiles` column order (not the raw
`score  mol_id  smiles` order some source files use — check with `head -3` and
reorder with `awk '{print $2, $1, $3}'` if needed before running anything).

---

## Stage 1 — Classify the top-hits pool

Computes per-compound descriptors (MW, clogP, ring counts, CF3/halogen/amide/biaryl
flags), the canonical Murcko scaffold, and a Butina cluster assignment at a given
Tanimoto similarity threshold (`--sim`, default 0.4).

```bash
python scripts/classify_seeds.py input/tc_merged.csv -o output/tc_merged_classified.csv
```
→ `output/tc_merged_classified.csv` (233 compounds, 184 distinct Murcko scaffolds)

Also classify the two seed pools the same way, so their descriptors are computed
identically and can be pooled later:

```bash
python scripts/classify_seeds.py input/seeds.smi -o output/seeds_classified.csv
```
→ `output/seeds_classified.csv` (16 compounds)

If you also want a generic run with no threshold suffix (e.g. as a default/sanity
copy):
```bash
python scripts/classify_seeds.py input/tc_merged.csv -o output/classified.csv
```
→ `output/classified.csv`

### Threshold sweep (`output/classified_0.15.csv` … `output/classified_0.40.csv`)

`multithreshold-classifier.sh` wraps `classify_seeds.py` across multiple `--sim`
values in one call, to inspect how Butina clustering behavior changes with the
threshold — this is how we diagnosed single-linkage chaining (cluster sizes
ballooning at low thresholds) earlier in the analysis; see the diversity
findings in Stage 4b for why exact-scaffold grouping ended up replacing Butina
clustering for anything diversity-related.

```bash
scripts/multithreshold-classifier.sh tc_merged.csv
```
→ `output/classified_0.40.csv` … `output/classified_0.15.csv`
→ `output/classified_sweep_summary.txt` (per-threshold compound/scaffold/cluster-size summary, captured from stderr)

Defaults: reads `input/tc_merged.csv`, writes to `output/`, sweeps
`0.40 0.35 0.30 0.25 0.20 0.15`, expects `ID`/`docking`/`SMILES` column headers.
All overridable via env var if your input differs, e.g.:
```bash
THRESHOLDS="0.5 0.4 0.3" ID_COL=cmpd_id scripts/multithreshold-classifier.sh other_file.csv
```

---

## Stage 2 — Merge in the seed baseline (fixes docking-score range restriction)

The 233-compound top-hits pool only spans a narrow docking-score band
(≈ −11.5 to −12.8 kcal/mol), which makes any score-vs-descriptor correlation
statistically unreliable (range restriction). Pooling in the 16 mediocre-scoring
hand-picked seeds (≈ −5 to −10 kcal/mol) restores enough variance for the
correlation analysis in Stage 4 to be meaningful.

```bash
python scripts/merge_with_source.py --label seed output/seeds_classified.csv \
                                     --label generated output/tc_merged_classified.csv \
                                     -o output/pooled.csv
```
→ `output/pooled.csv` (249 compounds, tagged `source=seed`/`source=generated`)

---

## Stage 3 — Group by exact scaffold identity

Butina similarity clustering turned out to be the wrong tool for grouping by
chemotype (chaining merges unrelated scaffolds at low thresholds; identical
scaffolds can even land in different clusters). `scaffold_groups.py` instead
groups compounds by **exact** canonical Murcko SMILES, which is deterministic
and chemically unambiguous.

```bash
python scripts/scaffold_groups.py output/tc_merged_classified.csv --top 6 \
    -o output/scaffold_labeled.csv 2> output/scaffold_groups.csv
```
→ `output/scaffold_labeled.csv` (original data + `scaffold_group` column:
`scaffold_1`…`scaffold_6` for the six scaffolds appearing ≥3 times, `other`
for the remaining 187 singleton/near-singleton compounds)
→ `output/scaffold_groups.csv` (the frequency table itself, captured from
stderr — group rank, count, %, and scaffold SMILES)

`--top 6` was chosen because that's where the natural break sits in the
frequency distribution (group 6 has n=3; group 7 would only have n=2, no
more meaningful than the singleton tail) — re-check this with `--top 10`
first on new data before assuming 6 is still the right cutoff.

---

## Stage 4 — Plots

All plotting scripts live in `scripts/plot/`.

### 4a. Correlation + regression (`plot_correlations.py`)

Computes Spearman rho, OLS regression (real kcal/mol units), a correlation
heatmap, per-descriptor scatter plots, and a leave-one-group-out sensitivity
check (tests whether one subgroup is driving a correlation). Ligand efficiency
uses the standard convention (`LE = -score / n_heavy_atoms`, positive = more
efficient).

**Pooled seed-vs-generated version** (`output/plots/pooled/`):
```bash
python scripts/plot/plot_correlations.py output/pooled.csv --group-col source \
    --out-prefix output/plots/pooled/pooled
```
→ `pooled_correlations.csv`, `pooled_regression.csv`, `pooled_heatmap.png`,
`pooled_scatter_{mw,clogp,n_rings,n_arom_rings,has_CF3,n_halogen,has_amide,biaryl}.png`

Sensitivity check here excludes `seed` or `generated` one at a time — excluding
`generated` leaves only n=16 and is too small to read as real evidence either
way; the informative direction is whether excluding `seed` (leaving the full
n=233 generated set) shifts anything, which it doesn't for most descriptors.

**Scaffold-family version** (`output/plots/fig_scaffold/`):
```bash
python scripts/plot/plot_correlations.py output/scaffold_labeled.csv \
    --group-col scaffold_group --out-prefix output/plots/fig_scaffold/fig_scaffold
```
→ `fig_scaffold_correlations.csv`, `fig_scaffold_regression.csv`,
`fig_scaffold_heatmap.png`, `fig_scaffold_scatter_*.png` (same 8 descriptors)

Sensitivity check here excludes each of `scaffold_1`…`scaffold_6` one at a
time — confirms no single recurring chemotype is driving the global
correlations (as of the last run, none were).

### 4b. Scaffold diversity comparison (`plot_scaffold_diversity.py`)

Compares scaffold diversity between DiffSBDD-inpaint and REINVENT4-mol2mol
(source inferred from the first letter of `cmpd`: `i`=inpaint, `m`=mol2mol).
Produces a rank-abundance curve, singleton/top-scaffold summary bars, and a
rarefaction check (subsamples the larger pool down to the smaller pool's size,
to rule out sample-size as the cause of any apparent diversity gap).

```bash
python scripts/plot/plot_scaffold_diversity.py output/tc_merged_classified.csv \
    --out-prefix output/plots/diversity/diversity
```
→ `diversity_rank_abundance.png`, `diversity_summary_bars.png`,
`diversity_summary.csv`, `diversity_rarefaction.png`

Result on this dataset: no real diversity difference between the two tools —
the rarefied inpaint value sits at the ~48th percentile of the rarefied
mol2mol distribution, i.e. right where you'd expect by chance.

### 4c. Training-corpus bias check (`compare_training_bias.py`)

Tests whether the ChEMBL augmentation step (`scripts/tl/03_augment_from_chembl.py`
in the `reinvent4-mol2mol` repo) shifted the mol2mol training data toward
bigger/greasier compounds relative to the original seeds, and whether the
generated output tracks that shift. Uses Mann-Whitney U (robust) + KS as a
cross-check, since this is a pure structure/descriptor comparison (training
compounds were never docked).

```bash
python scripts/plot/compare_training_bias.py --references input/references.smi \
    --train input/compounds_train.smi --val input/compounds_val.smi \
    --generated output/tc_merged_classified.csv \
    --out-prefix output/plots/training_bias/training_bias
```
→ `training_bias_stats.csv`, `training_bias_boxplots.png` (jittered points
overlaid on each box — important given the seed group is only n=18),
`training_bias_labeled.csv`

Result on this dataset: a real, monotonic drift — seed median MW 286 →
ChEMBL-analog median MW 327 → generated-mol2mol median MW 383 — confirming
augmentation introduced a size/lipophilicity bias that the model then learned
and amplified further. `n_arom_rings` is the one descriptor that does *not*
show this pattern (analog vs. generated: p=0.46, not significant).

---

## Stage 5 — Structural diversity diagnostic (source-split, table only)

A lighter-weight version of 4b without RDKit descriptor recomputation — reads
directly off a `classify_seeds.py` output's `murcko` column:

```bash
python scripts/scaffold_by_source.py input/tc_merged.csv
```
Prints singleton%/top-scaffold% per source directly to stderr; no file output.
Source is derived the same way as 4b (`i`/`m` prefix on the compound ID).

---

## Quick reference: which output came from which command

| output file | generating command |
|---|---|
| `tc_merged_classified.csv`, `seeds_classified.csv`, `classified.csv`, `classified_0.15–0.40.csv`, `classified_sweep_summary.txt` | Stage 1 |
| `pooled.csv` | Stage 2 |
| `scaffold_labeled.csv`, `scaffold_groups.csv` | Stage 3 |
| `plots/pooled/*` | Stage 4a (pooled) |
| `plots/fig_scaffold/*` | Stage 4a (scaffold_group) |
| `plots/diversity/*` | Stage 4b |
| `plots/training_bias/*` | Stage 4c |

`scripts/redock/` (Vina redocking + PyMOL review) sits upstream of all of the
above — it's what produces `tc_merged.csv` in the first place — and is
documented separately.
