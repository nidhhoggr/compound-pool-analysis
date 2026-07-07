#!/usr/bin/env bash
set -euo pipefail

# Sweep classify_seeds.py across multiple Butina --sim thresholds, so cluster-size
# distributions can be compared side by side. Low thresholds cause single-linkage
# chaining (structurally unrelated scaffolds merge into one giant cluster); this
# sweep is how that gets diagnosed before picking a --sim value to trust.
#
# Usage:
#   scripts/multithreshold-classifier.sh [input_filename] [output_dir]
#
# input_filename is resolved relative to $INPUT_DIR (default: input/) and
# defaults to tc_merged.csv. output_dir defaults to $OUTPUT_DIR (default:
# output/) -- both overridable via env var or positional argument.
#
# Thresholds swept default to 0.40 0.35 0.30 0.25 0.20 0.15; override with the
# THRESHOLDS env var (space-separated), e.g.:
#   THRESHOLDS="0.5 0.4 0.3" scripts/multithreshold-classifier.sh tc_merged.csv
#
# Column names default to ID/docking/SMILES (tc_merged.csv's header) --
# override with ID_COL/SCORE_COL/SMILES_COL env vars for a differently-headed
# CSV input.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_DIR="${INPUT_DIR:-input}"
OUTPUT_DIR="${OUTPUT_DIR:-output}"
THRESHOLDS="${THRESHOLDS:-0.40 0.35 0.30 0.25 0.20 0.15}"
ID_COL="${ID_COL:-ID}"
SCORE_COL="${SCORE_COL:-docking}"
SMILES_COL="${SMILES_COL:-SMILES}"

INPUT_FILE="${1:-tc_merged.csv}"
OUT_DIR="${2:-$OUTPUT_DIR}"

mkdir -p "$OUT_DIR"

SUMMARY_FILE="$OUT_DIR/classified_sweep_summary.txt"
: > "$SUMMARY_FILE"   # truncate/create fresh each run

for s in $THRESHOLDS; do
  echo "--- sim=$s ---" | tee -a "$SUMMARY_FILE"
  python "$SCRIPT_DIR/classify_seeds.py" \
      --sim "$s" \
      --id-col "$ID_COL" --score-col "$SCORE_COL" --smiles-col "$SMILES_COL" \
      -o "$OUT_DIR/classified_$s.csv" \
      "$INPUT_DIR/$INPUT_FILE" \
      2>> "$SUMMARY_FILE"
done

echo ""
echo "wrote $OUT_DIR/classified_<threshold>.csv for thresholds: $THRESHOLDS"
echo "per-threshold compound/scaffold/cluster-size summary: $SUMMARY_FILE"
