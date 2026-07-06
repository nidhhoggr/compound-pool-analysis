mkdir -p /workspace/output/classified

INPUT_FILE=${1:-merged}

for s in 0.40 0.35 0.30 0.25 0.20 0.15; do
  echo "--- sim=$s ---"
  python classify_seeds.py --sim $s --id-col ID --score-col docking --smiles-col SMILES -o /workspace/output/classified/classified_$s.csv /workspace/input/$INPUT_FILE
done
