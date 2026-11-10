#!/usr/bin/env bash
set -euo pipefail

# Run from any directory: bash scripts/run_vanilla_knockoff.sh
cd "$(dirname "$0")/.."

# Edit the parameter ranges here. Only reads existing Xk files.
sample_sizes=(100 400 581)
target_ids=({1..12})
amplitudes=(20 3.5)

for n in "${sample_sizes[@]}"; do
  B=10
  if [ "$n" -eq 581 ]; then
    B=1
  fi
  for id in "${target_ids[@]}"; do
    for amp in "${amplitudes[@]}"; do
      for ((b = 1; b <= B; b++)); do
        echo "Running n=$n B=$B id=$id amp=$amp b=$b"
        Rscript scripts/02_vanilla_knockoff_batch.R "$b" "$id" "$amp" "$n" "$B"
      done
    done
  done
done

echo "All vanilla knockoff experiments finished."
