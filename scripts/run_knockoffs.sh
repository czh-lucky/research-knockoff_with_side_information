#!/bin/bash

set -e

input_file="$1"
B="$2"

for b in $(seq 1 "$B"); do
  
  echo "================"
  echo "Running MC iteration: $b"
  echo "================"
  
  Rscript scripts/01_generate_Xk.R "$b" "$input_file"
  
done

echo "All knockoff simulations finished"