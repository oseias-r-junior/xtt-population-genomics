#!/bin/bash
#SBATCH --job-name=isescan
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=xtt_project/mobilome/insertion_sequences/isescan/logs/isescan_%A_%a.out
#SBATCH --error=xtt_project/mobilome/insertion_sequences/isescan/logs/isescan_%A_%a.err

set -euo pipefail

PROJECT="${PROJECT_DIR:-$PWD}"
SCRIPT="$PROJECT/xtt_project/mobilome/insertion_sequences/isescan/scripts/run_isescan.py"

eval "$(micromamba shell hook --shell bash)"
micromamba activate xtt_isescan

python \
"$SCRIPT" \
    --input  "$PROJECT/xtt_project/annotations/fna" \
    --output "$PROJECT/xtt_project/mobilome/insertion_sequences/isescan/results" \
    --threads "${SLURM_CPUS_PER_TASK:-4}"