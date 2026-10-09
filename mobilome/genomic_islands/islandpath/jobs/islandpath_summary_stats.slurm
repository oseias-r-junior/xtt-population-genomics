#!/bin/bash
#SBATCH --job-name=islandpath_stats
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --output=xtt_project/mobilome/genomic_islands/islandpath/logs/%x_%j.out
#SBATCH --error=xtt_project/mobilome/genomic_islands/islandpath/logs/%x_%j.err

set -euo pipefail

PROJECT=$(git rev-parse --show-toplevel)
cd "$PROJECT"

SCRIPT=$PROJECT/xtt_project/mobilome/genomic_islands/islandpath/scripts/islandpath_summary_stats.py

echo "Running IslandPath summary stats..."

PYTHONPATH="$PROJECT" micromamba run -n xtt_islandpath \
    python "${SCRIPT}"

echo "Finished IslandPath summary stats."