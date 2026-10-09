#!/bin/bash
#SBATCH --job-name=islandpath_matrices
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --output=xtt_project/mobilome/genomic_islands/islandpath/logs/%x_%j.out
#SBATCH --error=xtt_project/mobilome/genomic_islands/islandpath/logs/%x_%j.err

set -euo pipefail

PROJECT=$(git rev-parse --show-toplevel)
cd "$PROJECT"

SCRIPT=xtt_project/mobilome/genomic_islands/islandpath/scripts/build_islandpath_matrix.py

echo "Building IslandPath summary matrices..."

# Garantir que o Python veja a raiz do projeto (onde está 'workflow/')
PYTHONPATH="$PROJECT" micromamba run -n xtt_islandpath \
    python "${SCRIPT}"

echo "Finished building IslandPath matrices."
