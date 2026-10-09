#!/bin/bash
#SBATCH --job-name=islandpath
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00

#SBATCH --output=xtt_project/mobilome/genomic_islands/islandpath/logs/%x_%j.out
#SBATCH --error=xtt_project/mobilome/genomic_islands/islandpath/logs/%x_%j.err

set -euo pipefail

PROJECT=$(git rev-parse --show-toplevel)
cd "$PROJECT"

GBK_DIR=$PROJECT/xtt_project/annotations/gbk
OUTDIR=$PROJECT/xtt_project/mobilome/genomic_islands/islandpath/raw_outputs
SCRIPT=xtt_project/mobilome/genomic_islands/islandpath/scripts/run_islandpath.py

echo "Running IslandPath-DIMOB..."
echo "GenBank directory : ${GBK_DIR}"
echo "Output directory  : ${OUTDIR}"

for genome in "${GBK_DIR}"/*.gbk; do
    name=$(basename "${genome}")
    echo "Processing genome: ${name}"

    micromamba run -n xtt_islandpath \
        python "${SCRIPT}" \
            --genome "${genome}" \
            --outdir "${OUTDIR}"
done

echo "Finished IslandPath-DIMOB for all genomes."