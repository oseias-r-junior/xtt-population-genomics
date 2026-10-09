#!/bin/bash
#SBATCH --job-name=xtt_defensefinder
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=defense_systems/logs/defensefinder_%j.out
#SBATCH --error=defense_systems/logs/defensefinder_%j.err

set -euxo pipefail

PROJECT_ROOT="${PROJECT_DIR:-$PWD}"

echo "PROJECT_ROOT=${PROJECT_ROOT}"

cd "${PROJECT_ROOT}"

pwd

ls

eval "$(micromamba shell hook --shell bash)"
micromamba activate xtt_defense

# mkdir -p defense_systems/defensefinder/raw_outputs

echo "PWD = $(pwd)"
echo "Project root = ${PROJECT_ROOT}"

ls annotations/fna | head

for GENOME in annotations/fna/*.fna
do

    BASENAME=$(basename "${GENOME}" .fna)

    OUTDIR="defense_systems/defensefinder/raw_outputs/${BASENAME}"

    echo "======================================="
    echo "Running DefenseFinder on ${BASENAME}"
    echo "======================================="

    defense-finder run \
        "${GENOME}" \
        -o "${OUTDIR}" \
        -w ${SLURM_CPUS_PER_TASK}

done

echo "Finished DefenseFinder run."