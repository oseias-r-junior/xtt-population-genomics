#!/bin/bash
#SBATCH --job-name=padloc
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=xtt_project/defense_systems/padloc/logs/padloc_%j.out
#SBATCH --error=xtt_project/defense_systems/padloc/logs/padloc_%j.err

set -euo pipefail

###############################################################################
# Project paths
###############################################################################

PROJECT="${PROJECT_DIR:-$PWD}"

SCRIPT="$PROJECT/xtt_project/defense_systems/padloc/scripts/run_padloc.py"

INPUT="$PROJECT/xtt_project/annotations/fna"

OUTPUT="$PROJECT/xtt_project/defense_systems/padloc/results"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

###############################################################################
# Activate environment
###############################################################################

eval "$(micromamba shell hook --shell bash)"

micromamba activate xtt_padloc

###############################################################################
# Create directories
###############################################################################

mkdir -p \
    "$PROJECT/xtt_project/defense_systems/padloc/logs" \
    "$OUTPUT"

###############################################################################
# Run PADLOC
###############################################################################

python "$SCRIPT" \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --threads "$THREADS" \
    --max-genomes 5