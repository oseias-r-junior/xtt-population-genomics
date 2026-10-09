#!/bin/bash
#SBATCH --job-name=mobtest
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00

#SBATCH --output=xtt_project/mobilome/plasmids/mob_suite/logs/%x_%j_test.out
#SBATCH --error=xtt_project/mobilome/plasmids/mob_suite/logs/%x_%j_test.err

set -euo pipefail

export PATH=$HOME/bin:$PATH
eval "$(micromamba shell hook --shell=bash)"
micromamba activate xtt_mobsuite

PROJECT=${PROJECT_DIR:-$PWD}

cd $PROJECT

python mobilome/plasmids/mob_suite/scripts/run_mobsuite.py \
    --input annotations/fna \
    --output mobilome/plasmids/mob_suite/raw_outputs \
    --threads $SLURM_CPUS_PER_TASK \
    --limit 1