#!/bin/bash
#SBATCH --job-name=build_ssn_network
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=integration/ssn/logs/%x_%j.out
#SBATCH --error=integration/ssn/logs/%x_%j.err

set -euo pipefail

PROJECT=${PROJECT_DIR:-$PWD}

eval "$(micromamba shell hook --shell bash)"
micromamba activate xtt_blast

cd "$PROJECT"
python integration/ssn/scripts/build_ssn_network.py
