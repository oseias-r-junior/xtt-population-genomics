#!/bin/bash
#SBATCH --job-name=ssn_viz_v9
#SBATCH --output=integration/ssn/logs/%x_%j.out
#SBATCH --error=integration/ssn/logs/%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --partition=compute

source ~/.bashrc
cd ${PROJECT_DIR:-$PWD}
python3 -c "import matplotlib, networkx, scipy; from adjustText import adjust_text; print('env OK')" \
    || { echo "FATAL: pacote faltando ($(which python3))"; exit 1; }
python3 integration/ssn/scripts/visualize_ssn_network_v9.py
