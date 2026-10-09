#!/bin/bash
#SBATCH --job-name=mge_by_clade_genomesize
#SBATCH --output=integration/clade_correlation/logs/%x_%j.out
#SBATCH --error=integration/clade_correlation/logs/%x_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --partition=compute

source ~/.bashrc
cd ${PROJECT_DIR:-$PWD}
python3 -c "import matplotlib, pandas, numpy; print('env OK')" \
    || { echo "FATAL: pacote faltando ($(which python3))"; exit 1; }
python3 integration/clade_correlation/scripts/plot_mge_by_clade_genomesize.py
