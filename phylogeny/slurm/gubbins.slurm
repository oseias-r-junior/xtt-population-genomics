#!/bin/bash
# NOTE: historical job script kept for provenance; paths were generalised for publication (run from the project root or set PROJECT_DIR).
#SBATCH --job-name=GUBBINS
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=gubbins_%j.out
#SBATCH --error=gubbins_%j.err

cd ${PROJECT_DIR:-$PWD}/input_data_hrr_new_phylo_tree

# ==========================================
# MICROMAMBA
# ==========================================

eval "$(~/bin/micromamba shell hook -s bash)"

micromamba activate gubbins_env

# ==========================================
# RUN GUBBINS
# ==========================================

run_gubbins.py \
    --threads 8 \
    --prefix gubbins_core \
    core_gene_alignment.aln
