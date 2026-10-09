#!/bin/bash
# NOTE: historical job script kept for provenance; paths were generalised for publication (run from the project root or set PROJECT_DIR).
#SBATCH --job-name=IQTREE_XMFA
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=iqtree_xmfa_%j.out
#SBATCH --error=iqtree_xmfa_%j.err

# ==========================================
# WORK DIRECTORY
# ==========================================

cd ${PROJECT_DIR:-$PWD}/input_data_hrr_new_phylo_tree/cfml_xmfa_analysis

# ==========================================
# IQ-TREE PATH
# ==========================================


# ==========================================
# RUN IQ-TREE
# ==========================================

for aln in \
    xmfa_full_masked.aln \
    xmfa_conservative_masked.aln
do

    prefix=$(basename "$aln" .aln)

    echo "======================================"
    echo "Running IQ-TREE on: $aln"
    echo "Prefix: $prefix"
    echo "======================================"

    iqtree2 \
        -s "$aln" \
        -m GTR+G \
        -bb 1000 \
        -nt AUTO \
        -pre "$prefix"

done
