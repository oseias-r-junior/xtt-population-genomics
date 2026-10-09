#!/bin/bash
# NOTE: historical job script kept for provenance; paths were generalised for publication (run from the project root or set PROJECT_DIR).
#SBATCH --job-name=CFML_XMFA
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=cfml_xmfa_%j.out
#SBATCH --error=cfml_xmfa_%j.err

# ==========================================
# MOVE TO WORK DIRECTORY
# ==========================================

cd ${PROJECT_DIR:-$PWD}/cfml_xmfa_analysis

# ==========================================
# RUN CLONALFRAMEML
# ==========================================

${PROJECT_DIR:-$PWD}/third_party/ClonalFrameML/src/ClonalFrameML \
RAxML_bestTree.core_tree_best \
core_genes_corrected.xmfa \
cfml_xmfa \
-xmfa_file true
