#!/bin/bash
# NOTE: historical job script kept for provenance; paths were generalised for publication (run from the project root or set PROJECT_DIR).
#SBATCH --job-name=CFML
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=cfml_%j.out
#SBATCH --error=cfml_%j.err

cd ${PROJECT_DIR:-$PWD}/input_data_hrr_new_phylo_tree

${PROJECT_DIR:-$PWD}/third_party/ClonalFrameML/src/ClonalFrameML \
RAxML_bestTree.core_tree_best \
core_gene_alignment.aln \
cfml_output
