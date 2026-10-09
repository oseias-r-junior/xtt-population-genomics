#!/bin/bash
# NOTE: historical job script kept for provenance; paths were generalised for publication (run from the project root or set PROJECT_DIR).
#SBATCH --job-name=IQTREE
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=iqtree_%j.out
#SBATCH --error=iqtree_%j.err


which iqtree2
iqtree2 --version

cd ${PROJECT_DIR:-$PWD}/input_data_hrr_new_phylo_tree

bash ${PROJECT_DIR:-$PWD}/phylogeny/scripts/iqtree_loop.sh
