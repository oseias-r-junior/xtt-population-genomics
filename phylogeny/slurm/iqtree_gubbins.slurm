#!/bin/bash
# NOTE: historical job script kept for provenance; paths were generalised for publication (run from the project root or set PROJECT_DIR).
#SBATCH --job-name=IQTREE_GUBBINS
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=iqtree_gubbins_%j.out
#SBATCH --error=iqtree_gubbins_%j.err

cd ${PROJECT_DIR:-$PWD}/input_data_hrr_new_phylo_tree

iqtree2 \
-s core_gubbins_masked.aln \
-m GTR+G \
-B 1000 \
-T AUTO
