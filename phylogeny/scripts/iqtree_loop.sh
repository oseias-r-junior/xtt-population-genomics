#!/bin/bash

ALIGNMENTS=(
    "core_nonrecombinant.aln"
    "core_nonrecombinant_threshold.aln"
)

for aln in "${ALIGNMENTS[@]}"
do

    prefix=$(basename "$aln" .aln)

    echo "Running IQ-TREE for $aln"

    iqtree2 \
        -s "$aln" \
        -m GTR+G \
        -bb 1000 \
        -nt AUTO \
        --prefix "$prefix"

done
