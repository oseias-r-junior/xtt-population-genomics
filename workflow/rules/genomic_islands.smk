############################################################
# GENOMIC ISLANDS
############################################################
#
# 3-way consensus: AlienHunter x IslandPath x PanISLE (recalibrated).
#
# AlienHunter and PanISLE raw outputs are parsed into region-level matrices
# here (IslandPath already has one: islandpath/matrices/islandpath_coordinates.tsv,
# from the pre-existing IslandPath pipeline). PanISLE itself is run manually
# via SLURM array job (mobilome/genomic_islands/panisle/jobs/run_panisle_array.slurm)
# against the isolated genomic_islands_next code (see THIRD_PARTY_LICENSES.md --
# code inspired by GIPSy2 is not distributed, so it lives outside the repo and is
# not orchestrated by Snakemake); this rule set picks up its results/ directory
# as a plain input once that job has been run.

############################################################
# TARGETS
############################################################

GENOMIC_ISLAND_TARGETS = (
    "mobilome/genomic_islands/matrices/gi_consensus_regions.tsv",
    "mobilome/genomic_islands/stats/gi_consensus_genome_summary.tsv",
)

############################################################
# ALIENHUNTER
############################################################

rule parse_alienhunter_islands:
    input:
        "mobilome/genomic_islands/alienhunter/raw_outputs"

    output:
        "mobilome/genomic_islands/alienhunter/matrices/alienhunter_islands.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/genomic_islands/alienhunter/logs/parse_alienhunter_islands.log"

    message:
        "Parsing AlienHunter's native (already thresholded + merged) island calls."

    shell:
        "python mobilome/genomic_islands/alienhunter/scripts/parse_alienhunter_islands.py"

############################################################
# PANISLE (results/ produced manually via SLURM array job, outside Snakemake --
# see mobilome/genomic_islands/panisle/jobs/run_panisle_array.slurm)
############################################################

rule parse_panisle_islands:
    input:
        "mobilome/genomic_islands/panisle/results"

    output:
        "mobilome/genomic_islands/panisle/matrices/panisle_islands.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/genomic_islands/panisle/logs/parse_panisle_islands.log"

    message:
        "Parsing PanISLE (recalibrated z-score, threshold=1.0 SD) island calls with functional classification."

    shell:
        "python mobilome/genomic_islands/panisle/scripts/parse_panisle_islands.py"

############################################################
# CONSENSUS
############################################################

rule build_gi_consensus:
    input:
        alienhunter="mobilome/genomic_islands/alienhunter/matrices/alienhunter_islands.tsv",
        islandpath="mobilome/genomic_islands/islandpath/matrices/islandpath_coordinates.tsv",
        panisle="mobilome/genomic_islands/panisle/matrices/panisle_islands.tsv"

    output:
        "mobilome/genomic_islands/matrices/gi_consensus_regions.tsv",
        "mobilome/genomic_islands/stats/gi_consensus_genome_summary.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/genomic_islands/logs/build_gi_consensus.log"

    message:
        "Building 3-way (AlienHunter x IslandPath x PanISLE) genomic island consensus with functional classification."

    shell:
        "python mobilome/genomic_islands/scripts/build_gi_consensus.py"
