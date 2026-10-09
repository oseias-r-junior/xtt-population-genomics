############################################################
# PROPHAGES
############################################################

############################################################
# TARGETS
############################################################

PROPHAGE_TARGETS = (
    "mobilome/prophages/matrices/prophage_summary.tsv",
    "mobilome/prophages/stats/prophage_clade_presence.tsv",
    "mobilome/prophages/stats/prophage_kruskal.tsv",
    "mobilome/prophages/virsorter2/stats/virsorter2_summary.tsv",
    "mobilome/prophages/matrices/prophage_consensus_regions.tsv",
    "mobilome/prophages/stats/prophage_consensus_summary.tsv",
)

############################################################
# PHASTEST
############################################################

rule prophage_summary:
    input:
        "mobilome/prophages/matrices/prophage_family_presence_absence.tsv"

    output:
        "mobilome/prophages/matrices/prophage_summary.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/prophage_summary.log"

    message:
        "Building prophage summary statistics."

    shell:
        "python mobilome/prophages/scripts/prophage_summary_stats.py"


############################################################
# CLADE ANALYSES
############################################################

rule prophage_clade_stats:
    input:
        summary="mobilome/prophages/matrices/prophage_summary.tsv",
        clades="metadata/genome_clades_master.tsv"

    output:
        "mobilome/prophages/stats/prophage_clade_summary.tsv",
        "mobilome/prophages/stats/prophage_clade_presence.tsv",
        "mobilome/prophages/stats/prophage_kruskal.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/prophage_clade_stats.log"

    message:
        "Computing prophage clade statistics."

    shell:
        "python mobilome/prophages/scripts/stats_prophages_by_clade.py"


############################################################
# VIRSORTER2
############################################################

rule build_virsorter2_matrix:
    input:
        expand(
            "mobilome/prophages/virsorter2/raw_outputs/{genome}/final-viral-score.tsv",
            genome=GENOMES
        )

    output:
        "mobilome/prophages/virsorter2/matrices/virsorter2_summary_matrix.tsv",
        "mobilome/prophages/virsorter2/matrices/virsorter2_sequence_matrix.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/build_virsorter2_matrix.log"

    message:
        "Building VirSorter2 matrices."

    shell:
        "python mobilome/prophages/virsorter2/scripts/build_virsorter2_matrix.py"


rule virsorter2_summary:
    input:
        "mobilome/prophages/virsorter2/matrices/virsorter2_summary_matrix.tsv"

    output:
        "mobilome/prophages/virsorter2/stats/virsorter2_summary.tsv",
        "mobilome/prophages/virsorter2/plots/virsorter2_prevalence.png",
        "mobilome/prophages/virsorter2/plots/virsorter2_length_distribution.png",
        "mobilome/prophages/virsorter2/plots/virsorter2_score_distribution.png",
        "mobilome/prophages/virsorter2/plots/virsorter2_sequence_types.png"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/virsorter2_summary.log"

    message:
        "Summarizing VirSorter2 predictions."

    shell:
        "python mobilome/prophages/virsorter2/scripts/virsorter2_summary_stats.py"


rule virsorter2_clade_stats:
    input:
        summary="mobilome/prophages/virsorter2/matrices/virsorter2_summary_matrix.tsv",
        clades="metadata/genome_clades_master.tsv"

    output:
        "mobilome/prophages/virsorter2/stats/virsorter2_clade_summary.tsv",
        "mobilome/prophages/virsorter2/stats/virsorter2_clade_stats.tsv",
        "mobilome/prophages/virsorter2/plots/virsorter2_sequences_by_clade.png",
        "mobilome/prophages/virsorter2/plots/virsorter2_length_by_clade.png",
        "mobilome/prophages/virsorter2/plots/virsorter2_score_by_clade.png",
        "mobilome/prophages/virsorter2/plots/virsorter2_longest_by_clade.png"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/virsorter2_clade_stats.log"

    message:
        "Computing VirSorter2 clade statistics."

    shell:
        "python mobilome/prophages/virsorter2/scripts/virsorter2_clade_stats.py"


############################################################
# DETECTOR COMPARISON
############################################################

rule virsorter2_vs_phastest:
    input:
        phastest="mobilome/prophages/matrices/prophage_summary.tsv",
        virsorter="mobilome/prophages/virsorter2/matrices/virsorter2_summary_matrix.tsv"

    output:
        "mobilome/prophages/stats/virsorter2_vs_phastest.tsv",
        "mobilome/prophages/plots/virsorter2_vs_phastest_scatter.png",
        "mobilome/prophages/plots/virsorter2_vs_phastest_boxplot.png"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/virsorter2_vs_phastest.log"

    message:
        "Comparing VirSorter2 and PHASTEST."

    shell:
        "python mobilome/prophages/scripts/virsorter2_vs_phastest.py"


rule virsorter2_phastest_disagreement:
    input:
        phastest="mobilome/prophages/matrices/prophage_summary.tsv",
        virsorter="mobilome/prophages/virsorter2/matrices/virsorter2_summary_matrix.tsv"

    output:
        "mobilome/prophages/stats/virsorter2_phastest_disagreement.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/virsorter2_phastest_disagreement.log"

    message:
        "Identifying disagreements between VirSorter2 and PHASTEST."

    shell:
        "python mobilome/prophages/scripts/virsorter2_phastest_disagreement.py"


############################################################
# CONSENSUS
############################################################

# DEPRECATED 2026-09-25: replaced by the 3-way region-level consensus below
# (PHASTEST x VirSorter2 x geNomad, built from scratch with no 2-tool
# pre-consensus). Archived script:
# mobilome/prophages/scripts/_archive/build_prophage_consensus_matrix.py
#
# rule build_prophage_consensus_matrix:
#     input:
#         phastest="mobilome/prophages/matrices/prophage_summary.tsv",
#         virsorter="mobilome/prophages/virsorter2/matrices/virsorter2_summary_matrix.tsv"
#     output:
#         "mobilome/prophages/matrices/prophage_consensus_matrix.tsv",
#         "mobilome/prophages/stats/prophage_consensus_details.tsv"
#     threads:
#         1
#     conda:
#         "../envs/xtt_core.yml"
#     log:
#         "mobilome/prophages/logs/build_prophage_consensus_matrix.log"
#     message:
#         "Building consensus matrix from PHASTEST and VirSorter2."
#     shell:
#         "python mobilome/prophages/scripts/build_prophage_consensus_matrix.py"


rule parse_virsorter2_boundaries:
    input:
        "mobilome/prophages/virsorter2/raw_outputs"

    output:
        "mobilome/prophages/virsorter2/matrices/virsorter2_regions.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/parse_virsorter2_boundaries.log"

    message:
        "Parsing VirSorter2 per-region genomic boundaries."

    shell:
        "python mobilome/prophages/virsorter2/scripts/parse_virsorter2_boundaries.py"


rule build_prophage_consensus_regions:
    input:
        phastest_regions="mobilome/prophages/matrices/prophage_regions.tsv",
        virsorter2_regions="mobilome/prophages/virsorter2/matrices/virsorter2_regions.tsv",
        genomad_raw="mobilome/prophages/genomad/raw_outputs"

    output:
        "mobilome/prophages/matrices/prophage_consensus_regions.tsv",
        "mobilome/prophages/stats/prophage_consensus_genome_summary.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/build_prophage_consensus_regions.log"

    message:
        "Building 3-way (PHASTEST x VirSorter2 x geNomad) region-level prophage consensus."

    shell:
        "python mobilome/prophages/scripts/build_prophage_consensus_regions.py"


rule prophage_consensus_summary:
    input:
        "mobilome/prophages/matrices/prophage_consensus_regions.tsv",
        "mobilome/prophages/stats/prophage_consensus_genome_summary.tsv"

    output:
        "mobilome/prophages/stats/prophage_consensus_summary.tsv",
        "mobilome/prophages/plots/prophage_consensus_regions_prevalence.png",
        "mobilome/prophages/plots/prophage_consensus_support_level.png",
        "mobilome/prophages/plots/prophage_consensus_classification.png",
        "mobilome/prophages/plots/prophage_consensus_tools_combination.png"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/prophages/logs/prophage_consensus_summary.log"

    message:
        "Summarizing 3-way region-level consensus prophage predictions."

    shell:
        "python mobilome/prophages/scripts/prophage_consensus_summary_stats.py"
