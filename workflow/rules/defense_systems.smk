############################################################
# DEFENSE SYSTEMS
#
# Internal analyses only
############################################################

############################################################
# FINAL TARGETS
############################################################

DEFENSE_TARGETS = (
    "defense_systems/defensefinder/stats/defense_prevalence.tsv",
    "defense_systems/defensefinder/stats/defense_systems_per_genome.tsv",
    "defense_systems/defensefinder/plots/defense_prevalence.png",
    "defense_systems/defensefinder/plots/defense_clade_heatmap.png",
    "defense_systems/defensefinder/plots/defense_pca.png",
    "defense_systems/defensefinder/plots/defense_cooccurrence_network.png",
    "defense_systems/defensefinder/stats/cooccurrence/defense_cooccurrence_centrality.tsv",
)

############################################################
# MATRIX CONSTRUCTION
############################################################

rule build_defensefinder_matrix:
    output:
        "defense_systems/defensefinder/matrices/defensefinder_matrix.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/build_defensefinder_matrix.log"

    message:
        "Building DefenseFinder presence/absence matrix."

    shell:
        "python defense_systems/defensefinder/scripts/build_defensefinder_matrix.py"


############################################################
# SUMMARY ANALYSES
############################################################

rule defense_summary_stats:
    input:
        "defense_systems/defensefinder/matrices/defensefinder_matrix.tsv"

    output:
        "defense_systems/defensefinder/stats/defense_prevalence.tsv",
        "defense_systems/defensefinder/stats/defense_systems_per_genome.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_summary_stats.log"

    message:
        "Computing DefenseFinder summary statistics."

    shell:
        "python defense_systems/defensefinder/scripts/defense_summary_stats.py"


rule defense_prevalence_barplot:
    input:
        "defense_systems/defensefinder/stats/defense_prevalence.tsv"

    output:
        "defense_systems/defensefinder/plots/defense_prevalence.png"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_prevalence_barplot.log"

    message:
        "Plotting DefenseFinder prevalence barplot."

    shell:
        "python defense_systems/defensefinder/scripts/defense_prevalence_barplot.py"


############################################################
# PHYLOGENETIC ANALYSES
############################################################

rule defense_clade_stats:
    input:
        "defense_systems/defensefinder/matrices/defensefinder_matrix.tsv"

    output:
        "defense_systems/defensefinder/stats/defense_clade_stats/mean_defense_by_clade.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_clade_stats.log"

    message:
        "Computing DefenseFinder clade-level statistics."

    shell:
        "python defense_systems/defensefinder/scripts/defense_clade_stats.py"


rule defense_clade_heatmap:
    input:
        "defense_systems/defensefinder/stats/defense_clade_stats/mean_defense_by_clade.tsv"

    output:
        "defense_systems/defensefinder/plots/defense_clade_heatmap.png"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_clade_heatmap.log"

    message:
        "Plotting DefenseFinder clade heatmap."

    shell:
        "python defense_systems/defensefinder/scripts/defense_clade_heatmap.py"


############################################################
# PCA
############################################################

rule defense_pca:
    input:
        "defense_systems/defensefinder/matrices/defensefinder_matrix.tsv"

    output:
        "defense_systems/defensefinder/plots/defense_pca.png"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_pca.log"

    message:
        "Running PCA on DefenseFinder matrix."

    shell:
        "python defense_systems/defensefinder/scripts/defense_pca.py"


############################################################
# CORRELATION ANALYSES
############################################################

rule defense_correlation:
    input:
        "defense_systems/defensefinder/matrices/defensefinder_matrix.tsv"

    output:
        "defense_systems/defensefinder/stats/defense_correlation.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_correlation.log"

    message:
        "Computing DefenseFinder correlation matrix."

    shell:
        "python defense_systems/defensefinder/scripts/defense_correlation.py"


############################################################
# CO-OCCURRENCE ANALYSES
############################################################

rule defense_cooccurrence:
    input:
        "defense_systems/defensefinder/matrices/defensefinder_matrix.tsv"

    output:
        "defense_systems/defensefinder/stats/cooccurrence/defense_cooccurrence_significant.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_cooccurrence.log"

    message:
        "Computing DefenseFinder co-occurrence statistics."

    shell:
        "python defense_systems/defensefinder/scripts/defense_cooccurrence.py"


rule defense_cooccurrence_network:
    input:
        "defense_systems/defensefinder/stats/cooccurrence/defense_cooccurrence_significant.tsv"

    output:
        "defense_systems/defensefinder/plots/defense_cooccurrence_network.png",
        "defense_systems/defensefinder/stats/cooccurrence/defense_cooccurrence_edges.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_cooccurrence_network.log"

    message:
        "Plotting DefenseFinder co-occurrence network."

    shell:
        "python defense_systems/defensefinder/scripts/defense_cooccurrence_network.py"


rule defense_cooccurrence_centrality:
    input:
        "defense_systems/defensefinder/stats/cooccurrence/defense_cooccurrence_edges.tsv"

    output:
        "defense_systems/defensefinder/stats/cooccurrence/defense_cooccurrence_centrality.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "defense_systems/logs/defense_cooccurrence_centrality.log"

    message:
        "Computing centrality metrics for DefenseFinder co-occurrence network."

    shell:
        "python defense_systems/defensefinder/scripts/defense_cooccurrence_centrality.py"
