############################################################
# INSERTION SEQUENCES
############################################################

############################################################
# TARGETS
############################################################

IS_TARGETS = (
    "mobilome/insertion_sequences/plots/is_family_heatmap.png",
    "mobilome/insertion_sequences/plots/is_family_pca.png",
    "mobilome/insertion_sequences/plots/is_total_barplot.png",
    "mobilome/insertion_sequences/plots/is_family_clustered_heatmap.png",
    "mobilome/insertion_sequences/plots/pca_is_clades.png",
    "mobilome/insertion_sequences/plots/mantel_is_vs_phylogeny.png",
    "mobilome/insertion_sequences/stats/is_clade_stats/kruskal_results.tsv",
    "mobilome/insertion_sequences/stats/enrichment/is_family_enrichment.tsv",
    "mobilome/insertion_sequences/itol/itol_clades.txt",
    "mobilome/insertion_sequences/plots/genome_size_vs_is.png",
)

############################################################
# MATRIX CONSTRUCTION
############################################################

rule build_is_matrix:
    output:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/build_is_matrix.log"
    message: "Building IS family presence/absence matrix."
    shell:
        "python mobilome/insertion_sequences/scripts/build_is_matrix.py"

############################################################
# SUMMARY ANALYSES
############################################################

rule is_profile_plots:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/plots/is_family_heatmap.png",
        "mobilome/insertion_sequences/plots/is_family_pca.png",
        "mobilome/insertion_sequences/plots/is_total_barplot.png"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/is_profile_plots.log"
    message: "Generating IS family descriptive plots."
    shell:
        "python mobilome/insertion_sequences/scripts/plot_is_profiles.py"

############################################################
# HEATMAPS
############################################################

rule is_clustered_heatmap:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/plots/is_family_clustered_heatmap.png"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/is_clustered_heatmap.log"
    message: "Plotting clustered IS family heatmap."
    shell:
        "python mobilome/insertion_sequences/scripts/plot_clustered_is_heatmap.py"

############################################################
# MULTIVARIATE ANALYSES
############################################################

rule is_pca_clades:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/plots/pca_is_clades.png"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/is_pca_clades.log"
    message: "Running PCA on IS profiles with clade overlay."
    shell:
        "python mobilome/insertion_sequences/scripts/pca_overlay_clades.py"

############################################################
# PHYLOGENETIC ANALYSES
############################################################

rule mantel_is_vs_phylogeny:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/plots/mantel_is_vs_phylogeny.png"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/mantel_is_vs_phylogeny.log"
    message: "Computing Mantel test between IS profiles and phylogeny."
    shell:
        "python mobilome/insertion_sequences/scripts/mantel_is_vs_phylogeny.py"

rule is_clade_stats:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/stats/is_clade_stats/kruskal_results.tsv",
        "mobilome/insertion_sequences/stats/is_clade_stats/is_family_boxplots_by_clade.png",
        "mobilome/insertion_sequences/stats/is_clade_stats/mean_is_by_clade.tsv"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/is_clade_stats.log"
    message: "Computing IS family statistics by clade."
    shell:
        "python mobilome/insertion_sequences/scripts/is_family_stats_by_clade.py"

############################################################
# ENRICHMENT ANALYSES
############################################################

rule is_family_enrichment:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/stats/enrichment/is_family_enrichment.tsv",
        "mobilome/insertion_sequences/stats/enrichment/is_family_clade_means_heatmap.png"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/is_family_enrichment.log"
    message: "Testing IS family enrichment across clades."
    shell:
        "python mobilome/insertion_sequences/scripts/is_family_enrichment.py"

############################################################
# iTOL
############################################################

rule build_itol_is_tracks:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/itol/itol_clades.txt",
        "mobilome/insertion_sequences/itol/itol_total_is.txt",
        "mobilome/insertion_sequences/itol/itol_is_heatmap.txt"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/build_itol_is_tracks.log"
    message: "Building iTOL tracks for IS profiles."
    shell:
        "python mobilome/insertion_sequences/scripts/itol_is_tracks.py"

############################################################
# GENOME SIZE
############################################################

rule genome_size_vs_is:
    input:
        "mobilome/insertion_sequences/matrices/is_family_matrix.tsv"
    output:
        "mobilome/insertion_sequences/plots/genome_size_vs_is.png"
    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "mobilome/insertion_sequences/logs/genome_size_vs_is.log"
    message: "Analyzing relationship between genome size and IS content."
    shell:
        "python mobilome/insertion_sequences/scripts/genome_size_vs_is.py"
