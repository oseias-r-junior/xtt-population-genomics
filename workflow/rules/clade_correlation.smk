############################################################
# CLADE / CORRELATION MODULE (paper pipeline)
#
# Defense-system matrices, MGE x defense correlation, physical overlaps
# (defense systems and T3E/TSS inside MGEs) and the MGE distribution per
# clade. Scripts: defense_systems/scripts/ and
# integration/clade_correlation/scripts/.
#
# Raw DefenseFinder / PADLOC outputs and the SLURM-run BLAST/SSN steps are
# plain inputs (run outside Snakemake, same convention as ssn.smk).
############################################################

DF_RAW = "defense_systems/defensefinder/raw_outputs"
PADLOC_RAW = "defense_systems/padloc/results"
SSN_NODES = "integration/ssn/matrices/ssn_nodes_final.tsv"
T3E_TSS_XLSX = "metadata/Xtt_2021_2023_T3E_TSS (positions).xlsx"
CC = "integration/clade_correlation"

CLADE_CORR_TARGETS = (
    # defense systems
    "defense_systems/matrices/ds_consensus_gene_counts.tsv",
    "defense_systems/matrices/ds_consensus_presence.tsv",
    "defense_systems/matrices/ds_gene_level_matrix.tsv",
    "defense_systems/plots/ds_gene_heatmap_by_clade.png",
    # MGE per clade, MGE x DS correlation
    f"{CC}/plots/mge_distribution_by_clade.png",
    f"{CC}/matrices/mge_ds_correlation_long_presence.tsv",
    f"{CC}/plots/mge_ds_correlation_heatmap_v5_system_presence.png",
    # physical overlaps
    f"{CC}/matrices/ds_mge_physical_overlap_long_system.tsv",
    f"{CC}/plots/ds_mge_physical_overlap_heatmap_v2_system.png",
    f"{CC}/matrices/t3e_tss_mge_physical_overlap_long.tsv",
    f"{CC}/plots/t3e_tss_mge_overlap_heatmap.png",
)

############################################################
# DEFENSE SYSTEMS: consensus counts, presence, gene-level heatmap
############################################################

rule build_ds_consensus_matrix:
    input:
        df=DF_RAW,
        padloc=PADLOC_RAW
    output:
        wide="defense_systems/matrices/ds_consensus_gene_counts.tsv",
        long="defense_systems/matrices/ds_consensus_gene_counts_long.tsv"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: "defense_systems/logs/build_ds_consensus_matrix.log"
    message: "DefenseFinder + PADLOC consensus gene-count matrix (25 both-tool systems)."
    shell:
        "python defense_systems/scripts/build_ds_consensus_matrix.py > {log} 2>&1"

rule build_ds_presence_matrix:
    input:
        df=DF_RAW,
        padloc=PADLOC_RAW,
        counts="defense_systems/matrices/ds_consensus_gene_counts.tsv"
    output:
        wide="defense_systems/matrices/ds_consensus_presence.tsv",
        long="defense_systems/matrices/ds_consensus_instances_long.tsv"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: "defense_systems/logs/build_ds_presence_matrix.log"
    message: "Presence (0/1) of each defense system per genome, each system counted once."
    shell:
        "python defense_systems/scripts/build_ds_presence_matrix.py > {log} 2>&1"

rule plot_ds_gene_heatmap_by_clade:
    input:
        df=DF_RAW,
        padloc=PADLOC_RAW,
        clades=str(GENOME_CLADES),
        counts="defense_systems/matrices/ds_consensus_gene_counts.tsv"
    output:
        png="defense_systems/plots/ds_gene_heatmap_by_clade.png",
        matrix="defense_systems/matrices/ds_gene_level_matrix.tsv",
        completeness="defense_systems/matrices/ds_gene_level_system_completeness.tsv"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: "defense_systems/logs/plot_ds_gene_heatmap_by_clade.log"
    message: "Gene-by-gene defense-system heatmap grouped by clade."
    shell:
        "python defense_systems/scripts/plot_ds_gene_heatmap_by_clade.py > {log} 2>&1"

############################################################
# MGE DISTRIBUTION PER CLADE
############################################################

rule plot_mge_by_clade:
    input:
        nodes=SSN_NODES,
        clades=str(GENOME_CLADES)
    output:
        f"{CC}/plots/mge_distribution_by_clade.png"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: f"{CC}/logs/plot_mge_by_clade.log"
    message: "MGE composition per genome, grouped by clade (main figure)."
    shell:
        "python integration/clade_correlation/scripts/plot_mge_by_clade.py > {log} 2>&1"

############################################################
# MGE x DEFENSE CORRELATION
############################################################

rule build_mge_ds_correlation:
    input:
        nodes=SSN_NODES,
        presence="defense_systems/matrices/ds_consensus_presence.tsv",
        clades=str(GENOME_CLADES)
    output:
        f"{CC}/matrices/mge_ds_correlation_long_presence.tsv"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: f"{CC}/logs/build_mge_ds_correlation.log"
    message: "Spearman correlation MGE burden x defense burden (BH-FDR), pooled and per clade."
    shell:
        "python integration/clade_correlation/scripts/build_mge_ds_correlation.py > {log} 2>&1"

rule plot_mge_ds_correlation:
    input:
        f"{CC}/matrices/mge_ds_correlation_long_presence.tsv"
    output:
        f"{CC}/plots/mge_ds_correlation_heatmap_v5_system_presence.png"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: f"{CC}/logs/plot_mge_ds_correlation.log"
    shell:
        "python integration/clade_correlation/scripts/plot_mge_ds_correlation.py > {log} 2>&1"

############################################################
# PHYSICAL OVERLAP: defense-system loci inside MGEs
############################################################

rule build_ds_mge_physical_overlap:
    input:
        nodes=SSN_NODES,
        padloc=PADLOC_RAW,
        clades=str(GENOME_CLADES)
    output:
        f"{CC}/matrices/ds_mge_physical_overlap_long_system.tsv"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: f"{CC}/logs/build_ds_mge_physical_overlap.log"
    message: "Defense-system instances (PADLOC coordinates) overlapping MGE regions."
    shell:
        "python integration/clade_correlation/scripts/build_ds_mge_physical_overlap.py > {log} 2>&1"

rule plot_ds_mge_physical_overlap:
    input:
        f"{CC}/matrices/ds_mge_physical_overlap_long_system.tsv"
    output:
        f"{CC}/plots/ds_mge_physical_overlap_heatmap_v2_system.png"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: f"{CC}/logs/plot_ds_mge_physical_overlap.log"
    shell:
        "python integration/clade_correlation/scripts/plot_ds_mge_physical_overlap.py > {log} 2>&1"

############################################################
# PHYSICAL OVERLAP: T3E / secretion-system genes inside MGEs
# (needs openpyxl in the environment)
############################################################

rule build_t3e_tss_mge_overlap:
    input:
        xlsx=T3E_TSS_XLSX,
        nodes=SSN_NODES,
        clades=str(GENOME_CLADES)
    output:
        overlap=f"{CC}/matrices/t3e_tss_mge_physical_overlap_long.tsv",
        contig_map=f"{CC}/matrices/t3e_tss_chromosome_contig_map.tsv"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: f"{CC}/logs/build_t3e_tss_mge_overlap.log"
    message: "T3E and TSS gene hits (chromosome contig) overlapping MGE regions."
    shell:
        "python integration/clade_correlation/scripts/build_t3e_tss_mge_overlap.py > {log} 2>&1"

rule plot_t3e_tss_mge_overlap:
    input:
        f"{CC}/matrices/t3e_tss_mge_physical_overlap_long.tsv"
    output:
        f"{CC}/plots/t3e_tss_mge_overlap_heatmap.png"
    threads: 1
    conda: "../envs/xtt_blast.yml"
    log: f"{CC}/logs/plot_t3e_tss_mge_overlap.log"
    shell:
        "python integration/clade_correlation/scripts/plot_t3e_tss_mge_overlap.py > {log} 2>&1"