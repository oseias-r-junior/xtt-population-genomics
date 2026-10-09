############################################################
# SSN (Sequence Similarity Network)
#
# MGE nucleotide-level network mirroring Uceda-Campos et al. 2022
# (Microorganisms) Fig. 5: nodes = MGE instances (PPH/GI/IS/PLS), edges =
# BLASTN similarity >=50% identity / >=80% coverage.
#
# The all-vs-all BLASTN step itself is run manually via SLURM (see
# integration/ssn/jobs/run_blastn_ssn.slurm) rather than as a Snakemake
# rule -- same convention as PanISLE (mobilome/genomic_islands/panisle):
# a single heavy/long job better suited to its own SLURM allocation than
# to Snakemake's per-rule scheduling. This rule set picks up its filtered
# output (blast/all_vs_all_filtered.tsv) as a plain input once that job
# has been run.
############################################################

############################################################
# TARGETS
############################################################

SSN_TARGETS = (
    "integration/ssn/matrices/mge_ssn.graphml",
    "integration/ssn/matrices/ssn_nodes_final.tsv",
    "integration/ssn/matrices/ssn_edges_final.tsv",
)

############################################################
# EXTRACT MGE SEQUENCES
############################################################

rule extract_mge_sequences:
    input:
        prophages="mobilome/prophages/matrices/prophage_consensus_regions.tsv",
        gi="mobilome/genomic_islands/matrices/gi_consensus_regions.tsv",
        plasmids="mobilome/plasmids/matrices/plasmid_consensus.tsv"

    output:
        fasta="integration/ssn/fasta/all_mge_sequences.fna",
        nodes="integration/ssn/matrices/ssn_nodes.tsv"

    threads:
        1

    conda:
        "../envs/xtt_blast.yml"

    log:
        "integration/ssn/logs/extract_mge_sequences.log"

    message:
        "Extracting nucleotide sequences for all consensus MGE instances (PPH/GI/PLS/IS) for the SSN."

    shell:
        "python integration/ssn/scripts/extract_mge_sequences.py"

############################################################
# BUILD NETWORK (BLASTN all-vs-all run externally via SLURM, see
# integration/ssn/jobs/run_blastn_ssn.slurm -- picked up here as a plain
# input once available)
############################################################

rule build_ssn_network:
    input:
        nodes="integration/ssn/matrices/ssn_nodes.tsv",
        blast_filtered="integration/ssn/blast/all_vs_all_filtered.tsv"

    output:
        graphml="integration/ssn/matrices/mge_ssn.graphml",
        nodes_final="integration/ssn/matrices/ssn_nodes_final.tsv",
        edges_final="integration/ssn/matrices/ssn_edges_final.tsv"

    threads:
        4

    conda:
        "../envs/xtt_blast.yml"

    log:
        "integration/ssn/logs/build_ssn_network.log"

    message:
        "Building the MGE Sequence Similarity Network (nodes, edges, PPH-G/GI-G/IS-G/PLS-G groups)."

    shell:
        "python integration/ssn/scripts/build_ssn_network.py"
