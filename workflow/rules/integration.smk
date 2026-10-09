############################################################
# INTEGRATION ANALYSES
#
# Cross-module analyses
############################################################

############################################################
# FINAL TARGETS
############################################################

INTEGRATION_TARGETS = (

    "integration/stats/defense_vs_is_total.tsv",
    "integration/stats/is_association/defense_vs_is_family.tsv",
    "integration/stats/clade_corrected/defense_is_family_clade_corrected.tsv",
    "integration/plots/defense_vs_is_network.png",
    "integration/stats/defense_vs_prophage.tsv",
    "integration/stats/clade_corrected/defense_vs_prophage_clade_corrected.tsv",
    "integration/stats/clade_corrected/defense_vs_mobilome.tsv",
    "integration/stats/clade_corrected/defense_vs_mobilome_scaled.tsv",
    "integration/stats/clade_corrected/defense_vs_mobilome_centrality.tsv",
    "integration/plots/mantel_defense_vs_phylogeny.png",

)

############################################################
# DEFENSE × PROPHAGES
############################################################

rule defense_vs_prophage:
    input:
        defense="defense_systems/defensefinder/matrices/defensefinder_matrix.tsv",
        prophages="mobilome/prophages/matrices/prophage_consensus_matrix.tsv"

    output:
        "integration/stats/defense_vs_prophage.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_prophage.log"
    message: "Computing Defense × Prophage associations."
    shell:
        "python integration/scripts/defense_vs_prophage.py"


rule defense_vs_prophage_clade_corrected:
    input:
        defense="defense_systems/defensefinder/matrices/defensefinder_matrix.tsv",
        prophages="mobilome/prophages/matrices/prophage_consensus_matrix.tsv"

    output:
        "integration/stats/clade_corrected/defense_vs_prophage_clade_corrected.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_prophage_clade_corrected.log"
    message: "Computing clade-corrected Defense × Prophage associations."
    shell:
        "python integration/scripts/defense_vs_prophage_clade_corrected.py"


############################################################
# DEFENSE × INSERTION SEQUENCES
############################################################

rule defense_vs_is_total:
    input:
        defense="defense_systems/defensefinder/matrices/defensefinder_matrix.tsv",
        is_matrix="mobilome/insertion_sequences/matrices/is_family_matrix.tsv"

    output:
        "integration/stats/defense_vs_is_total.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_is_total.log"
    message: "Computing total Defense × IS associations."
    shell:
        "python integration/scripts/defense_vs_is_total.py"


rule defense_vs_is_family:
    input:
        defense="defense_systems/defensefinder/matrices/defensefinder_matrix.tsv",
        is_matrix="mobilome/insertion_sequences/matrices/is_family_matrix.tsv"

    output:
        "integration/stats/is_association/defense_vs_is_family.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_is_family.log"
    message: "Computing Defense × IS family associations."
    shell:
        "python integration/scripts/defense_vs_is_family.py"


rule defense_vs_is_family_clade_corrected:
    input:
        defense="defense_systems/defensefinder/matrices/defensefinder_matrix.tsv",
        is_matrix="mobilome/insertion_sequences/matrices/is_family_matrix.tsv"

    output:
        "integration/stats/clade_corrected/defense_is_family_clade_corrected.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_is_family_clade_corrected.log"
    message: "Computing clade-corrected Defense × IS family associations."
    shell:
        "python integration/scripts/defense_vs_is_family_clade_corrected.py"


rule defense_vs_is_network:
    input:
        "integration/stats/is_association/defense_vs_is_family.tsv"

    output:
        "integration/plots/defense_vs_is_network.png"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_is_network.log"
    message: "Plotting Defense × IS association network."
    shell:
        "python integration/scripts/defense_vs_is_network.py"


############################################################
# DEFENSE × MOBILOME
############################################################

rule defense_vs_mobilome:
    input:
        defense="defense_systems/defensefinder/matrices/defensefinder_matrix.tsv",
        is_matrix="mobilome/insertion_sequences/matrices/is_family_matrix.tsv",
        prophages="mobilome/prophages/matrices/prophage_consensus_matrix.tsv"

    output:
        "integration/stats/clade_corrected/defense_vs_mobilome.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_mobilome.log"
    message: "Computing Defense × Mobilome associations."
    shell:
        "python integration/scripts/defense_vs_mobilome.py"


rule defense_vs_mobilome_scaled:
    input:
        "integration/stats/clade_corrected/defense_vs_mobilome.tsv"

    output:
        "integration/stats/clade_corrected/defense_vs_mobilome_scaled.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_mobilome_scaled.log"
    message: "Scaling Defense × Mobilome associations."
    shell:
        "python integration/scripts/defense_vs_mobilome_scaled.py"


rule defense_vs_mobilome_centrality:
    input:
        mobilome="integration/stats/clade_corrected/defense_vs_mobilome_scaled.tsv",
        centrality="defense_systems/defensefinder/stats/cooccurrence/defense_cooccurrence_centrality.tsv"

    output:
        "integration/stats/clade_corrected/defense_vs_mobilome_centrality.tsv"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/defense_vs_mobilome_centrality.log"
    message: "Computing centrality-weighted Defense × Mobilome associations."
    shell:
        "python integration/scripts/defense_vs_mobilome_centrality.py"


############################################################
# DEFENSE × PHYLOGENY
############################################################

rule mantel_defense_vs_phylogeny:
    input:
        defense="defense_systems/defensefinder/matrices/defensefinder_matrix.tsv"

    output:
        "integration/plots/mantel_defense_vs_phylogeny.png"

    threads: 1
    conda: "../envs/xtt_core.yml"
    log: "integration/logs/mantel_defense_vs_phylogeny.log"
    message: "Computing Mantel test between DefenseFinder profiles and phylogeny."
    shell:
        "python integration/scripts/mantel_defense_vs_phylogeny.py"
