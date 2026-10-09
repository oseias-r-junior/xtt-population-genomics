############################################################
# PLASMIDS
############################################################

############################################################
# TARGETS
############################################################

PLASMID_TARGETS = (
    "mobilome/plasmids/matrices/plasmid_consensus.tsv",
    "mobilome/plasmids/stats/plasmid_consensus_genome_summary.tsv",
)

############################################################
# MOB-SUITE
############################################################

rule parse_mobsuite_plasmids:
    input:
        "mobilome/plasmids/mob_suite/raw_outputs"

    output:
        "mobilome/plasmids/mob_suite/matrices/mobsuite_plasmids.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/plasmids/logs/parse_mobsuite_plasmids.log"

    message:
        "Parsing MOB-suite plasmid contigs and typing."

    shell:
        "python mobilome/plasmids/mob_suite/scripts/parse_mobsuite_plasmids.py"

############################################################
# GENOMAD (shared run with prophages)
############################################################

rule parse_genomad_plasmids:
    input:
        "mobilome/prophages/genomad/raw_outputs"

    output:
        "mobilome/plasmids/genomad/matrices/genomad_plasmids.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/plasmids/logs/parse_genomad_plasmids.log"

    message:
        "Parsing geNomad plasmid calls (shared geNomad run with prophages)."

    shell:
        "python mobilome/plasmids/genomad/scripts/parse_genomad_plasmids.py"

############################################################
# CONSENSUS
############################################################

rule build_plasmid_consensus:
    input:
        mobsuite="mobilome/plasmids/mob_suite/matrices/mobsuite_plasmids.tsv",
        genomad="mobilome/plasmids/genomad/matrices/genomad_plasmids.tsv"

    output:
        "mobilome/plasmids/matrices/plasmid_consensus.tsv",
        "mobilome/plasmids/stats/plasmid_consensus_genome_summary.tsv"

    threads:
        1

    conda:
        "../envs/xtt_core.yml"

    log:
        "mobilome/plasmids/logs/build_plasmid_consensus.log"

    message:
        "Building 2-way (MOB-suite x geNomad) plasmid consensus with mobility typing."

    shell:
        "python mobilome/plasmids/scripts/build_plasmid_consensus.py"
