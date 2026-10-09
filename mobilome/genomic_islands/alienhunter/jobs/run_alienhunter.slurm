#!/bin/bash
#SBATCH --job-name=alienhunter
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00

#SBATCH --output=xtt_project/mobilome/genomic_islands/alienhunter/logs/%x_%j.out
#SBATCH --error=xtt_project/mobilome/genomic_islands/alienhunter/logs/%x_%j.err

set -euo pipefail

############################################################
# Project
############################################################

PROJECT="${PROJECT_DIR:-$PWD}"

GENOME=${GENOME:-}

FNA_DIR="$PROJECT/xtt_project/annotations/fna"

OUTDIR="$PROJECT/xtt_project/mobilome/genomic_islands/alienhunter/raw_outputs"

SCRIPT="$PROJECT/xtt_project/mobilome/genomic_islands/alienhunter/scripts/run_alienhunter.py"

############################################################
# Environment
############################################################

eval "$(micromamba shell hook --shell bash)"

############################################################
# Run
############################################################

echo "========================================="
echo "AlienHunter batch execution"
echo "Started: $(date)"
echo "========================================="

mkdir -p "${OUTDIR}"

N=0

# If GENOME is set, run only that genome
if [[ -n "${GENOME}" ]]; then
    GENOMES="${FNA_DIR}/${GENOME}.fna"
else
    GENOMES="${FNA_DIR}"/*.fna
fi

for fna in ${GENOMES}
do
    genome_name=$(basename "${fna}")
    # Remove a extensão .fna para pegar apenas o nome base (ex: AB22007)
    genome_base="${genome_name%.fna}" 

    # ---------------------------------------------------------
    # PONTO CHAVE: Defina como o arquivo de saída se chama.
    # Substitua ".sco" pela extensão real gerada pelo seu script Python.
    # ---------------------------------------------------------
    EXPECTED_OUTPUT="${OUTDIR}/${genome_base}.sco"

    # Se o arquivo de saída já existir (e não estiver vazio), pule este genoma.
    if [[ -s "${EXPECTED_OUTPUT}" ]]; then
        echo "Ignorando ${genome_name}: arquivo de saída já existe."
        continue
    fi

    # Se chegou aqui, o genoma ainda não foi processado
    ((N+=1))

    echo
    echo "[$N] Processing: ${genome_name}"

    micromamba run -n xtt_alienhunter \
        python "${SCRIPT}" \
        --genome "${fna}" \
        --outdir "${OUTDIR}"
done

echo
echo "========================================="
echo "Finished: $(date)"
echo "Processed genomes: ${N}"
echo "========================================="