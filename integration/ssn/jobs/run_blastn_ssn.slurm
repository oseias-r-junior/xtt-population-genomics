#!/bin/bash
#SBATCH --job-name=ssn_blastn
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

#SBATCH --output=integration/ssn/logs/%x_%j.out
#SBATCH --error=integration/ssn/logs/%x_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate xtt_blast

echo "=========================================="
echo "Environment"
echo "blastn : $(which blastn)"
blastn -version
echo "Threads: $SLURM_CPUS_PER_TASK"
echo "Started: $(date)"
echo "=========================================="

SSN=${PROJECT_DIR:-$PWD}/integration/ssn

cd "$SSN"

# makeblastdb already completed successfully in an earlier interactive run
# (db/mge_db.* present) -- skip rebuilding it unless missing.
if [[ ! -f "blast/db/mge_db.nsq" ]]; then
    echo "Rebuilding BLAST db (not found)..."
    makeblastdb -in fasta/all_mge_sequences.fna -dbtype nucl -out blast/db/mge_db
fi

blastn \
    -query fasta/all_mge_sequences.fna \
    -db blast/db/mge_db \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp" \
    -evalue 1e-10 \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    -out blast/all_vs_all_raw.tsv

echo
echo "Finished: $(date)"
wc -l blast/all_vs_all_raw.tsv
