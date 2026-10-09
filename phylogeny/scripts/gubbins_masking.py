from Bio import SeqIO
from Bio.Seq import Seq
import re

# ==========================================
# INPUT FILES
# ==========================================

alignment_file = "core_gene_alignment.aln"

gff_file = "gubbins_core.recombination_predictions.gff"

output_alignment = "core_gubbins_masked.aln"

# ==========================================
# LOAD ALIGNMENT
# ==========================================

records = list(
    SeqIO.parse(alignment_file, "fasta")
)

alignment_length = len(records[0].seq)

print(f"Alignment length: {alignment_length}")

# ==========================================
# LOAD RECOMBINATION REGIONS
# ==========================================

recomb_positions = set()

with open(gff_file) as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")

        start = int(fields[3])
        end = int(fields[4])

        for pos in range(start - 1, end):
            recomb_positions.add(pos)

print(f"Recombinant positions: {len(recomb_positions)}")

# ==========================================
# KEEP NONRECOMBINANT SITES
# ==========================================

keep_positions = [

    i for i in range(alignment_length)

    if i not in recomb_positions
]

print(f"Positions retained: {len(keep_positions)}")

# ==========================================
# BUILD FILTERED ALIGNMENT
# ==========================================

new_records = []

for record in records:

    seq = str(record.seq)

    filtered_seq = "".join(

        seq[i]

        for i in keep_positions
    )

    record.seq = Seq(filtered_seq)

    new_records.append(record)

# ==========================================
# SAVE OUTPUT
# ==========================================

with open(output_alignment, "w") as out_handle:

    SeqIO.write(
        new_records,
        out_handle,
        "fasta"
    )

print(f"Saved: {output_alignment}")
