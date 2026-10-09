from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import Counter

# =========================
# INPUT FILES
# =========================

alignment_file = "cfml_xmfa.ML_sequence.fasta"
recomb_file = "cfml_xmfa.importation_status.txt"

output_alignment = "xmfa_full_masked.aln"

# =========================
# LOAD ALIGNMENT
# =========================

records = list(SeqIO.parse(alignment_file, "fasta"))

seq_names = [r.id for r in records]
seqs = [list(str(r.seq)) for r in records]

alignment_length = len(seqs[0])

print(f"Alignment length: {alignment_length}")

# =========================
# LOAD RECOMBINATION REGIONS
# =========================

recomb_df = pd.read_csv(
    recomb_file,
    sep=r"\s+"
)

# Expected:
# column 0 = branch/node
# column 1 = start
# column 2 = end

recomb_positions = set()

for _, row in recomb_df.iterrows():

    start = int(row["Beg"])
    end = int(row["End"])

    for pos in range(start - 1, end):
        recomb_positions.add(pos)

print(f"Recombinant positions: {len(recomb_positions)}")

# =========================
# KEEP ONLY NONRECOMBINANT SITES
# =========================

keep_positions = [
    i for i in range(alignment_length)
    if i not in recomb_positions
]

print(f"Nonrecombinant positions kept: {len(keep_positions)}")

# =========================
# BUILD NEW ALIGNMENT
# =========================

new_records = []

for record in records:

    seq = str(record.seq)

    filtered_seq = "".join(
        seq[i] for i in keep_positions
    )

    record.seq = Seq(filtered_seq)

    new_records.append(record)

# =========================
# SAVE OUTPUT
# =========================

with open(output_alignment, "w") as out_handle:
    SeqIO.write(new_records, out_handle, "fasta")

print(f"Saved: {output_alignment}")



