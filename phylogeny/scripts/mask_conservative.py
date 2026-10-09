from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import Counter

# =========================
# RECOMBINATION THRESHOLD
# =========================

MIN_RECOMB_LENGTH = 1000
MIN_EVENT_COUNT = 2

# =========================
# INPUT FILES
# =========================

alignment_file = "cfml_xmfa.ML_sequence.fasta"
recomb_file = "cfml_xmfa.importation_status.txt"

output_alignment = "xmfa_conservative_masked.aln"

# =========================
# LOAD ALIGNMENT
# =========================

records = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(records[0].seq)

print(f"Alignment length: {alignment_length}")

# =========================
# LOAD RECOMBINATION REGIONS
# =========================

recomb_df = pd.read_csv(
    recomb_file,
    sep=r"\s+"
)

# =========================
# COUNT RECOMBINATION EVENTS
# =========================

recomb_counter = Counter()

for _, row in recomb_df.iterrows():

    start = int(row["Beg"])
    end = int(row["End"])

    length = end - start + 1

    # Skip short recombination events
    if length < MIN_RECOMB_LENGTH:
        continue

    # ClonalFrameML positions are usually 1-based
    for pos in range(start - 1, end):
        recomb_counter[pos] += 1

# =========================
# FILTER POSITIONS
# =========================

recomb_positions = {
    pos
    for pos, count in recomb_counter.items()
    if count >= MIN_EVENT_COUNT
}

print(f"Recombinant positions removed: {len(recomb_positions)}")

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
        seq[i]
        for i in keep_positions
    )

    record.seq = Seq(filtered_seq)

    new_records.append(record)

# =========================
# SAVE OUTPUT
# =========================

with open(output_alignment, "w") as out_handle:
    SeqIO.write(new_records, out_handle, "fasta")

print(f"Saved: {output_alignment}")
