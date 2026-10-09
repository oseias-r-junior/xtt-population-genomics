from Bio import SeqIO
from Bio.Seq import Seq

import pandas as pd
from collections import Counter
import numpy as np

# ==========================================
# INPUT FILES
# ==========================================

alignment_file = "core_gene_alignment.aln"

recomb_file = "cfml_output.importation_status.txt"

# ==========================================
# LOAD ALIGNMENT
# ==========================================

records = list(
    SeqIO.parse(alignment_file, "fasta")
)

alignment_length = len(records[0].seq)

print(f"\nAlignment length: {alignment_length}")

# ==========================================
# LOAD RECOMBINATION EVENTS
# ==========================================

recomb_df = pd.read_csv(
    recomb_file,
    sep=r"\s+"
)

# ==========================================
# EVENT LENGTHS
# ==========================================

recomb_df["Length"] = (
    recomb_df["End"] -
    recomb_df["Beg"] +
    1
)

# ==========================================
# DEFINE THRESHOLDS
# ==========================================

thresholds = {

    "mild": {
        "length": int(
            np.percentile(
                recomb_df["Length"],
                75
            )
        ),
        "recurrence": 3
    },

    "moderate": {
        "length": int(
            np.percentile(
                recomb_df["Length"],
                90
            )
        ),
        "recurrence": 5
    },

    "stringent": {
        "length": int(
            np.percentile(
                recomb_df["Length"],
                95
            )
        ),
        "recurrence": 6
    }
}

print("\nThresholds:")

for regime, values in thresholds.items():

    print(
        regime,
        values
    )

# ==========================================
# PROCESS EACH REGIME
# ==========================================

for regime, params in thresholds.items():

    print("\n======================================")
    print(f"Processing: {regime}")
    print("======================================")

    min_length = params["length"]

    min_recurrence = params["recurrence"]

    recomb_counter = Counter()

    # ======================================
    # COUNT EVENTS
    # ======================================

    for _, row in recomb_df.iterrows():

        start = int(row["Beg"])
        end = int(row["End"])

        length = end - start + 1

        # Skip short events
        if length < min_length:
            continue

        for pos in range(start - 1, end):

            recomb_counter[pos] += 1

    # ======================================
    # FILTER POSITIONS
    # ======================================

    recomb_positions = {

        pos
        for pos, count
        in recomb_counter.items()
        if count >= min_recurrence
    }

    print(
        f"Recombinant positions removed: "
        f"{len(recomb_positions)}"
    )

    # ======================================
    # KEEP POSITIONS
    # ======================================

    keep_positions = [

        i
        for i in range(alignment_length)
        if i not in recomb_positions
    ]

    print(
        f"Nonrecombinant positions kept: "
        f"{len(keep_positions)}"
    )

    # ======================================
    # BUILD NEW ALIGNMENT
    # ======================================

    new_records = []

    for record in records:

        seq = str(record.seq)

        filtered_seq = "".join(

            seq[i]
            for i in keep_positions
        )

        record_copy = record[:]

        record_copy.seq = Seq(
            filtered_seq
        )

        new_records.append(
            record_copy
        )

    # ======================================
    # SAVE OUTPUT
    # ======================================

    output_alignment = (
        f"core_masked_{regime}.aln"
    )

    with open(
        output_alignment,
        "w"
    ) as out_handle:

        SeqIO.write(
            new_records,
            out_handle,
            "fasta"
        )

    print(
        f"Saved: {output_alignment}"
    )
