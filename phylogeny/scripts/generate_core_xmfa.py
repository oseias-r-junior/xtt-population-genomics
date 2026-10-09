import pathlib
ROOT = str(pathlib.Path(__file__).resolve().parents[2]) + "/"
from pathlib import Path

# =========================
# SETTINGS
# =========================

input_dir = ROOT + "input_data_hrr_new_phylo_tree/pan_genome_sequences"

output_xmfa = "core_genes_only.xmfa"

expected_sequences = 113

# =========================
# FIND ALIGNMENTS
# =========================

alignment_files = sorted(
    Path(input_dir).glob("*.fa.aln")
)

kept = 0
removed = 0

# =========================
# BUILD XMFA
# =========================

with open(output_xmfa, "w") as outfile:

    for aln_file in alignment_files:

        # Count sequences
        nseq = 0

        with open(aln_file) as f:

            for line in f:
                if line.startswith(">"):
                    nseq += 1

        # Keep only core genes
        if nseq == expected_sequences:

            with open(aln_file) as f:
                outfile.write(f.read())

            outfile.write("\n=\n")

            kept += 1

        else:
            removed += 1

# =========================
# REPORT
# =========================

print(f"Core gene alignments kept: {kept}")
print(f"Non-core/incomplete alignments removed: {removed}")

print(f"XMFA saved: {output_xmfa}")
