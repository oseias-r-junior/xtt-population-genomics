import re

def normalize_genome_name(name):

    """
    Normalize genome/sample names across all analyses.
    """

    name = str(name)

    # remove common suffixes/extensions
    suffixes = [
        "_out",
        ".fasta",
        ".fa",
        ".fna",
        ".faa",
        ".gff"
    ]

    for s in suffixes:

        if name.endswith(s):

            name = name[:-len(s)]

    # remove whitespace
    name = name.strip()

    # harmonize common strain separators

    if re.match(r"^[A-Za-z]{2}\d+_\d+$", name):
        name = name.replace("_", "-")

    # ======================================
    # KNOWN GENOME ALIASES
    # ======================================

    aliases = {

        "MO22003_2": "MO22003-2",
        "L_G_5": "LG_5",
        "L_G_2": "LG_2",
        "Km8": "XtKm8",
        "Km9": "XtKm9",

    }

    if name in aliases:

        name = aliases[name]

    # replace spaces with underscore
    name = name.replace(" ", "_")

    # collapse repeated underscores
    name = re.sub(r"_+", "_", name)

    return name