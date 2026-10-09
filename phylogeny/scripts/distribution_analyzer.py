from collections import Counter
import pandas as pd
import numpy as np

# ==========================================
# INPUT FILE
# ==========================================

recomb_file = "cfml_output.importation_status.txt"

# ==========================================
# LOAD RECOMBINATION EVENTS
# ==========================================

recomb_df = pd.read_csv(
    recomb_file,
    sep=r"\s+"
)

print("\nLoaded recombination events:")
print(recomb_df.head())

# ==========================================
# EVENT LENGTHS
# ==========================================

recomb_df["Length"] = (
    recomb_df["End"] -
    recomb_df["Beg"] +
    1
)

lengths = recomb_df["Length"]

print("\n==========================================")
print("RECOMBINATION EVENT LENGTHS")
print("==========================================")

print(f"Total events: {len(lengths)}")

print(f"Mean length: {lengths.mean():.2f}")
print(f"Median length: {lengths.median():.2f}")

print(f"Min length: {lengths.min()}")
print(f"Max length: {lengths.max()}")

# ==========================================
# LENGTH PERCENTILES
# ==========================================

length_percentiles = [50, 75, 90, 95, 99]

print("\nLength percentiles:")

for p in length_percentiles:

    value = np.percentile(lengths, p)

    print(f"{p}th percentile: {value:.2f}")

# ==========================================
# POSITION RECURRENCE
# ==========================================

print("\n==========================================")
print("POSITION RECURRENCE")
print("==========================================")

recomb_counter = Counter()

for _, row in recomb_df.iterrows():

    start = int(row["Beg"])
    end = int(row["End"])

    for pos in range(start - 1, end):

        recomb_counter[pos] += 1

# Convert recurrence values to dataframe

recurrence_values = list(
    recomb_counter.values()
)

recurrence_series = pd.Series(
    recurrence_values
)

print(f"Positions affected: {len(recurrence_values)}")

print(f"Mean recurrence: {recurrence_series.mean():.2f}")
print(f"Median recurrence: {recurrence_series.median():.2f}")

print(f"Min recurrence: {recurrence_series.min()}")
print(f"Max recurrence: {recurrence_series.max()}")

# ==========================================
# RECURRENCE PERCENTILES
# ==========================================

recurrence_percentiles = [50, 75, 90, 95, 99]

print("\nRecurrence percentiles:")

for p in recurrence_percentiles:

    value = np.percentile(
        recurrence_values,
        p
    )

    print(f"{p}th percentile: {value:.2f}")

# ==========================================
# SUGGESTED THRESHOLDS
# ==========================================

print("\n==========================================")
print("SUGGESTED DATA-DRIVEN THRESHOLDS")
print("==========================================")

mild_length = int(
    np.percentile(lengths, 75)
)

moderate_length = int(
    np.percentile(lengths, 90)
)

stringent_length = int(
    np.percentile(lengths, 95)
)

mild_recurrence = int(
    np.percentile(recurrence_values, 75)
)

moderate_recurrence = int(
    np.percentile(recurrence_values, 90)
)

stringent_recurrence = int(
    np.percentile(recurrence_values, 95)
)

thresholds = pd.DataFrame({

    "Regime": [
        "Mild",
        "Moderate",
        "Stringent"
    ],

    "Length cutoff": [
        mild_length,
        moderate_length,
        stringent_length
    ],

    "Recurrence cutoff": [
        mild_recurrence,
        moderate_recurrence,
        stringent_recurrence
    ]
})

print(thresholds)

# ==========================================
# SAVE TABLE
# ==========================================

thresholds.to_excel(
    "data_driven_thresholds.xlsx",
    index=False
)

print("\nSaved:")
print("data_driven_thresholds.xlsx")
