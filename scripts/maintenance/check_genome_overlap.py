def report_genome_overlap(
    df1,
    df2,
    name1="df1",
    name2="df2"
):

    set1 = set(df1["Genome"])
    set2 = set(df2["Genome"])

    shared = set1 & set2

    only1 = set1 - set2

    only2 = set2 - set1

    print("\n===== Genome overlap report =====")

    print(f"{name1}: {len(set1)}")

    print(f"{name2}: {len(set2)}")

    print(f"Shared: {len(shared)}")

    print(f"\nOnly in {name1}: {len(only1)}")

    print(sorted(only1))

    print(f"\nOnly in {name2}: {len(only2)}")

    print(sorted(only2))