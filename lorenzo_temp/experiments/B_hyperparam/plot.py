#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt


CSV_FILE = "./lorenzo_temp/experiments/B_hyperparam/small/res_medians_score.csv"
OUTPUT_DIR = "plots_medians_score"

os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(CSV_FILE)

for (rho, epsilon), group in df.groupby(["rho", "epsilon"]):

    # Sort by rot ascending
    group = group.sort_values("rot")

    plt.figure(figsize=(12, 6))

    plt.bar(
        group["name"].astype(str),
        group["median_score"]
    )

    plt.xlabel("Name (ordered by rot)")
    plt.ylabel("Median score")
    plt.title(f"rho={rho}, epsilon={epsilon}")

    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Make a filesystem-safe filename
    filename = f"rho_{rho}_epsilon_{epsilon}.png"
    output_path = os.path.join(OUTPUT_DIR, filename)

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_path}")