#!/usr/bin/env python3

import os

import pandas as pd
import matplotlib.pyplot as plt


CSV_FILE = "./lorenzo_temp/experiments/B_hyperparam/small/res_medians_score.csv"
OUTPUT_DIR = "plots_medians_score"

os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(CSV_FILE)

# Create a configuration label
df["config"] = df.apply(
    lambda row: f"rho={row['rho']}, epsilon={row['epsilon']}",
    axis=1
)

# Determine the x-axis order from rot
# Assumes rot is the same for every configuration
order = (
    df[["name", "rot"]]
    .drop_duplicates()
    .sort_values("rot")
)

names = order["name"].astype(str).tolist()

plt.figure(figsize=(14, 7))

# Plot one line for each configuration
for config, group in df.groupby("config"):

    # Sort according to rot
    group = (
        group.set_index("name")
        .reindex(names)
        .reset_index()
    )

    plt.plot(
        group["name"],
        group["median_score"],
        marker="o",
        linewidth=2,
        label=config
    )

plt.xlabel("Name (ordered by rot)")
plt.ylabel("Median score")
# plt.yscale("log")
plt.title("Median score by hyperparameter configuration")

plt.xticks(rotation=45, ha="right")

plt.grid(axis="y", alpha=0.3)
plt.legend()

plt.tight_layout()

output_path = os.path.join(OUTPUT_DIR, "all_configs.png")

plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {output_path}")