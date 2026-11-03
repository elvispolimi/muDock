#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="Input CSV file")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    criteria = ["sr_rmsd", "sr_score_diff"]
    algorithms = df["algorithm"].unique()

    fig, ax = plt.subplots(figsize=(10, 6))

    markers = {
        "ga": "o",
        "lga": "s",
    }

    colors = {
        ("ga", "sr_rmsd"): "#90CAF9",        # Light blue
        ("lga", "sr_rmsd"): "#1565C0",       # Dark blue
        ("ga", "sr_score_diff"): "#FFCC80",  # Light orange
        ("lga", "sr_score_diff"): "#EF6C00", # Dark orange
    }

    linestyles = {
        "sr_rmsd": "-",
        "sr_score_diff": "--",
    }

    for algorithm in algorithms:
        for criterion in criteria:
            subset = df[df["algorithm"] == algorithm].sort_values("budget")

            ax.plot(
                subset["budget"],
                subset[criterion],
                color=colors[(algorithm, criterion)],
                marker=markers.get(algorithm, "o"),
                linestyle=linestyles[criterion],
                linewidth=2,
                markersize=7,
                label=f"{algorithm.upper()} - {criterion}",
            )

    ax.set_xlabel("Budget", fontsize=12)
    ax.set_ylabel("Success rate (%)", fontsize=12)
    ax.set_title("Success Rate vs Budget", fontsize=14)

    ax.set_xticks(sorted(df["budget"].unique()))
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.legend()
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig("success_rate_vs_budget.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()