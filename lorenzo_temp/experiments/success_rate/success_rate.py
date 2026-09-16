#!/usr/bin/env python3
"""
Compare docking success rates of GA vs LGA runs, per ligand.

A run is "successful" when score < 0.
GA runs are identified by missing (N/A) `rate` and `iterations` fields;
every other run is an LGA run, aggregated over all iterations/rate configs.

Usage:
    python success_rate.py results.csv [-o success_rate.png] [--csv summary.csv]
"""

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NA_VALUES = ["N/A", "NA", "n/a", "na", "", "-", "None", "null"]


def load(path):
    df = pd.read_csv(path, na_values=NA_VALUES, keep_default_na=True)

    required = {"name", "score", "rate", "iterations"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"Missing required column(s): {', '.join(sorted(missing))}")

    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["rate"] = pd.to_numeric(df["rate"], errors="coerce")
    df["iterations"] = pd.to_numeric(df["iterations"], errors="coerce")

    bad = df["score"].isna().sum()
    if bad:
        print(f"Warning: dropping {bad} row(s) with a non-numeric score.")
        df = df.dropna(subset=["score"])

    # GA: no local search parameters. LGA: everything else.
    df["algo"] = np.where(
        df["rate"].isna() & df["iterations"].isna(), "GA", "LGA"
    )
    return df


def summarize(df):
    g = df.groupby(["name", "algo"])["score"]
    stats = g.agg(total="size", successes=lambda s: (s < 0).sum())
    stats["success_rate"] = stats["successes"] / stats["total"]

    table = stats["success_rate"].unstack("algo")
    counts = stats["total"].unstack("algo")

    for algo in ("GA", "LGA"):
        if algo not in table.columns:
            table[algo] = np.nan
            counts[algo] = 0

    table = table[["GA", "LGA"]].sort_index()
    counts = counts[["GA", "LGA"]].reindex(table.index).fillna(0).astype(int)
    return table, counts


def plot(table, counts, out_path, annotate=True):
    ligands = list(table.index)
    x = np.arange(len(ligands))
    width = 0.38

    fig_w = max(8.0, 0.55 * len(ligands) + 3.0)
    fig, ax = plt.subplots(figsize=(fig_w, 6.0))

    ga = table["GA"].fillna(0.0).to_numpy() * 100
    lga = table["LGA"].fillna(0.0).to_numpy() * 100

    b1 = ax.bar(x - width / 2, ga, width, label="GA", color="#4C72B0")
    b2 = ax.bar(x + width / 2, lga, width, label="LGA", color="#DD8452")

    if annotate:
        for bars in (b1, b2):
            ax.bar_label(bars, fmt="%.0f", padding=2, fontsize=7)

    ax.set_xlabel("Ligand")
    ax.set_ylabel("Success rate (% of runs with score < 0)")
    ax.set_title("Docking success rate per ligand: GA vs LGA")
    ax.set_xticks(x)
    ax.set_xticklabels(ligands, rotation=45, ha="right")
    ax.set_ylim(0, 105)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    ax.legend(title="Algorithm")

    n_ga = counts["GA"].sum()
    n_lga = counts["LGA"].sum()
    fig.text(
        0.01, 0.01,
        f"{n_ga} GA runs, {n_lga} LGA runs (all iterations/rate configurations pooled)",
        fontsize=8, color="gray",
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="input CSV file")
    p.add_argument("-o", "--output", default="success_rate.png",
                   help="output image (default: success_rate.png)")
    p.add_argument("--csv-out", dest="csv_out",
                   help="also write the summary table to this CSV")
    p.add_argument("--no-labels", action="store_true",
                   help="don't annotate bars with their value")
    args = p.parse_args()

    df = load(args.csv)
    table, counts = summarize(df)

    report = pd.DataFrame({
        "GA_runs": counts["GA"],
        "GA_success_rate": table["GA"],
        "LGA_runs": counts["LGA"],
        "LGA_success_rate": table["LGA"],
    })
    print(report.to_string(float_format=lambda v: f"{v:.3f}"))

    if args.csv_out:
        report.to_csv(args.csv_out)
        print(f"Wrote {args.csv_out}")

    plot(table, counts, args.output, annotate=not args.no_labels)


if __name__ == "__main__":
    main()