#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare GA and LGA RMSD across seeds and plot the "
            "mean/median RMSD difference for each ligand."
        )
    )

    parser.add_argument(
        "csv",
        help="Path to the results CSV"
    )

    parser.add_argument(
        "--aggregation",
        choices=["mean", "median"],
        default="median",
        help="Aggregation across seeds (default: median)"
    )

    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Optional prefix for output image files"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------

    df = pd.read_csv(args.csv)

    # Remove possible whitespace from column names
    df.columns = df.columns.str.strip()

    # ------------------------------------------------------------------
    # Identify GA and LGA runs
    #
    # GA: rate and iterations are N/A
    # LGA: rate and iterations are present
    # ------------------------------------------------------------------
    # df['rate'] = pd.to_numeric(df['rate'], errors='coerce')
    # df['iterations'] = pd.to_numeric(df['iterations'], errors='coerce')
    
    ga = df[
        df["rate"].isna() &
        df["iterations"].isna()
    ].copy()

    lga = df[
        df["rate"].notna() &
        df["iterations"].notna()
    ].copy()

    print(f"GA runs:  {len(ga)}")
    print(f"LGA runs: {len(lga)}")



    # ------------------------------------------------------------------
    # Keep only the columns needed for the comparison
    # ------------------------------------------------------------------

    rmsd_columns = [
        "rmsd_best_scoring_pose",
        "rmsd_min",
    ]

    ga = ga[
        ["name", "seed"] + rmsd_columns
    ].rename(
        columns={
            column: f"{column}_ga"
            for column in rmsd_columns
        }
    )

    lga = lga[
        ["name", "seed"] + rmsd_columns
    ].rename(
        columns={
            column: f"{column}_lga"
            for column in rmsd_columns
        }
    )

    # ------------------------------------------------------------------
    # Match GA and LGA using ligand + seed
    # ------------------------------------------------------------------

    merged = pd.merge(
        ga,
        lga,
        on=["name", "seed"],
        how="inner"
    )

    print(f"Matched GA/LGA seed pairs: {len(merged)}")

    if merged.empty:
        raise RuntimeError(
            "No matching GA/LGA pairs were found. "
            "Check the rate/iterations columns and CSV format."
        )

    # ------------------------------------------------------------------
    # Compute differences:
    #
    # positive = LGA has lower RMSD
    # negative = LGA has higher RMSD
    # ------------------------------------------------------------------

    merged["difference_best_scoring_pose"] = (
        merged["rmsd_best_scoring_pose_ga"]
        - merged["rmsd_best_scoring_pose_lga"]
    )

    merged["difference_min"] = (
        merged["rmsd_min_ga"]
        - merged["rmsd_min_lga"]
    )

    # ------------------------------------------------------------------
    # Aggregate across seeds
    # ------------------------------------------------------------------

    aggregation_function = getattr(
        merged.groupby("name"),
        args.aggregation
    )

    aggregated = aggregation_function()[
        [
            "difference_best_scoring_pose",
            "difference_min",
        ]
    ]

    print("\nAggregated differences:")
    print(aggregated)

    # ------------------------------------------------------------------
    # Plot function
    # ------------------------------------------------------------------

    def plot_difference(column, ylabel, title, filename):
        values = aggregated[column]

        plt.figure(figsize=(12, 6))

        plt.bar(
            values.index,
            values.values
        )

        # Zero line makes interpretation easier
        plt.axhline(
            0,
            linewidth=1
        )

        plt.xlabel("Ligand")
        plt.ylabel(ylabel)
        plt.title(title)

        plt.xticks(
            rotation=45,
            ha="right"
        )

        plt.tight_layout()

        if filename:
            plt.savefig(
                filename,
                dpi=300,
                bbox_inches="tight"
            )

        plt.show()

    # ------------------------------------------------------------------
    # Plot 1: RMSD of best scoring pose
    # ------------------------------------------------------------------

    prefix = args.output_prefix

    output_best = (
        f"{prefix}_best_scoring_pose.png"
        if prefix
        else None
    )

    plot_difference(
        "difference_best_scoring_pose",
        "GA RMSD - LGA RMSD (Å)",
        f"{args.aggregation.capitalize()} improvement in RMSD of best scoring pose",
        output_best
    )

    # ------------------------------------------------------------------
    # Plot 2: minimum RMSD
    # ------------------------------------------------------------------

    output_min = (
        f"{prefix}_rmsd_min.png"
        if prefix
        else None
    )

    plot_difference(
        "difference_min",
        "GA RMSD - LGA RMSD (Å)",
        f"{args.aggregation.capitalize()} improvement in minimum RMSD",
        output_min
    )


if __name__ == "__main__":
    main()