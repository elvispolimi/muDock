import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


def main():

    # 1. Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Generate heatmaps for LGA vs GA docking scores."
    )

    parser.add_argument(
        "csv_file",
        type=str,
        help="Path to the results CSV file."
    )

    parser.add_argument(
        "--agg",
        choices=["mean", "median"],
        default="median",
        help="Aggregation method to combine seeds (default: median)."
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default=".",
        help="Directory to save the generated heatmaps (default: current directory)."
    )

    args = parser.parse_args()

    # 2. Read the CSV file
    # pandas automatically converts "N/A" strings to NaN by default
    df = pd.read_csv(args.csv_file)

    # Ensure rate and iteration are numeric
    # NaN for standard GA runs
    df["rate"] = pd.to_numeric(df["rate"], errors="coerce")
    df["iterations"] = pd.to_numeric(df["iterations"], errors="coerce")

    # 3. Separate GA runs (baseline) and LGA runs
    ga_runs = df[
        df["rate"].isna() &
        df["iterations"].isna()
    ]

    lga_runs = df[
        df["rate"].notna() &
        df["iterations"].notna()
    ]

    if ga_runs.empty or lga_runs.empty:
        print(
            "Error: Could not find both GA (N/A rate/iterations) "
            "and LGA runs in the dataset."
        )
        return

    # ------------------------------------------------------------------
    # 4. Remove invalid runs (score > 0), but keep track of them
    # ------------------------------------------------------------------

    # Positive scores are considered invalid/unreliable docking results.
    ga_positive = ga_runs["score"] > 0
    lga_positive = lga_runs["score"] > 0

    # Keep statistics BEFORE filtering
    ga_stats = (
        ga_runs
        .groupby("name")["score"]
        .agg(
            total_runs="count",
            positive_runs=lambda x: (x > 0).sum()
        )
        .reset_index()
    )

    ga_stats["valid_runs"] = (
        ga_stats["total_runs"] - ga_stats["positive_runs"]
    )

    lga_stats = (
        lga_runs
        .groupby(["name", "rate", "iterations"])["score"]
        .agg(
            total_runs="count",
            positive_runs=lambda x: (x > 0).sum()
        )
        .reset_index()
    )

    lga_stats["valid_runs"] = (
        lga_stats["total_runs"] - lga_stats["positive_runs"]
    )

    # Actually remove positive scores from the data used for aggregation
    ga_runs = ga_runs[~ga_positive]
    lga_runs = lga_runs[~lga_positive]

    # ------------------------------------------------------------------
    # 5. Aggregate valid data
    # ------------------------------------------------------------------

    agg_func = args.agg

    # Baseline GA scores per ligand
    ga_baseline = (
        ga_runs
        .groupby("name")["score"]
        .agg(agg_func)
        .reset_index()
    )

    ga_baseline.rename(
        columns={"score": "ga_score"},
        inplace=True
    )

    # LGA scores per ligand, rate, and iterations
    lga_agg = (
        lga_runs
        .groupby(["name", "rate", "iterations"])["score"]
        .agg(agg_func)
        .reset_index()
    )

    # ------------------------------------------------------------------
    # 6. Merge score data with reliability statistics
    # ------------------------------------------------------------------

    merged_data = pd.merge(
        lga_agg,
        ga_baseline,
        on="name"
    )

    # Add LGA reliability information
    merged_data = pd.merge(
        merged_data,
        lga_stats,
        on=["name", "rate", "iterations"],
        how="left"
    )

    # Add GA reliability information
    merged_data = pd.merge(
        merged_data,
        ga_stats[[
            "name",
            "total_runs",
            "positive_runs",
            "valid_runs"
        ]].rename(columns={
            "total_runs": "ga_total_runs",
            "positive_runs": "ga_positive_runs",
            "valid_runs": "ga_valid_runs"
        }),
        on="name",
        how="left"
    )

    # Difference: LGA - GA
    # Negative means LGA scored lower/better than GA
    merged_data["score_diff"] = (
        merged_data["score"] - merged_data["ga_score"]
    )

    # ------------------------------------------------------------------
    # 7. Create output directory
    # ------------------------------------------------------------------

    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # 8. Generate a heatmap for each ligand
    # ------------------------------------------------------------------

    ligands = merged_data["name"].unique()

    for ligand in ligands:

        # Filter data for the specific ligand
        ligand_data = merged_data[
            merged_data["name"] == ligand
        ]

        # Pivot the data to create a 2D matrix for the heatmap
        # Columns = rate (X-axis)
        # Index = iterations (Y-axis)
        pivot_table = ligand_data.pivot(
            index="iterations",
            columns="rate",
            values="score_diff"
        )

        pivot_table = pivot_table.sort_index(ascending=False)

        # --------------------------------------------------------------
        # Create annotations containing:
        #
        #   score difference
        #   valid runs / total runs
        #
        # Example:
        #
        #   -1.23
        #   9/10
        #
        # --------------------------------------------------------------

        annotation_table = ligand_data.copy()

        annotation_table["annotation"] = annotation_table.apply(
            lambda row: (
                f"{row['score_diff']:.2f}\n"
                f"{int(row['valid_runs'])}/{int(row['total_runs'])}"
            ),
            axis=1
        )

        annotation_table = annotation_table.pivot(
            index="iterations",
            columns="rate",
            values="annotation"
        )

        annotation_table = annotation_table.reindex(
            index=pivot_table.index,
            columns=pivot_table.columns
        )

        # --------------------------------------------------------------
        # Plotting
        # --------------------------------------------------------------

        plt.figure(figsize=(9, 7))

        sns.heatmap(
            pivot_table,
            annot=annotation_table,
            fmt="",
            cmap="coolwarm",
            center=0,
            cbar_kws={
                "label": f"Score Difference ({args.agg.capitalize()})"
            }
        )

        plt.title(
            f"{ligand} - LGA vs GA\n"
            f"(Score Difference: LGA - GA)\n"
            f"Annotation: valid runs / total runs"
        )

        plt.xlabel("Local Search Rate")
        plt.ylabel("Local Search Iterations")

        plt.tight_layout()

        # Save the figure
        out_path = os.path.join(
            args.outdir,
            f"heatmap_{ligand}_{args.agg}.png"
        )

        plt.savefig(out_path, dpi=300)
        plt.close()

        print(f"Saved heatmap for {ligand} to {out_path}")


if __name__ == "__main__":
    main()