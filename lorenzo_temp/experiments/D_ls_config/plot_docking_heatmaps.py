import argparse
import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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
    df = pd.read_csv(args.csv_file)

    # Ensure numeric columns
    for col in ["rate", "iter", "score", "rot"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "rot" not in df.columns:
        parser.error(
            "The CSV must contain a 'rot' column to create rotamer-bin heatmaps."
        )

    # 3. Separate GA runs (baseline) and LGA runs
    ga_runs = df[df["rate"].isna() & df["iter"].isna()]
    lga_runs = df[df["rate"].notna() & df["iter"].notna()]

    if ga_runs.empty or lga_runs.empty:
        print(
            "Error: Could not find both GA (N/A rate/iter) "
            "and LGA runs in the dataset."
        )
        return

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    # 4. Aggregate seeds using the selected metric
    agg_func = args.agg

    # Baseline GA scores per ligand
    ga_baseline = (
        ga_runs.groupby("name")["score"]
        .agg(agg_func)
        .reset_index()
        .rename(columns={"score": "ga_score"})
    )

    # LGA scores per ligand and configuration.
    # Keep rot so each complex can be assigned to a bin.
    lga_agg = (
        lga_runs.groupby(["name", "rate", "iter", "rot"])["score"]
        .agg(agg_func)
        .reset_index()
    )

    # 5. Merge and calculate score differences
    merged_data = pd.merge(
        lga_agg,
        ga_baseline,
        on="name",
        how="inner"
    )

    # Negative difference means LGA scored lower/better than GA
    merged_data["score_diff"] = (
        merged_data["score"] - merged_data["ga_score"]
    )

    # 6. Generate a heatmap for each ligand
    ligands = merged_data["name"].unique()

    for ligand in ligands:
        ligand_data = merged_data[merged_data["name"] == ligand]

        pivot_table = ligand_data.pivot(
            index="iter",
            columns="rate",
            values="score_diff"
        ).sort_index(ascending=False)

        plt.figure(figsize=(8, 6))
        sns.heatmap(
            pivot_table,
            annot=True,
            cmap="coolwarm",
            center=0,
            fmt=".2f",
            cbar_kws={
                "label": f"Score Difference ({args.agg.capitalize()})"
            }
        )

        plt.title(f"{ligand} - LGA vs GA\nScore Difference: LGA - GA")
        plt.xlabel("Local Search Rate")
        plt.ylabel("Local Search Iterations")
        plt.tight_layout()

        out_path = os.path.join(
            args.outdir,
            f"heatmap_{ligand}_{args.agg}.png"
        )
        plt.savefig(out_path, dpi=300)
        plt.close()

        print(f"Saved heatmap for {ligand} to {out_path}")

    # 7. Generate global heatmaps grouped by rotamer count
    # Bounds are inclusive: 1-9 and 10-50 rotamers.
    rotamer_bins = [
        (1, 7, "rot_1_7"),
        (8, 13, "rot_8_13"),
        (14, 50, "rot_14_50"),
    ]

    for rot_min, rot_max, bin_name in rotamer_bins:
        bin_data = merged_data[
            merged_data["rot"].between(rot_min, rot_max, inclusive="both")
        ]

        if bin_data.empty:
            print(
                f"Skipping global heatmap for {bin_name}: "
                f"no complexes with {rot_min}-{rot_max} rotamers."
            )
            continue

        # Average score differences across complexes for each configuration
        global_agg = (
            bin_data.groupby(["rate", "iter"])["score_diff"]
            .mean()
            .reset_index(name="mean_score_diff")
        )

        global_pivot = global_agg.pivot(
            index="iter",
            columns="rate",
            values="mean_score_diff"
        ).sort_index(ascending=False)

        num_complexes = bin_data["name"].nunique()

        plt.figure(figsize=(8, 6))
        sns.heatmap(
            global_pivot,
            annot=True,
            cmap="coolwarm",
            center=0,
            fmt=".2f",
            cbar_kws={
                "label": (
                    f"Mean Score Difference "
                    f"({args.agg.capitalize()} per complex)"
                )
            }
        )

        plt.title(
            f"Global LGA vs GA: {rot_min}-{rot_max} Rotamers\n"
            f"Mean Score Difference Across {num_complexes} Complexes "
            f"({args.agg.capitalize()} per complex)"
        )
        plt.xlabel("Local Search Rate")
        plt.ylabel("Local Search Iterations")
        plt.tight_layout()

        out_path = os.path.join(
            args.outdir,
            f"heatmap_global_{bin_name}_{args.agg}.png"
        )
        plt.savefig(out_path, dpi=300)
        plt.close()

        print(
            f"Saved global heatmap for {bin_name} "
            f"({num_complexes} complexes) to {out_path}"
        )


if __name__ == "__main__":
    main()
