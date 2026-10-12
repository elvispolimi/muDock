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
    df = pd.read_csv(args.csv_file)

    # Ensure numeric columns
    df["rate"] = pd.to_numeric(df["rate"], errors="coerce")
    df["iter"] = pd.to_numeric(df["iter"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")

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

    # LGA scores per ligand and configuration
    lga_agg = (
        lga_runs.groupby(["name", "rate", "iter"])["score"]
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

    # 7. Generate the global heatmap: mean across ligands
    global_agg = (
        merged_data.groupby(["rate", "iter"])["score_diff"]
        .mean()
        .reset_index(name="mean_score_diff")
    )

    global_pivot = global_agg.pivot(
        index="iter",
        columns="rate",
        values="mean_score_diff"
    ).sort_index(ascending=False)

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        global_pivot,
        annot=True,
        cmap="coolwarm",
        center=0,
        fmt=".2f",
        cbar_kws={
            "label": f"Mean Score Difference ({args.agg.capitalize()})"
        }
    )

    plt.title(
        f"Global LGA vs GA\n"
        f"Mean Score Difference Across Ligands ({args.agg.capitalize()} per Ligand)"
    )
    plt.xlabel("Local Search Rate")
    plt.ylabel("Local Search Iterations")
    plt.tight_layout()

    out_path = os.path.join(
        args.outdir,
        f"heatmap_global_mean_{args.agg}.png"
    )
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"Saved global heatmap to {out_path}")


if __name__ == "__main__":
    main()