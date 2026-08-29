#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Input format
# ---------------------------------------------------------------------------

COLUMNS = [
    "experiment",
    "ligand_id",
    "score",
    "generations",
    "seed",
    "num_rotamers",
    "num_atoms",
    "num_evaluations",
    "population_size",
    "ls_rate",
    "ls_iter",
]

NUMERIC_COLUMNS = [
    "score",
    "generations",
    "seed",
    "num_rotamers",
    "num_atoms",
    "num_evaluations",
    "population_size",
    "ls_rate",
    "ls_iter",
]


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_results(path: str | Path) -> pd.DataFrame:
    """
    Load experiment results from a whitespace-separated file.

    Expected format:

        Experiment ligand_id score generations seed num_rotamers
        num_atoms num_evaluations population_size ls_rate ls_iter

    The input file has no header.
    """

    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=COLUMNS,
    )

    # Remove exact duplicate rows.
    df = df.drop_duplicates().reset_index(drop=True)

    # Convert numerical columns.
    df[NUMERIC_COLUMNS] = df[NUMERIC_COLUMNS].apply(
        pd.to_numeric,
        errors="raise",
    )

    return df


# ---------------------------------------------------------------------------
# Statistics across random seeds
# ---------------------------------------------------------------------------

def summarize_runs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate different random-seed runs belonging to the same
    ligand/configuration.

    Each group is identified by:
        ligand_id
        ls_rate
        ls_iter

    The resulting dataframe contains statistics for both
    generations and score.
    """

    group_columns = [
        "ligand_id",
        "ls_rate",
        "ls_iter",
    ]

    summary = (
        df.groupby(group_columns)
        .agg(
            # Generations
            mean_generations=("generations", "mean"),
            std_generations=("generations", "std"),
            median_generations=("generations", "median"),
            min_generations=("generations", "min"),
            max_generations=("generations", "max"),

            # Score
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            median_score=("score", "median"),
            min_score=("score", "min"),
            max_score=("score", "max"),

            # Number of runs
            n_runs=("generations", "count"),
        )
        .reset_index()
    )

    return summary


# ---------------------------------------------------------------------------
# Heatmap data
# ---------------------------------------------------------------------------

def get_heatmap_data(
    summary: pd.DataFrame,
    ligand_id: str,
    value: str = "mean_generations",
) -> pd.DataFrame:
    """
    Return a matrix suitable for a heatmap.

    Rows    -> ls_iter
    Columns -> ls_rate
    Values  -> selected statistic
    """

    ligand_summary = summary[
        summary["ligand_id"] == ligand_id
    ]

    if ligand_summary.empty:
        raise ValueError(
            f"No data found for ligand '{ligand_id}'."
        )

    heatmap_data = ligand_summary.pivot(
        index="ls_iter",
        columns="ls_rate",
        values=value,
    )

    # Sort axes numerically.
    heatmap_data = heatmap_data.sort_index()
    heatmap_data = heatmap_data.sort_index(axis=1)

    return heatmap_data


# ---------------------------------------------------------------------------
# Speedup relative to baseline
# ---------------------------------------------------------------------------

def calculate_speedup(summary: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate convergence speedup relative to ls_rate = 0.

    speedup = baseline_generations / configuration_generations

    Therefore:
        1.0 -> same as baseline
        >1  -> faster convergence
        <1  -> slower convergence
    """

    baseline = (
        summary[summary["ls_rate"] == 0]
        [
            [
                "ligand_id",
                "ls_iter",
                "mean_generations",
            ]
        ]
        .rename(
            columns={
                "mean_generations": "baseline_generations"
            }
        )
    )

    result = summary.merge(
        baseline,
        on=["ligand_id", "ls_iter"],
        how="left",
    )

    result["speedup"] = (
        result["baseline_generations"]
        / result["mean_generations"]
    )

    return result


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_heatmap(
    heatmap_data: pd.DataFrame,
    ligand_id: str,
    title: str,
    colorbar_label: str,
    output_path: str | Path | None = None,
    value_format: str = ".1f",
):
    """
    Plot a heatmap from a dataframe.

    Rows    -> ls_iter
    Columns -> ls_rate
    Values  -> dataframe values
    """

    fig, ax = plt.subplots(
        figsize=(10, 7)
    )

    image = ax.imshow(
        heatmap_data.values,
        aspect="auto",
        origin="lower",
    )

    # X axis
    ax.set_xticks(
        range(len(heatmap_data.columns))
    )
    ax.set_xticklabels(
        heatmap_data.columns
    )

    # Y axis
    ax.set_yticks(
        range(len(heatmap_data.index))
    )
    ax.set_yticklabels(
        heatmap_data.index
    )

    ax.set_xlabel("LS rate (%)")
    ax.set_ylabel("LS iterations")

    ax.set_title(title)

    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(colorbar_label)

    # Write values inside cells.
    for i in range(len(heatmap_data.index)):
        for j in range(len(heatmap_data.columns)):
            value = heatmap_data.iloc[i, j]

            if pd.notna(value):
                ax.text(
                    j,
                    i,
                    f"{value:{value_format}}",
                    ha="center",
                    va="center",
                )

    fig.tight_layout()

    if output_path is not None:
        fig.savefig(
            output_path,
            dpi=300,
            bbox_inches="tight",
        )
        print(f"Saved heatmap to {output_path}")

    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyze muDock GA experiment results."
    )

    parser.add_argument(
        "input",
        help="Path to the experiment results file.",
    )

    parser.add_argument(
        "--ligand",
        help="Ligand to analyze. If omitted, all ligands are listed.",
    )

    parser.add_argument(
        "--output",
        "-o",
        help="Output path for the heatmap.",
    )

    args = parser.parse_args()

    # -----------------------------------------------------------------------
    # Load
    # -----------------------------------------------------------------------

    df = load_results(args.input)

    print("\nLoaded results:")
    print(df)

    print("\nNumber of rows:", len(df))

    # -----------------------------------------------------------------------
    # Show available ligands
    # -----------------------------------------------------------------------

    ligands = sorted(
        df["ligand_id"].unique()
    )

    print("\nLigands:")
    for ligand in ligands:
        print(f"  {ligand}")

    # -----------------------------------------------------------------------
    # Summarize stochastic runs
    # -----------------------------------------------------------------------

    summary = summarize_runs(df)

    print("\nSummary across seeds:")
    print(summary.to_string(index=False))

    # Save summary next to input file.
    input_path = Path(args.input)

    summary_path = (
        input_path.parent
        / f"{input_path.stem}_summary.csv"
    )

    summary.to_csv(
        summary_path,
        index=False,
    )

    print(
        f"\nSaved summary to {summary_path}"
    )

    # -----------------------------------------------------------------------
    # If no ligand was specified, stop after printing the summary.
    # -----------------------------------------------------------------------

    if args.ligand is None:
        print(
            "\nUse --ligand <ligand_id> to generate "
            "a heatmap."
        )
        return

    # -----------------------------------------------------------------------
    # Check ligand
    # -----------------------------------------------------------------------

    if args.ligand not in ligands:
        raise ValueError(
            f"Ligand '{args.ligand}' not found."
        )

    # -----------------------------------------------------------------------
    # Mean generations heatmap
    # -----------------------------------------------------------------------

    heatmap_data = get_heatmap_data(
        summary,
        args.ligand,
        value="mean_generations",
    )

    print(
        f"\nMean generations heatmap data "
        f"for {args.ligand}:"
    )
    print(heatmap_data)

    median_heatmap_data = get_heatmap_data(
        summary,
        args.ligand,
        value="median_generations",
    )

    print(
        f"\nMedian generations heatmap data "
        f"for {args.ligand}:"
    )
    print(median_heatmap_data)

    score_heatmap_data = get_heatmap_data(
        summary,
        args.ligand,
        value="mean_score",
    )

    print(
        f"\nMean score heatmap data "
        f"for {args.ligand}:"
    )
    print(score_heatmap_data)

    median_score_heatmap_data = get_heatmap_data(
        summary,
        args.ligand,
        value="median_score",
    )

    print(
        f"\nMedian score heatmap data "
        f"for {args.ligand}:"
    )
    print(median_score_heatmap_data)

    # -----------------------------------------------------------------------
    # Speedup
    # -----------------------------------------------------------------------

    speedup = calculate_speedup(summary)

    ligand_speedup = speedup[
        speedup["ligand_id"] == args.ligand
    ]

    print(
        f"\nSpeedup relative to ls_rate=0 "
        f"for {args.ligand}:"
    )

    print(
        ligand_speedup[
            [
                "ls_rate",
                "ls_iter",
                "mean_generations",
                "baseline_generations",
                "speedup",
            ]
        ].to_string(index=False)
    )

    # -----------------------------------------------------------------------
    # Plot mean generations heatmap
    # -----------------------------------------------------------------------

    if args.output is None:
        generations_output_path = (
            input_path.parent
            / f"{input_path.stem}_{args.ligand}_generations_heatmap.png"
        )
    else:
        generations_output_path = args.output

    plot_heatmap(
        heatmap_data,
        args.ligand,
        title=(
            f"{args.ligand} - "
            "Mean GA convergence generations"
        ),
        colorbar_label="Mean generations",
        output_path=generations_output_path,
    )

    # -----------------------------------------------------------------------
    # Plot median generations heatmap
    # -----------------------------------------------------------------------

    if args.output is None:
        median_generations_output_path = (
            input_path.parent
            / f"{input_path.stem}_{args.ligand}_median_generations_heatmap.png"
        )
    else:
        median_generations_output_path = args.output

    plot_heatmap(
        median_heatmap_data,
        args.ligand,
        title=(
            f"{args.ligand} - "
            "Median GA convergence generations"
        ),
        colorbar_label="Median generations",
        output_path=median_generations_output_path,
    )

    # -----------------------------------------------------------------------
    # Plot mean score heatmap
    # -----------------------------------------------------------------------

    score_output_path = (
        input_path.parent
        / f"{input_path.stem}_{args.ligand}_score_heatmap.png"
    )

    plot_heatmap(
        score_heatmap_data,
        args.ligand,
        title=(
            f"{args.ligand} - "
            "Mean final score"
        ),
        colorbar_label="Mean score",
        output_path=score_output_path,
    )

    # -----------------------------------------------------------------------
    # Plot median score heatmap
    # -----------------------------------------------------------------------

    median_score_output_path = (
        input_path.parent
        / f"{input_path.stem}_{args.ligand}_median_score_heatmap.png"
    )

    plot_heatmap(
        median_score_heatmap_data,
        args.ligand,
        title=(
            f"{args.ligand} - "
            "Median final score"
        ),
        colorbar_label="Median score",
        output_path=median_score_output_path,
    )


if __name__ == "__main__":
    main()