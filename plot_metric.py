#!/usr/bin/env python3
"""
Plot one curve for each dataset contained in a CSV file.

The CSV must contain exactly three columns:

    <id>, <x-axis>, <y-axis>

Examples:
    ligand,iteration,score
    ligand,iteration,gradient_norm
    ligand,generation,fitness
    id,time,distance

The first column identifies the dataset to plot, the second column is used
as the x-axis, and the third column is the measured value.

Optionally, the y-values can be centered so that the first recorded value of
each dataset becomes zero.

Usage:
    python plot_csv.py data.csv
    python plot_csv.py data.csv --center
    python plot_csv.py data.csv --out-dir plots
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def load_data(csv_path: Path):
    """Load the CSV and determine column names."""

    if not csv_path.exists():
        sys.exit(f"Error: file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    if len(df.columns) != 3:
        raise RuntimeError(
            f"Expected exactly 3 columns, found {len(df.columns)}."
        )

    id_col, x_col, y_col = df.columns

    return df, id_col, x_col, y_col


def center_values(df: pd.DataFrame, id_col: str, x_col: str, y_col: str):
    """Center each dataset on its first recorded value."""

    df = df.sort_values([id_col, x_col]).reset_index(drop=True).copy()

    first_values = (
        df.groupby(id_col)[y_col]
        .first()
    )

    df[y_col] = df[y_col] - df[id_col].map(first_values)

    return df


def plot_data(df, id_col, x_col, y_col, out_dir: Path, centered: bool, output_base: str | None):
    """Generate one plot per dataset."""

    out_dir.mkdir(parents=True, exist_ok=True)

    saved = []

    for dataset_id, group in df.groupby(id_col):

        fig, ax = plt.subplots(figsize=(8, 5))

        ax.plot(
            group[x_col],
            group[y_col],
            marker="o",
            markersize=3,
            linewidth=1,
        )

        if centered:
            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)

        ax.set_xlabel(x_col)

        ylabel = y_col
        if centered:
            ylabel += " (centered)"

        ax.set_ylabel(ylabel)
        ax.set_title(f"{dataset_id}")

        ax.grid(True, alpha=0.3)

        fig.tight_layout()

        if output_base is None:
            filename = f"{dataset_id}_{y_col}.png"
        else:
            filename = f"{output_base}_{dataset_id}.png"

        out_path = out_dir / filename

        fig.savefig(out_path, dpi=150)
        plt.close(fig)

        saved.append(out_path)

    return saved


def main():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "csv_path",
        help="CSV file to plot",
    )

    parser.add_argument(
        "--out-dir",
        default="plots",
        help="Output directory (default: plots)",
    )

    parser.add_argument(
        "--center",
        action="store_true",
        help="Center each curve on its first recorded value",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "Base name for output images. "
            "For example '-o score' produces score_<id>.png."
        ),
    )

    args = parser.parse_args()

    csv_path = Path(args.csv_path)
    out_dir = Path(args.out_dir)

    df, id_col, x_col, y_col = load_data(csv_path)

    if args.center:
        df = center_values(df, id_col, x_col, y_col)

    saved = plot_data(
        df,
        id_col,
        x_col,
        y_col,
        out_dir,
        args.center,
        args.output,
    )

    for path in saved:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()