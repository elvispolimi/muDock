#!/usr/bin/env python3
"""
Plot the evolution of a local-search metric for each ligand.

Reads a CSV file produced by muDock's CSV logger. The file must contain the
columns:

    iteration, ligand, <metric>

where <metric> is any single measured quantity (e.g. score, COM distance,
gradient norm, RMSD, etc.).

The script automatically detects the metric column, centers its values with
respect to the first recorded iteration for each ligand, and generates one
plot per ligand showing the metric variation over the course of the local
search.

Usage:
    python plot_metric.py [path/to/data.csv] [--out-dir DIR]

Examples:
    python plot_metric.py adadelta_scores.csv
    python plot_metric.py adadelta_com.csv
    python plot_metric.py gradient_norm.csv

If no CSV path is provided, "adadelta_scores.csv" in the current directory is
used.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def load_data(csv_path: Path) -> tuple[pd.DataFrame, str]:
    if not csv_path.exists():
        sys.exit(f"Error: file not found: {csv_path}")
    df = pd.read_csv(csv_path)
    metadata = [c for c in df.columns if c not in ("iteration", "ligand")]
    if len(metadata) != 1:
        raise RuntimeError("Expected exactly one metric column.")
    metric = metadata[0]
    expected_cols = {"iteration", "ligand", metric}
    if not expected_cols.issubset(df.columns):
        sys.exit(f"Error: CSV must contain columns {expected_cols}, found {set(df.columns)}")
    return df, metric


def center_on_first_iteration(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Subtract each ligand's data at iteration 0 from all its data."""
    df = df.sort_values(["ligand", "iteration"]).reset_index(drop=True).copy()

    # First recorded iteration's data, per ligand, broadcast back to every row.
    first_iter_idx = df.groupby("ligand")["iteration"].transform("min")
    is_first = df["iteration"] == first_iter_idx
    baseline_per_ligand = df.loc[is_first].groupby("ligand")[metric].first()

    df["centered_value"] = df[metric] - df["ligand"].map(baseline_per_ligand)    
    return df


def plot_per_ligand(df: pd.DataFrame, out_dir: Path, metric: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    saved = []
    for ligand_id, group in df.groupby("ligand"):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(group["iteration"], group["centered_value"], marker="o", markersize=3, linewidth=1)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Iteration")
        ax.set_ylabel(f"{metric} change from iteration 0")
        ax.set_title(f"Ligand {ligand_id}: local search {metric} convergence")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        out_path = out_dir / f"ligand_{ligand_id}_{metric}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        saved.append(out_path)
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", 
                        nargs="?", 
                        default="adadelta_scores.csv",
                        help="Path to the data CSV file (default: adadelta_scores.csv)")
    parser.add_argument("--out-dir", 
                        default="plots",
                        help="Directory to write plot images to (default: plots)")
    args = parser.parse_args()

    csv_path = Path(args.csv_path)
    out_dir = Path(args.out_dir)

    df, metric = load_data(csv_path)
    df = center_on_first_iteration(df, metric)

    saved = plot_per_ligand(df, out_dir, metric)
    for path in saved:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
