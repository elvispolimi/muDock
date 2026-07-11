#!/usr/bin/env python3
"""
Plot AdaDelta local-search score convergence per ligand.

Reads a CSV file with columns: iteration, ligand, score
(produced by the modified `run_standalone` in adadelta.hpp) and, for each
ligand, plots the score delta relative to its iteration-0 value against
the iteration number.

Usage:
    python plot_scores.py [path/to/adadelta_scores.csv] [--out-dir DIR]

If no path is given, "adadelta_scores.csv" in the current directory is used.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def load_scores(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        sys.exit(f"Error: file not found: {csv_path}")
    df = pd.read_csv(csv_path)
    expected_cols = {"iteration", "ligand", "score"}
    if not expected_cols.issubset(df.columns):
        sys.exit(f"Error: CSV must contain columns {expected_cols}, found {set(df.columns)}")
    return df


def center_on_first_iteration(df: pd.DataFrame) -> pd.DataFrame:
    """Subtract each ligand's score at iteration 0 from all its scores."""
    df = df.sort_values(["ligand", "iteration"]).reset_index(drop=True).copy()

    # First recorded iteration's score, per ligand, broadcast back to every row.
    first_iter_idx = df.groupby("ligand")["iteration"].transform("min")
    is_first = df["iteration"] == first_iter_idx
    baseline_per_ligand = df.loc[is_first].groupby("ligand")["score"].first()

    df["centered_score"] = df["score"] - df["ligand"].map(baseline_per_ligand)
    return df


def plot_per_ligand(df: pd.DataFrame, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    saved = []
    for ligand_id, group in df.groupby("ligand"):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(group["iteration"], group["centered_score"], marker="o", markersize=3, linewidth=1)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Score change from iteration 0")
        ax.set_title(f"Ligand {ligand_id}: local search score convergence")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        out_path = out_dir / f"ligand_{ligand_id}_score.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        saved.append(out_path)
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", nargs="?", default="adadelta_scores.csv",
                         help="Path to the scores CSV file (default: adadelta_scores.csv)")
    parser.add_argument("--out-dir", default="score_plots",
                         help="Directory to write plot images to (default: score_plots)")
    args = parser.parse_args()

    csv_path = Path(args.csv_path)
    out_dir = Path(args.out_dir)

    df = load_scores(csv_path)
    df = center_on_first_iteration(df)

    saved = plot_per_ligand(df, out_dir)
    for path in saved:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
