import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Compare GA and LGA success rates using two criteria."
    )
    parser.add_argument("csv", help="Input CSV containing GA and LGA runs")
    parser.add_argument(
        "--output-dir",
        default="success_rate_plots",
        help="Directory where plots will be saved (default: success_rate_plots)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read CSV and normalize column names.
    df = pd.read_csv(args.csv)
    df.columns = df.columns.str.strip().str.strip("*")

    required = {
        "rate",
        "iter",
        "score",
        "cry_score",
        "rmsd_best_scoring_pose",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required CSV columns: {sorted(missing)}")

    # Convert relevant columns to numeric.
    # Values such as N/A become NaN.
    numeric_cols = [
        "rate",
        "iter",
        "score",
        "cry_score",
        "rmsd_best_scoring_pose",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # GA: rate and iter are N/A.
    # LGA: rate and iter are numeric.
    ga_mask = df["rate"].isna() & df["iter"].isna()
    lga_mask = df["rate"].notna() & df["iter"].notna()

    runs = {
        "GA": df[ga_mask].copy(),
        "LGA": df[lga_mask].copy(),
    }

    unclassified = df[~(ga_mask | lga_mask)]
    if not unclassified.empty:
        print(
            f"Warning: ignoring {len(unclassified)} runs with "
            "inconsistent rate/iter values."
        )

    # Define the two success criteria.
    criteria = [
        {
            "column": "rmsd_best_scoring_pose",
            "threshold": 2.0,
            "operator": "<",
            "title": "Success Rate: RMSD < 2 Å",
            "filename": "success_rate_rmsd.png",
            "success": lambda data: data["rmsd_best_scoring_pose"] < 2,
        },
        {
            "column": "score - cry_score",
            "threshold": 1.0,
            "operator": "<",
            "title": "Success Rate: Score - Crystal Score < 1",
            "filename": "success_rate_score_difference.png",
            "success": lambda data: (
                data["score"] - data["cry_score"]
            ) < 1,
        },
    ]

    # Print summary and generate one plot per criterion.
    for criterion in criteria:
        labels = []
        rates = []

        print(f"\n{criterion['title']}")
        print("-" * len(criterion["title"]))

        for algorithm, data in runs.items():
            total = len(data)

            # Missing metric values count as unsuccessful runs.
            successful = int(criterion["success"](data).fillna(False).sum()) \
                if total else 0

            rate = 100 * successful / total if total else 0.0

            print(
                f"{algorithm}: {successful}/{total} successful "
                f"({rate:.2f}%)"
            )

            labels.append(algorithm)
            rates.append(rate)

        fig, ax = plt.subplots(figsize=(6, 5))

        bars = ax.bar(labels, rates, width=0.55)

        ax.set_title(criterion["title"])
        ax.set_ylabel("Success rate (%)")
        ax.set_ylim(0, 105)
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.set_axisbelow(True)

        # Annotate each bar with its percentage.
        for bar, rate in zip(bars, rates):
            ax.annotate(
                f"{rate:.2f}%",
                xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                va="bottom",
            )

        fig.tight_layout()

        output_path = output_dir / criterion["filename"]
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Plot saved to: {output_path}")

        plt.show()
        plt.close(fig)


if __name__ == "__main__":
    main()