
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_throughput(csv_file):
    df = pd.read_csv(csv_file)

    # Remove whitespace from column names
    df.columns = df.columns.str.strip()

    # Sort by average throughput (descending)
    df = df.sort_values("avg_throughput", ascending=False)

    plt.figure(figsize=(10, 6))

    bars = plt.bar(
        df["software"],
        df["avg_throughput"],
        width=0.6,
    )

    # Display throughput above each bar
    for bar in bars:
        height = bar.get_height()
        plt.annotate(
            f"{height:.2f}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
        )

    plt.xlabel("Application")
    plt.ylabel("Average throughput (lig/s)")
    plt.title("")
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.xticks()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", help="Path to CSV file")
    args = parser.parse_args()

    plot_throughput(args.csv_file)
