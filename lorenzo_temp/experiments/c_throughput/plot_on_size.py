
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_throughput(csv_file):
    df = pd.read_csv(csv_file)

    # Remove whitespace from column names
    df.columns = df.columns.str.strip()

    # Plot one curve per software
    plt.figure(figsize=(10, 6))

    for software, group in df.groupby("software"):
        group = group.sort_values("data_size")

        plt.plot(
            group["data_size"],
            group["avg_throughput"],
            marker="o",
            linewidth=2,
            label=software,
        )

    plt.xlabel("Data size")
    plt.ylabel("Average throughput")
    plt.title("Throughput vs Data Size")
    plt.legend(title="Software")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", help="Path to CSV file")
    args = parser.parse_args()

    plot_throughput(args.csv_file)
