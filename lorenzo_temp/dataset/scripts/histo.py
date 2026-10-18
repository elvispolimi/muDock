import pandas as pd
import matplotlib.pyplot as plt

# Input CSV file
CSV_FILE = "./lorenzo_temp/info/casf_2016_pocket.csv"

# Read CSV
df = pd.read_csv(CSV_FILE)

# Ensure rot contains numeric values
df["rot"] = pd.to_numeric(df["rot"], errors="coerce")
rot = df["rot"].dropna()

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(
    rot,
    bins=range(int(rot.min()), int(rot.max()) + 2),
    edgecolor="black",
    align="left"
)
plt.xticks(range(int(rot.min()), int(rot.max()) + 1))
plt.xlabel("Number of Rotatable Bonds (rot)")
plt.ylabel("Frequency")
plt.title("Histogram of Rotatable Bonds")
plt.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("rot_histogram.png", dpi=300)
plt.show()

print(f"Total valid entries: {len(rot)}")
print(f"Histogram saved to rot_histogram.png")
