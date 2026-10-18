
import pandas as pd

df = pd.read_csv("lorenzo_temp/info/casf2016_rotatable_bonds.csv")

# Adjust this column name if needed.
col = "rot"

bins = {
    "low (1-9)": df[col].between(1, 9),
    "medium (10-16)": df[col].between(10, 16),
    "high (17+)": df[col] >= 17,
}

print("Actual distribution:")
for name, mask in bins.items():
    n = int(mask.sum())
    print(f"{name:16s}: {n:3d} ({100*n/len(df):.1f}%)")

# Allocate 20 tuning ligands proportionally using
# the largest-remainder method.
counts = {name: int(mask.sum()) for name, mask in bins.items()}
quotas = {name: 20 * n / len(df) for name, n in counts.items()}
allocation = {name: int(q) for name, q in quotas.items()}

remaining = 20 - sum(allocation.values())
for name in sorted(
    quotas, key=lambda k: quotas[k] - allocation[k], reverse=True
)[:remaining]:
    allocation[name] += 1

print("\nProposed allocation:")
for name, n in counts.items():
    tune = allocation[name]
    print(
        f"{name:16s}: tuning={tune:3d}, "
        f"evaluation={n-tune:3d}"
    )
