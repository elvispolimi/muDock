import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

RESULTS_DIR = "/work/onedina/muDock_ON/scripts"
CSV_FILE = "macro_results_cpu_small_multi.csv"

# Read the data
csv_path = os.path.join(RESULTS_DIR, CSV_FILE)
if not os.path.exists(csv_path):
    print(f"[-] {CSV_FILE} not found in {RESULTS_DIR}")
    exit(1)

df = pd.read_csv(csv_path)

# Average across the runs
df_avg = df.groupby(['Backend', 'Population', 'Generations'])['Throughput (Evals/s)'].mean().reset_index()

sns.set_theme(style="whitegrid")

# ---------------------------------------------------------
# 1. HEATMAP: The Cost of Unrolling (Alpaka Serial No-Unroll vs Unroll)
# ---------------------------------------------------------
df_alp_unroll = df_avg[df_avg['Backend'] == 'Alpaka Serial (Unroll)'].set_index(['Population', 'Generations'])
df_alp_no = df_avg[df_avg['Backend'] == 'Alpaka Serial (No Unroll)'].set_index(['Population', 'Generations'])

# Speedup: (No-Unroll / Unroll)
speedup_unroll = (df_alp_no['Throughput (Evals/s)'] / df_alp_unroll['Throughput (Evals/s)']).reset_index()
heatmap_data_unroll = speedup_unroll.pivot(index="Population", columns="Generations", values="Throughput (Evals/s)")

plt.figure(figsize=(8, 6))
sns.heatmap(heatmap_data_unroll, annot=True, fmt=".2f", cmap="RdYlGn", cbar_kws={'label': 'Speedup Factor'})
plt.title("The Cost of Unrolling on CPU (Dataset: Small Multi)\nSpeedup: Alpaka (No-Unroll) / Alpaka (Unroll)")
plt.savefig(os.path.join(RESULTS_DIR, "cpu_advanced_01_unroll_penalty_heatmap_small_multi.pdf"), bbox_inches='tight')
plt.close()

# ---------------------------------------------------------
# 2. HEATMAP: The RNG Gap (CPP Serial vs Alpaka Serial No-Unroll)
# ---------------------------------------------------------
df_cpp = df_avg[df_avg['Backend'] == 'CPP Serial'].set_index(['Population', 'Generations'])

# Speedup: (Alpaka No-Unroll / CPP Serial)
speedup_rng = (df_alp_no['Throughput (Evals/s)'] / df_cpp['Throughput (Evals/s)']).reset_index()
heatmap_data_rng = speedup_rng.pivot(index="Population", columns="Generations", values="Throughput (Evals/s)")

plt.figure(figsize=(8, 6))
sns.heatmap(heatmap_data_rng, annot=True, fmt=".2f", cmap="Reds_r", cbar_kws={'label': 'Performance Ratio (1.0 = Parity)'})
plt.title("The RNG Gap (Philox vs MT) (Dataset: Small Multi)\nRatio: Alpaka (No-Unroll) / CPP Serial")
plt.savefig(os.path.join(RESULTS_DIR, "cpu_advanced_02_rng_gap_heatmap_small_multi.pdf"), bbox_inches='tight')
plt.close()

# ---------------------------------------------------------
# 3. OVERALL BAR PLOT: Average Throughput Collapse
# ---------------------------------------------------------
plt.figure(figsize=(10, 6))
avg_global = df_avg.groupby('Backend')['Throughput (Evals/s)'].mean().sort_values(ascending=False).reset_index()
sns.barplot(data=avg_global, x='Throughput (Evals/s)', y='Backend', palette="viridis")
plt.title("Global Average Throughput (Dataset: Small Multi)")
for index, value in enumerate(avg_global['Throughput (Evals/s)']):
    plt.text(value, index, f" {int(value)} Evals/s", va='center')
plt.xlim(0, avg_global['Throughput (Evals/s)'].max() * 1.15)
plt.savefig(os.path.join(RESULTS_DIR, "cpu_advanced_03_global_averages_small_multi.pdf"), bbox_inches='tight')
plt.close()

print("[+] Generated 3 advanced analytical plots for small_multi.")
