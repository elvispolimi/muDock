import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

sns.set_theme(style="ticks", context="paper")
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'axes.grid': True,
    'grid.alpha': 0.5,
    'legend.fontsize': 12,
    'figure.titlesize': 18,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

CSV_FILE = "profiling_macro_results.csv"

def generate_macro_plots():
    if not Path(CSV_FILE).exists():
        print(f"[!] ERROR: File {CSV_FILE} not found.")
        return

    df = pd.read_csv(CSV_FILE)
    if df.empty:
        print("[!] CSV file is empty.")
        return

    print("="*60)
    print(" GENERATING MACRO PLOTS (PAPER QUALITY)")
    print("="*60)

    print("[1/3] Generating Latency Scaling...")
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    
    sns.lineplot(
        data=df, 
        x="Generations", 
        y="Time (s)", 
        hue="Backend", 
        style="Population", 
        markers=True, 
        dashes=True,
        linewidth=2,
        markersize=8,
        ax=ax1,
        palette="Set1"
    )

    ax1.set_title("Performance Scaling: Execution Latency", pad=20, fontweight='bold')
    ax1.set_xlabel("Number of Generations")
    ax1.set_ylabel("Execution Time (Seconds) [Lower is Better]")
    sns.despine()
    fig1.savefig("macro_01_latency_scaling.pdf")
    
    print("[2/3] Generating Throughput Analysis...")
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    
    sns.lineplot(
        data=df,
        x="Total Evaluations",
        y="Throughput (Evals/s)",
        hue="Backend",
        style="Population",
        markers=["o", "s", "D", "X"],
        linewidth=2.5,
        markersize=9,
        ax=ax2,
        palette="Set1"
    )
    
    ax2.set_title("Hardware Utilization: Throughput Scaling", pad=20, fontweight='bold')
    ax2.set_xlabel("Total Workload (Population × Generations)")
    ax2.set_ylabel("Throughput (Evaluations / Second) [Higher is Better]")
    ax2.set_xscale("log")
    sns.despine()
    fig2.savefig("macro_02_throughput.pdf")

    print("[3/3] Generating Relative Speedup...")
    df_mean = df.groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"].mean().reset_index()
    
    df_cuda = df_mean[df_mean["Backend"] == "CUDA Native"].set_index(["Population", "Generations"])
    df_alpaka = df_mean[df_mean["Backend"] == "Alpaka CUDA"].set_index(["Population", "Generations"])
    
    if not df_cuda.empty and not df_alpaka.empty:
        df_speedup = (df_alpaka["Throughput (Evals/s)"] / df_cuda["Throughput (Evals/s)"]).reset_index()
        df_speedup.rename(columns={"Throughput (Evals/s)": "Speedup_Alpaka_vs_CUDA"}, inplace=True)
        
        heatmap_data = df_speedup.pivot(index="Population", columns="Generations", values="Speedup_Alpaka_vs_CUDA")
        
        fig3, ax3 = plt.subplots(figsize=(9, 7))
        sns.heatmap(
            heatmap_data, 
            annot=True, 
            fmt=".3f", 
            cmap="RdYlGn",
            center=1.0, 
            linewidths=.5, 
            cbar_kws={'label': 'Speedup (1.0 = Parity with CUDA)'},
            ax=ax3
        )
        
        ax3.set_title("Alpaka relative Speedup over CUDA Native", pad=20, fontweight='bold')
        ax3.invert_yaxis()
        fig3.savefig("macro_03_speedup_heatmap.pdf")
    else:
        print("[!] Cannot generate Speedup Heatmap (Incomplete data).")

    print("\n[+] MACRO plots generated successfully in the current directory.")

if __name__ == "__main__":
    generate_macro_plots()
