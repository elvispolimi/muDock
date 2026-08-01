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

DATASETS = ["single", "small_multi"]

def generate_macro_plots():
    for ds_name in DATASETS:
        csv_file = f"profiling_macro_results_{ds_name}.csv"
        
        if not Path(csv_file).exists():
            print(f"[!] ERROR: File {csv_file} not found. Skipping dataset {ds_name}.")
            continue

        df = pd.read_csv(csv_file)
        if df.empty:
            print(f"[!] CSV file {csv_file} is empty. Skipping dataset {ds_name}.")
            continue

        print("="*60)
        print(f" GENERATING MACRO PLOTS (PAPER QUALITY) - Dataset: {ds_name.upper()}")
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

        ax1.set_title(f"Performance Scaling: Execution Latency ({ds_name})", pad=20, fontweight='bold')
        ax1.set_xlabel("Number of Generations")
        ax1.set_ylabel("Execution Time (Seconds) [Lower is Better]")
        sns.despine()
        fig1.savefig(f"macro_01_latency_scaling_{ds_name}.pdf")
        plt.close(fig1)
        
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
        
        ax2.set_title(f"Hardware Utilization: Throughput Scaling ({ds_name})", pad=20, fontweight='bold')
        ax2.set_xlabel("Total Workload (Population × Generations)")
        ax2.set_ylabel("Throughput (Evaluations / Second) [Higher is Better]")
        ax2.set_xscale("log")
        sns.despine()
        fig2.savefig(f"macro_02_throughput_{ds_name}.pdf")
        plt.close(fig2)

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
            
            ax3.set_title(f"Alpaka relative Speedup over CUDA Native ({ds_name})", pad=20, fontweight='bold')
            ax3.invert_yaxis()
            fig3.savefig(f"macro_03_speedup_heatmap_{ds_name}.pdf")
            plt.close(fig3)
        else:
            print("[!] Cannot generate Speedup Heatmap (Incomplete data).")

        print("[4/4] Generating Statistical Table PDF (Mean & Variance)...")
        df_stats = df.groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"].agg(['mean', 'var']).reset_index()
        df_stats.rename(columns={'mean': 'Throughput_Mean', 'var': 'Throughput_Variance'}, inplace=True)
        
        # Save the table to CSV
        stats_csv = f"macro_04_statistics_table_{ds_name}.csv"
        df_stats.to_csv(stats_csv, index=False)
        
        # Arrotondiamo per estetica PDF
        df_stats['Throughput_Mean'] = df_stats['Throughput_Mean'].map('{:.2f}'.format)
        df_stats['Throughput_Variance'] = df_stats['Throughput_Variance'].map('{:.2f}'.format)
        
        fig4, ax4 = plt.subplots(figsize=(10, 12))
        ax4.axis('tight')
        ax4.axis('off')
        
        table = ax4.table(cellText=df_stats.values, colLabels=df_stats.columns, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        
        # Colora l'header
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_text_props(weight='bold')
                cell.set_facecolor('#d3d3d3')
                
        plt.title(f"Statistical Summary: Throughput Mean & Variance ({ds_name})", pad=20, fontweight='bold', fontsize=16)
        
        pdf_name = f"macro_04_statistics_table_{ds_name}.pdf"
        fig4.savefig(pdf_name, bbox_inches='tight')
        plt.close(fig4)
        
        print(f"-> Saved: {pdf_name}")
        print(f"-> Saved: {stats_csv}")

        print(f"\n[+] MACRO plots and statistics generated successfully for {ds_name}.")

if __name__ == "__main__":
    generate_macro_plots()
