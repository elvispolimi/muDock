import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

sns.set_theme(style="ticks", context="paper")
plt.rcParams.update({
    'font.size': 12, 'axes.labelsize': 14, 'axes.titlesize': 16,
    'axes.grid': True, 'grid.alpha': 0.5, 'legend.fontsize': 12,
    'figure.titlesize': 18, 'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

PROTEIN = "/work/onedina/muDock/data/1fkb/1fkb_pocket.pdbqt"
LIGAND = "/work/onedina/muDock/data/1fkb/1fkb_ligand.adtmol2"

BACKENDS = [
    ("CUDA Native", "/work/onedina/muDock/build/cuda/application/muDock", "CUDA:GPU:0"),
    ("Alpaka CUDA", "/work/onedina/muDock/build/alpaka-cuda/application/muDock", "ALPAKA:GPU:0")
]

# RAPID LOAD BENCHMARK
POPULATION = 100
GENERATIONS = 1000

def run_genetic_profiling():
    gpu_frames = []

    print("="*60)
    print(" MICRO PROFILING: GENETIC FOCUS (NSIGHT SYSTEMS)")
    print(f" Parameters: Pop={POPULATION}, Gen={GENERATIONS}")
    print("="*60)

    for name, bin_path, use_flag in BACKENDS:
        if not Path(bin_path).exists():
            print(f"[!] WARNING: Executable not found -> {bin_path}")
            continue

        safe_name = name.replace(" ", "_")
        report_base = f"nsys_genetic_{safe_name}"
        sqlite_file = f"{report_base}.sqlite"

        print(f"\n[+] Executing nsys profile for {name}...")
        
        cmd_profile = [
            "nsys", "profile",
            "--force-overwrite=true",
            "--stats=true",
            f"--output={report_base}",
            bin_path,
            "--protein", PROTEIN, "--ligand", LIGAND, "--use", use_flag,
            "--population", str(POPULATION), "--generations", str(GENERATIONS),
            "--observer", "1"
        ]
        
        try:
            subprocess.run(cmd_profile, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[!] nsys profile error for {name}: {e}")
            continue

        cmd_stats_kern = [
            "nsys", "stats", "--report", "gpukernsum", "--format", "csv", 
            "--force-overwrite", "true", "--output", report_base, sqlite_file
        ]
        subprocess.run(cmd_stats_kern, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        csv_kern = f"{report_base}_gpukernsum.csv"

        if Path(csv_kern).exists():
            df_k = pd.read_csv(csv_kern)
            
            # Filtriamo SOLO i kernel del genetic
            df_geom = df_k[df_k["Name"].str.contains("iterate|initialize|finalize|genetic", case=False, na=False)].copy()
            
            if df_geom.empty:
                print(f"[!] No genetic kernels found for {name}")
                continue
                
            df_geom["Kernel"] = "genetic"
            
            time_col = [c for c in df_geom.columns if "Total Time" in c or "Time (ns)" in c][0]
            inst_col = [c for c in df_geom.columns if "Instances" in c][0]
            avg_col = [c for c in df_geom.columns if "Avg" in c or "Average" in c][0]
            min_col = [c for c in df_geom.columns if "Min" in c][0]
            max_col = [c for c in df_geom.columns if "Max" in c][0]
            
            sum_k = df_geom.groupby("Kernel").agg({
                time_col: 'sum',
                inst_col: 'sum',
                avg_col: 'mean',
                min_col: 'min',
                max_col: 'max'
            }).reset_index()
            
            sum_k["Total Time (ms)"] = sum_k[time_col] / 1e6
            sum_k["Avg Time (us)"] = sum_k[avg_col] / 1e3
            sum_k["Min Time (us)"] = sum_k[min_col] / 1e3
            sum_k["Max Time (us)"] = sum_k[max_col] / 1e3
            sum_k["Backend"] = name
            
            print(f"  -> {name} genetic Total Time: {sum_k['Total Time (ms)'].iloc[0]:.2f} ms")
            print(f"  -> {name} genetic Avg Time:   {sum_k['Avg Time (us)'].iloc[0]:.2f} us")
            
            gpu_frames.append(sum_k)

    if gpu_frames:
        print("\n[+] Generating Genetic Plot...")
        df_gpu = pd.concat(gpu_frames, ignore_index=True)
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        fig.suptitle(f"Micro-Vision Analysis: GENETIC (Pop:{POPULATION}, Gen:{GENERATIONS})", fontweight='bold', fontsize=16)
        
        sns.barplot(data=df_gpu, x='Backend', y='Total Time (ms)', ax=axes[0], palette='Set1')
        axes[0].set_title("Total Execution Time")
        axes[0].set_ylabel("Time (ms)")
        
        sns.barplot(data=df_gpu, x='Backend', y='Avg Time (us)', ax=axes[1], palette='Set2')
        axes[1].set_title("Average Time per Launch")
        axes[1].set_ylabel("Time (us)")
        
        for ax in axes.flat:
            ax.set_xlabel("")
            sns.despine(ax=ax)
            for c in ax.containers:
                ax.bar_label(c, fmt='%.2f', label_type='edge', padding=3, fontsize=10)
                
        plt.tight_layout()
        out_name = f"focus_genetic_perf.pdf"
        fig.savefig(out_name)
        print(f"-> Saved: {out_name}")

if __name__ == "__main__":
    run_genetic_profiling()
