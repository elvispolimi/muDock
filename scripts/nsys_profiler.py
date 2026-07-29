import subprocess
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

PROTEIN = "/work/onedina/muDock/data/1fkb/1fkb_pocket.pdbqt"
LIGAND = "/work/onedina/muDock/data/1fkb/1fkb_ligand.adtmol2"

BACKENDS = [
    ("CUDA Native", "/work/onedina/muDock/build/cuda/application/muDock", "CUDA:GPU:0"),
    ("Alpaka CUDA", "/work/onedina/muDock/build/alpaka-cuda/application/muDock", "ALPAKA:GPU:0")
]

POPULATION = 100
GENERATIONS = 1000

def run_nsys_profiling():
    gpu_frames = []
    api_frames = []

    print("="*60)
    print(" MICRO PROFILING START (NSIGHT SYSTEMS)")
    print("="*60)

    for name, bin_path, use_flag in BACKENDS:
        if not Path(bin_path).exists():
            print(f"[!] WARNING: Executable not found -> {bin_path}")
            continue

        safe_name = name.replace(" ", "_")
        report_base = f"nsys_report_{safe_name}"
        sqlite_file = f"{report_base}.sqlite"

        print(f"\n[+] Executing nsys profile for {name}...")
        cmd_profile = [
            "nsys", "profile",
            "--force-overwrite=true",
            "--stats=true",
            f"--output={report_base}",
            bin_path,
            "--protein", PROTEIN, "--ligand", LIGAND, "--use", use_flag,
            "--population", str(POPULATION), "--generations", str(GENERATIONS)
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

        cmd_stats_api = [
            "nsys", "stats", "--report", "cudaapisum", "--format", "csv", 
            "--force-overwrite", "true", "--output", report_base, sqlite_file
        ]
        subprocess.run(cmd_stats_api, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        csv_api = f"{report_base}_cudaapisum.csv"

        if Path(csv_kern).exists():
            df_k = pd.read_csv(csv_kern)
            
            def clean_kname(kname):
                kname = str(kname).lower()
                if "calc_energy" in kname or "adt_score" in kname: return "adt_score"
                if "apply" in kname or "geom" in kname: return "geom_transform"
                if "iterate" in kname or "initialize" in kname or "finalize" in kname or "genetic" in kname: return "genetic"
                if "rand" in kname: return "rng_init"
                return "other_gpu"

            df_k["Kernel"] = df_k["Name"].apply(clean_kname)
            
            time_col = [c for c in df_k.columns if "Total Time" in c or "Time (ns)" in c][0]
            inst_col = [c for c in df_k.columns if "Instances" in c][0]
            avg_col = [c for c in df_k.columns if "Avg" in c or "Average" in c][0]
            min_col = [c for c in df_k.columns if "Min" in c][0]
            max_col = [c for c in df_k.columns if "Max" in c][0]
            
            sum_k = df_k.groupby("Kernel").agg({
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
            gpu_frames.append(sum_k)
        
        if Path(csv_api).exists():
            df_a = pd.read_csv(csv_api)
            
            def group_api(apiname):
                apiname = str(apiname)
                if "Launch" in apiname: return "Kernel Launch"
                if "Memcpy" in apiname: return "Memory Copy"
                if "Malloc" in apiname or "Free" in apiname: return "Memory Alloc"
                return "Other API"

            df_a["API_Category"] = df_a["Name"].apply(group_api)
            time_col_a = "Total Time (ns)" if "Total Time (ns)" in df_a.columns else df_a.columns[1]
            
            sum_a = df_a.groupby("API_Category")[time_col_a].sum().reset_index()
            sum_a["Time (ms)"] = sum_a[time_col_a] / 1e6
            sum_a["Backend"] = name
            api_frames.append(sum_a)

    if gpu_frames:
        print("\n[+] Generating Micro-Vision plots...")
        df_gpu = pd.concat(gpu_frames, ignore_index=True)
        unique_kernels = df_gpu['Kernel'].unique()
        
        for k in unique_kernels:
            k_df = df_gpu[df_gpu['Kernel'] == k]
            
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            fig.suptitle(f"Micro-Vision Analysis: {k.upper()}", fontweight='bold', fontsize=20)
            
            sns.barplot(data=k_df, x='Backend', y='Total Time (ms)', ax=axes[0,0], palette='Set1')
            axes[0,0].set_title("Total Execution Time")
            axes[0,0].set_ylabel("Time (ms)")
            
            sns.barplot(data=k_df, x='Backend', y='Avg Time (us)', ax=axes[0,1], palette='Set2')
            axes[0,1].set_title("Average Time per Launch")
            axes[0,1].set_ylabel("Time (us)")
            
            sns.barplot(data=k_df, x='Backend', y='Min Time (us)', ax=axes[1,0], palette='Set3')
            axes[1,0].set_title("Best (Min) Time per Launch")
            axes[1,0].set_ylabel("Time (us)")
            
            sns.barplot(data=k_df, x='Backend', y='Max Time (us)', ax=axes[1,1], palette='Dark2')
            axes[1,1].set_title("Worst (Max) Time per Launch")
            axes[1,1].set_ylabel("Time (us)")
            
            for ax in axes.flat:
                ax.set_xlabel("")
                sns.despine(ax=ax)
                for c in ax.containers:
                    ax.bar_label(c, fmt='%.2f', label_type='edge', padding=3, fontsize=10)
                    
            plt.tight_layout()
            out_name = f"micro_01_kernel_{k}.pdf"
            fig.savefig(out_name)
            print(f"-> Saved: {out_name}")

    if api_frames:
        df_api = pd.concat(api_frames, ignore_index=True)
        
        fig2, ax2 = plt.subplots(figsize=(10, 6))
        pivot_api = df_api.pivot(index='Backend', columns='API_Category', values='Time (ms)').fillna(0)
        pivot_api.plot(kind='bar', stacked=True, ax=ax2, colormap='Set2', width=0.6)
        
        ax2.set_title("Driver-Level Overhead: CPU Time spent in CUDA APIs", pad=20, fontweight='bold')
        ax2.set_ylabel("Total Host Time (ms)")
        ax2.set_xlabel("")
        plt.xticks(rotation=0)
        sns.despine()

        fig2.savefig("micro_02_cuda_apis.pdf")
        print("-> Saved: micro_02_cuda_apis.pdf")

    print("\n[+] ALL MICRO-PROFILING PLOTS GENERATED!")

if __name__ == "__main__":
    run_nsys_profiling()
