import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io

REPO_DIR = "/work/onedina/muDock_ON"
CUDA_BIN = os.path.join(REPO_DIR, "build/cuda/application/muDock")
ALPAKA_BIN = os.path.join(REPO_DIR, "build/alpaka-cuda/application/muDock")
PROTEIN = os.path.join(REPO_DIR, "data/1fkb/1fkb_pocket.pdbqt")
LIGAND = os.path.join(REPO_DIR, "data/1fkb/1fkb_ligand.adtmol2")

KERNEL_REGEX = "calc_energy|iterate|apply"
METRICS = "sm__warps_active.avg.pct_of_peak_sustained_active,sm__throughput.avg.pct_of_peak_sustained_elapsed,gpu__compute_memory_throughput.avg.pct_of_peak_sustained_elapsed,launch__registers_per_thread"

def run_ncu(binary, backend_name, use_flag):
    cmd = [
        "ncu", "--csv", "--page", "details",
        "--metrics", METRICS,
        binary,
        "--use", use_flag,
        "--protein", PROTEIN,
        "--ligand", LIGAND,
        "--population", "4000",
        "--generations", "1",
        "--seed", "42"
    ]
    print(f"Running NCU for {backend_name}... (this may take a minute or two)")
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return res.stdout

def parse_ncu_csv(output, backend):
    lines = output.split('\n')
    csv_lines = []
    in_csv = False
    
    for line in lines:
        if line.startswith('"ID"') or line.startswith('ID,'):
            in_csv = True
        if in_csv:
            csv_lines.append(line)
            
    if not csv_lines:
        print(f"[-] Warning: No CSV data found for {backend}")
        print("--- RAW OUTPUT START ---")
        print(output)
        print("--- RAW OUTPUT END ---")
        return []
        
    df = pd.read_csv(io.StringIO('\n'.join(csv_lines)))
    records = []
    
    for _, row in df.iterrows():
        kname = str(row.get("Kernel Name", ""))
        mname = str(row.get("Metric Name", ""))
        mval = str(row.get("Metric Value", "0"))
        
        if "adt_score" in kname or "calc_energy" in kname:
            k_clean = "adt_score"
        elif "genetic" in kname or "iterate" in kname or "finalize" in kname:
            k_clean = "genetic"
        elif "geom_transform" in kname or "apply" in kname:
            k_clean = "geom_transform"
        else:
            continue
            
        if mname == "sm__warps_active.avg.pct_of_peak_sustained_active":
            m_clean = "Achieved Occupancy (%)"
        elif mname == "sm__throughput.avg.pct_of_peak_sustained_elapsed":
            m_clean = "SM Compute Throughput (%)"
        elif mname == "gpu__compute_memory_throughput.avg.pct_of_peak_sustained_elapsed":
            m_clean = "Memory Bandwidth (%)"
        elif mname == "launch__registers_per_thread":
            m_clean = "Registers / Thread"
        else:
            continue
            
        try:
            val = float(mval.replace(',', ''))
        except:
            val = 0.0
            
        records.append({
            "Backend": backend,
            "Kernel": k_clean,
            "Metric": m_clean,
            "Value": val
        })
        
    return records

def main():
    print("=== NCU Runtime Hardware Profiler ===")
    
    out_cuda = run_ncu(CUDA_BIN, "CUDA Native", "CUDA:GPU:0")
    rec_cuda = parse_ncu_csv(out_cuda, "CUDA Native")
    
    out_alpaka = run_ncu(ALPAKA_BIN, "Alpaka CUDA", "ALPAKA:GPU:0")
    rec_alpaka = parse_ncu_csv(out_alpaka, "Alpaka CUDA")
    
    df = pd.DataFrame(rec_cuda + rec_alpaka)
    
    if not df.empty:
        # Average multiple invocations of the same kernel
        df = df.groupby(['Backend', 'Kernel', 'Metric']).mean().reset_index()
        
        csv_path = "ncu_metrics.csv"
        df.to_csv(csv_path, index=False)
        print(f"\n[+] Saved metrics to {csv_path}")
        print("\nExtracted Data:")
        print(df.to_string(index=False))
        
        # Plotting
        sns.set_theme(style="whitegrid")
        metrics_list = df['Metric'].unique()
        num_metrics = len(metrics_list)
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        for i, metric in enumerate(metrics_list):
            if i < 4:
                subset = df[df['Metric'] == metric]
                sns.barplot(data=subset, x="Kernel", y="Value", hue="Backend", ax=axes[i], palette="viridis")
                axes[i].set_title(f"{metric}")
                axes[i].set_ylabel("Value")
                # Do not clamp y-axis so small percentages are visible
                # Add value labels to the top of each bar
                for container in axes[i].containers:
                    axes[i].bar_label(container, fmt='%.2f', padding=3)
        plt.tight_layout()
        pdf_path = "ncu_analysis.pdf"
        plt.savefig(pdf_path, format="pdf")
        print(f"[+] Saved multi-plot to {pdf_path}")
        
    else:
        print("[-] Error: No data could be parsed. Did ncu run successfully?")

if __name__ == "__main__":
    main()
