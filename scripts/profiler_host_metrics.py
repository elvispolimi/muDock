import subprocess
import re
import csv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ─── CONFIG ──────────────────────────────────────────────────────────────────
PROTEIN = "/work/onedina/muDock_ON/data/1fkb/1fkb_pocket.pdbqt"
LIGAND = "/work/onedina/muDock_ON/data/1fkb/1fkb_ligand.adtmol2"

BACKENDS = [
    ("CUDA",   "/work/onedina/muDock_ON/build/cuda/application/muDock",        "CUDA:GPU:0"),
    ("Alpaka", "/work/onedina/muDock_ON/build/alpaka-cuda/application/muDock",  "ALPAKA:GPU:0"),
]

# We test how RAM and CPU scale with the population size (similar to the concurrent events in the paper)
POPULATIONS  = [50, 100, 200, 500, 1000]
GENERATIONS  = 1000
RUNS         = 3

def extract_host_metrics(time_output):
    """
    Parses the output of `/usr/bin/time -v` and returns (RSS_MB, CPU_Percent).
    """
    rss = None
    cpu = None
    
    for line in time_output.split('\n'):
        if "Maximum resident set size (kbytes):" in line:
            kb = int(re.search(r"(\d+)", line).group(1))
            rss = kb / 1024.0 # Convert to MB
        elif "Percent of CPU this job got:" in line:
            m = re.search(r"(\d+)%", line)
            if m:
                cpu = float(m.group(1))
                
    return rss, cpu

def run_host_profiling():
    results = []
    
    print("=== Host Metrics Profiler (RAM & CPU) ===")
    
    for name, bin_path, use_flag in BACKENDS:
        if not os.path.exists(bin_path):
            print(f"[-] Executable not found: {bin_path}")
            continue
            
        print(f"\n▶ Profiling Backend: {name}")
        for pop in POPULATIONS:
            cmd = [
                "/usr/bin/time", "-v",
                bin_path,
                "--protein", PROTEIN,
                "--ligand", LIGAND,
                "--use", use_flag,
                "--population", str(pop),
                "--generations", str(GENERATIONS)
            ]
            
            rss_runs = []
            cpu_runs = []
            
            print(f"  Pop: {pop:4d} | ", end="", flush=True)
            for r in range(RUNS):
                res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                # /usr/bin/time -v outputs to stderr
                rss, cpu = extract_host_metrics(res.stderr)
                
                if rss is not None and cpu is not None:
                    rss_runs.append(rss)
                    cpu_runs.append(cpu)
                    print(".", end="", flush=True)
                else:
                    print("x", end="", flush=True)
                    
            if rss_runs and cpu_runs:
                avg_rss = sum(rss_runs) / len(rss_runs)
                avg_cpu = sum(cpu_runs) / len(cpu_runs)
                print(f" | Avg RSS: {avg_rss:6.1f} MB, Avg CPU: {avg_cpu:5.1f}%")
                
                results.append({
                    "Backend": name,
                    "Population": pop,
                    "Peak RSS (MB)": avg_rss,
                    "CPU Utilization (%)": avg_cpu
                })

    if not results:
        print("[-] No data collected.")
        return

    # Save to CSV
    df = pd.DataFrame(results)
    csv_path = "host_metrics.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n[+] Saved metrics to {csv_path}")

    # Plotting
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Peak RSS Plot
    sns.lineplot(data=df, x="Population", y="Peak RSS (MB)", hue="Backend", marker="o", ax=axes[0])
    axes[0].set_title("Peak Host RSS (RAM Usage)")
    axes[0].set_ylabel("Memory (MB)")
    axes[0].set_xlabel("Population Size")

    # CPU Utilization Plot
    sns.lineplot(data=df, x="Population", y="CPU Utilization (%)", hue="Backend", marker="o", ax=axes[1])
    axes[1].set_title("Host CPU Utilization")
    axes[1].set_ylabel("CPU Utilization (%)")
    axes[1].set_xlabel("Population Size")
    
    # Set y-axis to start from 0 for better perspective
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylim(bottom=0, top=105)

    plt.tight_layout()
    pdf_path = "host_metrics_analysis.pdf"
    plt.savefig(pdf_path, format="pdf")
    print(f"[+] Saved plots to {pdf_path}")

if __name__ == "__main__":
    run_host_profiling()
