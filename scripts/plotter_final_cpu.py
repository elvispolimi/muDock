import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

RESULTS_DIR = "/work/onedina/muDock_ON/scripts"

def plot_cpu_results(csv_filename, output_filename, title, backend_filter=None, y_metric="Throughput (Evals/s)", y_label="Throughput (Evals/s)", format_int=True):
    csv_path = os.path.join(RESULTS_DIR, csv_filename)
    if not os.path.exists(csv_path):
        print(f"[-] {csv_filename} not found.")
        return

    df = pd.read_csv(csv_path)
    if df.empty:
        print(f"[-] {csv_filename} is empty.")
        return
        
    if backend_filter:
        df = df[df['Backend'].isin(backend_filter)]
        if df.empty:
            print(f"[-] No data for filtered backends in {csv_filename}.")
            return

    # Average the 3 runs for plotting
    df_avg = df.groupby(['Backend', 'Population', 'Generations'])[y_metric].mean().reset_index()

    sns.set_theme(style="whitegrid")
    
    # Create a grid of subplots: one for each Population
    populations = sorted(df_avg['Population'].unique())
    num_pops = len(populations)
    
    # Adjust width based on number of generations to fit nicely
    fig, axes = plt.subplots(1, num_pops, figsize=(5 * num_pops, 6), sharey=True)
    
    if num_pops == 1:
        axes = [axes]

    for i, pop in enumerate(populations):
        subset = df_avg[df_avg['Population'] == pop]
        sns.barplot(
            data=subset, 
            x="Generations", 
            y=y_metric, 
            hue="Backend", 
            ax=axes[i], 
            palette="Spectral"
        )
        axes[i].set_title(f"Population: {pop}")
        if i == 0:
            axes[i].set_ylabel(y_label)
        else:
            axes[i].set_ylabel("")
            
        # Add values on top of bars, rotate text slightly so it doesn't overlap
        for p in axes[i].patches:
            height = p.get_height()
            if height > 0:
                val_str = f'{int(height)}' if format_int else f'{height:.2f}'
                axes[i].annotate(val_str, 
                                 (p.get_x() + p.get_width() / 2., height), 
                                 ha='center', va='bottom', 
                                 fontsize=8, color='black', xytext=(0, 2), 
                                 textcoords='offset points', rotation=45)

    plt.suptitle(title, fontsize=16, y=1.05)
    plt.tight_layout()
    
    out_path = os.path.join(RESULTS_DIR, output_filename)
    plt.savefig(out_path, format="pdf", bbox_inches='tight')
    print(f"[+] Saved {output_filename}")
    plt.close()

if __name__ == "__main__":
    print("Generating Final CPU Macro Benchmark Plots (Separated by Serial / OMP)...")
    
    serial_backends = ['CPP Serial', 'Alpaka Serial (Unroll)', 'Alpaka Serial (No Unroll)']
    omp_backends = ['CPP OMP', 'Alpaka OMP (Unroll)', 'Alpaka OMP (No Unroll)']
    
    # --- SERIAL ---
    plot_cpu_results(
        "macro_results_cpu_single.csv", 
        "macro_plot_cpu_serial_throughput.pdf", 
        "CPU Serial Throughput (Dataset: Single)",
        backend_filter=serial_backends,
        y_metric="Throughput (Evals/s)",
        y_label="Throughput (Evals/s)",
        format_int=True
    )
    plot_cpu_results(
        "macro_results_cpu_single.csv", 
        "macro_plot_cpu_serial_latency.pdf", 
        "CPU Serial Latency (Dataset: Single)",
        backend_filter=serial_backends,
        y_metric="Time (s)",
        y_label="Time (Seconds)",
        format_int=False
    )
    
    # --- OMP ---
    plot_cpu_results(
        "macro_results_cpu_single.csv", 
        "macro_plot_cpu_omp_throughput.pdf", 
        "CPU OMP Throughput (Dataset: Single)",
        backend_filter=omp_backends,
        y_metric="Throughput (Evals/s)",
        y_label="Throughput (Evals/s)",
        format_int=True
    )
    plot_cpu_results(
        "macro_results_cpu_single.csv", 
        "macro_plot_cpu_omp_latency.pdf", 
        "CPU OMP Latency (Dataset: Single)",
        backend_filter=omp_backends,
        y_metric="Time (s)",
        y_label="Time (Seconds)",
        format_int=False
    )
