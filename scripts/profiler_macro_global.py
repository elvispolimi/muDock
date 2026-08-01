import subprocess
import re
import csv
import itertools
from pathlib import Path

DATASETS = {
    "single": {
        "PROTEIN": "/work/onedina/muDock/data/1fkb/1fkb_pocket.pdbqt",
        "LIGAND": "/work/onedina/muDock/data/1fkb/1fkb_ligand.adtmol2"
    },
    "small_multi": {
        "PROTEIN": "/work/onedina/muDock/data/1fkb/1fkb_pocket.pdbqt",
        "LIGAND": "/work/onedina/muDock/data/small 1/small.adtmol2"
    }
}

BACKENDS = [
    ("CUDA", "/work/onedina/muDock/build/cuda/application/muDock", "CUDA:GPU:0"),
    ("Alpaka", "/work/onedina/muDock/build/alpaka-cuda/application/muDock", "ALPAKA:GPU:0")
]

POPULATIONS = [50, 100, 200, 500]
GENERATIONS = [500, 1000, 2000, 5000]

WARMUP_RUNS = 1
RUNS = 3
TIMEOUT_SEC = 300  # 5 minuti massimi per run (evita stalli letali)



def extract_total_time(output_text):
    match = re.search(r'\[\s*([\d\.]+)\s*\]\s*INFO All Done!', output_text)
    return float(match.group(1)) if match else None

def run_macro_benchmark():
    total_tasks_per_ds = len(BACKENDS) * len(POPULATIONS) * len(GENERATIONS) * RUNS

    for ds_name, paths in DATASETS.items():
        results = []
        current_task = 0
        csv_output = f"profiling_macro_results_{ds_name}.csv"
        
        print("="*70)
        print(f" MACRO PROFILING START (GLOBAL VISION HPC) - Dataset: {ds_name.upper()}")
        print(f" Protocol: {WARMUP_RUNS} Warmups + {RUNS} Runs | Timeout: {TIMEOUT_SEC}s")
        print("="*70)

        fieldnames = ["Backend", "Population", "Generations", "Run", "Total Evaluations", "Time (s)", "Throughput (Evals/s)"]
        with open(csv_output, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

        for name, bin_path, use_flag in BACKENDS:
            if not Path(bin_path).exists():
                print(f"[!] WARNING: Executable not found -> {bin_path}. Skipping {name}.")
                continue
                
            for pop, gen in itertools.product(POPULATIONS, GENERATIONS):
                
                # --- WARM UP ---
                print(f"\n[>] {name} ({ds_name}) | Pop: {pop:<4} | Gen: {gen:<5} -> Executing {WARMUP_RUNS} Warmup runs...")
                cmd = [
                    bin_path, "--protein", paths["PROTEIN"], "--ligand", paths["LIGAND"],
                    "--use", use_flag, "--population", str(pop), "--generations", str(gen)
                ]
                
                skip_config = False
                for w in range(WARMUP_RUNS):
                    try:
                        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIMEOUT_SEC)
                    except subprocess.TimeoutExpired:
                        print(f"    [!] Warmup TIMED OUT (>10 mins). Skipping configuration Pop={pop}, Gen={gen}.")
                        skip_config = True
                        break
                
                if skip_config:
                    current_task += RUNS
                    continue

                # --- MISURAZIONI ---
                for run_idx in range(RUNS):
                    current_task += 1
                    progress = (current_task / total_tasks_per_ds) * 100
                    print(f"[{progress:5.1f}%] {name} | Pop: {pop:<4} | Gen: {gen:<5} | Run: {run_idx+1}/{RUNS}")
                    
                    try:
                        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True, timeout=TIMEOUT_SEC)
                        time_sec = extract_total_time(res.stdout)
                        
                        if time_sec is not None:
                            evaluations = pop * gen
                            throughput_evals_sec = evaluations / time_sec if time_sec > 0 else 0

                            row_data = {
                                "Backend": name,
                                "Population": pop,
                                "Generations": gen,
                                "Run": run_idx + 1,
                                "Total Evaluations": evaluations,
                                "Time (s)": time_sec,
                                "Throughput (Evals/s)": throughput_evals_sec
                            }
                            results.append(row_data)
                            
                            with open(csv_output, mode='a', newline='') as f:
                                writer = csv.DictWriter(f, fieldnames=fieldnames)
                                writer.writerow(row_data)
                        else:
                            print(f"    [!] ERROR: Execution time not found for {name}.")
                    except subprocess.TimeoutExpired:
                        print("    [!] Run TIMED OUT.")
                    except subprocess.CalledProcessError as e:
                        print(f"    [!] CRITICAL EXECUTION ERROR: {e}")
                        print(e.stdout)

        if results:
            print(f"\n[+] MACRO profiling completed for {ds_name}. Results saved to: '{csv_output}'.")
        else:
            print(f"\n[!] No valid results obtained for {ds_name}.")

if __name__ == "__main__":
    run_macro_benchmark()
