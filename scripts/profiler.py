import subprocess
import re
import csv
import itertools
from pathlib import Path

PROTEIN = "/work/onedina/muDock/data/1fkb/1fkb_pocket.pdbqt"
LIGAND = "/work/onedina/muDock/data/1fkb/1fkb_ligand.adtmol2"

BACKENDS = [
    ("CUDA Native", "/work/onedina/muDock/build/cuda/application/muDock", "CUDA:GPU:0"),
    ("Alpaka CUDA", "/work/onedina/muDock/build/alpaka-cuda/application/muDock", "ALPAKA:GPU:0")
]

POPULATIONS = [50, 100, 200, 500]
GENERATIONS = [500, 1000, 2000, 5000]
RUNS = 5

CSV_OUTPUT = "profiling_macro_results.csv"

def extract_total_time(output_text):
    match = re.search(r'\[\s*([\d\.]+)\s*\]\s*INFO All Done!', output_text)
    return float(match.group(1)) if match else None

def run_macro_benchmark():
    results = []
    total_tasks = len(BACKENDS) * len(POPULATIONS) * len(GENERATIONS) * RUNS
    current_task = 0

    print("="*60)
    print(" MACRO PROFILING START (GLOBAL VISION HPC)")
    print("="*60)

    fieldnames = ["Backend", "Population", "Generations", "Run", "Total Evaluations", "Time (s)", "Throughput (Evals/s)"]
    with open(CSV_OUTPUT, mode='w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

    for name, bin_path, use_flag in BACKENDS:
        if not Path(bin_path).exists():
            print(f"[!] WARNING: Executable not found -> {bin_path}. Skipping {name}.")
            continue
            
        for pop, gen in itertools.product(POPULATIONS, GENERATIONS):
            for run_idx in range(RUNS):
                current_task += 1
                progress = (current_task / total_tasks) * 100
                print(f"[{progress:5.1f}%] {name} | Pop: {pop:<4} | Gen: {gen:<5} | Run: {run_idx+1}/{RUNS}")
                
                cmd = [
                    bin_path,
                    "--protein", PROTEIN,
                    "--ligand", LIGAND,
                    "--use", use_flag,
                    "--population", str(pop),
                    "--generations", str(gen)
                ]
                
                try:
                    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
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
                        
                        with open(CSV_OUTPUT, mode='a', newline='') as f:
                            writer = csv.DictWriter(f, fieldnames=fieldnames)
                            writer.writerow(row_data)
                    else:
                        print(f"[!] ERROR: Execution time not found for {name}.")
                except subprocess.CalledProcessError as e:
                    print(f"[!] CRITICAL EXECUTION ERROR: {e}")
                    print(e.stdout)

    if results:
        print(f"\n[+] MACRO profiling completed. Results saved to: '{CSV_OUTPUT}'.")
    else:
        print("\n[!] No valid results obtained.")

if __name__ == "__main__":
    run_macro_benchmark()
