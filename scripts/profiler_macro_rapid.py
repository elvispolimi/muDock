import subprocess
import re
import csv
from pathlib import Path
import pandas as pd

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

# Configurazione Rapida (Solo 1 carico centrale)
CONFIGS = [
    (100, 1000)
]

WARMUP_RUNS = 1
RUNS = 3
TIMEOUT_SEC = 300  # 5 minuti


def extract_total_time(output_text):
    match = re.search(r'\[\s*([\d\.]+)\s*\]\s*INFO All Done!', output_text)
    return float(match.group(1)) if match else None

def run_rapid_benchmark():
    for ds_name, paths in DATASETS.items():
        results = []
        csv_output = f"rapid_macro_results_{ds_name}.csv"
        stats_output = f"rapid_macro_stats_{ds_name}.csv"
        
        print("="*60)
        print(f" RAPID MACRO PROFILING (3 CONFIGS) - {ds_name.upper()}")
        print("="*60)

        for name, bin_path, use_flag in BACKENDS:
            if not Path(bin_path).exists():
                continue
                
            for pop, gen in CONFIGS:
                cmd = [bin_path, "--protein", paths["PROTEIN"], "--ligand", paths["LIGAND"], "--use", use_flag, "--population", str(pop), "--generations", str(gen)]
                
                # WARMUP
                print(f"[>] {name} | Pop: {pop}, Gen: {gen} | Warmup...")
                try:
                    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIMEOUT_SEC)
                except:
                    pass
                
                # RUNS
                for r in range(RUNS):
                    print(f"    -> Run {r+1}/{RUNS}...")
                    try:
                        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=TIMEOUT_SEC)
                        time_sec = extract_total_time(res.stdout)
                        if time_sec:
                            results.append({
                                "Backend": name, "Population": pop, "Generations": gen, "Run": r+1,
                                "Time (s)": time_sec, "Throughput (Evals/s)": (pop*gen)/time_sec
                            })
                    except:
                        pass

        if results:
            df = pd.DataFrame(results)
            df.to_csv(csv_output, index=False)
            stats = df.groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"].agg(['mean', 'var']).reset_index()
            print(f"\n--- RAPID STATS ({ds_name}) ---")
            print(stats.to_string(index=False))
            stats.to_csv(stats_output, index=False)

if __name__ == "__main__":
    run_rapid_benchmark()
