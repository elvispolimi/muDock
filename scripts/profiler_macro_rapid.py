"""
profiler_macro_rapid.py
Rapid macro-profiling: 1 configuration (Pop=100, Gen=1000), 1 warmup + 3 runs.
Covers both 'single' and 'small_multi' datasets.
Produces: macro_results_rapid_{dataset}.csv  +  macro_stats_rapid_{dataset}.csv
"""

import subprocess
import re
import csv
import pandas as pd
from pathlib import Path
import sys

# ─── CONFIG ──────────────────────────────────────────────────────────────────
DATASETS = {
    "single": {
        "PROTEIN": "/work/onedina/muDock_ON/data/1fkb/1fkb_pocket.pdbqt",
        "LIGAND":  "/work/onedina/muDock_ON/data/1fkb/1fkb_ligand.adtmol2",
    },
    "small_multi": {
        "PROTEIN": "/work/onedina/muDock_ON/data/1fkb/1fkb_pocket.pdbqt",
        "LIGAND":  "/work/onedina/muDock_ON/data/small_1/small.adtmol2",
    },
}

BACKENDS = [
    ("CUDA",   "/work/onedina/muDock_ON/build/cuda/application/muDock",        "CUDA:GPU:0"),
    ("Alpaka", "/work/onedina/muDock_ON/build/alpaka-cuda/application/muDock",  "ALPAKA:GPU:0"),
]

CONFIGS      = [(100, 1000)]   # (population, generations)
WARMUP_RUNS  = 1
RUNS         = 3
TIMEOUT_SEC  = 300

FIELDNAMES   = ["Backend", "Population", "Generations", "Run", "Time (s)", "Throughput (Evals/s)"]

# ─── HELPERS ─────────────────────────────────────────────────────────────────
def banner(title: str, width: int = 62) -> None:
    print(f"\n╔{'═' * (width - 2)}╗")
    print(f"║  {title:<{width - 4}}║")
    print(f"╚{'═' * (width - 2)}╝")

def section(msg: str) -> None:
    print(f"\n  ▶  {msg}")

def ok(msg: str) -> None:
    print(f"  ✔  {msg}")

def warn(msg: str) -> None:
    print(f"  ⚠  {msg}", file=sys.stderr)

def extract_total_time(text: str) -> float | None:
    m = re.search(r'\[\s*([\d\.]+)\s*\]\s*INFO All Done!', text)
    return float(m.group(1)) if m else None

def print_run_row(run_idx: int, total_runs: int, name: str, pop: int,
                  gen: int, time_s: float, throughput: float) -> None:
    print(f"  [{run_idx:>2}/{total_runs}]  {name:<12}  Pop:{pop:<5}  Gen:{gen:<6}"
          f"  {time_s:>8.3f}s   {throughput:>10,.0f} Evals/s")

def print_stats_table(df: pd.DataFrame, ds_name: str) -> None:
    stats = (df
             .groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"]
             .agg(Mean="mean", Std="std", Min="min", Max="max")
             .reset_index())
    col_w = 14
    header = (f"  {'Backend':<{col_w}}  {'Pop':>6}  {'Gen':>6}"
              f"  {'Mean Evals/s':>13}  {'Std':>11}  {'Min':>11}  {'Max':>11}")
    sep    = "  " + "─" * (len(header) - 2)
    print(f"\n  ┌─ RAPID STATS: {ds_name.upper()} {'─' * max(0, 50 - len(ds_name))}┐")
    print(header)
    print(sep)
    for _, r in stats.iterrows():
        std_str = f"{r['Std']:>11,.0f}" if pd.notna(r["Std"]) else f"{'—':>11}"
        print(f"  {r['Backend']:<{col_w}}  {int(r['Population']):>6}  {int(r['Generations']):>6}"
              f"  {r['Mean']:>13,.0f}  {std_str}  {r['Min']:>11,.0f}  {r['Max']:>11,.0f}")
    print(f"  └{'─' * (len(header) - 2)}┘")


# ─── MAIN ────────────────────────────────────────────────────────────────────
def run_rapid_benchmark() -> None:
    banner("MACRO RAPID PROFILING  —  muDock  |  1 Warmup + 3 Runs")

    for ds_name, paths in DATASETS.items():
        section(f"Dataset: {ds_name.upper()}")
        results: list[dict] = []
        csv_output   = f"macro_results_rapid_{ds_name}.csv"
        stats_output = f"macro_stats_rapid_{ds_name}.csv"

        with open(csv_output, mode="w", newline="") as f:
            csv.DictWriter(f, fieldnames=FIELDNAMES).writeheader()

        for name, bin_path, use_flag in BACKENDS:
            if not Path(bin_path).exists():
                warn(f"Executable not found → {bin_path}  (skipping {name})")
                continue

            for pop, gen in CONFIGS:
                cmd = [
                    bin_path,
                    "--protein",     paths["PROTEIN"],
                    "--ligand",      paths["LIGAND"],
                    "--use",         use_flag,
                    "--population",  str(pop),
                    "--generations", str(gen),
                ]

                # ── Warmup ──
                print(f"\n  [warm-up]  {name}  Pop:{pop}  Gen:{gen} …", flush=True)
                try:
                    subprocess.run(cmd, stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL, timeout=TIMEOUT_SEC)
                except Exception:
                    pass

                # ── Measurements ──
                for r_idx in range(RUNS):
                    try:
                        res = subprocess.run(
                            cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            text=True,
                            timeout=TIMEOUT_SEC,
                        )
                        t = extract_total_time(res.stdout)
                        if t is None:
                            warn(f"Time not found in output (run {r_idx + 1})")
                            continue

                        evals_s = (pop * gen) / t
                        row = {
                            "Backend":             name,
                            "Population":          pop,
                            "Generations":         gen,
                            "Run":                 r_idx + 1,
                            "Time (s)":            t,
                            "Throughput (Evals/s)": evals_s,
                        }
                        results.append(row)
                        print_run_row(r_idx + 1, RUNS, name, pop, gen, t, evals_s)

                        with open(csv_output, mode="a", newline="") as f:
                            csv.DictWriter(f, fieldnames=FIELDNAMES).writerow(row)

                    except subprocess.TimeoutExpired:
                        warn(f"Run {r_idx + 1} timed out (>{TIMEOUT_SEC}s)")
                    except subprocess.CalledProcessError as e:
                        warn(f"Run {r_idx + 1} error: {e}")

        if results:
            df = pd.DataFrame(results)
            print_stats_table(df, ds_name)
            stats = (df
                     .groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"]
                     .agg(Mean="mean", Std="std", Min="min", Max="max")
                     .reset_index())
            stats.to_csv(stats_output, index=False)
            ok(f"Results → {csv_output}")
            ok(f"Stats   → {stats_output}")
        else:
            warn(f"No valid results for dataset {ds_name}")


if __name__ == "__main__":
    run_rapid_benchmark()
