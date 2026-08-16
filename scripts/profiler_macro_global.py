"""
profiler_macro_global.py
Full macro-profiling: 4x4 grid (Pop × Gen), 1 warmup + 3 runs.
Covers both 'single' and 'small_multi' datasets.
Produces: macro_results_global_{dataset}.csv
"""

import subprocess
import re
import csv
import itertools
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

POPULATIONS  = [50, 100, 200, 500]
GENERATIONS  = [500, 1000, 2000, 5000]
WARMUP_RUNS  = 1
RUNS         = 3
TIMEOUT_SEC  = 300   # timeout per single run (seconds)

FIELDNAMES   = ["Backend", "Population", "Generations", "Run",
                "Total Evaluations", "Time (s)", "Throughput (Evals/s)"]

# ─── HELPERS ─────────────────────────────────────────────────────────────────
def banner(title: str, width: int = 68) -> None:
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

def progress_bar(current: int, total: int, width: int = 30) -> str:
    filled = int(width * current / total)
    bar    = "█" * filled + "░" * (width - filled)
    pct    = 100.0 * current / total
    return f"[{bar}] {pct:5.1f}%  ({current}/{total})"


# ─── MAIN ────────────────────────────────────────────────────────────────────
def run_macro_benchmark() -> None:
    total_configs = len(BACKENDS) * len(POPULATIONS) * len(GENERATIONS)
    total_tasks   = total_configs * RUNS

    banner(
        f"MACRO GLOBAL PROFILING  —  muDock  |  "
        f"{len(POPULATIONS)}×{len(GENERATIONS)} grid  "
        f"|  {WARMUP_RUNS} Warmup + {RUNS} Runs"
    )

    for ds_name, paths in DATASETS.items():
        section(f"Dataset: {ds_name.upper()}")
        print(f"     Backends     : {', '.join(n for n,*_ in BACKENDS)}")
        print(f"     Populations  : {POPULATIONS}")
        print(f"     Generations  : {GENERATIONS}")
        print(f"     Total runs   : {total_tasks} per backend")

        results: list[dict] = []
        csv_output  = f"macro_results_global_{ds_name}.csv"
        task_done   = 0

        with open(csv_output, mode="w", newline="") as f:
            csv.DictWriter(f, fieldnames=FIELDNAMES).writeheader()

        for name, bin_path, use_flag in BACKENDS:
            if not Path(bin_path).exists():
                warn(f"Executable not found → {bin_path}  (skipping {name})")
                task_done += len(POPULATIONS) * len(GENERATIONS) * RUNS
                continue

            print()
            for pop, gen in itertools.product(POPULATIONS, GENERATIONS):
                cmd = [
                    bin_path,
                    "--protein",     paths["PROTEIN"],
                    "--ligand",      paths["LIGAND"],
                    "--use",         use_flag,
                    "--population",  str(pop),
                    "--generations", str(gen),
                ]

                # ── Warmup ──
                skip = False
                for _ in range(WARMUP_RUNS):
                    try:
                        subprocess.run(cmd, stdout=subprocess.DEVNULL,
                                       stderr=subprocess.DEVNULL, timeout=TIMEOUT_SEC)
                    except subprocess.TimeoutExpired:
                        warn(f"Warmup timed out → {name}  Pop:{pop}  Gen:{gen}  — config skipped")
                        skip = True
                        break

                if skip:
                    task_done += RUNS
                    continue

                # ── Measurements ──
                run_times: list[float] = []
                for r_idx in range(RUNS):
                    task_done += 1
                    pb = progress_bar(task_done, total_tasks)
                    print(f"  {pb}  {name:<6}  Pop:{pop:<5}  Gen:{gen:<6}  Run:{r_idx+1}/{RUNS}",
                          end="  ", flush=True)
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
                            print("→ ⚠ time not found")
                            continue

                        evals_s = (pop * gen) / t
                        run_times.append(t)
                        print(f"→ {t:.3f}s  ({evals_s:,.0f} Evals/s)")

                        row = {
                            "Backend":              name,
                            "Population":           pop,
                            "Generations":          gen,
                            "Run":                  r_idx + 1,
                            "Total Evaluations":    pop * gen,
                            "Time (s)":             t,
                            "Throughput (Evals/s)": evals_s,
                        }
                        results.append(row)
                        with open(csv_output, mode="a", newline="") as f:
                            csv.DictWriter(f, fieldnames=FIELDNAMES).writerow(row)

                    except subprocess.TimeoutExpired:
                        print(f"→ ⚠ timed out (>{TIMEOUT_SEC}s)")
                    except subprocess.CalledProcessError as e:
                        print(f"→ ✗ error: {e}")

                if run_times:
                    avg_t    = sum(run_times) / len(run_times)
                    avg_eval = (pop * gen) / avg_t
                    print(f"         └─ avg  {avg_t:.3f}s   {avg_eval:,.0f} Evals/s")

        if results:
            ok(f"Global profiling done for '{ds_name}'.")
            ok(f"Results → {csv_output}  ({len(results)} rows)")
        else:
            warn(f"No valid results for dataset '{ds_name}'.")


if __name__ == "__main__":
    run_macro_benchmark()
