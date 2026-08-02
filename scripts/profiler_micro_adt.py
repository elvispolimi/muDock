"""
profiler_micro_adt.py
Micro-profiling of the adt_score (calc_energy) kernel via NVIDIA Nsight Systems.
Produces: micro_adt_score.pdf
"""

import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path
import sys

# ─── CONFIG ──────────────────────────────────────────────────────────────────
PROTEIN    = "/work/onedina/muDock_ON/data/1fkb/1fkb_pocket.pdbqt"
LIGAND     = "/work/onedina/muDock_ON/data/1fkb/1fkb_ligand.adtmol2"
POPULATION = 100
GENERATIONS = 1000

BACKENDS = [
    ("CUDA Native",  "/work/onedina/muDock_ON/build/cuda/application/muDock",        "CUDA:GPU:0"),
    ("Alpaka CUDA",  "/work/onedina/muDock_ON/build/alpaka-cuda/application/muDock",  "ALPAKA:GPU:0"),
]

KERNEL_FILTER = r"calc_energy|adt_score"
KERNEL_LABEL  = "adt_score"
OUTPUT_PDF    = "micro_adt_score.pdf"

# Palette coerente (uguale in tutti i micro script)
PALETTE = {"CUDA Native": "#E74C3C", "Alpaka CUDA": "#3498DB"}

# ─── PLOT STYLE ──────────────────────────────────────────────────────────────
sns.set_theme(style="ticks", context="paper")
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.labelsize":    13,
    "axes.titlesize":    13,
    "axes.titleweight":  "bold",
    "axes.grid":         True,
    "grid.alpha":        0.35,
    "grid.linestyle":    "--",
    "legend.fontsize":   10,
    "figure.titlesize":  15,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

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

def find_col(cols, *keywords) -> str:
    for kw in keywords:
        for c in cols:
            if kw.lower() in c.lower():
                return c
    raise KeyError(f"Column not found: {keywords}")


def run_nsys(name: str, bin_path: str, use_flag: str) -> pd.DataFrame | None:
    safe_name   = name.replace(" ", "_")
    report_base = f"nsys_adt_{safe_name}"
    sqlite_file = f"{report_base}.sqlite"

    section(f"nsys profile  →  {name}")

    cmd_profile = [
        "nsys", "profile",
        "--force-overwrite=true",
        "--stats=true",
        f"--output={report_base}",
        bin_path,
        "--protein",     PROTEIN,
        "--ligand",      LIGAND,
        "--use",         use_flag,
        "--population",  str(POPULATION),
        "--generations", str(GENERATIONS),
        "--observer",    "1",
    ]
    try:
        subprocess.run(cmd_profile, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError as e:
        warn(f"nsys profile failed for {name}: {e}")
        return None

    cmd_stats = [
        "nsys", "stats",
        "--report", "gpukernsum",
        "--format", "csv",
        "--force-overwrite", "true",
        "--output", report_base,
        sqlite_file,
    ]
    subprocess.run(cmd_stats, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, check=True)

    csv_kern = f"{report_base}_gpukernsum.csv"
    if not Path(csv_kern).exists():
        warn(f"CSV not found: {csv_kern}")
        return None

    df_k  = pd.read_csv(csv_kern)
    df_kf = df_k[df_k["Name"].str.contains(KERNEL_FILTER, case=False, na=False)].copy()
    if df_kf.empty:
        warn(f"No {KERNEL_LABEL} kernels found for {name}")
        return None

    df_kf["Kernel"] = KERNEL_LABEL
    time_col = find_col(df_kf.columns, "Total Time", "Time (ns)")
    inst_col = find_col(df_kf.columns, "Instances", "Count")
    avg_col  = find_col(df_kf.columns, "Avg", "Average")
    min_col  = find_col(df_kf.columns, "Min")
    max_col  = find_col(df_kf.columns, "Max")

    agg = df_kf.groupby("Kernel").agg({
        time_col: "sum",
        inst_col: "sum",
        avg_col:  "mean",
        min_col:  "min",
        max_col:  "max",
    }).reset_index()

    agg["Total Time (ms)"] = agg[time_col] / 1e6
    agg["Avg Time (µs)"]   = agg[avg_col]  / 1e3
    agg["Min Time (µs)"]   = agg[min_col]  / 1e3
    agg["Max Time (µs)"]   = agg[max_col]  / 1e3
    agg["Backend"]         = name

    total_ms = agg["Total Time (ms)"].iloc[0]
    avg_us   = agg["Avg Time (µs)"].iloc[0]
    min_us   = agg["Min Time (µs)"].iloc[0]
    max_us   = agg["Max Time (µs)"].iloc[0]
    launches = int(agg[inst_col].iloc[0])

    ok(f"{name}")
    print(f"     {'Total':>10}  {'Avg':>10}  {'Min':>10}  {'Max':>10}  {'Launches':>10}")
    print(f"     {total_ms:>9.2f}ms  {avg_us:>9.2f}µs  {min_us:>9.2f}µs  {max_us:>9.2f}µs  {launches:>10,}")

    return agg


def make_barplot(ax, df, y_col, title, ylabel, fmt=".2f"):
    colors = [PALETTE.get(b, "#888") for b in df["Backend"]]
    bars = ax.bar(df["Backend"], df[y_col], color=colors,
                  edgecolor="white", linewidth=0.8, width=0.5)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    for bar, val in zip(bars, df[y_col]):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() * 1.02,
                f"{val:{fmt}}", ha="center", va="bottom",
                fontsize=10, fontweight="bold")
    sns.despine(ax=ax)


# ─── MAIN ────────────────────────────────────────────────────────────────────
def main():
    banner(f"MICRO PROFILING — adt_score kernel  |  Pop:{POPULATION}  Gen:{GENERATIONS}")

    frames = []
    for name, bin_path, use_flag in BACKENDS:
        if not Path(bin_path).exists():
            warn(f"Executable not found → {bin_path}  (skipping {name})")
            continue
        result = run_nsys(name, bin_path, use_flag)
        if result is not None:
            frames.append(result)

    if not frames:
        warn("No data collected. Exiting.")
        sys.exit(1)

    df = pd.concat(frames, ignore_index=True)

    section("Generating PDF plots…")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"Micro-Profiling: adt_score  —  Pop:{POPULATION}  Gen:{GENERATIONS}",
        fontsize=15, fontweight="bold", y=1.01,
    )

    make_barplot(axes[0, 0], df, "Total Time (ms)", "Total GPU Time",         "Time (ms)")
    make_barplot(axes[0, 1], df, "Avg Time (µs)",   "Average Time per Launch", "Time (µs)")
    make_barplot(axes[1, 0], df, "Min Time (µs)",   "Min Launch Latency",      "Time (µs)")
    make_barplot(axes[1, 1], df, "Max Time (µs)",   "Max Launch Latency",      "Time (µs)")

    plt.tight_layout()
    fig.savefig(OUTPUT_PDF)
    ok(f"Saved → {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
