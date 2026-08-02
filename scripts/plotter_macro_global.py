"""
plotter_macro_global.py
Generates paper-quality PDF plots from the CSV produced by profiler_macro_global.py.
Input : macro_results_global_{dataset}.csv
Output: macro_01_latency_{dataset}.pdf
        macro_02_throughput_{dataset}.pdf
        macro_03_speedup_heatmap_{dataset}.pdf
        macro_04_statistics_{dataset}.pdf / .csv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.cm as cm
import seaborn as sns
from pathlib import Path
import sys

# ─── CONFIG ──────────────────────────────────────────────────────────────────
DATASETS = ["single", "small_multi"]

# ─── PLOT STYLE ──────────────────────────────────────────────────────────────
sns.set_theme(style="ticks", context="paper")
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.labelsize":    13,
    "axes.titlesize":    13,
    "axes.titleweight":  "bold",
    "axes.grid":         True,
    "grid.alpha":        0.3,
    "grid.linestyle":    "--",
    "legend.fontsize":   9,
    "figure.titlesize":  14,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

# Color ramps: within each backend, darker = larger population
def _backend_colors(backend: str, n: int):
    if backend == "CUDA":
        return cm.Reds(np.linspace(0.40, 0.88, n))
    else:
        return cm.Blues(np.linspace(0.40, 0.88, n))

MARKERS = ["o", "s", "D", "^"]

# ─── HELPERS ─────────────────────────────────────────────────────────────────
def banner(title: str, width: int = 68) -> None:
    print(f"\n╔{'═' * (width - 2)}╗")
    print(f"║  {title:<{width - 4}}║")
    print(f"╚{'═' * (width - 2)}╝")

def section(msg: str) -> None:
    print(f"\n  ▶  {msg}")

def ok(msg: str) -> None:
    print(f"  ✔  Saved → {msg}")

def warn(msg: str) -> None:
    print(f"  ⚠  {msg}", file=sys.stderr)


# ─── PLOT FUNCTIONS ──────────────────────────────────────────────────────────
def plot_latency(df: pd.DataFrame, ds_name: str) -> str:
    mean_df = df.groupby(["Backend", "Population", "Generations"])["Time (s)"].mean().reset_index()
    populations = sorted(mean_df["Population"].unique())

    fig, ax = plt.subplots(figsize=(11, 6))

    for backend, grp in mean_df.groupby("Backend"):
        colors = _backend_colors(backend, len(populations))
        for idx, pop in enumerate(populations):
            pop_grp = grp[grp["Population"] == pop].sort_values("Generations")
            if pop_grp.empty:
                continue
            ax.plot(
                pop_grp["Generations"], pop_grp["Time (s)"],
                color=colors[idx], marker=MARKERS[idx],
                linewidth=2, markersize=7,
                label=f"{backend}  Pop={pop}",
                linestyle="-" if backend == "CUDA" else "--",
            )

    ax.set_title(f"Latency Scaling — {ds_name}")
    ax.set_xlabel("Number of Generations")
    ax.set_ylabel("Execution Time (s)")
    ax.legend(loc="upper left", framealpha=0.85, ncol=2, fontsize=9)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    sns.despine()

    out = f"macro_01_latency_{ds_name}.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_throughput(df: pd.DataFrame, ds_name: str) -> str:
    mean_df = df.groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"].mean().reset_index()
    populations = sorted(mean_df["Population"].unique())

    fig, ax = plt.subplots(figsize=(11, 6))

    for backend, grp in mean_df.groupby("Backend"):
        colors = _backend_colors(backend, len(populations))
        for idx, pop in enumerate(populations):
            pop_grp = grp[grp["Population"] == pop].sort_values("Generations")
            if pop_grp.empty:
                continue
            ax.plot(
                pop_grp["Generations"], pop_grp["Throughput (Evals/s)"],
                color=colors[idx], marker=MARKERS[idx],
                linewidth=2, markersize=7,
                label=f"{backend}  Pop={pop}",
                linestyle="-" if backend == "CUDA" else "--",
            )

    ax.set_title(f"Throughput Scaling — {ds_name}")
    ax.set_xlabel("Number of Generations")
    ax.set_ylabel("Throughput (Evals/s)")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.legend(loc="lower right", framealpha=0.85, ncol=2, fontsize=9)
    sns.despine()

    # Annotate gap at highest-pop / highest-gen point
    try:
        cuda_peak   = mean_df[(mean_df["Backend"]=="CUDA")   & (mean_df["Population"]==populations[-1])]["Throughput (Evals/s)"].max()
        alpaka_peak = mean_df[(mean_df["Backend"]=="Alpaka") & (mean_df["Population"]==populations[-1])]["Throughput (Evals/s)"].max()
        gap_pct = (alpaka_peak / cuda_peak - 1) * 100
        sign    = "+" if gap_pct >= 0 else ""
        color   = "#27AE60" if gap_pct >= 0 else "#C0392B"
        ax.annotate(
            f"Gap @ Pop={populations[-1]}: Alpaka {sign}{gap_pct:.1f}% vs CUDA",
            xy=(0.97, 0.05), xycoords="axes fraction", ha="right",
            fontsize=10, fontweight="bold", color=color,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7, edgecolor=color),
        )
    except Exception:
        pass

    out = f"macro_02_throughput_{ds_name}.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_speedup_heatmap(df: pd.DataFrame, ds_name: str) -> str | None:
    mean_df = df.groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"].mean().reset_index()
    df_cuda   = mean_df[mean_df["Backend"] == "CUDA"].set_index(["Population", "Generations"])
    df_alpaka = mean_df[mean_df["Backend"] == "Alpaka"].set_index(["Population", "Generations"])

    if df_cuda.empty or df_alpaka.empty:
        warn("Cannot generate speedup heatmap — missing one backend.")
        return None

    common = df_cuda.index.intersection(df_alpaka.index)
    if common.empty:
        warn("No common (Population, Generation) pairs between backends.")
        return None

    speedup = (df_alpaka.loc[common, "Throughput (Evals/s)"] /
               df_cuda.loc[common, "Throughput (Evals/s)"]).reset_index()
    speedup.rename(columns={"Throughput (Evals/s)": "Speedup"}, inplace=True)
    pivot = speedup.pivot(index="Population", columns="Generations", values="Speedup")

    fig, ax = plt.subplots(figsize=(9, 7))
    sns.heatmap(
        pivot, annot=True, fmt=".3f",
        cmap="RdYlGn", center=1.0,
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Speedup Alpaka / CUDA  (1.0 = parity)"},
        ax=ax,
    )
    ax.set_title(f"Relative Speedup: Alpaka vs CUDA — {ds_name}")
    ax.set_xlabel("Generations")
    ax.set_ylabel("Population")
    ax.invert_yaxis()

    # Add a 1.0 baseline annotation
    ax.annotate("< 1.0 = CUDA faster   |   > 1.0 = Alpaka faster",
                xy=(0.5, -0.08), xycoords="axes fraction",
                ha="center", fontsize=9, color="#555")

    out = f"macro_03_speedup_heatmap_{ds_name}.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_statistics_table(df: pd.DataFrame, ds_name: str) -> tuple[str, str]:
    stats = (df.groupby(["Backend", "Population", "Generations"])["Throughput (Evals/s)"]
             .agg(Mean="mean", Std="std", Min="min", Max="max")
             .reset_index())

    csv_out = f"macro_04_statistics_{ds_name}.csv"
    stats.to_csv(csv_out, index=False)

    # Format for display
    display = stats.copy()
    for col in ["Mean", "Std", "Min", "Max"]:
        display[col] = display[col].map(lambda x: f"{x:,.1f}")

    col_labels = ["Backend", "Pop", "Gen", "Mean (Evals/s)", "Std", "Min", "Max"]

    fig, ax = plt.subplots(figsize=(13, max(4, len(display) * 0.42 + 1.8)))
    ax.axis("off")

    table = ax.table(
        cellText=display.values,
        colLabels=col_labels,
        loc="center", cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.1, 1.55)

    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor("#2C3E50")
        else:
            backend = display.iloc[row - 1]["Backend"]
            base = "#FDEDEC" if backend == "CUDA" else "#EBF5FB"
            cell.set_facecolor(base if row % 2 == 1 else "white")
        cell.set_edgecolor("#D5D8DC")

    ax.set_title(
        f"Statistical Summary: Throughput (Evals/s) — {ds_name}",
        pad=16, fontweight="bold", fontsize=12, loc="left",
    )

    pdf_out = f"macro_04_statistics_{ds_name}.pdf"
    fig.savefig(pdf_out, bbox_inches="tight")
    plt.close(fig)
    return pdf_out, csv_out


# ─── MAIN ────────────────────────────────────────────────────────────────────
def generate_macro_plots() -> None:
    banner("MACRO PLOTTER  —  muDock  |  Paper-Quality PDF Generation")

    for ds_name in DATASETS:
        csv_file = f"macro_results_global_{ds_name}.csv"
        if not Path(csv_file).exists():
            warn(f"File not found: {csv_file}  (skipping '{ds_name}')")
            continue

        df = pd.read_csv(csv_file)
        if df.empty:
            warn(f"CSV empty: {csv_file}  (skipping)")
            continue

        section(f"Dataset: {ds_name.upper()}  ({len(df)} rows)")

        ok(plot_latency(df, ds_name))
        ok(plot_throughput(df, ds_name))

        heatmap_out = plot_speedup_heatmap(df, ds_name)
        if heatmap_out:
            ok(heatmap_out)

        pdf_out, csv_out = plot_statistics_table(df, ds_name)
        ok(pdf_out)
        ok(csv_out)

    print()


if __name__ == "__main__":
    generate_macro_plots()
