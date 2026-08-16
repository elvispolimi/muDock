"""
plotter_macro_rapid.py
Generates a PDF comparison plot from the CSV produced by profiler_macro_rapid.py.
Input : macro_results_rapid_{dataset}.csv
Output: macro_rapid_{dataset}.pdf
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path
import sys

# ─── CONFIG ──────────────────────────────────────────────────────────────────
DATASETS = ["single"]

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
    "legend.fontsize":   10,
    "figure.titlesize":  14,
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
    print(f"  ✔  Saved → {msg}")

def warn(msg: str) -> None:
    print(f"  ⚠  {msg}", file=sys.stderr)


# ─── PLOT ────────────────────────────────────────────────────────────────────
def plot_rapid(ds_name: str, df: pd.DataFrame) -> str:
    backends   = list(df["Backend"].unique())
    pop        = int(df["Population"].iloc[0])
    gen        = int(df["Generations"].iloc[0])
    n_runs     = df.groupby("Backend")["Run"].count().max()
    cfg_label  = f"Pop={pop}  ·  Gen={gen}  ·  {n_runs} runs"

    stats = (df.groupby("Backend")["Throughput (Evals/s)"]
             .agg(Mean="mean", Std="std", Min="min", Max="max")
             .reset_index())

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    fig.suptitle(
        f"Rapid Benchmark — {ds_name}  ·  {cfg_label}",
        fontsize=14, fontweight="bold", y=1.02,
    )

    # ── LEFT: grouped bar chart with error bars ───────────────────────────────
    ax = axes[0]
    x      = np.arange(len(backends))
    palette = sns.color_palette("muted", n_colors=len(backends))

    bars = ax.bar(
        x, stats.set_index("Backend").loc[backends, "Mean"],
        yerr=stats.set_index("Backend").loc[backends, "Std"],
        color=palette, edgecolor="white", linewidth=0.8,
        width=0.48, capsize=7,
        error_kw={"linewidth": 1.8, "capthick": 1.8, "ecolor": "#333"},
    )

    # Value labels
    for bar, b in zip(bars, backends):
        mean = stats.loc[stats["Backend"] == b, "Mean"].values[0]
        std  = stats.loc[stats["Backend"] == b, "Std"].values[0]
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            mean + std + stats["Mean"].max() * 0.012,
            f"{mean:,.0f}", ha="center", va="bottom",
            fontsize=11, fontweight="bold",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(backends, fontsize=11)
    ax.set_ylabel("Throughput (Evals/s)")
    ax.set_title("Mean Throughput ± Std Dev")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    # Y axis from 0 for honest comparison
    ax.set_ylim(bottom=0)
    sns.despine(ax=ax)

    # ── RIGHT: per-run dot plot (box + strip) ─────────────────────────────────
    ax2 = axes[1]
    rng = np.random.default_rng(42)

    for i, backend in enumerate(backends):
        color       = palette[i]
        color_light = sns.color_palette("pastel", n_colors=len(backends))[i]
        grp         = df[df["Backend"] == backend]["Throughput (Evals/s)"].values
        n           = len(grp)

        # IQR box (min/max range)
        lo, hi = grp.min(), grp.max()
        ax2.fill_between([i - 0.15, i + 0.15], [lo, lo], [hi, hi],
                         color=color_light, alpha=0.4, zorder=1)
        # Jittered dots
        jitter = rng.uniform(-0.10, 0.10, n)
        ax2.scatter(i + jitter, grp, color=color, s=70,
                    zorder=3, alpha=0.9, edgecolors="white", linewidths=0.6)
        # Mean line
        mean = grp.mean()
        ax2.hlines(mean, i - 0.2, i + 0.2, colors=color, linewidths=2.8, zorder=4)
        ax2.text(i + 0.22, mean, f"{mean:,.0f}", va="center",
                 fontsize=9, color=color, fontweight="bold")

    ax2.set_xticks(range(len(backends)))
    ax2.set_xticklabels(backends, fontsize=11)
    ax2.set_ylabel("Throughput (Evals/s)")
    ax2.set_title("Per-Run Scatter  (— = mean,  ▭ = range)")
    ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax2.set_ylim(bottom=0)
    sns.despine(ax=ax2)

    # Shared legend
    patches = [mpatches.Patch(color=palette[i], label=b)
               for i, b in enumerate(backends)]
    fig.legend(handles=patches, loc="lower center", ncol=len(backends),
               framealpha=0.85, bbox_to_anchor=(0.5, -0.04))

    plt.tight_layout()
    out = f"macro_rapid_cpu_{ds_name}.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


# ─── MAIN ────────────────────────────────────────────────────────────────────
def main() -> None:
    banner("RAPID PLOTTER  —  muDock")

    for ds_name in DATASETS:
        csv_file = Path(f"macro_results_rapid_cpu_{ds_name}.csv")
        if not csv_file.exists():
            warn(f"File not found: {csv_file}  (skipping '{ds_name}')")
            continue

        df = pd.read_csv(csv_file)
        if df.empty:
            warn(f"CSV empty: {csv_file}  (skipping)")
            continue

        section(f"Dataset: {ds_name.upper()}  ({len(df)} rows)")
        ok(plot_rapid(ds_name, df))

    print()


if __name__ == "__main__":
    main()
