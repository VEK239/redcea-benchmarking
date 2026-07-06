#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".tmp") / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path("results") / "vdjdb_method_benchmark" / "paper_table_vdjdb_trb_comparison.csv"
DEFAULT_OUTPUT_PNG = Path("results") / "vdjdb_method_benchmark" / "vdjdb_precision_recall_f1_contours.png"
DEFAULT_OUTPUT_SVG = Path("results") / "vdjdb_method_benchmark" / "vdjdb_precision_recall_f1_contours.svg"

EPITOPE_ORDER = ["GLC", "YLQ"]
METHOD_ORDER = [
    "1mm Hamming baseline",
    "GIANA",
    "GLIPH2",
    "clusTCR",
    "TCRdist3",
    "TCRnet",
    "RedCEA",
    "RedCEA+VJ",
]
DISPLAY_LABEL = {
    "1mm Hamming baseline": "Hamming",
}
F1_LEVELS = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
HIGH_F1_LEVELS = [0.80, 0.90]
HIGH_F1_LABEL_POSITIONS = {
    "GLC": {
        0.80: (0.945, 0.74),
        0.90: (0.905, 0.885),
    },
    "YLQ": {
        0.80: (0.945, 0.74),
        0.90: (0.905, 0.885),
    },
}
DEFAULT_LABEL_PLACEMENT = (6, -8, "left", "top")
LABEL_PLACEMENTS = {
    ("GLC", "1mm Hamming baseline"): (8, 0, "left", "center"),
    ("GLC", "GIANA"): (10, 0, "left", "center"),
    ("GLC", "GLIPH2"): (10, 0, "left", "center"),
    ("GLC", "clusTCR"): (-6, 0, "right", "center"),
    ("GLC", "TCRdist3"): (-10, 0, "right", "center"),
    ("GLC", "TCRnet"): (-10, 0, "right", "center"),
    ("GLC", "RedCEA"): (0, 8, "center", "bottom"),
    ("GLC", "RedCEA+VJ"): (8, 0, "left", "center"),
    ("YLQ", "1mm Hamming baseline"): (-12, 0, "right", "center"),
    ("YLQ", "GIANA"): (12, 0, "left", "center"),
    ("YLQ", "GLIPH2"): (12, 0, "left", "center"),
    ("YLQ", "clusTCR"): (12, 0, "left", "center"),
    ("YLQ", "TCRdist3"): (-12, 0, "right", "center"),
    ("YLQ", "TCRnet"): (12, 0, "left", "center"),
    ("YLQ", "RedCEA"): (12, 0, "left", "center"),
    ("YLQ", "RedCEA+VJ"): (12, 0, "left", "center"),
}
METHOD_COLORS = {
    "RedCEA": "#b2182b",
    "RedCEA+VJ": "#d6604d",
    "TCRnet": "#f4a582",
    "GLIPH2": "#92c5de",
    "TCRdist3": "#2166ac",
    "1mm Hamming baseline": "#4393c3",
    "GIANA": "#b8e186",
    "clusTCR": "#7fbc41",
}
CONTOUR_COLORS = [str(value) for value in np.linspace(0.78, 0.35, len(F1_LEVELS))]
PLOT_RCPARAMS = {
    "figure.figsize": (11.2, 5.8),
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333",
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "Times New Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot VDJdb precision-recall scatter with F1 contours.")
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Benchmark comparison table CSV.",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=DEFAULT_OUTPUT_PNG,
        help="Output PNG path.",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=DEFAULT_OUTPUT_SVG,
        help="Output SVG path.",
    )
    return parser.parse_args()


def f1_surface(precision: np.ndarray, recall: np.ndarray) -> np.ndarray:
    denom = precision + recall
    out = np.zeros_like(denom)
    mask = denom > 0
    out[mask] = 2 * precision[mask] * recall[mask] / denom[mask]
    return out


def add_f1_contours(ax: plt.Axes, epitope: str) -> None:
    p = np.linspace(0.001, 1.0, 900)
    r = np.linspace(0.001, 1.0, 900)
    pp, rr = np.meshgrid(p, r)
    ff = f1_surface(pp, rr)
    contour = ax.contour(
        pp,
        rr,
        ff,
        levels=F1_LEVELS,
        colors=CONTOUR_COLORS,
        linewidths=1.8,
        alpha=0.95,
    )
    ax.clabel(
        contour,
        inline=True,
        inline_spacing=8,
        fontsize=14,
        fmt=lambda value: f"F1 {value:.1f}",
        use_clabeltext=True,
        levels=[level for level in F1_LEVELS if level not in HIGH_F1_LEVELS],
    )
    for level in HIGH_F1_LEVELS:
        ax.clabel(
            contour,
            inline=True,
            inline_spacing=8,
            fontsize=14,
            fmt=lambda value: f"F1 {value:.1f}",
            use_clabeltext=True,
            levels=[level],
            manual=[HIGH_F1_LABEL_POSITIONS[epitope][level]],
        )


def plot_panel(ax: plt.Axes, frame: pd.DataFrame, epitope: str) -> None:
    add_f1_contours(ax, epitope)

    for _, row in frame.iterrows():
        method = str(row["Method"])
        x = float(row["Precision"])
        y = float(row["Recall"])
        ax.scatter(
            x,
            y,
            s=120,
            color=METHOD_COLORS.get(method, "#333333"),
            edgecolors="white",
            linewidths=0.8,
            zorder=3,
        )
        label = DISPLAY_LABEL.get(method, method)
        dx, dy, ha, va = LABEL_PLACEMENTS.get((epitope, method), DEFAULT_LABEL_PLACEMENT)
        ax.annotate(
            label,
            (x, y),
            textcoords="offset points",
            xytext=(dx, dy),
            ha=ha,
            va=va,
            fontsize=14,
            color="#1f1f1f",
        )

    ax.set_title(epitope, fontsize=20, fontweight="normal", pad=12)
    ax.set_xlabel("Precision", fontsize=16)
    ax.set_ylabel("Recall", fontsize=16)
    ax.set_xlim(0.0, 1.02)
    ax.set_ylim(0.0, 1.06)
    ax.set_xticks(np.linspace(0.0, 1.0, 6))
    ax.set_yticks(np.linspace(0.0, 1.0, 6))
    ax.grid(alpha=0.22, linestyle="-", linewidth=0.7, color="#b0b0b0")
    ax.set_axisbelow(True)
    ax.set_box_aspect(1)
def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    plt.rcParams.update(PLOT_RCPARAMS)

    fig, axes = plt.subplots(1, len(EPITOPE_ORDER), constrained_layout=False, sharex=True, sharey=True)
    if len(EPITOPE_ORDER) == 1:
        axes = [axes]

    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.12, wspace=0.08)

    for ax, epitope in zip(axes, EPITOPE_ORDER):
        frame = df[df["Epitope"].eq(epitope)].copy()
        plot_panel(ax, frame, epitope)

    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_png, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    fig.savefig(args.output_svg, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)

    print(f"Saved {args.output_png}")
    print(f"Saved {args.output_svg}")


if __name__ == "__main__":
    main()
