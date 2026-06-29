#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EPITOPE_ORDER = ["GLC", "YLQ"]
METHOD_STYLE = {
    "1mm Hamming baseline": {"color": "#c1121f"},
    "GIANA": {"color": "#43aa8b"},
    "GLIPH2": {"color": "#577590"},
    "clusTCR": {"color": "#6a994e"},
    "TCRdist3": {"color": "#f8961e"},
    "TCRnet": {"color": "#8c564b"},
    "RedCEA": {"color": "#277da1"},
    "RedCEA+VJ": {"color": "#264653"},
}
DISPLAY_LABEL = {
    "1mm Hamming baseline": "Hamming",
}
F1_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot VDJdb precision-recall scatter with F1 contours.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results") / "vdjdb_method_benchmark" / "paper_table_vdjdb_trb_comparison.csv",
        help="Benchmark comparison table CSV.",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=Path("results") / "vdjdb_method_benchmark" / "vdjdb_precision_recall_f1_contours.png",
        help="Output PNG path.",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=Path("results") / "vdjdb_method_benchmark" / "vdjdb_precision_recall_f1_contours.svg",
        help="Output SVG path.",
    )
    return parser.parse_args()


def f1_surface(precision: np.ndarray, recall: np.ndarray) -> np.ndarray:
    denom = precision + recall
    out = np.zeros_like(denom)
    mask = denom > 0
    out[mask] = 2 * precision[mask] * recall[mask] / denom[mask]
    return out


def add_f1_contours(ax: plt.Axes) -> None:
    p = np.linspace(0.001, 1.0, 900)
    r = np.linspace(0.001, 1.0, 900)
    pp, rr = np.meshgrid(p, r)
    ff = f1_surface(pp, rr)

    contour = ax.contour(
        pp,
        rr,
        ff,
        levels=F1_LEVELS,
        cmap="viridis",
        linewidths=1.4,
        alpha=0.85,
    )
    ax.clabel(contour, inline=True, fontsize=9, fmt=lambda value: f"F1={value:.1f}")


def plot_panel(ax: plt.Axes, frame: pd.DataFrame, epitope: str) -> None:
    add_f1_contours(ax)

    for _, row in frame.iterrows():
        method = str(row["Method"])
        style = METHOD_STYLE.get(method, {"color": "#333333"})
        x = float(row["Precision"])
        y = float(row["Recall"])
        ax.scatter(
            x,
            y,
            s=120,
            color=style["color"],
            edgecolors="white",
            linewidths=0.8,
            zorder=3,
        )
        label = DISPLAY_LABEL.get(method, method)
        xytext = (6, 6)
        ha = "left"
        if epitope == "GLC" and method == "clusTCR":
            xytext = (-8, 6)
            ha = "right"
        if epitope == "YLQ" and method in {"1mm Hamming baseline", "TCRdist3"}:
            xytext = (-8, 6)
            ha = "right"
        ax.annotate(
            label,
            (x, y),
            textcoords="offset points",
            xytext=xytext,
            ha=ha,
            fontsize=10,
            color="#222222",
        )

    ax.set_title(epitope, fontsize=14, fontweight="bold")
    ax.set_xlabel("Precision")
    ax.set_ylabel("Recall")
    ax.set_xlim(0.0, 1.02)
    ax.set_ylim(0.0, 1.02)
    ax.grid(alpha=0.2, linestyle="-", linewidth=0.8)
    ax.set_axisbelow(True)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)

    plt.rcParams.update(
        {
            "figure.figsize": (14, 6.5),
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 10,
        }
    )

    fig, axes = plt.subplots(1, len(EPITOPE_ORDER), constrained_layout=True)
    if len(EPITOPE_ORDER) == 1:
        axes = [axes]

    for ax, epitope in zip(axes, EPITOPE_ORDER):
        frame = df[df["Epitope"].eq(epitope)].copy()
        plot_panel(ax, frame, epitope)

    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_png, dpi=220, bbox_inches="tight", pad_inches=0.25)
    fig.savefig(args.output_svg, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)

    print(f"Saved {args.output_png}")
    print(f"Saved {args.output_svg}")


if __name__ == "__main__":
    main()
