from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DONOR_ORDER = ["P1", "P2", "Q1", "Q2", "S1", "S2"]
DONOR_COLORS = {
    "P1": "#1b9e77",
    "P2": "#66a61e",
    "Q1": "#d95f02",
    "Q2": "#e6ab02",
    "S1": "#7570b3",
    "S2": "#e7298a",
    "mixed": "#666666",
    "noise": "#c7c7c7",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run TCRemP on combined YFV clonotypes and build PCA50->UMAP scatter plots.",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("results/yfv_ms5_ek8_lr1/raw"),
        help="Directory with per-donor TCRemP result tables.",
    )
    parser.add_argument(
        "--input-pattern",
        default="*_enriched_clonotypes_tcremp.tsv",
        help="Glob pattern used inside --raw-dir to select per-donor clonotype tables.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search for input files recursively under --raw-dir.",
    )
    parser.add_argument(
        "--path-must-contain",
        nargs="*",
        default=(),
        help="Optional substrings that must all be present in the matched file path.",
    )
    parser.add_argument(
        "--donors",
        nargs="*",
        default=DONOR_ORDER,
        help="Donor IDs to include.",
    )
    parser.add_argument(
        "--allow-multiple-files-per-donor",
        action="store_true",
        help="Allow multiple matched input files per donor after filtering.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/yfv_ms5_ek8_lr1/article/tcremp_all_clonotypes"),
        help="Directory for the combined run, derived tables, and plots.",
    )
    parser.add_argument(
        "--prefix",
        default="yfv_all_donors_ms5_ek8_lr1",
        help="Output prefix for the combined TCRemP run.",
    )
    parser.add_argument(
        "--cluster-pc-components",
        type=int,
        default=50,
        help="Number of PCA components requested from TCRemP and used for UMAP.",
    )
    parser.add_argument(
        "--cluster-min-samples",
        type=int,
        default=5,
        help="DBSCAN min_samples passed to TCRemP.",
    )
    parser.add_argument(
        "--k-neighbors",
        type=int,
        default=20,
        help="k-th neighbor parameter passed to TCRemP.",
    )
    parser.add_argument(
        "--umap-neighbors",
        type=int,
        default=30,
        help="UMAP n_neighbors.",
    )
    parser.add_argument(
        "--umap-min-dist",
        type=float,
        default=0.15,
        help="UMAP min_dist.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for TCRemP and UMAP.",
    )
    parser.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Worker count for TCRemP.",
    )
    parser.add_argument(
        "--n-prototypes",
        type=int,
        default=512,
        help="Number of TCRemP prototypes to use for the combined run.",
    )
    parser.add_argument(
        "--skip-tcremp",
        action="store_true",
        help="Reuse existing combined TCRemP outputs if they are already present.",
    )
    return parser.parse_args()


def extract_donor_id(path: Path) -> str:
    parts = path.name.split("_")
    if len(parts) < 2:
        raise ValueError(f"Could not infer donor from {path.name}")
    return parts[1]


def discover_input_paths(
    raw_dir: Path,
    input_pattern: str,
    recursive: bool,
    path_must_contain: tuple[str, ...] | list[str],
    donors: tuple[str, ...] | list[str],
    allow_multiple_files_per_donor: bool,
) -> list[Path]:
    candidates = sorted(raw_dir.rglob(input_pattern) if recursive else raw_dir.glob(input_pattern))
    allowed_donors = set(donors)
    required_substrings = tuple(path_must_contain)

    filtered: list[Path] = []
    donor_to_paths: dict[str, list[Path]] = {}

    for path in candidates:
        donor_id = extract_donor_id(path)
        if donor_id not in allowed_donors:
            continue
        path_str = str(path).replace("\\", "/")
        if any(substr not in path_str for substr in required_substrings):
            continue
        filtered.append(path)
        donor_to_paths.setdefault(donor_id, []).append(path)

    if not filtered:
        details = f"pattern={input_pattern!r}, recursive={recursive}, required_substrings={required_substrings!r}"
        raise FileNotFoundError(f"No files matched under {raw_dir} with {details}")

    missing_donors = [donor for donor in donors if donor not in donor_to_paths]
    if missing_donors:
        raise FileNotFoundError(f"No matched files found for donor(s): {', '.join(missing_donors)}")

    if not allow_multiple_files_per_donor:
        duplicate_donors = {donor: paths for donor, paths in donor_to_paths.items() if len(paths) > 1}
        if duplicate_donors:
            lines = ["Multiple files matched for the same donor; refine --path-must-contain or --input-pattern:"]
            for donor, paths in sorted(duplicate_donors.items()):
                lines.append(f"{donor}:")
                lines.extend(f"  {path}" for path in paths)
            raise ValueError("\n".join(lines))

    return filtered


def load_combined_input(
    raw_dir: Path,
    input_pattern: str,
    recursive: bool,
    path_must_contain: tuple[str, ...] | list[str],
    donors: tuple[str, ...] | list[str],
    allow_multiple_files_per_donor: bool,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    input_paths = discover_input_paths(
        raw_dir=raw_dir,
        input_pattern=input_pattern,
        recursive=recursive,
        path_must_contain=path_must_contain,
        donors=donors,
        allow_multiple_files_per_donor=allow_multiple_files_per_donor,
    )
    for path in input_paths:
        donor_id = extract_donor_id(path)
        df = pd.read_csv(path, sep="\t", usecols=["clone_id", "cdr3aa_beta", "v_beta", "j_beta"])
        df = df.rename(
            columns={
                "cdr3aa_beta": "junction_aa",
                "v_beta": "v_call",
                "j_beta": "j_call",
            }
        ).copy()
        df["clone_id"] = donor_id + "__" + df["clone_id"].astype(str)
        df["donor_id"] = donor_id
        df["locus"] = "beta"
        df["source_path"] = str(path)
        rows.append(df)

    combined = pd.concat(rows, ignore_index=True)
    combined = combined.dropna(subset=["clone_id", "junction_aa", "v_call", "j_call", "locus"]).copy()
    combined["junction_aa"] = combined["junction_aa"].astype(str).str.strip()
    combined["v_call"] = combined["v_call"].astype(str).str.strip()
    combined["j_call"] = combined["j_call"].astype(str).str.strip()
    combined = combined.loc[combined["junction_aa"].ne("")].copy()
    return combined


def run_tcremp(args: argparse.Namespace, input_path: Path, output_dir: Path) -> None:
    from tcremp.tcremp_run import main as tcremp_main

    mplconfigdir = Path(tempfile.gettempdir()) / "mplconfig-redcea-icml-2026"
    mplconfigdir.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mplconfigdir)

    cli_args = [
        "tcremp-run",
        "--input",
        str(input_path),
        "--output",
        str(output_dir),
        "--prefix",
        args.prefix,
        "--index-col",
        "clone_id",
        "--labels-col",
        "donor_id",
        "--chain",
        "TRB",
        "--n-prototypes",
        str(args.n_prototypes),
        "--cluster-pc-components",
        str(args.cluster_pc_components),
        "--cluster-min-samples",
        str(args.cluster_min_samples),
        "--k-neighbors",
        str(args.k_neighbors),
        "--random-seed",
        str(args.random_state),
        "--nproc",
        str(args.nproc),
        "--no-save-dists",
        "--tsne",
    ]

    original_argv = sys.argv[:]
    try:
        sys.argv = cli_args
        tcremp_main()
    finally:
        sys.argv = original_argv


def composition_label(donors: list[str], is_noise: bool) -> str:
    if is_noise:
        return "noise"
    if not donors:
        return "unassigned"
    if len(donors) == 1:
        return donors[0]
    if len(donors) == 2:
        return "+".join(donors)
    return "mixed"


def dominant_label(counts: pd.Series, is_noise: bool) -> str:
    if is_noise:
        return "noise"
    counts = counts.sort_values(ascending=False)
    if counts.empty:
        return "unassigned"
    top_donor = str(counts.index[0])
    if len(counts) == 1:
        return top_donor
    if counts.iloc[0] > counts.iloc[1]:
        return top_donor
    return "mixed"


def build_cluster_annotations(cluster_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    membership = (
        cluster_df.loc[:, ["cluster_id", "clone_id", "donor_id"]]
        .drop_duplicates()
        .copy()
    )
    cluster_counts = (
        membership.groupby(["cluster_id", "donor_id"])
        .size()
        .rename("n_clonotypes")
        .reset_index()
    )
    pivot = (
        cluster_counts.pivot(index="cluster_id", columns="donor_id", values="n_clonotypes")
        .fillna(0)
        .astype(int)
        .reindex(columns=DONOR_ORDER, fill_value=0)
    )
    cluster_meta = pivot.reset_index()
    cluster_meta["is_noise"] = cluster_meta["cluster_id"].eq(-1)
    donor_sets = [
        [donor for donor in DONOR_ORDER if row[donor] > 0]
        for _, row in cluster_meta.iterrows()
    ]
    cluster_meta["donor_composition"] = [
        composition_label(donors, is_noise)
        for donors, is_noise in zip(donor_sets, cluster_meta["is_noise"])
    ]
    cluster_meta["donor_presence"] = [";".join(donors) if donors else "" for donors in donor_sets]
    cluster_meta["dominant_donor"] = [
        dominant_label(row[DONOR_ORDER], bool(row["is_noise"]))
        for _, row in cluster_meta.iterrows()
    ]
    return cluster_meta, cluster_counts


def run_umap_embedding(
    pca_df: pd.DataFrame,
    n_neighbors: int,
    min_dist: float,
    random_state: int,
) -> pd.DataFrame:
    import umap

    pc_cols = [col for col in pca_df.columns if col.startswith("PC")]
    if not pc_cols:
        raise ValueError("No PCA columns starting with 'PC' were found in the TCRemP PCA output.")

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric="euclidean",
        random_state=random_state,
    )
    coords = reducer.fit_transform(pca_df.loc[:, pc_cols].to_numpy(dtype=float))
    out = pca_df.copy()
    out["umap_x"] = coords[:, 0]
    out["umap_y"] = coords[:, 1]
    return out


def make_composition_palette(labels: list[str]) -> dict[str, str]:
    pair_palette = {
        "P1+P2": "#1f78b4",
        "Q1+Q2": "#ff7f00",
        "S1+S2": "#6a3d9a",
    }
    palette: dict[str, str] = {}
    for label in labels:
        if label in DONOR_COLORS:
            palette[label] = DONOR_COLORS[label]
        elif label in pair_palette:
            palette[label] = pair_palette[label]
        elif "+" in label:
            palette[label] = "#a6761d"
        else:
            palette[label] = "#666666"
    return palette


def plot_composition_scatter(df: pd.DataFrame, out_path: Path) -> None:
    labels = sorted(df["donor_composition"].dropna().astype(str).unique().tolist())
    palette = make_composition_palette(labels)

    fig, ax = plt.subplots(figsize=(10, 8))
    for label in labels:
        sub = df.loc[df["donor_composition"].eq(label)]
        ax.scatter(
            sub["umap_x"],
            sub["umap_y"],
            s=8,
            alpha=0.72 if label != "noise" else 0.25,
            c=palette[label],
            linewidths=0,
            label=label,
        )

    ax.set_title("TCRemP UMAP colored by donor composition of each cluster")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(
        title="Cluster composition",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=2.5,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_presence_panels(df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
    base = df.loc[:, ["umap_x", "umap_y"]]

    for ax, donor in zip(axes.ravel(), DONOR_ORDER):
        present_mask = df[donor].gt(0)
        ax.scatter(
            base["umap_x"],
            base["umap_y"],
            s=5,
            c="#d9d9d9",
            alpha=0.22,
            linewidths=0,
        )
        ax.scatter(
            df.loc[present_mask, "umap_x"],
            df.loc[present_mask, "umap_y"],
            s=8,
            c=DONOR_COLORS[donor],
            alpha=0.85,
            linewidths=0,
        )
        ax.set_title(f"{donor} present in cluster")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")

    fig.suptitle("TCRemP UMAP with per-donor cluster presence", fontsize=14)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_dominant_scatter(df: pd.DataFrame, out_path: Path) -> None:
    labels = [label for label in DONOR_ORDER if label in set(df["dominant_donor"])] + ["mixed", "noise"]
    labels = [label for label in labels if label in set(df["dominant_donor"])]

    fig, ax = plt.subplots(figsize=(10, 8))
    for label in labels:
        sub = df.loc[df["dominant_donor"].eq(label)]
        ax.scatter(
            sub["umap_x"],
            sub["umap_y"],
            s=8,
            alpha=0.72 if label != "noise" else 0.25,
            c=DONOR_COLORS.get(label, "#666666"),
            linewidths=0,
            label=label,
        )

    ax.set_title("TCRemP UMAP colored by dominant donor within cluster")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(
        title="Dominant donor",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=2.5,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_presence_summary(cluster_meta: pd.DataFrame, cluster_counts: pd.DataFrame, output_dir: Path, prefix: str) -> None:
    cluster_meta.to_csv(output_dir / f"{prefix}_cluster_donor_presence.tsv", sep="\t", index=False)
    cluster_counts.to_csv(output_dir / f"{prefix}_cluster_donor_counts_long.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    mplconfigdir = Path(tempfile.gettempdir()) / "mplconfig-redcea-icml-2026"
    mplconfigdir.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mplconfigdir)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    combined_input = load_combined_input(
        raw_dir=args.raw_dir,
        input_pattern=args.input_pattern,
        recursive=args.recursive,
        path_must_contain=args.path_must_contain,
        donors=args.donors,
        allow_multiple_files_per_donor=args.allow_multiple_files_per_donor,
    )
    combined_input_path = output_dir / f"{args.prefix}_input.tsv"
    combined_input.to_csv(combined_input_path, sep="\t", index=False)

    pca_path = output_dir / f"{args.prefix}_tcremp_pca.tsv"
    cluster_path = output_dir / f"{args.prefix}_tcremp_clusters.tsv"

    if not args.skip_tcremp:
        run_tcremp(args, combined_input_path, output_dir)

    if not pca_path.exists():
        raise FileNotFoundError(f"Missing PCA output: {pca_path}")
    if not cluster_path.exists():
        raise FileNotFoundError(f"Missing cluster output: {cluster_path}")

    pca_df = pd.read_csv(pca_path, sep="\t")
    cluster_df = pd.read_csv(cluster_path, sep="\t")

    cluster_meta, cluster_counts = build_cluster_annotations(cluster_df)
    umap_df = run_umap_embedding(
        pca_df,
        n_neighbors=args.umap_neighbors,
        min_dist=args.umap_min_dist,
        random_state=args.random_state,
    )
    scatter_df = umap_df.merge(cluster_meta, on="cluster_id", how="left")

    scatter_path = output_dir / f"{args.prefix}_pca50_umap_points.tsv"
    scatter_df.to_csv(scatter_path, sep="\t", index=False)
    save_presence_summary(cluster_meta, cluster_counts, output_dir, args.prefix)

    plot_composition_scatter(scatter_df, output_dir / f"{args.prefix}_pca50_umap_cluster_composition.png")
    plot_composition_scatter(scatter_df, output_dir / f"{args.prefix}_pca50_umap_cluster_composition.svg")
    plot_presence_panels(scatter_df, output_dir / f"{args.prefix}_pca50_umap_cluster_presence_panels.png")
    plot_presence_panels(scatter_df, output_dir / f"{args.prefix}_pca50_umap_cluster_presence_panels.svg")
    plot_dominant_scatter(scatter_df, output_dir / f"{args.prefix}_pca50_umap_dominant_donor.png")
    plot_dominant_scatter(scatter_df, output_dir / f"{args.prefix}_pca50_umap_dominant_donor.svg")


if __name__ == "__main__":
    main()
