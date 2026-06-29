#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


RUN_CONFIG = {
    "glc_0127": {
        "epitope": "GLCTLVAML",
        "dataset_name": "GLC",
        "clusters": "vdjdb_glc_vdbscan_leiden_grid_0127_tcremp_clusters.tsv",
    },
    "ylq_0371": {
        "epitope": "YLQPRTFLL",
        "dataset_name": "YLQ",
        "clusters": "vdjdb_ylq_vdbscan_leiden_grid_0371_tcremp_clusters.tsv",
    },
}


def convert_run(run_dir: Path, *, epitope: str, dataset_name: str, filename: str) -> pd.DataFrame:
    path = run_dir / filename
    df = pd.read_csv(path, sep="\t").copy()
    df = df[df["cluster_id"].ne(-1)].copy()
    out = pd.DataFrame(
        {
            "gene": "TRB",
            "cdr3aa": df["cdr3aa_beta"].astype(str),
            "v.segm": df["v_beta"].astype(str),
            "j.segm": df["j_beta"].astype(str),
            "cid": "H.B." + epitope + "." + df["cluster_id"].astype(str),
            "antigen.epitope": epitope,
            "method": "redcea",
            "dataset_name": dataset_name,
        }
    )
    return out.drop_duplicates().reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert imported upstream RedCEA runs into benchmark cluster_members_TRB format.")
    parser.add_argument(
        "--runs-root",
        type=Path,
        default=Path("results") / "vdjdb_upstream_runs",
        help="Root directory containing imported upstream run folders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results") / "vdjdb_method_benchmark" / "redcea" / "cluster_members_TRB.txt",
        help="Output cluster-members TSV for benchmark evaluation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frames = []
    for run_id, cfg in RUN_CONFIG.items():
        frames.append(
            convert_run(
                args.runs_root / run_id,
                epitope=cfg["epitope"],
                dataset_name=cfg["dataset_name"],
                filename=cfg["clusters"],
            )
        )
    out = pd.concat(frames, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Saved {args.output}")
    print(out.groupby('antigen.epitope')["cid"].nunique().to_string())


if __name__ == "__main__":
    main()
