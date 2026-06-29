#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

ENV_NAME="${1:-clustcr-bench}"
VDJDB_SLIM="${VDJDB_SLIM:-$PWD/data/vdjdb/vdjdb.slim.txt}"
TRUTH_CSV="${TRUTH_CSV:-$PWD/data/private/01_05_2025_TCRvdb.csv}"
WORK_DIR="${WORK_DIR:-$PWD/work/vdjdb_method_benchmark_clustcr}"
RESULTS_DIR="${RESULTS_DIR:-$PWD/results/vdjdb_method_benchmark_clustcr}"
DATASET_NAME="${DATASET_NAME:-ALL}"
CLUSTCR_METHOD="${CLUSTCR_METHOD:-MCL}"
CLUSTCR_CPUS="${CLUSTCR_CPUS:-1}"
CLUSTCR_MIN_CLUSTER_SIZE="${CLUSTCR_MIN_CLUSTER_SIZE:-3}"
MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/.tmp/mplconfig}"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: '$1' not found in PATH" >&2
    exit 1
  }
}

require_path() {
  [[ -e "$1" ]] || {
    echo "Error: required path not found: $1" >&2
    exit 1
  }
}

require_cmd conda
require_cmd python
require_path "$VDJDB_SLIM"
require_path "$TRUTH_CSV"

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

export PYTHONNOUSERSITE=1
export MPLCONFIGDIR

mkdir -p "$WORK_DIR" "$RESULTS_DIR" "$MPLCONFIGDIR"

python scripts/benchmark_vdjdb_methods.py prepare \
  --vdjdb "$VDJDB_SLIM" \
  --truth "$TRUTH_CSV" \
  --work-dir "$WORK_DIR"

mkdir -p "$RESULTS_DIR/clustcr"
python scripts/benchmark_vdjdb_methods.py run-clustcr \
  --input "$WORK_DIR/inputs/generic/${DATASET_NAME}.tsv" \
  --output "$RESULTS_DIR/clustcr/cluster_members_TRB.txt" \
  --dataset-name "$DATASET_NAME" \
  --method "$CLUSTCR_METHOD" \
  --n-cpus "$CLUSTCR_CPUS" \
  --min-cluster-size "$CLUSTCR_MIN_CLUSTER_SIZE"

python scripts/benchmark_vdjdb_methods.py evaluate \
  --truth "$WORK_DIR/truth_glc_ylq.csv" \
  --results-root "$RESULTS_DIR" \
  --output "$RESULTS_DIR/metrics.csv"

echo
echo "Done."
echo "Results: $RESULTS_DIR"
