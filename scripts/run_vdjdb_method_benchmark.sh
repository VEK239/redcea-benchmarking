#!/usr/bin/env bash
set -eo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

ENV_NAME="${1:-redcea-benchmark}"
METHODS="${METHODS:-tcrdist3 giana gliph2}"

VDJDB_SLIM="${VDJDB_SLIM:-$PWD/data/vdjdb/vdjdb.slim.txt}"
TRUTH_CSV="${TRUTH_CSV:-$PWD/data/private/01_05_2025_TCRvdb.csv}"
TOOLS_DIR="${TOOLS_DIR:-$PWD/tools}"
WORK_DIR="${WORK_DIR:-$PWD/work/vdjdb_method_benchmark}"
RESULTS_DIR="${RESULTS_DIR:-$PWD/results/vdjdb_method_benchmark}"

GIANA_HOME="${GIANA_HOME:-$TOOLS_DIR/GIANA}"
GLIPH2_RSCRIPT="${GLIPH2_RSCRIPT:-Rscript}"

TCRDIST3_RADIUS="${TCRDIST3_RADIUS:-24}"
TCRDIST3_MIN_CLUSTER_SIZE="${TCRDIST3_MIN_CLUSTER_SIZE:-3}"
TCRDIST3_CPUS="${TCRDIST3_CPUS:-4}"
TCRDIST3_CHUNK_SIZE="${TCRDIST3_CHUNK_SIZE:-100}"
GIANA_MIN_CLUSTER_SIZE="${GIANA_MIN_CLUSTER_SIZE:-3}"
GLIPH2_MIN_CLUSTER_SIZE="${GLIPH2_MIN_CLUSTER_SIZE:-3}"
DATASET_NAME="${DATASET_NAME:-ALL}"

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
require_cmd git

require_path "$VDJDB_SLIM"
require_path "$TRUTH_CSV"

if [[ " $METHODS " == *" tcrnet "* ]]; then
  echo "Error: tcrnet was removed from this benchmark workflow." >&2
  echo "Use the vdjdb-motifs workflow for TCRNet runs instead." >&2
  exit 1
fi

if [[ " $METHODS " == *" giana "* ]]; then
  require_path "$GIANA_HOME"
fi

if [[ " $METHODS " == *" gliph2 "* ]]; then
  require_cmd "$GLIPH2_RSCRIPT"
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

mkdir -p "$WORK_DIR" "$RESULTS_DIR"

python scripts/benchmark_vdjdb_methods.py prepare \
  --vdjdb "$VDJDB_SLIM" \
  --truth "$TRUTH_CSV" \
  --work-dir "$WORK_DIR"

if [[ " $METHODS " == *" tcrdist3 "* ]]; then
  mkdir -p "$RESULTS_DIR/tcrdist3"
  python scripts/benchmark_vdjdb_methods.py run-tcrdist3 \
    --input "$WORK_DIR/inputs/generic/${DATASET_NAME}.tsv" \
    --output "$RESULTS_DIR/tcrdist3/cluster_members_TRB.txt" \
    --dataset-name "$DATASET_NAME" \
    --radius "$TCRDIST3_RADIUS" \
    --cpus "$TCRDIST3_CPUS" \
    --chunk-size "$TCRDIST3_CHUNK_SIZE" \
    --min-cluster-size "$TCRDIST3_MIN_CLUSTER_SIZE"
fi

if [[ " $METHODS " == *" giana "* ]]; then
  mkdir -p "$RESULTS_DIR/giana" "$WORK_DIR/giana_raw"
  input_path="$WORK_DIR/inputs/giana/${DATASET_NAME}.tsv"
  raw_output="$WORK_DIR/giana_raw/${DATASET_NAME}.txt"
  python "$GIANA_HOME/GIANA4.1.py" -f "$input_path" -O "$raw_output"
  python scripts/benchmark_vdjdb_methods.py convert-giana \
    --input "$WORK_DIR/inputs/generic/${DATASET_NAME}.tsv" \
    --giana-output "$raw_output" \
    --output "$RESULTS_DIR/giana/cluster_members_TRB.txt" \
    --dataset-name "$DATASET_NAME" \
    --min-cluster-size "$GIANA_MIN_CLUSTER_SIZE"
fi

if [[ " $METHODS " == *" gliph2 "* ]]; then
  mkdir -p "$RESULTS_DIR/gliph2" "$WORK_DIR/gliph2_raw"
  gliph2_outputs=()
  for short_name in GLC YLQ; do
    input_path="$WORK_DIR/inputs/generic/${short_name}.tsv"
    raw_output="$WORK_DIR/gliph2_raw/${short_name}_cluster_members_TRB.txt"
    converted_output="$RESULTS_DIR/gliph2/cluster_members_${short_name}.txt"
    method_work_dir="$WORK_DIR/gliph2_raw/${short_name}"

    "$GLIPH2_RSCRIPT" scripts/run_gliph2_benchmark.R \
      "$input_path" \
      "$raw_output" \
      "$short_name" \
      "$GLIPH2_MIN_CLUSTER_SIZE" \
      "$method_work_dir"

    python scripts/benchmark_vdjdb_methods.py convert-gliph2 \
      --gliph2-output "$raw_output" \
      --output "$converted_output" \
      --dataset-name "$short_name"

    gliph2_outputs+=("$converted_output")
  done

  python scripts/benchmark_vdjdb_methods.py merge-cluster-members \
    --inputs "${gliph2_outputs[@]}" \
    --output "$RESULTS_DIR/gliph2/cluster_members_TRB.txt"
fi

python scripts/benchmark_vdjdb_methods.py evaluate \
  --truth "$WORK_DIR/truth_glc_ylq.csv" \
  --results-root "$RESULTS_DIR" \
  --output "$RESULTS_DIR/metrics.csv"

echo
echo "Done."
echo "Results: $RESULTS_DIR"
