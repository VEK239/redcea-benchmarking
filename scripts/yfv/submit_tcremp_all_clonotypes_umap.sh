#!/bin/bash

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/projects/immunestatus/pogorelyy/redcea-icml-2026}"
RAW_DIR="${RAW_DIR:-/results/redcea_runs}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/results/yfv_ms5_ek8_lr1/article/tcremp_all_clonotypes}"
LOG_DIR="${LOG_DIR:-$ROOT_DIR/logs/tcremp_all_clonotypes}"

PYTHON_BIN="${PYTHON_BIN:-python}"
CPUS_PER_TASK="${CPUS_PER_TASK:-8}"
MEMORY="${MEMORY:-64gb}"
TIME_LIMIT="${TIME_LIMIT:-24:00:00}"
PARTITION="${PARTITION:-medium}"
CONSTRAINT="${CONSTRAINT:-hpc}"

DONORS="${DONORS:-P1 P2 Q1 Q2 S1 S2}"
PREFIX="${PREFIX:-yfv_all_donors_ms5_ek8_lr1_np512_n8}"
RUN_PATH_FILTER="${RUN_PATH_FILTER:-grid_0011_ms5_kn20_ek8_ebs_srasym_lr1_pc50}"
INPUT_PATTERN="${INPUT_PATTERN:-*_enriched_clonotypes_tcremp.tsv}"
NPROC="${NPROC:-8}"
N_PROTOTYPES="${N_PROTOTYPES:-512}"
CLUSTER_PC_COMPONENTS="${CLUSTER_PC_COMPONENTS:-50}"
CLUSTER_MIN_SAMPLES="${CLUSTER_MIN_SAMPLES:-5}"
K_NEIGHBORS="${K_NEIGHBORS:-20}"
UMAP_NEIGHBORS="${UMAP_NEIGHBORS:-30}"
UMAP_MIN_DIST="${UMAP_MIN_DIST:-0.15}"
RANDOM_STATE="${RANDOM_STATE:-42}"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
read -r -a DONOR_ARRAY <<< "$DONORS"

SBATCH_SCRIPT="$(mktemp)"
cat > "$SBATCH_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=${PREFIX}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOG_DIR}/${PREFIX}.%j.log
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --partition=${PARTITION}

set -euo pipefail

cd $(printf '%q' "$ROOT_DIR")
mkdir -p $(printf '%q' "$OUTPUT_DIR")
DONORS=($(printf '%q ' "${DONOR_ARRAY[@]}"))

$(printf '%q' "$PYTHON_BIN") scripts/yfv/run_tcremp_all_clonotypes_umap.py \\
  --raw-dir $(printf '%q' "$RAW_DIR") \\
  --recursive \\
  --input-pattern $(printf '%q' "$INPUT_PATTERN") \\
  --path-must-contain $(printf '%q' "$RUN_PATH_FILTER") \\
  --output-dir $(printf '%q' "$OUTPUT_DIR") \\
  --prefix $(printf '%q' "$PREFIX") \\
  --donors "\${DONORS[@]}" \\
  --nproc $(printf '%q' "$NPROC") \\
  --n-prototypes $(printf '%q' "$N_PROTOTYPES") \\
  --cluster-pc-components $(printf '%q' "$CLUSTER_PC_COMPONENTS") \\
  --cluster-min-samples $(printf '%q' "$CLUSTER_MIN_SAMPLES") \\
  --k-neighbors $(printf '%q' "$K_NEIGHBORS") \\
  --umap-neighbors $(printf '%q' "$UMAP_NEIGHBORS") \\
  --umap-min-dist $(printf '%q' "$UMAP_MIN_DIST") \\
  --random-state $(printf '%q' "$RANDOM_STATE")
EOF

sbatch "$SBATCH_SCRIPT"
rm -f "$SBATCH_SCRIPT"
