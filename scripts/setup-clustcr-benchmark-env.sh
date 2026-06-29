#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-clustcr-bench}"
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: '$1' not found in PATH" >&2
    exit 1
  }
}

require_cmd conda
require_cmd git

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  conda create -y -n "$ENV_NAME" "python=${PYTHON_VERSION}" pip
fi

conda activate "$ENV_NAME"

python -m pip install --upgrade pip setuptools wheel

# Keep the environment intentionally minimal to avoid heavy conda solving.
python -m pip install \
  numpy \
  pandas \
  scipy \
  scikit-learn \
  networkx \
  matplotlib \
  faiss-cpu \
  markov_clustering \
  parmap \
  python-louvain

python -m pip install --no-deps "git+https://github.com/svalkiers/clusTCR"

cat <<EOF

Done.

Activate environment:
  conda activate $ENV_NAME

Installed:
  clustcr                 -> GitHub package (no-deps)
  numpy/pandas/scipy      -> pip
  scikit-learn            -> pip
  networkx/matplotlib     -> pip
  faiss-cpu               -> pip
  markov_clustering       -> pip
  parmap                  -> pip
  python-louvain          -> pip

Quick check:
  python -c "import clustcr; print('clustcr ok')"
EOF
