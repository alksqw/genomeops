#!/usr/bin/env bash

set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "${PROJECT_ROOT}"

mkdir -p "${HOME}/.local/bin"

PATH_LINE='export PATH="$HOME/.local/bin:$PATH"'

if ! grep -qxF "${PATH_LINE}" "${HOME}/.bashrc"; then
    echo "${PATH_LINE}" >> "${HOME}/.bashrc"
fi

export PATH="${HOME}/.local/bin:${PATH}"

echo
echo "Starting GenomeOps environment setup..."
echo

echo "Creating the GenomeOps Python environment..."

python -m venv "${PROJECT_ROOT}/.venv"

source "${PROJECT_ROOT}/.venv/bin/activate"

python -m pip install --upgrade pip

python -m pip install -e ".[dev]"

bash scripts/install_nextflow.sh
bash scripts/install_nf_test.sh
bash scripts/install_nf_core.sh
bash scripts/verify_environment.sh

echo
echo "GenomeOps Codespace setup completed successfully."