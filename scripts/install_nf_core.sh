#!/usr/bin/env bash

set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${PROJECT_ROOT}/.tool-versions.env"

INSTALL_DIR="${HOME}/.local/bin"
VENV_DIR="${HOME}/.venvs/nf-core"

mkdir -p "${INSTALL_DIR}"
mkdir -p "$(dirname "${VENV_DIR}")"

echo "Installing nf-core/tools ${NFCORE_VERSION}..."

rm -rf "${VENV_DIR}"

python -m venv "${VENV_DIR}"

"${VENV_DIR}/bin/python" -m pip install --upgrade pip

"${VENV_DIR}/bin/python" -m pip install \
    "nf-core==${NFCORE_VERSION}"

ln -sfn \
    "${VENV_DIR}/bin/nf-core" \
    "${INSTALL_DIR}/nf-core"

echo "nf-core/tools installed:"
"${INSTALL_DIR}/nf-core" --version