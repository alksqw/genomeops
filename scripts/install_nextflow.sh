#!/usr/bin/env bash

set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${PROJECT_ROOT}/.tool-versions.env"

INSTALL_DIR="${HOME}/.local/bin"
TEMP_DIR="$(mktemp -d)"

cleanup() {
    rm -rf "${TEMP_DIR}"
}

trap cleanup EXIT

mkdir -p "${INSTALL_DIR}"

echo "Installing Nextflow ${NEXTFLOW_VERSION}..."

cd "${TEMP_DIR}"

export NXF_VER="${NEXTFLOW_VERSION}"

curl -fsSL https://get.nextflow.io | bash

install -m 0755 nextflow "${INSTALL_DIR}/nextflow"

echo "Nextflow installed:"
"${INSTALL_DIR}/nextflow" -version