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

echo "Installing nf-test ${NFTEST_VERSION}..."

cd "${TEMP_DIR}"

curl -fsSL https://get.nf-test.com | bash -s "${NFTEST_VERSION}"

install -m 0755 nf-test "${INSTALL_DIR}/nf-test"

echo "nf-test installed:"
"${INSTALL_DIR}/nf-test" version