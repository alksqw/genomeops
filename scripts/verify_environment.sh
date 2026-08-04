#!/usr/bin/env bash

set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${PROJECT_ROOT}/.tool-versions.env"

export PATH="${HOME}/.local/bin:${PATH}"

echo
echo "========================================"
echo " GenomeOps environment verification"
echo "========================================"
echo

check_command() {
    local command_name="$1"

    if command -v "${command_name}" >/dev/null 2>&1; then
        echo "[OK] ${command_name}: $(command -v "${command_name}")"
    else
        echo "[ERROR] Command not found: ${command_name}" >&2
        exit 1
    fi
}

echo "Checking commands..."
echo

check_command git
check_command curl
check_command python
check_command java
check_command docker
check_command nextflow
check_command nf-test
check_command nf-core

echo
echo "Checking versions..."
echo

echo "--- Git ---"
git --version

echo
echo "--- Python ---"
python --version

echo
echo "--- Java ---"
java -version

echo
echo "--- Docker client ---"
docker --version

echo
echo "--- Nextflow ---"
nextflow -version

echo
echo "--- nf-test ---"
nf-test version

echo
echo "--- nf-core/tools ---"
nf-core --version

echo
echo "Checking Docker daemon..."

DOCKER_READY=false

for attempt in {1..15}; do
    if docker info >/dev/null 2>&1; then
        DOCKER_READY=true
        break
    fi

    echo "Docker is not ready yet. Attempt ${attempt}/15..."
    sleep 2
done

if [[ "${DOCKER_READY}" != "true" ]]; then
    echo "[ERROR] Docker command exists, but Docker daemon is unavailable." >&2
    exit 1
fi

echo "[OK] Docker daemon is available."

echo
echo "Expected project tool versions:"
echo "  Nextflow: ${NEXTFLOW_VERSION}"
echo "  nf-test:  ${NFTEST_VERSION}"
echo "  nf-core:  ${NFCORE_VERSION}"

echo
echo "========================================"
echo " Environment verification completed"
echo "========================================"