#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_DIR="${REPO_ROOT}/data/benchmark/DF5120"
RELEASE_TAG="data-assets"
ASSET_NAME="df5120_benchmark_v1.tar.gz"
REPO="pgonzale60/gfapdefinisher"

mkdir -p "${TARGET_DIR}"

TMP_ARCHIVE="$(mktemp "${TMPDIR:-/tmp}/df5120_benchmark_XXXXXX.tar.gz")"
trap 'rm -f "${TMP_ARCHIVE}"' EXIT

echo "Fetching ${ASSET_NAME} from release ${RELEASE_TAG}..."
DOWNLOAD_SUCCESS=0
if command -v gh >/dev/null 2>&1; then
    echo "Attempting download with GitHub CLI..."
    if gh release download "${RELEASE_TAG}" -R "${REPO}" -p "${ASSET_NAME}" -O "${TMP_ARCHIVE}" --clobber; then
        DOWNLOAD_SUCCESS=1
    else
        echo "gh release download failed, falling back to curl..."
    fi
fi

if [ "${DOWNLOAD_SUCCESS}" -eq 0 ]; then
    echo "Downloading via curl..."
    curl -fsSL "https://github.com/${REPO}/releases/download/${RELEASE_TAG}/${ASSET_NAME}" -o "${TMP_ARCHIVE}"
fi

echo "Extracting ${ASSET_NAME} into ${TARGET_DIR}..."
tar -xzf "${TMP_ARCHIVE}" -C "${TARGET_DIR}"
echo "Successfully staged DF5120 benchmark dataset into ${TARGET_DIR}."
