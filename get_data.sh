#!/usr/bin/env bash
# get_data.sh — Download TCGA KIRC PanCancer Atlas 2018 data files
#
# Downloads data_mutations.txt and data_clinical_sample.txt from the
# cBioPortal datahub git-LFS store into data/.
#
# Requirements: git, git-lfs
#   macOS:  brew install git-lfs
#   Linux:  sudo apt install git-lfs   (Debian/Ubuntu)
#           sudo yum install git-lfs   (RHEL/Amazon Linux)

set -euo pipefail

STUDY="kirc_tcga_pan_can_atlas_2018"
REPO="https://github.com/cBioPortal/datahub.git"
TMPDIR_CLONE="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_CLONE"' EXIT
# Absolute path: the script later cd's into the temp clone
DATADIR="$(cd "$(dirname "$0")" && pwd)/data"

if [[ -s "$DATADIR/data_mutations.txt" && -s "$DATADIR/data_clinical_sample.txt" ]]; then
    echo "Data already present in $DATADIR; nothing to do."
    exit 0
fi

echo "==> Checking git-lfs..."
if ! command -v git-lfs &>/dev/null; then
    echo "ERROR: git-lfs not found. Install it first:" >&2
    echo "  macOS:  brew install git-lfs" >&2
    echo "  Ubuntu: sudo apt install git-lfs" >&2
    echo "  RHEL:   sudo yum install git-lfs" >&2
    exit 1
fi

git lfs install --skip-smudge 2>/dev/null || true

echo "==> Cloning datahub (sparse, no blobs)..."
git clone \
    --no-checkout \
    --depth=1 \
    --filter=blob:none \
    "$REPO" \
    "$TMPDIR_CLONE"

echo "==> Checking out study directory (LFS pointers only)..."
cd "$TMPDIR_CLONE"
git sparse-checkout init --cone
git sparse-checkout set "public/$STUDY"
GIT_LFS_SKIP_SMUDGE=1 git checkout master

echo "==> Fetching LFS objects for the two target files..."
git lfs pull --include \
    "public/$STUDY/data_mutations.txt,public/$STUDY/data_clinical_sample.txt"

echo "==> Copying files to data/..."
mkdir -p "$DATADIR"
cp "public/$STUDY/data_mutations.txt"      "$DATADIR/data_mutations.txt"
cp "public/$STUDY/data_clinical_sample.txt" "$DATADIR/data_clinical_sample.txt"

if head -c 100 "$DATADIR/data_mutations.txt" | grep -q "git-lfs"; then
    echo "ERROR: got a git-lfs pointer instead of the real file." >&2
    exit 1
fi

echo ""
echo "Done. Files written to $DATADIR:"
ls -lh "$DATADIR/data_mutations.txt" "$DATADIR/data_clinical_sample.txt"
