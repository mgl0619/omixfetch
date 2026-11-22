#!/usr/bin/env bash
# Build everything omixfetch needs, in one go.
#
#   bash scripts/install.sh
#
# Creates the main environment (aria2 + omixfetch), installs the package, then
# builds the separate sra-tools environment that --convert fastq uses. Two
# environments because sra-tools and aria2 cannot be solved into one on
# osx-arm64. After this, --convert fastq finds the toolkit automatically.
set -euo pipefail

cd "$(dirname "$0")/.."

for runner in micromamba mamba conda; do
    if command -v "$runner" >/dev/null 2>&1; then RUNNER="$runner"; break; fi
done
if [ -z "${RUNNER:-}" ]; then
    echo "error: need micromamba, mamba or conda on PATH" >&2
    exit 1
fi
echo "using $RUNNER"

MAIN_ENV="${OMIXFETCH_ENV:-omixfetch}"

if "$RUNNER" env list | awk '{print $1}' | grep -qx "$MAIN_ENV"; then
    echo "environment '$MAIN_ENV' exists; updating"
    "$RUNNER" env update -y -n "$MAIN_ENV" -f environment.yml
else
    "$RUNNER" env create -y -n "$MAIN_ENV" -f environment.yml
fi

"$RUNNER" run -n "$MAIN_ENV" pip install -e ".[dev]"

# Builds omixfetch-convert, or does nothing if fasterq-dump is already usable.
"$RUNNER" run -n "$MAIN_ENV" omixfetch setup

echo
"$RUNNER" run -n "$MAIN_ENV" omixfetch --version
echo
echo "done. activate with:  $RUNNER activate $MAIN_ENV"
