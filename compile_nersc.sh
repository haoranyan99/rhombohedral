#!/bin/bash
set -euo pipefail

# Compile current production executables on NERSC Perlmutter CPU.
# Usage:
#   ./compile_nersc.sh                 # build all
#   ./compile_nersc.sh calculate_fermi # build one Makefile target

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

module purge
module load craype-x86-milan
module load craype-network-ofi
module load libfabric/1.22.0
module load PrgEnv-gnu/8.5.0
module load cray-mpich/8.1.30

echo "[INFO] Compile root: $ROOT_DIR"
echo "[INFO] Modules:"
module list 2>&1

TARGET="${1:-all}"
JOBS="${JOBS:-8}"

if [ "$TARGET" = "clean" ]; then
  make clean
  exit 0
fi

make -j "$JOBS" "$TARGET"

echo "[DONE] Built target: $TARGET"
echo "[DONE] Executables:"
ls -lh bin || true
