#!/bin/sh
# piv_pressure/cfd/run_cases.sh
# Generate, mesh and run the coarse / medium / fine cylinder cases, then
# build the grid-independence table and export the PIV window.
#
# Usage (from piv_pressure/, with the OpenFOAM environment loaded):
#     sh cfd/run_cases.sh [nProcs]            # nProcs > 1 runs each case in parallel
#     sh cfd/run_cases.sh [nProcs] --steady   # also run the optional steady Re = 20 case
#
# Expected run time grows ~8x per level (4x cells, ~2x smaller time step).
set -e
cd "$(dirname "$0")/.." || exit 1
NPROCS=${1:-1}
STEADY=${2:-}

command -v blockMesh >/dev/null 2>&1 || {
    echo "OpenFOAM not found: source its etc/bashrc (or open the OpenFOAM app/shell) first."; exit 1; }

python3 cfd/generate_cases.py
for level in coarse medium fine; do
    echo "=== running $level ==="
    sh "cfd/cases/$level/Allrun" "$NPROCS"
done

if [ "$STEADY" = "--steady" ]; then
    echo "=== running steady_re20 (optional; stops automatically at steady state) ==="
    sh "cfd/cases/steady_re20/Allrun" "$NPROCS"
fi

python3 cfd/grid_independence.py
python3 cfd/export_window.py
python3 cfd/parse_log.py
python3 cfd/change_norm.py
echo "Done. Now run:  python3 run_all.py"
