#!/usr/bin/env bash
# Run the preliminary pipeline inside WSL (Ubuntu-24.04, ~/edaenv).
#   wsl.exe -e bash -lc 'bash /mnt/c/Users/HUAWEI/Desktop/hbond-stretch-eda/scripts/run_wsl.sh'
# Environment variables: NPROC (threads, default 8), STAGES (default "opt scan fc density").
set -euo pipefail
source ~/edaenv/bin/activate
export OMP_NUM_THREADS=${NPROC:-8} OPENBLAS_NUM_THREADS=${NPROC:-8} MKL_NUM_THREADS=${NPROC:-8}
export TMPDIR=/tmp PYSCF_TMPDIR=/tmp
cd "$(dirname "$0")/.."

XC=${XC:-revpbe0}
DISP=${DISP:-d3bj}
BASIS=${BASIS:-def2-tzvp}
TAG=${TAG:-${XC}-${DISP}_${BASIS#def2-}}
GEOM=${GEOM:-geometries/dimer_revpbe0-d3bj_tzvp.xyz}
STAGES=${STAGES:-"opt scan fc density"}

for stage in $STAGES; do
  case "$stage" in
    opt)
      python scripts/00_optimize_dimer.py --xc "$XC" --disp "$DISP" --basis "$BASIS" --tag "$TAG" \
        2>&1 | tee "logs/opt_${TAG}.out" ;;
    scan)
      python scripts/01_scan_dimer.py --geom "$GEOM" --xc "$XC" --disp "$DISP" --basis "$BASIS" --tag "$TAG" \
        --rmin "${RMIN:-2.55}" --rmax "${RMAX:-3.60}" --step "${STEP:-0.05}" \
        --fine-halfwidth "${FINE_HW:-0.16}" --fine-step "${FINE_STEP:-0.02}" \
        ${OH_ELONG:+--oh-elongation "$OH_ELONG"} --restart 2>&1 | tee -a "logs/scan_${TAG}.out" ;;
    fc)
      python scripts/02_force_constants.py --scan "results/scan_${TAG}.json" 2>&1 | tee "logs/fc_${TAG}.out" ;;
    density)
      python scripts/03_density_profile.py --geom "$GEOM" --xc "$XC" --disp "$DISP" --basis "$BASIS" --tag "$TAG" \
        --R ${DENSITY_R:-2.70 2.80 2.90 3.00 3.10 3.20} 2>&1 | tee "logs/density_${TAG}.out" ;;
    *) echo "unknown stage $stage" >&2; exit 2 ;;
  esac
done
