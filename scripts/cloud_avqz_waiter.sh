#!/usr/bin/env bash
# Runs ON the VM.  Waits until the aug-cc-pVTZ CCSD(T) scan and the def2-QZVP
# DFT scan have finished, then runs the reduced-grid aug-cc-pVQZ CCSD(T) scan
# (9 points, 2.80-3.08 A, step 0.04, plus the reference R) on all 32 threads.
cd "$(dirname "$0")/.."
source ~/edaenv/bin/activate
export TMPDIR=/scratch/tmp PYSCF_TMPDIR=/scratch/tmp
# bracket trick: the pattern must not match this script's own command line
while pgrep -f "basis aug-cc-pvt[z]" >/dev/null; do sleep 60; done
while pgrep -f "01_scan_dime[r].py" >/dev/null; do sleep 60; done
echo "# $(date) starting aVQZ"
OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32 \
python scripts/08_ccsdt_scan.py --geom geometries/dimer_revpbe0-d3bj_tzvp.xyz --basis aug-cc-pvqz --tag avqz \
  --rmin 2.80 --rmax 3.08 --step 0.04 --fine-halfwidth 0 --fine-step 0.04 --max-memory 150000 --restart \
  > logs/ccsdt_avqz.out 2>&1
echo "# $(date) aVQZ done"
