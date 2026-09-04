#!/usr/bin/env bash
# Runs ON the VM from ~/hbond-stretch-eda (via cloud_vm.sh run 'bash scripts/cloud_jobs_wp3.sh').
# Launches the WP3 jobs in the background:
#   job 1 (12 threads): revPBE0-D3(BJ)/def2-QZVP rigid scan + force-constant partition
#   job 2 (20 threads): CP-corrected CCSD(T) scans, aug-cc-pVTZ then aug-cc-pVQZ
set -euo pipefail
cd "$(dirname "$0")/.."
source ~/edaenv/bin/activate
export TMPDIR=/scratch/tmp PYSCF_TMPDIR=/scratch/tmp
mkdir -p logs results figures /scratch/tmp

if ! pgrep -f "01_scan_dimer.py.*qzvp" >/dev/null; then
  XC=revpbe0 DISP=d3bj BASIS=def2-qzvp TAG=revpbe0-d3bj_qzvp STAGES="scan fc" NPROC=12 \
    nohup bash scripts/run_wsl.sh > logs/pipeline_qzvp.out 2>&1 < /dev/null &
  echo "started QZVP scan (pid $!)"
fi

if ! pgrep -f "08_ccsdt_scan.py" >/dev/null; then
  OMP_NUM_THREADS=20 OPENBLAS_NUM_THREADS=20 MKL_NUM_THREADS=20 nohup bash -c '
    python scripts/08_ccsdt_scan.py --geom geometries/dimer_revpbe0-d3bj_tzvp.xyz --basis aug-cc-pvtz --tag avtz \
      --rmin 2.70 --rmax 3.30 --step 0.05 --fine-halfwidth 0.16 --fine-step 0.02 --max-memory 200000 --restart \
      > logs/ccsdt_avtz.out 2>&1
    python scripts/08_ccsdt_scan.py --geom geometries/dimer_revpbe0-d3bj_tzvp.xyz --basis aug-cc-pvqz --tag avqz \
      --rmin 2.70 --rmax 3.30 --step 0.05 --fine-halfwidth 0.16 --fine-step 0.02 --max-memory 200000 --restart \
      > logs/ccsdt_avqz.out 2>&1
  ' > logs/ccsdt_chain.out 2>&1 < /dev/null &
  echo "started CCSD(T) chain (pid $!)"
fi
sleep 2
pgrep -af "python scripts" | cut -c1-140
