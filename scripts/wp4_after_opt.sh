#!/usr/bin/env bash
# Wait for a cluster optimization (06) to finish, then scan every requested
# fragment set (07) and partition the force constant (02).  Run inside WSL:
#   NPROC=7 setsid nohup bash scripts/wp4_after_opt.sh acid acid > logs/wp4_acid.out 2>&1 &
#   NPROC=7 setsid nohup bash scripts/wp4_after_opt.sh trimer acid_control base_control > logs/wp4_trimer.out 2>&1 &
set -euo pipefail
system=$1; shift
sets=("$@")
TAG=${TAG:-revpbe0-d3bj_tzvp}
XC=${XC:-revpbe0}; DISP=${DISP:-d3bj}; BASIS=${BASIS:-def2-tzvp}
cd "$(dirname "$0")/.."
source ~/edaenv/bin/activate
export OMP_NUM_THREADS=${NPROC:-7} OPENBLAS_NUM_THREADS=${NPROC:-7} MKL_NUM_THREADS=${NPROC:-7}
export TMPDIR=/tmp PYSCF_TMPDIR=/tmp

cluster="results/cluster_${system}_${TAG}.json"
echo "# waiting for $cluster"
until [ -f "$cluster" ] && ! pgrep -f "06_build_clusters.py --system $system " >/dev/null; do
  if grep -q Traceback "logs/opt_cluster_${system}.out" 2>/dev/null && ! pgrep -f "06_build_clusters.py --system $system " >/dev/null; then
    echo "# optimization of $system failed; see logs/opt_cluster_${system}.out"; exit 1
  fi
  sleep 30
done
echo "# $(date) optimization done, starting scans: ${sets[*]}"
for set in "${sets[@]}"; do
  tag="${system}_${set}_${TAG}"
  [ "$set" = "$system" ] && tag="${system}_${TAG}"
  python scripts/07_scan_cluster.py --cluster "$cluster" --set "$set" --xc "$XC" --disp "$DISP" --basis "$BASIS" \
    --tag "$tag" --rmin "${RMIN:--0.30}" --rmax "${RMAX:-0.50}" --step "${STEP:-0.05}" \
    --fine-halfwidth "${FINE_HW:-0.16}" --fine-step "${FINE_STEP:-0.02}" --restart 2>&1 | tee "logs/scan_${tag}.out"
  python scripts/02_force_constants.py --scan "results/scan_${tag}.json" 2>&1 | tee "logs/fc_${tag}.out"
done
echo "# $(date) all done for $system"
