#!/usr/bin/env bash
# Wait for the cluster optimization (06), then scan the target bond with only
# its two waters kept at the frozen cluster geometry (07 --pair) and partition.
#   NPROC=4 setsid nohup bash scripts/wp4_pair_scan.sh acid acid > logs/wp4_pair_acid.out 2>&1 &
set -euo pipefail
system=$1; shift
sets=("$@")
TAG=${TAG:-revpbe0-d3bj_tzvp}
XC=${XC:-revpbe0}; DISP=${DISP:-d3bj}; BASIS=${BASIS:-def2-tzvp}
cd "$(dirname "$0")/.."
source ~/edaenv/bin/activate
export OMP_NUM_THREADS=${NPROC:-4} OPENBLAS_NUM_THREADS=${NPROC:-4} MKL_NUM_THREADS=${NPROC:-4}
export TMPDIR=/tmp PYSCF_TMPDIR=/tmp
cluster="results/cluster_${system}_${TAG}.json"
until [ -f "$cluster" ] && ! pgrep -f "06_build_clusters.py --system $system " >/dev/null; do sleep 60; done
echo "# $(date) pair scans for $system: ${sets[*]}"
for set in "${sets[@]}"; do
  tag="${system}_${set}_pair_${TAG}"
  [ "$set" = "$system" ] && tag="${system}_pair_${TAG}"
  python scripts/07_scan_cluster.py --cluster "$cluster" --set "$set" --pair --xc "$XC" --disp "$DISP" --basis "$BASIS" \
    --tag "$tag" --rmin "${RMIN:--0.30}" --rmax "${RMAX:-0.50}" --step "${STEP:-0.05}" \
    --fine-halfwidth "${FINE_HW:-0.16}" --fine-step "${FINE_STEP:-0.02}" --restart 2>&1 | tee "logs/scan_${tag}.out"
  python scripts/02_force_constants.py --scan "results/scan_${tag}.json" 2>&1 | tee "logs/fc_${tag}.out"
done
echo "# $(date) pair scans done for $system"
