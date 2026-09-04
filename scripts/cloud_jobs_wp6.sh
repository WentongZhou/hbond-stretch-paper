#!/usr/bin/env bash
# Runs ON the VM from ~/hbond-stretch-eda:
#   nohup bash scripts/cloud_jobs_wp6.sh liquid 6 5 > logs/wp6_master.out 2>&1 &
#   nohup bash scripts/cloud_jobs_wp6.sh liquid_h3o 6 5 --density-fit > logs/wp6_h3o_master.out 2>&1 &
# <tag> <n_workers> <threads_per_worker> [extra args for the full scan, e.g. --density-fit]
# Workers claim clusters from results/liquid/index_<tag>.json with mkdir locks and run,
# for each: pair scan (two waters, frozen geometry) + fc, then full scan + fc.
# The fragment-set name is read from the first cluster JSON of the index.
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/edaenv/bin/activate
export TMPDIR=/scratch/tmp PYSCF_TMPDIR=/scratch/tmp
export EDA_MAX_MEMORY_MB=${EDA_MAX_MEMORY_MB:-20000}
tag=${1:-liquid}; nw=${2:-4}; nt=${3:-8}; shift 3 || true
EXTRA="$*"
index="results/liquid/index_${tag}.json"
names=$(python -c "import json; print(' '.join(e['name'] for e in json.load(open('$index'))))")
first=$(echo $names | cut -d' ' -f1)
fset=$(python -c "import json; print(list(json.load(open('results/liquid/cluster_${first}.json'))['fragment_sets'])[0])")
mkdir -p logs/wp6 results/locks /scratch/tmp
GRID="--rmin -0.30 --rmax 0.30 --step 0.05 --fine-halfwidth 0"
PAIR_GRID="--rmin -0.30 --rmax 0.50 --step 0.05 --fine-halfwidth 0"
echo "$(date) tag=$tag set=$fset workers=$nw threads=$nt extra=[$EXTRA] clusters: $names"

worker() {
  export OMP_NUM_THREADS=$nt OPENBLAS_NUM_THREADS=$nt MKL_NUM_THREADS=$nt
  for name in $names; do
    mkdir "results/locks/$name" 2>/dev/null || continue
    cl="results/liquid/cluster_${name}.json"
    t0=$(date +%s)
    python scripts/07_scan_cluster.py --cluster "$cl" --set "$fset" --pair --tag "${name}_pair" $PAIR_GRID --restart \
      > "logs/wp6/scan_${name}_pair.out" 2>&1
    python scripts/02_force_constants.py --scan "results/scan_${name}_pair.json" --no-plot > "logs/wp6/fc_${name}_pair.out" 2>&1
    python scripts/07_scan_cluster.py --cluster "$cl" --set "$fset" --tag "${name}_full" $GRID --restart $EXTRA \
      > "logs/wp6/scan_${name}_full.out" 2>&1
    python scripts/02_force_constants.py --scan "results/scan_${name}_full.json" --no-plot > "logs/wp6/fc_${name}_full.out" 2>&1
    echo "$(date) done $name in $(( $(date +%s) - t0 )) s"
  done
}
for i in $(seq 1 "$nw"); do worker > "logs/wp6/worker_${tag}_$i.out" 2>&1 & done
wait
echo "$(date) all workers finished"
