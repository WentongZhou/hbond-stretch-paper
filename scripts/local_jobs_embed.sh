#!/usr/bin/env bash
# Embedded pair scans (07 --embed) + force-constant fits for every cluster of the given tags.
#   NPROC=2 setsid nohup bash scripts/local_jobs_embed.sh liquid liquid_h3o liquid_oh > logs/wp6_embed_local.out 2>&1 &
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/edaenv/bin/activate
export TMPDIR=/tmp PYSCF_TMPDIR=/tmp
export OMP_NUM_THREADS=${NPROC:-2} OPENBLAS_NUM_THREADS=${NPROC:-2} MKL_NUM_THREADS=${NPROC:-2}
mkdir -p logs/wp6 results/locks_embed
GRID="--rmin -0.30 --rmax 0.50 --step 0.05 --fine-halfwidth 0"
for tag in "$@"; do
  names=$(python -c "import json; print(' '.join(e['name'] for e in json.load(open('results/liquid/index_${tag}.json'))))")
  first=$(echo $names | cut -d' ' -f1)
  fset=$(python -c "import json; print(list(json.load(open('results/liquid/cluster_${first}.json'))['fragment_sets'])[0])")
  for name in $names; do
    mkdir "results/locks_embed/$name" 2>/dev/null || continue
    t0=$(date +%s)
    python scripts/07_scan_cluster.py --cluster "results/liquid/cluster_${name}.json" --set "$fset" --embed \
      --tag "${name}_embed" $GRID --restart > "logs/wp6/scan_${name}_embed.out" 2>&1
    python scripts/02_force_constants.py --scan "results/scan_${name}_embed.json" --no-plot > "logs/wp6/fc_${name}_embed.out" 2>&1
    echo "$(date) done embed $name in $(( $(date +%s) - t0 )) s"
  done
done
echo "$(date) embed batch finished"
