#!/usr/bin/env bash
# Collect every WP6 scan (three VMs + local), refit all force constants with the
# current 02 (R_min and R0 blocks), and run the liquid statistics for the three tags.
#   wsl.exe -e bash -lc 'bash /mnt/c/Users/HUAWEI/Desktop/hbond-stretch-eda/scripts/wp6_finalize.sh [pull]'
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/edaenv/bin/activate
export OMP_NUM_THREADS=2
if [ "${1:-}" = "pull" ]; then
  PROJECT=<gcp-project-1>  VM=eda-wp3  bash scripts/cloud_vm.sh pull "scan_liquid*" 2>&1 | grep -E "rror" || true
  PROJECT=<gcp-project-2>  VM=eda-wp6b bash scripts/cloud_vm.sh pull "scan_liquid*" 2>&1 | grep -E "rror" || true
  PROJECT=<gcp-project-3>  VM=eda-wp6c bash scripts/cloud_vm.sh pull "scan_liquid*" 2>&1 | grep -E "rror" || true
fi
n=0
for f in results/scan_liquid*_pair.json results/scan_liquid*_embed.json results/scan_liquid*_full.json; do
  [ -f "$f" ] || continue
  python scripts/02_force_constants.py --scan "$f" --no-plot > /dev/null 2>&1 || echo "fit failed: $f"
  n=$((n + 1))
done
echo "# refitted $n scans"
for tag in liquid liquid_h3o liquid_oh; do
  python scripts/12_liquid_analysis.py --tag "$tag" 2>&1 | grep -E "^# wrote|^# figure|target H-bonds"
done
