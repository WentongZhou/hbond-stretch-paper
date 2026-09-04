#!/usr/bin/env bash
# Runs ON the VM: wait for the neutral WP6 batch, then run the acid and base batches.
#   nohup bash scripts/cloud_jobs_wp6_queue.sh > logs/wp6_queue.out 2>&1 &
cd "$(dirname "$0")/.."
while ! grep -q "all workers finished" logs/wp6_master.out 2>/dev/null; do sleep 120; done
echo "$(date) neutral batch finished; starting h3o"
bash scripts/cloud_jobs_wp6.sh liquid_h3o 6 5 --density-fit > logs/wp6_h3o_master.out 2>&1
echo "$(date) h3o finished; starting oh"
bash scripts/cloud_jobs_wp6.sh liquid_oh 6 5 --density-fit > logs/wp6_oh_master.out 2>&1
echo "$(date) queue finished"
