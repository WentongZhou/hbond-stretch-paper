#!/usr/bin/env bash
# Manage the Google Compute Engine VM used for the heavy WP3 jobs (def2-QZVP
# scan, CCSD(T)/aug-cc-pVTZ benchmark).  Run inside WSL, where gcloud lives:
#   wsl.exe -e bash -lc 'bash /mnt/c/Users/HUAWEI/Desktop/hbond-stretch-eda/scripts/cloud_vm.sh <cmd> [args]'
# Environment: VM (name), ZONE, MTYPE, PROJECT (GCP project id; use another project to get another 32-vCPU quota).
# Commands:
#   create            create the VM (n2-highmem-32, 16 local NVMe SSDs, 400 GB boot; project quota CPUS_ALL_REGIONS=32)
#   setup             copy cloud_setup_vm.sh to the VM and run it (RAID0 /scratch, venv, pyscf)
#   push              tar D:\pyscf-dm-eda (no .git/examples) and scripts/geometries/results, install editable
#   run  '<cmd>'      run a shell command in ~/hbond-stretch-eda inside the venv
#   pull '<glob>'     copy results/<glob> from the VM into the local results/
#   logs '<glob>'     copy logs/<glob> from the VM into the local logs/
#   status            uptime, load, running python jobs
#   delete            delete the VM (stops billing)
set -euo pipefail
VM=${VM:-eda-wp3}
ZONE=${ZONE:-us-east4-c}
MTYPE=${MTYPE:-n2-highmem-32}
PROJECT=${PROJECT:-your-gcp-project-id}     # each project has its own CPUS_ALL_REGIONS quota (32)
PROJ_DIR=/mnt/c/Users/HUAWEI/Desktop/hbond-stretch-eda
REPO_PARENT=/mnt/d
gcloud() { command gcloud --project "$PROJECT" "$@"; }
ssh_cmd() { gcloud compute ssh "$VM" --zone "$ZONE" --quiet --command "$1"; }

case "${1:-}" in
  create)
    ssd=""; for i in $(seq 1 16); do ssd="$ssd --local-ssd interface=nvme"; done
    gcloud compute instances create "$VM" --zone "$ZONE" --machine-type "$MTYPE"       --image-family ubuntu-2404-lts-amd64 --image-project ubuntu-os-cloud       --boot-disk-size 400GB --boot-disk-type pd-balanced $ssd       --labels purpose=hbond-eda ;;
  setup)
    gcloud compute scp "$PROJ_DIR/scripts/cloud_setup_vm.sh" "$VM":~ --zone "$ZONE" --quiet
    ssh_cmd 'bash ~/cloud_setup_vm.sh' ;;
  push)
    tar czf /tmp/dmeda.tgz --exclude='__pycache__' --exclude='*.pyc' -C "$REPO_PARENT" pyscf-dm-eda/src pyscf-dm-eda/pyproject.toml pyscf-dm-eda/README.md pyscf-dm-eda/LICENSE
    tar czf /tmp/proj.tgz --exclude='__pycache__' -C "$(dirname "$PROJ_DIR")" hbond-stretch-eda/scripts hbond-stretch-eda/geometries hbond-stretch-eda/results
    ls -la /tmp/dmeda.tgz /tmp/proj.tgz
    gcloud compute scp /tmp/dmeda.tgz /tmp/proj.tgz "$VM":~ --zone "$ZONE" --quiet
    ssh_cmd 'tar xzf dmeda.tgz && tar xzf proj.tgz && source ~/edaenv/bin/activate && pip install -q -e ~/pyscf-dm-eda && mkdir -p ~/hbond-stretch-eda/logs ~/hbond-stretch-eda/figures && python -c "import pyscf_dm_eda; print(pyscf_dm_eda.__file__)"' ;;
  run)
    shift
    ssh_cmd "cd ~/hbond-stretch-eda && source ~/edaenv/bin/activate && export TMPDIR=/scratch/tmp PYSCF_TMPDIR=/scratch/tmp && $*" ;;
  pull)
    gcloud compute scp "$VM:~/hbond-stretch-eda/results/${2:-*}" "$PROJ_DIR/results/" --zone "$ZONE" --quiet ;;
  logs)
    gcloud compute scp "$VM:~/hbond-stretch-eda/logs/${2:-*}" "$PROJ_DIR/logs/" --zone "$ZONE" --quiet ;;
  status)
    ssh_cmd 'uptime; df -h /scratch | tail -1; echo; ps -eo pid,etime,pcpu,rss,args --sort=-pcpu | grep -E "python" | grep -v grep | cut -c1-160; echo; for f in ~/hbond-stretch-eda/logs/*.out; do echo "-- $f"; tail -2 "$f"; done' ;;
  delete)
    gcloud compute instances delete "$VM" --zone "$ZONE" --quiet ;;
  *)
    echo "usage: $0 create|setup|push|run '<cmd>'|pull '<glob>'|logs '<glob>'|status|delete" >&2
    exit 2 ;;
esac
