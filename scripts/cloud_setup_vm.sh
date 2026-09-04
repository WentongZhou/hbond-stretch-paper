#!/usr/bin/env bash
# Runs ON the GCE VM (copied there by cloud_vm.sh setup).  Idempotent.
#  - RAID0 of all local NVMe SSDs mounted at /scratch (PySCF scratch space)
#  - Python venv ~/edaenv with pyscf, pyscf-dispersion, ase, matplotlib
set -euo pipefail
sudo apt-get update -qq >/dev/null
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq python3-venv python3-dev build-essential mdadm >/dev/null

if ! mountpoint -q /scratch; then
  disks=$(ls /dev/disk/by-id/google-local-nvme-ssd-* 2>/dev/null | tr '\n' ' ' || true)
  n=$(echo $disks | wc -w)
  if [ "$n" -gt 0 ]; then
    echo "building RAID0 from $n local NVMe SSDs"
    sudo mdadm --create /dev/md0 --level=0 --raid-devices="$n" $disks --force --run
    sudo mkfs.ext4 -F -q /dev/md0
    sudo mkdir -p /scratch
    sudo mount -o discard,defaults /dev/md0 /scratch
    sudo chmod 1777 /scratch
  else
    echo "no local SSDs found; using boot disk for /scratch"
    sudo mkdir -p /scratch && sudo chmod 1777 /scratch
  fi
fi
mkdir -p /scratch/tmp
df -h /scratch | tail -1

if [ ! -d ~/edaenv ]; then python3 -m venv ~/edaenv; fi
source ~/edaenv/bin/activate
pip install -q --upgrade pip
pip install -q pyscf pyscf-dispersion ase numpy scipy matplotlib
python -c "import pyscf, ase; print('pyscf', pyscf.__version__, 'ase', ase.__version__)"
nproc
free -g | head -2
