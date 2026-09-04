"""Peek at the geometry of a running PySCF/ASE optimization through the
temporary SCF chkfiles PySCF leaves in $PYSCF_TMPDIR, or analyze an xyz file.

Usage:
    python scripts/peek_geometry.py                 # scan /tmp/tmp* chkfiles, write geometries/cluster_<sys>_snapshot.xyz
    python scripts/peek_geometry.py file.xyz [...]  # H-bond table of the given xyz files
"""
from __future__ import annotations

import glob
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

B = importlib.import_module("06_build_clusters")


def report(symbols, xyz, title):
    print(f"\n# {title}: {len(symbols)} atoms")
    for r in B.hbond_table(symbols, xyz):
        print(f"   O{r['donor_O']}-H{r['H']}...O{r['acceptor_O']}  R={r['R_OO']:.3f}  r(OH)={r['r_OH']:.3f}  angle={r['angle_OHO']:.1f}")
    o_idx = [i for i, s in enumerate(symbols) if s == "O"]
    print("   O..O from O0: " + " ".join(f"O{j}:{np.linalg.norm(xyz[j] - xyz[0]):.2f}" for j in o_idx[1:]))
    # nearest O for every H (covalent partner) to spot proton transfer
    for i, s in enumerate(symbols):
        if s == "H":
            d = [(np.linalg.norm(xyz[i] - xyz[j]), j) for j in o_idx]
            d.sort()
            if d[0][0] > 1.05 or (len(d) > 1 and d[1][0] < 1.35):
                print(f"   H{i}: nearest O{d[0][1]} {d[0][0]:.3f}, next O{d[1][1]} {d[1][0]:.3f}  <-- shared/transferring proton?")


def main() -> int:
    if len(sys.argv) > 1:
        for f in sys.argv[1:]:
            sym, xyz = C.read_xyz(f)
            report(sym, xyz, f)
        return 0
    from pyscf import lib

    for f in sorted(glob.glob("/tmp/tmp*")):
        try:
            mol = lib.chkfile.load_mol(f)
        except Exception:
            continue
        sym = [mol.atom_symbol(i) for i in range(mol.natm)]
        xyz = mol.atom_coords(unit="Angstrom")
        tag = {1: "acid", -1: "base", 0: "neutral"}.get(mol.charge, "q")
        report(sym, xyz, f"{f} charge={mol.charge} basis={mol.basis} -> {tag}")
        C.write_xyz(C.GEOM_DIR / f"cluster_{tag}_snapshot.xyz", sym, xyz, f"{tag} snapshot from running optimization")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
