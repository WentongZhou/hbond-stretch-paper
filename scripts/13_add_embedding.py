"""Add an electrostatic-embedding charge set to liquid cluster specs (WP6).

For every results/liquid/cluster_<name>.json the source snapshot is re-read
and every molecule (other than the two target waters) with an atom within
--embed-cut of O_d or O_a is turned into point charges: TIP4P-Ew sites for
water (H +0.52422, M -1.04844 at 0.125 A from O along the bisector), the
model charges of H3O+ (O -0.32, H +0.44), OH- (O -1.32, H +0.32), Na+ (+1),
Cl- (-1).  Coordinates are unwrapped into the cluster frame (donor O at the
origin, same minimum-image convention as the extraction scripts).  The list
is stored under spec["embedding"] = {"coords": [...], "charges": [...]} and is
used by 07_scan_cluster.py --embed.

Usage:
    python scripts/13_add_embedding.py --tag liquid --embed-cut 8.0
    python scripts/13_add_embedding.py --tag liquid_h3o
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

Q_H, Q_M, D_OM = 0.52422, -1.04844, 0.125
ION_Q = {"H3O": {"O": -0.32, "H": 0.44}, "OHM": {"O": -1.32, "H": 0.32}, "NA": {"Na": 1.0}, "CL": {"Cl": -1.0}}


def min_image(v, L):
    return v - L * np.round(v / L)


def read_snapshot(xyz_path):
    from ase.io import read

    atoms = read(xyz_path, format="extxyz")
    symbols, pos, L = atoms.get_chemical_symbols(), atoms.get_positions(), float(atoms.cell.lengths()[0])
    meta = Path(xyz_path).with_suffix(".json")
    if meta.exists():
        residues = C.load_json(meta)["residues"]
    else:  # neutral water: O,H,H triples
        residues = [{"name": "HOH", "atoms": [i, i + 1, i + 2]} for i in range(0, len(symbols), 3)]
    return symbols, pos, L, residues


def water_sites(o, h1, h2):
    bis = (h1 - o) / np.linalg.norm(h1 - o) + (h2 - o) / np.linalg.norm(h2 - o)
    bis /= np.linalg.norm(bis)
    m = o + D_OM * bis
    return [(h1, Q_H), (h2, Q_H), (m, Q_M)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="liquid")
    parser.add_argument("--embed-cut", type=float, default=8.0, help="include molecules with any atom within this distance (A) of O_d or O_a")
    args = parser.parse_args()
    index = C.load_json(C.RESULT_DIR / "liquid" / f"index_{args.tag}.json")
    for entry in index:
        spec_path = C.RESULT_DIR / "liquid" / f"cluster_{entry['name']}.json"
        spec = C.load_json(spec_path)
        src = spec["source"]
        symbols, pos, L, residues = read_snapshot(src["snapshot"])
        od, oa = src["donor_O"], src["acceptor_O"]
        ref = pos[od]
        # the cluster frame: cluster coords = pos + shift - pos[od] with shift the min-image lattice vector
        coords, charges, kinds = [], [], []
        for res in residues:
            idx = res["atoms"]
            if od in idx or oa in idx:
                continue
            shift = min_image(pos[idx[0]] - ref, L) - (pos[idx[0]] - ref)
            atoms_here = pos[idx] + shift - ref
            d = min(np.linalg.norm(atoms_here - 0.0, axis=1).min(),
                    np.linalg.norm(atoms_here - (min_image(pos[oa] - ref, L)), axis=1).min())
            if d > args.embed_cut:
                continue
            if res["name"] == "HOH":
                o, h1, h2 = atoms_here[0], atoms_here[1], atoms_here[2]
                for c, q in water_sites(o, h1, h2):
                    coords.append(c.tolist()); charges.append(q); kinds.append("HOH")
            else:
                qmap = ION_Q[res["name"]]
                for i, c in zip(idx, atoms_here):
                    coords.append(c.tolist()); charges.append(qmap[symbols[i]]); kinds.append(res["name"])
        spec["embedding"] = {"coords": coords, "charges": charges, "cutoff_A": args.embed_cut,
                             "n_sites": len(charges), "n_molecules": len(charges) // 3 if all(k == "HOH" for k in kinds) else None,
                             "total_charge": float(np.sum(charges)),
                             "model": "TIP4P-Ew water sites; H3O+ O-0.32/H+0.44; OH- O-1.32/H+0.32; Na+ +1; Cl- -1"}
        C.dump_json(spec_path, spec)
        print(f"{entry['name']}: {len(charges)} charge sites, total charge {np.sum(charges):+.3f}, "
              f"kinds {sorted(set(kinds))}")
    print(f"# updated {len(index)} cluster specs for tag {args.tag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
