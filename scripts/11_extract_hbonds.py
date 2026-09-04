"""Pick water-water H-bonds from liquid snapshots and cut clusters for the
DM-EDA force-constant scans (WP6).

For every snapshot (extended xyz with a periodic cell, molecules stored as
O,H,H triples) all H-bonds satisfying R(O..O) <= r_oo_max and
angle(O_d-H...O_a) >= angle_min are listed; ``--per-snapshot`` of them are
drawn at random (seeded).  For each target the cluster contains the donor
water D, the acceptor water A and every other water whose O lies within
``--env-cut`` of O_d or O_a (minimum image, molecules kept whole).  The
cluster xyz and a spec JSON in the layout of 06_build_clusters.py are
written so 07_scan_cluster.py can scan the target bond:

  fragment set "liquid":  [D + environment] / [A],  A moves along O_d -> O_a.

Usage:
    python scripts/11_extract_hbonds.py --snapshots geometries/liquid/snap_*.xyz --per-snapshot 2 --seed 11
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def read_extxyz(path):
    from ase.io import read

    atoms = read(path, format="extxyz")
    return atoms.get_chemical_symbols(), atoms.get_positions(), np.array(atoms.cell.lengths())


def min_image(v, L):
    return v - L * np.round(v / L)


def find_hbonds(symbols, pos, L, r_oo_max, angle_min):
    o_idx = [i for i, s in enumerate(symbols) if s == "O"]
    mol = {o: [o, o + 1, o + 2] for o in o_idx}  # O,H,H triples
    hbonds = []
    for od in o_idx:
        for h in mol[od][1:]:
            for oa in o_idx:
                if oa == od:
                    continue
                d_oo = min_image(pos[oa] - pos[od], L)
                r = float(np.linalg.norm(d_oo))
                if r > r_oo_max:
                    continue
                v_ho = min_image(pos[od] - pos[h], L)
                v_ha = min_image(pos[oa] - pos[h], L)
                cosang = np.dot(v_ho, v_ha) / np.linalg.norm(v_ho) / np.linalg.norm(v_ha)
                ang = float(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
                if ang >= angle_min:
                    hbonds.append({"donor_O": od, "H": h, "acceptor_O": oa, "R_OO": r, "angle": ang,
                                   "r_HA": float(np.linalg.norm(v_ha))})
    return hbonds, mol


def count_hbonds(hbonds, o):
    donated = sum(1 for hb in hbonds if hb["donor_O"] == o)
    accepted = sum(1 for hb in hbonds if hb["acceptor_O"] == o)
    return donated, accepted


def cut_cluster(symbols, pos, L, mol, target, env_cut):
    od, oa = target["donor_O"], target["acceptor_O"]
    ref = pos[od]
    env = []
    for o in mol:
        if o in (od, oa):
            continue
        d1 = np.linalg.norm(min_image(pos[o] - pos[od], L))
        d2 = np.linalg.norm(min_image(pos[o] - pos[oa], L))
        if min(d1, d2) <= env_cut:
            env.append((min(d1, d2), o))
    env.sort()
    order = [od, oa] + [o for _, o in env]
    coords, syms = [], []
    for o in order:
        shift = min_image(pos[o] - ref, L) - (pos[o] - ref)  # lattice shift keeping the molecule whole
        for i in mol[o]:
            coords.append(pos[i] + shift)
            syms.append(symbols[i])
    coords = np.asarray(coords)
    coords -= coords[0]
    # bonded H of the donor: the one closer to the acceptor O
    hb_local = 1 if np.linalg.norm(coords[1] - coords[3]) < np.linalg.norm(coords[2] - coords[3]) else 2
    return syms, coords, hb_local, len(env)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshots", nargs="+", default=None)
    parser.add_argument("--per-snapshot", type=int, default=2)
    parser.add_argument("--r-oo-max", type=float, default=3.3)
    parser.add_argument("--angle-min", type=float, default=140.0)
    parser.add_argument("--env-cut", type=float, default=3.4, help="O..O cutoff for environment waters")
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--tag", default="liquid")
    args = parser.parse_args()
    files = args.snapshots or sorted(glob.glob(str(C.GEOM_DIR / "liquid" / "snap_*.xyz")))
    rng = np.random.default_rng(args.seed)
    geom_dir = C.GEOM_DIR / "liquid"
    res_dir = C.RESULT_DIR / "liquid"
    geom_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)

    index = []
    for f in files:
        symbols, pos, L = read_extxyz(f)
        hbonds, mol = find_hbonds(symbols, pos, L[0], args.r_oo_max, args.angle_min)
        n_o = sum(1 for s in symbols if s == "O")
        print(f"# {Path(f).name}: {len(hbonds)} H-bonds, {len(hbonds) / n_o:.2f} per water")
        picks = rng.choice(len(hbonds), size=min(args.per_snapshot, len(hbonds)), replace=False)
        snap = Path(f).stem.split("_")[-1]
        for j, p in enumerate(sorted(picks)):
            hb = hbonds[p]
            syms, coords, hb_local, n_env = cut_cluster(symbols, pos, L[0], mol, hb, args.env_cut)
            name = f"{args.tag}_{snap}_{j}"
            xyz_path = geom_dir / f"cluster_{name}.xyz"
            C.write_xyz(xyz_path, syms, coords, f"{name}: target O0-H{hb_local}...O3 from {Path(f).name}, "
                                                 f"R={hb['R_OO']:.3f} angle={hb['angle']:.1f}, {n_env} env waters")
            n_atoms = len(syms)
            dd, da = count_hbonds(hbonds, hb["donor_O"])
            ad, aa = count_hbonds(hbonds, hb["acceptor_O"])
            spec = {
                "system": name, "charge": 0, "tag": args.tag, "xyz": str(xyz_path.relative_to(C.ROOT)).replace("\\", "/"),
                "symbols": syms,
                "source": {"snapshot": str(f), "donor_O": hb["donor_O"], "H": hb["H"], "acceptor_O": hb["acceptor_O"]},
                "descriptors": {"R_OO": hb["R_OO"], "angle_OHO": hb["angle"], "r_HA": hb["r_HA"], "n_env_waters": n_env,
                                "donor_nHB_donated": dd, "donor_nHB_accepted": da,
                                "acceptor_nHB_donated": ad, "acceptor_nHB_accepted": aa},
                "fragment_sets": {
                    "liquid": {
                        "fragments": [{"label": "donor_env", "atoms": [0, 1, 2] + list(range(6, n_atoms)), "charge": 0},
                                      {"label": "water_a", "atoms": [3, 4, 5], "charge": 0}],
                        "moving": "water_a", "donor_O": 0, "bonded_H": hb_local, "acceptor_O": 3, "ct_label": "water_a",
                        "description": "liquid H-bond: donor + first shells of both waters / acceptor water; acceptor moves",
                    }
                },
            }
            C.dump_json(res_dir / f"cluster_{name}.json", spec)
            index.append({"name": name, "cluster_json": str(res_dir / f"cluster_{name}.json"), **spec["descriptors"]})
            print(f"   {name}: R={hb['R_OO']:.3f} angle={hb['angle']:.1f} env={n_env} waters, {n_atoms} atoms, "
                  f"D(d/a)={dd}/{da} A(d/a)={ad}/{aa}")
    C.dump_json(res_dir / f"index_{args.tag}.json", index)
    print(f"# {len(index)} clusters, index at {res_dir / f'index_{args.tag}.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
