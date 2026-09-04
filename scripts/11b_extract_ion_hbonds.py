"""Pick water-water H-bonds between the first and second solvation shell of a
dissolved H3O+ or OH- from the liquid snapshots of 10b_ion_md.py and cut
clusters for the DM-EDA scans (WP6, acid/base in the liquid).

  acid (H3O+): D = first-shell water (O within --shell of the ion O), A = a
               second-shell water accepting from D.  A moves.
  base (OH-):  A = first-shell water, D = a second-shell water donating to A.
               D moves.

Cluster = ion + all first-shell waters + D + A + waters with O within
--env-cut of O_D or O_A (+ the counter-ion if within --env-cut + 1 A of them).
Fragment set "ionliquid": [everything except the moving water] / [moving water].
ct_label is the fragment containing the acceptor O of the target bond.

Usage:
    python scripts/11b_extract_ion_hbonds.py --ion h3o --per-snapshot 1 --seed 13
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

ION_CHARGE = {"H3O": 1, "OHM": -1, "NA": 1, "CL": -1}


def read_snapshot(xyz_path):
    from ase.io import read

    atoms = read(xyz_path, format="extxyz")
    meta = C.load_json(Path(xyz_path).with_suffix(".json"))
    return atoms.get_chemical_symbols(), atoms.get_positions(), float(atoms.cell.lengths()[0]), meta["residues"]


def min_image(v, L):
    return v - L * np.round(v / L)


def water_hbonds(pos, L, waters, r_oo_max, angle_min):
    """H-bonds among water residues; waters = list of [O,H,H] index triples."""
    out = []
    for wd in waters:
        od = wd[0]
        for h in wd[1:]:
            for wa in waters:
                oa = wa[0]
                if oa == od:
                    continue
                d = min_image(pos[oa] - pos[od], L)
                r = float(np.linalg.norm(d))
                if r > r_oo_max:
                    continue
                v1 = min_image(pos[od] - pos[h], L)
                v2 = min_image(pos[oa] - pos[h], L)
                ang = float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))))
                if ang >= angle_min:
                    out.append({"donor_O": od, "H": h, "acceptor_O": oa, "R_OO": r, "angle": ang})
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ion", required=True, choices=["h3o", "oh"])
    parser.add_argument("--snapshots", nargs="+", default=None)
    parser.add_argument("--per-snapshot", type=int, default=1)
    parser.add_argument("--shell", type=float, default=None, help="ion O..O_w cutoff for the first shell (default 3.0 / 3.1)")
    parser.add_argument("--r-oo-max", type=float, default=3.3)
    parser.add_argument("--angle-min", type=float, default=140.0)
    parser.add_argument("--env-cut", type=float, default=3.2)
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()
    ion_res = "H3O" if args.ion == "h3o" else "OHM"
    shell = args.shell or (3.0 if args.ion == "h3o" else 3.1)
    tag = f"liquid_{args.ion}"
    files = args.snapshots or sorted(glob.glob(str(C.GEOM_DIR / tag / "snap_*.xyz")))
    rng = np.random.default_rng(args.seed)
    geom_dir = C.GEOM_DIR / tag
    res_dir = C.RESULT_DIR / "liquid"
    res_dir.mkdir(parents=True, exist_ok=True)

    index = []
    for f in files:
        symbols, pos, L, residues = read_snapshot(f)
        waters = [r["atoms"] for r in residues if r["name"] == "HOH"]
        ion = next(r for r in residues if r["name"] == ion_res)
        counter = [r for r in residues if r["name"] in ("NA", "CL")]
        ion_o = pos[ion["atoms"][0]]
        first = [w for w in waters if np.linalg.norm(min_image(pos[w[0]] - ion_o, L)) <= shell]
        first_o = {w[0] for w in first}
        hbonds = water_hbonds(pos, L, waters, args.r_oo_max, args.angle_min)
        if args.ion == "h3o":
            cands = [hb for hb in hbonds if hb["donor_O"] in first_o and hb["acceptor_O"] not in first_o]
        else:
            cands = [hb for hb in hbonds if hb["acceptor_O"] in first_o and hb["donor_O"] not in first_o]
        snap = Path(f).stem.split("_")[-1]
        print(f"# {Path(f).name}: {len(first)} first-shell waters, {len(cands)} candidate shell-1 => shell-2 bonds")
        if not cands:
            continue
        picks = rng.choice(len(cands), size=min(args.per_snapshot, len(cands)), replace=False)
        for j, p in enumerate(sorted(picks)):
            hb = cands[p]
            od, oa = hb["donor_O"], hb["acceptor_O"]
            wat = {w[0]: w for w in waters}
            D, A = wat[od], wat[oa]
            # environment of D and A
            env = []
            for w in waters:
                if w[0] in (od, oa) or w[0] in first_o:
                    continue
                dmin = min(np.linalg.norm(min_image(pos[w[0]] - pos[od], L)), np.linalg.norm(min_image(pos[w[0]] - pos[oa], L)))
                if dmin <= args.env_cut:
                    env.append((dmin, w))
            env.sort(key=lambda x: x[0])
            shell_waters = [w for w in first if w[0] not in (od, oa)]
            counter_in = []
            for r in counter:
                c = pos[r["atoms"][0]]
                dmin = min(np.linalg.norm(min_image(c - pos[od], L)), np.linalg.norm(min_image(c - pos[oa], L)))
                if dmin <= args.env_cut + 1.0:
                    counter_in.append(r)
            groups = [("D", D), ("A", A), ("ion", ion["atoms"])] + [("shell", w) for w in shell_waters] + \
                     [("env", w) for _, w in env] + [("counter", r["atoms"]) for r in counter_in]
            ref = pos[od]
            coords, syms, labels = [], [], []
            for lab, idx in groups:
                shift = min_image(pos[idx[0]] - ref, L) - (pos[idx[0]] - ref)
                for i in idx:
                    coords.append(pos[i] + shift)
                    syms.append(symbols[i])
                    labels.append(lab)
            coords = np.asarray(coords) - coords[0]
            hb_local = 1 if np.linalg.norm(coords[1] - coords[3]) < np.linalg.norm(coords[2] - coords[3]) else 2
            n_atoms = len(syms)
            moving_atoms = [3, 4, 5] if args.ion == "h3o" else [0, 1, 2]
            rest = [i for i in range(n_atoms) if i not in moving_atoms]
            big_charge = ION_CHARGE[ion_res] + sum(ION_CHARGE[r["name"]] for r in counter_in)
            name = f"{tag}_{snap}_{j}"
            xyz_path = geom_dir / f"cluster_{name}.xyz"
            C.write_xyz(xyz_path, syms, coords, f"{name}: target O0-H{hb_local}...O3, R={hb['R_OO']:.3f} angle={hb['angle']:.1f}, "
                                                 f"{len(shell_waters)} shell + {len(env)} env waters, counter-ion {len(counter_in)}")
            acceptor_label = "water_2nd" if args.ion == "h3o" else "ion_env"
            spec = {
                "system": name, "charge": big_charge, "tag": tag, "xyz": str(xyz_path.relative_to(C.ROOT)).replace("\\", "/"),
                "symbols": syms, "atom_groups": labels,
                "source": {"snapshot": str(f), "donor_O": od, "H": hb["H"], "acceptor_O": oa},
                "descriptors": {"R_OO": hb["R_OO"], "angle_OHO": hb["angle"],
                                "ion_to_donor_O": float(np.linalg.norm(min_image(pos[od] - ion_o, L))),
                                "ion_to_acceptor_O": float(np.linalg.norm(min_image(pos[oa] - ion_o, L))),
                                "n_first_shell": len(first), "n_env_waters": len(env), "n_atoms": n_atoms,
                                "counterion_included": len(counter_in),
                                "counterion_distance": float(min(np.linalg.norm(min_image(pos[r["atoms"][0]] - pos[od], L)) for r in counter)) if counter else None},
                "fragment_sets": {
                    "ionliquid": {
                        "fragments": [{"label": "ion_env", "atoms": rest, "charge": big_charge},
                                      {"label": "water_2nd", "atoms": moving_atoms, "charge": 0}],
                        "moving": "water_2nd", "donor_O": 0, "bonded_H": hb_local, "acceptor_O": 3, "ct_label": acceptor_label,
                        "description": f"{args.ion}: shell-1 => shell-2 water-water H-bond; the second-shell water moves",
                    }
                },
            }
            C.dump_json(res_dir / f"cluster_{name}.json", spec)
            index.append({"name": name, "cluster_json": str(res_dir / f"cluster_{name}.json"), **spec["descriptors"]})
            print(f"   {name}: R={hb['R_OO']:.3f} angle={hb['angle']:.1f} ion-D {spec['descriptors']['ion_to_donor_O']:.2f} "
                  f"ion-A {spec['descriptors']['ion_to_acceptor_O']:.2f}  {n_atoms} atoms, charge {big_charge}")
    C.dump_json(res_dir / f"index_{tag}.json", index)
    print(f"# {len(index)} clusters, index at {res_dir / f'index_{tag}.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
