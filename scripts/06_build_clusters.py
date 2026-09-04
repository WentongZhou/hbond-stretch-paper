"""Build and optimize the WP4 clusters: ion-perturbed water-water H-bonds and a
chain-trimer control.

Systems (fragment labels in brackets, target H-bond marked =>):

  trimer : W0 -> Wd => Wa      open chain built on the optimized dimer; the
                               O0-Od-Oa angle is held (ASE FixInternals) so the
                               chain cannot collapse into the cyclic trimer.
                               Two scan sets: acid_control ([W0+Wd] / Wa, Wa moves)
                               and base_control (W0 / [Wd+Wa], W0 moves).
  acid   : H3O+(W1a,W1b,W1c); W1a => W2     [ion_shell]+ / [water2]; W2 (acceptor) moves
  base   : OH-(W1a..W1d);     W2 => W1a     [ion_shell]- / [water2]; W2 (donor) moves

Usage (inside the WSL venv):
    python scripts/06_build_clusters.py --system acid --xc revpbe0 --disp d3bj \
        --basis def2-tzvp --preopt-basis def2-svp --tag revpbe0-d3bj_tzvp

Writes geometries/cluster_<system>_<tag>.xyz and results/cluster_<system>_<tag>.json
(fragment/scan specification consumed by scripts/07_scan_cluster.py).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

R_OH = 0.965
HOH = 104.5


# ----------------------------------------------------------------------------
# geometry primitives
# ----------------------------------------------------------------------------
def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def perp(v):
    v = unit(v)
    a = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    return unit(np.cross(v, a))


def water_atoms(o, bis, in_plane, r_oh=R_OH, hoh=HOH):
    """O at ``o``; H atoms at +-hoh/2 from the bisector ``bis`` in the plane (bis, in_plane)."""
    o = np.asarray(o, dtype=float)
    bis = unit(bis)
    p = unit(np.asarray(in_plane, dtype=float) - np.dot(in_plane, bis) * bis)
    half = np.radians(hoh / 2.0)
    return [o, o + r_oh * (np.cos(half) * bis + np.sin(half) * p),
            o + r_oh * (np.cos(half) * bis - np.sin(half) * p)]


def acceptor_water(o_d, d, R, tilt_deg=57.0, n=None):
    """Water accepting an H-bond from a donor O at ``o_d`` whose O-H points along ``d``.

    The acceptor bisector is tilted by ``tilt_deg`` from the H-bond axis toward
    ``n`` and its molecular plane is perpendicular to the (d, n) plane, as in
    the trans-linear water dimer.
    """
    d = unit(d)
    n = perp(d) if n is None else unit(np.asarray(n, dtype=float) - np.dot(n, d) * d)
    o = np.asarray(o_d, dtype=float) + R * d
    t = np.radians(tilt_deg)
    return water_atoms(o, np.cos(t) * d + np.sin(t) * n, np.cross(d, n))


def donor_water(o_a, lp, R, p=None, hoh=HOH):
    """Water donating an H-bond to an acceptor O at ``o_a`` from direction ``lp``.

    The bonded H points straight at ``o_a``; the free H lies in the plane
    (lp, p) on the ``p`` side.
    """
    lp = unit(lp)
    p = perp(lp) if p is None else unit(np.asarray(p, dtype=float) - np.dot(p, lp) * lp)
    o = np.asarray(o_a, dtype=float) + R * lp
    h1 = -lp
    h2 = np.cos(np.radians(hoh)) * h1 + np.sin(np.radians(hoh)) * p
    return water_atoms(o, h1 + h2, p, hoh=hoh)


def lone_pair_dirs(o, h1, h2, tilt_deg=55.0):
    b = unit(unit(h1 - o) + unit(h2 - o))
    m = unit(np.cross(h1 - o, h2 - o))
    t = np.radians(tilt_deg)
    return [-np.cos(t) * b + np.sin(t) * m, -np.cos(t) * b - np.sin(t) * m]


def min_distance(coords, new_atoms):
    if len(coords) == 0:
        return np.inf
    d = np.linalg.norm(np.asarray(coords)[:, None, :] - np.asarray(new_atoms)[None, :, :], axis=-1)
    return float(d.min())


def build_trimer():
    symbols, xyz = C.read_xyz(C.GEOM_DIR / "dimer_revpbe0-d3bj_tzvp.xyz")
    o_d, h_b, h_f = xyz[0], xyz[1], xyz[2]
    # third water donating to the donor's lone pair on the +z side of the dimer plane
    lp = lone_pair_dirs(o_d, h_b, h_f)[0]
    w0 = donor_water(o_d, lp, 2.85)
    coords = np.vstack([xyz, w0])
    symbols = list(symbols) + ["O", "H", "H"]
    spec = {
        "system": "trimer",
        "charge": 0,
        "description": "open chain W0 -> Wd -> Wa built on the optimized dimer; O0-Od-Oa angle held fixed",
        "angle_constraint": [6, 0, 3],
        "fragment_sets": {
            "acid_control": {
                "fragments": [{"label": "chain", "atoms": [0, 1, 2, 6, 7, 8], "charge": 0},
                              {"label": "water_a", "atoms": [3, 4, 5], "charge": 0}],
                "moving": "water_a", "donor_O": 0, "bonded_H": 1, "acceptor_O": 3, "ct_label": "water_a",
                "description": "Wd => Wa with W0 donating to Wd (acceptor of the target bond moves)",
            },
            "base_control": {
                "fragments": [{"label": "water_0", "atoms": [6, 7, 8], "charge": 0},
                              {"label": "chain", "atoms": [0, 1, 2, 3, 4, 5], "charge": 0}],
                "moving": "water_0", "donor_O": 6, "bonded_H": 8, "acceptor_O": 0, "ct_label": "chain",
                "description": "W0 => Wd with Wd donating on to Wa (donor of the target bond moves)",
            },
        },
    }
    return symbols, coords, spec


def build_acid():
    o = np.zeros(3)
    z = np.array([0.0, 0.0, 1.0])
    symbols, coords = ["O"], [o]
    us = []
    for phi in (0.0, 120.0, 240.0):
        u = unit([np.cos(np.radians(phi)), np.sin(np.radians(phi)), -0.25])
        us.append(u)
        symbols.append("H")
        coords.append(o + 0.98 * u)
    waters = [acceptor_water(o, u, 2.55, tilt_deg=15.0, n=z) for u in us]
    for w in waters:
        symbols += ["O", "H", "H"]
        coords += w
    w1a = waters[0]
    d = unit(w1a[1] - w1a[0])
    w2 = acceptor_water(w1a[0], d, 2.75, tilt_deg=57.0)
    symbols += ["O", "H", "H"]
    coords += w2
    spec = {
        "system": "acid",
        "charge": 1,
        "description": "Eigen cation H3O+(H2O)3 with one second-shell water accepting from W1a",
        "fragment_sets": {
            "acid": {
                "fragments": [{"label": "ion_shell", "atoms": list(range(0, 13)), "charge": 1},
                              {"label": "water2", "atoms": [13, 14, 15], "charge": 0}],
                "moving": "water2", "donor_O": 4, "bonded_H": 5, "acceptor_O": 13, "ct_label": "water2",
                "description": "W1a => W2, W1a accepts from H3O+ (acceptor of the target bond moves)",
            },
        },
    }
    return symbols, np.asarray(coords), spec


def build_base():
    o = np.zeros(3)
    z = np.array([0.0, 0.0, 1.0])
    symbols, coords = ["O", "H"], [o, o + 0.97 * z]
    waters = []
    for phi in (0.0, 90.0, 180.0, 270.0):
        u = unit([np.cos(np.radians(phi)), np.sin(np.radians(phi)), -0.30])
        p = -z - np.dot(-z, u) * u  # free H points away from the hydroxide H
        waters.append(donor_water(o, u, 2.65, p=p))
    for w in waters:
        symbols += ["O", "H", "H"]
        coords += w
    w1a = waters[0]
    v = unit(w1a[0] - o)                      # ion O -> W1a O
    zp = unit(z - np.dot(z, v) * v)           # "up" component perpendicular to v
    th = np.radians(65.0)                     # O_ion-O_W1a-O_W2 angle = 180 - 65 = 115 deg
    d = np.cos(th) * v + np.sin(th) * zp
    w2 = donor_water(w1a[0], d, 2.85, p=np.cross(v, zp))
    symbols += ["O", "H", "H"]
    coords += w2
    spec = {
        "system": "base",
        "charge": -1,
        "description": "OH-(H2O)4 with one second-shell water donating to W1a; W2 held in the (O_ion, O_W1a, z) plane "
                       "by an O0-O2-O14 angle and H1-O0-O2-O14 dihedral constraint so it cannot bridge to W1b/W1d",
        "angle_constraint": [0, 2, 14],
        "dihedral_constraint": [1, 0, 2, 14],
        "fragment_sets": {
            "base": {
                "fragments": [{"label": "ion_shell", "atoms": list(range(0, 14)), "charge": -1},
                              {"label": "water2", "atoms": [14, 15, 16], "charge": 0}],
                "moving": "water2", "donor_O": 14, "bonded_H": 16, "acceptor_O": 2, "ct_label": "ion_shell",
                "description": "W2 => W1a, W1a donates to OH- (donor of the target bond moves)",
            },
        },
    }
    return symbols, np.asarray(coords), spec


BUILDERS = {"trimer": build_trimer, "acid": build_acid, "base": build_base}


# ----------------------------------------------------------------------------
# PySCF energy/gradient and ASE optimization
# ----------------------------------------------------------------------------
def energy_and_gradient(symbols, coords, charge, xc, basis, disp, grid_level):
    from pyscf import dft, gto

    mol = gto.M(atom=[(s, tuple(p)) for s, p in zip(symbols, coords)], basis=basis, charge=charge,
                unit="Angstrom", verbose=0, max_memory=8000)
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = grid_level
    mf.conv_tol = 1e-10
    if disp:
        mf.disp = disp
    e = mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        e = mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")
    g = mf.nuc_grad_method().kernel()
    return float(e), np.asarray(g)


def optimize(symbols, coords, charge, xc, basis, disp, grid_level, fmax, steps, logfile, angle_constraint=None,
             dihedral_constraint=None):
    from ase import Atoms, units
    from ase.calculators.calculator import Calculator, all_changes
    from ase.constraints import FixInternals
    from ase.optimize import BFGS

    class PySCFCalc(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            e, g = energy_and_gradient(atoms.get_chemical_symbols(), atoms.get_positions(), charge,
                                       xc, basis, disp, grid_level)
            self.results["energy"] = e * units.Hartree
            self.results["forces"] = -g * units.Hartree / units.Bohr

    atoms = Atoms(symbols=symbols, positions=coords)
    angles, dihedrals = [], []
    if angle_constraint:
        i, j, k = angle_constraint
        angles.append([atoms.get_angle(i, j, k), [i, j, k]])
    if dihedral_constraint:
        i, j, k, l = dihedral_constraint
        dihedrals.append([atoms.get_dihedral(i, j, k, l), [i, j, k, l]])
    if angles or dihedrals:
        atoms.set_constraint(FixInternals(angles_deg=angles or None, dihedrals_deg=dihedrals or None))
    atoms.calc = PySCFCalc()
    opt = BFGS(atoms, logfile=str(logfile), maxstep=0.1, trajectory=str(Path(logfile).with_suffix(".traj")))
    opt.run(fmax=fmax, steps=steps)
    e = atoms.get_potential_energy() / units.Hartree
    fmax_final = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1)).max())
    return atoms.get_positions(), e, opt.get_number_of_steps(), fmax_final


def hbond_table(symbols, coords, cutoff=3.3):
    """All O-H...O contacts with R(O..O) < cutoff and angle > 120 deg."""
    o_idx = [i for i, s in enumerate(symbols) if s == "O"]
    h_idx = [i for i, s in enumerate(symbols) if s == "H"]
    rows = []
    for od in o_idx:
        for h in h_idx:
            if np.linalg.norm(coords[h] - coords[od]) > 1.15:
                continue
            for oa in o_idx:
                if oa == od:
                    continue
                r = float(np.linalg.norm(coords[oa] - coords[od]))
                if r > cutoff:
                    continue
                v1, v2 = coords[od] - coords[h], coords[oa] - coords[h]
                ang = float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2), -1.0, 1.0))))
                if ang > 120:
                    rows.append({"donor_O": od, "H": h, "acceptor_O": oa, "R_OO": r,
                                 "r_OH": float(np.linalg.norm(coords[h] - coords[od])), "angle_OHO": ang})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", required=True, choices=sorted(BUILDERS))
    parser.add_argument("--xc", default="revpbe0")
    parser.add_argument("--disp", default="d3bj")
    parser.add_argument("--basis", default="def2-tzvp")
    parser.add_argument("--preopt-basis", default="def2-svp", help="empty string skips the pre-optimization")
    parser.add_argument("--grid-level", type=int, default=4)
    parser.add_argument("--tag", default=None)
    parser.add_argument("--fmax", type=float, default=0.005, help="eV/Angstrom, final basis")
    parser.add_argument("--preopt-fmax", type=float, default=0.03)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--build-only", action="store_true")
    args = parser.parse_args()
    disp = None if args.disp in ("", "none", "None") else args.disp
    tag = args.tag or f"{args.xc}{'-' + disp if disp else ''}_{args.basis.replace('def2-', '')}"

    symbols, coords, spec = BUILDERS[args.system]()
    charge = spec["charge"]
    C.GEOM_DIR.mkdir(exist_ok=True)
    C.RESULT_DIR.mkdir(exist_ok=True)
    (C.ROOT / "logs").mkdir(exist_ok=True)
    guess_path = C.GEOM_DIR / f"cluster_{args.system}_guess.xyz"
    C.write_xyz(guess_path, symbols, coords, f"{args.system} initial guess, charge {charge}")
    print(f"# {args.system}: {len(symbols)} atoms, charge {charge}, guess written to {guess_path}")
    print("# initial H-bonds:")
    for row in hbond_table(symbols, coords):
        print(f"#   O{row['donor_O']}-H{row['H']}...O{row['acceptor_O']}  R={row['R_OO']:.3f}  angle={row['angle_OHO']:.1f}")
    dmin = min(np.linalg.norm(coords[i] - coords[j]) for i in range(len(symbols)) for j in range(i + 1, len(symbols))
               if not (symbols[i] == "H" and symbols[j] == "H"))
    print(f"# minimum non-HH interatomic distance: {dmin:.3f} A")
    if args.build_only:
        return 0

    timer = C.Timer()
    history = []
    stages = []
    if args.preopt_basis:
        stages.append((args.preopt_basis, args.preopt_fmax))
    stages.append((args.basis, args.fmax))
    for basis, fmax in stages:
        log = C.ROOT / "logs" / f"opt_cluster_{args.system}_{tag}_{basis.replace('def2-', '')}.log"
        coords, e, nsteps, fmax_final = optimize(symbols, coords, charge, args.xc, basis, disp, args.grid_level,
                                                  fmax, args.steps, log, spec.get("angle_constraint"),
                                                  spec.get("dihedral_constraint"))
        dt = timer.lap()
        history.append({"basis": basis, "energy_hartree": e, "bfgs_steps": nsteps, "fmax_final_eV_A": fmax_final,
                        "wall_seconds": dt})
        C.write_xyz(C.GEOM_DIR / f"cluster_{args.system}_{tag}_{basis.replace('def2-', '')}.xyz", symbols, coords,
                    f"{args.system} {basis} stage, E={e:.10f} Eh, fmax={fmax_final:.4f} eV/A")
        print(f"# {basis}: E={e:.8f} Eh after {nsteps} steps, fmax={fmax_final:.4f} eV/A, {dt:.0f} s", flush=True)

    xyz_path = C.GEOM_DIR / f"cluster_{args.system}_{tag}.xyz"
    C.write_xyz(xyz_path, symbols, coords,
                f"{args.system} {args.xc}{'-' + disp if disp else ''}/{args.basis} charge {charge} E={history[-1]['energy_hartree']:.10f} Eh")
    hb = hbond_table(symbols, coords)
    targets = {}
    for name, fs in spec["fragment_sets"].items():
        od, h, oa = fs["donor_O"], fs["bonded_H"], fs["acceptor_O"]
        v1, v2 = coords[od] - coords[h], coords[oa] - coords[h]
        targets[name] = {
            "R_OO": float(np.linalg.norm(coords[oa] - coords[od])),
            "r_OH_bonded": float(np.linalg.norm(coords[h] - coords[od])),
            "angle_OHO": float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2), -1.0, 1.0)))),
        }
    out = {
        "tag": tag, "xc": args.xc, "dispersion": disp, "basis": args.basis, "grid_level": args.grid_level,
        "xyz": str(xyz_path), "symbols": symbols, "optimization": history, "hbonds": hb, "targets": targets,
        **spec,
    }
    C.dump_json(C.RESULT_DIR / f"cluster_{args.system}_{tag}.json", out)
    print("# final H-bonds:")
    for row in hb:
        print(f"#   O{row['donor_O']}-H{row['H']}...O{row['acceptor_O']}  R={row['R_OO']:.3f}  r(OH)={row['r_OH']:.4f}  angle={row['angle_OHO']:.1f}")
    print(json.dumps(targets, indent=2))
    print(f"# wrote {xyz_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
