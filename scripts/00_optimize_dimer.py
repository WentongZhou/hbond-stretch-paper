"""Optimize the water dimer with PySCF (RKS + D3/D4) through an ASE BFGS driver.

Usage (inside the WSL venv):
    python scripts/00_optimize_dimer.py --xc revpbe0 --disp d3bj --basis def2-tzvp --tag revpbe0-d3bj_tzvp

Writes geometries/dimer_<tag>.xyz and results/opt_<tag>.json, and prints a
finite-difference check of the analytic gradient along the O...O axis so the
D3 gradient contribution is verified rather than assumed.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def make_mf(symbols, coords, xc, basis, disp, grid_level):
    from pyscf import dft, gto

    mol = gto.M(atom=[(s, tuple(p)) for s, p in zip(symbols, coords)], basis=basis,
                unit="Angstrom", verbose=0)
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = grid_level
    mf.conv_tol = 1e-10
    if disp:
        mf.disp = disp
    return mf


def energy_and_gradient(symbols, coords, xc, basis, disp, grid_level):
    mf = make_mf(symbols, coords, xc, basis, disp, grid_level)
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF did not converge")
    g = mf.nuc_grad_method().kernel()  # Eh/bohr, includes dispersion when mf.disp is set
    return float(e), np.asarray(g)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xc", default="revpbe0")
    parser.add_argument("--disp", default="d3bj")
    parser.add_argument("--basis", default="def2-tzvp")
    parser.add_argument("--grid-level", type=int, default=4)
    parser.add_argument("--tag", default=None)
    parser.add_argument("--fmax", type=float, default=0.003, help="eV/Angstrom")
    parser.add_argument("--steps", type=int, default=200)
    args = parser.parse_args()
    disp = None if args.disp in ("", "none", "None") else args.disp
    tag = args.tag or f"{args.xc}{'-' + disp if disp else ''}_{args.basis.replace('def2-', '')}"

    from ase import Atoms, units
    from ase.calculators.calculator import Calculator, all_changes
    from ase.optimize import BFGS

    class PySCFCalc(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            e, g = energy_and_gradient(atoms.get_chemical_symbols(), atoms.get_positions(),
                                       args.xc, args.basis, disp, args.grid_level)
            self.results["energy"] = e * units.Hartree
            self.results["forces"] = -g * units.Hartree / units.Bohr

    symbols, coords = C.dimer_guess()
    atoms = Atoms(symbols=symbols, positions=coords)
    atoms.calc = PySCFCalc()
    timer = C.Timer()
    opt = BFGS(atoms, logfile=str(C.ROOT / "logs" / f"opt_{tag}.log"))
    opt.run(fmax=args.fmax, steps=args.steps)
    coords = atoms.get_positions()
    e_final = atoms.get_potential_energy() / units.Hartree
    elapsed = timer.lap()

    # Gradient sanity check along the O...O coordinate (central differences).
    h = 0.005
    e_plus, _ = energy_and_gradient(symbols, C.set_oo_distance(coords, C.oo_distance(coords) + h),
                                    args.xc, args.basis, disp, args.grid_level)
    e_minus, _ = energy_and_gradient(symbols, C.set_oo_distance(coords, C.oo_distance(coords) - h),
                                     args.xc, args.basis, disp, args.grid_level)
    fd_grad = (e_plus - e_minus) / (2 * h)  # Eh/Angstrom along R_OO
    _, g = energy_and_gradient(symbols, coords, args.xc, args.basis, disp, args.grid_level)
    axis = coords[C.ACCEPTOR[0]] - coords[C.DONOR[0]]
    axis /= np.linalg.norm(axis)
    an_grad = sum(float(g[i] @ axis) for i in C.ACCEPTOR) / C.BOHR_ANGSTROM  # Eh/Angstrom

    C.GEOM_DIR.mkdir(exist_ok=True)
    C.RESULT_DIR.mkdir(exist_ok=True)
    xyz_path = C.GEOM_DIR / f"dimer_{tag}.xyz"
    C.write_xyz(xyz_path, symbols, coords, f"water dimer {args.xc}{'-' + disp if disp else ''}/{args.basis} E={e_final:.10f} Eh")
    summary = {
        "tag": tag,
        "xc": args.xc,
        "dispersion": disp,
        "basis": args.basis,
        "grid_level": args.grid_level,
        "energy_hartree": e_final,
        "R_OO_angstrom": C.oo_distance(coords),
        "r_OH_bonded_angstrom": float(np.linalg.norm(coords[1] - coords[0])),
        "r_OH_free_angstrom": float(np.linalg.norm(coords[2] - coords[0])),
        "r_OH_acceptor_angstrom": [float(np.linalg.norm(coords[i] - coords[3])) for i in (4, 5)],
        "angle_OHO_deg": float(np.degrees(np.arccos(
            np.dot(coords[0] - coords[1], coords[3] - coords[1])
            / np.linalg.norm(coords[0] - coords[1]) / np.linalg.norm(coords[3] - coords[1])))),
        "bfgs_steps": opt.get_number_of_steps(),
        "fmax_eV_per_angstrom": args.fmax,
        "wall_seconds": elapsed,
        "gradient_check": {
            "analytic_dE_dR_Eh_per_angstrom": an_grad,
            "finite_difference_dE_dR_Eh_per_angstrom": fd_grad,
        },
        "xyz": str(xyz_path),
    }
    C.dump_json(C.RESULT_DIR / f"opt_{tag}.json", summary)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    import json

    raise SystemExit(main())
