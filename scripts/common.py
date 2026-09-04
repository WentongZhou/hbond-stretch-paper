"""Shared helpers for the H-bond stretch DM-EDA project.

All scripts import this module.  It wraps ``pyscf_dm_eda`` so that the
intermediate density matrices (promolecule P0, Pauli P_Pauli, final P_S) can
be captured for real-space analysis, and adds IAO fragment charges as a
basis-set-robust complement to the Mulliken charge transfer reported by the
package.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Sequence

import numpy as np

# ----------------------------------------------------------------------------
# Constants (CODATA 2018)
# ----------------------------------------------------------------------------
BOHR_ANGSTROM = 0.529177210903
HARTREE_J = 4.3597447222071e-18
AMU_KG = 1.66053906660e-27
C_CM_S = 2.99792458e10
HARTREE_KCAL = 627.5094740631
M_H = 1.00782503223
M_D = 2.01410177812
M_O = 15.99491461957
M_H2O = M_O + 2 * M_H
M_D2O = M_O + 2 * M_D
MU_H2O = M_H2O / 2.0  # reduced mass of two identical water molecules
MU_D2O = M_D2O / 2.0

ROOT = Path(__file__).resolve().parent.parent
GEOM_DIR = ROOT / "geometries"
RESULT_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures"

# Atom bookkeeping for the dimer: donor O, bonded H, free H, acceptor O, H, H
DONOR = (0, 1, 2)
ACCEPTOR = (3, 4, 5)

# Energy keys in the order used for tables and plots
PRIMARY_KEYS = (
    "Electrostatic Interaction",
    "Exchange Int.",
    "Repulsion",
    "Orbital Relaxation",
    "Correlation Interaction",
    "Dispersion Interaction",
)
DERIVED_KEYS = (
    "Exchange-Repulsion",
    "Corr_Disp",
    "Total Interaction energy",
)
SHORT = {
    "Electrostatic Interaction": "Elec",
    "Exchange Int.": "Exch",
    "Repulsion": "Rep",
    "Exchange-Repulsion": "ExRep",
    "Orbital Relaxation": "OrbRel",
    "Correlation Interaction": "Corr",
    "Dispersion Interaction": "Disp",
    "Corr_Disp": "CorrDisp",
    "Total Interaction energy": "Total",
    "Nuc---Nuc": "NucNuc",
    "1-electron": "OneE",
    "2-electron": "TwoE",
    "Steric": "Steric",
    "Closure Error": "Closure",
}
PRIMARY_SHORT = ("Elec", "Exch", "Rep", "OrbRel", "Corr", "Disp")


# ----------------------------------------------------------------------------
# Units
# ----------------------------------------------------------------------------
def k_to_wavenumber(k_hartree_per_angstrom2: float, mu_amu: float) -> float:
    """Harmonic wavenumber (cm^-1) from a force constant in Eh/Angstrom^2.

    Negative force constants return a negative wavenumber (imaginary mode).
    """
    k_si = k_hartree_per_angstrom2 * HARTREE_J / (1e-10) ** 2
    omega = math.sqrt(abs(k_si) / (mu_amu * AMU_KG)) / (2.0 * math.pi * C_CM_S)
    return math.copysign(omega, k_si)


def k_to_newton_per_metre(k_hartree_per_angstrom2: float) -> float:
    return k_hartree_per_angstrom2 * HARTREE_J / (1e-10) ** 2


# ----------------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------------
def read_xyz(path: str | Path) -> tuple[list[str], np.ndarray]:
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    symbols, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return symbols, np.asarray(coords, dtype=float)


def write_xyz(path: str | Path, symbols: Sequence[str], coords: np.ndarray, comment: str = "") -> None:
    lines = [str(len(symbols)), comment]
    for s, (x, y, z) in zip(symbols, coords):
        lines.append(f"{s:<2s} {x: .8f} {y: .8f} {z: .8f}")
    Path(path).write_text("\n".join(lines) + "\n")


def dimer_guess(r_oh: float = 0.9572, hoh_deg: float = 104.52, r_oo: float = 2.91,
                alpha_deg: float = 5.0, tilt_deg: float = 56.0) -> tuple[list[str], np.ndarray]:
    """Cs trans-linear water dimer guess.

    Donor O at the origin, O...O along +x.  The donor plane is the xy plane;
    the acceptor bisector lies in that plane tilted by ``tilt_deg`` from +x
    on the side opposite to the donor free H (trans-linear arrangement).
    """
    a = math.radians(alpha_deg)
    hoh = math.radians(hoh_deg)
    tilt = math.radians(tilt_deg)
    half = hoh / 2.0
    o_d = np.zeros(3)
    h_b = r_oh * np.array([math.cos(a), -math.sin(a), 0.0])
    h_f = r_oh * np.array([math.cos(hoh - a), math.sin(hoh - a), 0.0])
    o_a = np.array([r_oo, 0.0, 0.0])
    b = np.array([math.cos(tilt), -math.sin(tilt), 0.0])
    z = np.array([0.0, 0.0, 1.0])
    h_a1 = o_a + r_oh * (math.cos(half) * b + math.sin(half) * z)
    h_a2 = o_a + r_oh * (math.cos(half) * b - math.sin(half) * z)
    coords = np.vstack([o_d, h_b, h_f, o_a, h_a1, h_a2])
    return ["O", "H", "H", "O", "H", "H"], coords


def oo_distance(coords: np.ndarray) -> float:
    return float(np.linalg.norm(coords[ACCEPTOR[0]] - coords[DONOR[0]]))


def set_oo_distance(coords: np.ndarray, r_oo: float) -> np.ndarray:
    """Rigidly translate the acceptor water along the O...O axis."""
    axis = coords[ACCEPTOR[0]] - coords[DONOR[0]]
    axis = axis / np.linalg.norm(axis)
    shift = (r_oo - oo_distance(coords)) * axis
    new = coords.copy()
    for i in ACCEPTOR:
        new[i] += shift
    return new


def shift_fragment(coords: np.ndarray, fixed_o: int, moving_o: int, moving_atoms, r_target: float) -> np.ndarray:
    """Rigidly translate ``moving_atoms`` along the fixed_o -> moving_o axis so the O...O distance is ``r_target``."""
    axis = coords[moving_o] - coords[fixed_o]
    r_now = float(np.linalg.norm(axis))
    axis = axis / r_now
    new = coords.copy()
    for i in moving_atoms:
        new[i] += (r_target - r_now) * axis
    return new


def target_from_spec(fs: dict):
    """(fragments, moving_atoms, fixed_o, moving_o, donor_O, bonded_H, acceptor_O, i_donor_frag, i_acceptor_frag)."""
    from pyscf_dm_eda import FragmentSpec

    fragments = [FragmentSpec(tuple(f["atoms"]), int(f["charge"]), 0, f["label"]) for f in fs["fragments"]]
    moving = next(f for f in fs["fragments"] if f["label"] == fs["moving"])
    moving_atoms = list(moving["atoms"])
    od, hb, oa = fs["donor_O"], fs["bonded_H"], fs["acceptor_O"]
    moving_o, fixed_o = (oa, od) if oa in moving_atoms else (od, oa)
    i_don = next(i for i, f in enumerate(fs["fragments"]) if od in f["atoms"])
    i_acc = next(i for i, f in enumerate(fs["fragments"]) if oa in f["atoms"])
    return fragments, moving_atoms, fixed_o, moving_o, od, hb, oa, i_don, i_acc


def water_of(symbols: Sequence[str], coords: np.ndarray, i_o: int) -> list[int]:
    """Atom indices of the water molecule whose oxygen is ``i_o`` (O plus its two nearest H)."""
    d = np.linalg.norm(coords - coords[i_o], axis=1)
    hs = [i for i in np.argsort(d) if symbols[i] == "H"][:2]
    return sorted([int(i_o)] + [int(i) for i in hs])


def pair_fragment_set(fs: dict, symbols: Sequence[str], coords: np.ndarray) -> dict:
    """Fragment set restricted to the two neutral waters of the target bond (same atom numbering)."""
    donor_w = water_of(symbols, coords, fs["donor_O"])
    acc_w = water_of(symbols, coords, fs["acceptor_O"])
    moving_is_acceptor = fs["acceptor_O"] in next(f for f in fs["fragments"] if f["label"] == fs["moving"])["atoms"]
    out = dict(fs)
    out["fragments"] = [{"label": "donor_water", "atoms": donor_w, "charge": 0},
                        {"label": "acceptor_water", "atoms": acc_w, "charge": 0}]
    out["moving"] = "acceptor_water" if moving_is_acceptor else "donor_water"
    out["ct_label"] = "acceptor_water"
    out["description"] = fs.get("description", "") + " [pair only, frozen cluster geometry]"
    return out


def elongate_bond(coords: np.ndarray, i_center: int, i_moving: int, delta: float) -> np.ndarray:
    """Move atom ``i_moving`` away from ``i_center`` along the bond by ``delta`` (Angstrom)."""
    v = coords[i_moving] - coords[i_center]
    v = v / np.linalg.norm(v)
    new = coords.copy()
    new[i_moving] += delta * v
    return new


# ----------------------------------------------------------------------------
# DM-EDA wrapper
# ----------------------------------------------------------------------------
def make_config(xc: str, basis: str, disp: str | None, grid_level: int = 4, **kw):
    from pyscf_dm_eda import SCFConfig

    import os

    max_memory = float(os.environ.get("EDA_MAX_MEMORY_MB", 6000.0))
    return SCFConfig(method=xc, basis=basis, dispersion=disp, grid_level=grid_level,
                     conv_tol=1e-10, max_memory=max_memory, **kw)


def _capturing_eda_class():
    from pyscf_dm_eda import PySCFEDA

    class CapturingEDA(PySCFEDA):
        """PySCFEDA that keeps the SCF states so densities can be re-used."""

        def _decompose(self, super_state, fragment_states):
            self.super_state = super_state
            self.fragment_states = list(fragment_states)
            return super()._decompose(super_state, fragment_states)

    return CapturingEDA


def run_eda(symbols: Sequence[str], coords: np.ndarray, config, fragments=None,
            charge: int | None = None, spin: int | None = None):
    """Run DM-EDA and return (result, eda_object).  Default fragments: water dimer."""
    from pyscf_dm_eda import FragmentSpec

    if fragments is None:
        fragments = [FragmentSpec(DONOR, 0, 0, "donor"), FragmentSpec(ACCEPTOR, 0, 0, "acceptor")]
    atoms = [(s, tuple(float(x) for x in xyz)) for s, xyz in zip(symbols, coords)]
    eda = _capturing_eda_class()(atoms, fragments, config, charge=charge, spin=spin)
    result = eda.run()
    return result, eda


def total_density_matrices(eda) -> dict:
    """Spin-summed P0, P_Pauli and P_S in the common AO basis."""
    from pyscf_dm_eda import eda as _eda

    overlap = eda.super_state.mol.intor_symmetric("int1e_ovlp")
    p0 = np.sum([s.spin_dm for s in eda.fragment_states], axis=0)
    pauli, _ = _eda._build_pauli_density(eda.fragment_states, overlap, eda.config.linear_dep_threshold)
    ps = eda.super_state.spin_dm
    return {
        "P0": np.sum(p0, axis=0),
        "Pauli": np.sum(pauli, axis=0),
        "S": np.sum(ps, axis=0),
        "fragments": [np.sum(s.spin_dm, axis=0) for s in eda.fragment_states],
    }


def rho_on_points(mol, dm_total: np.ndarray, points_angstrom: np.ndarray) -> np.ndarray:
    """Electron density (e/bohr^3) of an AO density matrix at Cartesian points."""
    from pyscf import dft

    coords_bohr = np.asarray(points_angstrom, dtype=float) / BOHR_ANGSTROM
    ao = dft.numint.eval_ao(mol, coords_bohr)
    return dft.numint.eval_rho(mol, ao, dm_total)


def iao_fragment_charges(eda) -> dict:
    """IAO fragment charge transfer with the package's Mulliken sign convention
    (positive = the fragment lost electrons relative to its formal charge)."""
    from pyscf import lo

    mf = eda.super_state.mf
    mol = eda.super_state.mol
    mo_coeff = np.asarray(mf.mo_coeff)
    mo_occ = np.asarray(mf.mo_occ)
    s = mf.get_ovlp()
    if mo_coeff.ndim == 3:  # UKS
        dm = mf.make_rdm1()
        dm_total = dm[0] + dm[1]
        occ_orbs = np.hstack([mo_coeff[0][:, mo_occ[0] > 0], mo_coeff[1][:, mo_occ[1] > 0]])
    else:
        dm_total = mf.make_rdm1()
        occ_orbs = mo_coeff[:, mo_occ > 0]
    iao = lo.iao.iao(mol, occ_orbs)
    iao = lo.vec_lowdin(iao, s)
    pop = np.diag(iao.T @ s @ dm_total @ s @ iao)
    pmol = lo.iao.reference_mol(mol)
    aoslice = pmol.aoslice_by_atom()
    atom_pop = np.array([pop[p0:p1].sum() for (_, _, p0, p1) in aoslice])
    charges = {}
    for frag in eda.fragments:
        z = sum(mol.atom_charge(i) for i in frag.atom_indices)
        n = sum(atom_pop[i] for i in frag.atom_indices)
        charges[frag.label] = float(z - n - frag.charge)
    return charges


def eda_row(result, eda, extra: dict | None = None) -> dict:
    comp = result.components("kcal/mol")
    row = {SHORT.get(k, k): float(v) for k, v in comp.items()}
    row["mulliken_ct"] = {k: float(v) for k, v in result.fragment_charge_transfer.items()}
    try:
        row["iao_ct"] = iao_fragment_charges(eda)
    except Exception as exc:  # IAO is optional
        row["iao_ct"] = {"error": repr(exc)}
    row["E_super_hartree"] = float(eda.super_state.total_energy)
    if extra:
        row.update(extra)
    return row


def dump_json(path: str | Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")


def load_json(path: str | Path):
    return json.loads(Path(path).read_text())


class Timer:
    def __init__(self):
        self.t0 = time.time()

    def lap(self) -> float:
        now = time.time()
        dt = now - self.t0
        self.t0 = now
        return dt
