"""PySCF implementation of density-matrix energy decomposition analysis.

This module implements the DM-EDA equations of Zhang et al., using exact J/K
by default and an explicitly labelled PySCF density-fitting approximation when
requested.
(J. Chem. Phys. 160, 174101, 2024; DOI: 10.1063/5.0202787).  It is a standalone
implementation built directly on PySCF density matrices and J/K builders; it
does not parse output from another quantum-chemistry program.

Conventions
-----------
* Energies are stored internally in Hartree.
* Every fragment SCF uses the complete supermolecule AO basis.  Atoms outside
  a fragment are represented as PySCF ghost atoms (counterpoise convention).
* ``spin`` is the signed value ``N_alpha - N_beta``.  For an antiferromagnetic
  broken-symmetry product, use opposite signs for the two fragments.
* Alpha and beta occupied spaces are antisymmetrized independently to form
  the Pauli density.
* DFT contributions not present in the HF-like ``E0[P]`` functional are
  collected in ``Correlation Interaction``.  PySCF can evaluate this residual
  for pure, hybrid, meta-GGA, and range-separated functionals; it closes the
  reported identity algebraically, while its chemical interpretation still
  depends on the chosen functional.
* An empirical D3/D4 correction, when requested and available, is removed
  from the correlation residual and reported exactly once as dispersion.
* Implicit solvation (PySCF C-PCM/IEF-PCM/SS(V)PE/COSMO or SMD) can be
  attached to every SCF.  All gas-phase-operator terms are then evaluated with
  the solvated densities, and the change in the reaction-field energy is
  reported exactly once as ``Desolvation``.  Fragment cavities are built from
  the fragment's real atoms only; ghost atoms never carry a cavity sphere.

PySCF does not provide NBO.  The grid-compatible row therefore reports a
clearly named Mulliken fragment charge transfer instead of mislabelling it as
``NBO_charge``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
import warnings
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
from pyscf import __version__ as pyscf_version
from pyscf import dft, gto, scf
from pyscf.solvent import pcm as pyscf_pcm
from pyscf.solvent import smd as pyscf_smd


HARTREE_TO_KCAL_MOL = 627.5094740631
HARTREE_TO_KJ_MOL = 2625.4996394799
HARTREE_TO_EV = 27.211386245988

_ENERGY_UNIT_FACTORS = {
    "hartree": 1.0,
    "eh": 1.0,
    "au": 1.0,
    "a.u.": 1.0,
    "kcal/mol": HARTREE_TO_KCAL_MOL,
    "kcalmol-1": HARTREE_TO_KCAL_MOL,
    "kcalmol^-1": HARTREE_TO_KCAL_MOL,
    "kj/mol": HARTREE_TO_KJ_MOL,
    "kjmol-1": HARTREE_TO_KJ_MOL,
    "kjmol^-1": HARTREE_TO_KJ_MOL,
    "ev": HARTREE_TO_EV,
}

ENERGY_COMPONENTS = (
    "Total Interaction energy",
    "Electrostatic Interaction",
    "Nuc---Nuc",
    "1-electron",
    "2-electron",
    "Exchange-Repulsion",
    "Exchange Int.",
    "Repulsion",
    "Orbital Relaxation",
    "Correlation Interaction",
    "Dispersion Interaction",
    "Desolvation",
)

GRID_EDA_COLUMNS = (
    "Grid_Index",
    "Tot",
    "Electro",
    "Nuc_Nuc",
    "1e",
    "2e",
    "Exc_Rep",
    "Exc",
    "Rep",
    "Orb_Relax",
    "Corr",
    "Disp",
    "Desolv",
    "Corr_Disp",
    "Steric",
    "Mulliken_CT",
    "Closure_Error",
)

_HF_METHODS = {"hf", "rhf", "uhf"}
# Implicit-solvent models built on PySCF's SWIG-discretised PCM surface.  The
# canonical value is the ``method`` string understood by ``pyscf.solvent.pcm``;
# ``"SMD"`` is handled by ``pyscf.solvent.smd`` (PCM electrostatics + CDS).
_SOLVENT_MODEL_ALIASES = {
    "pcm": "C-PCM",
    "cpcm": "C-PCM",
    "c-pcm": "C-PCM",
    "iefpcm": "IEF-PCM",
    "ief-pcm": "IEF-PCM",
    "cosmo": "COSMO",
    "ssvpe": "SS(V)PE",
    "ss(v)pe": "SS(V)PE",
    "smd": "SMD",
}
# ddCOSMO/ddPCM discretise the cavity per atom of ``mol`` and evaluate the
# density on atom-centred DFT grids, so their cavity cannot be restricted to
# the real atoms of a ghost-padded fragment molecule.
_UNSUPPORTED_SOLVENT_MODELS = frozenset({"ddcosmo", "ddpcm", "pe", "cosmors"})
# Attributes of the PySCF solvent object that this module manages itself or
# that are internal state; ``solvent_options`` must not override them.
_RESERVED_SOLVENT_OPTIONS = frozenset(
    {"mol", "method", "eps", "solvent", "e", "v", "surface", "v_grids_n",
     "cavity_mol", "frozen", "state_id", "stdout", "verbose", "max_memory"}
)
# Libxc/PySCF expose only the XC part of these functionals; their empirical
# dispersion tail is part of the parameterisation and cannot be supplied later
# via ``dispersion=...`` with the same parameters.  Reject them instead of
# silently running "half a functional" with ``mf.disp`` disabled.
_BUILTIN_DISPERSION_METHODS = frozenset({"wb97x-d", "b97-d"})
_CHECKPOINT_PREFIX = "# pyscf-dm-eda-grid-v1 fingerprint="


class EDAError(RuntimeError):
    """Base error for invalid or unsuccessful EDA calculations."""


class SCFConvergenceError(EDAError):
    """Raised when a fragment or supermolecule SCF does not converge."""


class IncompatibleFragmentError(EDAError):
    """Raised when fragment AO or spin spaces cannot form the supermolecule."""


class EDAValidationError(EDAError):
    """Raised when numerical validation fails; carries the full result.

    ``result`` is the complete :class:`EDAResult` (components, charge transfer,
    diagnostics and metadata) so callers can inspect every quantity instead of
    having to re-run an expensive calculation.
    """

    def __init__(self, message: str, result: EDAResult) -> None:
        super().__init__(message)
        self.result = result


@dataclass(frozen=True)
class Atom:
    """An XYZ atom in Cartesian coordinates."""

    symbol: str
    xyz: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not self.symbol or self.symbol.lower().startswith("ghost"):
            raise ValueError("Input atoms must be real element labels, not ghosts")
        try:
            xyz = tuple(float(value) for value in self.xyz)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid Cartesian coordinate for {self.symbol}: {self.xyz}"
            ) from exc
        if len(xyz) != 3 or not np.all(np.isfinite(xyz)):
            raise ValueError(f"Invalid Cartesian coordinate for {self.symbol}: {self.xyz}")
        object.__setattr__(self, "xyz", xyz)


@dataclass(frozen=True)
class FragmentSpec:
    """Atoms, charge, and signed spin projection of one fragment."""

    atom_indices: tuple[int, ...]
    charge: int = 0
    spin: int = 0
    label: str = ""

    def __post_init__(self) -> None:
        if any(
            isinstance(index, (bool, np.bool_))
            or not isinstance(index, (int, np.integer))
            for index in self.atom_indices
        ):
            raise TypeError("Fragment atom indices must be integers")
        indices = tuple(int(index) for index in self.atom_indices)
        object.__setattr__(self, "atom_indices", indices)
        if not indices:
            raise ValueError("A fragment must contain at least one atom")
        if len(set(indices)) != len(indices):
            raise ValueError(f"A fragment contains duplicate atom indices: {indices}")
        if (
            isinstance(self.charge, (bool, np.bool_))
            or isinstance(self.spin, (bool, np.bool_))
            or not isinstance(self.charge, (int, np.integer))
            or not isinstance(self.spin, (int, np.integer))
        ):
            raise TypeError("Fragment charge and spin must be integers")
        object.__setattr__(self, "charge", int(self.charge))
        object.__setattr__(self, "spin", int(self.spin))


@dataclass(frozen=True)
class SCFConfig:
    """Electronic-structure controls shared by all CP calculations."""

    method: str = "r2scan"
    basis: Any = "def2-svp"
    ecp: Any = None
    dispersion: str | None = None
    grid_level: int = 4
    conv_tol: float = 1e-9
    max_cycle: int = 100
    max_memory: float = 4000.0
    density_fit: bool = False
    auxbasis: Any = None
    unrestricted: bool | None = None
    level_shift: float = 0.0
    damp: float = 0.0
    init_guess: str = "minao"
    newton_fallback: bool = True
    linear_dep_threshold: float = 1e-9
    validation_tol: float = 1e-7
    strict_validation: bool = True
    verbose: int = 0
    unit: str = "Angstrom"
    # Implicit solvation through PySCF.  ``solvent`` selects the model
    # (C-PCM, IEF-PCM, SS(V)PE, COSMO, or SMD; ``"pcm"`` means C-PCM).
    # ``solvent_eps`` overrides the dielectric constant; ``solvent_name`` picks
    # a solvent from PySCF's SMD database (required descriptors for SMD, or the
    # dielectric constant for PCM models).  ``solvent_options`` are extra
    # attributes set on the PySCF solvent object (e.g. ``lebedev_order``).
    solvent: str | None = None
    solvent_eps: float | None = None
    solvent_name: str | None = None
    solvent_options: Mapping[str, Any] | None = None
    # Electrostatic embedding: ``(coords_angstrom, charges)`` of external point
    # charges felt by every SCF (supermolecule and fragments).  The
    # decomposition keeps the fragment-fragment cross terms free of the
    # embedding potential, while the frozen, Pauli and relaxed densities are
    # evaluated with it, so closure and the frozen identity still hold exactly.
    point_charges: Any = None

    def __post_init__(self) -> None:
        if not self.method.strip():
            raise ValueError("SCF method/XC functional cannot be empty")
        if self.point_charges is not None:
            try:
                coords, charges = self.point_charges
                coords = np.asarray(coords, dtype=float).reshape(-1, 3)
                charges = np.asarray(charges, dtype=float).reshape(-1)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "point_charges must be (coords_angstrom[n,3], charges[n])"
                ) from exc
            if len(coords) != len(charges):
                raise ValueError("point_charges: one charge per coordinate is required")
            object.__setattr__(self, "point_charges", (coords, charges))
            if self.solvent is not None:
                raise ValueError("point_charges cannot be combined with an implicit solvent")
        for name in ("grid_level", "max_cycle", "verbose"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(
                value, (int, np.integer)
            ):
                raise TypeError(f"{name} must be an integer")
        if not isinstance(self.density_fit, (bool, np.bool_)):
            raise TypeError("density_fit must be a boolean")
        if not isinstance(self.strict_validation, (bool, np.bool_)):
            raise TypeError("strict_validation must be a boolean")
        if self.unrestricted is not None and not isinstance(
            self.unrestricted, (bool, np.bool_)
        ):
            raise TypeError("unrestricted must be True, False, or None")
        if self.grid_level < 0:
            raise ValueError("grid_level must be non-negative")
        if self.conv_tol <= 0 or self.max_cycle <= 0:
            raise ValueError("conv_tol and max_cycle must be positive")
        if self.linear_dep_threshold <= 0:
            raise ValueError("linear_dep_threshold must be positive")
        if self.validation_tol <= 0:
            raise ValueError("validation_tol must be positive")
        method_lower = self.method.strip().lower()
        if "-d3" in method_lower or "-d4" in method_lower:
            raise ValueError(
                "Pass the base XC functional in method and the D3/D4 model "
                "separately through dispersion"
            )
        if method_lower in _BUILTIN_DISPERSION_METHODS:
            raise ValueError(
                f"{self.method!r} has empirical dispersion built into its "
                "parameterisation. PySCF/libxc provides only the XC part, and "
                "this implementation keeps mf.disp disabled, so running it "
                "would silently drop the dispersion tail. Use the underlying "
                "XC functional together with an explicit dispersion model "
                "instead."
            )
        if method_lower == "uhf" and self.unrestricted is False:
            raise ValueError("method='uhf' conflicts with unrestricted=False")
        if method_lower == "rhf" and self.unrestricted is True:
            raise ValueError("method='rhf' conflicts with unrestricted=True")
        _validate_solvent_config(self)


@dataclass
class SCFState:
    """Numerical state needed by the density-matrix decomposition."""

    label: str
    mol: Any
    mf: Any
    spin_dm: np.ndarray
    occupied_coeff: tuple[np.ndarray, np.ndarray]
    total_energy: float
    dispersion_energy: float
    # Reaction-field energy of the converged density in this state's own
    # cavity (electrostatic PCM part plus CDS for SMD); 0 without solvent.
    solvent_energy: float = 0.0
    solvent_cds_energy: float = 0.0
    # PySCF's own ``scf_summary`` value, kept only as a consistency diagnostic.
    solvent_scf_summary_energy: float | None = None


@dataclass(frozen=True)
class EDAResult:
    """A complete, closed DM-EDA result in atomic units."""

    components_hartree: Mapping[str, float]
    fragment_charge_transfer: Mapping[str, float]
    diagnostics: Mapping[str, Any]
    metadata: Mapping[str, Any]

    def components(self, unit: str = "kcal/mol") -> dict[str, float]:
        """Return all primary and derived energy components in ``unit``."""

        normalized = unit.strip().lower().replace(" ", "")
        try:
            factor = _ENERGY_UNIT_FACTORS[normalized]
        except KeyError as exc:
            raise ValueError(f"Unsupported energy unit: {unit!r}") from exc

        values = {
            key: float(self.components_hartree[key]) * factor
            for key in ENERGY_COMPONENTS
        }
        values["Corr_Disp"] = (
            values["Correlation Interaction"] + values["Dispersion Interaction"]
        )
        values["Steric"] = values["Exchange-Repulsion"] + values["Corr_Disp"]
        values["Closure Error"] = float(self.diagnostics["closure_error_hartree"]) * factor
        return values

    def grid_row(
        self,
        grid_index: int,
        probe_fragment: str | None = None,
        unit: str = "kcal/mol",
    ) -> dict[str, float | int]:
        """Map the result to the compact batch-grid column order."""

        values = self.components(unit)
        if probe_fragment is None:
            probe_fragment = list(self.fragment_charge_transfer)[-1]
        try:
            charge_transfer = float(self.fragment_charge_transfer[probe_fragment])
        except KeyError as exc:
            raise KeyError(f"Unknown probe fragment label: {probe_fragment!r}") from exc

        return {
            "Grid_Index": int(grid_index),
            "Tot": values["Total Interaction energy"],
            "Electro": values["Electrostatic Interaction"],
            "Nuc_Nuc": values["Nuc---Nuc"],
            "1e": values["1-electron"],
            "2e": values["2-electron"],
            "Exc_Rep": values["Exchange-Repulsion"],
            "Exc": values["Exchange Int."],
            "Rep": values["Repulsion"],
            "Orb_Relax": values["Orbital Relaxation"],
            "Corr": values["Correlation Interaction"],
            "Disp": values["Dispersion Interaction"],
            "Desolv": values["Desolvation"],
            "Corr_Disp": values["Corr_Disp"],
            "Steric": values["Steric"],
            "Mulliken_CT": charge_transfer,
            "Closure_Error": values["Closure Error"],
        }

    def as_dict(self, unit: str = "kcal/mol") -> dict[str, Any]:
        """Return a JSON-serializable representation."""

        return {
            "energy_unit": unit,
            "components": self.components(unit),
            "fragment_charge_transfer_e": dict(self.fragment_charge_transfer),
            "diagnostics": _jsonable(self.diagnostics),
            "metadata": _jsonable(self.metadata),
        }

    def write_json(self, path: str | os.PathLike[str], unit: str = "kcal/mol") -> Path:
        """Write a stable, atomically replaced machine-readable result file."""

        output = Path(path)
        _atomic_write_text(
            output,
            json.dumps(self.as_dict(unit), indent=2, sort_keys=True) + "\n",
        )
        return output


@dataclass(frozen=True)
class GridRunSummary:
    """Paths and counts produced by :class:`PySCFGridRunner`."""

    energy_table: Path
    extended_xyz: Path
    metadata_file: Path
    completed_points: int


class PySCFEDA:
    """Run counterpoise-corrected density-matrix EDA with PySCF."""

    def __init__(
        self,
        atoms: Sequence[Atom | Sequence[Any]],
        fragments: Sequence[FragmentSpec],
        config: SCFConfig | None = None,
        *,
        charge: int | None = None,
        spin: int | None = None,
    ) -> None:
        self.atoms = tuple(_coerce_atom(atom) for atom in atoms)
        self.fragments = _normalize_fragments(fragments, len(self.atoms))
        self.config = config or SCFConfig()

        fragment_charge = sum(fragment.charge for fragment in self.fragments)
        fragment_spin = sum(fragment.spin for fragment in self.fragments)
        self.charge = (
            fragment_charge if charge is None else _strict_int(charge, "charge")
        )
        self.spin = fragment_spin if spin is None else _strict_int(spin, "spin")

        if self.charge != fragment_charge:
            raise IncompatibleFragmentError(
                "Supermolecule charge must equal the sum of fragment charges "
                f"({self.charge} != {fragment_charge})"
            )
        if self.spin != fragment_spin:
            raise IncompatibleFragmentError(
                "The determinant implementation requires supermolecule spin "
                "N_alpha-N_beta to equal the signed sum of fragment spins. "
                "For broken-symmetry coupling, give one fragment a negative spin "
                f"({self.spin} != {fragment_spin})."
            )

    @classmethod
    def from_xyz(
        cls,
        path: str | os.PathLike[str],
        fragments: Sequence[FragmentSpec],
        config: SCFConfig | None = None,
        *,
        charge: int | None = None,
        spin: int | None = None,
    ) -> "PySCFEDA":
        """Construct an EDA calculation from an XYZ file."""

        return cls(read_xyz(path), fragments, config, charge=charge, spin=spin)

    def run(self) -> EDAResult:
        """Run all SCFs and evaluate the DM-EDA equations."""

        super_mol = _build_molecule(
            self.atoms,
            self.config,
            charge=self.charge,
            spin=self.spin,
            real_indices=range(len(self.atoms)),
        )
        _validate_ao_overlap(super_mol, self.config, "supermolecule")

        fragment_states: list[SCFState] = []
        for fragment in self.fragments:
            fragment_mol = _build_molecule(
                self.atoms,
                self.config,
                charge=fragment.charge,
                spin=fragment.spin,
                real_indices=fragment.atom_indices,
            )
            physical_fragment_mol = _build_molecule(
                tuple(self.atoms[index] for index in fragment.atom_indices),
                self.config,
                charge=fragment.charge,
                spin=fragment.spin,
                real_indices=range(len(fragment.atom_indices)),
            )
            fragment_states.append(
                _run_scf(
                    fragment_mol,
                    self.config,
                    fragment.label,
                    physical_mol=physical_fragment_mol,
                    force_unrestricted=self.config.unrestricted is True,
                )
            )

        method_lower = self.config.method.lower()
        if method_lower == "uhf":
            use_unrestricted_super = True
        elif method_lower == "rhf":
            use_unrestricted_super = False
        else:
            use_unrestricted_super = (
                self.config.unrestricted
                if self.config.unrestricted is not None
                else self.spin != 0
                or any(fragment.spin != 0 for fragment in self.fragments)
            )
        if not use_unrestricted_super and self.spin != 0:
            raise IncompatibleFragmentError(
                "A restricted reference cannot represent a nonzero spin"
            )
        promolecule_guess = np.sum(
            [state.spin_dm for state in fragment_states], axis=0
        )
        super_state = _run_scf(
            super_mol,
            self.config,
            "supermolecule",
            dm0=(
                promolecule_guess
                if use_unrestricted_super
                else np.sum(promolecule_guess, axis=0)
            ),
            physical_mol=super_mol,
            force_unrestricted=bool(use_unrestricted_super),
        )

        return self._decompose(super_state, fragment_states)

    def _decompose(
        self,
        super_state: SCFState,
        fragment_states: Sequence[SCFState],
    ) -> EDAResult:
        super_mol = super_state.mol
        overlap = super_mol.intor_symmetric("int1e_ovlp")
        kinetic = super_mol.intor_symmetric("int1e_kin")
        # ``super_core``/``fragment_core`` include the embedding potential when
        # point charges are present (they enter E0[P] for every density); the
        # gas-phase cores below define the fragment-fragment cross terms.
        super_core = np.asarray(super_state.mf.get_hcore(super_mol))
        embedded = self.config.point_charges is not None
        super_core_gas = (
            np.asarray(scf.hf.get_hcore(super_mol)) if embedded else super_core
        )

        for state in fragment_states:
            _assert_same_ao_space(super_mol, state.mol, overlap, kinetic, state.label)

        fragment_core = [
            np.asarray(state.mf.get_hcore(state.mol)) for state in fragment_states
        ]
        fragment_core_gas = (
            [np.asarray(scf.hf.get_hcore(state.mol)) for state in fragment_states]
            if embedded
            else fragment_core
        )
        fragment_potential = [core - kinetic for core in fragment_core_gas]
        core_partition_error = _max_abs(
            super_core_gas - kinetic - np.sum(fragment_potential, axis=0)
        )
        if core_partition_error > self.config.validation_tol:
            raise IncompatibleFragmentError(
                "Fragment one-electron operators do not reconstruct the "
                f"supermolecule core Hamiltonian (max error {core_partition_error:.3e})."
            )

        p0 = np.sum([state.spin_dm for state in fragment_states], axis=0)
        pauli_dm, pauli_diagnostics = _build_pauli_density(
            fragment_states,
            overlap,
            self.config.linear_dep_threshold,
        )

        # One batched J/K call shares the same two-electron integral build for
        # every fragment and intermediate density.
        spin_densities: list[np.ndarray] = []
        fragment_jk_indices: list[tuple[int, int]] = []
        for state in fragment_states:
            first = len(spin_densities)
            spin_densities.extend((state.spin_dm[0], state.spin_dm[1]))
            fragment_jk_indices.append((first, first + 1))
        pauli_indices = (len(spin_densities), len(spin_densities) + 1)
        spin_densities.extend((pauli_dm[0], pauli_dm[1]))
        super_indices = (len(spin_densities), len(spin_densities) + 1)
        spin_densities.extend((super_state.spin_dm[0], super_state.spin_dm[1]))
        vj, vk = _batch_jk(super_state, spin_densities)

        fragment_e0: list[float] = []
        for state, core, indices in zip(
            fragment_states, fragment_core, fragment_jk_indices
        ):
            fragment_e0.append(
                _e0(state.spin_dm, core, vj[list(indices)], vk[list(indices)])
            )

        p0_vj = np.stack(
            [
                np.sum([vj[indices[spin]] for indices in fragment_jk_indices], axis=0)
                for spin in (0, 1)
            ]
        )
        p0_vk = np.stack(
            [
                np.sum([vk[indices[spin]] for indices in fragment_jk_indices], axis=0)
                for spin in (0, 1)
            ]
        )
        e0_p0 = _e0(p0, super_core, p0_vj, p0_vk)
        e0_pauli = _e0(
            pauli_dm,
            super_core,
            vj[list(pauli_indices)],
            vk[list(pauli_indices)],
        )
        e0_super = _e0(
            super_state.spin_dm,
            super_core,
            vj[list(super_indices)],
            vk[list(super_indices)],
        )

        nuclear = float(super_mol.energy_nuc()) - sum(
            float(state.mol.energy_nuc()) for state in fragment_states
        )
        electron_nuclear = 0.0
        for density_index, density_state in enumerate(fragment_states):
            density = np.sum(density_state.spin_dm, axis=0)
            for potential_index, potential in enumerate(fragment_potential):
                if density_index != potential_index:
                    electron_nuclear += _trace_product(density, potential)

        coulomb = 0.0
        exchange = 0.0
        for left in range(len(fragment_states)):
            left_state = fragment_states[left]
            left_total = np.sum(left_state.spin_dm, axis=0)
            for right in range(left + 1, len(fragment_states)):
                right_indices = fragment_jk_indices[right]
                right_j = vj[right_indices[0]] + vj[right_indices[1]]
                coulomb += _trace_product(left_total, right_j)
                exchange -= sum(
                    _trace_product(left_state.spin_dm[sigma], vk[right_indices[sigma]])
                    for sigma in (0, 1)
                )

        electrostatic = nuclear + electron_nuclear + coulomb
        repulsion = e0_pauli - e0_p0
        polarization = e0_super - e0_pauli

        dispersion = super_state.dispersion_energy - sum(
            state.dispersion_energy for state in fragment_states
        )
        # The reaction-field energy is part of ``mf.e_tot`` but not of the
        # gas-phase functional ``E0[P]``; take it out of the residual and
        # report its change exactly once as desolvation.
        corr_super = (
            super_state.total_energy
            - float(super_mol.energy_nuc())
            - super_state.dispersion_energy
            - super_state.solvent_energy
            - e0_super
        )
        corr_fragments = [
            state.total_energy
            - float(state.mol.energy_nuc())
            - state.dispersion_energy
            - state.solvent_energy
            - e0
            for state, e0 in zip(fragment_states, fragment_e0)
        ]
        correlation = corr_super - sum(corr_fragments)
        desolvation = super_state.solvent_energy - sum(
            state.solvent_energy for state in fragment_states
        )

        total_interaction = super_state.total_energy - sum(
            state.total_energy for state in fragment_states
        )
        exchange_repulsion = exchange + repulsion
        component_sum = (
            electrostatic
            + exchange
            + repulsion
            + polarization
            + correlation
            + dispersion
            + desolvation
        )
        closure_error = total_interaction - component_sum
        frozen_identity_error = (
            e0_p0
            + float(super_mol.energy_nuc())
            - sum(
                e0 + float(state.mol.energy_nuc())
                for e0, state in zip(fragment_e0, fragment_states)
            )
            - electrostatic
            - exchange
        )

        n_point_charges = (
            int(len(self.config.point_charges[1])) if embedded else 0
        )
        components = {
            "Total Interaction energy": total_interaction,
            "Electrostatic Interaction": electrostatic,
            "Nuc---Nuc": nuclear,
            "1-electron": electron_nuclear,
            "2-electron": coulomb,
            "Exchange-Repulsion": exchange_repulsion,
            "Exchange Int.": exchange,
            "Repulsion": repulsion,
            "Orbital Relaxation": polarization,
            "Correlation Interaction": correlation,
            "Dispersion Interaction": dispersion,
            "Desolvation": desolvation,
        }

        solvent_diagnostics = _solvent_diagnostics(
            super_state, fragment_states, p0, pauli_dm, desolvation
        )
        fragment_charge_transfer = _mulliken_fragment_charge_transfer(
            super_state.spin_dm,
            super_mol,
            overlap,
            self.fragments,
        )
        super_electrons = _electron_count(super_state.spin_dm, overlap)
        promolecule_electrons = _electron_count(p0, overlap)
        pauli_electrons = _electron_count(pauli_dm, overlap)
        fragment_idempotency = {
            state.label: _idempotency_error(state.spin_dm, overlap)
            for state in fragment_states
        }
        fragment_hermiticity = {
            state.label: _hermiticity_error(state.spin_dm)
            for state in fragment_states
        }
        electron_count_error = max(
            abs(super_electrons["alpha"] - super_mol.nelec[0]),
            abs(super_electrons["beta"] - super_mol.nelec[1]),
            abs(
                promolecule_electrons["alpha"]
                - sum(state.mol.nelec[0] for state in fragment_states)
            ),
            abs(
                promolecule_electrons["beta"]
                - sum(state.mol.nelec[1] for state in fragment_states)
            ),
            abs(
                pauli_electrons["alpha"]
                - sum(state.mol.nelec[0] for state in fragment_states)
            ),
            abs(
                pauli_electrons["beta"]
                - sum(state.mol.nelec[1] for state in fragment_states)
            ),
        )
        diagnostics = {
            "closure_error_hartree": closure_error,
            "frozen_identity_error_hartree": frozen_identity_error,
            "n_point_charges": n_point_charges,
            "core_partition_max_abs": core_partition_error,
            "ao_overlap_condition": float(np.linalg.cond(overlap)),
            "electron_count_max_abs_error": electron_count_error,
            "supermolecule_electrons_from_density": super_electrons,
            "promolecule_electrons_from_density": promolecule_electrons,
            "pauli_electrons_from_density": pauli_electrons,
            "fragment_scf_converged": {
                state.label: bool(state.mf.converged) for state in fragment_states
            },
            "supermolecule_scf_converged": bool(super_state.mf.converged),
            "supermolecule_idempotency_max_abs": _idempotency_error(
                super_state.spin_dm, overlap
            ),
            "pauli_idempotency_max_abs": _idempotency_error(pauli_dm, overlap),
            "promolecule_hermiticity_max_abs": _hermiticity_error(p0),
            "supermolecule_hermiticity_max_abs": _hermiticity_error(
                super_state.spin_dm
            ),
            "pauli_hermiticity_max_abs": _hermiticity_error(pauli_dm),
            "fragment_idempotency_max_abs": fragment_idempotency,
            "fragment_hermiticity_max_abs": fragment_hermiticity,
            "solvent": solvent_diagnostics,
            **pauli_diagnostics,
        }
        validation_errors = {
            "closure": abs(closure_error),
            "frozen identity": abs(frozen_identity_error),
            "core partition": core_partition_error,
            "electron count": electron_count_error,
            "supermolecule idempotency": diagnostics[
                "supermolecule_idempotency_max_abs"
            ],
            "Pauli idempotency": diagnostics["pauli_idempotency_max_abs"],
            "density hermiticity": max(
                diagnostics["promolecule_hermiticity_max_abs"],
                diagnostics["supermolecule_hermiticity_max_abs"],
                diagnostics["pauli_hermiticity_max_abs"],
                *fragment_hermiticity.values(),
            ),
            "fragment idempotency": max(fragment_idempotency.values()),
        }
        diagnostics["validation_errors"] = validation_errors
        worst_name = max(validation_errors, key=validation_errors.get)
        worst_error = validation_errors[worst_name]
        diagnostics["validation_worst_case"] = worst_name
        diagnostics["validation_tol"] = self.config.validation_tol
        metadata = {
            "scheme": (
                "DM-EDA (PySCF density-fitted J/K approximation)"
                if self.config.density_fit
                else "DM-EDA (exact J/K)"
            ),
            "jk_evaluation": "density-fitted" if self.config.density_fit else "exact",
            "reference_doi": "10.1063/5.0202787",
            "pyscf_version": pyscf_version,
            "method": self.config.method,
            "basis": repr(self.config.basis),
            "ecp": repr(self.config.ecp),
            "dispersion": self.config.dispersion,
            "density_fit": self.config.density_fit,
            "solvent": _solvent_metadata(self.config, super_state),
            "strict_validation": bool(self.config.strict_validation),
            "validation_tol": self.config.validation_tol,
            "unrestricted_supermolecule": np.asarray(
                super_state.mf.mo_coeff
            ).ndim == 3,
            "supermolecule_charge": self.charge,
            "supermolecule_spin": self.spin,
            "fragments": [
                {
                    "label": fragment.label,
                    "atom_indices_zero_based": list(fragment.atom_indices),
                    "charge": fragment.charge,
                    "spin": fragment.spin,
                }
                for fragment in self.fragments
            ],
            "charge_analysis": "Mulliken",
            "limitations": [
                (
                    "Implicit solvent: Desolvation is the change of the PySCF "
                    "reaction-field energy (plus CDS for SMD) between the "
                    "supermolecule cavity and the isolated real-atom fragment "
                    "cavities; all other terms use gas-phase operators on the "
                    "solvated densities. No explicit solvent, no thermal "
                    "corrections."
                    if self.config.solvent
                    else "Gas phase only; no PCM/desolvation term."
                ),
                "No NBO analysis; Mulliken charge transfer is reported separately.",
                "DM-EDA definition; not a numerical clone of other EDA schemes.",
                "Fragments use the full ghost basis (counterpoise convention).",
                "Corr excludes empirical D3/D4; Corr_Disp is the paper's "
                "combined correlation term.",
            ],
        }
        result = EDAResult(
            components, fragment_charge_transfer, diagnostics, metadata
        )
        if worst_error > self.config.validation_tol:
            message = (
                f"DM-EDA numerical validation failed: {worst_name} error "
                f"{worst_error:.3e} exceeds {self.config.validation_tol:.3e}"
            )
            if not self.config.strict_validation:
                warnings.warn(message, RuntimeWarning, stacklevel=2)
                return result
            raise EDAValidationError(message, result)
        return result


@dataclass
class PySCFGridRunner:
    """Evaluate a directory of already positioned probe geometries.

    The runner consumes ``<name>_filtered.xyz`` and
    ``<name>_probe/mol_probe_<i>.xyz`` without depending on the program that
    generated them. Results are checkpointed one grid point at a time as TSV;
    failed calculations raise and never write a synthetic zero-energy row.

    The checkpoint writer currently rewrites all completed rows for every
    finished point (atomic replace with ``fsync``), so I/O grows quadratically
    with the number of points. This is fine up to a few thousand points; for
    grids of 10⁴ points or more, either use ``indices=...`` to shard the job
    into several outputs or expect a longer checkpoint commit tail.

    ``indices`` selects a zero-based subset of grid points for one shard; the
    main process and each shard must write to different output paths (the lock
    is per output file). ``progress_callback`` is called as
    ``callback(index, total_selected)`` after each successfully checkpointed
    point.

    For linear polyatomic probes whose template origin lies off the molecular
    axis, the Kabsch rotation has an unphysical free rotation around the axis;
    use ``probe_anchor_atom`` for such cases (see ``_validate_grid_geometry``).
    """

    molecule_xyz: str | os.PathLike[str]
    molecule_charge: int
    molecule_spin: int
    probe_charge: int
    probe_spin: int
    config: SCFConfig = field(default_factory=SCFConfig)
    supermolecule_charge: int | None = None
    supermolecule_spin: int | None = None
    probe_directory: str | os.PathLike[str] | None = None
    combined_pattern: str = "mol_probe_{index}.xyz"
    grid_xyz: str | os.PathLike[str] | None = None
    energy_output: str | os.PathLike[str] | None = None
    xyz_output: str | os.PathLike[str] | None = None
    geometry_tolerance: float = 1e-4
    probe_anchor_atom: int | None = None
    probe_template_xyz: str | os.PathLike[str] | None = None
    indices: Sequence[int] | None = None
    progress_callback: Callable[[int, int], None] | None = field(
        default=None, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        for name in (
            "molecule_charge",
            "molecule_spin",
            "probe_charge",
            "probe_spin",
        ):
            setattr(self, name, _strict_int(getattr(self, name), name))
        for name in ("supermolecule_charge", "supermolecule_spin"):
            value = getattr(self, name)
            if value is not None:
                setattr(self, name, _strict_int(value, name))
        if not isinstance(self.config, SCFConfig):
            raise TypeError("config must be an SCFConfig")
        if self.geometry_tolerance <= 0:
            raise ValueError("geometry_tolerance must be positive")
        if self.probe_anchor_atom is not None:
            self.probe_anchor_atom = _strict_int(
                self.probe_anchor_atom, "probe_anchor_atom"
            )
            if self.probe_anchor_atom < 0:
                raise ValueError("probe_anchor_atom must be non-negative")
        if self.probe_anchor_atom is not None and self.probe_template_xyz is not None:
            raise ValueError(
                "Use either probe_anchor_atom or probe_template_xyz, not both"
            )
        if self.indices is not None:
            indices = tuple(_strict_int(index, "indices") for index in self.indices)
            if not indices:
                raise ValueError("indices must contain at least one grid point")
            if len(set(indices)) != len(indices):
                raise ValueError("indices must not contain duplicates")
            if any(index < 0 for index in indices):
                raise ValueError("indices must be non-negative")
            self.indices = indices
        if self.progress_callback is not None and not callable(self.progress_callback):
            raise TypeError("progress_callback must be callable or None")
        try:
            sample_combined_path = Path(self.combined_pattern.format(index=0))
        except (IndexError, KeyError, ValueError) as exc:
            raise ValueError(
                "combined_pattern must be format-compatible and contain {index}"
            ) from exc
        if "{index" not in self.combined_pattern:
            raise ValueError("combined_pattern must contain an {index} field")
        if sample_combined_path.is_absolute() or ".." in sample_combined_path.parts:
            raise ValueError("combined_pattern must stay inside probe_directory")

    def _output_paths(self) -> tuple[Path, Path]:
        """Return the resolved energy and extended-XYZ output paths."""

        molecule_path = Path(self.molecule_xyz)
        base_directory = molecule_path.resolve().parent
        stem = molecule_path.stem
        energy_path = Path(
            self.energy_output or base_directory / f"energy_values_pyscf_{stem}.tsv"
        )
        xyz_path = Path(self.xyz_output or base_directory / f"{stem}_pyscf.xyz")
        return energy_path, xyz_path

    def _selected_indices(self, grid_count: int) -> tuple[int, ...]:
        if self.indices is None:
            return tuple(range(grid_count))
        selected = tuple(self.indices)
        if selected and max(selected) >= grid_count:
            raise ValueError(
                f"indices contains {max(selected)}, but the filtered grid has "
                f"only {grid_count} points (valid range 0..{grid_count - 1})"
            )
        return selected

    def run(self, *, restart: bool = False) -> GridRunSummary:
        """Run one grid job while holding an exclusive output lock."""

        energy_path, _ = self._output_paths()
        lock_path = energy_path.with_name(f"{energy_path.name}.lock")
        with _exclusive_output_lock(lock_path):
            return self._run_unlocked(restart=restart)

    def _run_unlocked(self, *, restart: bool) -> GridRunSummary:
        molecule_path = Path(self.molecule_xyz)
        molecule_atoms = read_xyz(molecule_path)
        base_directory = molecule_path.resolve().parent
        stem = molecule_path.stem
        probe_directory = Path(
            self.probe_directory or base_directory / f"{stem}_probe"
        )
        grid_path = Path(self.grid_xyz or base_directory / f"{stem}_filtered.xyz")
        energy_path, xyz_path = self._output_paths()
        metadata_path = energy_path.with_suffix(".json")
        probe_template_path = (
            None if self.probe_template_xyz is None else Path(self.probe_template_xyz)
        )
        energy_path.parent.mkdir(parents=True, exist_ok=True)
        xyz_path.parent.mkdir(parents=True, exist_ok=True)

        grid_rows = _read_xyz_rows(grid_path)
        selected_indices = self._selected_indices(len(grid_rows))
        settings = self._checkpoint_settings(
            molecule_path,
            probe_directory,
            grid_path,
            probe_template_path,
            len(grid_rows),
        )
        fingerprint = hashlib.sha256(
            json.dumps(settings, sort_keys=True).encode("utf-8")
        ).hexdigest()
        settings["fingerprint"] = fingerprint

        probe_template_atoms = (
            None if probe_template_path is None else read_xyz(probe_template_path)
        )
        if probe_template_path is not None and _sha256_file(
            probe_template_path
        ) != settings["input_sha256"]["probe_template_xyz"]:
            raise EDAError(
                f"Probe template changed during validation: {probe_template_path}"
            )
        combined_hashes = settings["input_sha256"]["combined_geometries"]
        for index in selected_indices:
            combined_path = self._combined_path(probe_directory, index)
            combined_atoms = read_xyz(combined_path)
            if _sha256_file(combined_path) != combined_hashes[str(index)]:
                raise EDAError(
                    f"Combined geometry changed during validation: {combined_path}"
                )
            _validate_grid_geometry(
                molecule_atoms,
                combined_atoms,
                grid_rows[index],
                combined_path,
                self.geometry_tolerance,
                self.probe_anchor_atom,
                probe_template_atoms,
            )

        completed: dict[int, dict[str, str]] = {}
        if restart and energy_path.exists():
            if not metadata_path.exists():
                raise EDAError(
                    f"Cannot safely restart without metadata file: {metadata_path}"
                )
            previous = json.loads(metadata_path.read_text(encoding="utf-8"))
            if previous.get("fingerprint") != fingerprint:
                raise EDAError("Restart settings differ from the existing checkpoint")
            completed = _read_checkpoint(energy_path, fingerprint)
        else:
            # Commit the self-identifying empty checkpoint first.  If the
            # metadata update then fails, a later restart sees an old manifest
            # and refuses rather than attaching old rows to new inputs.
            _write_checkpoint_atomic(energy_path, completed, fingerprint)
            _atomic_write_text(
                metadata_path,
                json.dumps(settings, indent=2, sort_keys=True) + "\n",
            )

        super_charge = (
            self.molecule_charge + self.probe_charge
            if self.supermolecule_charge is None
            else self.supermolecule_charge
        )
        super_spin = (
            self.molecule_spin + self.probe_spin
            if self.supermolecule_spin is None
            else self.supermolecule_spin
        )

        for index in selected_indices:
            if index in completed:
                continue
            combined_path = self._combined_path(probe_directory, index)
            if _sha256_file(combined_path) != combined_hashes[str(index)]:
                raise EDAError(f"Combined geometry changed during run: {combined_path}")
            combined_atoms = read_xyz(combined_path)
            _validate_grid_geometry(
                molecule_atoms,
                combined_atoms,
                grid_rows[index],
                combined_path,
                self.geometry_tolerance,
                self.probe_anchor_atom,
                probe_template_atoms,
            )
            fragments = (
                FragmentSpec(
                    tuple(range(len(molecule_atoms))),
                    self.molecule_charge,
                    self.molecule_spin,
                    "molecule",
                ),
                FragmentSpec(
                    tuple(range(len(molecule_atoms), len(combined_atoms))),
                    self.probe_charge,
                    self.probe_spin,
                    "probe",
                ),
            )
            result = PySCFEDA(
                combined_atoms,
                fragments,
                self.config,
                charge=super_charge,
                spin=super_spin,
            ).run()
            row = result.grid_row(index, probe_fragment="probe")
            completed[index] = {
                key: _format_table_value(row[key]) for key in GRID_EDA_COLUMNS
            }
            _write_checkpoint_atomic(energy_path, completed, fingerprint)
            if self.progress_callback is not None:
                self.progress_callback(index, len(selected_indices))

        completed = _read_checkpoint(energy_path, fingerprint)
        expected_indices = set(selected_indices)
        actual_indices = set(completed)
        if actual_indices != expected_indices:
            missing = sorted(expected_indices - actual_indices)
            unexpected = sorted(actual_indices - expected_indices)
            raise EDAError(
                "Energy checkpoint grid indices do not match the filtered grid: "
                f"missing={missing[:10]}, unexpected={unexpected[:10]}"
            )
        _write_extended_xyz(xyz_path, grid_rows, completed)
        return GridRunSummary(
            energy_path,
            xyz_path,
            metadata_path,
            len(completed),
        )

    def _combined_path(self, probe_directory: Path, index: int) -> Path:
        return probe_directory / self.combined_pattern.format(index=index)

    def _checkpoint_settings(
        self,
        molecule_path: Path,
        probe_directory: Path,
        grid_path: Path,
        probe_template_path: Path | None,
        grid_count: int,
    ) -> dict[str, Any]:
        combined_hashes = {
            str(index): _sha256_file(
                self._combined_path(probe_directory, index)
            )
            for index in range(grid_count)
        }
        return {
            "schema": "pyscf-dm-eda-grid-v1",
            "pyscf_version": pyscf_version,
            "implementation_sha256": _sha256_file(Path(__file__)),
            "molecule_xyz": str(molecule_path.resolve()),
            "probe_directory": str(probe_directory.resolve()),
            "combined_pattern": self.combined_pattern,
            "grid_xyz": str(grid_path.resolve()),
            "molecule_charge": self.molecule_charge,
            "molecule_spin": self.molecule_spin,
            "probe_charge": self.probe_charge,
            "probe_spin": self.probe_spin,
            "supermolecule_charge": self.supermolecule_charge,
            "supermolecule_spin": self.supermolecule_spin,
            "geometry_tolerance": self.geometry_tolerance,
            "probe_anchor_atom": self.probe_anchor_atom,
            "probe_template_xyz": (
                None
                if probe_template_path is None
                else str(probe_template_path.resolve())
            ),
            "indices": (
                None if self.indices is None else list(self.indices)
            ),
            "scf_config": _jsonable(asdict(self.config)),
            "input_sha256": {
                "molecule_xyz": _sha256_file(molecule_path),
                "grid_xyz": _sha256_file(grid_path),
                "probe_template_xyz": (
                    None
                    if probe_template_path is None
                    else _sha256_file(probe_template_path)
                ),
                "combined_geometries": combined_hashes,
            },
        }


def read_xyz(path: str | os.PathLike[str]) -> tuple[Atom, ...]:
    """Read a conventional XYZ file and validate its atom count."""

    input_path = Path(path)
    try:
        lines = input_path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        lines = input_path.read_text(encoding="latin-1").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Not an XYZ file: {input_path}")
    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid XYZ atom count in {input_path}") from exc
    atom_lines = [line for line in lines[2:] if line.strip()]
    if len(atom_lines) != atom_count:
        raise ValueError(
            f"XYZ atom count mismatch in {input_path}: expected {atom_count}, "
            f"found {len(atom_lines)}"
        )
    atoms: list[Atom] = []
    for line_number, line in enumerate(atom_lines, start=3):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"Invalid XYZ row {line_number} in {input_path}")
        try:
            xyz = tuple(float(value) for value in fields[1:4])
        except ValueError as exc:
            raise ValueError(
                f"Invalid coordinate on XYZ row {line_number} in {input_path}"
            ) from exc
        atoms.append(Atom(fields[0], xyz))
    return tuple(atoms)


def _coerce_atom(atom: Atom | Sequence[Any]) -> Atom:
    if isinstance(atom, Atom):
        return atom
    if len(atom) == 2:
        symbol, xyz = atom
    elif len(atom) == 4:
        symbol, *xyz = atom
    else:
        raise ValueError(f"Atom must be (symbol, xyz) or (symbol, x, y, z): {atom!r}")
    return Atom(str(symbol), tuple(float(value) for value in xyz))


def _strict_int(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, np.integer)
    ):
        raise TypeError(f"{name} must be an integer")
    return int(value)


def _normalize_fragments(
    fragments: Sequence[FragmentSpec], atom_count: int
) -> tuple[FragmentSpec, ...]:
    if len(fragments) < 2:
        raise ValueError("EDA requires at least two fragments")
    normalized: list[FragmentSpec] = []
    occupied: set[int] = set()
    labels: set[str] = set()
    for number, fragment in enumerate(fragments, start=1):
        if not isinstance(fragment, FragmentSpec):
            raise TypeError("fragments must contain FragmentSpec instances")
        label = fragment.label or f"fragment_{number}"
        if label in labels:
            raise ValueError(f"Duplicate fragment label: {label!r}")
        labels.add(label)
        indices = set(fragment.atom_indices)
        if min(indices) < 0 or max(indices) >= atom_count:
            raise ValueError(f"Fragment {label!r} has an out-of-range atom index")
        overlap = occupied.intersection(indices)
        if overlap:
            raise ValueError(f"Atoms assigned to multiple fragments: {sorted(overlap)}")
        occupied.update(indices)
        normalized.append(
            FragmentSpec(fragment.atom_indices, fragment.charge, fragment.spin, label)
        )
    missing = sorted(set(range(atom_count)) - occupied)
    if missing:
        raise ValueError(f"Atoms not assigned to a fragment: {missing}")
    return tuple(normalized)


def _build_molecule(
    atoms: Sequence[Atom],
    config: SCFConfig,
    *,
    charge: int,
    spin: int,
    real_indices: Iterable[int],
) -> Any:
    real = set(real_indices)
    atom_spec = [
        (atom.symbol if index in real else f"ghost-{atom.symbol}", atom.xyz)
        for index, atom in enumerate(atoms)
    ]
    kwargs: dict[str, Any] = {
        "atom": atom_spec,
        "unit": config.unit,
        "basis": config.basis,
        "charge": charge,
        "spin": spin,
        "symmetry": False,
        "verbose": config.verbose,
        "max_memory": config.max_memory,
    }
    if config.ecp is not None:
        kwargs["ecp"] = config.ecp
    try:
        return gto.M(**kwargs)
    except Exception as exc:
        raise IncompatibleFragmentError(
            f"Could not build molecule with charge={charge}, spin={spin}: {exc}"
        ) from exc


def _validate_ao_overlap(mol: Any, config: SCFConfig, label: str) -> None:
    overlap = mol.intor_symmetric("int1e_ovlp")
    eigenvalues = np.linalg.eigvalsh(overlap)
    largest = float(np.max(eigenvalues))
    smallest = float(np.min(eigenvalues))
    if smallest <= config.linear_dep_threshold * largest:
        raise IncompatibleFragmentError(
            f"AO basis is linearly dependent for {label!r}: "
            f"min/max overlap eigenvalue={smallest:.3e}/{largest:.3e}"
        )


def _run_scf(
    mol: Any,
    config: SCFConfig,
    label: str,
    *,
    dm0: np.ndarray | None = None,
    physical_mol: Any | None = None,
    force_unrestricted: bool = False,
) -> SCFState:
    """Run one SCF in the full ghost basis.

    ``physical_mol`` holds the real atoms of this system only.  It is used for
    the atom-pair dispersion correction and for the solvent cavity, both of
    which must not see ghost atoms.
    """

    method = config.method.strip()
    method_lower = method.lower()
    if method_lower == "uhf":
        unrestricted = True
    elif method_lower == "rhf":
        if force_unrestricted or mol.spin != 0:
            raise IncompatibleFragmentError(
                f"RHF was requested for open-shell state {label!r}"
            )
        unrestricted = False
    else:
        unrestricted = force_unrestricted or mol.spin != 0
    if method.lower() in _HF_METHODS:
        mf = scf.UHF(mol) if unrestricted else scf.RHF(mol)
    else:
        mf = dft.UKS(mol) if unrestricted else dft.RKS(mol)
        mf.xc = method
        mf.grids.level = config.grid_level

    if config.density_fit:
        mf = mf.density_fit(auxbasis=config.auxbasis)
    if config.point_charges is not None:
        from pyscf import qmmm

        coords, charges = config.point_charges
        mf = qmmm.mm_charge(mf, coords, charges, unit="Angstrom")
    if config.solvent is not None:
        # Attach the solvent after density fitting so the DF J/K builder stays
        # the inner class (the order PySCF documents for ``mol.RKS().density_fit().PCM()``).
        solvent = _solvent_model(mol, physical_mol, config)
        mf = mf.SMD(solvent) if config.solvent == "SMD" else mf.PCM(solvent)
    mf.conv_tol = config.conv_tol
    mf.max_cycle = config.max_cycle
    mf.max_memory = config.max_memory
    mf.level_shift = config.level_shift
    mf.damp = config.damp
    mf.init_guess = config.init_guess
    # Keep D3/D4 out of the SCF energy here.  Fragment SCFs contain ghost
    # atoms, while atom-pair dispersion must be evaluated on real atoms only.
    # ``None`` is PySCF's documented "disabled" value; relying on False being
    # falsy would couple this code to an implementation detail.
    mf.disp = None

    try:
        mf.kernel(dm0=dm0)
        if not mf.converged and config.newton_fallback:
            retry_dm = mf.make_rdm1()
            mf = mf.newton()
            mf.max_cycle = config.max_cycle
            mf.conv_tol = config.conv_tol
            mf.kernel(dm0=retry_dm)
    except Exception as exc:
        raise EDAError(f"SCF failed for {label!r}: {exc}") from exc
    if not mf.converged:
        raise SCFConvergenceError(
            f"SCF did not converge for {label!r} in {config.max_cycle} cycles"
        )

    spin_dm, occupied = _spin_density_and_occupied(mf, label)
    dispersion_energy = _dispersion_energy(
        physical_mol if physical_mol is not None else mol,
        config,
    )
    solvent_energy, cds_energy, summary_energy = _scf_solvent_energies(mf, label)
    return SCFState(
        label=label,
        mol=mol,
        mf=mf,
        spin_dm=spin_dm,
        occupied_coeff=occupied,
        total_energy=float(np.real_if_close(mf.e_tot)) + dispersion_energy,
        dispersion_energy=dispersion_energy,
        solvent_energy=solvent_energy,
        solvent_cds_energy=cds_energy,
        solvent_scf_summary_energy=summary_energy,
    )


def _dispersion_energy(mol: Any, config: SCFConfig) -> float:
    """Evaluate D3/D4 on a molecule containing real atoms only."""

    if not config.dispersion:
        return 0.0
    method = config.method.strip()
    if method.lower() in _HF_METHODS:
        evaluator = scf.UHF(mol) if mol.spin != 0 else scf.RHF(mol)
    else:
        evaluator = dft.UKS(mol) if mol.spin != 0 else dft.RKS(mol)
        evaluator.xc = method
    evaluator.verbose = config.verbose
    try:
        return float(
            np.real_if_close(evaluator.get_dispersion(disp=config.dispersion))
        )
    except Exception as exc:
        raise EDAError(
            f"Dispersion calculation {config.dispersion!r} failed: {exc}"
        ) from exc


class _CavityMixin:
    """Build the solvent cavity from ``cavity_mol`` but integrate over ``mol``.

    PySCF assigns every atom of ``mol`` a cavity sphere, and a ghost atom gets
    the table's Z=0 radius (2.0 Å), so a counterpoise fragment would otherwise
    be solvated in a spurious supermolecule-sized cavity.  While the surface,
    the nuclear potential on the surface, and (for SMD) the CDS term are built
    the object temporarily points at the real-atom molecule; the AO integrals
    always use the full ghost basis.
    """

    _keys = {"cavity_mol"}
    cavity_mol: Any = None

    @contextmanager
    def _cavity_scope(self):
        if self.cavity_mol is None:
            yield
            return
        mol = self.mol
        self.mol = self.cavity_mol
        try:
            yield
        finally:
            self.mol = mol

    def build(self, ng=None):
        with self._cavity_scope():
            return super().build(ng)


class _CavityPCM(_CavityMixin, pyscf_pcm.PCM):
    def __init__(self, mol: Any, cavity_mol: Any | None = None) -> None:
        super().__init__(mol)
        self.cavity_mol = cavity_mol


class _CavitySMD(_CavityMixin, pyscf_smd.SMD):
    def __init__(self, mol: Any, cavity_mol: Any | None = None) -> None:
        super().__init__(mol)
        self.cavity_mol = cavity_mol

    def get_cds(self):
        with self._cavity_scope():
            return super().get_cds()


def _solvent_option_keys(model: str) -> frozenset[str]:
    """Attributes of the PySCF solvent object that ``solvent_options`` may set."""

    cls = _CavitySMD if model == "SMD" else _CavityPCM
    keys: set[str] = set()
    for klass in cls.__mro__:
        keys.update(getattr(klass, "_keys", ()))
    return frozenset(keys - _RESERVED_SOLVENT_OPTIONS)


def _validate_solvent_config(config: SCFConfig) -> None:
    """Normalise and validate the solvent fields of a frozen ``SCFConfig``."""

    model = config.solvent
    if model is None:
        for name in ("solvent_eps", "solvent_name", "solvent_options"):
            if getattr(config, name) is not None:
                raise ValueError(f"{name} requires a solvent model (solvent=...)")
        return
    if not isinstance(model, str) or not model.strip():
        raise ValueError("solvent must be a non-empty model name or None")
    key = model.strip().lower()
    if key in _UNSUPPORTED_SOLVENT_MODELS:
        raise ValueError(
            f"Solvent model {model!r} is not supported: its cavity is built "
            "per atom of the ghost-padded fragment molecule and cannot be "
            "restricted to the real atoms. Use one of "
            f"{sorted(set(_SOLVENT_MODEL_ALIASES.values()))}."
        )
    try:
        canonical = _SOLVENT_MODEL_ALIASES[key]
    except KeyError as exc:
        raise ValueError(
            f"Unknown solvent model {model!r}; choose from "
            f"{sorted(_SOLVENT_MODEL_ALIASES)}"
        ) from exc
    object.__setattr__(config, "solvent", canonical)

    eps = config.solvent_eps
    if eps is not None:
        if isinstance(eps, (bool, np.bool_)) or not isinstance(
            eps, (int, float, np.integer, np.floating)
        ):
            raise TypeError("solvent_eps must be a number or None")
        eps = float(eps)
        if not eps >= 1.0:  # also rejects NaN; inf is the conductor limit
            raise ValueError("solvent_eps must be >= 1 (1 corresponds to vacuum)")
        object.__setattr__(config, "solvent_eps", eps)

    name = config.solvent_name
    if name is not None:
        if not isinstance(name, str) or not name.strip():
            raise ValueError("solvent_name must be a non-empty string or None")
        name = name.strip().lower()
        if name not in pyscf_smd.solvent_db:
            raise ValueError(
                f"Unknown solvent name {config.solvent_name!r}; it is not in "
                "PySCF's SMD database (pyscf.solvent.smd.solvent_db)"
            )
        object.__setattr__(config, "solvent_name", name)
    if canonical == "SMD" and name is None and eps is not None:
        raise ValueError(
            "SMD needs solvent_name for its CDS descriptors; solvent_eps can "
            "only override the dielectric constant of a named solvent"
        )

    options = config.solvent_options
    if options is not None:
        if not isinstance(options, Mapping):
            raise TypeError(
                "solvent_options must be a mapping of PySCF solvent attributes"
            )
        allowed = _solvent_option_keys(canonical)
        normalized: dict[str, Any] = {}
        for option, value in options.items():
            if not isinstance(option, str):
                raise TypeError("solvent_options keys must be strings")
            if option in _RESERVED_SOLVENT_OPTIONS or option not in allowed:
                raise ValueError(
                    f"solvent option {option!r} is not settable for "
                    f"{canonical}; allowed: {sorted(allowed)}"
                )
            normalized[option] = value
        object.__setattr__(config, "solvent_options", normalized)


def _solvent_model(mol: Any, physical_mol: Any | None, config: SCFConfig) -> Any:
    """Create the PySCF solvent object for ``mol`` with a real-atom cavity."""

    cavity_mol = None if physical_mol is None or physical_mol is mol else physical_mol
    if config.solvent == "SMD":
        solvent = _CavitySMD(mol, cavity_mol=cavity_mol)
        solvent.solvent = config.solvent_name or "water"
        if config.solvent_eps is not None:
            solvent.eps = config.solvent_eps
    else:
        solvent = _CavityPCM(mol, cavity_mol=cavity_mol)
        solvent.method = config.solvent
        if config.solvent_eps is not None:
            solvent.eps = config.solvent_eps
        elif config.solvent_name is not None:
            solvent.eps = float(pyscf_smd.solvent_db[config.solvent_name][5])
    solvent.verbose = config.verbose
    solvent.max_memory = config.max_memory
    for option, value in (config.solvent_options or {}).items():
        setattr(solvent, option, value)
    return solvent


def _has_solvent(mf: Any) -> bool:
    return getattr(mf, "with_solvent", None) is not None


def _scf_solvent_energies(mf: Any, label: str) -> tuple[float, float, float | None]:
    """Reaction-field energy of the converged SCF density.

    Returns ``(total, cds, scf_summary)`` in Hartree.  ``total`` is evaluated
    from the final density so that it is consistent with ``mf.e_tot`` and the
    density used in the decomposition; PySCF's ``scf_summary`` value is kept
    only as a diagnostic.
    """

    if not _has_solvent(mf):
        return 0.0, 0.0, None
    solvent = mf.with_solvent
    dm = np.asarray(mf.make_rdm1())
    if dm.ndim == 3:
        dm = dm[0] + dm[1]
    try:
        electrostatic = float(np.real_if_close(solvent.kernel(dm)[0]))
        cds = (
            float(np.real_if_close(solvent.get_cds()))
            if isinstance(solvent, pyscf_smd.SMD)
            else 0.0
        )
    except Exception as exc:
        raise EDAError(f"Solvent energy evaluation failed for {label!r}: {exc}") from exc
    summary = getattr(mf, "scf_summary", {}).get("e_solvent")
    summary_total: float | None = None
    if summary is not None:
        summary_total = float(np.real_if_close(summary)) + float(
            np.real_if_close(mf.scf_summary.get("e_cds", 0.0))
        )
    return electrostatic + cds, cds, summary_total


def _solvent_energy_for_density(state: SCFState, total_dm: np.ndarray) -> float:
    """Reaction-field energy of an arbitrary spin-summed density in ``state``'s cavity.

    The CDS term (SMD) depends only on the cavity geometry and is therefore
    the same for every density evaluated in this cavity.
    """

    energy = state.mf.with_solvent.kernel(np.asarray(total_dm))[0]
    return float(np.real_if_close(energy)) + state.solvent_cds_energy


def _solvent_diagnostics(
    super_state: SCFState,
    fragment_states: Sequence[SCFState],
    p0: np.ndarray,
    pauli_dm: np.ndarray,
    desolvation: float,
) -> dict[str, Any] | None:
    """Split the desolvation term along the DM-EDA density sequence.

    ``E_solv`` of the promolecule and Pauli densities is evaluated in the
    supermolecule cavity, which turns the single desolvation term into a
    frozen (cavity + mutual screening) part, a Pauli response and an orbital
    relaxation response.  The three parts sum to ``desolvation`` exactly.
    """

    if not _has_solvent(super_state.mf):
        return None
    promolecule = _solvent_energy_for_density(super_state, np.sum(p0, axis=0))
    pauli = _solvent_energy_for_density(super_state, np.sum(pauli_dm, axis=0))
    fragment_total = sum(state.solvent_energy for state in fragment_states)
    summary_deviation = max(
        abs(state.solvent_energy - state.solvent_scf_summary_energy)
        for state in (super_state, *fragment_states)
        if state.solvent_scf_summary_energy is not None
    ) if any(
        state.solvent_scf_summary_energy is not None
        for state in (super_state, *fragment_states)
    ) else None
    surface_points = {
        state.label: int(state.mf.with_solvent.surface["grid_coords"].shape[0])
        for state in (super_state, *fragment_states)
    }
    return {
        "supermolecule_solvent_energy_hartree": super_state.solvent_energy,
        "fragment_solvent_energy_hartree": {
            state.label: state.solvent_energy for state in fragment_states
        },
        "supermolecule_cds_energy_hartree": super_state.solvent_cds_energy,
        "fragment_cds_energy_hartree": {
            state.label: state.solvent_cds_energy for state in fragment_states
        },
        "promolecule_in_supermolecule_cavity_hartree": promolecule,
        "pauli_in_supermolecule_cavity_hartree": pauli,
        "desolvation_hartree": desolvation,
        "desolvation_frozen_hartree": promolecule - fragment_total,
        "desolvation_pauli_response_hartree": pauli - promolecule,
        "desolvation_polarization_response_hartree": (
            super_state.solvent_energy - pauli
        ),
        "scf_summary_max_abs_deviation_hartree": summary_deviation,
        "cavity_surface_points": surface_points,
    }


def _solvent_metadata(config: SCFConfig, super_state: SCFState) -> dict[str, Any] | None:
    if config.solvent is None:
        return None
    solvent = super_state.mf.with_solvent
    if config.solvent == "SMD":
        descriptors = solvent.solvent_descriptors or pyscf_smd.solvent_db[solvent.solvent]
        eps = float(solvent.eps or descriptors[5])
        name = solvent.solvent
    else:
        eps = float(solvent.eps)
        name = config.solvent_name
    return {
        "model": config.solvent,
        "implementation": (
            "pyscf.solvent.smd" if config.solvent == "SMD" else "pyscf.solvent.pcm"
        ),
        "eps": eps,
        "solvent_name": name,
        "options": _jsonable(dict(config.solvent_options or {})),
        "cavity_convention": (
            "fragments: real atoms only (ghost atoms carry no sphere); "
            "supermolecule: all atoms"
        ),
        "desolvation_definition": (
            "E_solv[supermolecule] - sum_A E_solv[fragment A], each evaluated "
            "with its own converged density in its own cavity"
        ),
    }


def _spin_density_and_occupied(
    mf: Any, label: str
) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray]]:
    coeff = np.asarray(mf.mo_coeff)
    occupation = np.asarray(mf.mo_occ)
    density = np.asarray(mf.make_rdm1())
    tolerance = 1e-7

    if coeff.ndim == 2:
        if np.any(
            ~np.logical_or(
                np.isclose(occupation, 0.0, atol=tolerance),
                np.isclose(occupation, 2.0, atol=tolerance),
            )
        ):
            raise IncompatibleFragmentError(
                f"Fractional restricted occupations are not supported for {label!r}"
            )
        occupied = coeff[:, occupation > tolerance]
        spin_dm = np.stack((density * 0.5, density * 0.5))
        return spin_dm, (occupied, occupied.copy())

    if coeff.ndim != 3 or coeff.shape[0] != 2:
        raise IncompatibleFragmentError(
            f"Unsupported orbital coefficient shape for {label!r}: {coeff.shape}"
        )
    if np.any(
        ~np.logical_or(
            np.isclose(occupation, 0.0, atol=tolerance),
            np.isclose(occupation, 1.0, atol=tolerance),
        )
    ):
        raise IncompatibleFragmentError(
            f"Fractional unrestricted occupations are not supported for {label!r}"
        )
    occupied = tuple(
        coeff[spin][:, occupation[spin] > tolerance] for spin in (0, 1)
    )
    return np.asarray(density), occupied  # type: ignore[return-value]


def _assert_same_ao_space(
    super_mol: Any,
    fragment_mol: Any,
    super_overlap: np.ndarray,
    super_kinetic: np.ndarray,
    label: str,
) -> None:
    if super_mol.nao_nr() != fragment_mol.nao_nr():
        raise IncompatibleFragmentError(
            f"Fragment {label!r} has {fragment_mol.nao_nr()} AOs; "
            f"supermolecule has {super_mol.nao_nr()}"
        )
    fragment_overlap = fragment_mol.intor_symmetric("int1e_ovlp")
    fragment_kinetic = fragment_mol.intor_symmetric("int1e_kin")
    if not np.allclose(super_overlap, fragment_overlap, atol=1e-10, rtol=1e-10):
        raise IncompatibleFragmentError(f"AO overlap mismatch for fragment {label!r}")
    if not np.allclose(super_kinetic, fragment_kinetic, atol=1e-10, rtol=1e-10):
        raise IncompatibleFragmentError(f"AO kinetic mismatch for fragment {label!r}")


def _build_pauli_density(
    states: Sequence[SCFState],
    overlap: np.ndarray,
    threshold: float,
) -> tuple[np.ndarray, dict[str, float]]:
    pauli_spin: list[np.ndarray] = []
    minimum_eigenvalue = np.inf
    maximum_condition = 1.0
    for spin in (0, 1):
        occupied = np.concatenate(
            [state.occupied_coeff[spin] for state in states], axis=1
        )
        if occupied.shape[1] == 0:
            pauli_spin.append(np.zeros_like(overlap))
            continue
        gram = occupied.conj().T @ overlap @ occupied
        gram = 0.5 * (gram + gram.conj().T)
        eigenvalues, eigenvectors = np.linalg.eigh(gram)
        largest = float(np.max(eigenvalues))
        smallest = float(np.min(eigenvalues))
        minimum_eigenvalue = min(minimum_eigenvalue, smallest)
        if smallest <= threshold * largest:
            raise IncompatibleFragmentError(
                "Occupied fragment orbitals are linearly dependent in the "
                f"supermolecule AO metric for spin channel {spin}: "
                f"min/max={smallest:.3e}/{largest:.3e}"
            )
        maximum_condition = max(maximum_condition, largest / smallest)
        inverse = (eigenvectors / eigenvalues) @ eigenvectors.conj().T
        density = occupied @ inverse @ occupied.conj().T
        pauli_spin.append(0.5 * (density + density.conj().T))
    if not np.isfinite(minimum_eigenvalue):
        minimum_eigenvalue = 0.0
    return np.asarray(pauli_spin), {
        "pauli_occupied_overlap_min_eigenvalue": minimum_eigenvalue,
        "pauli_occupied_overlap_condition": maximum_condition,
    }


def _batch_jk(
    super_state: SCFState, spin_densities: Sequence[np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    dm = np.asarray(spin_densities)
    vj, vk = super_state.mf.get_jk(super_state.mol, dm, hermi=1)
    return np.asarray(vj), np.asarray(vk)


def _e0(
    spin_dm: np.ndarray,
    hcore: np.ndarray,
    spin_j: np.ndarray,
    spin_k: np.ndarray,
) -> float:
    total_dm = np.sum(spin_dm, axis=0)
    total_j = np.sum(spin_j, axis=0)
    return (
        _trace_product(total_dm, hcore)
        + 0.5 * _trace_product(total_dm, total_j)
        - 0.5
        * sum(_trace_product(spin_dm[spin], spin_k[spin]) for spin in (0, 1))
    )


def _mulliken_fragment_charge_transfer(
    spin_dm: np.ndarray,
    mol: Any,
    overlap: np.ndarray,
    fragments: Sequence[FragmentSpec],
) -> dict[str, float]:
    total_dm = np.sum(spin_dm, axis=0)
    ao_population = np.einsum("ij,ji->i", total_dm, overlap).real
    ao_slices = mol.aoslice_by_atom()
    charges = np.asarray(mol.atom_charges(), dtype=float)
    result: dict[str, float] = {}
    for fragment in fragments:
        electron_population = sum(
            float(np.sum(ao_population[ao_slices[index, 2] : ao_slices[index, 3]]))
            for index in fragment.atom_indices
        )
        fragment_charge = (
            float(np.sum(charges[list(fragment.atom_indices)])) - electron_population
        )
        result[fragment.label] = fragment_charge - fragment.charge
    return result


def _electron_count(spin_dm: np.ndarray, overlap: np.ndarray) -> dict[str, float]:
    alpha = _trace_product(spin_dm[0], overlap)
    beta = _trace_product(spin_dm[1], overlap)
    return {"alpha": alpha, "beta": beta, "total": alpha + beta}


def _idempotency_error(spin_dm: np.ndarray, overlap: np.ndarray) -> float:
    return max(
        _max_abs(spin_dm[spin] @ overlap @ spin_dm[spin] - spin_dm[spin])
        for spin in (0, 1)
    )


def _hermiticity_error(spin_dm: np.ndarray) -> float:
    return max(
        _max_abs(spin_dm[spin] - spin_dm[spin].conj().T) for spin in (0, 1)
    )


def _trace_product(left: np.ndarray, right: np.ndarray) -> float:
    value = np.einsum("ij,ji->", left, right)
    return float(np.real_if_close(value))


def _max_abs(array: np.ndarray) -> float:
    return float(np.max(np.abs(array))) if array.size else 0.0


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return repr(value)


def _read_xyz_rows(path: Path) -> list[list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Not an XYZ file: {path}")
    count = int(lines[0].strip())
    rows = [line.split() for line in lines[2:] if line.strip()]
    if len(rows) != count:
        raise ValueError(
            f"XYZ row count mismatch in {path}: expected {count}, found {len(rows)}"
        )
    return rows


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_grid_geometry(
    molecule_atoms: Sequence[Atom],
    combined_atoms: Sequence[Atom],
    grid_row: Sequence[str],
    combined_path: Path,
    tolerance: float,
    probe_anchor_atom: int | None,
    probe_template_atoms: Sequence[Atom] | None,
) -> None:
    """Validate that a combined geometry is an unmodified host + placed probe.

    For a linear polyatomic probe the Kabsch rotation has an unphysical
    free rotation around the molecular axis.  If the template origin is not
    exactly on that axis, the recovered translation becomes ambiguous and can
    fail the grid-point check for a perfectly valid placement.  Use
    ``probe_anchor_atom`` for linear probes with an off-axis template origin.
    """

    if len(combined_atoms) <= len(molecule_atoms):
        raise EDAError(f"No probe atoms found in {combined_path}")
    combined_molecule = combined_atoms[: len(molecule_atoms)]
    for index, (expected, actual) in enumerate(
        zip(molecule_atoms, combined_molecule)
    ):
        if expected.symbol.lower() != actual.symbol.lower() or not np.allclose(
            expected.xyz, actual.xyz, atol=tolerance, rtol=0.0
        ):
            raise EDAError(
                f"Host atom {index} in {combined_path} does not match molecule_xyz"
            )
    if len(grid_row) < 4:
        raise EDAError("A filtered grid row must contain atom label and X/Y/Z")
    try:
        grid_coordinate = np.asarray(grid_row[1:4], dtype=float)
    except ValueError as exc:
        raise EDAError(f"Invalid filtered-grid coordinate: {grid_row[:4]}") from exc
    probe_coordinates = np.asarray(
        [atom.xyz for atom in combined_atoms[len(molecule_atoms) :]], dtype=float
    )
    if probe_anchor_atom is not None:
        if probe_anchor_atom >= len(probe_coordinates):
            raise EDAError(
                f"Probe anchor atom {probe_anchor_atom} is out of range in "
                f"{combined_path}"
            )
        anchor = probe_coordinates[probe_anchor_atom]
        if not np.allclose(anchor, grid_coordinate, atol=tolerance, rtol=0.0):
            raise EDAError(
                f"Probe anchor atom in {combined_path} does not match its "
                "filtered grid point"
            )
        return

    if len(probe_coordinates) == 1:
        if not np.allclose(
            probe_coordinates[0], grid_coordinate, atol=tolerance, rtol=0.0
        ):
            raise EDAError(
                f"Monatomic probe in {combined_path} does not match its "
                "filtered grid point"
            )
        return

    if probe_template_atoms is None:
        raise EDAError(
            "Polyatomic probe geometry requires probe_anchor_atom or "
            "probe_template_xyz for an unambiguous grid-point check"
        )
    probe_atoms = combined_atoms[len(molecule_atoms) :]
    if len(probe_template_atoms) != len(probe_atoms):
        raise EDAError(
            f"Probe template atom count does not match {combined_path}"
        )
    template_symbols = tuple(atom.symbol.lower() for atom in probe_template_atoms)
    probe_symbols = tuple(atom.symbol.lower() for atom in probe_atoms)
    if template_symbols != probe_symbols:
        raise EDAError(f"Probe template atom order does not match {combined_path}")

    template_coordinates = np.asarray(
        [atom.xyz for atom in probe_template_atoms], dtype=float
    )
    _, translation, rmsd = _rigid_transform(
        template_coordinates, probe_coordinates
    )
    if rmsd > tolerance:
        raise EDAError(
            f"Probe in {combined_path} is not a rigid copy of its template "
            f"(RMSD={rmsd:.3e})"
        )
    if not np.allclose(translation, grid_coordinate, atol=tolerance, rtol=0.0):
        raise EDAError(
            f"Probe-template origin in {combined_path} does not match its "
            "filtered grid point"
        )


def _rigid_transform(
    source: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return the proper row-vector rigid transform ``source @ R + t``."""

    source_centroid = np.mean(source, axis=0)
    target_centroid = np.mean(target, axis=0)
    source_centered = source - source_centroid
    target_centered = target - target_centroid
    left, _, right_t = np.linalg.svd(source_centered.T @ target_centered)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    translation = target_centroid - source_centroid @ rotation
    residual = source @ rotation + translation - target
    rmsd = float(np.sqrt(np.mean(np.sum(residual * residual, axis=1))))
    return rotation, translation, rmsd


def _read_checkpoint(
    path: Path, expected_fingerprint: str
) -> dict[int, dict[str, str]]:
    completed: dict[int, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        preamble = stream.readline().rstrip("\r\n")
        if preamble != f"{_CHECKPOINT_PREFIX}{expected_fingerprint}":
            raise EDAError(
                f"Checkpoint fingerprint does not match its metadata: {path}"
            )
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GRID_EDA_COLUMNS:
            raise EDAError(f"Unexpected checkpoint columns in {path}")
        for row in reader:
            if None in row or any(
                row.get(column) in (None, "") for column in GRID_EDA_COLUMNS
            ):
                raise EDAError(f"Incomplete checkpoint row in {path}")
            try:
                index_value = float(row["Grid_Index"])
                numeric_values = [float(row[column]) for column in GRID_EDA_COLUMNS[1:]]
            except (TypeError, ValueError) as exc:
                raise EDAError(f"Non-numeric checkpoint row in {path}") from exc
            if not index_value.is_integer() or index_value < 0:
                raise EDAError(f"Invalid grid index {row['Grid_Index']!r} in {path}")
            if not np.all(np.isfinite(numeric_values)):
                raise EDAError(f"Non-finite checkpoint value in {path}")
            index = int(index_value)
            if index in completed:
                raise EDAError(f"Duplicate grid index {index} in {path}")
            completed[index] = row
    return completed


def _write_checkpoint_atomic(
    path: Path,
    completed: Mapping[int, Mapping[str, str]],
    fingerprint: str,
) -> None:
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
            newline="",
            encoding="utf-8",
        ) as stream:
            temporary_name = stream.name
            stream.write(f"{_CHECKPOINT_PREFIX}{fingerprint}\n")
            writer = csv.DictWriter(
                stream,
                fieldnames=GRID_EDA_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for index in sorted(completed):
                writer.writerow(completed[index])
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


@contextmanager
def _exclusive_output_lock(path: Path):
    """Hold a non-blocking advisory lock for one grid output path."""

    path.parent.mkdir(parents=True, exist_ok=True)
    stream = path.open("a+b")
    if os.name == "nt":
        import msvcrt

        stream.seek(0, os.SEEK_END)
        if stream.tell() == 0:
            stream.write(b"\0")
            stream.flush()
        stream.seek(0)
        try:
            msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as exc:
            stream.close()
            raise EDAError(f"Another grid runner holds the output lock: {path}") from exc

        try:
            yield
        finally:
            stream.seek(0)
            try:
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            finally:
                stream.close()
        return

    import fcntl

    try:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        stream.close()
        raise EDAError(f"Another grid runner holds the output lock: {path}") from exc
    try:
        yield
    finally:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        finally:
            stream.close()


def _atomic_write_text(path: Path, content: str) -> None:
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
            newline="\n",
            encoding="utf-8",
        ) as stream:
            temporary_name = stream.name
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def _write_extended_xyz(
    path: Path,
    grid_rows: Sequence[Sequence[str]],
    completed: Mapping[int, Mapping[str, str]],
) -> None:
    energy_columns = GRID_EDA_COLUMNS[1:]
    written_indices = sorted(completed)
    lines = [
        str(len(written_indices)),
        "Grid with PySCF DM-EDA values (kcal/mol)",
    ]
    for index in written_indices:
        values = [completed[index][column] for column in energy_columns]
        lines.append(" ".join([*grid_rows[index], *values]))
    _atomic_write_text(path, "\n".join(lines) + "\n")


def _format_table_value(value: float | int) -> str:
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):.12f}"


def _parse_fragment_indices(text: str, atom_count: int) -> tuple[int, ...]:
    """Parse one-based CLI selections such as ``1-3,7``."""

    indices: list[int] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            if "-" in item:
                first_text, last_text = item.split("-", maxsplit=1)
                first, last = int(first_text), int(last_text)
                if last < first:
                    raise argparse.ArgumentTypeError(f"Invalid atom range: {item}")
                indices.extend(range(first - 1, last))
            else:
                indices.append(int(item) - 1)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid atom selection: {item}") from exc
    if not indices:
        raise argparse.ArgumentTypeError(f"Empty fragment selection: {text!r}")
    if min(indices) < 0 or max(indices) >= atom_count:
        raise argparse.ArgumentTypeError(
            f"Fragment selection {text!r} is outside 1..{atom_count}"
        )
    if len(set(indices)) != len(indices):
        raise argparse.ArgumentTypeError(
            f"Fragment selection {text!r} contains duplicate atoms"
        )
    return tuple(indices)


def _parse_index_subset(text: str) -> tuple[int, ...]:
    """Parse a zero-based grid shard selection such as ``0-99,150``."""

    indices: list[int] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            if "-" in item:
                first_text, last_text = item.split("-", maxsplit=1)
                first, last = int(first_text), int(last_text)
                if last < first:
                    raise argparse.ArgumentTypeError(f"Invalid index range: {item}")
                indices.extend(range(first, last + 1))
            else:
                indices.append(int(item))
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid index selection: {item}") from exc
    if not indices or min(indices) < 0:
        raise argparse.ArgumentTypeError(f"Invalid index selection: {text!r}")
    if len(set(indices)) != len(indices):
        raise argparse.ArgumentTypeError(
            f"Index selection {text!r} contains duplicate entries"
        )
    return tuple(indices)


def _add_scf_config_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--method", default="r2scan")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--ecp", default=None)
    parser.add_argument("--dispersion", default=None)
    parser.add_argument("--grid-level", type=int, default=4)
    parser.add_argument("--conv-tol", type=float, default=1e-9)
    parser.add_argument("--max-cycle", type=int, default=100)
    parser.add_argument("--max-memory", type=float, default=4000.0)
    parser.add_argument("--density-fit", action="store_true")
    parser.add_argument("--auxbasis", default=None)
    parser.add_argument("--level-shift", type=float, default=0.0)
    parser.add_argument("--damp", type=float, default=0.0)
    parser.add_argument("--init-guess", default="minao")
    parser.add_argument("--linear-dep-threshold", type=float, default=1e-9)
    reference_group = parser.add_mutually_exclusive_group()
    reference_group.add_argument(
        "--unrestricted",
        dest="unrestricted",
        action="store_const",
        const=True,
        default=None,
        help="Force UHF/UKS, including for a closed-shell supermolecule",
    )
    reference_group.add_argument(
        "--restricted",
        dest="unrestricted",
        action="store_const",
        const=False,
        help="Force RHF/RKS for the supermolecule; rejected for nonzero total spin",
    )
    parser.add_argument("--validation-tol", type=float, default=1e-7)
    parser.add_argument(
        "--strict-validation",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Raise an error carrying the full result when validation fails",
    )
    parser.add_argument("--verbose", type=int, default=0)
    solvent_group = parser.add_argument_group(
        "implicit solvent (PySCF PCM/SMD)",
        "Fragments are solvated in their own real-atom cavities and the change "
        "of the reaction-field energy is reported as Desolvation.",
    )
    solvent_group.add_argument(
        "--solvent",
        default=None,
        help="Solvent model: cpcm (alias pcm), iefpcm, ssvpe, cosmo, or smd",
    )
    solvent_group.add_argument(
        "--solvent-eps",
        type=float,
        default=None,
        help="Dielectric constant for PCM models (PySCF default 78.3553, water)",
    )
    solvent_group.add_argument(
        "--solvent-name",
        default=None,
        help=(
            "Solvent from PySCF's SMD database, e.g. water or acetonitrile; "
            "supplies the SMD descriptors or the PCM dielectric constant"
        ),
    )
    solvent_group.add_argument(
        "--solvent-option",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "Extra attribute for the PySCF solvent object, e.g. "
            "lebedev_order=29 or vdw_scale=1.1 (repeatable; values are "
            "parsed as JSON when possible)"
        ),
    )


def _parse_solvent_options(items: Sequence[str]) -> dict[str, Any] | None:
    if not items:
        return None
    options: dict[str, Any] = {}
    for item in items:
        key, separator, raw = item.partition("=")
        key = key.strip()
        if not separator or not key:
            raise ValueError(f"--solvent-option expects KEY=VALUE, got {item!r}")
        raw = raw.strip()
        try:
            value: Any = json.loads(raw)
        except ValueError:
            value = raw
        options[key] = value
    return options


def _scf_config_from_args(args: argparse.Namespace) -> SCFConfig:
    return SCFConfig(
        method=args.method,
        basis=args.basis,
        ecp=args.ecp,
        dispersion=args.dispersion,
        grid_level=args.grid_level,
        conv_tol=args.conv_tol,
        max_cycle=args.max_cycle,
        max_memory=args.max_memory,
        density_fit=args.density_fit,
        auxbasis=args.auxbasis,
        unrestricted=args.unrestricted,
        level_shift=args.level_shift,
        damp=args.damp,
        init_guess=args.init_guess,
        linear_dep_threshold=args.linear_dep_threshold,
        validation_tol=args.validation_tol,
        strict_validation=args.strict_validation,
        verbose=args.verbose,
        solvent=args.solvent,
        solvent_eps=args.solvent_eps,
        solvent_name=args.solvent_name,
        solvent_options=_parse_solvent_options(args.solvent_option),
    )


def _normalize_energy_unit(text: str) -> str:
    return text.strip().lower().replace(" ", "")


def _build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Counterpoise-corrected PySCF density-matrix EDA. "
            "Run 'pyscf-dm-eda grid --help' for the checkpointed batch-grid runner."
        )
    )
    parser.add_argument("xyz", help="Supermolecule XYZ file")
    parser.add_argument(
        "--fragment",
        action="append",
        required=True,
        help="One-based atom selection; repeat for every fragment (for example 1-3)",
    )
    parser.add_argument("--fragment-charge", action="append", type=int, default=[])
    parser.add_argument(
        "--fragment-spin",
        action="append",
        type=int,
        default=[],
        help="Signed N_alpha-N_beta; repeat in fragment order",
    )
    parser.add_argument("--fragment-label", action="append", default=[])
    parser.add_argument("--charge", type=int, default=None)
    parser.add_argument("--spin", type=int, default=None)
    _add_scf_config_arguments(parser)
    parser.add_argument("--output", help="Optional JSON output path")
    parser.add_argument(
        "--unit",
        type=_normalize_energy_unit,
        default="kcal/mol",
        metavar="{hartree,kcal/mol,kj/mol,ev}",
    )
    return parser


def _build_grid_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pyscf-dm-eda grid",
        description="Checkpointed batch-grid DM-EDA runner",
    )
    parser.add_argument("molecule_xyz", help="Host molecule XYZ file")
    parser.add_argument("--molecule-charge", type=int, required=True)
    parser.add_argument("--molecule-spin", type=int, required=True)
    parser.add_argument("--probe-charge", type=int, required=True)
    parser.add_argument("--probe-spin", type=int, required=True)
    parser.add_argument("--supermolecule-charge", type=int, default=None)
    parser.add_argument("--supermolecule-spin", type=int, default=None)
    parser.add_argument("--probe-directory", default=None)
    parser.add_argument("--combined-pattern", default="mol_probe_{index}.xyz")
    parser.add_argument("--grid-xyz", default=None)
    parser.add_argument("--energy-output", default=None)
    parser.add_argument("--xyz-output", default=None)
    parser.add_argument("--geometry-tolerance", type=float, default=1e-4)
    parser.add_argument("--probe-anchor-atom", type=int, default=None)
    parser.add_argument("--probe-template-xyz", default=None)
    parser.add_argument(
        "--indices",
        default=None,
        help="Optional zero-based shard selection, e.g. 0-99 or 0-49,120",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Print one stderr progress line per completed grid point",
    )
    parser.add_argument("--restart", action="store_true")
    _add_scf_config_arguments(parser)
    return parser


def _report_runtime_error(exc: Exception) -> int:
    print(f"pyscf-dm-eda: error: {exc}", file=sys.stderr)
    return 1


def _single_geometry_main(parser: argparse.ArgumentParser, argv: Sequence[str]) -> int:
    args = parser.parse_args(argv)
    fragment_count = len(args.fragment)
    for option, values in (
        ("--fragment-charge", args.fragment_charge),
        ("--fragment-spin", args.fragment_spin),
        ("--fragment-label", args.fragment_label),
    ):
        if values and len(values) != fragment_count:
            parser.error(f"{option} must be omitted or repeated {fragment_count} times")
    try:
        atoms = read_xyz(args.xyz)
        charges = args.fragment_charge or [0] * fragment_count
        spins = args.fragment_spin or [0] * fragment_count
        labels = args.fragment_label or [
            f"fragment_{i + 1}" for i in range(fragment_count)
        ]
        fragments = [
            FragmentSpec(
                _parse_fragment_indices(selection, len(atoms)),
                charges[index],
                spins[index],
                labels[index],
            )
            for index, selection in enumerate(args.fragment)
        ]
        config = _scf_config_from_args(args)
        result = PySCFEDA(
            atoms,
            fragments,
            config,
            charge=args.charge,
            spin=args.spin,
        ).run()
        payload = json.dumps(result.as_dict(args.unit), indent=2, sort_keys=True)
        if args.output:
            output_path = Path(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            _atomic_write_text(output_path, payload + "\n")
        else:
            print(payload)
    except (argparse.ArgumentTypeError, TypeError, ValueError) as exc:
        parser.error(str(exc))
    except (EDAError, OSError) as exc:
        return _report_runtime_error(exc)
    return 0


def grid_main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point for the checkpointed batch-grid runner."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "grid":
        arguments = arguments[1:]
    parser = _build_grid_cli_parser()
    args = parser.parse_args(arguments)
    try:
        indices = (
            None if args.indices is None else _parse_index_subset(args.indices)
        )
        progress_callback = None
        if args.progress:
            def progress_callback(index: int, total: int) -> None:
                print(
                    f"grid point {index} checkpointed ({index + 1}/{total})",
                    file=sys.stderr,
                )

        runner = PySCFGridRunner(
            molecule_xyz=args.molecule_xyz,
            molecule_charge=args.molecule_charge,
            molecule_spin=args.molecule_spin,
            probe_charge=args.probe_charge,
            probe_spin=args.probe_spin,
            config=_scf_config_from_args(args),
            supermolecule_charge=args.supermolecule_charge,
            supermolecule_spin=args.supermolecule_spin,
            probe_directory=args.probe_directory,
            combined_pattern=args.combined_pattern,
            grid_xyz=args.grid_xyz,
            energy_output=args.energy_output,
            xyz_output=args.xyz_output,
            geometry_tolerance=args.geometry_tolerance,
            probe_anchor_atom=args.probe_anchor_atom,
            probe_template_xyz=args.probe_template_xyz,
            indices=indices,
            progress_callback=progress_callback,
        )
        summary = runner.run(restart=args.restart)
        print(
            json.dumps(
                {
                    "energy_table": str(summary.energy_table),
                    "extended_xyz": str(summary.extended_xyz),
                    "metadata_file": str(summary.metadata_file),
                    "completed_points": summary.completed_points,
                },
                sort_keys=True,
            )
        )
    except (argparse.ArgumentTypeError, TypeError, ValueError) as exc:
        parser.error(str(exc))
    except (EDAError, OSError) as exc:
        return _report_runtime_error(exc)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point for both single and grid calculations."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "grid":
        return grid_main(arguments[1:])
    return _single_geometry_main(_build_cli_parser(), arguments)


__all__ = [
    "Atom",
    "EDAError",
    "EDAResult",
    "EDAValidationError",
    "FragmentSpec",
    "GridRunSummary",
    "HARTREE_TO_KCAL_MOL",
    "HARTREE_TO_EV",
    "HARTREE_TO_KJ_MOL",
    "IncompatibleFragmentError",
    "PySCFEDA",
    "PySCFGridRunner",
    "SCFConfig",
    "SCFConvergenceError",
    "grid_main",
    "read_xyz",
]


if __name__ == "__main__":
    raise SystemExit(main())
