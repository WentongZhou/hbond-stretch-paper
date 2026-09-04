"""Water-dimer DM-EDA in the gas phase and in PySCF implicit solvent.

The same hydrogen-bonded water dimer is decomposed three times:

* gas phase,
* C-PCM water (eps = 78.3553),
* SMD water (PCM electrostatics + CDS).

Every fragment SCF and the supermolecule SCF carry the solvent model.  The
fragments are solvated in their own real-atom cavities (ghost atoms never get
a cavity sphere), and the change of the reaction-field energy is reported as
``Desolvation``.  The other terms keep their gas-phase operator definitions
but are evaluated with the solvated densities, so the printed table shows how
solvation screens electrostatics and shifts polarisation into the solvent
response.

Run with the edaenv interpreter:

    PYTHONPATH=src ~/edaenv/bin/python examples/water_dimer_solvent_eda.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from pyscf_dm_eda import Atom, FragmentSpec, PySCFEDA, SCFConfig

ATOMS = (
    Atom("O", (-1.551007, -0.114520, 0.000000)),
    Atom("H", (-1.934259, 0.762503, 0.000000)),
    Atom("H", (-0.599677, 0.040712, 0.000000)),
    Atom("O", (1.350625, 0.111469, 0.000000)),
    Atom("H", (1.680398, -0.373741, -0.758561)),
    Atom("H", (1.680398, -0.373741, 0.758561)),
)
FRAGMENTS = (
    FragmentSpec((0, 1, 2), label="donor"),
    FragmentSpec((3, 4, 5), label="acceptor"),
)
COMPONENTS = (
    "Total Interaction energy",
    "Electrostatic Interaction",
    "Exchange Int.",
    "Repulsion",
    "Orbital Relaxation",
    "Correlation Interaction",
    "Dispersion Interaction",
    "Desolvation",
    "Closure Error",
)
SOLVENT_SPLIT = (
    ("frozen (cavity + screening)", "desolvation_frozen_hartree"),
    ("Pauli response", "desolvation_pauli_response_hartree"),
    ("orbital-relaxation response", "desolvation_polarization_response_hartree"),
)
HARTREE_TO_KCAL_MOL = 627.5094740631


def main(argv: list[str]) -> int:
    method = argv[1] if len(argv) > 1 else "pbe"
    basis = argv[2] if len(argv) > 2 else "def2-svp"
    common = dict(method=method, basis=basis, grid_level=3)
    runs = {
        "gas": SCFConfig(**common),
        "C-PCM water": SCFConfig(**common, solvent="cpcm"),
        "SMD water": SCFConfig(**common, solvent="smd", solvent_name="water"),
    }

    results = {}
    for label, config in runs.items():
        results[label] = PySCFEDA(ATOMS, FRAGMENTS, config).run()

    width = max(len(name) for name in COMPONENTS)
    header = f"{'component (kcal/mol)':<{width}}" + "".join(
        f"{label:>16}" for label in runs
    )
    print(header)
    print("-" * len(header))
    for name in COMPONENTS:
        row = f"{name:<{width}}"
        for label in runs:
            row += f"{results[label].components('kcal/mol')[name]:16.3f}"
        print(row)

    print()
    print("desolvation split along the DM-EDA density sequence (kcal/mol):")
    for label, result in results.items():
        solvent = result.diagnostics["solvent"]
        if solvent is None:
            continue
        parts = ", ".join(
            f"{name} = {solvent[key] * HARTREE_TO_KCAL_MOL:.3f}"
            for name, key in SOLVENT_SPLIT
        )
        print(f"  {label}: {parts}")

    output_directory = Path(__file__).resolve().parent / "water_dimer_solvent"
    output_directory.mkdir(exist_ok=True)
    for label, result in results.items():
        stem = label.lower().replace(" ", "_").replace("-", "")
        result.write_json(output_directory / f"water_dimer_{stem}.json")
    summary = {
        label: {name: result.components("kcal/mol")[name] for name in COMPONENTS}
        for label, result in results.items()
    }
    (output_directory / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"\nJSON results written to {output_directory}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
