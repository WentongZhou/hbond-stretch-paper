import argparse
import io
import json
import os
import tempfile
import unittest
import warnings
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from pyscf_dm_eda import (  # noqa: E402
    Atom,
    EDAError,
    EDAResult,
    EDAValidationError,
    FragmentSpec,
    IncompatibleFragmentError,
    PySCFEDA,
    PySCFGridRunner,
    SCFConfig,
    read_xyz,
)
from pyscf_dm_eda.eda import (  # noqa: E402
    ENERGY_COMPONENTS,
    GRID_EDA_COLUMNS,
    _build_molecule,
    _exclusive_output_lock,
    _parse_fragment_indices,
    _parse_index_subset,
    _parse_solvent_options,
    _solvent_model,
    _validate_grid_geometry,
    _write_extended_xyz,
    main,
)

import numpy as np  # noqa: E402
from pyscf.solvent import pcm as pyscf_pcm  # noqa: E402


class PySCFEDATest(unittest.TestCase):
    def setUp(self):
        self.atoms = (Atom("He", (0.0, 0.0, 0.0)), Atom("He", (0.0, 0.0, 3.0)))
        self.fragments = (
            FragmentSpec((0,), label="A"),
            FragmentSpec((1,), label="B"),
        )

    def test_config_rejects_builtin_dispersion_functionals(self):
        for method in ("wb97x-d", "b97-d"):
            with self.subTest(method=method):
                with self.assertRaisesRegex(ValueError, "dispersion"):
                    SCFConfig(method=method)

    def test_atom_rejects_non_numeric_coordinates_with_value_error(self):
        with self.assertRaisesRegex(ValueError, "Invalid Cartesian"):
            Atom("He", ("a", "b", "c"))

    def test_density_fit_path_runs_and_reports_df_scheme(self):
        result = PySCFEDA(
            self.atoms,
            self.fragments,
            SCFConfig(
                method="hf",
                basis="sto-3g",
                density_fit=True,
                conv_tol=1e-11,
            ),
        ).run()
        values = result.components("hartree")

        self.assertIn("density-fitted", result.metadata["scheme"])
        self.assertEqual(result.metadata["jk_evaluation"], "density-fitted")
        self.assertLess(abs(values["Correlation Interaction"]), 1e-8)
        self.assertLess(abs(values["Closure Error"]), 1e-11)

    def test_dispersion_uses_real_atom_geometries_and_closes(self):
        def fake_dispersion(mol, config):
            return 0.25 if mol.natm == 1 else 1.0

        with patch(
            "pyscf_dm_eda.eda._dispersion_energy", side_effect=fake_dispersion
        ):
            result = PySCFEDA(
                self.atoms,
                self.fragments,
                SCFConfig(method="hf", basis="sto-3g", dispersion="d3"),
            ).run()
        values = result.components("hartree")

        self.assertAlmostEqual(values["Dispersion Interaction"], 0.5, places=12)
        self.assertAlmostEqual(values["Closure Error"], 0.0, places=11)

    def test_uks_broken_symmetry_path_closes(self):
        atoms = (Atom("H", (0.0, 0.0, -1.5)), Atom("H", (0.0, 0.0, 1.5)))
        fragments = (
            FragmentSpec((0,), spin=1, label="alpha"),
            FragmentSpec((1,), spin=-1, label="beta"),
        )
        result = PySCFEDA(
            atoms,
            fragments,
            SCFConfig(
                method="lda,vwn",
                basis="sto-3g",
                grid_level=1,
                conv_tol=1e-10,
            ),
            spin=0,
        ).run()

        self.assertTrue(result.metadata["unrestricted_supermolecule"])
        self.assertAlmostEqual(
            result.diagnostics["closure_error_hartree"], 0.0, places=11
        )

    def test_validation_failure_carries_result_and_can_downgrade(self):
        config = SCFConfig(
            method="hf",
            basis="sto-3g",
            conv_tol=1e-11,
            validation_tol=1e-15,
        )
        with self.assertRaises(EDAValidationError) as caught:
            PySCFEDA(self.atoms, self.fragments, config).run()
        exception = caught.exception
        self.assertIn("Total Interaction energy", exception.result.components())
        worst = exception.result.diagnostics["validation_worst_case"]
        self.assertGreater(
            exception.result.diagnostics["validation_errors"][worst],
            config.validation_tol,
        )
        self.assertIn("validation_errors", exception.result.diagnostics)

        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            result = PySCFEDA(
                self.atoms,
                self.fragments,
                SCFConfig(
                    method="hf",
                    basis="sto-3g",
                    conv_tol=1e-11,
                    validation_tol=1e-15,
                    strict_validation=False,
                ),
            ).run()
        self.assertEqual(len(caught_warnings), 1)
        self.assertIn("validation failed", str(caught_warnings[0].message))
        self.assertIn("Total Interaction energy", result.components())

    def test_units_include_ev_and_kj_per_mol_and_reject_unknown(self):
        result = EDAResult(
            {key: 1.0 for key in ENERGY_COMPONENTS},
            {"A": 0.25},
            {"closure_error_hartree": 2.0},
            {},
        )

        self.assertAlmostEqual(
            result.components("eV")["Total Interaction energy"], 27.211386245988
        )
        self.assertAlmostEqual(
            result.components("kJ/mol")["Total Interaction energy"],
            2625.4996394799,
        )
        with self.assertRaisesRegex(ValueError, "Unsupported energy unit"):
            result.components("furlong")
        self.assertEqual(result.grid_row(0, probe_fragment="A")["Mulliken_CT"], 0.25)
        with self.assertRaisesRegex(KeyError, "Unknown probe fragment"):
            result.grid_row(0, probe_fragment="missing")

    def test_write_json_uses_atomic_writer(self):
        result = EDAResult(
            {key: 1.0 for key in ENERGY_COMPONENTS},
            {"A": 0.25},
            {"closure_error_hartree": 0.0},
            {},
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "result.json"
            with patch("pyscf_dm_eda.eda._atomic_write_text") as atomic_write:
                result.write_json(output)
            atomic_write.assert_called_once()

    def test_fragment_index_parser_boundaries(self):
        self.assertEqual(
            _parse_fragment_indices("1-3,7", 7), (0, 1, 2, 6)
        )
        for selection in ("3-1", "0", "1,1", ""):
            with self.subTest(selection=selection):
                with self.assertRaises(argparse.ArgumentTypeError):
                    _parse_fragment_indices(selection, 3)
        self.assertEqual(_parse_index_subset("0-2,4"), (0, 1, 2, 4))
        with self.assertRaises(argparse.ArgumentTypeError):
            _parse_index_subset("2-1")

    def test_cli_single_geometry_end_to_end(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            xyz = Path(temporary_directory) / "he2.xyz"
            output = Path(temporary_directory) / "result.json"
            xyz.write_text("2\ncomment\nHe 0 0 0\nHe 0 0 3\n", encoding="utf-8")
            return_code = main(
                [
                    str(xyz),
                    "--fragment",
                    "1",
                    "--fragment",
                    "2",
                    "--method",
                    "hf",
                    "--basis",
                    "sto-3g",
                    "--unit",
                    "ev",
                    "--output",
                    str(output),
                ]
            )
            self.assertEqual(return_code, 0)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["energy_unit"], "ev")
            self.assertIn("Total Interaction energy", payload["components"])

    def test_cli_runtime_errors_exit_with_code_1(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            xyz = Path(temporary_directory) / "he2.xyz"
            xyz.write_text("2\ncomment\nHe 0 0 0\nHe 0 0 3\n", encoding="utf-8")
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                return_code = main(
                    [
                        str(xyz),
                        "--fragment",
                        "1",
                        "--fragment",
                        "2",
                        "--method",
                        "hf",
                        "--basis",
                        "not-a-basis",
                    ]
                )
            self.assertEqual(return_code, 1)
            self.assertIn("pyscf-dm-eda: error", stderr.getvalue())

    def test_cli_grid_subcommand_end_to_end(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            molecule = directory / "he.xyz"
            grid = directory / "he_filtered.xyz"
            probe_directory = directory / "he_probe"
            combined = probe_directory / "mol_probe_0.xyz"
            molecule.write_text("1\nhost\nHe 0 0 0\n", encoding="utf-8")
            grid.write_text("1\ngrid\nX 0 0 3\n", encoding="utf-8")
            probe_directory.mkdir()
            combined.write_text("2\ncombined\nHe 0 0 0\nHe 0 0 3\n", encoding="utf-8")
            stdout = io.StringIO()
            stderr = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                return_code = main(
                    [
                        "grid",
                        str(molecule),
                        "--molecule-charge",
                        "0",
                        "--molecule-spin",
                        "0",
                        "--probe-charge",
                        "0",
                        "--probe-spin",
                        "0",
                        "--method",
                        "hf",
                        "--basis",
                        "sto-3g",
                        "--indices",
                        "0",
                        "--progress",
                        "--energy-output",
                        str(directory / "eda.tsv"),
                    ]
                )
            self.assertEqual(return_code, 0)
            summary = json.loads(stdout.getvalue())
            self.assertEqual(summary["completed_points"], 1)
            self.assertIn("grid point 0 checkpointed", stderr.getvalue())
            self.assertTrue(Path(summary["energy_table"]).exists())


    def test_hf_decomposition_closes_and_has_zero_correlation(self):
        result = PySCFEDA(
            self.atoms,
            self.fragments,
            SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-11),
        ).run()
        values = result.components("hartree")

        self.assertAlmostEqual(
            values["Electrostatic Interaction"],
            values["Nuc---Nuc"] + values["1-electron"] + values["2-electron"],
            places=12,
        )
        self.assertAlmostEqual(
            values["Exchange-Repulsion"],
            values["Exchange Int."] + values["Repulsion"],
            places=12,
        )
        self.assertAlmostEqual(values["Correlation Interaction"], 0.0, places=12)
        self.assertAlmostEqual(values["Closure Error"], 0.0, places=12)
        self.assertAlmostEqual(sum(result.fragment_charge_transfer.values()), 0.0, places=12)
        self.assertLess(result.diagnostics["pauli_idempotency_max_abs"], 1e-12)

    def test_solvent_config_validation_and_normalisation(self):
        config = SCFConfig(
            method="hf", basis="sto-3g", solvent="pcm", solvent_name="Acetonitrile"
        )
        self.assertEqual(config.solvent, "C-PCM")
        self.assertEqual(config.solvent_name, "acetonitrile")
        self.assertEqual(SCFConfig(method="hf", solvent="SS(V)PE").solvent, "SS(V)PE")
        self.assertEqual(SCFConfig(method="hf", solvent="IEFPCM").solvent, "IEF-PCM")
        options = SCFConfig(
            method="hf", solvent="cpcm", solvent_options={"lebedev_order": 17}
        ).solvent_options
        self.assertEqual(options, {"lebedev_order": 17})
        with self.assertRaisesRegex(ValueError, "Unknown solvent model"):
            SCFConfig(method="hf", solvent="magic")
        with self.assertRaisesRegex(ValueError, "not supported"):
            SCFConfig(method="hf", solvent="ddcosmo")
        with self.assertRaisesRegex(ValueError, "requires a solvent model"):
            SCFConfig(method="hf", solvent_eps=80.0)
        with self.assertRaisesRegex(ValueError, ">= 1"):
            SCFConfig(method="hf", solvent="cpcm", solvent_eps=0.5)
        with self.assertRaises(TypeError):
            SCFConfig(method="hf", solvent="cpcm", solvent_eps="80")
        with self.assertRaisesRegex(ValueError, "Unknown solvent name"):
            SCFConfig(method="hf", solvent="smd", solvent_name="unobtainium")
        with self.assertRaisesRegex(ValueError, "SMD needs solvent_name"):
            SCFConfig(method="hf", solvent="smd", solvent_eps=40.0)
        with self.assertRaisesRegex(ValueError, "not settable"):
            SCFConfig(method="hf", solvent="cpcm", solvent_options={"eps": 20.0})
        with self.assertRaisesRegex(ValueError, "not settable"):
            SCFConfig(method="hf", solvent="cpcm", solvent_options={"bogus": 1})
        with self.assertRaises(TypeError):
            SCFConfig(method="hf", solvent="cpcm", solvent_options=["lebedev_order"])
        self.assertEqual(
            _parse_solvent_options(
                ["lebedev_order=17", "surface_discretization_method=ISWIG"]
            ),
            {"lebedev_order": 17, "surface_discretization_method": "ISWIG"},
        )
        self.assertIsNone(_parse_solvent_options([]))
        with self.assertRaisesRegex(ValueError, "KEY=VALUE"):
            _parse_solvent_options(["lebedev_order"])

    def test_fragment_cavity_excludes_ghost_atoms(self):
        # PySCF gives a ghost atom the Z=0 radius of its table (2.0 Angstrom),
        # so a naive fragment PCM would be solvated in an oversized cavity.
        config = SCFConfig(method="hf", basis="sto-3g", solvent="cpcm")
        ghost_mol = _build_molecule(
            self.atoms, config, charge=0, spin=0, real_indices=(0,)
        )
        physical_mol = _build_molecule(
            self.atoms[:1], config, charge=0, spin=0, real_indices=(0,)
        )
        solvent = _solvent_model(ghost_mol, physical_mol, config)
        solvent.build()
        reference = pyscf_pcm.PCM(physical_mol)
        reference.build()
        naive = pyscf_pcm.PCM(ghost_mol)
        naive.build()

        points = solvent.surface["grid_coords"].shape[0]
        self.assertEqual(points, reference.surface["grid_coords"].shape[0])
        self.assertGreater(naive.surface["grid_coords"].shape[0], points)
        np.testing.assert_allclose(solvent.v_grids_n, reference.v_grids_n)
        # AO integrals still run in the full ghost basis.
        self.assertIs(solvent.mol, ghost_mol)
        supermolecule = _solvent_model(ghost_mol, ghost_mol, config)
        self.assertIsNone(supermolecule.cavity_mol)

    def test_pcm_decomposition_closes_and_isolates_desolvation(self):
        gas = PySCFEDA(
            self.atoms,
            self.fragments,
            SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-11),
        ).run()
        solvated = PySCFEDA(
            self.atoms,
            self.fragments,
            SCFConfig(
                method="hf",
                basis="sto-3g",
                conv_tol=1e-11,
                solvent="cpcm",
                solvent_eps=78.3553,
            ),
        ).run()
        gas_values = gas.components("hartree")
        values = solvated.components("hartree")

        self.assertIn("Desolvation", ENERGY_COMPONENTS)
        self.assertIn("Desolv", GRID_EDA_COLUMNS)
        self.assertEqual(gas_values["Desolvation"], 0.0)
        self.assertIsNone(gas.diagnostics["solvent"])
        self.assertIsNone(gas.metadata["solvent"])
        self.assertEqual(gas.metadata["limitations"][0], "Gas phase only; no PCM/desolvation term.")

        solvent = solvated.diagnostics["solvent"]
        self.assertNotEqual(solvent["supermolecule_solvent_energy_hartree"], 0.0)
        self.assertNotEqual(values["Desolvation"], 0.0)
        self.assertNotAlmostEqual(
            values["Total Interaction energy"],
            gas_values["Total Interaction energy"],
            places=9,
        )
        # HF: E0 + E_solv accounts for the whole SCF energy, so the residual
        # must stay zero once the reaction field is removed from it.
        self.assertAlmostEqual(values["Correlation Interaction"], 0.0, places=11)
        self.assertAlmostEqual(values["Closure Error"], 0.0, places=11)
        self.assertAlmostEqual(
            values["Total Interaction energy"],
            values["Electrostatic Interaction"]
            + values["Exchange Int."]
            + values["Repulsion"]
            + values["Orbital Relaxation"]
            + values["Correlation Interaction"]
            + values["Dispersion Interaction"]
            + values["Desolvation"],
            places=11,
        )
        self.assertAlmostEqual(
            solvent["desolvation_frozen_hartree"]
            + solvent["desolvation_pauli_response_hartree"]
            + solvent["desolvation_polarization_response_hartree"],
            values["Desolvation"],
            places=12,
        )
        self.assertEqual(solvent["supermolecule_cds_energy_hartree"], 0.0)
        points = solvent["cavity_surface_points"]
        self.assertEqual(points["A"], points["B"])
        self.assertLess(points["A"], points["supermolecule"])
        self.assertEqual(solvated.metadata["solvent"]["model"], "C-PCM")
        self.assertEqual(solvated.metadata["solvent"]["eps"], 78.3553)
        self.assertIn("Implicit solvent", solvated.metadata["limitations"][0])

        row = solvated.grid_row(0, probe_fragment="B")
        self.assertAlmostEqual(
            row["Desolv"], solvated.components("kcal/mol")["Desolvation"], places=12
        )
        self.assertEqual(set(row), set(GRID_EDA_COLUMNS))
        payload = solvated.as_dict("kcal/mol")
        json.dumps(payload)
        self.assertEqual(payload["metadata"]["solvent"]["options"], {})

    def test_smd_reports_cds_and_closes(self):
        result = PySCFEDA(
            self.atoms,
            self.fragments,
            SCFConfig(
                method="hf",
                basis="sto-3g",
                conv_tol=1e-11,
                solvent="smd",
                solvent_name="water",
            ),
        ).run()
        values = result.components("hartree")
        solvent = result.diagnostics["solvent"]

        self.assertAlmostEqual(values["Correlation Interaction"], 0.0, places=11)
        self.assertAlmostEqual(values["Closure Error"], 0.0, places=11)
        self.assertNotEqual(values["Desolvation"], 0.0)
        self.assertIn("supermolecule_cds_energy_hartree", solvent)
        self.assertIn("A", solvent["fragment_cds_energy_hartree"])
        self.assertEqual(result.metadata["solvent"]["model"], "SMD")
        self.assertEqual(result.metadata["solvent"]["solvent_name"], "water")
        self.assertEqual(
            result.metadata["solvent"]["implementation"], "pyscf.solvent.smd"
        )
        self.assertAlmostEqual(result.metadata["solvent"]["eps"], 78.355, places=6)

    def test_cli_accepts_solvent_options(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            xyz = Path(temporary_directory) / "he2.xyz"
            output = Path(temporary_directory) / "result.json"
            xyz.write_text("2\ncomment\nHe 0 0 0\nHe 0 0 3\n", encoding="utf-8")
            with redirect_stdout(io.StringIO()):
                return_code = main(
                    [
                        str(xyz),
                        "--fragment",
                        "1",
                        "--fragment",
                        "2",
                        "--method",
                        "hf",
                        "--basis",
                        "sto-3g",
                        "--solvent",
                        "iefpcm",
                        "--solvent-name",
                        "acetonitrile",
                        "--solvent-option",
                        "lebedev_order=17",
                        "--output",
                        str(output),
                    ]
                )
            self.assertEqual(return_code, 0)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["metadata"]["solvent"]["model"], "IEF-PCM")
            self.assertAlmostEqual(payload["metadata"]["solvent"]["eps"], 35.688)
            self.assertEqual(
                payload["metadata"]["solvent"]["options"], {"lebedev_order": 17}
            )
            self.assertIn("Desolvation", payload["components"])

            with redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as context:
                    main(
                        [
                            str(xyz),
                            "--fragment",
                            "1",
                            "--fragment",
                            "2",
                            "--solvent",
                            "ddcosmo",
                        ]
                    )
            self.assertEqual(context.exception.code, 2)

    def test_dft_residual_is_reported_once_and_closes(self):
        result = PySCFEDA(
            self.atoms,
            self.fragments,
            SCFConfig(
                method="lda,vwn",
                basis="sto-3g",
                grid_level=1,
                conv_tol=1e-10,
            ),
        ).run()
        values = result.components("hartree")

        self.assertNotAlmostEqual(values["Correlation Interaction"], 0.0, places=10)
        self.assertAlmostEqual(values["Dispersion Interaction"], 0.0, places=12)
        self.assertAlmostEqual(values["Closure Error"], 0.0, places=11)

    def test_spin_projection_must_match_fragment_sum(self):
        fragments = (
            FragmentSpec((0,), spin=1, label="alpha"),
            FragmentSpec((1,), spin=1, label="also_alpha"),
        )
        with self.assertRaisesRegex(IncompatibleFragmentError, "negative spin"):
            PySCFEDA(self.atoms, fragments, spin=0)

    def test_opposite_fragment_spins_select_unrestricted_supermolecule(self):
        atoms = (Atom("H", (0.0, 0.0, -1.5)), Atom("H", (0.0, 0.0, 1.5)))
        fragments = (
            FragmentSpec((0,), spin=1, label="alpha"),
            FragmentSpec((1,), spin=-1, label="beta"),
        )
        result = PySCFEDA(
            atoms,
            fragments,
            SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-10),
            spin=0,
        ).run()

        self.assertTrue(result.metadata["unrestricted_supermolecule"])
        self.assertAlmostEqual(
            result.diagnostics["supermolecule_electrons_from_density"]["alpha"],
            1.0,
            places=10,
        )
        self.assertAlmostEqual(
            result.diagnostics["supermolecule_electrons_from_density"]["beta"],
            1.0,
            places=10,
        )
        self.assertAlmostEqual(
            result.diagnostics["closure_error_hartree"], 0.0, places=11
        )

    def test_xyz_reader_and_json_result(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            xyz = Path(temporary_directory) / "he2.xyz"
            xyz.write_text("2\ncomment\nHe 0 0 0\nHe 0 0 3\n", encoding="utf-8")
            self.assertEqual(read_xyz(xyz), self.atoms)

            result = PySCFEDA.from_xyz(
                xyz,
                self.fragments,
                SCFConfig(method="hf", basis="sto-3g"),
            ).run()
            output = Path(temporary_directory) / "eda.json"
            result.write_json(output)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["metadata"]["scheme"], "DM-EDA (exact J/K)")
            self.assertIn("Closure Error", payload["components"])

    def test_grid_checkpoint_restart_and_input_fingerprint(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            molecule = directory / "he.xyz"
            grid = directory / "he_filtered.xyz"
            probe_directory = directory / "he_probe"
            combined = probe_directory / "complex_000.xyz"
            molecule.write_text("1\nhost\nHe 0 0 0\n", encoding="utf-8")
            grid.write_text("1\ngrid\nX 0 0 3\n", encoding="utf-8")
            probe_directory.mkdir()
            combined.write_text(
                "2\ncombined\nHe 0 0 0\nHe 0 0 3\n", encoding="utf-8"
            )
            progress = []
            runner = PySCFGridRunner(
                molecule_xyz=molecule,
                molecule_charge=0,
                molecule_spin=0,
                probe_charge=0,
                probe_spin=0,
                combined_pattern="complex_{index:03d}.xyz",
                config=SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-11),
                indices=(0,),
                progress_callback=lambda index, total: progress.append((index, total)),
            )

            first = runner.run()
            restarted = runner.run(restart=True)
            self.assertEqual(first.completed_points, 1)
            self.assertEqual(restarted.completed_points, 1)
            self.assertEqual(progress, [(0, 1)])
            checkpoint_text = first.energy_table.read_text(encoding="utf-8")
            self.assertEqual(len(checkpoint_text.splitlines()), 3)
            self.assertNotIn("\r\n", checkpoint_text)
            metadata = json.loads(first.metadata_file.read_text(encoding="utf-8"))
            self.assertEqual(metadata["schema"], "pyscf-dm-eda-grid-v1")
            self.assertEqual(metadata["indices"], [0])
            self.assertIn("implementation_sha256", metadata)
            self.assertIn("combined_geometries", metadata["input_sha256"])

            checkpoint_lines = checkpoint_text.splitlines()
            checkpoint_lines[0] = (
                "# pyscf-dm-eda-grid-v1 fingerprint=invalid"
            )
            first.energy_table.write_text(
                "\n".join(checkpoint_lines) + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(EDAError, "fingerprint"):
                runner.run(restart=True)
            first.energy_table.write_text(checkpoint_text, encoding="utf-8")

            combined.write_text(
                "2\nchanged\nHe 0 0 0\nHe 0 0 4\n", encoding="utf-8"
            )
            grid.write_text("1\ngrid\nX 0 0 4\n", encoding="utf-8")
            with patch(
                "pyscf_dm_eda.eda._atomic_write_text",
                side_effect=OSError("simulated metadata failure"),
            ):
                with self.assertRaisesRegex(OSError, "simulated metadata failure"):
                    runner.run(restart=False)
            with self.assertRaisesRegex(EDAError, "Restart settings differ"):
                runner.run(restart=True)

    def test_polyatomic_probe_requires_and_uses_an_unambiguous_origin(self):
        molecule_atoms = (Atom("He", (0.0, 0.0, 0.0)),)
        template_atoms = (
            Atom("C", (0.0, 0.0, 0.0)),
            Atom("O", (0.0, 0.0, -1.4)),
        )
        combined_atoms = molecule_atoms + (
            Atom("C", (0.0, 0.0, 3.0)),
            Atom("O", (0.0, 0.0, 1.6)),
        )
        combined_path = Path("mol_probe_0.xyz")

        with self.assertRaisesRegex(EDAError, "requires probe_anchor_atom"):
            _validate_grid_geometry(
                molecule_atoms,
                combined_atoms,
                ("X", "0", "0", "3"),
                combined_path,
                1e-6,
                None,
                None,
            )
        _validate_grid_geometry(
            molecule_atoms,
            combined_atoms,
            ("X", "0", "0", "3"),
            combined_path,
            1e-6,
            None,
            template_atoms,
        )
        with self.assertRaisesRegex(EDAError, "template origin"):
            _validate_grid_geometry(
                molecule_atoms,
                combined_atoms,
                ("X", "0", "0", "2.5"),
                combined_path,
                1e-6,
                None,
                template_atoms,
            )

    def test_grid_output_lock_rejects_a_second_writer(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            lock_path = Path(temporary_directory) / "grid.tsv.lock"
            with _exclusive_output_lock(lock_path):
                with self.assertRaisesRegex(EDAError, "holds the output lock"):
                    with _exclusive_output_lock(lock_path):
                        self.fail("A second writer acquired the same output lock")

    def test_grid_shard_indices_are_validated(self):
        runner = PySCFGridRunner(
            molecule_xyz="molecule.xyz",
            molecule_charge=0,
            molecule_spin=0,
            probe_charge=0,
            probe_spin=0,
            indices=(2,),
        )
        with self.assertRaisesRegex(ValueError, "only 2 points"):
            runner._selected_indices(2)
        with self.assertRaisesRegex(ValueError, "duplicates"):
            PySCFGridRunner(
                molecule_xyz="molecule.xyz",
                molecule_charge=0,
                molecule_spin=0,
                probe_charge=0,
                probe_spin=0,
                indices=(0, 0),
            )

    def test_extended_xyz_writes_only_checkpointed_shard_rows(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "shard.xyz"
            rows = [("X", "0", "0", "3"), ("X", "0", "0", "4")]
            values = {column: "0.0" for column in GRID_EDA_COLUMNS}
            _write_extended_xyz(path, rows, {1: values})
            lines = path.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[0], "1")
            self.assertTrue(lines[-1].startswith("X 0 0 4"))

    def test_combined_pattern_cannot_escape_probe_directory(self):
        with self.assertRaisesRegex(ValueError, "inside probe_directory"):
            PySCFGridRunner(
                molecule_xyz="molecule.xyz",
                molecule_charge=0,
                molecule_spin=0,
                probe_charge=0,
                probe_spin=0,
                combined_pattern="../complex_{index}.xyz",
            )


if __name__ == "__main__":
    unittest.main()


class PointChargeEmbeddingTest(unittest.TestCase):
    """Electrostatic embedding keeps the algebraic identities and responds to the field."""

    def setUp(self):
        self.atoms = [
            Atom("O", (0.0, 0.0, 0.0)), Atom("H", (0.96, 0.0, 0.0)), Atom("H", (-0.24, 0.93, 0.0)),
            Atom("O", (2.9, 0.0, 0.0)), Atom("H", (3.2, -0.5, 0.76)), Atom("H", (3.2, -0.5, -0.76)),
        ]
        self.fragments = [FragmentSpec((0, 1, 2), 0, 0, "donor"), FragmentSpec((3, 4, 5), 0, 0, "acceptor")]
        self.charges = ([(1.5, 3.5, 0.5), (1.5, 3.9, 0.5), (1.5, -3.2, 0.2)], [-0.8, 0.4, 0.4])

    def test_closure_and_frozen_identity_hold_with_point_charges(self):
        result = PySCFEDA(
            self.atoms, self.fragments,
            SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-11, point_charges=self.charges),
        ).run()
        self.assertEqual(result.diagnostics["n_point_charges"], 3)
        self.assertLess(abs(result.diagnostics["closure_error_hartree"]), 1e-10)
        self.assertLess(abs(result.diagnostics["frozen_identity_error_hartree"]), 1e-10)
        self.assertLess(result.diagnostics["core_partition_max_abs"], 1e-10)

    def test_field_changes_the_interaction_but_not_the_gas_phase_limit(self):
        gas = PySCFEDA(self.atoms, self.fragments, SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-11)).run()
        far = PySCFEDA(
            self.atoms, self.fragments,
            SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-11,
                      point_charges=([(500.0, 0.0, 0.0)], [1.0])),
        ).run()
        near = PySCFEDA(
            self.atoms, self.fragments,
            SCFConfig(method="hf", basis="sto-3g", conv_tol=1e-11, point_charges=self.charges),
        ).run()
        e_gas = gas.components("hartree")["Total Interaction energy"]
        e_far = far.components("hartree")["Total Interaction energy"]
        e_near = near.components("hartree")["Total Interaction energy"]
        self.assertAlmostEqual(e_gas, e_far, places=6)
        self.assertGreater(abs(e_near - e_gas), 1e-5)

    def test_point_charges_are_validated(self):
        with self.assertRaisesRegex(ValueError, "one charge per coordinate"):
            SCFConfig(method="hf", point_charges=([(0.0, 0.0, 0.0)], [1.0, 2.0]))
        with self.assertRaisesRegex(ValueError, "implicit solvent"):
            SCFConfig(method="hf", solvent="pcm", point_charges=([(0.0, 0.0, 0.0)], [1.0]))
