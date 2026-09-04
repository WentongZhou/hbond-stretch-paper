# hbond-stretch-eda

Code, scripts, geometries and data for

> Wentong Zhou, *A Local-Coordinate Analysis of Water's Hydrogen-Bond Stretch: Tracking the Exchange–Repulsion Wall*,
> submitted to *J. Phys. Chem. Lett.* (2026).

The paper partitions the curvature of the intermolecular hydrogen-bond stretch of water along R(O···O) into
density-matrix energy-decomposition (DM-EDA) channels — electrostatics, exchange–repulsion, orbital relaxation,
and correlation plus dispersion — for the water dimer, small ion-containing clusters, and 40 hydrogen bonds cut
from classical molecular-dynamics snapshots of liquid water, acidic water (H₃O⁺/Cl⁻) and basic water (OH⁻/Na⁺).

## Layout

| Path | Content |
|---|---|
| `pyscf-dm-eda/` | The DM-EDA code (PySCF-based), snapshot of commit `63424fa` **plus** the point-charge electrostatic-embedding extension written for this work (`src/pyscf_dm_eda/eda.py`, tests in `tests/test_eda.py`). MIT licence (see `pyscf-dm-eda/LICENSE`). |
| `pyscf-dm-eda_point_charges.patch` | The embedding extension as a diff against commit `63424fa`, for reference; it is already applied in the snapshot. |
| `scripts/` | The calculation pipeline (numbered stages, see below) and the shell drivers used to run it locally (WSL) and on Google Compute Engine VMs. |
| `geometries/` | Optimized dimer and cluster geometries (`*.xyz`), and the liquid-snapshot clusters (`liquid/`, `liquid_h3o/`, `liquid_oh/`; `*.json` = cluster specifications with fragment definitions and embedding charges). |
| `results/` | Every calculation output in JSON (`scan_*.json` DM-EDA scans, `fc_*.json` force-constant fits, `ccsdt_*.json`, `density_*.json`, `cluster_*.json`, `wp4_*.json`, `liquid_*_analysis.json`, `results/liquid/` cluster specs) and the Markdown tables generated from them. |
| `logs/` | Stdout/stderr of every job (`*.out`, `*.log`). |
| `figures/` | Diagnostic figures written by the pipeline (one per scan/force-constant fit, density profiles, environment comparisons). |
| `paper/` | Analysis and figure scripts for the manuscript: `revision_stats.py` (all derived statistics quoted in the text and the Supporting Information), `make_paper_figures_v4.py` (Figures 1–3), `make_si_figures_v4.py` (Figures S1–S4), `paper/data/` (the subset of `results/` used by the paper, grouped by type; produced by `collect_materials.py`), `paper/tables/`, and the final figure files in `paper/figures_paper/v4/`. |
| `requirements.txt` | `pip freeze` of the environment used for every calculation. |

## Environment

All calculations were run on Ubuntu 24.04 (WSL2 and Google Compute Engine n2-highmem-32 VMs) with
Python 3.12.3, PySCF 2.14.0, pyscf-dispersion 1.5.0 (D3(BJ)/D4), NumPy 2.5.2, SciPy 1.18.1, OpenMM 8.6.0
(liquid MD, rigid TIP4P-Ew) and ASE 3.29.0 (geometry optimization driver). Exact versions are in `requirements.txt`.

```bash
python -m venv edaenv && source edaenv/bin/activate
pip install -r requirements.txt
pip install -e ./pyscf-dm-eda          # installs the DM-EDA package with the embedding extension
cd pyscf-dm-eda && python -m pytest tests -q && cd ..
```

## Pipeline

Each stage is a script in `scripts/`; `scripts/common.py` holds the shared helpers (fragment handling,
DM-EDA driver call, JSON I/O). Paths are relative to the repository root.

| Stage | Script | What it does |
|---|---|---|
| 00 | `00_optimize_dimer.py` | Optimize the water dimer (revPBE0-D3(BJ)/def2-TZVP, ASE BFGS driver). |
| 01 | `01_scan_dimer.py` | Rigid scan of the acceptor along O···O with DM-EDA at every point (`results/scan_<tag>.json`). |
| 02 | `02_force_constants.py` | Fit each DM-EDA component around the minimum and partition k = d²E/dR² (`results/fc_<tag>.json`, `figures/fc_<tag>.png`). |
| 03 | `03_density_profile.py` | Stage-resolved electron density along the O–H···O axis (`results/density_<tag>.json`). |
| 04 | `04_summary_tables.py` | Markdown summary tables from all `fc_*.json` / `density_*.json`. |
| 05 | `05_overlap_proxies.py` | Decay lengths of charge transfer, exchange–repulsion and the axial density minimum. |
| 06 | `06_build_clusters.py` | Build and optimize the gas-phase clusters (chain trimer, H₃O⁺(H₂O)₄, OH⁻(H₂O)₅). |
| 07 | `07_scan_cluster.py` | Rigid scan of one target hydrogen bond inside a cluster (full cluster, bare pair, or point-charge-embedded pair). |
| 08 | `08_ccsdt_scan.py` | Counterpoise-corrected HF/MP2/CCSD/CCSD(T) scan of the dimer (aug-cc-pVTZ/QZ). |
| 09 | `09_compare_environments.py` | Compare the force-constant partition across environments (`results/wp4_environments.json`). |
| 10, 10b | `10_liquid_md.py`, `10b_ion_md.py` | OpenMM MD of neutral, acidic and basic TIP4P-Ew water (NPT, 1 bar). |
| 11, 11b | `11_extract_hbonds.py`, `11b_extract_ion_hbonds.py` | Select water–water hydrogen bonds from snapshots and cut clusters (`geometries/liquid*/`). |
| 12 | `12_liquid_analysis.py` | Statistics of the partition over the liquid hydrogen bonds (`results/liquid_*_analysis.json`). |
| 13 | `13_add_embedding.py` | Add the electrostatic-embedding charge set to the liquid cluster specifications. |
| 14 | `14_liquid_compare.py` | Compare neutral, acidic and basic liquid statistics (`results/liquid_compare.md`). |

`scripts/run_wsl.sh` runs stages 00–03 for one method (`XC`, `DISP`, `BASIS`, `TAG`, `STAGES`, `RMIN`/`RMAX`/`STEP`,
`OH_ELONG` environment variables); `wp4_*.sh`, `local_jobs_*.sh` and `cloud_jobs_*.sh` are the batch drivers that
were used for the cluster and liquid work. The `cloud_*.sh` scripts manage the Compute Engine VMs; the Google Cloud
project ids have been replaced by placeholders (`PROJECT=…`), and the comments still show the original Windows/WSL
paths of the author's machine.

## Data formats

* `results/scan_<tag>.json`: metadata (`xc`, `dispersion`, `basis`, `grid_level`, `oh_elongation_angstrom`,
  `reference_R_OO`) and `points`, one entry per R with the DM-EDA terms in kcal mol⁻¹ (`Total`, `Elec`, `ExRep` =
  `Exch` + `Rep`, `OrbRel`, `Corr`, `Disp`, `CorrDisp`, `Steric`, `Closure`), `mulliken_ct` and `iao_ct`
  (charge-transfer proxies, electrons), `R_OO` (Å), `E_super_hartree` and timing.
* `results/fc_<tag>.json`: the fit window, `R_min_A`, `E_int_min_kcal`, `k_total_kcal_per_A2`,
  `omega_H2O_cm-1` (harmonic wavenumber for the reduced mass of two water molecules), `charge_transfer_acceptor`,
  `max_abs_closure_hartree`, and `components[<channel>]` with `value_kcal`, `slope_kcal_per_A`, `k_kcal_per_A2`,
  `cubic_kcal_per_A3`, `fraction_of_k` (= k_X / k_total).
* `results/liquid/cluster_<name>.json` and `geometries/liquid*/*.json`: cluster specification (`system`, `charge`,
  `xyz`, `fragment_sets`, `descriptors`, `embedding` charges).
* `results/liquid_<system>_analysis.json`: `rows` (one hydrogen bond each, with `full`, `pair` and `embed` fits)
  and `fits` (pooled log–log fits).

## Reproducing the paper's numbers and figures

```bash
cd paper
python revision_stats.py             # statistics quoted in the text / SI
python make_paper_figures_v4.py      # Figures 1–3 and the data-only TOC graphic -> figures_paper/v4/
python make_si_figures_v4.py         # Figures S1–S4
```

`paper/data/` is a copy of the relevant `results/` files grouped by type (`scans/`, `fc/`, `ccsdt/`, `density/`,
`analysis/`, `liquid_clusters/`, `optimization/`); `make_paper_figures_v4.py` checks that the copy it reads is
identical to the original output when the source tree is present.

## Licence and citation

The DM-EDA code is released under the MIT licence (`pyscf-dm-eda/LICENSE`). The licence for the scripts and data in
this repository will be stated here by the author. Please cite the paper above when using these data.
