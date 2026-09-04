"""Classical MD of liquid water with OpenMM (rigid TIP4P-Ew, SETTLE, PME, NPT at 1 bar) to
generate liquid snapshots for WP6.

A cubic water box of the requested edge is built with Modeller.addSolvent,
minimized, thermalized with a Langevin integrator (300 K, 2 fs) and, after
equilibration, snapshots are written every ``--save-every`` ps to
geometries/liquid/snap_<index>.xyz as extended xyz (O,H,H per molecule, the
virtual site dropped, unwrapped coordinates, cell).  A final O-O radial
distribution function is written for sanity checking.

Usage:
    python scripts/10_liquid_md.py --box 15.5 --equil-ps 50 --prod-ps 40 --save-every 4 --seed 7
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def oo_rdf(frames, L, rmax=6.0, nbins=120):
    """O-O RDF; L may be a scalar or one box length per frame."""
    edges = np.linspace(0, rmax, nbins + 1)
    hist = np.zeros(nbins)
    n_o = 0
    Ls = [L] * len(frames) if np.isscalar(L) else list(L)
    for pos, L in zip(frames, Ls):
        o = pos[0::3]
        n_o = len(o)
        d = o[:, None, :] - o[None, :, :]
        d -= L * np.round(d / L)
        r = np.linalg.norm(d, axis=-1)[np.triu_indices(n_o, 1)]
        hist += np.histogram(r, bins=edges)[0]
    rho = n_o / float(np.mean(Ls)) ** 3
    shell = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = hist / len(frames) / (0.5 * n_o * rho * shell)
    return 0.5 * (edges[1:] + edges[:-1]), g


def write_extxyz(path, symbols, pos, L, comment):
    lines = [str(len(symbols)),
             f'Lattice="{L:.6f} 0.0 0.0 0.0 {L:.6f} 0.0 0.0 0.0 {L:.6f}" Properties=species:S:1:pos:R:3 pbc="T T T" {comment}']
    for s, (x, y, z) in zip(symbols, pos):
        lines.append(f"{s} {x:.6f} {y:.6f} {z:.6f}")
    Path(path).write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--box", type=float, default=15.5, help="cubic box edge in Angstrom")
    parser.add_argument("--model", default="tip4pew", choices=["tip4pew", "tip3p"])
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--dt-fs", type=float, default=2.0)
    parser.add_argument("--equil-ps", type=float, default=50.0)
    parser.add_argument("--prod-ps", type=float, default=40.0)
    parser.add_argument("--save-every", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--outdir", default=str(C.GEOM_DIR / "liquid"))
    args = parser.parse_args()

    import openmm
    from openmm import app, unit

    ff = app.ForceField(f"{args.model}.xml")
    modeller = app.Modeller(app.Topology(), [])
    L_nm = args.box / 10.0
    modeller.addSolvent(ff, model=args.model, boxSize=openmm.Vec3(L_nm, L_nm, L_nm) * unit.nanometer)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.7 * unit.nanometer,
                             constraints=app.HBonds, rigidWater=True)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.bar, args.temperature * unit.kelvin, 25))
    integrator = openmm.LangevinMiddleIntegrator(args.temperature * unit.kelvin, 1.0 / unit.picosecond,
                                                 args.dt_fs * unit.femtosecond)
    integrator.setRandomNumberSeed(args.seed)
    platform = openmm.Platform.getPlatformByName("CPU")
    sim = app.Simulation(modeller.topology, system, integrator, platform, {"Threads": str(args.threads)})
    sim.context.setPositions(modeller.positions)
    sim.context.setPeriodicBoxVectors(*[openmm.Vec3(*v) for v in np.eye(3) * L_nm])
    n_water = sum(1 for r in modeller.topology.residues())
    # real (non-virtual) atoms in O,H,H order per residue
    real_idx, symbols = [], []
    for res in modeller.topology.residues():
        atoms = list(res.atoms())
        o = [a for a in atoms if a.element is not None and a.element.symbol == "O"]
        h = [a for a in atoms if a.element is not None and a.element.symbol == "H"]
        for a in o + h:
            real_idx.append(a.index)
            symbols.append(a.element.symbol)
    real_idx = np.array(real_idx)
    print(f"# {n_water} waters ({args.model}), box {args.box:.2f} A, {len(real_idx)} real atoms, "
          f"density {n_water * 18.01528 / 0.602214076 / args.box ** 3:.3f} g/cm3", flush=True)

    t0 = time.time()
    sim.minimizeEnergy(maxIterations=500)
    sim.context.setVelocitiesToTemperature(args.temperature * unit.kelvin, args.seed)
    steps_equil = int(round(args.equil_ps * 1000 / args.dt_fs))
    steps_save = int(round(args.save_every * 1000 / args.dt_fs))
    n_snap = int(round(args.prod_ps / args.save_every))
    log = C.ROOT / "logs" / "liquid_md.log"
    (C.ROOT / "logs").mkdir(exist_ok=True)
    sim.reporters.append(app.StateDataReporter(str(log), max(steps_save // 4, 100), step=True, time=True,
                                               potentialEnergy=True, temperature=True, density=True))
    sim.step(steps_equil)
    state = sim.context.getState(getEnergy=True)
    print(f"# equilibrated {args.equil_ps} ps in {time.time() - t0:.0f} s, "
          f"E_pot={state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / n_water:.2f} kJ/mol per water", flush=True)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    frames, boxes = [], []
    for i in range(n_snap):
        sim.step(steps_save)
        st = sim.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
        pos = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)[real_idx]
        L = float(st.getPeriodicBoxVectors()[0][0].value_in_unit(unit.angstrom))
        frames.append(pos.copy())
        boxes.append(L)
        write_extxyz(outdir / f"snap_{i:02d}.xyz", symbols, pos, L,
                     f"model={args.model} T={args.temperature} t_ps={(i + 1) * args.save_every:.1f}")
        print(f"# snapshot {i:02d} at {(i + 1) * args.save_every:.1f} ps, L={L:.3f} A, "
              f"density={n_water * 18.01528 / 0.602214076 / L ** 3:.3f}, "
              f"E_pot={st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / n_water:.2f} kJ/mol/water, "
              f"{time.time() - t0:.0f} s", flush=True)

    L_mean = float(np.mean(boxes))
    r, g = oo_rdf(frames, boxes)
    sel = (r > 2.4) & (r < 3.5)
    i_peak = int(np.flatnonzero(sel)[np.argmax(g[sel])])
    rho = n_water / L_mean ** 3
    mask = r <= 3.3
    n_coord = float(np.sum(4 * np.pi * rho * g[mask] * r[mask] ** 2 * (r[1] - r[0])))
    summary = {"nwater": n_water, "L_A": L_mean, "L_per_snapshot": boxes, "initial_box_A": args.box, "model": args.model, "temperature_K": args.temperature,
               "dt_fs": args.dt_fs, "equil_ps": args.equil_ps, "prod_ps": args.prod_ps, "save_every_ps": args.save_every,
               "n_snapshots": n_snap, "goo_first_peak_A": float(r[i_peak]), "goo_first_peak_height": float(g[i_peak]),
               "coordination_below_3.3A": n_coord, "rdf_r": r.tolist(), "rdf_g": g.tolist(),
               "wall_seconds": time.time() - t0}
    C.RESULT_DIR.mkdir(exist_ok=True)
    C.dump_json(C.RESULT_DIR / "liquid_md.json", summary)
    print(f"# g_OO first peak at {r[i_peak]:.2f} A (height {g[i_peak]:.2f}), n_OO(<3.3 A)={n_coord:.2f}")
    print(f"# wrote {n_snap} snapshots to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
