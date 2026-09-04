"""OpenMM MD of one H3O+ (with Cl-) or one OH- (with Na+) in TIP4P-Ew water
to generate liquid acid/base snapshots for WP6.

The ions are simple non-polarizable rigid models (charges: H3O+ O -0.32 / H
+0.44, OH- O -1.32 / H +0.32, LJ of the TIP4P-Ew oxygen; counter-ions from
Joung-Cheatham TIP4P-Ew) that only serve to generate solvation structures;
the electronic structure is done later by DM-EDA on cut-out clusters.
Snapshots are written as extended xyz (real atoms only, unwrapped) plus a
JSON listing the residues and their atom indices.

Usage:
    python scripts/10b_ion_md.py --ion h3o --box 15.5 --equil-ps 60 --prod-ps 40 --save-every 4 --seed 7
    python scripts/10b_ion_md.py --ion oh  ...
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

ION_XML = """<ForceField>
 <AtomTypes>
  <Type name="h3o-O" class="OH3" element="O" mass="15.99943"/>
  <Type name="h3o-H" class="HO3" element="H" mass="1.007947"/>
  <Type name="ohm-O" class="OHM" element="O" mass="15.99943"/>
  <Type name="ohm-H" class="HOM" element="H" mass="1.007947"/>
  <Type name="na-ion" class="Na" element="Na" mass="22.98977"/>
  <Type name="cl-ion" class="Cl" element="Cl" mass="35.453"/>
 </AtomTypes>
 <Residues>
  <Residue name="H3O">
   <Atom name="O" type="h3o-O" charge="-0.32"/>
   <Atom name="H1" type="h3o-H" charge="0.44"/>
   <Atom name="H2" type="h3o-H" charge="0.44"/>
   <Atom name="H3" type="h3o-H" charge="0.44"/>
   <Bond atomName1="O" atomName2="H1"/>
   <Bond atomName1="O" atomName2="H2"/>
   <Bond atomName1="O" atomName2="H3"/>
  </Residue>
  <Residue name="OHM">
   <Atom name="O" type="ohm-O" charge="-1.32"/>
   <Atom name="H" type="ohm-H" charge="0.32"/>
   <Bond atomName1="O" atomName2="H"/>
  </Residue>
  <Residue name="NA"><Atom name="Na" type="na-ion" charge="1"/></Residue>
  <Residue name="CL"><Atom name="Cl" type="cl-ion" charge="-1"/></Residue>
 </Residues>
 <HarmonicBondForce>
  <Bond class1="OH3" class2="HO3" length="0.098" k="462750"/>
  <Bond class1="OHM" class2="HOM" length="0.096" k="462750"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle class1="HO3" class2="OH3" class3="HO3" angle="1.948" k="836"/>
 </HarmonicAngleForce>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="h3o-O" charge="-0.32" sigma="0.316435" epsilon="0.680946"/>
  <Atom type="h3o-H" charge="0.44" sigma="1" epsilon="0"/>
  <Atom type="ohm-O" charge="-1.32" sigma="0.316435" epsilon="0.680946"/>
  <Atom type="ohm-H" charge="0.32" sigma="1" epsilon="0"/>
  <Atom type="na-ion" charge="1" sigma="0.2184" epsilon="0.7047"/>
  <Atom type="cl-ion" charge="-1" sigma="0.4918" epsilon="0.04895"/>
 </NonbondedForce>
</ForceField>
"""


def ion_topology(ion: str, centre_nm: float):
    from openmm import Vec3, app, unit

    top = app.Topology()
    chain = top.addChain()
    c = centre_nm
    if ion == "h3o":
        res = top.addResidue("H3O", chain)
        o = top.addAtom("O", app.element.oxygen, res)
        hs = [top.addAtom(f"H{i}", app.element.hydrogen, res) for i in (1, 2, 3)]
        for h in hs:
            top.addBond(o, h)
        r, dz = 0.098, 0.03
        rho = np.sqrt(r ** 2 - dz ** 2)
        pos = [Vec3(c, c, c)] + [Vec3(c + rho * np.cos(a), c + rho * np.sin(a), c - dz)
                                 for a in (0.0, 2 * np.pi / 3, 4 * np.pi / 3)]
    elif ion == "oh":
        res = top.addResidue("OHM", chain)
        o = top.addAtom("O", app.element.oxygen, res)
        h = top.addAtom("H", app.element.hydrogen, res)
        top.addBond(o, h)
        pos = [Vec3(c, c, c), Vec3(c, c, c + 0.096)]
    else:
        raise ValueError(ion)
    return top, pos * unit.nanometer


def write_extxyz(path, symbols, pos, L, comment):
    lines = [str(len(symbols)),
             f'Lattice="{L:.6f} 0.0 0.0 0.0 {L:.6f} 0.0 0.0 0.0 {L:.6f}" Properties=species:S:1:pos:R:3 pbc="T T T" {comment}']
    for s, (x, y, z) in zip(symbols, pos):
        lines.append(f"{s} {x:.6f} {y:.6f} {z:.6f}")
    Path(path).write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ion", required=True, choices=["h3o", "oh"])
    parser.add_argument("--box", type=float, default=15.5)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--dt-fs", type=float, default=2.0)
    parser.add_argument("--equil-ps", type=float, default=60.0)
    parser.add_argument("--prod-ps", type=float, default=40.0)
    parser.add_argument("--save-every", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--no-counterion", action="store_true")
    args = parser.parse_args()

    import openmm
    from openmm import app, unit

    xml_path = C.ROOT / "scripts" / "ions_tip4pew.xml"
    xml_path.write_text(ION_XML)
    ff = app.ForceField("tip4pew.xml", str(xml_path))
    L_nm = args.box / 10.0
    top, pos = ion_topology(args.ion, L_nm / 2)
    modeller = app.Modeller(top, pos)
    modeller.addSolvent(ff, model="tip4pew", boxSize=openmm.Vec3(L_nm, L_nm, L_nm) * unit.nanometer,
                        neutralize=not args.no_counterion, positiveIon="Na+", negativeIon="Cl-")
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.7 * unit.nanometer,
                             constraints=app.HBonds, rigidWater=True)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.bar, args.temperature * unit.kelvin, 25))
    integrator = openmm.LangevinMiddleIntegrator(args.temperature * unit.kelvin, 1.0 / unit.picosecond,
                                                 args.dt_fs * unit.femtosecond)
    integrator.setRandomNumberSeed(args.seed)
    sim = app.Simulation(modeller.topology, system, integrator, openmm.Platform.getPlatformByName("CPU"),
                         {"Threads": str(args.threads)})
    sim.context.setPositions(modeller.positions)

    residues, real_idx, symbols = [], [], []
    for res in modeller.topology.residues():
        idx = []
        for a in res.atoms():
            if a.element is None:
                continue
            idx.append(len(real_idx))
            real_idx.append(a.index)
            symbols.append(a.element.symbol)
        residues.append({"name": res.name, "atoms": idx})
    real_idx = np.array(real_idx)
    names = [r["name"] for r in residues]
    n_water = names.count("HOH")
    print(f"# {args.ion}: {n_water} waters, residues {sorted(set(names))}, {len(real_idx)} real atoms, "
          f"{'no counter-ion' if args.no_counterion else 'neutralized'}", flush=True)

    t0 = time.time()
    sim.minimizeEnergy(maxIterations=500)
    sim.context.setVelocitiesToTemperature(args.temperature * unit.kelvin, args.seed)
    steps_equil = int(round(args.equil_ps * 1000 / args.dt_fs))
    steps_save = int(round(args.save_every * 1000 / args.dt_fs))
    n_snap = int(round(args.prod_ps / args.save_every))
    (C.ROOT / "logs").mkdir(exist_ok=True)
    sim.reporters.append(app.StateDataReporter(str(C.ROOT / "logs" / f"liquid_md_{args.ion}.log"), max(steps_save // 4, 100),
                                               step=True, time=True, potentialEnergy=True, temperature=True, density=True))
    sim.step(steps_equil)
    print(f"# equilibrated {args.equil_ps} ps in {time.time() - t0:.0f} s", flush=True)

    outdir = C.GEOM_DIR / f"liquid_{args.ion}"
    outdir.mkdir(parents=True, exist_ok=True)
    boxes = []
    for i in range(n_snap):
        sim.step(steps_save)
        st = sim.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
        pos = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)[real_idx]
        L = float(st.getPeriodicBoxVectors()[0][0].value_in_unit(unit.angstrom))
        boxes.append(L)
        write_extxyz(outdir / f"snap_{i:02d}.xyz", symbols, pos, L,
                     f"ion={args.ion} T={args.temperature} t_ps={(i + 1) * args.save_every:.1f}")
        C.dump_json(outdir / f"snap_{i:02d}.json", {"ion": args.ion, "L_A": L, "residues": residues})
        # ion first-shell report
        ion_o = pos[residues[0]["atoms"][0]]
        d = []
        for r in residues:
            if r["name"] == "HOH":
                v = pos[r["atoms"][0]] - ion_o
                v -= L * np.round(v / L)
                d.append(np.linalg.norm(v))
        d = np.sort(d)
        print(f"# snapshot {i:02d} at {(i + 1) * args.save_every:.1f} ps, L={L:.3f}, "
              f"E_pot={st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / n_water:.2f} kJ/mol/water, "
              f"ion O..O_w: {' '.join(f'{x:.2f}' for x in d[:5])}, {time.time() - t0:.0f} s", flush=True)
    C.dump_json(C.RESULT_DIR / f"liquid_md_{args.ion}.json",
                {"ion": args.ion, "nwater": n_water, "L_per_snapshot": boxes, "n_snapshots": n_snap,
                 "counterion": not args.no_counterion, "wall_seconds": time.time() - t0})
    print(f"# wrote {n_snap} snapshots to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
