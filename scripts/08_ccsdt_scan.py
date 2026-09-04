"""Counterpoise-corrected HF / MP2 / CCSD / CCSD(T) rigid O...O scan of the water
dimer: the total-curvature benchmark for the DFT force constants (WP3).

For every R three calculations in the full dimer basis (dimer, donor with
ghost acceptor, acceptor with ghost donor) give E_int at each level.  Points
are appended to results/ccsdt_<tag>.json as they finish (restartable).  A
polynomial fit around the minimum then gives R_min, k_total and the harmonic
wavenumber at each level, written to results/ccsdt_<tag>.md.

Usage:
    python scripts/08_ccsdt_scan.py --geom geometries/dimer_revpbe0-d3bj_tzvp.xyz \
        --basis aug-cc-pvtz --tag ccsdt_avtz --rmin 2.70 --rmax 3.30 --step 0.05 \
        --fine-halfwidth 0.16 --fine-step 0.02 --restart
    python scripts/08_ccsdt_scan.py --analyze results/ccsdt_ccsdt_avtz.json
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

LEVELS = ("HF", "MP2", "CCSD", "CCSD(T)")


def energies(symbols, coords, real, basis, max_memory):
    """Return dict of total energies (Eh) with atoms where real[i] is False as ghosts."""
    from pyscf import cc, gto, mp, scf

    atom = [(s if r else f"ghost-{s}", tuple(float(x) for x in c)) for s, c, r in zip(symbols, coords, real)]
    mol = gto.M(atom=atom, basis=basis, unit="Angstrom", verbose=0, max_memory=max_memory)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10
    e_hf = mf.kernel()
    if not mf.converged:
        raise RuntimeError("HF did not converge")
    n_core = sum(1 for s, r in zip(symbols, real) if r and s == "O")  # freeze O 1s
    pt2 = mp.MP2(mf).set(frozen=n_core).run()
    mycc = cc.CCSD(mf).set(frozen=n_core)
    mycc.conv_tol = 1e-8
    mycc.conv_tol_normt = 1e-6
    mycc.kernel()
    if not mycc.converged:
        raise RuntimeError("CCSD did not converge")
    e_t = mycc.ccsd_t()
    return {"HF": float(e_hf), "MP2": float(pt2.e_tot), "CCSD": float(mycc.e_tot), "CCSD(T)": float(mycc.e_tot + e_t)}


def analyze(path, window, degree):
    data = C.load_json(path)
    pts = sorted(data["points"], key=lambda p: p["R_OO"])
    r = np.array([p["R_OO"] for p in pts])
    lines = [f"# CP-corrected total-curvature benchmark: {data['tag']} ({data['basis']}, frozen core)", "",
             f"{len(pts)} points, fit window +-{window} A, degree {degree}.  omega with mu = m(H2O)/2.", "",
             "| level | R_min A | E_int kcal/mol | k kcal/mol/A^2 | k N/m | omega(H2O) cm-1 | omega(D2O) | cubic kcal/mol/A^3 |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    summary = {}
    for lev in LEVELS:
        y = np.array([p["E_int_kcal"][lev] for p in pts])
        i0 = int(np.argmin(y))
        mask = np.abs(r - r[i0]) <= window + 1e-9
        if mask.sum() < degree + 2:
            lines.append(f"| {lev} | too few points in window ({mask.sum()}) |")
            continue
        poly = np.polynomial.polynomial.Polynomial.fit(r[mask], y[mask], degree)
        roots = poly.deriv(1).roots()
        real = roots[np.isreal(roots)].real
        real = real[np.abs(real - r[i0]) <= window]
        if len(real) == 0:
            lines.append(f"| {lev} | no stationary point in window |")
            continue
        r_min = float(real[np.argmin(np.abs(real - r[i0]))])
        poly = np.polynomial.polynomial.Polynomial.fit(r[np.abs(r - r_min) <= window + 1e-9],
                                                       y[np.abs(r - r_min) <= window + 1e-9], degree)
        k = float(poly.deriv(2)(r_min))
        k_h = k / C.HARTREE_KCAL
        summary[lev] = {"R_min_A": r_min, "E_int_min_kcal": float(poly(r_min)), "k_kcal_per_A2": k,
                        "k_N_per_m": C.k_to_newton_per_metre(k_h),
                        "omega_H2O_cm-1": C.k_to_wavenumber(k_h, C.MU_H2O),
                        "omega_D2O_cm-1": C.k_to_wavenumber(k_h, C.MU_D2O),
                        "cubic_kcal_per_A3": float(poly.deriv(3)(r_min)), "points_in_window": int(mask.sum())}
        s = summary[lev]
        lines.append(f"| {lev} | {r_min:.4f} | {s['E_int_min_kcal']:.3f} | {k:.3f} | {s['k_N_per_m']:.2f} | "
                     f"{s['omega_H2O_cm-1']:.1f} | {s['omega_D2O_cm-1']:.1f} | {s['cubic_kcal_per_A3']:.2f} |")
    lines += ["", "| R A | " + " | ".join(f"E_int {lev}" for lev in LEVELS) + " |",
              "|---:|" + "---:|" * len(LEVELS)]
    for p in pts:
        lines.append(f"| {p['R_OO']:.3f} | " + " | ".join(f"{p['E_int_kcal'][lev]:.4f}" for lev in LEVELS) + " |")
    text = "\n".join(lines) + "\n"
    data["fit"] = {"window_A": window, "degree": degree, "summary": summary}
    C.dump_json(path, data)
    md = Path(path).with_suffix(".md")
    md.write_text(text)
    print(text)
    print(f"# wrote {md}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", default=None)
    parser.add_argument("--basis", default="aug-cc-pvtz")
    parser.add_argument("--tag", default=None)
    parser.add_argument("--rmin", type=float, default=2.70)
    parser.add_argument("--rmax", type=float, default=3.30)
    parser.add_argument("--step", type=float, default=0.05)
    parser.add_argument("--center", type=float, default=None)
    parser.add_argument("--fine-halfwidth", type=float, default=0.16)
    parser.add_argument("--fine-step", type=float, default=0.02)
    parser.add_argument("--max-memory", type=float, default=60000.0, help="MB for PySCF")
    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--analyze", default=None, help="only analyze this JSON")
    parser.add_argument("--window", type=float, default=0.16)
    parser.add_argument("--degree", type=int, default=4)
    args = parser.parse_args()
    if args.analyze:
        analyze(args.analyze, args.window, args.degree)
        return 0

    tag = args.tag or f"ccsdt_{args.basis.replace('aug-cc-p', 'a').replace('cc-p', '').lower()}"
    symbols, coords0 = C.read_xyz(args.geom)
    center = args.center if args.center is not None else C.oo_distance(coords0)
    coarse = np.arange(args.rmin, args.rmax + 1e-9, args.step)
    fine = np.arange(center - args.fine_halfwidth, center + args.fine_halfwidth + 1e-9, args.fine_step)
    grid = np.unique(np.round(np.concatenate([coarse, fine]), 4))

    C.RESULT_DIR.mkdir(exist_ok=True)
    out = C.RESULT_DIR / f"ccsdt_{tag}.json"
    payload = {"tag": tag, "geometry": str(args.geom), "basis": args.basis, "frozen_core": True,
               "counterpoise": True, "reference_R_OO": center, "points": []}
    if args.restart and out.exists():
        payload = C.load_json(out)
    done = {round(p["R_OO"], 4) for p in payload["points"]}
    real_d = [i in C.DONOR for i in range(len(symbols))]
    real_a = [i in C.ACCEPTOR for i in range(len(symbols))]
    print(f"# {tag}: {len(grid)} points, {len(done)} done", flush=True)
    for r in grid:
        if round(float(r), 4) in done:
            continue
        t0 = time.time()
        coords = C.set_oo_distance(coords0, float(r))
        e_ab = energies(symbols, coords, [True] * len(symbols), args.basis, args.max_memory)
        e_a = energies(symbols, coords, real_d, args.basis, args.max_memory)
        e_b = energies(symbols, coords, real_a, args.basis, args.max_memory)
        row = {"R_OO": float(r), "E_dimer": e_ab, "E_donor_cp": e_a, "E_acceptor_cp": e_b,
               "E_int_kcal": {lev: (e_ab[lev] - e_a[lev] - e_b[lev]) * C.HARTREE_KCAL for lev in LEVELS},
               "wall_seconds": time.time() - t0}
        payload["points"].append(row)
        payload["points"].sort(key=lambda p: p["R_OO"])
        C.dump_json(out, payload)
        print(f"R={r:6.3f}  " + "  ".join(f"{lev}={row['E_int_kcal'][lev]:8.4f}" for lev in LEVELS)
              + f"  {row['wall_seconds']:.0f}s", flush=True)
    analyze(out, args.window, args.degree)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
