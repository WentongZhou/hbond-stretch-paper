"""Rigid O...O scan of the water dimer with DM-EDA at every point.

The acceptor water is translated rigidly along the O...O axis of an optimized
dimer; monomer geometries and the relative orientation are frozen.  Every
point yields the full DM-EDA partition, Mulliken and IAO charge transfer and
the closure error.  Results are appended to results/scan_<tag>.json after each
point so a crash loses at most one point.

Usage:
    python scripts/01_scan_dimer.py --geom geometries/dimer_revpbe0-d3bj_tzvp.xyz \
        --xc revpbe0 --disp d3bj --basis def2-tzvp --tag revpbe0-d3bj_tzvp \
        --rmin 2.55 --rmax 3.60 --step 0.05 --fine-halfwidth 0.16 --fine-step 0.02
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def build_grid(rmin, rmax, step, center, fine_halfwidth, fine_step):
    coarse = np.arange(rmin, rmax + 1e-9, step)
    if fine_halfwidth > 0 and center is not None:
        fine = np.arange(center - fine_halfwidth, center + fine_halfwidth + 1e-9, fine_step)
        grid = np.concatenate([coarse, fine])
    else:
        grid = coarse
    grid = np.unique(np.round(grid, 4))
    return grid


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", required=True)
    parser.add_argument("--xc", default="revpbe0")
    parser.add_argument("--disp", default="d3bj")
    parser.add_argument("--basis", default="def2-tzvp")
    parser.add_argument("--grid-level", type=int, default=4)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--rmin", type=float, default=2.55)
    parser.add_argument("--rmax", type=float, default=3.60)
    parser.add_argument("--step", type=float, default=0.05)
    parser.add_argument("--center", type=float, default=None, help="centre of the fine window; default: R_OO of --geom")
    parser.add_argument("--fine-halfwidth", type=float, default=0.16)
    parser.add_argument("--fine-step", type=float, default=0.02)
    parser.add_argument("--oh-elongation", type=float, default=0.0,
                        help="stretch the donor bonded O-H by this many Angstrom before scanning")
    parser.add_argument("--restart", action="store_true", help="skip points already in the JSON")
    args = parser.parse_args()
    disp = None if args.disp in ("", "none", "None") else args.disp

    symbols, coords0 = C.read_xyz(args.geom)
    if args.oh_elongation:
        coords0 = C.elongate_bond(coords0, C.DONOR[0], C.DONOR[1], args.oh_elongation)
    center = args.center if args.center is not None else C.oo_distance(coords0)
    grid = build_grid(args.rmin, args.rmax, args.step, center, args.fine_halfwidth, args.fine_step)

    C.RESULT_DIR.mkdir(exist_ok=True)
    out = C.RESULT_DIR / f"scan_{args.tag}.json"
    payload = {
        "tag": args.tag,
        "geometry": str(args.geom),
        "xc": args.xc,
        "dispersion": disp,
        "basis": args.basis,
        "grid_level": args.grid_level,
        "oh_elongation_angstrom": args.oh_elongation,
        "reference_R_OO": center,
        "points": [],
    }
    if args.restart and out.exists():
        payload = C.load_json(out)
    done = {round(p["R_OO"], 4) for p in payload["points"]}

    config = C.make_config(args.xc, args.basis, disp, args.grid_level)
    print(f"# {args.tag}: {len(grid)} points, {len(done)} already done", flush=True)
    for r in grid:
        if round(float(r), 4) in done:
            continue
        t0 = time.time()
        coords = C.set_oo_distance(coords0, float(r))
        result, eda = C.run_eda(symbols, coords, config)
        row = C.eda_row(result, eda, extra={
            "R_OO": float(r),
            "wall_seconds": time.time() - t0,
            "closure_hartree": float(result.diagnostics["closure_error_hartree"]),
            "frozen_identity_hartree": float(result.diagnostics["frozen_identity_error_hartree"]),
        })
        payload["points"].append(row)
        payload["points"].sort(key=lambda p: p["R_OO"])
        C.dump_json(out, payload)
        ct = row["mulliken_ct"].get("acceptor", float("nan"))
        iao = row["iao_ct"].get("acceptor", float("nan")) if "error" not in row["iao_ct"] else float("nan")
        print(f"R={r:6.3f}  Total={row['Total']:9.4f}  Elec={row['Elec']:9.4f}  ExRep={row['ExRep']:9.4f}  "
              f"OrbRel={row['OrbRel']:9.4f}  CorrDisp={row['CorrDisp']:9.4f}  "
              f"CT(Mull)={ct:+.4f}  CT(IAO)={iao:+.4f}  closure={row['closure_hartree']:.1e}  "
              f"{row['wall_seconds']:.1f}s", flush=True)
    print(f"# wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
