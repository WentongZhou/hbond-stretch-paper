"""Rigid scan of one target H-bond inside a cluster with DM-EDA at every point.

Generalization of 01_scan_dimer.py: the fragment definition, the moving
fragment and the target O...O pair come from the cluster JSON written by
06_build_clusters.py.  The moving fragment is translated rigidly along the
axis joining the two target oxygens; everything else is frozen.  The output
JSON has the same layout as scan_<tag>.json from 01 so 02_force_constants.py
and 05_overlap_proxies.py work unchanged (they read ``ct_label``).

Usage:
    python scripts/07_scan_cluster.py --cluster results/cluster_acid_revpbe0-d3bj_tzvp.json \
        --set acid --xc revpbe0 --disp d3bj --basis def2-tzvp --tag acid_revpbe0-d3bj_tzvp \
        --rmin -0.30 --rmax 0.50 --step 0.05 --fine-halfwidth 0.16 --fine-step 0.02 --restart
``--rmin/--rmax`` are offsets from the optimized target distance (Angstrom).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def build_grid(center, rmin_off, rmax_off, step, fine_halfwidth, fine_step):
    coarse = center + np.arange(rmin_off, rmax_off + 1e-9, step)
    fine = center + np.arange(-fine_halfwidth, fine_halfwidth + 1e-9, fine_step) if fine_halfwidth > 0 else []
    return np.unique(np.round(np.concatenate([coarse, fine]), 4))


def shift_fragment(coords, fixed_o, moving_o, moving_atoms, r_target):
    axis = coords[moving_o] - coords[fixed_o]
    r_now = float(np.linalg.norm(axis))
    axis /= r_now
    new = coords.copy()
    for i in moving_atoms:
        new[i] += (r_target - r_now) * axis
    return new


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster", required=True, help="results/cluster_<system>_<tag>.json from 06")
    parser.add_argument("--set", required=True, help="fragment set name inside the cluster JSON")
    parser.add_argument("--geom", default=None, help="override the xyz recorded in the cluster JSON")
    parser.add_argument("--xc", default="revpbe0")
    parser.add_argument("--disp", default="d3bj")
    parser.add_argument("--basis", default="def2-tzvp")
    parser.add_argument("--grid-level", type=int, default=4)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--rmin", type=float, default=-0.30, help="offset from R0 (Angstrom)")
    parser.add_argument("--rmax", type=float, default=0.50)
    parser.add_argument("--step", type=float, default=0.05)
    parser.add_argument("--fine-halfwidth", type=float, default=0.16)
    parser.add_argument("--fine-step", type=float, default=0.02)
    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--density-fit", action="store_true", help="PySCF density fitting in every SCF (labelled approximation)")
    parser.add_argument("--embed", action="store_true",
                        help="pair scan with the rest of the environment as point charges (spec['embedding'] from 13_add_embedding.py)")
    parser.add_argument("--pair", action="store_true",
                        help="keep only the two water molecules of the target bond (frozen cluster geometry): "
                             "isolates the geometric distortion from the electronic effect of the rest of the cluster")
    args = parser.parse_args()
    disp = None if args.disp in ("", "none", "None") else args.disp

    from pyscf_dm_eda import FragmentSpec

    cl = C.load_json(args.cluster)
    fs = cl["fragment_sets"][args.set]
    geom = args.geom or cl["xyz"]
    if not Path(geom).is_absolute():
        geom = C.ROOT / geom
    symbols, coords0 = C.read_xyz(geom)
    fs = dict(fs)
    point_charges = None
    if args.embed:
        if "embedding" not in cl:
            raise SystemExit("cluster spec has no 'embedding'; run 13_add_embedding.py first")
        args.pair = True
        point_charges = (np.asarray(cl["embedding"]["coords"], dtype=float), np.asarray(cl["embedding"]["charges"], dtype=float))
    if args.pair:
        fs = C.pair_fragment_set(fs, symbols, coords0)
        keep = sorted(set(fs["fragments"][0]["atoms"]) | set(fs["fragments"][1]["atoms"]))
        remap = {old: new for new, old in enumerate(keep)}
        symbols = [symbols[i] for i in keep]
        coords0 = coords0[keep]
        for f in fs["fragments"]:
            f["atoms"] = [remap[i] for i in f["atoms"]]
        for key in ("donor_O", "bonded_H", "acceptor_O"):
            fs[key] = remap[fs[key]]
    fragments = [FragmentSpec(tuple(f["atoms"]), int(f["charge"]), 0, f["label"]) for f in fs["fragments"]]
    moving = next(f for f in fs["fragments"] if f["label"] == fs["moving"])
    moving_atoms = list(moving["atoms"])
    od, oa = fs["donor_O"], fs["acceptor_O"]
    moving_o, fixed_o = (oa, od) if oa in moving_atoms else (od, oa)
    r0 = float(np.linalg.norm(coords0[oa] - coords0[od]))
    grid = build_grid(r0, args.rmin, args.rmax, args.step, args.fine_halfwidth, args.fine_step)

    C.RESULT_DIR.mkdir(exist_ok=True)
    out = C.RESULT_DIR / f"scan_{args.tag}.json"
    payload = {
        "tag": args.tag,
        "geometry": str(geom),
        "cluster": str(args.cluster),
        "fragment_set": args.set + ("_embed" if args.embed else "_pair" if args.pair else ""),
        "pair_only": bool(args.pair),
        "symbols": list(symbols),
        "system": cl.get("system"),
        "xc": args.xc,
        "dispersion": disp,
        "basis": args.basis,
        "grid_level": args.grid_level,
        "oh_elongation_angstrom": 0.0,
        "reference_R_OO": r0,
        "fragments": fs["fragments"],
        "moving": fs["moving"],
        "target": {"donor_O": od, "bonded_H": fs["bonded_H"], "acceptor_O": oa},
        "ct_label": fs["ct_label"],
        "points": [],
    }
    if args.restart and out.exists():
        old = C.load_json(out)
        payload["points"] = old.get("points", [])
    done = {round(p["R_OO"], 4) for p in payload["points"]}

    config = C.make_config(args.xc, args.basis, disp, args.grid_level, density_fit=args.density_fit,
                           point_charges=point_charges)
    payload["density_fit"] = bool(args.density_fit)
    payload["embedded"] = bool(args.embed)
    if args.embed:
        payload["embedding"] = {k: v for k, v in cl["embedding"].items() if k not in ("coords", "charges")}
    print(f"# {args.tag}: set={args.set} R0={r0:.4f} moving={fs['moving']} ct_label={fs['ct_label']} "
          f"{len(grid)} points, {len(done)} done", flush=True)
    for r in grid:
        if round(float(r), 4) in done:
            continue
        t0 = time.time()
        coords = shift_fragment(coords0, fixed_o, moving_o, moving_atoms, float(r))
        result, eda = C.run_eda(symbols, coords, config, fragments=fragments)
        row = C.eda_row(result, eda, extra={
            "R_OO": float(r),
            "wall_seconds": time.time() - t0,
            "closure_hartree": float(result.diagnostics["closure_error_hartree"]),
            "frozen_identity_hartree": float(result.diagnostics["frozen_identity_error_hartree"]),
        })
        payload["points"].append(row)
        payload["points"].sort(key=lambda p: p["R_OO"])
        C.dump_json(out, payload)
        lab = fs["ct_label"]
        ct = row["mulliken_ct"].get(lab, float("nan"))
        iao = row["iao_ct"].get(lab, float("nan")) if "error" not in row["iao_ct"] else float("nan")
        print(f"R={r:6.3f}  Total={row['Total']:9.4f}  Elec={row['Elec']:9.4f}  ExRep={row['ExRep']:9.4f}  "
              f"OrbRel={row['OrbRel']:9.4f}  CorrDisp={row['CorrDisp']:9.4f}  "
              f"CT(Mull)={ct:+.4f}  CT(IAO)={iao:+.4f}  closure={row['closure_hartree']:.1e}  "
              f"{row['wall_seconds']:.1f}s", flush=True)
    print(f"# wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
