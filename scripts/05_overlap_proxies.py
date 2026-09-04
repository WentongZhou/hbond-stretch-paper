"""Do charge transfer, exchange repulsion and the saddle-point density all decay
with the same exponential length along the O...O coordinate?

If yes, any of them can be used as an empirical proxy for "how far up the
repulsive wall the H-bond sits", which is what the stretch frequency measures.
Fits ln(y) = a - b R over a window and reports b (1/Angstrom) and the
effective decay length 1/b for every candidate.

Usage:
    python scripts/05_overlap_proxies.py --scan results/scan_revpbe0-d3bj_tzvp.json \
        --density results/density_revpbe0-d3bj_tzvp.json --rmin 2.7 --rmax 3.3
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def logfit(r, y):
    mask = np.isfinite(y) & (y > 0)
    coef = np.polyfit(r[mask], np.log(y[mask]), 1)
    pred = np.polyval(coef, r[mask])
    ss_res = np.sum((np.log(y[mask]) - pred) ** 2)
    ss_tot = np.sum((np.log(y[mask]) - np.log(y[mask]).mean()) ** 2)
    return -coef[0], 1 - ss_res / ss_tot, int(mask.sum())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan", required=True)
    parser.add_argument("--density", default=None)
    parser.add_argument("--rmin", type=float, default=2.7)
    parser.add_argument("--rmax", type=float, default=3.3)
    args = parser.parse_args()

    data = C.load_json(args.scan)
    pts = [p for p in sorted(data["points"], key=lambda p: p["R_OO"]) if args.rmin - 1e-9 <= p["R_OO"] <= args.rmax + 1e-9]
    r = np.array([p["R_OO"] for p in pts])
    cands = {
        "ExRep (exchange + Pauli repulsion)": np.array([p["ExRep"] for p in pts]),
        "Rep (Pauli repulsion only)": np.array([p["Rep"] for p in pts]),
        "-Exch (interfragment exchange)": np.array([-p["Exch"] for p in pts]),
        "-OrbRel (orbital relaxation)": np.array([-p["OrbRel"] for p in pts]),
        "-Elec (electrostatics)": np.array([-p["Elec"] for p in pts]),
        "-CorrDisp": np.array([-p["CorrDisp"] for p in pts]),
        "CT Mulliken (acceptor -> donor)": np.array([p["mulliken_ct"]["acceptor"] for p in pts]),
        "CT IAO (acceptor -> donor)": np.array([p["iao_ct"].get("acceptor", np.nan) for p in pts]),
    }
    rows = []
    for name, y in cands.items():
        b, r2, n = logfit(r, y)
        rows.append((name, b, 1 / b, r2, n))
    if args.density:
        dd = C.load_json(args.density)
        recs = [x for x in dd["records"] if args.rmin - 1e-9 <= x["R_OO"] <= args.rmax + 1e-9]
        rr = np.array([x["R_OO"] for x in recs])
        for name, key in (("rho_final at saddle", "rho_final"), ("rho_0 at saddle (promolecule)", "rho_0"),
                          ("d rho_relax at saddle", "d_rho_relax"), ("-d rho_Pauli at saddle", "d_rho_pauli")):
            y = np.array([x["saddle"][key] for x in recs])
            if key == "d_rho_pauli":
                y = -y
            b, r2, n = logfit(rr, y)
            rows.append((name, b, 1 / b, r2, n))

    lines = [f"# Exponential decay constants along R(O···O), window {args.rmin}-{args.rmax} Å, {data['tag']}", "",
             "ln y = a − b·R.  Similar b means the quantities are interchangeable proxies for the position on the overlap wall.", "",
             "| quantity | b / Å⁻¹ | decay length 1/b / Å | R² of log-linear fit | n |", "|---|---:|---:|---:|---:|"]
    for name, b, L, r2, n in rows:
        lines.append(f"| {name} | {b:.3f} | {L:.3f} | {r2:.4f} | {n} |")
    text = "\n".join(lines) + "\n"
    out = C.RESULT_DIR / f"proxies_{data['tag']}.md"
    out.write_text(text)
    print(text)
    print(f"# wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
