"""Force-constant partition of the H-bond stretch from a DM-EDA scan.

For every energy component X(R) the script fits a polynomial in a window
around the minimum of the total interaction energy and reports

    k_X = d^2 X / dR^2 |_{R_min},      f_X = k_X / k_total,
    g_X = dX / dR |_{R_min}            (the forces balance: sum g_X = 0),

the harmonic wavenumber of the intermolecular stretch for H2O and D2O
reduced masses, the cubic coefficient (anharmonicity), and the charge
transfer and its slope at R_min.

Usage:
    python scripts/02_force_constants.py --scan results/scan_revpbe0-d3bj_tzvp.json
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

KEYS = ("Total", "Elec", "Exch", "Rep", "ExRep", "OrbRel", "Corr", "Disp", "CorrDisp", "Steric")


def series(points, key):
    return np.array([p[key] for p in points], dtype=float)


def ct_series(points, kind, label="acceptor"):
    out = []
    for p in points:
        block = p.get(kind, {})
        out.append(block.get(label, np.nan) if "error" not in block else np.nan)
    return np.array(out, dtype=float)


def fit_window(r, y, center, halfwidth, degree):
    mask = np.abs(r - center) <= halfwidth + 1e-9
    if mask.sum() < degree + 2:
        # near a grid edge: take the degree+3 points nearest the centre instead (asymmetric window)
        order = np.argsort(np.abs(r - center))[: degree + 3]
        mask = np.zeros(len(r), dtype=bool)
        mask[order] = True
        if mask.sum() < degree + 2:
            raise ValueError(f"only {mask.sum()} points available; widen the scan or lower --degree")
    return np.polynomial.polynomial.Polynomial.fit(r[mask], y[mask], degree), int(mask.sum())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan", required=True)
    parser.add_argument("--window", type=float, default=0.16, help="half-width (Angstrom) of the fit window")
    parser.add_argument("--degree", type=int, default=4)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    data = C.load_json(args.scan)
    tag = data["tag"]
    ct_label = data.get("ct_label", "acceptor")
    pts = sorted(data["points"], key=lambda p: p["R_OO"])
    r = series(pts, "R_OO")
    total = series(pts, "Total")

    # Locate the minimum of the total interaction energy.
    i0 = int(np.argmin(total))
    poly_total, npts = fit_window(r, total, r[i0], args.window, args.degree)
    d1 = poly_total.deriv(1)
    roots = d1.roots()
    real = roots[np.isreal(roots)].real
    real = real[np.abs(real - r[i0]) <= args.window]
    r0 = data.get("reference_R_OO")
    if len(real) == 0:
        if r0 is None:
            raise RuntimeError("no stationary point inside the window")
        print(f"# WARNING: no stationary point in the scanned window; reporting the curvature at R0={r0:.4f} only")
        r_min = None
    else:
        r_min = float(real[np.argmin(np.abs(real - r[i0]))])
    edge_note = ""
    if r_min is None:
        edge_note = " (no minimum in the scanned range)"
    elif r_min < r.min() or r_min > r.max():
        edge_note = f" (WARNING: R_min {r_min:.3f} outside the scanned range {r.min():.3f}-{r.max():.3f}; extrapolated)"
    elif min(r_min - r.min(), r.max() - r_min) < args.window:
        edge_note = f" (note: R_min within {args.window} A of the scan edge; asymmetric fit window)"

    def component_table(at):
        tab = {}
        for key in KEYS:
            y = series(pts, key)
            poly, _ = fit_window(r, y, at, args.window, args.degree)
            tab[key] = {
                "value_kcal": float(poly(at)),
                "slope_kcal_per_A": float(poly.deriv(1)(at)),
                "k_kcal_per_A2": float(poly.deriv(2)(at)),
                "cubic_kcal_per_A3": float(poly.deriv(3)(at)),
            }
        kt = tab["Total"]["k_kcal_per_A2"]
        for key in KEYS:
            tab[key]["fraction_of_k"] = tab[key]["k_kcal_per_A2"] / kt if kt else float("nan")
        return tab

    def ct_table(at):
        out = {}
        for name, y in (("mulliken", ct_m), ("iao", ct_i)):
            if np.all(np.isfinite(y)):
                poly, _ = fit_window(r, y, at, args.window, min(args.degree, 3))
                out[name] = {"value_e": float(poly(at)), "slope_e_per_A": float(poly.deriv(1)(at))}
            else:
                out[name] = None
        return out

    ct_m = ct_series(pts, "mulliken_ct", ct_label)
    ct_i = ct_series(pts, "iao_ct", ct_label)

    # instantaneous curvature at the reference geometry (always available)
    at_r0 = None
    if r0 is not None and r.min() - 1e-6 <= r0 <= r.max() + 1e-6:
        t0 = component_table(float(r0))
        k0 = t0["Total"]["k_kcal_per_A2"]
        at_r0 = {"R0_A": float(r0), "components": t0, "k_total_kcal_per_A2": k0,
                 "omega_H2O_cm-1": C.k_to_wavenumber(k0 / C.HARTREE_KCAL, C.MU_H2O),
                 "residual_force_kcal_per_A": t0["Total"]["slope_kcal_per_A"],
                 "charge_transfer_acceptor": ct_table(float(r0))}

    if r_min is None:
        # no minimum: use the R0 block as the reported table so downstream tools still work
        r_min = float(r0)
    table = component_table(r_min)
    k_total = table["Total"]["k_kcal_per_A2"]
    k_sum_primary = sum(table[k]["k_kcal_per_A2"] for k in C.PRIMARY_SHORT)
    slope_sum = sum(table[k]["slope_kcal_per_A"] for k in C.PRIMARY_SHORT)

    k_hartree_A2 = k_total / C.HARTREE_KCAL
    omega_h2o = C.k_to_wavenumber(k_hartree_A2, C.MU_H2O)
    omega_d2o = C.k_to_wavenumber(k_hartree_A2, C.MU_D2O)
    ct_fit = ct_table(r_min)

    closure_max = max(abs(p.get("closure_hartree", 0.0)) for p in pts)
    summary = {
        "tag": tag,
        "xc": data.get("xc"),
        "dispersion": data.get("dispersion"),
        "basis": data.get("basis"),
        "oh_elongation_angstrom": data.get("oh_elongation_angstrom", 0.0),
        "n_points": len(pts),
        "fit_window_halfwidth_A": args.window,
        "fit_degree": args.degree,
        "points_in_window": npts,
        "R_min_A": r_min,
        "minimum_found": "no minimum" not in edge_note,
        "edge_note": edge_note,
        "at_R0": at_r0,
        "E_int_min_kcal": table["Total"]["value_kcal"],
        "k_total_kcal_per_A2": k_total,
        "k_total_N_per_m": C.k_to_newton_per_metre(k_hartree_A2),
        "k_sum_primary_minus_total": k_sum_primary - k_total,
        "force_balance_sum_of_slopes": slope_sum,
        "omega_H2O_cm-1": omega_h2o,
        "omega_D2O_cm-1": omega_d2o,
        "omega_ratio_D2O_over_H2O_classical": omega_d2o / omega_h2o,
        "ct_label": ct_label,
        "charge_transfer_acceptor": ct_fit,
        "max_abs_closure_hartree": closure_max,
        "components": table,
    }
    out_json = C.RESULT_DIR / f"fc_{tag}.json"
    C.dump_json(out_json, summary)

    # Markdown table
    lines = [
        f"# Force-constant partition: {tag}",
        "",
        f"R_min = {r_min:.4f} Å{edge_note}, E_int = {table['Total']['value_kcal']:.3f} kcal/mol, "
        f"k_total = {k_total:.3f} kcal/mol/Å² = {summary['k_total_N_per_m']:.2f} N/m, "
        f"ω(H2O) = {omega_h2o:.1f} cm⁻¹, ω(D2O) = {omega_d2o:.1f} cm⁻¹ (mass effect only)",
        "",
        "| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key in KEYS:
        t = table[key]
        lines.append(f"| {key} | {t['value_kcal']:.3f} | {t['slope_kcal_per_A']:.3f} | {t['k_kcal_per_A2']:.3f} | "
                     f"{t['fraction_of_k']:+.3f} | {t['cubic_kcal_per_A3']:.2f} |")
    lines += [
        "",
        f"Sum of primary slopes at R_min (should be ~0): {slope_sum:.4f} kcal/mol/Å",
        f"Sum of primary k minus k_total (closure): {k_sum_primary - k_total:.2e} kcal/mol/Å²",
        f"Electrons lost by fragment '{ct_label}' (acceptor of the target bond) at R_min: Mulliken {ct_fit['mulliken']}, IAO {ct_fit['iao']}",
        f"Max |closure| over scan: {closure_max:.2e} Eh",
    ]
    if at_r0 is not None:
        c0 = at_r0["components"]
        lines += ["", f"At the reference geometry R0 = {at_r0['R0_A']:.4f} Å (residual force {at_r0['residual_force_kcal_per_A']:+.3f} kcal/mol/Å): "
                  f"k = {at_r0['k_total_kcal_per_A2']:.3f}, ω = {at_r0['omega_H2O_cm-1']:.1f} cm⁻¹, "
                  + ", ".join(f"f_{k} = {c0[k]['fraction_of_k']:+.3f}" for k in ("Elec", "ExRep", "OrbRel", "CorrDisp"))]
    (C.RESULT_DIR / f"fc_{tag}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    if not args.no_plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        ax = axes[0]
        for key, color in zip(("Total", "Elec", "ExRep", "OrbRel", "CorrDisp"),
                              ("k", "tab:blue", "tab:red", "tab:green", "tab:purple")):
            ax.plot(r, series(pts, key), "o-", ms=3, color=color, label=key)
        ax.axvline(r_min, color="gray", ls="--", lw=0.8)
        ax.set_xlabel("R(O···O) / Å")
        ax.set_ylabel("kcal/mol")
        ax.set_title(f"DM-EDA components, {tag}")
        ax.legend(fontsize=8)

        ax = axes[1]
        names = list(C.PRIMARY_SHORT) + ["Total"]
        vals = [table[k]["k_kcal_per_A2"] for k in names]
        ax.bar(names, vals, color=["tab:blue", "tab:orange", "tab:red", "tab:green", "tab:purple", "tab:brown", "k"])
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylabel("k_X = d²X/dR² / kcal mol⁻¹ Å⁻²")
        ax.set_title(f"Force-constant partition at R_min={r_min:.3f} Å")
        for i, v in enumerate(vals):
            ax.text(i, v, f"{v / k_total:+.2f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8)

        ax = axes[2]
        ax.plot(r, ct_m, "o-", ms=3, label="Mulliken")
        if np.all(np.isfinite(ct_i)):
            ax.plot(r, ct_i, "s-", ms=3, label="IAO")
        ax.axvline(r_min, color="gray", ls="--", lw=0.8)
        ax.set_xlabel("R(O···O) / Å")
        ax.set_ylabel(f"electrons lost by '{ct_label}' / e")
        ax.set_title("CT acceptor → donor of the target bond")
        ax.legend(fontsize=8)
        fig.tight_layout()
        C.FIG_DIR.mkdir(exist_ok=True)
        fig.savefig(C.FIG_DIR / f"fc_{tag}.png", dpi=160)
        print(f"# figure: {C.FIG_DIR / f'fc_{tag}.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
