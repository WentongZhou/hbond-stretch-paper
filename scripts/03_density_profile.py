"""Decompose the electron density along the O-H···O axis into DM-EDA stages.

rho_0      = rho[P_donor] + rho[P_acceptor]        (promolecule overlap)
d rho_P    = rho[P_Pauli] - rho_0                  (antisymmetrization)
d rho_relax= rho[P_S]     - rho[P_Pauli]           (polarization + charge transfer)

The saddle point is taken as the minimum of the final density along the
O_d -> O_a line between the bonded H and the acceptor O; this is the quantity
the CVS paper's k ~ dq^2 model calls "charge density at the saddle point".

Usage:
    python scripts/03_density_profile.py --geom geometries/dimer_revpbe0-d3bj_tzvp.xyz \
        --xc revpbe0 --disp d3bj --basis def2-tzvp --tag revpbe0-d3bj_tzvp --R 2.7 2.8 2.9 3.0 3.1 3.2
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", default=None)
    parser.add_argument("--cluster", default=None, help="results/cluster_<system>_<tag>.json from 06 (target bond from --set)")
    parser.add_argument("--set", default=None, help="fragment set inside the cluster JSON")
    parser.add_argument("--xc", default="revpbe0")
    parser.add_argument("--disp", default="d3bj")
    parser.add_argument("--basis", default="def2-tzvp")
    parser.add_argument("--grid-level", type=int, default=4)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--R", type=float, nargs="*", default=None)
    parser.add_argument("--npts", type=int, default=600)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    disp = None if args.disp in ("", "none", "None") else args.disp

    if args.cluster:
        cl = C.load_json(args.cluster)
        fs = cl["fragment_sets"][args.set]
        geom = args.geom or cl["xyz"]
        if not Path(geom).is_absolute():
            geom = C.ROOT / geom
        symbols, coords0 = C.read_xyz(geom)
        fragments, moving_atoms, fixed_o, moving_o, i_od, i_hb, i_oa, i_don, i_acc = C.target_from_spec(fs)
        ct_label = fs["ct_label"]
        r_ref = float(np.linalg.norm(coords0[i_oa] - coords0[i_od]))

        def place(r):
            return C.shift_fragment(coords0, fixed_o, moving_o, moving_atoms, float(r))
    else:
        geom = args.geom
        symbols, coords0 = C.read_xyz(geom)
        fragments = None
        i_od, i_hb, i_oa, i_don, i_acc = C.DONOR[0], C.DONOR[1], C.ACCEPTOR[0], 0, 1
        ct_label = "acceptor"
        r_ref = C.oo_distance(coords0)

        def place(r):
            return C.set_oo_distance(coords0, float(r))

    radii = args.R or [r_ref]
    config = C.make_config(args.xc, args.basis, disp, args.grid_level)

    records = []
    profiles = {}
    for r_oo in radii:
        coords = place(r_oo)
        result, eda = C.run_eda(symbols, coords, config, fragments=fragments)
        dms = C.total_density_matrices(eda)
        mol = eda.super_state.mol
        o_d, h_b, o_a = coords[i_od], coords[i_hb], coords[i_oa]
        t = np.linspace(0.0, 1.0, args.npts)
        line = o_d[None, :] + t[:, None] * (o_a - o_d)[None, :]
        rho_d = C.rho_on_points(mol, dms["fragments"][i_don], line)
        rho_a = C.rho_on_points(mol, dms["fragments"][i_acc], line)
        rho_0 = C.rho_on_points(mol, dms["P0"], line)
        rho_p = C.rho_on_points(mol, dms["Pauli"], line)
        rho_s = C.rho_on_points(mol, dms["S"], line)
        # projection of the bonded H on the axis
        t_h = float(np.dot(h_b - o_d, o_a - o_d) / np.dot(o_a - o_d, o_a - o_d))
        window = (t > t_h + 0.05) & (t < 0.97)
        i_saddle = int(np.flatnonzero(window)[np.argmin(rho_s[window])])
        i_mid = int(np.argmin(np.abs(t - 0.5)))
        comp = result.components("kcal/mol")

        def at(i):
            return {
                "t": float(t[i]),
                "distance_from_Od_A": float(t[i] * r_oo),
                "distance_from_Hb_A": float(np.linalg.norm(line[i] - h_b)),
                "rho_donor": float(rho_d[i]),
                "rho_acceptor": float(rho_a[i]),
                "rho_0": float(rho_0[i]),
                "rho_pauli": float(rho_p[i]),
                "rho_final": float(rho_s[i]),
                "d_rho_pauli": float(rho_p[i] - rho_0[i]),
                "d_rho_relax": float(rho_s[i] - rho_p[i]),
                "d_rho_total": float(rho_s[i] - rho_0[i]),
                "fraction_promolecule": float(rho_0[i] / rho_s[i]),
                "fraction_pauli": float((rho_p[i] - rho_0[i]) / rho_s[i]),
                "fraction_relax": float((rho_s[i] - rho_p[i]) / rho_s[i]),
            }

        rec = {
            "R_OO": float(r_oo),
            "saddle": at(i_saddle),
            "midpoint": at(i_mid),
            "components_kcal": {C.SHORT.get(k, k): float(v) for k, v in comp.items()},
            "mulliken_ct": {k: float(v) for k, v in result.fragment_charge_transfer.items()},
            "iao_ct": C.iao_fragment_charges(eda),
        }
        records.append(rec)
        profiles[f"{r_oo:.3f}"] = {
            "t": t.tolist(), "rho_0": rho_0.tolist(), "rho_pauli": rho_p.tolist(), "rho_final": rho_s.tolist(),
            "rho_donor": rho_d.tolist(), "rho_acceptor": rho_a.tolist(), "t_Hb": t_h,
        }
        s = rec["saddle"]
        print(f"R={r_oo:.3f}  saddle at {s['distance_from_Hb_A']:.3f} Å from H_b: rho_final={s['rho_final']:.5f} "
              f"rho_0={s['rho_0']:.5f} ({100 * s['fraction_promolecule']:.1f}%)  dPauli={s['d_rho_pauli']:+.5f} "
              f"({100 * s['fraction_pauli']:+.1f}%)  drelax={s['d_rho_relax']:+.5f} ({100 * s['fraction_relax']:+.1f}%)  "
              f"Total={comp['Total Interaction energy']:.3f}", flush=True)

    C.RESULT_DIR.mkdir(exist_ok=True)
    out = C.RESULT_DIR / f"density_{args.tag}.json"
    C.dump_json(out, {"tag": args.tag, "xc": args.xc, "dispersion": disp, "basis": args.basis,
                      "geometry": str(geom), "cluster": args.cluster, "fragment_set": args.set, "ct_label": ct_label,
                      "reference_R_OO": r_ref, "records": records, "profiles": profiles})
    print(f"# wrote {out}")

    if not args.no_plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        key = f"{radii[len(radii) // 2]:.3f}"
        prof = profiles[key]
        tt = np.array(prof["t"]) * float(key)
        ax = axes[0]
        ax.semilogy(tt, prof["rho_0"], label="ρ₀ (promolecule)")
        ax.semilogy(tt, prof["rho_pauli"], label="ρ_Pauli")
        ax.semilogy(tt, prof["rho_final"], label="ρ_final")
        ax.axvline(prof["t_Hb"] * float(key), color="gray", ls=":", lw=0.8, label="H_b projection")
        ax.set_xlabel("distance from O_d along O···O / Å")
        ax.set_ylabel("ρ / e bohr⁻³")
        ax.set_title(f"Density along the H-bond axis, R = {key} Å")
        ax.legend(fontsize=8)
        ax = axes[1]
        ax.plot(tt, np.array(prof["rho_pauli"]) - np.array(prof["rho_0"]), label="Δρ_Pauli")
        ax.plot(tt, np.array(prof["rho_final"]) - np.array(prof["rho_pauli"]), label="Δρ_relax")
        ax.plot(tt, np.array(prof["rho_final"]) - np.array(prof["rho_0"]), "k--", label="Δρ_total")
        ax.axhline(0, color="k", lw=0.6)
        ax.axvline(prof["t_Hb"] * float(key), color="gray", ls=":", lw=0.8)
        ax.set_xlim(prof["t_Hb"] * float(key) - 0.3, float(key))
        ax.set_ylim(-0.02, 0.02)
        ax.set_xlabel("distance from O_d along O···O / Å")
        ax.set_ylabel("Δρ / e bohr⁻³")
        ax.set_title("Stage-resolved density change")
        ax.legend(fontsize=8)
        fig.tight_layout()
        C.FIG_DIR.mkdir(exist_ok=True)
        fig.savefig(C.FIG_DIR / f"density_{args.tag}.png", dpi=160)

        if len(records) > 1:
            fig, ax = plt.subplots(figsize=(5, 4))
            rr = [x["R_OO"] for x in records]
            ax.plot(rr, [x["saddle"]["rho_0"] for x in records], "o-", label="ρ₀ at saddle")
            ax.plot(rr, [x["saddle"]["d_rho_pauli"] for x in records], "s-", label="Δρ_Pauli at saddle")
            ax.plot(rr, [x["saddle"]["d_rho_relax"] for x in records], "^-", label="Δρ_relax at saddle")
            ax.plot(rr, [x["saddle"]["rho_final"] for x in records], "k--", label="ρ_final at saddle")
            ax.axhline(0, color="k", lw=0.6)
            ax.set_xlabel("R(O···O) / Å")
            ax.set_ylabel("e bohr⁻³")
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(C.FIG_DIR / f"density_saddle_vs_R_{args.tag}.png", dpi=160)
        print(f"# figures in {C.FIG_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
