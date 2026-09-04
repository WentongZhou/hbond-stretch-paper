"""Compare the force-constant partition of the target water-water H-bond across
environments (WP4, hypotheses H3 and H2b).

  * table of k_X, f_X, omega and charge transfer for every environment
  * change relative to a reference (dimer, or the matching chain-trimer
    control): which DM-EDA channel carries Delta k?
  * regression ln k = a + n ln(delta q) across environments for the Mulliken
    and IAO charge-transfer proxies: is the CVS paper's exponent 2 recovered?

Usage:
    python scripts/09_compare_environments.py            (auto-detects results/fc_*.json)
    python scripts/09_compare_environments.py --env dimer=results/fc_revpbe0-d3bj_tzvp.json --env acid=results/fc_acid_revpbe0-d3bj_tzvp.json ...
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

CHANNELS = ("Elec", "ExRep", "OrbRel", "CorrDisp")
DEFAULT_ENVS = [
    ("dimer", "fc_revpbe0-d3bj_tzvp.json"),
    ("dimer OH+0.005", "fc_revpbe0-d3bj_tzvp_oh+0.005.json"),
    ("dimer OH+0.010", "fc_revpbe0-d3bj_tzvp_oh+0.010.json"),
    ("dimer OH+0.020", "fc_revpbe0-d3bj_tzvp_oh+0.020.json"),
    ("trimer acid-ctrl", "fc_trimer_acid_control_revpbe0-d3bj_tzvp.json"),
    ("trimer base-ctrl", "fc_trimer_base_control_revpbe0-d3bj_tzvp.json"),
    ("acid H3O+", "fc_acid_revpbe0-d3bj_tzvp.json"),
    ("base OH-", "fc_base_revpbe0-d3bj_tzvp.json"),
]
CONTROLS = {"acid H3O+": "trimer acid-ctrl", "base OH-": "trimer base-ctrl"}


def load_envs(pairs):
    envs = []
    for label, path in pairs:
        p = Path(path)
        if not p.exists():
            p = C.RESULT_DIR / Path(path).name
        if not p.exists():
            print(f"# skipping {label}: {p} not found")
            continue
        d = C.load_json(p)
        ct = d["charge_transfer_acceptor"]
        envs.append({
            "label": label, "tag": d["tag"], "R_min": d["R_min_A"], "E_int": d["E_int_min_kcal"],
            "k_total": d["k_total_kcal_per_A2"], "omega": d["omega_H2O_cm-1"],
            "k": {ch: d["components"][ch]["k_kcal_per_A2"] for ch in CHANNELS},
            "f": {ch: d["components"][ch]["fraction_of_k"] for ch in CHANNELS},
            "slope": {ch: d["components"][ch]["slope_kcal_per_A"] for ch in CHANNELS},
            "ct_mull": ct["mulliken"]["value_e"] if ct.get("mulliken") else np.nan,
            "ct_iao": ct["iao"]["value_e"] if ct.get("iao") else np.nan,
            "ct_label": d.get("ct_label", "acceptor"),
        })
    return envs


def loglog_fit(x, y):
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    lx, ly = np.log(x[m]), np.log(y[m])
    n, a = np.polyfit(lx, ly, 1)
    pred = a + n * lx
    r2 = 1 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return float(n), float(r2), int(m.sum())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--env", action="append", default=None, help="label=path, repeatable")
    parser.add_argument("--out", default="wp4_environments")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    pairs = [tuple(e.split("=", 1)) for e in args.env] if args.env else DEFAULT_ENVS
    envs = load_envs(pairs)
    if not envs:
        print("no environments found")
        return 1
    by_label = {e["label"]: e for e in envs}

    lines = ["# Target H-bond force constant across environments (WP4)", "",
             "k_X = d²X/dR² at R_min of the target O···O (kcal/mol/Å²), f_X = k_X/k_total. "
             "CT = electrons lost by the acceptor fragment of the target bond at R_min.", "",
             "| environment | R_min Å | E_int | k_total | ω(H₂O) | k_Elec | k_ExRep | k_OrbRel | k_CorrDisp | f_Elec | f_ExRep | f_OrbRel | f_CorrDisp | CT Mull | CT IAO |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for e in envs:
        lines.append(f"| {e['label']} | {e['R_min']:.4f} | {e['E_int']:.3f} | {e['k_total']:.3f} | {e['omega']:.1f} | "
                     + " | ".join(f"{e['k'][ch]:.2f}" for ch in CHANNELS) + " | "
                     + " | ".join(f"{e['f'][ch]:+.3f}" for ch in CHANNELS)
                     + f" | {e['ct_mull']:+.4f} | {e['ct_iao']:+.4f} |")

    # Channel attribution of Delta k relative to a reference
    lines += ["", "## Where does Δk come from?", "",
              "Δk_X = k_X(env) − k_X(ref); share = Δk_X / Δk_total. A share > 1 means the channel over-explains the change and is compensated by the others.", "",
              "| environment | reference | Δk_total | Δk/k | Δω cm⁻¹ | ΔR_min Å | Δk_Elec (share) | Δk_ExRep (share) | Δk_OrbRel (share) | Δk_CorrDisp (share) | ΔCT IAO |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    delta_rows = []
    for e in envs:
        refs = []
        if e["label"] != "dimer" and "dimer" in by_label:
            refs.append("dimer")
        if e["label"] in CONTROLS and CONTROLS[e["label"]] in by_label:
            refs.append(CONTROLS[e["label"]])
        for ref_label in refs:
            r = by_label[ref_label]
            dk = e["k_total"] - r["k_total"]
            cells = []
            for ch in CHANNELS:
                d = e["k"][ch] - r["k"][ch]
                share = d / dk if abs(dk) > 1e-6 else np.nan
                cells.append(f"{d:+.2f} ({share:+.2f})")
            delta_rows.append({"env": e["label"], "ref": ref_label, "dk": dk,
                               "dk_channels": {ch: e["k"][ch] - r["k"][ch] for ch in CHANNELS}})
            lines.append(f"| {e['label']} | {ref_label} | {dk:+.3f} | {dk / r['k_total']:+.3f} | {e['omega'] - r['omega']:+.1f} | "
                         f"{e['R_min'] - r['R_min']:+.4f} | " + " | ".join(cells) + f" | {e['ct_iao'] - r['ct_iao']:+.4f} |")

    # log-log regression of k against the CT proxies
    k = np.array([e["k_total"] for e in envs])
    fits = {}
    lines += ["", "## Apparent exponent n in k ∝ δq^n across environments", "",
              "| proxy | n | R² | points |", "|---|---:|---:|---:|"]
    for name, key in (("Mulliken CT", "ct_mull"), ("IAO CT", "ct_iao")):
        q = np.array([e[key] for e in envs])
        n, r2, npts = loglog_fit(q, k)
        fits[name] = {"n": n, "R2": r2, "n_points": npts}
        lines.append(f"| {name} | {n:.2f} | {r2:.3f} | {npts} |")
    lines += ["", "Attractive channels at R_min (share of the restoring force −dX/dR among channels with dX/dR > 0):", "",
              "| environment | Elec | OrbRel | CorrDisp |", "|---|---:|---:|---:|"]
    for e in envs:
        pos = {ch: max(e["slope"][ch], 0.0) for ch in ("Elec", "OrbRel", "CorrDisp")}
        tot = sum(pos.values()) or np.nan
        lines.append(f"| {e['label']} | " + " | ".join(f"{pos[ch] / tot:.3f}" for ch in ("Elec", "OrbRel", "CorrDisp")) + " |")

    text = "\n".join(lines) + "\n"
    out_md = C.RESULT_DIR / f"{args.out}.md"
    out_md.write_text(text)
    C.dump_json(C.RESULT_DIR / f"{args.out}.json", {"environments": envs, "deltas": delta_rows, "loglog_fits": fits})
    print(text)
    print(f"# wrote {out_md}")

    if not args.no_plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
        ax = axes[0]
        labels = [e["label"] for e in envs]
        x = np.arange(len(envs))
        width = 0.18
        colors = {"Elec": "tab:blue", "ExRep": "tab:red", "OrbRel": "tab:green", "CorrDisp": "tab:purple"}
        for i, ch in enumerate(CHANNELS):
            ax.bar(x + (i - 1.5) * width, [e["k"][ch] for e in envs], width, color=colors[ch], label=ch)
        ax.plot(x, [e["k_total"] for e in envs], "k_", ms=18, mew=2, label="k_total")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
        ax.set_ylabel("k_X / kcal mol⁻¹ Å⁻²")
        ax.set_title("Force-constant partition of the target H-bond")
        ax.legend(fontsize=8, ncol=3)

        ax = axes[1]
        for name, key, marker in (("Mulliken", "ct_mull", "o"), ("IAO", "ct_iao", "s")):
            q = np.array([e[key] for e in envs])
            ax.plot(q, k, marker, label=f"{name} (n = {fits[name + ' CT']['n']:.2f})")
            n, r2, _ = loglog_fit(q, k)
            if np.isfinite(n):
                m = np.isfinite(q) & (q > 0)
                a = np.mean(np.log(k[m]) - n * np.log(q[m]))
                qq = np.linspace(q[m].min() * 0.9, q[m].max() * 1.1, 50)
                ax.plot(qq, np.exp(a) * qq ** n, "-", lw=0.8, color="gray")
        for e in envs:
            ax.annotate(e["label"], (e["ct_iao"], e["k_total"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("δq = electrons lost by the acceptor / e")
        ax.set_ylabel("k_total / kcal mol⁻¹ Å⁻²")
        ax.set_title("k against charge transfer across environments")
        ax.legend(fontsize=8)
        fig.tight_layout()
        C.FIG_DIR.mkdir(exist_ok=True)
        fig.savefig(C.FIG_DIR / f"{args.out}.png", dpi=160)
        print(f"# figure: {C.FIG_DIR / f'{args.out}.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
