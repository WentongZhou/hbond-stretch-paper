"""Statistics of the force-constant partition over liquid-snapshot H-bonds (WP6).

Reads results/liquid/index_<tag>.json and, for every cluster, the force-constant
files results/fc_<name>_full.json (donor + environment / acceptor) and
results/fc_<name>_pair.json (the two waters alone at the same geometry) if present.

Outputs results/liquid_<tag>_analysis.md/.json and figures/liquid_<tag>.png:
  * distribution (mean, sd, min, max) of k_total, omega, f_X, CT for full and pair
  * per-bond Delta k (full - pair) channel shares: what the environment does
  * regression ln k = a + n ln(delta q) over all bonds, IAO and Mulliken, full and pair
  * correlation of k with the snapshot geometry (R_OO, angle) and coordination

Usage:
    python scripts/12_liquid_analysis.py --tag liquid
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

CHANNELS = ("Elec", "ExRep", "OrbRel", "CorrDisp")


def load_fc(path):
    if not Path(path).exists():
        return None
    d = C.load_json(path)
    ct = d["charge_transfer_acceptor"]
    return {
        "R_min": d["R_min_A"], "E_int": d["E_int_min_kcal"], "k": d["k_total_kcal_per_A2"], "omega": d["omega_H2O_cm-1"],
        "kX": {ch: d["components"][ch]["k_kcal_per_A2"] for ch in CHANNELS},
        "fX": {ch: d["components"][ch]["fraction_of_k"] for ch in CHANNELS},
        "ct_mull": ct["mulliken"]["value_e"] if ct.get("mulliken") else np.nan,
        "ct_iao": ct["iao"]["value_e"] if ct.get("iao") else np.nan,
        "n_points": d["n_points"],
        "minimum_found": d.get("minimum_found", True),
        "at_R0": d.get("at_R0"),
    }


def loglog(x, y):
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if m.sum() < 3:
        return np.nan, np.nan, np.nan, int(m.sum())
    lx, ly = np.log(x[m]), np.log(y[m])
    n, a = np.polyfit(lx, ly, 1)
    pred = a + n * lx
    r2 = 1 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return float(n), float(a), float(r2), int(m.sum())


def stats(values):
    v = np.array([x for x in values if np.isfinite(x)], dtype=float)
    if len(v) == 0:
        return "n/a"
    return f"{v.mean():.3f} ± {v.std(ddof=1) if len(v) > 1 else 0:.3f} [{v.min():.3f}, {v.max():.3f}] (n={len(v)})"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="liquid")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    index = C.load_json(C.RESULT_DIR / "liquid" / f"index_{args.tag}.json")
    rows = []
    for entry in index:
        name = entry["name"]
        full = load_fc(C.RESULT_DIR / f"fc_{name}_full.json")
        pair = load_fc(C.RESULT_DIR / f"fc_{name}_pair.json")
        embed = load_fc(C.RESULT_DIR / f"fc_{name}_embed.json")
        if full is None and pair is None and embed is None:
            continue
        rows.append({"name": name, "desc": entry, "full": full, "pair": pair, "embed": embed})
    if not rows:
        print("no force-constant files found yet")
        return 1

    lines = [f"# Liquid-snapshot H-bonds: force-constant partition statistics ({args.tag})", "",
             f"{len(rows)} target H-bonds with results ({sum(r['full'] is not None for r in rows)} full, "
             f"{sum(r['embed'] is not None for r in rows)} embedded, {sum(r['pair'] is not None for r in rows)} pair).  "
             "k in kcal/mol/Å², ω in cm⁻¹, CT = electrons lost by the acceptor fragment of the target bond (IAO / Mulliken).", "",
             "Three levels: pair = the two waters alone at the liquid geometry; embed = the same two waters inside the point charges "
             "of everything within 8 Å (TIP4P-Ew sites, model ions); full = donor + first shells as QM fragment / acceptor water "
             "(includes the cage of the acceptor's other neighbours).", ""]
    KIND_DESC = {"full": "donor + first shells (QM) / acceptor water", "embed": "two waters in the point-charge environment",
                 "pair": "two waters only, frozen liquid geometry"}
    for kind in ("pair", "embed", "full"):
        sel = [r[kind] for r in rows if r[kind] is not None]
        if not sel:
            continue
        lines += [f"## {kind}: {KIND_DESC[kind]}", "",
                  "| quantity | mean ± sd [min, max] |", "|---|---|",
                  f"| R_min Å | {stats([s['R_min'] for s in sel])} |",
                  f"| k_total | {stats([s['k'] for s in sel])} |",
                  f"| ω(H₂O) | {stats([s['omega'] for s in sel])} |"]
        for ch in CHANNELS:
            lines.append(f"| f_{ch} | {stats([s['fX'][ch] for s in sel])} |")
        lines += [f"| CT IAO | {stats([s['ct_iao'] for s in sel])} |",
                  f"| CT Mulliken | {stats([s['ct_mull'] for s in sel])} |", ""]

    def delta_block(kind_hi, kind_lo, title):
        both = [r for r in rows if r[kind_hi] is not None and r[kind_lo] is not None]
        if not both:
            return {}
        lines.extend(["", f"## {title}: Δk = k({kind_hi}) − k({kind_lo}), channel shares Δk_X/Δk", "",
                      "| bond | R_snap Å | angle | env waters | k lo | k hi | Δk | Elec | ExRep | OrbRel | CorrDisp | CT IAO lo→hi |",
                      "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"])
        shares = {ch: [] for ch in CHANNELS}
        for r in both:
            dk = r[kind_hi]["k"] - r[kind_lo]["k"]
            cells = []
            for ch in CHANNELS:
                d = r[kind_hi]["kX"][ch] - r[kind_lo]["kX"][ch]
                sh = d / dk if abs(dk) > 0.5 else np.nan
                shares[ch].append(sh)
                cells.append(f"{sh:+.2f}" if np.isfinite(sh) else "—")
            de = r["desc"]
            lines.append(f"| {r['name']} | {de['R_OO']:.3f} | {de['angle_OHO']:.0f} | {de.get('n_env_waters', '')} | "
                         f"{r[kind_lo]['k']:.2f} | {r[kind_hi]['k']:.2f} | {dk:+.2f} | " + " | ".join(cells) +
                         f" | {r[kind_lo]['ct_iao']:.3f}→{r[kind_hi]['ct_iao']:.3f} |")
        lines.extend(["", "Mean share of Δk over bonds with |Δk| > 0.5: " +
                      ", ".join(f"{ch} {np.nanmean(shares[ch]):+.2f}" for ch in CHANNELS),
                      "Mean Δk/k(lo): " + f"{np.mean([(r[kind_hi]['k'] - r[kind_lo]['k']) / r[kind_lo]['k'] for r in both]):+.2f}"])
        return {ch: [x for x in shares[ch] if np.isfinite(x)] for ch in CHANNELS}

    shares_field = delta_block("embed", "pair", "Electrostatic environment effect")
    shares_cage = delta_block("full", "embed", "QM environment + cage on top of the field")
    shares_total = delta_block("full", "pair", "Total environment effect")

    lines += ["## Apparent exponent n in k ∝ δq^n across liquid H-bonds", "",
              "| set | proxy | n | R² | points |", "|---|---|---:|---:|---:|"]
    fits = {}
    for kind in ("pair", "embed", "full"):
        sel = [r[kind] for r in rows if r[kind] is not None]
        if len(sel) < 3:
            continue
        k = np.array([s["k"] for s in sel])
        for pname, key in (("IAO", "ct_iao"), ("Mulliken", "ct_mull")):
            q = np.array([s[key] for s in sel])
            n, a, r2, npts = loglog(q, k)
            fits[f"{kind}_{pname}"] = {"n": n, "a": a, "R2": r2, "n_points": npts}
            lines.append(f"| {kind} | {pname} | {n:.2f} | {r2:.3f} | {npts} |")
    # k(full) against the target-bond CT taken from the embedded pair (the full-cluster CT is not the target bond's)
    sel = [r for r in rows if r["full"] is not None and r["embed"] is not None]
    if len(sel) >= 3:
        k = np.array([r["full"]["k"] for r in sel])
        q = np.array([r["embed"]["ct_iao"] for r in sel])
        n, a, r2, npts = loglog(q, k)
        fits["full_vs_embedCT_IAO"] = {"n": n, "a": a, "R2": r2, "n_points": npts}
        lines.append(f"| full k vs embed δq | IAO | {n:.2f} | {r2:.3f} | {npts} |")
    # geometry correlations
    lines += ["", "## Correlation of k with snapshot geometry", "",
              "| set | corr(k, R_OO snapshot) | corr(k, angle) | corr(k, n_env) | corr(k, R_min) |", "|---|---:|---:|---:|---:|"]
    for kind in ("pair", "embed", "full"):
        sel = [r for r in rows if r[kind] is not None]
        if len(sel) < 3:
            continue
        k = np.array([r[kind]["k"] for r in sel])
        def corr(x):
            x = np.array(x, dtype=float)
            return float(np.corrcoef(x, k)[0, 1]) if np.std(x) > 0 else np.nan
        lines.append(f"| {kind} | {corr([r['desc']['R_OO'] for r in sel]):+.2f} | {corr([r['desc']['angle_OHO'] for r in sel]):+.2f} | "
                     f"{corr([r['desc']['n_env_waters'] for r in sel]):+.2f} | {corr([r[kind]['R_min'] for r in sel]):+.2f} |")

    # instantaneous curvature at the snapshot geometry (all bonds, no minimum needed)
    lines += ["", "## Curvature at the snapshot geometry R0 (instantaneous, all bonds)", "",
              "| set | n | k(R0) | ω(R0) | f_Elec | f_ExRep | f_OrbRel | f_CorrDisp | bonds without a minimum in range |",
              "|---|---:|---|---|---|---|---|---|---:|"]
    for kind in ("pair", "embed", "full"):
        sel = [r[kind] for r in rows if r[kind] is not None and r[kind].get("at_R0")]
        if not sel:
            continue
        a = [x["at_R0"] for x in sel]
        nomin = sum(1 for x in sel if not x.get("minimum_found", True))
        lines.append(f"| {kind} | {len(a)} | {stats([x['k_total_kcal_per_A2'] for x in a])} | {stats([x['omega_H2O_cm-1'] for x in a])} | "
                     + " | ".join(stats([x["components"][ch]["fraction_of_k"] for x in a]) for ch in CHANNELS)
                     + f" | {nomin} |")
    k0 = {}
    for kind in ("pair", "embed", "full"):
        sel = [r for r in rows if r[kind] is not None and r[kind].get("at_R0")]
        if len(sel) >= 3:
            k0[kind] = (np.array([r[kind]["at_R0"]["k_total_kcal_per_A2"] for r in sel]),
                        np.array([r[kind]["at_R0"]["charge_transfer_acceptor"]["iao"]["value_e"] if r[kind]["at_R0"]["charge_transfer_acceptor"].get("iao") else np.nan for r in sel]))
    if k0:
        lines += ["", "ln k(R0) against ln δq(R0), IAO: " + ", ".join(
            f"{kind} n = {loglog(q, k)[0]:.2f} (R² {loglog(q, k)[2]:.2f}, {loglog(q, k)[3]} pts)" for kind, (k, q) in k0.items())]

    text = "\n".join(lines) + "\n"
    out = C.RESULT_DIR / f"liquid_{args.tag}_analysis.md"
    out.write_text(text)
    C.dump_json(C.RESULT_DIR / f"liquid_{args.tag}_analysis.json", {"rows": rows, "fits": fits})
    print(text)
    print(f"# wrote {out}")

    if not args.no_plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
        ax = axes[0]
        for kind, marker, alpha in (("pair", "s", 0.35), ("embed", "^", 0.6), ("full", "o", 0.9)):
            sel = [r[kind] for r in rows if r[kind] is not None]
            if not sel:
                continue
            for ch, color in zip(CHANNELS, ("tab:blue", "tab:red", "tab:green", "tab:purple")):
                ax.plot([s["k"] for s in sel], [s["fX"][ch] for s in sel], marker, ms=4, color=color,
                        label=f"{ch} ({kind})", alpha=alpha)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel("k_total / kcal mol⁻¹ Å⁻²")
        ax.set_ylabel("f_X = k_X / k_total")
        ax.set_title("Channel fractions across liquid H-bonds")
        ax.legend(fontsize=6, ncol=2)

        ax = axes[1]
        for kind, marker in (("pair", "s"), ("embed", "^"), ("full", "o")):
            sel = [r[kind] for r in rows if r[kind] is not None]
            if len(sel) < 3:
                continue
            q = np.array([s["ct_iao"] for s in sel])
            k = np.array([s["k"] for s in sel])
            f = fits.get(f"{kind}_IAO")
            lab = f"{kind} (n = {f['n']:.2f}, R² = {f['R2']:.2f})" if f else kind
            ax.plot(q, k, marker, ms=4, label=lab)
            if f and np.isfinite(f["n"]):
                qq = np.linspace(np.nanmin(q) * 0.9, np.nanmax(q) * 1.1, 50)
                ax.plot(qq, np.exp(f["a"]) * qq ** f["n"], "-", lw=0.8, color="gray")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("δq (IAO, electrons lost by the acceptor) / e")
        ax.set_ylabel("k_total / kcal mol⁻¹ Å⁻²")
        ax.set_title("k against charge transfer, liquid H-bonds")
        ax.legend(fontsize=8)

        ax = axes[2]
        pos = 0
        ticks, labels = [], []
        for (hi, lo, lab, color) in (("embed", "pair", "field", "tab:orange"), ("full", "embed", "QM env + cage", "tab:gray")):
            sel = [r for r in rows if r[hi] is not None and r[lo] is not None]
            if not sel:
                continue
            data = [[r[hi]["kX"][ch] - r[lo]["kX"][ch] for r in sel] for ch in CHANNELS]
            bp = ax.boxplot(data, positions=[pos + i for i in range(len(CHANNELS))], widths=0.5, patch_artist=True)
            for b in bp["boxes"]:
                b.set_facecolor(color)
                b.set_alpha(0.5)
            ticks += [pos + i for i in range(len(CHANNELS))]
            labels += [ch + "\n" + lab for ch in CHANNELS]
            pos += len(CHANNELS) + 1
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=7)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_ylabel("Δk_X / kcal mol⁻¹ Å⁻²")
        ax.set_title("What the environment adds, by channel")
        fig.tight_layout()
        C.FIG_DIR.mkdir(exist_ok=True)
        fig.savefig(C.FIG_DIR / f"liquid_{args.tag}.png", dpi=160)
        print(f"# figure: {C.FIG_DIR / f'liquid_{args.tag}.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
