"""Compare the liquid statistics of neutral water, acid and base (WP6): one table
of means (with standard errors) per level, plus the environment Δk shares and
the exponents.  Reads results/liquid_<tag>_analysis.json for the given tags.

Usage:
    python scripts/14_liquid_compare.py --tags liquid liquid_h3o liquid_oh
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

CHANNELS = ("Elec", "ExRep", "OrbRel", "CorrDisp")
LABEL = {"liquid": "neutral water", "liquid_h3o": "next to H3O+ (shell 1 => shell 2)", "liquid_oh": "next to OH- (shell 2 => shell 1)"}


def mse(v):
    v = np.array([x for x in v if x is not None and np.isfinite(x)], dtype=float)
    if len(v) == 0:
        return "n/a"
    se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    d = 4 if np.abs(v).max() < 0.2 else 2
    return f"{v.mean():.{d}f} ± {se:.{d}f}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tags", nargs="+", default=["liquid", "liquid_h3o", "liquid_oh"])
    parser.add_argument("--out", default="liquid_compare")
    args = parser.parse_args()
    data = {}
    for tag in args.tags:
        p = C.RESULT_DIR / f"liquid_{tag}_analysis.json"
        if p.exists():
            data[tag] = C.load_json(p)
    lines = ["# Liquid H-bonds: neutral water vs acid vs base (means ± standard error)", "",
             "k in kcal/mol/Å², ω in cm⁻¹ (μ = m(H₂O)/2), CT = IAO electrons lost by the acceptor of the target bond.", ""]
    for kind in ("pair", "embed", "full"):
        lines += [f"## {kind}", "",
                  "| environment | n | R_snap Å | R_min Å | k(min) | ω(min) | k(R0) | ω(R0) | CT(min) | f_Elec | f_ExRep | f_OrbRel | f_CorrDisp |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for tag, d in data.items():
            rows = [r for r in d["rows"] if r.get(kind)]
            if not rows:
                continue
            a = [r[kind]["at_R0"] for r in rows if r[kind].get("at_R0")]
            lines.append(f"| {LABEL.get(tag, tag)} | {len(rows)} | {mse([r['desc']['R_OO'] for r in rows])} | "
                         f"{mse([r[kind]['R_min'] for r in rows])} | {mse([r[kind]['k'] for r in rows])} | "
                         f"{mse([r[kind]['omega'] for r in rows])} | {mse([x['k_total_kcal_per_A2'] for x in a])} | "
                         f"{mse([x['omega_H2O_cm-1'] for x in a])} | {mse([r[kind]['ct_iao'] for r in rows])} | "
                         + " | ".join(mse([r[kind]["fX"][ch] for r in rows]) for ch in CHANNELS) + " |")
        lines.append("")
    lines += ["## Environment effects (mean over bonds; shares Δk_X/Δk for |Δk| > 0.5)", "",
              "| environment | Δk/k field (embed − pair) | shares Elec / ExRep / OrbRel / CorrDisp | Δk/k cage+QM (full − embed) | shares | ΔCT field |",
              "|---|---:|---|---:|---|---:|"]
    for tag, d in data.items():
        rows = d["rows"]
        fe = [r for r in rows if r.get("embed") and r.get("pair")]
        cq = [r for r in rows if r.get("full") and r.get("embed")]

        def shares(sel, hi, lo):
            out = []
            for ch in CHANNELS:
                v = [(r[hi]["kX"][ch] - r[lo]["kX"][ch]) / (r[hi]["k"] - r[lo]["k"]) for r in sel if abs(r[hi]["k"] - r[lo]["k"]) > 0.5]
                out.append(f"{np.mean(v):+.2f}" if v else "—")
            return " / ".join(out)

        lines.append(f"| {LABEL.get(tag, tag)} | {mse([(r['embed']['k'] - r['pair']['k']) / r['pair']['k'] for r in fe])} | {shares(fe, 'embed', 'pair')} | "
                     f"{mse([(r['full']['k'] - r['embed']['k']) / r['embed']['k'] for r in cq])} | {shares(cq, 'full', 'embed')} | "
                     f"{mse([r['embed']['ct_iao'] - r['pair']['ct_iao'] for r in fe])} |")
    lines += ["", "## Exponent n in k ∝ δq^n (IAO) across bonds", "", "| environment | pair | embed | full vs embed δq |", "|---|---:|---:|---:|"]
    for tag, d in data.items():
        f = d["fits"]
        def fmt(key):
            x = f.get(key)
            return f"{x['n']:.2f} (R² {x['R2']:.2f}, n={x['n_points']})" if x and np.isfinite(x["n"]) else "—"
        lines.append(f"| {LABEL.get(tag, tag)} | {fmt('pair_IAO')} | {fmt('embed_IAO')} | {fmt('full_vs_embedCT_IAO')} |")
    # pooled regression across all environments (embed level)
    ks, qs = [], []
    for d in data.values():
        for r in d["rows"]:
            if r.get("embed"):
                ks.append(r["embed"]["k"]); qs.append(r["embed"]["ct_iao"])
    ks, qs = np.array(ks), np.array(qs)
    m = np.isfinite(qs) & (qs > 0) & (ks > 0)
    if m.sum() >= 3:
        n, a = np.polyfit(np.log(qs[m]), np.log(ks[m]), 1)
        pred = a + n * np.log(qs[m])
        r2 = 1 - np.sum((np.log(ks[m]) - pred) ** 2) / np.sum((np.log(ks[m]) - np.log(ks[m]).mean()) ** 2)
        lines += ["", f"Pooled embed-level regression over all {m.sum()} liquid bonds: n = {n:.2f}, R² = {r2:.2f}."]
    # figure: k(embed) vs CT(embed) for all bonds by environment, and the field-effect shares
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        ax = axes[0]
        colors = {"liquid": "tab:blue", "liquid_h3o": "tab:red", "liquid_oh": "tab:green"}
        for tag, d in data.items():
            rows = [r for r in d["rows"] if r.get("embed")]
            ax.plot([r["embed"]["ct_iao"] for r in rows], [r["embed"]["k"] for r in rows], "o", ms=5,
                    color=colors.get(tag, "k"), label=LABEL.get(tag, tag), alpha=0.8)
            rows = [r for r in d["rows"] if r.get("pair")]
            ax.plot([r["pair"]["ct_iao"] for r in rows], [r["pair"]["k"] for r in rows], "s", ms=3,
                    color=colors.get(tag, "k"), alpha=0.3)
        if m.sum() >= 3:
            qq = np.linspace(qs[m].min() * 0.9, qs[m].max() * 1.1, 50)
            ax.plot(qq, np.exp(a) * qq ** n, "-", color="gray", lw=1, label=f"pooled embed fit, n = {n:.2f}")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("δq (IAO electrons lost by the acceptor) / e")
        ax.set_ylabel("k_total / kcal mol⁻¹ Å⁻²")
        ax.set_title("Liquid H-bonds: k against charge transfer (embed: circles, pair: squares)", fontsize=9)
        ax.legend(fontsize=7)
        ax = axes[1]
        width = 0.25
        for i, (tag, d) in enumerate(data.items()):
            fe = [r for r in d["rows"] if r.get("embed") and r.get("pair") and abs(r["embed"]["k"] - r["pair"]["k"]) > 0.5]
            vals = [np.mean([(r["embed"]["kX"][ch] - r["pair"]["kX"][ch]) / (r["embed"]["k"] - r["pair"]["k"]) for r in fe]) for ch in CHANNELS]
            ax.bar(np.arange(len(CHANNELS)) + (i - 1) * width, vals, width, color=colors.get(tag, "k"), label=LABEL.get(tag, tag))
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(np.arange(len(CHANNELS)))
        ax.set_xticklabels(CHANNELS)
        ax.set_ylabel("share of Δk (embed − pair)")
        ax.set_title("Which channel carries the field-induced stiffening", fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout()
        C.FIG_DIR.mkdir(exist_ok=True)
        fig.savefig(C.FIG_DIR / f"{args.out}.png", dpi=160)
        lines.append(f"\nFigure: figures/{args.out}.png")
    except Exception as exc:  # plotting is optional
        lines.append(f"\n(figure skipped: {exc!r})")

    text = "\n".join(lines) + "\n"
    out = C.RESULT_DIR / f"{args.out}.md"
    out.write_text(text)
    print(text)
    print(f"# wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
