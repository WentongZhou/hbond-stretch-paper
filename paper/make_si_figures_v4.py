"""SI figures, revision v4 (labels "control" for the geometry-matched bare-pair controls; output figures_paper/v4).

Supporting-Information figures with clean legends (no legend over data).

    figures_paper/figS1_gas.png          gas-phase environments: k_X bars and k against charge transfer
    figures_paper/figS2_liquid.png       neutral liquid, 20 bonds
    figures_paper/figS3_liquid_h3o.png   next to H3O+, 10 bonds
    figures_paper/figS4_liquid_oh.png    next to OH-, 10 bonds
Usage:  python make_si_figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
D = HERE / "data"
OUT = HERE / "figures_paper" / "v4"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "legend.fontsize": 6.5, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "axes.linewidth": 0.6,
    "lines.linewidth": 1.0, "figure.dpi": 100, "savefig.dpi": 300, "mathtext.fontset": "dejavusans",
    "legend.frameon": False, "legend.handlelength": 1.4,
})
COL = {"Elec": "#1f77b4", "ExRep": "#d62728", "OrbRel": "#2ca02c", "CorrDisp": "#9467bd"}
NAME = {"Elec": "electrostatics", "ExRep": "exchange–repulsion", "OrbRel": "orbital relaxation", "CorrDisp": "correlation + dispersion"}
CH = ("Elec", "ExRep", "OrbRel", "CorrDisp")


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def panel_label(ax, s, x=-0.2):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def loglog_fit(q, k):
    m = np.isfinite(q) & np.isfinite(k) & (q > 0) & (k > 0)
    n, a = np.polyfit(np.log(q[m]), np.log(k[m]), 1)
    pred = a + n * np.log(q[m])
    r2 = 1 - np.sum((np.log(k[m]) - pred) ** 2) / np.sum((np.log(k[m]) - np.log(k[m]).mean()) ** 2)
    return n, a, r2


def log_axes(ax, xt, yt):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(FixedLocator(xt))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_major_locator(FixedLocator(yt))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


# ---------------------------------------------------------------------------
def figS1():
    wp4 = load(D / "analysis" / "wp4_environments.json")
    envs = wp4["environments"]
    labels = {"dimer": "dimer", "dimer OH+0.010": "dimer, O–H +0.01 Å", "trimer acid-ctrl": "trimer Wd⇒Wa",
              "trimer base-ctrl": "trimer W0⇒Wd", "acid H3O+": "H₃O⁺ cluster", "acid pair": "H₃O⁺ control",
              "base OH-": "OH⁻ cluster", "base pair": "OH⁻ control"}
    order = [e for e in envs if e["label"] in labels]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.35, "left": 0.08, "right": 0.98, "top": 0.82, "bottom": 0.27})
    ax = axes[0]
    x = np.arange(len(order))
    w = 0.19
    for i, ch in enumerate(CH):
        ax.bar(x + (i - 1.5) * w, [e["k"][ch] for e in order], w, color=COL[ch], label=NAME[ch])
    ax.plot(x, [e["k_total"] for e in order], "_", color="k", ms=12, mew=1.5, label=r"$k_{\mathrm{total}}$")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[e["label"]] for e in order], fontsize=6, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_ylabel(r"$k_X$ / kcal mol⁻¹ Å⁻²")
    ax.set_ylim(-95, 215)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=6, columnspacing=1.0)
    panel_label(ax, "a", x=-0.14)

    ax = axes[1]
    k = np.array([e["k_total"] for e in order])
    for key, marker, color, name in (("ct_iao", "o", "#ff7f0e", "IAO"), ("ct_mull", "s", "0.4", "Mulliken")):
        q = np.array([e[key] for e in order])
        n, a, r2 = loglog_fit(q, k)
        ax.plot(q, k, marker, color=color, ms=4, label=f"{name}: n = {n:.2f}, R² = {r2:.2f}")
        qq = np.linspace(q.min() * 0.85, q.max() * 1.15, 40)
        ax.plot(qq, np.exp(a) * qq ** n, "-", color=color, lw=0.8)
    log_axes(ax, [0.02, 0.03, 0.05, 0.08], [20, 30, 50])
    ax.set_xlim(0.014, 0.12)
    ax.set_ylim(12, 70)
    ax.set_xlabel("δq, electrons lost by the acceptor / e")
    ax.set_ylabel(r"$k_{\mathrm{total}}$ / kcal mol⁻¹ Å⁻²")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1, title="k ∝ δqⁿ over the eight environments", title_fontsize=6.3, fontsize=6.3)
    panel_label(ax, "b", x=-0.2)
    fig.savefig(OUT / "figS1_gas.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
def fig_liquid(tag, fname, title):
    d = load(D / "analysis" / f"liquid_{tag}_analysis.json")
    rows = d["rows"]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2), gridspec_kw={"hspace": 0.55, "wspace": 0.32, "left": 0.09, "right": 0.98, "top": 0.86, "bottom": 0.09})
    fig.suptitle(title, fontsize=8.5, y=0.985)

    # (a) channel fractions against k, pair and embed
    ax = axes[0, 0]
    for kind, marker, alpha in (("pair", "s", 0.45), ("embed", "o", 0.9)):
        sel = [r[kind] for r in rows if r.get(kind)]
        for ch in CH:
            ax.plot([s["k"] for s in sel], [s["fX"][ch] for s in sel], marker, color=COL[ch], ms=3.5, alpha=alpha, mew=0)
    ax.axhline(0, color="k", lw=0.5)
    handles = [Patch(color=COL[ch], label=NAME[ch]) for ch in CH] + [
        Line2D([], [], marker="s", color="0.4", ls="", ms=3.5, alpha=0.6, label="bare pair"),
        Line2D([], [], marker="o", color="0.4", ls="", ms=3.5, label="embedded pair")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=5.8, columnspacing=0.9)
    ax.set_xlabel(r"$k_{\mathrm{total}}$ / kcal mol⁻¹ Å⁻²")
    ax.set_ylabel(r"$f_X = k_X / k_{\mathrm{total}}$")
    ax.set_ylim(-2.4, 5.6)
    panel_label(ax, "a")

    # (b) k against dq, pair and embed, with fits; legend in the empty lower right
    ax = axes[0, 1]
    for kind, marker, color, name in (("pair", "s", "0.45", "bare pair"), ("embed", "o", "#ff7f0e", "embedded pair")):
        sel = [r[kind] for r in rows if r.get(kind)]
        q = np.array([s["ct_iao"] for s in sel])
        k = np.array([s["k"] for s in sel])
        n, a, r2 = loglog_fit(q, k)
        ax.plot(q, k, marker, color=color, ms=3.8, label=f"{name}: n = {n:.2f}, R² = {r2:.2f}", mew=0)
        qq = np.linspace(q.min() * 0.85, q.max() * 1.15, 40)
        ax.plot(qq, np.exp(a) * qq ** n, "-", color=color, lw=0.8)
    log_axes(ax, [0.015, 0.02, 0.03, 0.04, 0.06], [3, 5, 10, 20, 40])
    ax.set_xlim(0.011, 0.09)
    ax.set_ylim(2.5, 60)
    ax.set_xlabel("δq, electrons lost by the acceptor (IAO) / e")
    ax.set_ylabel(r"$k_{\mathrm{total}}$ / kcal mol⁻¹ Å⁻²")
    ax.legend(loc="lower right", fontsize=6)
    panel_label(ax, "b")

    # (c) field effect: dk_X = embed - pair ; (d) cage + QM environment: full - embed
    for ax, (hi, lo, lab, letter) in zip(axes[1], (("embed", "pair", "field: embedded − bare pair", "c"),
                                                    ("full", "embed", "QM environment + cage: full − embedded", "d"))):
        sel = [r for r in rows if r.get(hi) and r.get(lo)]
        data = [[r[hi]["kX"][ch] - r[lo]["kX"][ch] for r in sel] for ch in CH]
        bp = ax.boxplot(data, positions=range(len(CH)), widths=0.55, patch_artist=True,
                        medianprops={"color": "k", "lw": 1.0}, flierprops={"marker": ".", "ms": 3})
        for b, ch in zip(bp["boxes"], CH):
            b.set_facecolor(COL[ch])
            b.set_alpha(0.55)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(CH)))
        ax.set_xticklabels(["electro-\nstatics", "exchange–\nrepulsion", "orbital\nrelaxation", "corr. +\ndisp."], fontsize=6.3)
        ax.set_ylabel(r"$\Delta k_X$ / kcal mol⁻¹ Å⁻²")
        ax.set_title(f"{lab} (n = {len(sel)})", fontsize=7)
        panel_label(ax, letter)
    fig.savefig(OUT / fname)
    plt.close(fig)


if __name__ == "__main__":
    figS1()
    fig_liquid("liquid", "figS2_liquid.png", "Neutral liquid water, 20 hydrogen bonds")
    fig_liquid("liquid_h3o", "figS3_liquid_h3o.png", "Next to H₃O⁺, 10 shell-1 ⇒ shell-2 hydrogen bonds")
    fig_liquid("liquid_oh", "figS4_liquid_oh.png", "Next to OH⁻, 10 shell-2 ⇒ shell-1 hydrogen bonds")
    print("SI figures written to", OUT)
