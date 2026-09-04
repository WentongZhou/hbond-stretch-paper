"""Paper figures for the JPC Letters manuscript, revision v4 (drawn from the collected JSON data).

Identical to make_paper_figures_v3.py except for Fig. 1a, which is redrawn from the original scan output
(hbond-stretch-eda/results/scan_revpbe0-d3bj_qzvp.json; the copy in data/scans/ is checked to be identical):
all 39 computed points are shown as markers joined by lines, as in the pipeline figure fc_revpbe0-d3bj_qzvp.png;
the x-range is the scan range (2.50-3.65 Å) so the curves fill the panel; the five direct curve labels of v3
(which collided with the zero line and each other, and "correlation + dispersion" ran past the axes frame)
are replaced by a legend with the channel abbreviations used in panels (b) and (c) (Elec, ExRep, OrbRel,
Corr+Disp, total; defined in the caption) in the empty upper-right region; the R_min label sits at the foot of
the dotted line. Figs. 2, 3 and toc_plain are unchanged.

Outputs (300 dpi PNG, double-column width 7.0 in) in figures_paper/v4/:
    fig1_dimer.png, fig2_density.png, fig3_environments.png, toc_plain.png (data-only TOC alternative)
Usage:  python make_paper_figures_v4.py            # all figures
        python make_paper_figures_v4.py fig1       # only the named figure(s): fig1 fig2 fig3 toc_plain
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

import revision_stats as rs

HERE = Path(__file__).resolve().parent
D = HERE / "data"
OUT = HERE / "figures_paper" / "v4"
OUT.mkdir(parents=True, exist_ok=True)
SRC_RESULTS = HERE.parent / "hbond-stretch-eda" / "results"  # original calculation output (collect_materials.py copies it to data/)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "legend.fontsize": 6.5, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "axes.linewidth": 0.6,
    "lines.linewidth": 1.1, "lines.markersize": 3.5, "figure.dpi": 100, "savefig.dpi": 300,
    "mathtext.fontset": "dejavusans", "legend.frameon": False, "legend.handlelength": 1.6,
})
COL = {"Elec": "#1f77b4", "ExRep": "#d62728", "OrbRel": "#2ca02c", "CorrDisp": "#9467bd", "Total": "k"}
NAME = {"Elec": "electrostatics", "ExRep": "exchange–repulsion", "OrbRel": "orbital relaxation",
        "CorrDisp": "correlation + dispersion", "Total": "total"}
CH = ("Elec", "ExRep", "OrbRel", "CorrDisp")
R0 = r"$\rho_0$"
RP = r"$\rho_{\mathrm{Pauli}}$"
RF = r"$\rho_{\mathrm{final}}$"
DP = r"$\Delta\rho_{\mathrm{Pauli}}$"
DR = r"$\Delta\rho_{\mathrm{relax}}$"
XLAB_R = "$R$(O···O) / Å"
XLAB_OD = "distance from O$_\\mathrm{d}$ along O···O / Å"


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def panel_label(ax, s, x=-0.2):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def load_original(sub, name):
    """Load a result file from the original hbond-stretch-eda output; the copy under data/ must be identical."""
    local = D / sub / name
    src = SRC_RESULTS / name
    if src.exists():
        if src.read_bytes() != local.read_bytes():
            raise RuntimeError(f"{local} differs from the original {src}; rerun collect_materials.py")
        return load(src)
    return load(local)


# ---------------------------------------------------------------------------
def fig1():
    scan = load_original("scans", "scan_revpbe0-d3bj_qzvp.json")
    fc = load_original("fc", "fc_revpbe0-d3bj_qzvp.json")
    cc = load(D / "ccsdt" / "ccsdt_avqz.json")["fit"]["summary"]
    pts = sorted(scan["points"], key=lambda p: p["R_OO"])
    r = np.array([p["R_OO"] for p in pts])
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5), gridspec_kw={"width_ratios": [1.0, 1.15, 1.0]})

    # (a) the 39 computed scan points of each component (markers) joined by lines over the full scan range;
    # abbreviated legend (as in panels b, c) in the empty upper-right region (ExRep < 3.5 kcal/mol for R > 3.1 Å),
    # R_min label at the foot of the dotted line
    ax = axes[0]
    lines = {}
    for key in ("Elec", "ExRep", "OrbRel", "CorrDisp", "Total"):
        y = np.array([p[key] for p in pts])
        lines[key], = ax.plot(r, y, "o-", color=COL[key], ms=2.0, mew=0, lw=1.3 if key == "Total" else 0.9)
    ax.axvline(fc["R_min_A"], color="gray", ls=":", lw=0.7)
    ax.text(fc["R_min_A"] + 0.03, -19.2, r"$R_{\min}$", fontsize=6.5, color="gray", va="bottom", ha="left")
    ax.axhline(0, color="k", lw=0.5)
    order = ("ExRep", "CorrDisp", "OrbRel", "Total", "Elec")  # top-to-bottom order of the curves at large R
    abbr = {"Elec": "Elec", "ExRep": "ExRep", "OrbRel": "OrbRel", "CorrDisp": "Corr+Disp", "Total": "total"}
    ax.legend([lines[k] for k in order], [abbr[k] for k in order], loc="upper right", fontsize=6.2,
              labelspacing=0.35, handlelength=1.4, handletextpad=0.5, borderaxespad=0.3)
    ax.set_xlim(2.50, 3.65)
    ax.set_ylim(-20, 30)
    ax.set_xticks([2.6, 3.0, 3.4])
    ax.set_xlabel(XLAB_R)
    ax.set_ylabel("energy component / kcal mol⁻¹")
    panel_label(ax, "a")

    # (b) curvature partition; CCSD(T) as a short reference mark above the total bar
    ax = axes[1]
    comp = fc["components"]
    keys = ("Elec", "ExRep", "OrbRel", "CorrDisp")
    vals = [comp[k]["k_kcal_per_A2"] for k in keys]
    ktot = fc["k_total_kcal_per_A2"]
    ax.bar(range(len(keys)), vals, color=[COL[k] for k in keys], width=0.62)
    ax.bar(len(keys), ktot, color="0.35", width=0.62)
    kcc = cc["CCSD(T)"]["k_kcal_per_A2"]
    ax.plot([len(keys) - 0.42, len(keys) + 0.42], [kcc, kcc], "k-", lw=1.2)
    ax.text(len(keys), kcc + 3, "CCSD(T)", ha="center", va="bottom", fontsize=6)
    for i, v in enumerate(vals):
        ax.text(i, v + (2.5 if v > 0 else -2.5), f"{v / ktot:+.2f}", ha="center", va="bottom" if v > 0 else "top", fontsize=6.5)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(len(keys) + 1))
    ax.set_xticklabels(["Elec", "ExRep", "OrbRel", "Corr\n+Disp", "total"], fontsize=6.2)
    ax.set_xlim(-0.6, len(keys) + 0.6)
    ax.tick_params(axis="x", pad=2)
    ax.set_ylabel(r"$k_X = \mathrm{d}^2E_X/\mathrm{d}R^2$ / kcal mol⁻¹ Å⁻²")
    ax.set_ylim(-45, 100)
    panel_label(ax, "b", x=-0.18)

    # (c) decay along R; legend in the empty upper-right region created by the extended x-range
    ax = axes[2]
    mask = (r >= 2.7) & (r <= 3.3)
    series = {
        "ExRep": np.array([p["ExRep"] for p in pts]),
        "−OrbRel": np.array([-p["OrbRel"] for p in pts]),
        "−Elec": np.array([-p["Elec"] for p in pts]),
        "δ$q$ (IAO)": np.array([p["iao_ct"]["acceptor"] for p in pts]),
    }
    colors = [COL["ExRep"], COL["OrbRel"], COL["Elec"], "#ff7f0e"]
    for (name, y), c in zip(series.items(), colors):
        y0 = y / y[np.argmin(np.abs(r - 2.9))]
        b = -np.polyfit(r[mask], np.log(y[mask]), 1)[0]
        ax.semilogy(r[mask], y0[mask], "o-", color=c, ms=2.3, label=f"{name} ({1 / b:.2f} Å)")
    ax.set_xlim(2.65, 3.45)
    ax.set_ylim(0.18, 3.0)
    ax.set_xticks([2.7, 2.9, 3.1, 3.3])
    ax.set_xlabel(XLAB_R)
    ax.set_ylabel("value / value at 2.9 Å")
    ax.legend(loc="lower left", fontsize=5.6, title="decay length", title_fontsize=6, labelspacing=0.4, handletextpad=0.5)
    panel_label(ax, "c")
    fig.tight_layout(w_pad=1.0)
    fig.savefig(OUT / "fig1_dimer.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
def fig2():
    dens = load(D / "density" / "density_revpbe0-d3bj_tzvp.json")
    prof = dens["profiles"]["2.900"]
    t = np.array(prof["t"]) * 2.9
    rec = next(x for x in dens["records"] if abs(x["R_OO"] - 2.9) < 1e-6)
    xs = rec["saddle"]["distance_from_Od_A"]
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5))

    # (a) densities; legend in the empty upper-left region (y-range extended)
    ax = axes[0]
    ax.semilogy(t, prof["rho_0"], color="0.4", label=R0 + " (superposition)")
    ax.semilogy(t, prof["rho_pauli"], color=COL["ExRep"], ls="--", label=RP)
    ax.semilogy(t, prof["rho_final"], color="k", label=RF)
    ax.axvline(prof["t_Hb"] * 2.9, color="gray", ls=":", lw=0.7)
    ax.axvline(xs, color="#ff7f0e", ls=":", lw=0.9)
    ax.text(xs + 0.03, 1.2e-2, "axial min.", fontsize=6, color="#ff7f0e")
    ax.text(prof["t_Hb"] * 2.9 - 0.03, 1.3e-2, "H", fontsize=6, color="gray", ha="right")
    ax.set_xlim(0.5, 2.9)
    ax.set_ylim(5e-3, 200)
    ax.set_xlabel(XLAB_OD)
    ax.set_ylabel("ρ / e bohr⁻³")
    ax.legend(loc="upper left")
    panel_label(ax, "a")

    # (b) stage densities; legend in the empty top band (y-range extended)
    ax = axes[1]
    dp = np.array(prof["rho_pauli"]) - np.array(prof["rho_0"])
    dr = np.array(prof["rho_final"]) - np.array(prof["rho_pauli"])
    ax.plot(t, dp, color=COL["ExRep"], label=DP + " (antisymmetrization)")
    ax.plot(t, dr, color=COL["OrbRel"], label=DR + " (polarization + CT)")
    ax.plot(t, dp + dr, "k--", lw=0.8, label=r"$\Delta\rho_{\mathrm{total}}$")
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(xs, color="#ff7f0e", ls=":", lw=0.9)
    ax.text(xs + 0.02, -0.0072, "axial min.", fontsize=6, color="#ff7f0e")
    ax.set_xlim(1.2, 2.55)
    ax.set_ylim(-0.008, 0.0125)
    ax.set_yticks([-0.008, -0.004, 0, 0.004, 0.008])
    ax.set_xlabel(XLAB_OD)
    ax.set_ylabel(r"$\Delta\rho$ / e bohr⁻³")
    ax.legend(loc="upper left")
    panel_label(ax, "b")

    # (c) axial-minimum fractions vs R; system legend on the right in the empty band 0.35-0.75
    ax = axes[2]
    sets = [("dimer", "density_revpbe0-d3bj_tzvp.json", "k", "o"),
            ("trimer", "density_trimer_acid_control_revpbe0-d3bj_tzvp.json", "0.5", "s"),
            ("H₃O⁺ cluster", "density_acid_revpbe0-d3bj_tzvp.json", COL["ExRep"], "^"),
            ("OH⁻ cluster", "density_base_revpbe0-d3bj_tzvp.json", COL["Elec"], "v")]
    for lab, f, c, m in sets:
        p = D / "density" / f
        if not p.exists():
            continue
        recs = sorted(load(p)["records"], key=lambda x: x["R_OO"])
        rr = [x["R_OO"] for x in recs]
        ax.plot(rr, [x["saddle"]["fraction_promolecule"] for x in recs], m + "-", color=c, label=lab, lw=0.9, ms=3.2)
        ax.plot(rr, [x["saddle"]["fraction_relax"] for x in recs], m + "--", color=c, lw=0.7, ms=3.2, mfc="none")
        ax.plot(rr, [x["saddle"]["fraction_pauli"] for x in recs], m + ":", color=c, lw=0.7, ms=3.2, mfc="none")
    ax.axhline(0, color="k", lw=0.5)
    ax.text(2.47, 1.06, R0 + "/" + RF, fontsize=6.5, ha="left", va="bottom")
    ax.text(2.47, 0.30, DR + "/" + RF, fontsize=6.5, ha="left", va="bottom")
    ax.text(2.47, -0.28, DP + "/" + RF, fontsize=6.5, ha="left", va="top")
    ax.set_xlabel(XLAB_R)
    ax.set_ylabel("fraction of ρ at the axial minimum")
    ax.set_xlim(2.45, 3.3)
    ax.set_ylim(-0.35, 1.2)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.55), fontsize=6)
    panel_label(ax, "c")
    fig.tight_layout(w_pad=1.0)
    fig.savefig(OUT / "fig2_density.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
def fig3():
    wp4 = load(D / "analysis" / "wp4_environments.json")
    envs = {e["label"]: e for e in wp4["environments"]}
    liq = {tag: load(D / "analysis" / f"liquid_{tag}_analysis.json") for tag in ("liquid", "liquid_h3o", "liquid_oh")}
    fig = plt.figure(figsize=(7.0, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.62, wspace=0.32, top=0.90, bottom=0.12, left=0.09, right=0.98)

    # (a) gas-phase clusters: k_X bars; legend above the panel
    ax = fig.add_subplot(gs[0, 0])
    order = [("dimer", "dimer"), ("trimer acid-ctrl", "trimer"), ("acid pair", "H₃O⁺\ncontrol"),
             ("acid H3O+", "H₃O⁺\ncluster"), ("base pair", "OH⁻\ncontrol"), ("base OH-", "OH⁻\ncluster")]
    order = [(k, l) for k, l in order if k in envs]
    x = np.arange(len(order))
    w = 0.19
    for i, ch in enumerate(CH):
        ax.bar(x + (i - 1.5) * w, [envs[k]["k"][ch] for k, _ in order], w, color=COL[ch], label=NAME[ch])
    ax.plot(x, [envs[k]["k_total"] for k, _ in order], "_", color="k", ms=14, mew=1.6, label=r"$k_{\mathrm{total}}$")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in order], fontsize=6.3)
    ax.set_ylabel(r"$k_X$ / kcal mol⁻¹ Å⁻²")
    ax.set_ylim(-95, 215)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=5.8, columnspacing=1.0, handlelength=1.3)
    panel_label(ax, "a", x=-0.22)

    # (b) liquid: k(embed) vs dq(embed); legend in the empty lower-right region
    ax = fig.add_subplot(gs[0, 1])
    lab = {"liquid": ("neutral water", COL["Elec"]), "liquid_h3o": ("next to H₃O⁺", COL["ExRep"]),
           "liquid_oh": ("next to OH⁻", COL["OrbRel"])}
    ks, qs = [], []
    for tag, d in liq.items():
        rows = [r for r in d["rows"] if r.get("embed")]
        q = np.array([r["embed"]["ct_iao"] for r in rows])
        k = np.array([r["embed"]["k"] for r in rows])
        ax.plot(q, k, "o", color=lab[tag][1], label=lab[tag][0] + " (embedded pair)", ms=4, alpha=0.85)
        rows = [r for r in d["rows"] if r.get("pair")]
        ax.plot([r["pair"]["ct_iao"] for r in rows], [r["pair"]["k"] for r in rows], "s", color=lab[tag][1],
                ms=2.2, alpha=0.3, mew=0)
        ks += list(k)
        qs += list(q)
    ax.plot([], [], "s", color="0.5", ms=2.2, alpha=0.5, mew=0, label="bare pair")
    ks, qs = np.array(ks), np.array(qs)
    n, a = np.polyfit(np.log(qs), np.log(ks), 1)
    pred = a + n * np.log(qs)
    r2 = 1 - np.sum((np.log(ks) - pred) ** 2) / np.sum((np.log(ks) - np.log(ks).mean()) ** 2)
    qq = np.linspace(0.012, 0.075, 50)
    ax.plot(qq, np.exp(a) * qq ** n, "-", color="0.3", lw=1, label=f"fit: $k$ ∝ δ$q^n$, $n$ = {n:.2f}, $R^2$ = {r2:.2f}")
    ax.plot(qq, np.exp(a) * qs.mean() ** (n - 2) * qq ** 2, ":", color="0.5", lw=1,
            label="$n$ = 2 reference (normalized at mean δ$q$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.012, 0.09)
    ax.set_ylim(2.5, 60)
    ax.xaxis.set_major_locator(FixedLocator([0.015, 0.02, 0.03, 0.04, 0.06, 0.08]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_major_locator(FixedLocator([3, 5, 10, 20, 40]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("δ$q$, electrons lost by the acceptor (IAO) / e")
    ax.set_ylabel("$k$ / kcal mol⁻¹ Å⁻²")
    ax.legend(loc="lower right", fontsize=5.4, labelspacing=0.5)
    panel_label(ax, "b", x=-0.2)

    # (c) channel shares; legend above the panel, three-line tick labels
    ax = fig.add_subplot(gs[1, :])
    groups = []

    def gas_share(hi, lo, label):
        e, rr = envs[hi], envs[lo]
        dk = e["k_total"] - rr["k_total"]
        groups.append((label, [(e["k"][ch] - rr["k"][ch]) / dk for ch in CH], f"Δ$k$/$k$ = {dk / rr['k_total']:+.0%}", None))

    if "trimer acid-ctrl" in envs:
        gas_share("trimer acid-ctrl", "dimer", "gas phase\ntrimer\n− dimer")
    if "acid H3O+" in envs and "acid pair" in envs:
        gas_share("acid H3O+", "acid pair", "gas phase\nH₃O⁺ cluster\n− bare-pair control")
    if "base OH-" in envs and "base pair" in envs:
        gas_share("base OH-", "base pair", "gas phase\nOH⁻ cluster\n− bare-pair control")
    # liquid: one prespecified estimator (sum dk_X / sum dk over all bonds), 95% bootstrap percentile intervals
    LS = rs.liquid_shares()
    for tag, (l, c) in lab.items():
        v = LS[tag]
        sh = [v["shares"][ch] for ch in CH]
        err = ([v["shares"][ch] - v["ci"][ch][0] for ch in CH], [v["ci"][ch][1] - v["shares"][ch] for ch in CH])
        m, se = v["dkk"]
        groups.append((f"liquid\n{l}\nfield − bare pair", sh, f"Δ$k$/$k$ = {m:+.0%} ± {se:.0%}\n($n$ = {v['N']})", err))
    x = np.arange(len(groups))
    w = 0.19
    for i, ch in enumerate(CH):
        ax.bar(x + (i - 1.5) * w, [g[1][i] for g in groups], w, color=COL[ch], label=NAME[ch])
        xe = [xi + (i - 1.5) * w for xi, g in zip(x, groups) if g[3] is not None]
        ye = [g[1][i] for g in groups if g[3] is not None]
        ee = np.array([[g[3][0][i] for g in groups if g[3] is not None], [g[3][1][i] for g in groups if g[3] is not None]])
        ax.errorbar(xe, ye, yerr=ee, fmt="none", ecolor="k", elinewidth=0.6, capsize=1.5)
    for xi, g in zip(x, groups):
        top = max(g[1]) + (g[3][1][1] if g[3] is not None else 0.0)
        ax.text(xi, top + 0.2, g[2], ha="center", va="bottom", fontsize=6.0, linespacing=1.1)
    ax.axhline(0, color="k", lw=0.5)
    ax.axhline(1, color="k", lw=0.5, ls=":")
    ax.text(len(groups) - 0.5, 1.08, "Δ$k$", fontsize=6, ha="right", va="bottom", color="0.3")
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups], fontsize=6.2)
    ax.set_ylabel(r"share of $\Delta k$, $\Delta k_X/\Delta k$")
    ax.set_ylim(-3.0, 8.6)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=4, fontsize=6.3, columnspacing=1.5, handlelength=1.4)
    panel_label(ax, "c", x=-0.09)
    fig.savefig(OUT / "fig3_environments.png")
    plt.close(fig)
    return n, r2


# ---------------------------------------------------------------------------
def toc_plain():
    """Data-only TOC alternative (3.25 x 1.75 in) with the revised wording."""
    fc = load(D / "fc" / "fc_revpbe0-d3bj_qzvp.json")
    comp = fc["components"]
    fig, ax = plt.subplots(figsize=(3.25, 1.75))
    keys = ("Elec", "ExRep", "OrbRel", "CorrDisp")
    vals = [comp[k]["k_kcal_per_A2"] / fc["k_total_kcal_per_A2"] for k in keys]
    ax.bar(range(4), vals, color=[COL[k] for k in keys], width=0.6)
    ax.axhline(0, color="k", lw=0.5)
    ax.axhline(1, color="k", lw=0.6, ls="--")
    ax.text(3.35, 1.1, r"$k_{\mathrm{total}}$", fontsize=6.5, ha="right")
    for i, v in enumerate(vals):
        ax.text(i, v + (0.15 if v > 0 else -0.15), f"{v:+.1f}".replace("-", "−"), ha="center",
                va="bottom" if v > 0 else "top", fontsize=5.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["electrostatics", "exchange–\nrepulsion", "orbital\nrelaxation", "correlation\n+ dispersion"], fontsize=5.4)
    ax.set_ylabel(r"$k_X / k_{\mathrm{total}}$", fontsize=7)
    ax.set_ylim(-2.7, 5.4)
    ax.set_title("Local-coordinate stiffness tracks the exchange–repulsion wall;\n"
                 "orbital relaxation softens the coordinate", fontsize=6.4)
    fig.tight_layout()
    fig.savefig(OUT / "toc_plain.png", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    import sys

    wanted = set(sys.argv[1:]) or {"fig1", "fig2", "fig3", "toc_plain"}
    if "fig1" in wanted:
        fig1()
    if "fig2" in wanted:
        fig2()
    if "fig3" in wanted:
        n, r2 = fig3()
        print(f"pooled liquid exponent n = {n:.2f}, R2 = {r2:.2f}")
    if "toc_plain" in wanted:
        toc_plain()
    print(f"figures written to {OUT}: {', '.join(sorted(wanted))}")
