"""Fit statistics used by the revised manuscript and Supporting Information.

Every number is recomputed from data/ so that the text, Table S4 and Table S5 stay consistent.
    python revision_stats.py      prints the summary
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
D = HERE / "data"
SCRATCH = Path(r"C:\Users\HUAWEI\AppData\Local\Temp\claude\C--Users-HUAWEI-Desktop-hbond-stretch-paper\905032d0-211c-4841-9551-d6940fa6a092\scratchpad")
CH = ("Elec", "ExRep", "OrbRel", "CorrDisp")


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def loglog_fit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    res = stats.linregress(np.log(x), np.log(y))
    n = len(x)
    t = stats.t.ppf(0.975, n - 2)
    return {"n": res.slope, "se": res.stderr, "lo": res.slope - t * res.stderr, "hi": res.slope + t * res.stderr,
            "r2": res.rvalue ** 2, "N": n, "a": res.intercept}


def linfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    res = stats.linregress(x, y)
    t = stats.t.ppf(0.975, len(x) - 2)
    return {"slope": res.slope, "se": res.stderr, "lo": res.slope - t * res.stderr, "hi": res.slope + t * res.stderr,
            "r2": res.rvalue ** 2, "N": len(x)}


def decay_constants(tag):
    """b (ln y = a - bR) of the proxies from tables/proxies_<tag>.md."""
    out = {}
    for line in (HERE / "tables" / f"proxies_{tag}.md").read_text(encoding="utf-8").splitlines():
        if line.startswith("| ") and "---" not in line and "quantity" not in line:
            c = [x.strip() for x in line.strip("|").split("|")]
            out[c[0]] = {"b": float(c[1]), "length": float(c[2]), "r2": float(c[3]), "N": int(c[4])}
    return out


def liquid_rows():
    return {t: load(D / "analysis" / f"liquid_{t}_analysis.json")["rows"] for t in ("liquid", "liquid_h3o", "liquid_oh")}


def liquid_fits():
    liq = liquid_rows()
    out = {}
    q, k, q0, k0, qp, kp = [], [], [], [], [], []
    for t, rows in liq.items():
        e = [r["embed"] for r in rows if r.get("embed")]
        p = [r["pair"] for r in rows if r.get("pair")]
        qq = [x["ct_iao"] for x in e]
        kk = [x["k"] for x in e]
        qq0 = [x["at_R0"]["charge_transfer_acceptor"]["iao"]["value_e"] for x in e]
        kk0 = [x["at_R0"]["k_total_kcal_per_A2"] for x in e]
        out[f"{t}_embed_min"] = loglog_fit(qq, kk)
        out[f"{t}_embed_R0"] = loglog_fit(qq0, kk0)
        q += qq
        k += kk
        q0 += qq0
        k0 += kk0
        qp += [x["ct_iao"] for x in p]
        kp += [x["k"] for x in p]
    out["pooled_embed_min"] = loglog_fit(q, k)
    out["pooled_pair_min"] = loglog_fit(qp, kp)
    out["pooled_embed_R0"] = loglog_fit(q0, k0)
    return out


def gas_fits():
    envs = load(D / "analysis" / "wp4_environments.json")["environments"]
    k = [e["k_total"] for e in envs]
    return {"IAO": loglog_fit([e["ct_iao"] for e in envs], k), "Mulliken": loglog_fit([e["ct_mull"] for e in envs], k)}


def across_method():
    meth = [("hf_tzvp", "HF"), ("revpbe0-d3bj_tzvp", "revPBE0-D3(BJ)"), ("pbe0-d3bj_tzvp", "PBE0-D3(BJ)"),
            ("r2scan-d4_tzvp", "r²SCAN-D4"), ("scan_tzvp", "SCAN")]
    pts = []
    for tag, name in meth:
        f = load(D / "fc" / f"fc_{tag}.json")
        pts.append((name, f["R_min_A"], f["k_total_kcal_per_A2"]))
    fit = linfit([p[1] for p in pts], np.log([p[2] for p in pts]))
    b = decay_constants("revpbe0-d3bj_tzvp")
    proxies = {"axial density minimum": b["rho_final at saddle"]["b"], "IAO charge": b["CT IAO (acceptor -> donor)"]["b"],
               "Mulliken charge": b["CT Mulliken (acceptor -> donor)"]["b"]}
    n = {}
    for name, bb in proxies.items():
        n[name] = {"n": -fit["slope"] / bb, "lo": -fit["hi"] / bb, "hi": -fit["lo"] / bb, "b": bb}
    return {"points": pts, "fit": fit, "n": n}


def single_pes(tag, window=0.16, grid=np.arange(2.80, 3.101, 0.05)):
    """Local curvature k(R) = d2E/dR2 from a quartic fit centred at each R of the same rigid dimer scan,
    against the proxies evaluated at the same R (interpolated on the scan grid)."""
    scan = load(D / "scans" / f"scan_{tag}.json")
    pts = sorted(scan["points"], key=lambda p: p["R_OO"])
    r = np.array([p["R_OO"] for p in pts])
    E = np.array([p["Total"] for p in pts])
    qi = np.array([p["iao_ct"]["acceptor"] for p in pts])
    qm = np.array([p["mulliken_ct"]["acceptor"] for p in pts])
    rows = []
    for Rc in grid:
        m = np.abs(r - Rc) <= window + 1e-3
        c = np.polyfit(r[m] - Rc, E[m], 4)
        rows.append({"R": float(Rc), "npts": int(m.sum()), "k": 2 * c[-3], "iao": float(np.interp(Rc, r, qi)),
                     "mull": float(np.interp(Rc, r, qm))})
    dens = D / "density" / f"density_{tag}.json"
    if dens.exists():
        recs = {round(x["R_OO"], 3): x["saddle"]["rho_final"] for x in load(dens)["records"]}
        for row in rows:
            row["rho"] = recs.get(round(row["R"], 3))
    K = np.array([x["k"] for x in rows])
    out = {"rows": rows, "dlnk_dR": linfit([x["R"] for x in rows], np.log(K)),
           "IAO": loglog_fit([x["iao"] for x in rows], K), "Mulliken": loglog_fit([x["mull"] for x in rows], K)}
    rho_rows = [x for x in rows if x.get("rho")]
    if len(rho_rows) >= 3:
        out["rho"] = loglog_fit([x["rho"] for x in rho_rows], [x["k"] for x in rho_rows])
    return out


def block_bootstrap(n_boot=5000, seed=0):
    liq = liquid_rows()
    rng = np.random.default_rng(seed)
    rows = [r for r in liq["liquid"] if r.get("embed") and r.get("pair")]
    frames = {}
    for r in rows:
        frames.setdefault(r["name"].split("_")[1], []).append(r)
    frames = list(frames.values())
    dkk = np.array([(r["embed"]["k"] - r["pair"]["k"]) / r["pair"]["k"] for r in rows])
    wmin = np.array([r["embed"]["omega"] for r in rows])
    b_dkk, b_w, b_n = [], [], []
    for _ in range(n_boot):
        sel = [r for i in rng.integers(0, len(frames), len(frames)) for r in frames[i]]
        b_dkk.append(np.mean([(r["embed"]["k"] - r["pair"]["k"]) / r["pair"]["k"] for r in sel]))
        b_w.append(np.mean([r["embed"]["omega"] for r in sel]))
        b_n.append(np.polyfit(np.log([r["embed"]["ct_iao"] for r in sel]), np.log([r["embed"]["k"] for r in sel]), 1)[0])
    allframes = list(frames) + [[r] for t in ("liquid_h3o", "liquid_oh") for r in liq[t] if r.get("embed")]
    b_pool = []
    for _ in range(n_boot):
        sel = [r for i in rng.integers(0, len(allframes), len(allframes)) for r in allframes[i]]
        b_pool.append(np.polyfit(np.log([r["embed"]["ct_iao"] for r in sel]), np.log([r["embed"]["k"] for r in sel]), 1)[0])
    return {
        "n_frames": len(frames), "n_bonds": len(rows),
        "dkk_mean": dkk.mean(), "dkk_se_naive": dkk.std(ddof=1) / np.sqrt(len(dkk)), "dkk_se_block": float(np.std(b_dkk)),
        "w_mean": wmin.mean(), "w_se_naive": wmin.std(ddof=1) / np.sqrt(len(wmin)), "w_se_block": float(np.std(b_w)),
        "n_neutral_ci": tuple(np.percentile(b_n, [2.5, 97.5])), "n_pooled_ci": tuple(np.percentile(b_pool, [2.5, 97.5])),
        "n_pooled_se": float(np.std(b_pool)),
    }


def liquid_env_summary():
    liq = liquid_rows()
    out = {}
    for t, rows in liq.items():
        allr = [r for r in rows if r.get("embed") and r.get("pair")]
        rr = [r for r in allr if abs(r["embed"]["k"] - r["pair"]["k"]) > 0.5]
        se = lambda v: float(np.std(v, ddof=1) / np.sqrt(len(v)))
        dk = np.array([(r["embed"]["k"] - r["pair"]["k"]) / r["pair"]["k"] for r in allr])
        out[t] = {
            "N": len(allr),
            "shares": {ch: float(np.mean([(r["embed"]["kX"][ch] - r["pair"]["kX"][ch]) / (r["embed"]["k"] - r["pair"]["k"]) for r in rr])) for ch in CH},
            "dkk": (dk.mean(), se(dk)),
            "dR": float(np.mean([r["embed"]["R_min"] - r["pair"]["R_min"] for r in allr])),
            "w_min": (np.mean([r["embed"]["omega"] for r in allr]), se(np.array([r["embed"]["omega"] for r in allr]))),
            "w_R0": (np.mean([r["embed"]["at_R0"]["omega_H2O_cm-1"] for r in allr]), se(np.array([r["embed"]["at_R0"]["omega_H2O_cm-1"] for r in allr]))),
            "k_R0_embed": float(np.mean([r["embed"]["at_R0"]["k_total_kcal_per_A2"] for r in allr])),
            "k_R0_pair": float(np.mean([r["pair"]["at_R0"]["k_total_kcal_per_A2"] for r in allr])),
            "R0_minus_Rmin": float(np.mean([r["embed"]["at_R0"]["R0_A"] - r["embed"]["R_min"] for r in allr])),
        }
    return out


def max_closure():
    return max(load(f)["max_abs_closure_hartree"] for f in (D / "fc").glob("fc_*.json"))


SOURCE_REPO = Path(r"C:\Users\HUAWEI\Desktop\hbond-stretch-eda")


def cutoff_test():
    """Results of the 12 Å embedding-cutoff test.  The archived copy data/analysis/cutoff_test_12A.json
    (written 2026-09-04, includes the charge-site counts) is preferred; the scratchpad JSONs of the
    original run are the fallback.  Charge-site counts are taken from the cluster specifications
    (8 Å: data/liquid_clusters; 12 Å: source repo) when missing."""
    archived = D / "analysis" / "cutoff_test_12A.json"
    if archived.exists():
        return load(archived)
    out = {}
    for f in sorted(SCRATCH.glob("cutoff_*.json")):
        out.update(load(f))
    for name, v in out.items():
        for lab, path in (("cut8", D / "liquid_clusters" / f"cluster_{name}.json"),
                          ("cut12", SOURCE_REPO / "results" / "liquid" / f"cluster_{name}_cut12.json")):
            if lab in v and path.exists():
                emb = load(path).get("embedding", {})
                v[lab]["n_sites"] = emb.get("n_sites")
                v[lab]["n_molecules"] = emb.get("n_molecules")
    return out


def decay_lengths(tag="revpbe0-d3bj_qzvp"):
    """Decay lengths 1/b (Å) of the quantities plotted in Figure 1c, from tables/proxies_<tag>.md."""
    b = decay_constants(tag)
    return {"ExRep": b["ExRep (exchange + Pauli repulsion)"]["length"], "Rep": b["Rep (Pauli repulsion only)"]["length"],
            "OrbRel": b["-OrbRel (orbital relaxation)"]["length"], "Elec": b["-Elec (electrostatics)"]["length"],
            "CorrDisp": b["-CorrDisp"]["length"], "IAO": b["CT IAO (acceptor -> donor)"]["length"]}


def gas_shares():
    """Channel shares dk_X/dk of the three gas-phase perturbations of Figure 3c (Table S7)."""
    envs = {e["label"]: e for e in load(D / "analysis" / "wp4_environments.json")["environments"]}
    out = {}
    for key, hi, lo in (("trimer", "trimer acid-ctrl", "dimer"), ("acid", "acid H3O+", "acid pair"), ("base", "base OH-", "base pair")):
        e, r = envs[hi], envs[lo]
        dk = e["k_total"] - r["k_total"]
        out[key] = {"dkk": dk / r["k_total"], "shares": {ch: (e["k"][ch] - r["k"][ch]) / dk for ch in CH}}
    return out


def liquid_shares(n_boot=5000, seed=0):
    """Figure 3c liquid estimator: shares sum(dk_X)/sum(dk) over ALL bonds of a sample (dk = k_embed - k_bare at the
    respective minima), 95% percentile intervals from a bootstrap over bonds, the mean per-bond dk/k_bare with its
    standard error over the same bonds, and (for the sensitivity note) the mean of per-bond ratios after excluding
    |dk| <= 0.5 kcal/mol/A^2, which was the construction used in the previous version of Figure 3c."""
    rng = np.random.default_rng(seed)
    out = {}
    for t, rows in liquid_rows().items():
        rr = [r for r in rows if r.get("embed") and r.get("pair")]
        dk = np.array([r["embed"]["k"] - r["pair"]["k"] for r in rr])
        kb = np.array([r["pair"]["k"] for r in rr])
        dkX = {ch: np.array([r["embed"]["kX"][ch] - r["pair"]["kX"][ch] for r in rr]) for ch in CH}
        boots = {ch: [] for ch in CH}
        for _ in range(n_boot):
            i = rng.integers(0, len(rr), len(rr))
            s = dk[i].sum()
            for ch in CH:
                boots[ch].append(dkX[ch][i].sum() / s)
        keep = np.abs(dk) > 0.5
        out[t] = {
            "N": len(rr), "N_filtered": int(keep.sum()),
            "shares": {ch: float(dkX[ch].sum() / dk.sum()) for ch in CH},
            "ci": {ch: tuple(float(x) for x in np.percentile(boots[ch], [2.5, 97.5])) for ch in CH},
            "shares_filtered_mean_of_ratios": {ch: float(np.mean(dkX[ch][keep] / dk[keep])) for ch in CH},
            "dkk": (float(np.mean(dk / kb)), float(np.std(dk / kb, ddof=1) / np.sqrt(len(rr)))),
        }
    return out


def oh_rows():
    """Donor O-H elongation/shortening scans available in data/fc (Table S6), sorted by the O-H change."""
    rows = []
    for f in sorted((D / "fc").glob("fc_revpbe0-d3bj_tzvp_oh*.json")) + [D / "fc" / "fc_revpbe0-d3bj_tzvp.json"]:
        d = load(f)
        rows.append((d.get("oh_elongation_angstrom", 0.0), f.stem[3:]))
    return [tag for _, tag in sorted(rows)]


def force_balance_summary():
    """Largest |sum_X dE_X/dR| at R_min over the fitted scans that have an interior minimum, and the scans that
    have none (their fc JSON carries minimum_found = False; the reported values are then taken at R0)."""
    ok, no_minimum = [], []
    for f in sorted((D / "fc").glob("fc_*.json")):
        d = load(f)
        if d.get("minimum_found", True):
            ok.append((abs(d["force_balance_sum_of_slopes"]), f.stem[3:]))
        else:
            no_minimum.append(f.stem[3:])
    return {"max_ok": max(ok)[0], "worst_ok": max(ok)[1], "n_ok": len(ok), "n_total": len(ok) + len(no_minimum),
            "no_minimum": no_minimum}


def fmt_n(f, nd=2):
    return f"{f['n']:.{nd}f} ({f['lo']:.{nd}f}–{f['hi']:.{nd}f})"


if __name__ == "__main__":
    lf = liquid_fits()
    for k, v in lf.items():
        print(f"{k:24s} n = {fmt_n(v)}  se {v['se']:.3f}  R² {v['r2']:.3f}  N {v['N']}")
    for k, v in gas_fits().items():
        print(f"gas {k:20s} n = {fmt_n(v)}  R² {v['r2']:.3f}")
    am = across_method()
    print("across-method d ln k/dR = %.2f (%.2f–%.2f), R² %.3f" % (am["fit"]["slope"], am["fit"]["lo"], am["fit"]["hi"], am["fit"]["r2"]))
    for k, v in am["n"].items():
        print(f"   n {k}: {v['n']:.2f} ({v['lo']:.2f}–{v['hi']:.2f}), b = {v['b']}")
    for tag in ("revpbe0-d3bj_tzvp", "revpbe0-d3bj_qzvp"):
        sp = single_pes(tag)
        print(tag, "d ln k_loc/dR = %.2f ± %.2f" % (sp["dlnk_dR"]["slope"], sp["dlnk_dR"]["se"]),
              "IAO", fmt_n(sp["IAO"]), "Mull", fmt_n(sp["Mulliken"]), "rho", fmt_n(sp["rho"]) if "rho" in sp else "—")
    bb = block_bootstrap()
    print(bb)
    for t, v in liquid_env_summary().items():
        print(t, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in v.items()})
    print("max closure", max_closure())
    print("cutoff test", cutoff_test())
