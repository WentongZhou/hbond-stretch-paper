"""Collect everything the Letter needs from hbond-stretch-eda into this directory.

Re-runnable: copies are overwritten, nothing in the source tree is touched.
    python collect_materials.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "hbond-stretch-eda"
REPO = Path("D:/pyscf-dm-eda")

LAYOUT = {
    "figures": [("figures", "*.png")],
    "tables": [("results", "fc_*.md"), ("results", "summary_tables.md"), ("results", "proxies_*.md"),
               ("results", "ccsdt_*.md"), ("results", "wp4_*.md"), ("results", "liquid_*_analysis.md"),
               ("results", "liquid_compare.md")],
    "docs": [("results", "PRELIMINARY_RESULTS.md"), ("results", "RESULTS_2026-09-03.md"),
             ("results", "RESULTS_WP6_liquid.md"), (".", "RESEARCH_PLAN.md"), (".", "README.md")],
    "data/fc": [("results", "fc_*.json")],
    "data/scans": [("results", "scan_*.json")],
    "data/ccsdt": [("results", "ccsdt_*.json")],
    "data/density": [("results", "density_*.json")],
    "data/optimization": [("results", "opt_*.json"), ("results", "cluster_*.json")],
    "data/analysis": [("results", "wp4_*.json"), ("results", "liquid_*_analysis.json"), ("results", "liquid_md*.json")],
    "data/liquid_clusters": [("results/liquid", "*.json")],
    "geometries": [("geometries", "*.xyz")],
    "geometries/liquid": [("geometries/liquid", "*.xyz")],
    "geometries/liquid_h3o": [("geometries/liquid_h3o", "*.xyz"), ("geometries/liquid_h3o", "*.json")],
    "geometries/liquid_oh": [("geometries/liquid_oh", "*.xyz"), ("geometries/liquid_oh", "*.json")],
    "methods/scripts": [("scripts", "*.py"), ("scripts", "*.sh"), ("scripts", "*.xml")],
    "logs": [("logs", "*.out"), ("logs", "*.log")],
}


def main() -> int:
    manifest = [f"# Manifest (collected {date.today().isoformat()} from {SRC})", ""]
    total = 0
    for target, sources in LAYOUT.items():
        dst = HERE / target
        dst.mkdir(parents=True, exist_ok=True)
        names = []
        for sub, pattern in sources:
            for f in sorted((SRC / sub).glob(pattern)):
                if f.is_file():
                    shutil.copy2(f, dst / f.name)
                    names.append(f.name)
        total += len(names)
        manifest.append(f"## {target}/  ({len(names)} files)")
        manifest.append("")
        manifest.extend(f"- {n}" for n in names)
        manifest.append("")
    # package patch and environment
    meth = HERE / "methods"
    meth.mkdir(exist_ok=True)
    try:
        diff = subprocess.run(["git", "-C", str(REPO), "diff", "--", "src/pyscf_dm_eda/eda.py", "tests/test_eda.py"],
                              capture_output=True, text=True, check=True).stdout
        (meth / "pyscf-dm-eda_point_charges.patch").write_text(diff, encoding="utf-8")
        head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        manifest.append(f"## methods/pyscf-dm-eda_point_charges.patch\n\nUncommitted diff of pyscf-dm-eda (HEAD {head}): SCFConfig.point_charges embedding + tests.\n")
    except Exception as exc:  # git may be unavailable on the Windows side
        manifest.append(f"## methods/pyscf-dm-eda_point_charges.patch\n\n(not exported: {exc!r}; run `git -C D:/pyscf-dm-eda diff` in WSL)\n")
    shutil.copy2(REPO / "src/pyscf_dm_eda/eda.py", meth / "eda.py") if (REPO / "src/pyscf_dm_eda/eda.py").exists() else None
    (HERE / "MANIFEST.md").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(f"copied {total} files; manifest at {HERE / 'MANIFEST.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
