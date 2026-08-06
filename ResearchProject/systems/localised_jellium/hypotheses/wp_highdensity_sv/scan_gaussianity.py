"""
scan_gaussianity — read the launch-distance injection scan and report, for every
trial, (a) how much weight orthogonalisation removed and (b) how Gaussian the
surviving packet still is in momentum and in real space.

CONTEXT. The effective-sigma hypothesis (docs/plans/effective-sigma-near-launch.md)
says the wavepacket's ARRIVAL width, not its launch sigma, sets the interaction
with the slab. Testing it means launching close to the slab — inside the
electronic spill-out, where the packet overlaps the occupied bath and inqkit's
injector projects that component out. The scan finds the smallest standoff at
which that projection removes < 3 % of the packet (user criterion, 2026-08-01).

WHAT DECIDES WHAT. `removed_weight` alone decides accept/reject — that gate lives
in the C++ scan program and is already applied by the time this script runs.
Gaussianity is REPORTED here and never vetoes (user decision, 2026-08-01). This
script exists so the reported half is quantitative rather than eyeballed.

THE RIGHT GAUSSIANITY OBSERVABLE. For a Gaussian wavepacket the k_z MARGINAL
P(k_z) = sum_{kx,ky} |psi~(k)|^2 is exactly N(k0, sigma_p^2) with
sigma_p = 1/(sqrt2 sigma). The RADIAL distribution n(|k|) that inqkit's
MomentumDistribution computes is NOT Gaussian for a drifting packet, so it cannot
answer this question. See inqview.visualisation.field_io.kz_marginal.

VALIDATION OF THIS PIPELINE (2026-08-01, on the far-launch v2p0 t=0 packet, which
is undeformed by construction at 11.5 Bohr standoff): recovered <k_z> = 1.999791
vs 2.0 (-0.01 %), sigma_kz = 1.414473 vs 1.414214 (+0.018 %), skewness -0.005,
excess kurtosis +0.022, R^2 vs the ANALYTIC Gaussian = 1.000000. The extraction is
therefore known-good before it is applied to any near-launch trial.

Usage:
    python scan_gaussianity.py            # analyse every trial, write CSV + figure
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "inq-stack" / "python"))

from inqview.visualisation.field_io import (  # noqa: E402
    load_complex_vti, kz_marginal, gaussian_fit_quality)

SCAN_DIR = (REPO / "ResearchProject/systems/localised_jellium/scripts"
            / "wp_highdensity_sv/inject_scan/results/scan")
OUT_DIR = Path(__file__).resolve().parent
SLAB_FACE_Z = -12.5      # Cfg::SLAB_CENTER_BOHR - Cfg::SLAB_HALF_WIDTH
EDGE_WIDTH = 1.0         # Cfg::EDGE_WIDTH_BOHR (erfc softening)
CRITERION_PC = 3.0       # user criterion


_KV = re.compile(r"(\w+)\s*=\s*([^\s]+)")


def read_trial(path: Path) -> dict:
    """Parse a trial.txt into {key: float|str}.

    These files (like run_summary.txt across the project) put SEVERAL
    `key = value` pairs on one line, e.g.
        wp_k0 = 2  launch_z = -14  standoff_bohr = 1.5
    so a naive split on the first '=' swallows the rest of the line into the
    first value. Match every pair instead.
    """
    out: dict = {}
    for line in path.read_text().splitlines():
        for k, v in _KV.findall(line):
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def analyse_trial(trial_dir: Path) -> dict | None:
    tfile = trial_dir / "trial.txt"
    if not tfile.exists():
        return None
    rec = read_trial(tfile)

    sigma = rec["wp_sigma_bohr"]
    k0 = rec["wp_k0"]
    sigma_p = 1.0 / (np.sqrt(2.0) * sigma)

    row = {
        "trial": trial_dir.name,
        "launch_z": rec["launch_z"],
        "standoff_bohr": rec["standoff_bohr"],
        "k0": k0,
        "sigma_wp": sigma,
        "removed_percent": rec["removed_percent"],
        "sum_overlap_sq": rec["sum_overlap_sq"],
        "max_overlap": rec["max_overlap"],
        "closure_residual_rel": rec["closure_residual_rel"],
        "accept": rec.get("accept") == "true",
        # C++-side moments (node convention, the authoritative ones)
        "centroid_z_cpp": rec["centroid_z"],
        "density_std_z_cpp": rec["density_std_z"],
        "mean_pz_cpp": rec["mean_pz"],
        "sigma_pz2_cpp": rec["sigma_pz2"],
    }

    wf = trial_dir / "raw/vti/wavefunction_wp/wavefunction_t000000.vti"
    if wf.exists():
        field = load_complex_vti(wf)
        kz, p = kz_marginal(field)
        row.update(gaussian_fit_quality(kz, p, k0=k0, sigma_p=sigma_p))
        row["_kz"] = kz
        row["_p"] = p
        # real-space |psi|^2 profile along z (planar-summed)
        n = np.abs(field.data) ** 2
        prof = n.sum(axis=(0, 1))
        prof = prof / prof.sum()
        row["_z"] = field.z
        row["_nz"] = prof
    return row


def main() -> int:
    if not SCAN_DIR.exists():
        print(f"no scan directory yet: {SCAN_DIR}")
        return 1

    rows = []
    for d in sorted(SCAN_DIR.iterdir()):
        if not d.is_dir():
            continue
        r = analyse_trial(d)
        if r is not None:
            rows.append(r)
    if not rows:
        print("no completed trials found")
        return 1

    # Production trials only for the summary table; the z=-24 regression trial is
    # a build check, not a scan point, but IS kept for the figure as the
    # undeformed reference.
    df = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")}
                       for r in rows]).sort_values("standoff_bohr")

    csv = OUT_DIR / "scan_gaussianity.csv"
    df.to_csv(csv, index=False)

    pd.set_option("display.width", 200)
    cols = ["trial", "launch_z", "standoff_bohr", "removed_percent", "accept",
            "mean_kz", "std_kz", "skewness", "excess_kurtosis", "r2_analytic"]
    cols = [c for c in cols if c in df.columns]
    print(f"\n=== injection scan: {len(df)} trials  (criterion removed < {CRITERION_PC} %) ===")
    print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.5g}"))
    print(f"\nwrote {csv}")

    _figure(rows)
    return 0


def _figure(rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from inqview.visualisation.style import apply_theme
    apply_theme()

    scan = sorted([r for r in rows if not r["trial"].startswith("regression")],
                  key=lambda r: r["standoff_bohr"])
    ref = [r for r in rows if r["trial"].startswith("regression")]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))

    # (a) the decision: removed weight vs standoff
    ax = axes[0]
    if scan:
        s = [r["standoff_bohr"] for r in scan]
        w = [r["removed_percent"] for r in scan]
        ax.plot(s, w, "o-", label="scan trials")
    for r in ref:
        ax.plot(r["standoff_bohr"], r["removed_percent"], "s",
                color="0.4", label="far launch (z=-24)")
    ax.axhline(CRITERION_PC, ls="--", color="crimson", lw=1.2)
    ax.text(0.98, CRITERION_PC, f" {CRITERION_PC:g} % criterion", color="crimson",
            va="bottom", ha="right", transform=ax.get_yaxis_transform(), fontsize=8)
    ax.axvspan(0, EDGE_WIDTH, color="0.85", zorder=0)
    ax.text(EDGE_WIDTH / 2, ax.get_ylim()[1], "erfc\nsoftening", ha="center",
            va="top", fontsize=7, color="0.35")
    ax.set_yscale("log")
    ax.set_xlabel("standoff from slab face (Bohr)")
    ax.set_ylabel("weight removed by orthogonalisation (%)")
    ax.set_title("(a) the accept/reject decision")
    ax.legend(frameon=False, fontsize=8)

    # (b) k_z marginal vs the ANALYTIC Gaussian
    ax = axes[1]
    for r in scan + ref:
        if "_kz" not in r:
            continue
        ax.plot(r["_kz"], r["_p"], lw=1.2,
                label=f"z={r['launch_z']:.1f} ({r['removed_percent']:.2f} %)")
    if scan or ref:
        r0 = (scan + ref)[0]
        kz = r0["_kz"]
        sp = 1.0 / (np.sqrt(2.0) * r0["sigma_wp"])
        analytic = np.exp(-0.5 * ((kz - r0["k0"]) / sp) ** 2) / (sp * np.sqrt(2 * np.pi))
        ax.plot(kz, analytic, "k--", lw=1.0, label="analytic N(k$_0$, $\\sigma_p^2$)")
        ax.set_xlim(r0["k0"] - 5 * sp, r0["k0"] + 5 * sp)
    ax.set_xlabel("$k_z$ (Bohr$^{-1}$)")
    ax.set_ylabel("$P(k_z)$")
    ax.set_title("(b) momentum profile (reported, not a veto)")
    ax.legend(frameon=False, fontsize=7)

    # (c) real-space packet vs the slab
    ax = axes[2]
    for r in scan + ref:
        if "_z" not in r:
            continue
        ax.plot(r["_z"], r["_nz"], lw=1.2, label=f"z={r['launch_z']:.1f}")
    ax.axvline(SLAB_FACE_Z, color="crimson", ls="--", lw=1.2)
    ax.axvspan(SLAB_FACE_Z, 0, color="0.85", zorder=0)
    ax.text(SLAB_FACE_Z, ax.get_ylim()[1], " slab", color="crimson",
            va="top", ha="left", fontsize=8)
    ax.set_xlim(-26, -8)
    ax.set_xlabel("z (Bohr)")
    ax.set_ylabel("normalised $|\\psi(z)|^2$")
    ax.set_title("(c) where the packet starts")
    ax.legend(frameon=False, fontsize=7)

    fig.tight_layout()
    out = OUT_DIR / "scan_gaussianity.png"
    fig.savefig(out, dpi=160)
    print(f"wrote {out}")


if __name__ == "__main__":
    raise SystemExit(main())
