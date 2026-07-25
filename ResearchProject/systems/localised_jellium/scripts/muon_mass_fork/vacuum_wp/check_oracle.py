#!/usr/bin/env python3
"""Phase-2 analytic oracle checker for a vacuum free-particle WP run.

Validates the per-state mass fork against the EXACT free-Gaussian spreading law
(hbar = 1):

    sigma_rho_z(t)^2 = sigma_rho0^2 + b * t^2 ,   b = 1 / (4 m^2 sigma_rho0^2)
    sigma_rho0 = sigma_WP / sqrt(2)

so a linear fit of sigma_z2 vs t^2 gives intercept sigma_rho0^2 and slope b, from
which m_fit = 1 / (2 sigma_rho0 sqrt(b)). Also checks:
  - t=0 density variance  sigma_z2(0) == sigma_WP^2/2         (WP injection sanity)
  - group velocity        d<z>/dt == k0 / m                   (if k0 > 0)
  - KE conservation       <T> drift < tol                     (free particle)

Usage:
  check_oracle.py <run_dir> --sigma_wp 0.5 --mass 1.0 [--k0 0.0] [--json out.json]

Exit code 0 iff all applicable oracles pass. Called by orchestrate.py phase2.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np


def _read_csv(path: Path) -> dict[str, np.ndarray]:
    # inqkit observable CSVs may lead with a '#' provenance comment line before
    # the real header row (step,time_au,...); skip any leading '#' lines.
    lines = [ln for ln in path.read_text().strip().splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    hdr = [h.strip() for h in lines[0].split(",")]
    rows = np.array([[float(x) for x in ln.split(",")] for ln in lines[1:]])
    return {h: rows[:, i] for i, h in enumerate(hdr)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--sigma_wp", type=float, required=True)
    ap.add_argument("--mass", type=float, required=True)
    ap.add_argument("--k0", type=float, default=0.0)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    run = Path(a.run_dir)
    rs_csv = run / "raw/observables/wp_real_space_stats.csv"
    ms_csv = run / "raw/observables/wp_momentum_stats.csv"
    if not rs_csv.exists():
        print(f"FAIL: missing {rs_csv}"); return 2

    rs = _read_csv(rs_csv)
    t = rs["time_au"]
    sig_z2 = rs["sigma_z2"]
    z_mean = rs["z_mean"]
    norm = rs.get("norm_check", np.ones_like(t))

    sigma_rho0 = a.sigma_wp / np.sqrt(2.0)
    sig_rho0_sq = sigma_rho0 ** 2                       # expected intercept = sigma_WP^2/2
    b_expected = 1.0 / (4.0 * a.mass**2 * sig_rho0_sq)  # expected slope vs t^2

    results = {}
    ok = True

    # --- (1) t=0 density variance -----------------------------------------
    v0 = float(sig_z2[0])
    rel0 = abs(v0 - sig_rho0_sq) / sig_rho0_sq
    p0 = rel0 < 5e-3
    ok &= p0
    results["sigma_z2(0)"] = dict(got=v0, expect=sig_rho0_sq, rel=rel0, pass_=bool(p0), tol=5e-3)

    # --- (2) parabola fit over the clean (in-box) window ------------------
    # keep points while the packet is well inside the box: use the recorded
    # variance itself (4*sigma_rho < L/2 proxy) AND norm still ~1.
    good = (norm > 0.99) & (norm < 1.01)
    # drop points once variance stops being monotonic (boundary wrap) --
    # keep the leading monotonic-increasing stretch.
    mono = np.ones_like(sig_z2, dtype=bool)
    for i in range(1, len(sig_z2)):
        if sig_z2[i] < sig_z2[i-1] - 1e-9:
            mono[i:] = False; break
    mask = good & mono
    n_fit = int(mask.sum())
    if n_fit >= 5:
        A = np.vstack([np.ones(n_fit), (t[mask])**2]).T
        coef, *_ = np.linalg.lstsq(A, sig_z2[mask], rcond=None)
        intercept, slope = float(coef[0]), float(coef[1])
        slope = max(slope, 1e-300)
        m_fit = 1.0 / (2.0 * sigma_rho0 * np.sqrt(slope))
        rel_m = abs(m_fit - a.mass) / a.mass
        # tolerance: 5% for a well-spreading packet; a barely-spreading muon
        # over a short window has a tiny, noise-dominated slope -> relax to 25%.
        spread = (sig_z2[mask][-1] - v0) / v0
        tol_m = 0.05 if spread > 0.05 else 0.25
        pm = rel_m < tol_m
        ok &= pm
        results["mass_fit"] = dict(m_fit=m_fit, m_true=a.mass, rel=rel_m,
                                   slope=slope, slope_expect=b_expected,
                                   intercept=intercept, n_fit=n_fit,
                                   spread_frac=spread, tol=tol_m, pass_=bool(pm))
    else:
        ok = False
        results["mass_fit"] = dict(error=f"only {n_fit} clean points", pass_=False)

    # --- (3) group velocity  d<z>/dt = k0/m -------------------------------
    if a.k0 > 0.0 and n_fit >= 5:
        v_expect = a.k0 / a.mass
        A = np.vstack([np.ones(n_fit), t[mask]]).T
        coef, *_ = np.linalg.lstsq(A, z_mean[mask], rcond=None)
        v_fit = float(coef[1])
        rel_v = abs(v_fit - v_expect) / max(abs(v_expect), 1e-12)
        pv = rel_v < 0.02
        ok &= pv
        results["v_group"] = dict(v_fit=v_fit, v_expect=v_expect, rel=rel_v,
                                  tol=0.02, pass_=bool(pv))

    # --- (4) norm + KE conservation ---------------------------------------
    nd = float(np.max(np.abs(norm - 1.0)))
    pn = nd < 1e-3
    ok &= pn
    results["norm_drift"] = dict(max_abs=nd, tol=1e-3, pass_=bool(pn))

    # KE oracle from the WP's OWN momentum stats (<k^2> is exactly conserved for
    # a free particle, and is mass-independent so it isolates the dynamics from
    # the spectator electron that the total-energy column would include).
    if ms_csv.exists():
        ms = _read_csv(ms_csv)
        if all(c in ms for c in ("px2_mean", "py2_mean", "pz2_mean")):
            k2 = ms["px2_mean"] + ms["py2_mean"] + ms["pz2_mean"]
            k20 = k2[0] if abs(k2[0]) > 1e-12 else 1.0
            drift = float(np.max(np.abs(k2 - k2[0])) / abs(k20))
            pke = drift < 5e-3
            ok &= pke
            results["k2_drift"] = dict(max_rel=drift, k2_0=float(k2[0]),
                                       tol=5e-3, pass_=bool(pke))
        # absolute momentum check: <k_z> should equal k0
        if "pz_mean" in ms and a.k0 > 0.0:
            pz0 = float(ms["pz_mean"][0])
            rel = abs(pz0 - a.k0) / a.k0
            pkz = rel < 0.05
            ok &= pkz
            results["kz_mean(0)"] = dict(got=pz0, expect=a.k0, rel=rel,
                                         tol=0.05, pass_=bool(pkz))

    verdict = "PASS" if ok else "FAIL"
    print(f"=== vacuum oracle {verdict} : {run.name} "
          f"(sigma_WP={a.sigma_wp}, mass={a.mass}, k0={a.k0}) ===")
    for k, v in results.items():
        mark = "ok " if v.get("pass_") else "XX "
        print(f"  [{mark}] {k}: {v}")
    if a.json:
        Path(a.json).write_text(json.dumps(dict(verdict=verdict, results=results), indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
