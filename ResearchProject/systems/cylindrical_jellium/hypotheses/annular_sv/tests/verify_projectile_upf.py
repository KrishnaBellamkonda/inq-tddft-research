#!/usr/bin/env python3
"""Phase-1a(ii) gate: VERIFY the classical-electron Gaussian projectile UPF.

The cylindrical-jellium campaign locks a classical electron projectile carried by
the fictitious "H" species with an erf-smoothed GAUSSIAN local potential of width
sigma_pot = sigma_WP/sqrt(2) = 0.5/sqrt(2) ~ 0.354 Bohr (the sqrt(2) rule).

The candidate electron_gaussian_*.upf files carry STALE is_coulomb="T" headers
(and one a stale "bare Coulomb" comment) — they must NOT be trusted by header.
This script verifies by DATA: it parses PP_R + PP_LOCAL, confirms

  (1) a FINITE, REPULSIVE core  V(0) > 0  (no -Z/r singularity),
  (2) the tabulated potential fits the erf-Gaussian form  V(r) = Z * erf(r/(sqrt(2) sigma_pot)) / r
      with Z ~ 1 (unit charge magnitude) and sigma_pot ~ 0.354 Bohr,
  (3) monotone decay to ~0 at large r.

INQ uses the PP_LOCAL table verbatim (upf2.hpp:235 reads it, *0.5 Ry->Ha; the
is_coulomb flag is referenced nowhere in inq/src, pseudopod, or upf2.hpp), so a
clean Gaussian PP_LOCAL is exactly what the projectile feels. Verified 2026-06-28.

Run:  venv/bin/python3 verify_projectile_upf.py
"""
from __future__ import annotations
import sys, math, re
from pathlib import Path
import numpy as np

PSP = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials")
CANDIDATES = {
    "wpsigma0p5": PSP / "electron_gaussian_wpsigma0p5.upf",  # labelled by sigma_WP=0.5 -> sigma_pot=0.354
    "sigma0p35":  PSP / "electron_gaussian_sigma0p35.upf",   # labelled by sigma_pot~0.35
}
TARGET_SIGMA_POT = 0.5 / math.sqrt(2.0)  # 0.35355...
OUT = Path(__file__).parent


def _read_block(text: str, tag: str) -> np.ndarray:
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", text, re.DOTALL)
    if not m:
        raise ValueError(f"{tag} not found")
    return np.array([float(x) for x in m.group(1).split()])


def verify(name: str, path: Path) -> dict:
    text = path.read_text()
    r = _read_block(text, "PP_R")
    vloc_ry = _read_block(text, "PP_LOCAL")
    n = min(len(r), len(vloc_ry))
    r, vloc_ry = r[:n], vloc_ry[:n]
    v_ha = 0.5 * vloc_ry  # Ry -> Ha (matches upf2.hpp:246)

    v0 = float(v_ha[0])
    repulsive_core = v0 > 0 and np.isfinite(v0)

    # Fit erf-Gaussian: V(r) = Z * erf(r/(sqrt(2) sigma)) / r over a window away from r=0.
    # Equivalent: r*V(r) = Z * erf(r/(sqrt(2) sigma)), a clean 1-parameter (sigma) shape
    # once Z is fixed by the large-r asymptote (r*V -> Z).
    mask = (r > 0.05) & (r < 6.0)
    rr, vv = r[mask], v_ha[mask]
    rv = rr * vv
    Z = float(np.mean(rv[rr > 4.0]))  # asymptotic r*V -> Z (charge magnitude)

    # grid-search sigma_pot minimising residual of r*V vs Z*erf(r/(sqrt(2)sigma))
    sig_grid = np.linspace(0.20, 0.60, 4001)
    best_sig, best_res = None, np.inf
    for s in sig_grid:
        model = Z * np.array([math.erf(x / (math.sqrt(2.0) * s)) for x in rr])
        res = float(np.sum((rv - model) ** 2))
        if res < best_res:
            best_res, best_sig = res, s
    rms = math.sqrt(best_res / len(rr))

    # closed-form cross-check from V(0): V(0) = Z * sqrt(2/pi) / sigma
    sigma_from_v0 = Z * math.sqrt(2.0 / math.pi) / v0 if v0 > 0 else float("nan")

    monotone = bool(np.all(np.diff(v_ha[r < 3.0]) <= 1e-9))

    return dict(name=name, v0_ha=v0, repulsive_core=repulsive_core, Z=Z,
                sigma_fit=best_sig, sigma_from_v0=sigma_from_v0, fit_rms=rms,
                monotone=monotone, r=r, v_ha=v_ha)


def main() -> int:
    results = {}
    print(f"TARGET sigma_pot = sigma_WP/sqrt(2) = {TARGET_SIGMA_POT:.4f} Bohr (sigma_WP=0.5)\n")
    for name, path in CANDIDATES.items():
        if not path.exists():
            print(f"  {name}: MISSING {path}"); continue
        res = verify(name, path)
        results[name] = res
        ok_sig = abs(res["sigma_fit"] - TARGET_SIGMA_POT) < 0.02
        print(f"[{name}]  ({path.name})")
        print(f"  V(0)               = {res['v0_ha']:+.3f} Ha   "
              f"(finite repulsive core: {'YES' if res['repulsive_core'] else 'NO'})")
        print(f"  charge magnitude Z = {res['Z']:.3f}        (target ~1.0)")
        print(f"  sigma_pot (fit)    = {res['sigma_fit']:.4f} Bohr (erf-Gaussian fit, rms={res['fit_rms']:.2e} Ha)")
        print(f"  sigma_pot (V0)     = {res['sigma_from_v0']:.4f} Bohr (closed-form cross-check)")
        print(f"  monotone decay     = {'YES' if res['monotone'] else 'NO'}")
        print(f"  --> sigma matches target 0.354 within 0.02: {'PASS' if ok_sig else 'FAIL'}\n")

    # plot V(r) for the candidates (no preview; user-only per feedback_no_image_preview)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        for name, res in results.items():
            m = res["r"] < 4.0
            ax.plot(res["r"][m], res["v_ha"][m], label=f"{name} (sigma_pot={res['sigma_fit']:.3f})")
            # overlay the erf-Gaussian model
            rr = res["r"][(res["r"] > 0) & (res["r"] < 4.0)]
            model = res["Z"] * np.array([math.erf(x/(math.sqrt(2)*res["sigma_fit"]))/x for x in rr])
            ax.plot(rr, model, "--", lw=1, alpha=0.6)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel("r (Bohr)"); ax.set_ylabel("V(r) (Ha)")
        ax.set_title("Classical-electron projectile: Gaussian local potential V(r)\n"
                     "(solid = UPF PP_LOCAL, dashed = erf-Gaussian fit)")
        ax.legend(); fig.tight_layout()
        png = OUT / "projectile_upf_Vr.png"
        fig.savefig(png, dpi=130); print(f"wrote {png}")
    except Exception as e:
        print(f"(plot skipped: {e})")

    # gate verdict
    chosen = "wpsigma0p5"
    res = results.get(chosen)
    ok = res and res["repulsive_core"] and abs(res["sigma_fit"] - TARGET_SIGMA_POT) < 0.02 and 0.8 < res["Z"] < 1.2
    print(f"\nGATE [{chosen}]: {'PASS — verified -1 Gaussian projectile, sigma_pot~0.354' if ok else 'FAIL — regenerate'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
