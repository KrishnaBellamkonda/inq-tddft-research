#!/usr/bin/env python3
"""Mandatory pre-run grid-cutoff (aliasing) guard for every TDDFT projectile run.

A real-space grid of spacing h can only represent momenta up to the Nyquist
wavevector k_Nyq = π/h, i.e. single-particle kinetic energies up to
E_cut = ½(π/h)². A projectile carrying momentum beyond k_Nyq ALIASES: high-k
content wraps to spurious momenta and injects energy, so the dynamics — and any
stopping power — are garbage (positive late-time energy slope, runaway S). This
guard BLOCKS such a run *before* it launches.

POLICY (two-tier, user decision 2026-06-27)
  Wavepacket projectile (Gaussian, NOT monochromatic):
      σ_p = 1/(√2·σ_WP)   ← momentum-space std; ψ ∝ exp(−r²/2σ²) (wavepacket.hpp:254)
                            σ_WP is the .sigma() ENVELOPE width (density std = σ_WP/√2);
                            σ_WP=0.5 ⇒ σ_p=1.41.  [NOT 1/(2σ_WP) — that was a √2 error.]
      p0  = √(2·E_drift)  ← mean (drift) momentum
      aliased tail fraction f = 1 − Φ((k_Nyq − p0)/σ_p)
      • HARD BLOCK if f > block_tail_frac (default 2%)            — genuinely aliased
      • WARN       if p0 + n_sigma·σ_p > k_Nyq but f ≤ 2%          — marginal, runs
      • PASS       otherwise
      (Empirical calibration: f≈1% ran clean, ≈5% borderline, ≈18% destroyed.)

  Classical projectile (monochromatic, p0 = √(2·E_kin)):
      • HARD BLOCK if E_cut < margin·E_kin (margin default 1.10, i.e. "10% higher").
      [response_factor>1 accounts for induced-electron pickup to ~response_factor·v;
       a projectile at v can forward-scatter electrons to ~2v, so response_factor≈2 is
       the physically complete bound — default 1.0 follows the projectile-momentum
       convention; set it if you want the stricter check.]

MANDATORY for every projectile run; a BLOCK is a HARD STOP. The ONLY exception is an
explicit user instruction to override (--override / override=True). A WARN proceeds.

CLI:
  python3 cutoff_guard.py --spacing 0.5 --kind wp        --energy-ev 340 --sigma-wp 0.5
  python3 cutoff_guard.py --spacing 0.5 --kind classical --energy-ev 122
Exit 0 = PASS/WARN/override, 3 = BLOCK, 2 = bad args. Importable: check_run(...) -> dict.
"""
from __future__ import annotations
import argparse, math, sys

HA_EV = 27.211386245988


def grid_cutoff(spacing_bohr: float):
    """Return (k_Nyq [a.u.], E_cut [Ha]) for a real-space grid of the given spacing."""
    k_nyq = math.pi / spacing_bohr
    return k_nyq, 0.5 * k_nyq * k_nyq


def _phi(z: float) -> float:
    """Standard-normal CDF via erf (no scipy dependency)."""
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def check_run(spacing_bohr: float, kind: str, energy_ev: float, *,
              sigma_wp_bohr: float | None = None, n_sigma: float = 3.0,
              block_tail_frac: float = 0.02, classical_margin: float = 1.10,
              response_factor: float = 1.0, override: bool = False) -> dict:
    """Verdict dict: status(pass|warn|block), ok, block, reason + all numbers (eV/au)."""
    k_nyq, e_cut_ha = grid_cutoff(spacing_bohr)
    e_cut_ev = e_cut_ha * HA_EV
    e_ha = energy_ev / HA_EV
    kind = kind.lower()
    r: dict = dict(kind=kind, spacing_bohr=spacing_bohr, k_nyq=k_nyq,
                   e_cut_ha=e_cut_ha, e_cut_ev=e_cut_ev, energy_ev=energy_ev)

    if kind == "classical":
        p0 = math.sqrt(2.0 * e_ha)
        p_eff = response_factor * p0
        need_ev = classical_margin * 0.5 * p_eff * p_eff * HA_EV
        status = "pass" if e_cut_ev >= need_ev else "block"
        r.update(p0=p0, p_eff=p_eff, required_ecut_ev=need_ev, headroom=e_cut_ev / need_ev,
                 reason=(f"E_cut={e_cut_ev:.0f} eV {'≥' if status=='pass' else '<'} "
                         f"{classical_margin:.2f}×E_kin"
                         f"{f'(×{response_factor:g} response)' if response_factor != 1 else ''}"
                         f"={need_ev:.0f} eV"))

    elif kind in ("wp", "wavepacket"):
        if sigma_wp_bohr is None or sigma_wp_bohr <= 0:
            return dict(status="block", ok=False, block=True, reason="WP guard needs --sigma-wp > 0")
        sigma_p = 1.0 / (math.sqrt(2.0) * sigma_wp_bohr)
        p0 = math.sqrt(2.0 * e_ha)
        p_max = p0 + n_sigma * sigma_p
        f = 1.0 - _phi((k_nyq - p0) / sigma_p)            # aliased tail fraction
        strict_ok = k_nyq >= p_max
        if f > block_tail_frac:
            status = "block"
        elif not strict_ok:
            status = "warn"
        else:
            status = "pass"
        r.update(sigma_wp_bohr=sigma_wp_bohr, sigma_p=sigma_p, p0=p0, p_max=p_max,
                 aliased_tail_frac=f, strict_3sigma_ok=strict_ok, block_tail_frac=block_tail_frac,
                 reason=(f"aliased tail = {f*100:.2f}% "
                         f"({'>' if f>block_tail_frac else '≤'} {block_tail_frac*100:g}% block); "
                         f"k_Nyq={k_nyq:.2f} {'≥' if strict_ok else '<'} "
                         f"p0+{n_sigma:g}σ_p={p_max:.2f}"))
    else:
        return dict(status="block", ok=False, block=True, reason=f"unknown kind {kind!r} (classical|wp)")

    r["status"] = status
    r["ok"] = (status == "pass")
    r["warn"] = (status == "warn")
    r["override"] = override
    r["block"] = (status == "block" and not override)
    return r


def _print(r: dict) -> int:
    st = r.get("status", "block")
    tag = {"pass": "PASS", "warn": "WARN (marginal — proceeds)",
           "block": "OVERRIDE (user)" if r.get("override") else "BLOCK"}[st]
    print(f"[cutoff-guard] {tag}: {r.get('reason','')}")
    if "e_cut_ev" in r:
        print(f"  spacing={r['spacing_bohr']} Bohr  k_Nyq={r['k_nyq']:.3f} a.u.  E_cut={r['e_cut_ev']:.0f} eV")
    return 3 if r.get("block") else 0


def _selftest():
    def e(v): return 0.5 * v * v * HA_EV
    # two-tier: v5@0.5 (18%) blocks; v5@0.35 (0.25%) warns (runs); v3@0.5 (1%) warns; v4@0.5 (5.4%) blocks
    assert check_run(0.5,  "wp", e(5), sigma_wp_bohr=0.5)["status"] == "block"
    assert check_run(0.35, "wp", e(5), sigma_wp_bohr=0.5)["status"] == "warn"
    assert check_run(0.5,  "wp", e(3), sigma_wp_bohr=0.5)["status"] == "warn"
    assert check_run(0.5,  "wp", e(4), sigma_wp_bohr=0.5)["status"] == "block"
    assert check_run(0.5,  "classical", 122.0)["status"] == "pass"
    o = check_run(0.5, "wp", e(5), sigma_wp_bohr=0.5, override=True)
    assert o["status"] == "block" and not o["block"], o     # override clears the hard stop
    print("[cutoff-guard] selftest OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--spacing", type=float, help="grid spacing h (Bohr)")
    ap.add_argument("--kind", choices=("classical", "wp", "wavepacket"))
    ap.add_argument("--energy-ev", type=float, help="projectile kinetic (classical) / drift (WP) energy, eV")
    ap.add_argument("--sigma-wp", type=float, default=None, help="WP envelope σ (Bohr); density std = σ/√2")
    ap.add_argument("--n-sigma", type=float, default=3.0, help="WARN if p0+n_sigma·σ_p > k_Nyq")
    ap.add_argument("--block-tail-frac", type=float, default=0.02, help="HARD BLOCK if aliased tail > this")
    ap.add_argument("--classical-margin", type=float, default=1.10)
    ap.add_argument("--response-factor", type=float, default=1.0,
                    help="classical induced-electron momentum factor (1=projectile, ~2=electron pickup)")
    ap.add_argument("--override", action="store_true", help="explicit user override — bypass a BLOCK")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        _selftest(); sys.exit(0)
    if a.spacing is None or a.kind is None or a.energy_ev is None:
        ap.error("need --spacing, --kind, --energy-ev")
    sys.exit(_print(check_run(
        a.spacing, a.kind, a.energy_ev, sigma_wp_bohr=a.sigma_wp, n_sigma=a.n_sigma,
        block_tail_frac=a.block_tail_frac, classical_margin=a.classical_margin,
        response_factor=a.response_factor, override=a.override)))
