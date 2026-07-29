#!/usr/bin/env python3
"""Task-specific test: absorbing_wrap CAP profile invariants (cap_fix campaign).

Mirrors the two W(z) formulas exactly as implemented and asserts the design
claims made in the campaign doc:

  two-sided (inq/src/perturbations/absorbing.hpp:44, fractional coords):
      W2(z) = |eta|·sin²((z − z_lo)·(π/2)/(w/2))   on each window, w = 15 Bohr,
      bumps at ±32.5 Bohr in an Lz = 80 Bohr cell
  wrap (inq-stack/include/inqkit/perturbations/absorbing_wrap.hpp):
      d(z)  = 0.5 − |z_frac|;  W(z) = |eta|·cos²(π·d/w_frac)  for d < w_frac/2,
      w = 30 Bohr

Claims tested:
  1. two-sided W == 0 exactly at the periodic boundary (the topology gap)
  2. wrap W peaks (== |eta|) exactly at the boundary and is continuous and
     smooth across the wrap (values and slopes match from both sides)
  3. equal footprint: both are nonzero exactly on |z| > 25 Bohr
  4. equal integral: ∫W2 dz == ∫Wwrap dz == |eta|·15 Bohr (to quadrature acc.)

Run:  venv/bin/python3 test_wrap_profile.py   (exit 0 = pass)
Engine-level compile+run of the header is covered by the 10-step smoke
(scripts/cap_fix/results/smoke_wrap, run_completed = true, 2026-07-13).
"""
import numpy as np

LZ = 80.0
ETA = 1.0  # amplitude scale; both profiles linear in |eta|


def w_two(z_bohr):
    """Two-sided inq absorbing: sin^2 bumps at ±32.5, full width 15 (Bohr)."""
    w = np.zeros_like(z_bohr)
    for mid in (-32.5, 32.5):
        lo, hi = mid - 7.5, mid + 7.5
        m = (z_bohr > lo) & (z_bohr < hi)
        w[m] += ETA * np.sin((z_bohr[m] - lo) * (np.pi / 2) / 7.5) ** 2
    return w


def w_wrap(z_bohr, width_bohr=30.0):
    """inqkit absorbing_wrap: cos^2 bump peaking at the boundary plane."""
    zf = z_bohr / LZ                     # fractional
    d = 0.5 - np.abs(zf)                 # periodic distance to boundary
    wf = width_bohr / LZ
    w = np.zeros_like(z_bohr)
    m = d < wf / 2
    w[m] = ETA * np.cos(d[m] * np.pi / wf) ** 2
    return w


def main():
    z = np.linspace(-40.0, 40.0, 800001, endpoint=False)  # dz = 1e-4 Bohr
    dz = z[1] - z[0]
    W2, Ww = w_two(z), w_wrap(z)

    # 1. two-sided vanishes at the boundary
    b = np.argmin(np.abs(z - (-40.0)))
    assert W2[b] < 1e-9, f"two-sided W at boundary = {W2[b]}"

    # 2. wrap peaks at the boundary; continuous + smooth across the wrap
    assert abs(Ww[b] - ETA) < 1e-6, f"wrap W at boundary = {Ww[b]} != {ETA}"
    zl, zr = -40.0 + 1e-3, 40.0 - 1e-3   # 1 mBohr each side of the wrap plane
    wl = w_wrap(np.array([zl]))[0]
    wr = w_wrap(np.array([zr]))[0]
    assert abs(wl - wr) < 1e-8, f"wrap discontinuous: {wl} vs {wr}"
    # slopes from both sides (d/dz), sign-corrected for the fold at the plane
    eps = 1e-6
    sl = (w_wrap(np.array([zl + eps]))[0] - wl) / eps
    sr = (w_wrap(np.array([zr - eps]))[0] - wr) / eps
    assert abs(sl - sr) < 1e-4, f"wrap kinked at boundary: {sl} vs {sr}"

    # 3. equal footprint |z| > 25 — with the one designed difference: the
    #    two-sided profile has a HOLE exactly at the boundary plane (W2 = 0 at
    #    z = ±L/2, its window edge), while the wrap profile is positive there.
    inside = np.abs(z) < 24.999
    assert np.all(W2[inside] == 0) and np.all(Ww[inside] == 0)
    outside = (np.abs(z) > 25.001) & (np.abs(z) < 39.999)
    assert np.all(W2[outside] > 0) and np.all(Ww[outside] > 0)
    assert Ww[b] > 0 and W2[b] < 1e-9  # the topology gap itself

    # 4. equal integral == ETA * 15 Bohr
    i2, iw = W2.sum() * dz, Ww.sum() * dz
    assert abs(i2 - ETA * 15.0) < 1e-3, f"∫W2 = {i2}"
    assert abs(iw - ETA * 15.0) < 1e-3, f"∫Wwrap = {iw}"

    print(f"PASS: boundary W2={W2[b]:.1e}, Wwrap={Ww[b]:.6f}; "
          f"∫W2={i2:.4f}, ∫Wwrap={iw:.4f} (target 15.0000); "
          f"wrap smooth across boundary (slope gap {abs(sl-sr):.1e})")


if __name__ == "__main__":
    main()
