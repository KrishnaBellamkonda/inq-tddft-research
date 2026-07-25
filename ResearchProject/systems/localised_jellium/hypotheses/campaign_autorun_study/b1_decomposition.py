#!/usr/bin/env python3
"""B1 — exact t=0 decomposition of d(H+E) for the WP-vs-classical insertion pair.

Campaign: docs/campaigns/localised_jellium_parameter_study_2/ (Energy book-keeping).
Advisor-gated route (ruling 2026-07-11): parse UPF -> exact grid decomposition at
r=12 and r=4 -> known-case gate (sum of terms reproduces measured d(H+E)) ->
only then extend radii.

Identity being tested (exact at t=0 because both runs' baths ARE the GS):
  d(H+E)(0) = [H(n_b+n_w) - H(n_b)] + [v_bg.(n_b+n_w) - (v_bg+v_ghost).n_b]
            = E_wb + E_selfH + E_bgw - E_ghb
  E_wb    = int n_w P[n_b]        (WP<->bath Hartree cross term)
  E_selfH = 1/2 int n_w P[n_w]    (WP self-Hartree, in the run's own solver)
  E_bgw   = -int phi_plus n_w     (WP<->background, phi_plus = P[n_plus]; v_bg=-phi_plus)
  E_ghb   = int v_ghost n_b       (ghost<->bath, v_ghost from the PARSED UPF, truncated
                                   at its mesh end r_max=50, lateral periodic images)
P = the run's Poisson convention: periodicity 2 = periodic x,y + OPEN z.
Implemented per lateral G: G!=0 -> (2 pi/G) exp(-G|z-z'|) kernel; G=0 -> -2 pi |z-z'|.
(The G=0 open-sheet gauge is the documented convention; constants cancel in the
d(H+E) IDENTITY because every term uses the same P.)

UPF (parsed by data 2026-07-11): V(r) = +erf(r/0.5)/r Ha on mesh r in [0,50] Bohr
(pure +1/r tail up to the mesh end; z_valence=0 so no analytic continuation beyond).
Ambiguity to test: what INQ does beyond the mesh (V=0 assumed; the gate arbitrates).
"""
import sys, glob
import numpy as np
from scipy.special import erf

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview import load_vti

HA = 27.211386
LJ = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
SW = LJ + "/scripts/campaign_autorun/runs/screening_wp"
UPF = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf"
A_HALF, N0 = 12.5, 1.312e-3


def load_fields(run):
    base = f"{SW}/{run}/results/{run}"
    d_w = load_vti(glob.glob(base + "/density_wp/*.vti")[0])
    d_b = load_vti(glob.glob(base + "/density_bath/*.vti")[0])
    return d_w, d_b


def poisson_p2(rho, x, z):
    """phi = P[rho] with periodic x,y + open z (per-G exponential kernels)."""
    nx = len(x); nz = len(z)
    dx = x[1] - x[0]; dz = z[1] - z[0]
    rk = np.fft.fftn(rho, axes=(0, 1))                      # lateral FFT only
    k = 2 * np.pi * np.fft.fftfreq(nx, dx)
    G = np.sqrt(k[:, None] ** 2 + k[None, :] ** 2)
    absdz = np.abs(z[:, None] - z[None, :])
    phik = np.empty_like(rk)
    # group identical |G| to reuse kernels
    Gq = np.round(G, 12)
    for g in np.unique(Gq):
        m = Gq == g
        cols = rk[m, :]                                     # (nm, nz)
        if g == 0.0:
            K = -2 * np.pi * absdz * dz
        else:
            K = (2 * np.pi / g) * np.exp(-g * absdz) * dz
        phik[m, :] = cols @ K.T
    return np.real(np.fft.ifftn(phik, axes=(0, 1)))


def v_ghost_grid(x, z, zp, images=1, rmax=None):
    """Parsed-UPF ghost potential at (0,0,zp), lateral periodic images."""
    import re
    t = open(UPF).read()
    rr = np.array([float(v) for v in re.search(r"<PP_R[^>]*>(.*?)</PP_R>", t, re.S).group(1).split()])
    vv = np.array([float(v) for v in re.search(r"<PP_LOCAL[^>]*>(.*?)</PP_LOCAL>", t, re.S).group(1).split()]) * 0.5  # Ry->Ha
    if rmax is None:
        rmax = rr[-1]
    L = x[-1] - x[0] + (x[1] - x[0])
    X = x[:, None, None]; Y = x[None, :, None]; Z = z[None, None, :]
    v = np.zeros((len(x), len(x), len(z)))
    for mx in range(-images, images + 1):
        for my in range(-images, images + 1):
            s = np.sqrt((X - mx * L) ** 2 + (Y - my * L) ** 2 + (Z - zp) ** 2)
            v += np.where(s <= rmax, np.interp(s, rr, vv, right=0.0), 0.0)
    return v


def decompose(run, zp, rmax=None, images=1):
    d_w, d_b = load_fields(run)
    x, z = d_w.x, d_w.z
    dv = (x[1] - x[0]) ** 2 * (z[1] - z[0])
    n_w, n_b = d_w.data, d_b.data
    n_plus = np.where(np.abs(z)[None, None, :] <= A_HALF, N0, 0.0) * np.ones_like(n_w)
    phi_b = poisson_p2(n_b, x, z)
    phi_w = poisson_p2(n_w, x, z)
    phi_p = poisson_p2(n_plus, x, z)
    vg = v_ghost_grid(x, z, zp, images=images, rmax=rmax)
    E_wb = np.sum(n_w * phi_b) * dv
    E_selfH = 0.5 * np.sum(n_w * phi_w) * dv
    E_bgw = -np.sum(phi_p * n_w) * dv
    E_ghb = np.sum(vg * n_b) * dv
    return dict(E_wb=E_wb, E_selfH=E_selfH, E_bgw=E_bgw, E_ghb=E_ghb,
                total=E_wb + E_selfH + E_bgw - E_ghb,
                norms=(np.sum(n_w) * dv, np.sum(n_b) * dv))


if __name__ == "__main__":
    MEASURED = {"wp_r12_p2": -125.8 / HA, "wp_r4_p2": -169.4 / HA}   # p2 ledger d(H+E)
    ZP = {"wp_r12_p2": -(12.5 + 12), "wp_r4_p2": -(12.5 + 4)}
    for run in ("wp_r12_p2", "wp_r4_p2"):
        r = decompose(run, ZP[run])
        print(f"--- {run}  (norms: n_w={r['norms'][0]:.4f}, n_b={r['norms'][1]:.4f})")
        for k in ("E_wb", "E_selfH", "E_bgw", "E_ghb"):
            print(f"  {k:8s} = {r[k]*HA:+10.2f} eV")
        print(f"  SUM(E_wb+E_selfH+E_bgw-E_ghb) = {r['total']*HA:+.2f} eV   "
              f"measured d(H+E) = {MEASURED[run]*HA:+.2f} eV   "
              f"gap = {(r['total']-MEASURED[run])*HA:+.2f} eV")
