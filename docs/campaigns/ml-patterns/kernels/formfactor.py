"""Projectile form factors and q-space induced-density ratio.

Campaign-local kernel (ml-patterns). PRE-GATED in T1.

F_WP(q): a Gaussian charge cloud of std sigma_pot has form factor
    F_WP(q) = exp(-q^2 sigma_pot^2 / 2)
(Fourier transform of a normalised 3D Gaussian; standard, e.g. Jackson
"Classical Electrodynamics" Ch. 3 form-factor of a charge distribution).
sigma_WP convention: the WP density std is sigma_WP; the classical-potential std
is sigma_pot = sigma_WP/sqrt(2) (project rule sigma-wp-convention). The form
factor that multiplies the q-space coupling is exp(-q^2 sigma_pot^2/2) with the
sigma_pot read from the run database column `sigma_pot_bohr`.

F_ONCV(q): the point-classical projectile imposes an EXTERNAL local potential
V_ext(r) = PP_LOCAL read verbatim (INQ ignores the is_coulomb flag — project
memory reference_inq_ignores_is_coulomb_upf_flag). A bare point charge Z gives
V(q) = -4 pi Z / q^2, i.e. form factor 1. The pseudised ONCV local potential is
regularised at small r, so its effective charge density rho = -(1/4pi) Lap V is
smeared and F_ONCV(q) = V_ext(q) q^2 / (-4 pi Z) rolls off at high q.
We compute V_ext(q) by the radial (s=r V) sine transform of the tabulated
PP_LOCAL, then normalise F_ONCV(q->0)=1. The q-range where F_ONCV ~= 1 (within a
stated tolerance) is the window in which the T2 prediction reduces to
exp(-q^2 sigma_pot^2/2); T1 establishes it.

R(q): the headline T2 metric is the azimuthally + temporally reduced ratio of
induced bath densities  R(q) = |dn_WP(q)| / |dn_classical(q)|, compared to
F_WP(q)/F_ONCV(q). Radial q-binning of a 3D FFT power spectrum is provided.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass


# ----------------------------------------------------------------------------
# Analytic form factors
# ----------------------------------------------------------------------------
def F_WP(q: np.ndarray, sigma_pot: float) -> np.ndarray:
    """Gaussian form factor exp(-q^2 sigma_pot^2 / 2)."""
    q = np.asarray(q, dtype=np.float64)
    return np.exp(-0.5 * (q * sigma_pot) ** 2)


# ----------------------------------------------------------------------------
# UPF parsing + F_ONCV from the local potential
# ----------------------------------------------------------------------------
def _read_upf_block(path: str, tag: str) -> np.ndarray:
    """Read the whitespace-separated floats inside <tag ...> ... </tag>."""
    with open(path) as fh:
        text = fh.read()
    start = text.find(f"<{tag}")
    if start < 0:
        raise KeyError(f"{tag} not found in {path}")
    gt = text.find(">", start)
    end = text.find(f"</{tag}>", gt)
    body = text[gt + 1:end]
    return np.fromstring(body, sep=" ", dtype=np.float64)


def read_upf_local(path: str):
    """Return (r, V_loc) from a UPF: PP_R mesh and PP_LOCAL (Rydberg units in UPF).

    UPF PP_LOCAL is in Rydberg; convert to Hartree (V_Ha = V_Ry / 2).
    """
    r = _read_upf_block(path, "PP_R")
    vloc_ry = _read_upf_block(path, "PP_LOCAL")
    n = min(len(r), len(vloc_ry))
    return r[:n], 0.5 * vloc_ry[:n]


def F_ONCV_from_upf(path: str, q: np.ndarray, z: float | None = None):
    """Effective form factor of the ONCV projectile's local potential.

    The local potential has a pure-Coulomb long-range tail V -> Z_s/r (signed;
    for the electron-in-jellium projectile it is REPULSIVE, V -> +1/r). The UPF
    extends it to 50 Bohr. Transforming r*V directly rings because the Coulomb
    tail is truncated. Split V = Z_s/r + dV(r) with dV short-ranged (dV -> 0 at
    large r):
        V_ext(q) = 4 pi Z_s / q^2 + dV(q),  dV(q) = (4 pi / q) int [r dV] sin(q r) dr
        F(q) = V_ext(q) / (4 pi Z_s / q^2) = 1 + (q / Z_s) int_0^inf (r V - Z_s) sin(q r) dr
    with r dV = r V - Z_s localised near the origin -> a clean, smooth F(q->0)=1.
    F(q) is the projectile form factor relative to a point charge of the same
    asymptotic strength: ~1 at low q, rolling off where pseudisation regularises
    the singularity.

    Z_s is read from the Coulomb tail itself: Z_s = median(r V) over the outer
    mesh (signed; robust to the header — INQ ignores is_coulomb anyway).
    """
    r, V = read_upf_local(path)
    q = np.asarray(q, dtype=np.float64)
    tail = slice(int(0.6 * len(r)), len(r))
    if z is None:
        z = float(np.median(r[tail] * V[tail]))   # signed asymptotic charge
    resid = r * V - z          # localised short-range remainder (-> 0 at large r)
    F = np.empty_like(q)
    for i, qi in enumerate(q):
        if qi == 0:
            F[i] = 1.0
            continue
        integ = np.trapezoid(resid * np.sin(qi * r), r)
        F[i] = 1.0 + (qi / z) * integ
    return F


def foncv_unity_range(path: str, q: np.ndarray, z: float = 1.0, tol: float = 0.05):
    """Largest q_max such that |F_ONCV(q)-1| <= tol for all q in (0, q_max].

    Returns (q_max, F_ONCV array aligned to q).
    """
    q = np.asarray(q, dtype=np.float64)
    F = F_ONCV_from_upf(path, q, z=z)
    qpos = q > 0
    qq = q[qpos]
    FF = F[qpos]
    order = np.argsort(qq)
    qq, FF = qq[order], FF[order]
    qmax = 0.0
    for qi, Fi in zip(qq, FF):
        if abs(Fi - 1.0) <= tol:
            qmax = qi
        else:
            break
    return qmax, F


# ----------------------------------------------------------------------------
# Radial (azimuthal) q-binning of a 3D field
# ----------------------------------------------------------------------------
def radial_power_spectrum(field: np.ndarray, dx: float, nbins: int = 60,
                          qmax: float | None = None):
    """Azimuthally averaged |FFT(field)| vs |q|.

    field : 3D real array in PHYSICAL order (do NOT fftshift; we use np.fft.fftn
            and the matching np.fft.fftfreq grid which is self-consistent).
    Returns (q_centres, amp) where amp = sqrt(<|F(q)|^2>) per radial shell, with
    the q=0 (DC) bin dropped.
    """
    field = np.asarray(field, dtype=np.float64)
    nx, ny, nz = field.shape
    Fk = np.fft.fftn(field)
    power = np.abs(Fk) ** 2
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=dx)
    kz = 2 * np.pi * np.fft.fftfreq(nz, d=dx)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    Q = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
    if qmax is None:
        qmax = float(Q.max()) / np.sqrt(3)  # avoid corner aliasing
    edges = np.linspace(0, qmax, nbins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    qflat = Q.ravel()
    pflat = power.ravel()
    sel = qflat <= qmax
    idx = np.digitize(qflat[sel], edges) - 1
    idx = np.clip(idx, 0, nbins - 1)
    sums = np.bincount(idx, weights=pflat[sel], minlength=nbins)
    cnts = np.bincount(idx, minlength=nbins)
    mean_power = np.where(cnts > 0, sums / np.maximum(cnts, 1), 0.0)
    amp = np.sqrt(mean_power)
    # drop DC bin
    return centres[1:], amp[1:]


def q_ratio(field_wp: np.ndarray, field_cl: np.ndarray, dx: float,
            nbins: int = 60, qmax: float | None = None, floor: float = 1e-12):
    """R(q) = amp_WP(q) / amp_classical(q) on a shared radial-q grid."""
    qc, a_wp = radial_power_spectrum(field_wp, dx, nbins, qmax)
    _, a_cl = radial_power_spectrum(field_cl, dx, nbins, qmax)
    R = a_wp / np.maximum(a_cl, floor)
    return qc, R
