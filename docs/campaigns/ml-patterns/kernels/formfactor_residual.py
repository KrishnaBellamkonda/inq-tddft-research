"""Linear-response residual / form-factor test — classical vs WP induced density.

Campaign-local kernel (ml-patterns). Chosen by the scientific panel (2026-07-06)
over Floquet/Koopman, HAVOK, optimal transport, wavelet+transfer-entropy.

THE PHYSICS (the null this test challenges)
-------------------------------------------
Linear response:  n_ind(q,w) = chi(q,w) * V_ext(q).  chi is a property of the
medium (the HEG) -> IDENTICAL for both projectiles.  Only the drive differs:
    point charge : V_ext(q) = 4*pi/q^2
    Gaussian WP  : V_ext(q) = 4*pi/q^2 * exp(-q^2 sigma^2 / 2)
so the WP is a LOW-PASS-FILTERED point charge with form factor
    F(q) = exp(-q^2 sigma^2 / 2).
In linear response, FRAME BY FRAME:  n_WP(q,t) = F(q) * n_cl(q,t).  The ratio
cancels chi in the TIME domain -- no omega-binning (the matched runs are too
short in time for any frequency-resolved technique, Dw >> w_p).

Two discriminants, both d'Alembert-safe (dividing by V_ext annihilates any rigid
f(z - v t), whose q-space translation is a pure phase that cancels in |.|):
  (1) |R(q,t)| = |n_WP(q,t)| / |n_cl(q,t)| must equal F(q).
  (2) |R(q,t)| must be FLAT in t.  t-drift => the equal-trajectory / linear
      premise breaks (deceleration mismatch, WP spreading, nonlinearity).

Fork A (the sqrt(2) trap) is resolved EMPIRICALLY, never hardcoded: fit the
exponent a in |R(q)| ~ exp(-a q^2), then test a(sigma) across sigma=0.5/3/8.
Slope 0.5 in sigma^2 => physical width is sigma_WP; slope 0.25 => sigma_pot =
sigma_WP/sqrt(2).  (See .claude/rules/sigma-wp-convention.md.)

Never np.fft.fftshift a VTI (physical order); magnitudes are shift-invariant so
we never need to.  float64 spectra; input frames stay float32 upstream.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass


HA_EV = 27.211386


# ---------------------------------------------------------------------------
# Spectra
# ---------------------------------------------------------------------------
def radial_spectrum(field3d: np.ndarray, dx: float, nbin: int | None = None):
    """3-D isotropic |q| spectrum of a real field.

    Returns (q_shells, amp, noise, count):
      amp[m]   = mean |FFT| over the m-th |q| shell   (signal proxy)
      noise[m] = std |FFT| over the shell / sqrt(count)  (SNR proxy; the reason
                 the 3-D lens beats the 1-D one -- each shell averages many modes)
      count[m] = number of Fourier modes in the shell
    Uses the FULL complex FFT (not rfft) so shell statistics are unbiased; |.|
    makes the result invariant to the projectile's rigid position (d'Alembert).
    """
    nx, ny, nz = field3d.shape
    F = np.fft.fftn(field3d)
    mag = np.abs(F)
    kx = np.fft.fftfreq(nx, d=dx) * 2.0 * np.pi
    ky = np.fft.fftfreq(ny, d=dx) * 2.0 * np.pi
    kz = np.fft.fftfreq(nz, d=dx) * 2.0 * np.pi
    QX, QY, QZ = np.meshgrid(kx, ky, kz, indexing="ij")
    qmag = np.sqrt(QX**2 + QY**2 + QZ**2)
    dq = 2.0 * np.pi / (nx * dx)                 # radial bin = grid fundamental
    qmax = qmag.max()
    nbin = nbin or int(np.floor(qmax / dq)) + 1
    edges = (np.arange(nbin + 1)) * dq
    which = np.clip((qmag / dq).astype(int), 0, nbin - 1)
    amp = np.zeros(nbin); noise = np.zeros(nbin); count = np.zeros(nbin)
    flat_m = mag.ravel(); flat_w = which.ravel()
    for m in range(nbin):
        sel = flat_m[flat_w == m]
        count[m] = sel.size
        if sel.size:
            amp[m] = sel.mean()
            noise[m] = sel.std() / np.sqrt(sel.size)
    q = 0.5 * (edges[:-1] + edges[1:])
    return q, amp, noise, count


def axial_spectrum(field3d: np.ndarray, dx: float):
    """1-D axial |q_z| spectrum on the q_perp=0 line (transverse mean then FFT).

    Returns (q_z, amp) with amp = |rfft_z(mean_xy field)|.  F(q_z) = exp(-q_z^2
    sigma^2/2) still holds on this line, so it is a valid cheap cross-check.
    """
    u = field3d.mean(axis=(0, 1))                # (nz,)
    nz = u.size
    A = np.abs(np.fft.rfft(u))
    qz = np.fft.rfftfreq(nz, d=dx) * 2.0 * np.pi
    return qz, A


# ---------------------------------------------------------------------------
# Time resampling onto a common grid
# ---------------------------------------------------------------------------
def resample_time(times: np.ndarray, arr: np.ndarray, tcommon: np.ndarray):
    """Linear-interp arr (T, Nq) from `times` onto `tcommon` (M,), per q column."""
    out = np.empty((tcommon.size, arr.shape[1]), dtype=np.float64)
    for j in range(arr.shape[1]):
        out[:, j] = np.interp(tcommon, times, arr[:, j])
    return out


# ---------------------------------------------------------------------------
# Form factor + exponent fit
# ---------------------------------------------------------------------------
def form_factor(q: np.ndarray, sigma: float) -> np.ndarray:
    """F(q) = exp(-q^2 sigma^2 / 2)  (Gaussian charge-density form factor)."""
    return np.exp(-0.5 * (q * sigma) ** 2)


def fit_gaussian_exponent(q: np.ndarray, ratio: np.ndarray, weights=None):
    """Fit log(ratio) = c - a*q^2  by weighted least squares over q>0.

    Returns dict(a, sigma_fit, c, r2).  sigma_fit = sqrt(2a) is the width implied
    by |R(q)| = exp(-a q^2); compare it to sigma_WP and sigma_pot=sigma_WP/sqrt2.
    Only finite, positive ratios inside the SNR band should be passed in.
    """
    m = np.isfinite(ratio) & (ratio > 0) & (q > 0)
    if m.sum() < 3:
        return dict(a=np.nan, sigma_fit=np.nan, c=np.nan, r2=np.nan, n=int(m.sum()))
    x = (q[m] ** 2)
    y = np.log(ratio[m])
    w = np.ones_like(x) if weights is None else np.asarray(weights)[m]
    # weighted linear fit y = c - a x
    W = w / w.sum()
    xb = (W * x).sum(); yb = (W * y).sum()
    sxx = (W * (x - xb) ** 2).sum()
    sxy = (W * (x - xb) * (y - yb)).sum()
    slope = sxy / sxx if sxx > 0 else np.nan
    a = -slope
    c = yb - slope * xb
    yhat = c + slope * x
    ss_res = (W * (y - yhat) ** 2).sum(); ss_tot = (W * (y - yb) ** 2).sum()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    sigma_fit = np.sqrt(2.0 * a) if a > 0 else np.nan
    return dict(a=float(a), sigma_fit=float(sigma_fit), c=float(c), r2=float(r2),
                n=int(m.sum()))


# ---------------------------------------------------------------------------
# The residual test (per matched pair)
# ---------------------------------------------------------------------------
@dataclass
class PairResult:
    sigma: float
    q: np.ndarray            # |q| shells (a.u.)
    ncl_t: np.ndarray        # (M, Nq) classical shell amplitude vs common time
    nwp_t: np.ndarray        # (M, Nq) WP shell amplitude vs common time
    noise_cl: np.ndarray     # (Nq,) per-shell noise floor (median over t)
    tcommon: np.ndarray      # (M,) common time grid (a.u.)
    noise_wp: np.ndarray | None = None   # (Nq,) WP noise floor; defaults to noise_cl


def residual_test(res: PairResult, snr: float = 3.0):
    """Run the form-factor / residual test on one matched pair.

    Returns a JSON-able dict: empirical exponent + implied width (vs sigma_WP and
    sigma_pot), t-flatness of |R|, and the normalized high-q residual after the
    best-fit F.  All quantities restricted to the SNR band |n_cl| > snr*noise.
    """
    q = res.q
    ncl = res.ncl_t; nwp = res.nwp_t                    # (M, Nq)
    noise = res.noise_cl                                # (Nq,)
    noise_wp = res.noise_wp if res.noise_wp is not None else noise
    ncl_bar = np.median(ncl, axis=0)                    # (Nq,) time-robust signal
    nwp_bar = np.median(nwp, axis=0)
    # SNR band: BOTH signals must clear their own noise floor. Gating on the
    # classical signal alone lets the band run past where the WP signal exists
    # (the point charge couples to all q), pulling noise into the F(q) fit.
    band = ((ncl_bar > snr * np.maximum(noise, 1e-30)) &
            (nwp_bar > snr * np.maximum(noise_wp, 1e-30)) & (q > 0))

    # frame-by-frame ratio in the band
    both = (ncl > snr * noise[None, :]) & (nwp > snr * noise_wp[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(both, nwp / ncl, np.nan)           # (M, Nq)
    with np.errstate(invalid="ignore"):
        allnan = ~np.isfinite(R).any(axis=0)
        R_bar = np.full(q.shape, np.nan)
        R_bar[~allnan] = np.nanmedian(R[:, ~allnan], axis=0)
        mad = np.full(q.shape, np.nan)
        mad[~allnan] = np.nanmedian(np.abs(R[:, ~allnan] - R_bar[~allnan][None, :]), axis=0)
        flat = mad / np.abs(R_bar)
    flatness_band = float(np.nanmedian(flat[band])) if band.any() else np.nan

    # empirical exponent: fit F(q) ONLY on the CONTIGUOUS DESCENDING ARM of |R|,
    # from q=0 up to where |R| bottoms out. Beyond that turnaround |R| rises onto
    # the numerical noise-ratio floor (the total-minus-wp blob-subtraction residual
    # leaks ~5-10% high-q power, so |R| plateaus rather than -> 0) and would bend
    # the log-linear fit and bias sigma_fit high. Parameter-free; the wider `band`
    # is still used for the excess / t-flatness diagnostics.
    fit_band = np.zeros_like(band)
    idx = np.where(band & np.isfinite(R_bar) & (R_bar > 0.0))[0]
    if idx.size >= 3:
        turn = idx[int(np.argmin(R_bar[idx]))]          # first bottom of |R|
        arm = idx[idx <= turn]
        fit_band[arm] = (R_bar[arm] < 1.5)              # drop any q with |R|>1.5 (noise)
    fit = fit_gaussian_exponent(q[fit_band], R_bar[fit_band], weights=ncl_bar[fit_band])

    # normalized additive residual after best-fit F(q): delta = n_WP - F*n_cl
    if np.isfinite(fit["a"]):
        Ffit = np.exp(fit["c"]) * np.exp(-fit["a"] * q**2)   # includes amplitude c
        delta = nwp_bar - Ffit * ncl_bar
        resid_norm = delta / np.maximum(noise, 1e-30)
    else:
        resid_norm = np.full_like(q, np.nan)

    # high-q excess: median normalized residual in the upper half of the band
    band_q = q[band]
    excess = np.nan
    if band_q.size >= 4:
        qhi = np.median(band_q)
        hi = band & (q >= qhi)
        excess = float(np.nanmedian(resid_norm[hi])) if hi.any() else np.nan

    sigma = res.sigma
    # SNR-adequacy: a FLOORED fit (the WP form factor e-folds within a few shells,
    # so the descending arm hits the ~5-10% blob-subtraction plateau) grossly
    # UNDER-estimates the exponent. In linear response a = 0.5*sigma_phys^2, with
    # the two candidate widths giving a = 0.5*sigma^2 (sigma_WP) or 0.25*sigma^2
    # (sigma_pot). Require the measured a to reach at least ~30% of the smaller
    # (sigma_pot) expectation: a >= 0.15*sigma^2. Below that the floor has eaten the
    # signal (panel: sigma>=3 is SNR-dead). Synthetic clean data passes trivially.
    q_efold = np.sqrt(2.0) / sigma if sigma > 0 else np.inf
    snr_adequate = bool(sigma > 0 and np.isfinite(fit["a"])
                        and fit["a"] >= 0.15 * sigma**2 and fit["n"] >= 5)
    return dict(
        sigma=sigma,
        sigma_pot=sigma / np.sqrt(2.0),
        q_efold=float(q_efold),
        snr_adequate=snr_adequate,
        n_common=int(res.tcommon.size),
        t_overlap_au=float(res.tcommon[-1]),
        band_qmax=float(band_q.max()) if band_q.size else 0.0,
        band_nq=int(band.sum()),
        fit_nq=int(fit_band.sum()),
        fit_qmax=float(q[fit_band].max()) if fit_band.any() else 0.0,
        exponent_a=fit["a"],
        sigma_fit=fit["sigma_fit"],
        fit_r2=fit["r2"],
        matches_sigma_wp=bool(np.isfinite(fit["sigma_fit"]) and
                              abs(fit["sigma_fit"] - sigma) < abs(fit["sigma_fit"] - sigma / np.sqrt(2.0))),
        t_flatness=flatness_band,
        highq_excess_over_noise=excess,
        # arrays for plotting / downstream (lists)
        q=q.tolist(),
        R_median=[float(x) for x in R_bar],
        F_at_sigma_wp=[float(x) for x in form_factor(q, sigma)],
        F_at_sigma_pot=[float(x) for x in form_factor(q, sigma / np.sqrt(2.0))],
        resid_over_noise=[float(x) for x in resid_norm],
        band_mask=[bool(x) for x in band],
    )


def collapse_fork_a(per_pair: list[dict]):
    """Resolve Fork A across sigma: fit exponent_a(sigma) vs sigma^2.

    In linear response a = 0.5*sigma_phys^2, so a-vs-sigma^2 is a line through the
    origin with slope 0.5 (physical width = sigma_WP) or 0.25 (sigma_pot). ONLY
    SNR-adequate pairs are used (SNR-dead broad-sigma pairs fit the noise floor,
    not F(q), and would poison the slope). With <2 adequate points the collapse is
    INCONCLUSIVE -- report the single-point lean instead of a bogus slope.
    """
    adq = [p for p in per_pair if p.get("snr_adequate")
           and np.isfinite(p.get("exponent_a", np.nan)) and p["sigma"] > 0]
    dead = [p["sigma"] for p in per_pair if not p.get("snr_adequate")]
    if len(adq) < 2:
        # single-point lean: compare implied width to sigma_WP vs sigma_pot
        if adq:
            p = adq[0]
            lean = ("sigma_WP" if abs(p["sigma_fit"] - p["sigma"])
                    < abs(p["sigma_fit"] - p["sigma_pot"]) else "sigma_pot")
            return dict(ok=False, inconclusive=True, snr_dead_sigmas=dead,
                        reason=f"only {len(adq)} SNR-adequate sigma (need >=2)",
                        single_point=dict(sigma=p["sigma"], sigma_fit=p["sigma_fit"],
                                          leans=lean))
        return dict(ok=False, inconclusive=True, snr_dead_sigmas=dead,
                    reason="no SNR-adequate sigma pairs")
    s2 = np.array([p["sigma"]**2 for p in adq]); a = np.array([p["exponent_a"] for p in adq])
    slope = float((s2 * a).sum() / (s2 * s2).sum())     # through-origin LS
    return dict(ok=True, slope=slope, slope_sigma_wp=0.5, slope_sigma_pot=0.25,
                selects="sigma_WP" if abs(slope - 0.5) < abs(slope - 0.25) else "sigma_pot",
                snr_dead_sigmas=dead,
                points=[dict(sigma=p["sigma"], a=p["exponent_a"]) for p in adq])
