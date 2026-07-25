"""High-level analysis pipeline for the ml-patterns campaign (T2 form-factor, T3 wake).

Ties the kernels (celldb, normaliser, formfactor, pod, dmd) into the two headline
analyses with a single tunable CONFIG (ADR 0011: tuned on calibration, frozen,
verdict read on held-out). float32, <=max_frames per series.
"""
from __future__ import annotations
import os, glob
import numpy as np
from dataclasses import dataclass, asdict, field

from . import celldb, normaliser as N, formfactor as FF, dmd as D, pod as P
from inqview import load_vti

ONCV_UPF = ("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
            "shared/pseudopotentials/electron-ONCV-1.2.upf")
DX3 = None  # set per-cell


@dataclass
class FormFactorConfig:
    nbins: int = 32             # wider shells -> more voxels/shell -> less ratio noise
    qmax: float = 1.9           # 1/Bohr; <= F_ONCV unity range (~1.89)
    q_fit_lo: float = 0.2       # low-q window for shape comparison / normalisation
    q_fit_hi: float = 1.6       # absolute high-q cutoff (stay where F_ONCV ~ 1)
    pred_floor: float = 0.2     # truncate fit where prediction drops below this
    smooth: int = 3             # moving-average window over q-shells (odd)
    max_frames: int = 90
    window_frac: tuple = (0.05, 0.6)  # temporal window (early-mid, constant-v-ish)
    normalise_at_qlo: bool = True

    def to_dict(self):
        return asdict(self)


def _smooth(a, w):
    if w is None or w < 3:
        return a
    w = int(w) | 1
    k = np.ones(w) / w
    return np.convolve(a, k, mode="same")


@dataclass
class WakeConfig:
    rank: int = 12
    max_frames: int = 250
    window_frac: tuple = (0.0, 0.5)   # early near-constant-velocity stretch
    project_pod_rank: int = 30        # POD pre-compression before DMD
    line_axis: int = 2                # z = projectile direction

    def to_dict(self):
        return asdict(self)


# ----------------------------------------------------------------------------
# Bath induced-density resolution (handles _wf vs legacy density_system)
# ----------------------------------------------------------------------------
def _electron_count(field, dx):
    return float(field.sum()) * dx ** 3


def load_bath_induced(bath_dir, gs_path, dx, max_frames, n_bath_expected,
                      wp_dir="", total_dir=""):
    """Return (delta (T,nx,ny,nz), axes, method) with the GS-subtracted bath series.

    Resolution: if density_system integrates to ~N_bath -> use it (bath-only,
    the _wf convention). Else if a WP series exists -> total - wp (frame-min
    aligned). Else fall back to density_total with a logged WP-included caveat.
    """
    series, axes, idx = N.load_series(bath_dir, max_frames)
    gs = N.load_gs(gs_path)
    ne = _electron_count(series[len(series) // 2], dx)
    method = "system_bathonly"
    if abs(ne - (n_bath_expected + 1)) < 0.5 and total_dir and wp_dir and os.path.isdir(wp_dir):
        # WP-included; subtract the WP electron series (frame-aligned by index)
        wp_series, wp_axes, _ = N.load_series(wp_dir, max_frames)
        t = min(len(series), len(wp_series))
        series = series[:t] - N.cogrid(wp_series[:t], wp_axes, axes)
        method = "total_minus_wp"
    elif abs(ne - (n_bath_expected + 1)) < 0.5:
        method = "total_wp_included_CAVEAT"
    delta = (series - gs[None]).astype(np.float32)
    return delta, axes, method


# ----------------------------------------------------------------------------
# T2 form-factor: R(q) per cell + agreement vs F_WP/F_ONCV
# ----------------------------------------------------------------------------
def _temporal_window(T, frac):
    a = int(frac[0] * T); b = max(a + 2, int(frac[1] * T))
    return a, min(b, T)


def cell_Rq(cell, cfg: FormFactorConfig, n_bath=162):
    dx = cell["dx"]
    wp_delta, wp_axes, wp_method = load_bath_induced(
        cell["wp_bath_dir"], cell["wp_gs"], dx, cfg.max_frames, n_bath,
        wp_dir=cell.get("wp_wp_dir", ""), total_dir=cell.get("wp_total_dir", ""))
    cl_delta, cl_axes, cl_method = load_bath_induced(
        cell["cl_bath_dir"], cell["cl_gs"], dx, cfg.max_frames, n_bath)
    # co-grid classical onto WP grid if needed
    if not N.grids_match(cl_axes, wp_axes):
        cl_delta = N.cogrid(cl_delta, cl_axes, wp_axes)
    # temporal window
    aw, bw = _temporal_window(len(wp_delta), cfg.window_frac)
    ac, bc = _temporal_window(len(cl_delta), cfg.window_frac)
    # per-frame radial spectra, per-q median (robust temporal reduction)
    def med_spec(delta, a, b):
        amps = []
        q = None
        for k in range(a, b):
            qq, amp = FF.radial_power_spectrum(delta[k], dx, nbins=cfg.nbins, qmax=cfg.qmax)
            q = qq; amps.append(amp)
        return q, np.median(np.asarray(amps), axis=0)
    q, amp_wp = med_spec(wp_delta, aw, bw)
    _, amp_cl = med_spec(cl_delta, ac, bc)
    amp_wp_s = _smooth(amp_wp, cfg.smooth)
    amp_cl_s = _smooth(amp_cl, cfg.smooth)
    R = amp_wp_s / np.maximum(amp_cl_s, 1e-12)
    return {"q": q, "R": R, "wp_method": wp_method, "cl_method": cl_method,
            "amp_wp": amp_wp_s, "amp_cl": amp_cl_s}


def predict_FF(q, sigma_pot, oncv_upf=ONCV_UPF):
    fwp = FF.F_WP(q, sigma_pot)
    foncv = FF.F_ONCV_from_upf(oncv_upf, q)
    return fwp / np.maximum(foncv, 1e-6)


def cell_agreement(cell, cfg: FormFactorConfig, n_bath=162):
    res = cell_Rq(cell, cfg, n_bath)
    q, R = res["q"], res["R"]
    pred = predict_FF(q, cell["sigma_pot"])
    # window: above q_fit_lo, below q_fit_hi, and where the prediction is still
    # appreciable (>= pred_floor) so the ratio is not dominated by division noise.
    sel = (q >= cfg.q_fit_lo) & (q <= cfg.q_fit_hi) & (pred >= cfg.pred_floor)
    qs, Rs, ps = q[sel], R[sel].copy(), pred[sel]
    if cfg.normalise_at_qlo and len(Rs) > 0 and Rs[0] != 0:
        Rs = Rs / Rs[0] * ps[0]   # anchor measured to prediction at window start
    rel = np.abs(Rs - ps) / np.maximum(np.abs(ps), 1e-6)
    frac_within = float(np.mean(rel <= 0.20)) if len(rel) else 0.0
    # robust slope cross-check: log R ~ -q^2 sigma_eff^2/2 in the F_ONCV~1 window
    sigma_eff = np.nan
    if len(qs) >= 3:
        pos = Rs > 0
        if pos.sum() >= 3:
            slope = np.polyfit(qs[pos] ** 2, np.log(Rs[pos]), 1)[0]
            if slope < 0:
                sigma_eff = float(np.sqrt(-2 * slope))
    sig_rel = (abs(sigma_eff - cell["sigma_pot"]) / cell["sigma_pot"]
               if np.isfinite(sigma_eff) else np.nan)
    res.update({"pred": pred, "q_sel": qs, "R_sel": Rs, "pred_sel": ps,
                "frac_within20": frac_within,
                "median_rel": float(np.median(rel)) if len(rel) else np.nan,
                "sigma_eff": sigma_eff, "sigma_pot": cell["sigma_pot"],
                "sigma_eff_rel": sig_rel, "sigma_wp": cell["sigma_wp"],
                "n_window": int(len(qs))})
    return res


# ----------------------------------------------------------------------------
# T3 wake: DMD dominant angular frequency vs omega_p; wavelength vs 2pi v/omega_p
# ----------------------------------------------------------------------------
HA_EV = 27.211386


def cell_wake_dmd(cell, cfg: WakeConfig, n_bath=162):
    dx = cell["dx"]
    delta, axes, method = load_bath_induced(
        cell["wp_bath_dir"], cell["wp_gs"], dx, cfg.max_frames, n_bath,
        wp_dir=cell.get("wp_wp_dir", ""), total_dir=cell.get("wp_total_dir", ""))
    fdt = cell["frame_dt_au_wp"]
    if fdt is None or not np.isfinite(fdt):
        return {"method_invalid": "no frame_dt", "wp_method": method}
    T = len(delta)
    a = int(cfg.window_frac[0] * T); b = max(a + 4, int(cfg.window_frac[1] * T))
    b = min(b, T)
    # Nyquist guard: dt < pi/omega_p
    omega_p_au = cell["omega_p_ev"] / HA_EV
    nyquist_ok = fdt < np.pi / omega_p_au
    # flatten, POD pre-compress, DMD
    sub = delta[a:b].reshape(b - a, -1).T.astype(np.float32)   # (n_features, t)
    pod = P.pod(sub, rank=min(cfg.project_pod_rank, b - a - 1))
    coeffs = pod.coeffs   # (r, t) latent dynamics
    try:
        res = D.dmd(coeffs, dt=fdt, rank=min(cfg.rank, coeffs.shape[0]))
    except Exception as e:
        return {"method_invalid": f"dmd_fail:{e}", "wp_method": method,
                "nyquist_ok": nyquist_ok}
    i, w_au, g, amp = res.dominant()
    w_ev = w_au * HA_EV
    v = cell["velocity_au"]
    lam_dmd = 2 * np.pi * v / w_au if w_au > 0 else np.nan
    lam_theory = 2 * np.pi * v / omega_p_au
    return {"wp_method": method, "nyquist_ok": bool(nyquist_ok), "fdt": fdt,
            "omega_dmd_ev": w_ev, "omega_p_ev": cell["omega_p_ev"],
            "growth": g, "amp": amp, "lambda_dmd": lam_dmd,
            "lambda_theory": lam_theory, "velocity_au": v,
            "window": (a, b), "n_pod_modes": coeffs.shape[0]}
