"""Physical-anchor heuristics for jellium-slab projectile runs (groups A–I).

Pure numpy/pandas (deps-clean — no matplotlib/VTK). Composes the existing
analysis kernels (``lindhard_elf`` scales, ``wp_integrity`` KL) and adds the
electron-gas scales, projectile timescales, wavepacket zero-point energy,
spreading factor, and norm/absorption diagnostics used across the localised
jellium campaign. All energies in Hartree unless a name ends ``_ev``.

References: the first jellium-slab test campaign (qa_jellium_slab_baselines) and
standard HEG relations (Giuliani & Vignale, *Quantum Theory of the Electron
Liquid*; Ashcroft & Mermin for k_TF, ω_p).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from . import lindhard_elf as _L

HA_EV = 27.211386245988


# --------------------------------------------------------------- A. e-gas scales
def electron_gas_scales(rs: float) -> dict:
    """HEG anchors for Wigner–Seitz radius ``rs`` (a.u., m_e=1)."""
    kF = _L.kF_from_rs(rs)                     # (9π/4)^(1/3)/rs
    n0 = kF**3 / (3.0 * np.pi**2)
    return {
        "rs": rs,
        "n0": n0,                              # density (a0^-3)
        "kF": kF,                              # Fermi wavevector (a0^-1)
        "vF": kF,                              # Fermi velocity (a.u.)
        "EF_ha": 0.5 * kF**2,                  # Fermi energy (Ha)
        "EF_ev": 0.5 * kF**2 * HA_EV,
        "lambda_F_friedel": np.pi / kF,        # Friedel wavelength π/kF (Bohr)
        "omega_p_ha": _L.omega_p(kF),          # plasmon frequency sqrt(4π n0) (Ha)
        "omega_p_ev": _L.omega_p(kF) * HA_EV,
        "T_plasmon_au": 2.0 * np.pi / _L.omega_p(kF),   # plasmon period (a.u.)
        "k_TF": _L.k_TF(kF),                   # Thomas–Fermi screening (a0^-1)
        "t_heg_ha_per_e": 1.104954 / rs**2,    # HEG kinetic energy per electron
    }


# ------------------------------------------------------------- B. timescales
def projectile_timescales(z0: float, v: float, slab_half: float,
                          box_half: float) -> dict:
    """Constant-mean-velocity crossing anchors (Bohr, a.u.). ``z0`` launch z<0,
    projectile moving +z at speed ``v``."""
    return {
        "v": v,
        "t_enter_slab_au": (-slab_half - z0) / v,    # reach near face −slab_half
        "t_exit_slab_au": (slab_half - z0) / v,      # reach far face +slab_half (END)
        "t_cross_au": (2.0 * slab_half) / v,         # near→far transit
        "t_reach_box_edge_au": (box_half - z0) / v,  # reach +box edge (wrap onset)
    }


# ------------------------------------------------------- C. wavepacket kinetics
def wp_zero_point(sigma_wp: float) -> dict:
    """Gaussian wavepacket zero-point KE = 3/(4σ²) (a.u.); charge std = σ/√2."""
    zp = 3.0 / (4.0 * sigma_wp**2)
    return {
        "sigma_wp": sigma_wp,
        "sigma_charge": sigma_wp / np.sqrt(2.0),
        "zero_point_ke_ha": zp,
        "zero_point_ke_ev": zp * HA_EV,
    }


# ------------------------------------------------------ D. norm / absorption
def norm_absorption(N_total: np.ndarray, N_wp: Optional[np.ndarray] = None) -> dict:
    """Absorbed norm from time series. ``N_total`` includes the projectile; if
    ``N_wp`` (WP orbital norm) is given, split the WP-orbital vs bath overflow."""
    out = {
        "N_total_0": float(N_total[0]),
        "N_total_f": float(N_total[-1]),
        "total_absorbed": float(N_total[0] - N_total[-1]),
    }
    if N_wp is not None and len(N_wp):
        wp_abs = float(N_wp[0] - N_wp[-1])
        out.update({
            "N_wp_0": float(N_wp[0]),
            "N_wp_f": float(N_wp[-1]),
            "wp_orbital_absorbed": wp_abs,
            "wp_fraction_absorbed": wp_abs / float(N_wp[0]) if N_wp[0] else float("nan"),
            "bath_overflow_absorbed": out["total_absorbed"] - wp_abs,
        })
    return out


# ------------------------------------------------------------- H. spreading
def spreading(sigma_z: np.ndarray) -> dict:
    """Longitudinal spreading factor σ_z(t)/σ_z(0)."""
    s0 = float(sigma_z[0]); sf = float(sigma_z[-1])
    return {"sigma_z_0": s0, "sigma_z_f": sf,
            "spread_factor": sf / s0 if s0 else float("nan"),
            "sigma_z_max": float(np.nanmax(sigma_z))}


# ------------------------------------------------------------- E. stopping refs
def lindhard_references(rs: float, v: float) -> dict:
    """Point- and (handled by caller for finite-σ) Lindhard stopping references."""
    kF = _L.kF_from_rs(rs)
    return {
        "S_point_ha_per_bohr": _L.stopping_power_point(v, kF),
        "S_point_ev_per_bohr": _L.stopping_power_point(v, kF) * HA_EV,
    }


# ------------------------------------------------------ high-level orchestrator
@dataclass
class Heuristics:
    run_type: str                      # "wp" | "classical" | "baseline"
    eg_scales: dict = field(default_factory=dict)
    timescales: dict = field(default_factory=dict)
    wp_kinetics: dict = field(default_factory=dict)
    norms: dict = field(default_factory=dict)
    spreading: dict = field(default_factory=dict)
    stopping_refs: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)

    def flat(self) -> dict:
        d = {}
        for grp in ("eg_scales", "timescales", "wp_kinetics", "norms",
                    "spreading", "stopping_refs", "extra"):
            for k, v in getattr(self, grp).items():
                d[f"{grp}.{k}"] = v
        return d


def _read_csv(path, **kw):
    return pd.read_csv(path, **kw) if os.path.exists(path) else None


def compute_heuristics(results_dir: str, *, rs: float, v0: float, z0: float,
                       slab_half: float, box_half: float,
                       sigma_wp: Optional[float] = None) -> Heuristics:
    """Read a run's CSVs and assemble the campaign heuristics (auto WP/classical)."""
    obs = os.path.join(results_dir, "raw", "observables")
    has_wp = os.path.exists(os.path.join(obs, "momentum_distribution.csv"))
    has_cl = os.path.exists(os.path.join(obs, "electron_track.csv"))
    rtype = "wp" if has_wp else ("classical" if has_cl else "baseline")

    H = Heuristics(run_type=rtype)
    H.eg_scales = electron_gas_scales(rs)
    H.timescales = projectile_timescales(z0, v0, slab_half, box_half)
    H.stopping_refs = lindhard_references(rs, v0)
    if sigma_wp is not None:
        H.wp_kinetics = wp_zero_point(sigma_wp)

    en = _read_csv(os.path.join(obs, "electron_number.csv"))
    N_wp = None
    if rtype == "wp":
        rs_stats = _read_csv(os.path.join(obs, "wp_real_space_stats.csv"), comment="#")
        if rs_stats is not None and "norm_check" in rs_stats:
            N_wp = rs_stats["norm_check"].values
            H.spreading = spreading(np.sqrt(rs_stats["sigma_z2"].values))
            H.extra["wp_zmean_max"] = float(rs_stats["z_mean"].max())
    if en is not None:
        H.norms = norm_absorption(en["N_total"].values, N_wp)
    return H
