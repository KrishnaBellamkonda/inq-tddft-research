"""Extract electronic stopping power S(v) from a free-Ehrenfest projectile track.

For a decelerating classical projectile (mass m, here m_e=1) the instantaneous
stopping power is the local energy loss per unit path:

    S(v(t)) = - d(KE_proj)/ds ,   KE_proj = 1/2 m |v|^2 ,   ds = |v| dt

Because v changes along the trajectory, S is a *local* derivative binned by the
instantaneous speed v(t) — NOT a single linear slope. A cross-check uses the
electronic energy gain dE_electrons/ds (equal and opposite, up to drift).

Input: electron_track.csv with columns step,time_au,x,y,z,vx,vy,vz
(written every step by run_sv_sigma0p5). Motion is along +z here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Track", "load_track", "stopping_vs_v"]


@dataclass(frozen=True)
class Track:
    t: np.ndarray
    s: np.ndarray        # path length along motion (Bohr)
    v: np.ndarray        # instantaneous speed (a.u.)
    ke: np.ndarray       # projectile kinetic energy (Ha)


def load_track(path: str, mass: float = 1.0, axis: str = "z") -> Track:
    import csv

    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh)
        for row in r:
            # tolerate a partial last line from a live-being-written CSV
            if any(row.get(k) in (None, "") for k in
                   ("time_au", "x", "y", "z", "vx", "vy", "vz")):
                continue
            rows.append(row)
    # dedupe the duplicated t=0 header row some runs write
    seen = set()
    uniq = []
    for row in rows:
        key = (row["step"], row["time_au"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(row)
    t = np.array([float(x["time_au"]) for x in uniq])
    comp = {"x": "x", "y": "y", "z": "z"}[axis]
    pos = np.array([float(x[comp]) for x in uniq])
    vx = np.array([float(x["vx"]) for x in uniq])
    vy = np.array([float(x["vy"]) for x in uniq])
    vz = np.array([float(x["vz"]) for x in uniq])
    v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    s = np.abs(pos - pos[0])
    ke = 0.5 * mass * v ** 2
    order = np.argsort(t)
    return Track(t[order], s[order], v[order], ke[order])


def stopping_vs_v(
    track: Track, *, transient_bohr: float = 3.0, window: int = 11,
) -> tuple[np.ndarray, np.ndarray]:
    """Local S(v) = -dKE/ds, after discarding an initial transient.

    Uses a Savitzky-Golay-like local linear slope (odd `window`); returns
    (v_mid, S) sampled at interior points.
    """
    mask = track.s >= transient_bohr
    s = track.s[mask]
    ke = track.ke[mask]
    v = track.v[mask]
    if s.size < window + 2:
        # too short: single global slope
        if s.size < 3:
            return np.array([]), np.array([])
        A = np.vstack([s, np.ones_like(s)]).T
        slope = np.linalg.lstsq(A, ke, rcond=None)[0][0]
        return np.array([v.mean()]), np.array([-slope])

    half = window // 2
    v_out, S_out = [], []
    for i in range(half, s.size - half):
        ss = s[i - half : i + half + 1]
        kk = ke[i - half : i + half + 1]
        A = np.vstack([ss, np.ones_like(ss)]).T
        slope = np.linalg.lstsq(A, kk, rcond=None)[0][0]
        v_out.append(v[i])
        S_out.append(-slope)
    return np.asarray(v_out), np.asarray(S_out)
