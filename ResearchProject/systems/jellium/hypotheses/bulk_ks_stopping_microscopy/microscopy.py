"""Microscopic diagnostics for the classical-vs-wavepacket stopping gap.

Plan: docs/plans/bulk-jellium-ks-stopping.md
Handover: docs/handovers/bulk-jellium-ks-stopping.md

THE QUESTION
------------
Both completed twin pairs give S_classical / S_WP ~ 5-6x, and the density lever
moved it only 13 % (6.49x at r_s = 5.702 -> 5.65x at r_s = 3.987). The stopping
numbers say THAT the projectiles differ; they do not say WHY. This module builds
the real-space picture: what the bath actually does around each projectile.

WHAT DRAG IS, MICROSCOPICALLY
-----------------------------
A charged projectile polarises the electron gas. The induced density does not sit
symmetrically about the projectile -- it LAGS, because the gas responds at a
finite rate (set by omega_p) while the projectile keeps moving. That lag leaves
more induced charge BEHIND the projectile than in front, and the resulting
induced field pulls backwards. Stopping power IS that asymmetry:

    S  =  -q * E_induced(at the projectile)   ~   the front/back imbalance of dn

So the microscopic question "why does the classical projectile slow faster?"
becomes two measurable sub-questions:
  (a) does the classical projectile induce a LARGER dn than the wavepacket?
  (b) is the classical dn more ASYMMETRIC (more lag) than the wavepacket's?

(a) is about coupling strength -- a compact charge couples to shorter wavelengths
than a smeared one. (b) is about the response being able to keep up.

THE DECOMPOSITION THAT MATTERS
------------------------------
`density_delta` = n_total(t) - n_total(0). For the CLASSICAL run that is exactly
the bath response, because the projectile is an external potential and is not
part of n at all.

For the WAVEPACKET run it is NOT: the WP is itself an occupied orbital, so
density_delta contains the packet moving (a big positive blob at its current
position, a negative one where it started) on top of the bath response. Comparing
that to the classical density_delta would be comparing different things.

So for the WP run this module forms the BATH explicitly:

    n_bath(t)   = n_total(t) - n_wp(t)
    dn_bath(t)  = n_bath(t) - n_bath(0)

and it is `dn_bath` that is compared with the classical `density_delta`. This is
the only apples-to-apples comparison available, and getting it wrong would make
the wavepacket look like it induces an enormous response when most of that is
just the packet itself.

PERIODICITY
-----------
Every profile is taken on the box's own z axis and then re-expressed in the
PROJECTILE FRAME (zeta = z - z_proj) with wrapping into (-L/2, +L/2]. The
projectile position comes from the classical track (`electron_track.csv`) or,
for the WP, the circular centroid (`wp_real_space_stats.csv`, z_mean_circ) --
never the naive centroid, which slides to a wrong value near a cell face.
"""
from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HA_TO_EV = 27.211386


# ---------------------------------------------------------------------------
# Frame discovery
# ---------------------------------------------------------------------------

def frame_paths(run_dir: str | Path, kind: str) -> list[Path]:
    """Sorted VTI paths for one field kind, e.g. 'density_total'."""
    p = Path(run_dir) / "results" / "raw" / "vti" / kind
    return [Path(f) for f in sorted(glob.glob(str(p / "*.vti")))]


def frame_step(path: str | Path) -> int:
    """Extract the timestep index from a `*_t000123.vti` filename."""
    m = re.search(r"_t(\d+)\.vti$", str(path))
    if not m:
        raise ValueError(f"cannot parse a step index from {path}")
    return int(m.group(1))


def stride_for(run_dir: str | Path, kind: str, target_frames: int = 40) -> int:
    """Frame stride that yields ~``target_frames`` regardless of write cadence.

    WHY THIS EXISTS. A fixed frame stride silently assumes every run wrote its
    density at the same cadence. That stopped being true on 2026-07-31 when
    ``WRITE_EVERY`` went 2 -> 8: runs made before have 301-347 frames, runs made
    after have ~87. A blanket ``stride=8`` then samples the new runs 4x more
    sparsely than the old ones — which is exactly how the r_s=3.987 sigma=1 pair
    produced a spurious residual of 3.32 against ~2.1 for every other pair, and
    briefly looked like a real width dependence.

    Striding to a TARGET COUNT instead makes the sampling comparable across runs
    whatever their cadence. Always prefer this over a hard-coded stride when
    comparing runs against each other.
    """
    n = len(frame_paths(run_dir, kind))
    if n == 0:
        raise FileNotFoundError(f"no {kind} frames under {run_dir}")
    return max(1, n // max(1, target_frames))


# ---------------------------------------------------------------------------
# Axial profiles
# ---------------------------------------------------------------------------

@dataclass
class AxialSeries:
    """On-axis z profiles of one field over time.

    profiles[i, :] is the profile at time t[i] on the axis `z` (Bohr).
    `radius_bohr` records how the transverse average was taken: 0 means a single
    on-axis line, >0 means averaged over a cylinder of that radius (which is far
    less noisy and is what the drag actually integrates over).
    """
    t: np.ndarray
    z: np.ndarray
    profiles: np.ndarray
    radius_bohr: float
    kind: str


def axial_series(run_dir: str | Path, kind: str, *, dt_au: float,
                 stride: int = 1, radius_bohr: float = 3.0,
                 loader=None) -> AxialSeries:
    """Load `kind` frames and reduce each to a 1-D z profile.

    The transverse reduction is a MEAN over the cylinder x^2 + y^2 <= radius^2,
    not a single grid line: a one-line cut through a 3-D density is noisy and
    grid-orientation dependent, while the quantity that produces drag is the
    charge in a tube around the trajectory.
    """
    if loader is None:                       # imported lazily: VTK is slow
        from inqview import load_vti as loader
    paths = frame_paths(run_dir, kind)[::stride]
    if not paths:
        raise FileNotFoundError(f"no {kind} frames under {run_dir}")
    first = loader(paths[0])
    X, Y = np.meshgrid(first.x, first.y, indexing="ij")
    mask = (X**2 + Y**2) <= radius_bohr**2 if radius_bohr > 0 else None
    if mask is not None and not mask.any():
        raise ValueError(f"radius {radius_bohr} Bohr encloses no grid points")

    prof, times = [], []
    for p in paths:
        f = loader(p)
        if mask is None:
            ix = int(np.argmin(np.abs(f.x))); iy = int(np.argmin(np.abs(f.y)))
            prof.append(f.data[ix, iy, :])
        else:
            prof.append(f.data[mask, :].mean(axis=0))
        times.append(frame_step(p) * dt_au)
    return AxialSeries(np.asarray(times), first.z, np.asarray(prof),
                       radius_bohr, kind)


def induced_series(run_dir: str | Path, *, dt_au: float, stride: int = 1,
                   radius_bohr: float = 3.0, loader=None) -> AxialSeries:
    """Induced density n(t) - n(0), DERIVED from ``density_total``.

    The runs used to also write a full-resolution ``density_delta`` field, which
    is exactly n(t) - n(t0) by construction (see
    ``inqkit/observables/density_delta.hpp``: "defined relative to a
    user-supplied reference density n(r, t0)"). That made it a byte-for-byte
    duplicate of something already derivable from ``density_total``, at ~10 GB
    per run — it was the single largest redundancy on disk and is no longer
    written (``emit_raw_vti = false``, 2026-07-31).

    This function reproduces it exactly by subtracting the first frame. Prefer it
    over reading ``density_delta``; that directory is absent for any run made
    after the change and has been pruned from the earlier ones.
    """
    tot = axial_series(run_dir, "density_total", dt_au=dt_au, stride=stride,
                       radius_bohr=radius_bohr, loader=loader)
    return AxialSeries(tot.t, tot.z, tot.profiles - tot.profiles[0],
                       radius_bohr, "induced_total")


def bath_series(run_dir: str | Path, *, dt_au: float, stride: int = 1,
                radius_bohr: float = 3.0, loader=None) -> AxialSeries:
    """Induced BATH density for a wavepacket run: (n_total - n_wp) minus its t=0 value.

    See the module docstring: for a WP run `density_delta` is contaminated by the
    packet itself, so it must not be compared with the classical `density_delta`.
    """
    tot = axial_series(run_dir, "density_total", dt_au=dt_au, stride=stride,
                       radius_bohr=radius_bohr, loader=loader)
    wp = axial_series(run_dir, "density_wp", dt_au=dt_au, stride=stride,
                      radius_bohr=radius_bohr, loader=loader)
    n = min(len(tot.t), len(wp.t))
    if not np.allclose(tot.t[:n], wp.t[:n]):
        raise ValueError("density_total and density_wp frames are not time-aligned")
    bath = tot.profiles[:n] - wp.profiles[:n]
    return AxialSeries(tot.t[:n], tot.z, bath - bath[0], radius_bohr,
                       "induced_bath")


# ---------------------------------------------------------------------------
# Projectile frame
# ---------------------------------------------------------------------------

def to_projectile_frame(series: AxialSeries, z_proj: np.ndarray,
                        box_length_z: float,
                        zeta: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Re-express each profile against zeta = z - z_proj(t), wrapped periodically.

    `z_proj` must already be sampled at `series.t`. Returns (zeta, values) with
    values[i, :] the profile at time i interpolated onto the common zeta grid.
    """
    if zeta is None:
        zeta = np.linspace(-box_length_z / 2, box_length_z / 2, 401)
    out = np.empty((len(series.t), len(zeta)))
    for i, zp in enumerate(z_proj):
        d = series.z - zp
        d = (d + box_length_z / 2) % box_length_z - box_length_z / 2   # wrap
        order = np.argsort(d)
        out[i] = np.interp(zeta, d[order], series.profiles[i][order])
    return zeta, out


def wake_asymmetry(zeta: np.ndarray, values: np.ndarray,
                   window_bohr: float = 12.0) -> dict[str, np.ndarray]:
    """Front/back imbalance of the induced density — the drag signature.

    Returns integrals of dn over [0, +w] (ahead) and [-w, 0] (behind) and their
    difference. A projectile that drags a lagging polarisation cloud has MORE
    induced charge behind it than ahead; the difference is what the induced field
    acts on, so its sign and size track the stopping power.
    """
    ahead = (zeta >= 0) & (zeta <= window_bohr)
    behind = (zeta <= 0) & (zeta >= -window_bohr)
    dz = float(np.mean(np.diff(zeta)))
    a = values[:, ahead].sum(axis=1) * dz
    b = values[:, behind].sum(axis=1) * dz
    return {"ahead": a, "behind": b, "asymmetry": b - a,
            "peak_depletion": values.min(axis=1),
            "peak_pileup": values.max(axis=1),
            "at_projectile": values[:, int(np.argmin(np.abs(zeta)))]}


# ---------------------------------------------------------------------------
# Energy channels
# ---------------------------------------------------------------------------

ENERGY_COLS = ["energy_kinetic", "energy_hartree", "energy_xc",
               "energy_external", "energy_nonlocal", "energy_ion",
               "energy_ion_kinetic", "energy_exact_exchange"]


def energy_channels(df, t0: float | None = None, t1: float | None = None):
    """Per-component energy CHANGE (eV) from the window start.

    The point of this table is to say WHERE the deposited energy goes, not how
    much there is. For the classical run `energy_total` is electronic-only (INQ
    leaves `energy_ion_kinetic` at zero — verified 2026-07-30), so the components
    should sum to a rise equal to the projectile's KE loss.
    """
    m = np.ones(len(df), dtype=bool)
    if t0 is not None:
        m &= df["time_au"].to_numpy() >= t0
    if t1 is not None:
        m &= df["time_au"].to_numpy() <= t1
    sub = df[m]
    out = {}
    for c in ENERGY_COLS + ["energy_total"]:
        if c in sub.columns:
            v = sub[c].to_numpy() * HA_TO_EV
            if np.ptp(v) == 0.0 and np.all(v == 0.0):
                continue                     # channel not populated by INQ
            out[c] = v - v[0]
    out["time_au"] = sub["time_au"].to_numpy()
    return out
