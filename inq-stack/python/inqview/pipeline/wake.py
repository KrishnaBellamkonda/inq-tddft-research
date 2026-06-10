#!/usr/bin/env python3
"""inqview.postprocess.wake — canonical bath-density wake extraction + movies.

CANONICAL BATH DENSITY (use this everywhere; do not improvise)
--------------------------------------------------------------
The jellium "system"/bath density is

    n_system(r,t) = n_total(r,t) - n_wp(r,t)

i.e. the full electronic density MINUS the injected wave-packet orbital
(occupation 1.0).  The *induced* bath wake is the t0-subtracted bath:

    dn_system(r,t) = n_system(r,t) - n_system(r,t0)

WHY total - wp AND NOT the saved `density_system` field
-------------------------------------------------------
The saved `density_system` VTI field is NOT consistent across run
generations (verified 2026-06-01 by integration):
  * old runs (e.g. run_wp_n162_L50_E100_sigma1_v2): density_system = 163 e
    -> it INCLUDES the WP orbital (NOT the bath).
  * new `_wf` runs (run_template fix): density_system = 162 e (bath only).
Both have n_total = 163 e (162 bath + 1 WP) and n_wp = 1.000 e at all times.
Therefore the only run-independent bath density is `total - wp`, computed in
post-processing.  This module always does that.

CLASSICAL runs have no WP orbital in the electron density (the projectile is
an external ion potential), so for them n_system = n_total directly (no
density_wp field is written).  `bath_*` functions below detect this and skip
the subtraction.

SHARED-COLORBAR RULE (project rule, user-mandated 2026-06-01)
-------------------------------------------------------------
Panels that are meant to be compared DIRECTLY (e.g. WP wake vs classical wake)
MUST share an identical colour scale.  Use `shared_clim()` to compute one
(vmin,vmax) over all such panels and pass it to every one of them.  A
difference panel (WP - classical) is evaluated on its OWN scale (independent
quantity).  Provide BOTH linear and log-scale views.

Verified known-cases (2026-06-01, dev-feedback-loop rule):
  * integral(n_total)=163, integral(n_wp)=1.000, total-wp=162=N_e (t0/mid/late).
  * dn_system at t0 == 0 by construction.

RECIPE — induced-wake difference plots (the canonical 3 steps)
--------------------------------------------------------------
Step 1  n_system at each timestep, EXACT same-timeframe subtraction:
        n_system(r,t) = n_total(r,t) - n_wp(r,t), both at the SAME exact step
        (`bath_volume` snaps to a frame with an exact density_wp partner and
        raises otherwise). NEVER subtract n_wp from a nearest/different step —
        the WP is a moving Gaussian, so a mismatched subtraction leaves a
        charge-neutral MOVING-WP DIPOLE residual (looks like the WP is still
        present; was ~85% of the sigma=1 signal, 2026-06-01). For sparse
        density_wp (v2 runs) sample at `wp_frame_times(run)`. Classical run:
        no WP orbital -> n_system = n_total. Verify ∫n_system dV = N_electrons.
Step 2  induced wake: dn_system(r,t) = n_system(r,t) - n_system(r,t0)
        (`bath_line_z`/`bath_slice_xz` then subtract the t0 result;
        dn_system(t0)=0 by construction).
Step 3  cross-system metric difference at a MATCHED time t:
        D(r,t) = dn_system^A(r,t) - dn_system^B(r,t)   (e.g. A=WP, B=classical),
        BOTH evaluated at the same physical t. Show A and B in panels sharing
        ONE colour scale (`shared_clim` over both); the difference D gets its
        OWN scale. Mark the WP centroid (`wp_centroid_z`). Linear AND symlog;
        fix the scale ONCE over all animation frames (never per-frame).
Reference impl: ResearchProject/systems/jellium/scripts/wake_movie_driver.py
"""
from __future__ import annotations
import glob, re
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

HA2EV = 27.211386245988

# TODO: postprocess/ submodule, I don't know if its making sense. I need to come up with
# a simple and effective organisational structure for the python library. 

# TODO: Eventually, when this behaviour of total - wp = system the same as teh system 
# output, in the inqkit side of the library, then, we should can reduce the
# complexity of this file. 

# TODO: It would be a good idea to use the center of density calculated in the inqkit
# library and visualise it. This way we can compare the wave packet centre of density,
# the total electronic system, and the jellium bath system. 

# ---------------------------------------------------------------- VTI loading
def _read_vti(path):
    """Return (array3d[nz,ny,nx], origin(x,y,z), spacing(x,y,z))."""
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(str(path)); r.Update()
    img = r.GetOutput()
    nx, ny, nz = img.GetDimensions()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return a, img.GetOrigin(), img.GetSpacing()


def frames(run, sub):
    """Sorted [(step, path), ...] for results/raw/vti/<sub>/*.vti."""
    out = []
    base = Path(run) / "results/raw/vti" / sub
    if not base.is_dir():
        base = Path(run) / "results/vti" / sub
    for f in glob.glob(f"{base}/*.vti"):
        m = re.search(r"_t(\d+)\.vti$", f)
        if m:
            out.append((int(m.group(1)), f))
    return sorted(out)


def dt_of(run):
    txt = (Path(run) / "results/run_summary.txt").read_text()
    m = re.search(r"dt_au\s*=\s*([\d.]+)", txt)
    return float(m.group(1)) if m else 0.02


def has_wp(run):
    return len(frames(run, "density_wp")) > 0


# ------------------------------------------------------- bath density (t-pick)
def _pick(steps, dt, t):
    return int(np.argmin(np.abs(np.asarray(steps) * dt - t)))


def _aligned_wp(run, target_step, exact=True):
    """Path of the density_wp frame at target_step. If exact and there is no
    EXACT-step match, return None (do NOT subtract a temporally-displaced WP —
    a moving Gaussian subtracted at the wrong step leaves a dipole residual)."""
    fw = frames(run, "density_wp")
    if not fw:
        return None
    ws = np.array([s for s, _ in fw])
    j = int(np.argmin(np.abs(ws - target_step)))
    if exact and ws[j] != target_step:
        return None
    return fw[j][1]


def wp_frame_times(run, dt=None):
    """Times (a.u.) at which density_wp EXISTS — use these to sample WP runs so
    the n_total - n_wp subtraction is exact (no moving-WP residual). Empty for
    classical runs."""
    dt = dt or dt_of(run)
    return [s * dt for s, _ in frames(run, "density_wp")]


def bath_volume(run, t, dt=None, require_exact_wp=True):
    """3D bath density n_system = n_total - n_wp (or n_total for classical).

    The density_total frame nearest t is chosen; the WP is subtracted ONLY from
    an EXACT-step density_wp frame. If require_exact_wp and the total frame has
    no exact WP partner, raises — callers must sample at `wp_frame_times`.

    Returns (n_system[nz,ny,nx], origin, spacing, actual_time_au, step).
    """
    dt = dt or dt_of(run)
    if has_wp(run):
        # snap to the nearest density_total frame that HAS an exact WP partner
        wpt = set(s for s, _ in frames(run, "density_wp"))
        ft = [(s, p) for s, p in frames(run, "density_total") if s in wpt]
        if not ft:
            ft = frames(run, "density_total")
        steps = [s for s, _ in ft]
        i = _pick(steps, dt, t); st = steps[i]
        ntot, origin, spacing = _read_vti(ft[i][1])
        wp_path = _aligned_wp(run, st, exact=require_exact_wp)
        if wp_path is None:
            if require_exact_wp:
                raise ValueError(f"{run}: no exact density_wp at step {st}")
            nsys = ntot
        else:
            nwp, _, _ = _read_vti(wp_path)
            nsys = ntot - nwp
    else:
        ft = frames(run, "density_total")
        steps = [s for s, _ in ft]
        i = _pick(steps, dt, t); st = steps[i]
        ntot, origin, spacing = _read_vti(ft[i][1])
        nsys = ntot                       # classical: no WP orbital in e-density
    return nsys, origin, spacing, st * dt, st


def bath_line_z(run, t, dt=None):
    """1D bath z-profile (sum over x,y * dx*dy), e/Bohr.  Returns (z, line, t_au)."""
    nsys, origin, spacing, t_au, _ = bath_volume(run, t, dt)
    nz = nsys.shape[0]
    z = origin[2] + spacing[2] * np.arange(nz)
    line = nsys.sum(axis=(1, 2)) * spacing[0] * spacing[1]
    return z, line, t_au


def bath_slice_xz(run, t, dt=None):
    """2D xz bath slice at the central y plane, e/Bohr^3.  Returns (x, z, slab, t_au)."""
    nsys, origin, spacing, t_au, _ = bath_volume(run, t, dt)
    nz, ny, nx = nsys.shape
    slab = nsys[:, ny // 2, :]                       # [nz, nx]
    x = origin[0] + spacing[0] * np.arange(nx)
    z = origin[2] + spacing[2] * np.arange(nz)
    return x, z, slab, t_au


# ----------------------------------------------------------------- WP centroid
def wp_centroid_z(run, t, dt=None):
    """WP centroid z(t) in Bohr from the density_wp first moment (None if no WP).
    Picks the nearest density_wp frame directly (exact by construction)."""
    if not has_wp(run):
        return None
    dt = dt or dt_of(run)
    fw = frames(run, "density_wp")
    steps = [s for s, _ in fw]
    nwp, origin, spacing = _read_vti(fw[_pick(steps, dt, t)][1])
    nz = nwp.shape[0]
    z = origin[2] + spacing[2] * np.arange(nz)
    wz = nwp.sum(axis=(1, 2))
    wz = np.clip(wz, 0, None)
    s = wz.sum()
    return float((z * wz).sum() / s) if s > 0 else None


# ------------------------------------------------------------ colorbar helpers
def shared_clim(*arrays, symmetric=True, pct=100.0):
    """One (vmin,vmax) over ALL arrays — for directly-compared panels.

    symmetric=True -> (-m, m) about zero (signed Δn).  pct<100 clips to a
    percentile to suppress lone spikes (e.g. WP self-spike near boundary).
    """
    m = 0.0
    for a in arrays:
        a = np.asarray(a)
        v = np.percentile(np.abs(a), pct) if pct < 100 else np.abs(a).max()
        m = max(m, float(v))
    return (-m, m) if symmetric else (0.0, m)
