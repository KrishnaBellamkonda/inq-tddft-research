"""Induced-density normaliser + subtraction ladder + VTI series loading.

Campaign-local kernel (ml-patterns). PRE-GATED in T1.

Subtraction ladder (research Q5; mandatory before ANY ML):
  1. GS          : dn = n_bath(t) - n_bath^GS         (remove static background)
  2. rigid-motion: |FFT| metric is translation-invariant, so for the q-space
                   ratio R(q) the common rigid projectile translation cancels in
                   magnitude; for real-space POD we additionally mean-subtract.
  3. Lindhard    : linear-response reference (reused from inqview.analysis.
                   lindhard_elf) — for the wake NONLINEAR residual (T3). Cancels
                   in the matched-v ratio R(q) (same chi), so T2 does not need it.
  4. vacuum-WP   : subtract no-bath WP run at matched (v,sigma) for n_wp-based
                   signatures (T6); requires a free_wp partner (skip+log if absent).

Bath density (project memory reference_canonical_bath_density):
  n_bath = n_total - n_wp.  In `_wf` runs density_system is ALREADY bath-only
  (sum ~ N_bath electrons); in legacy runs density_system == density_total
  (WP-included). The resolver below returns bath-only by preferring
  density_system when it integrates to N_bath, else density_total - density_wp.

Never np.fft.fftshift a VTI (physical order; load via inqview.load_vti).
float32; subsample to <= max_frames per series (campaign memory budget).
"""
from __future__ import annotations
import os, glob
import numpy as np
from dataclasses import dataclass, field
from inqview import load_vti


# ----------------------------------------------------------------------------
# VTI series loading
# ----------------------------------------------------------------------------
def list_vti(series_dir: str) -> list[str]:
    files = sorted(glob.glob(os.path.join(series_dir, "density_t*.vti")))
    return files


def subsample_indices(n: int, max_frames: int) -> np.ndarray:
    """Uniform-STRIDE subsample (constant spacing).

    A constant stride keeps the time spacing uniform so finite-difference time
    derivatives are artefact-free; linspace-round gives mixed 1/2 spacing that
    injects a sawtooth into u_t/u_tt (bulk PDE-discovery redo, 2026-07-03).
    """
    if n <= max_frames:
        return np.arange(n)
    stride = int(np.ceil(n / max_frames))
    return np.arange(0, n, stride)


def load_series(series_dir: str, max_frames: int = 300):
    """Load a VTI series -> (data (T, nx, ny, nz) float32, axes (x,y,z), times_idx).

    Loads via inqview.load_vti (physical order). Subsamples to <= max_frames.
    """
    files = list_vti(series_dir)
    if not files:
        raise FileNotFoundError(f"no VTI frames in {series_dir}")
    idx = subsample_indices(len(files), max_frames)
    sel = [files[i] for i in idx]
    first = load_vti(sel[0])
    nx, ny, nz = first.data.shape
    out = np.empty((len(sel), nx, ny, nz), dtype=np.float32)
    out[0] = first.data.astype(np.float32)
    for k, f in enumerate(sel[1:], start=1):
        out[k] = load_vti(f).data.astype(np.float32)
    axes = (np.asarray(first.x), np.asarray(first.y), np.asarray(first.z))
    return out, axes, idx


def load_gs(gs_path: str) -> np.ndarray:
    return load_vti(gs_path).data.astype(np.float32)


# ----------------------------------------------------------------------------
# Co-gridding (resample two fields to a common grid before subtraction)
# ----------------------------------------------------------------------------
def grids_match(ax_a, ax_b, atol=1e-6) -> bool:
    return all(
        a.shape == b.shape and np.allclose(a, b, atol=atol)
        for a, b in zip(ax_a, ax_b)
    )


def cogrid(data, axes_src, axes_dst):
    """Trilinear-resample `data` (..., nx,ny,nz) from axes_src onto axes_dst.

    Identity (returns input) when grids already match. Uses scipy
    RegularGridInterpolator per frame. Returns resampled array.
    """
    if grids_match(axes_src, axes_dst):
        return data
    from scipy.interpolate import RegularGridInterpolator
    xs, ys, zs = axes_src
    xd, yd, zd = axes_dst
    XD, YD, ZD = np.meshgrid(xd, yd, zd, indexing="ij")
    pts = np.stack([XD.ravel(), YD.ravel(), ZD.ravel()], axis=-1)
    single = data.ndim == 3
    frames = data[None] if single else data
    out = np.empty((frames.shape[0], len(xd), len(yd), len(zd)), dtype=np.float32)
    for k in range(frames.shape[0]):
        itp = RegularGridInterpolator((xs, ys, zs), frames[k], bounds_error=False,
                                      fill_value=0.0)
        out[k] = itp(pts).reshape(len(xd), len(yd), len(zd)).astype(np.float32)
    return out[0] if single else out


# ----------------------------------------------------------------------------
# Subtraction ladder
# ----------------------------------------------------------------------------
@dataclass
class InducedSeries:
    delta: np.ndarray              # (T, nx, ny, nz) dn = n_bath - n_bath_GS
    axes: tuple
    dx: float
    rungs_applied: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def induced_from_bath(bath_series: np.ndarray, gs_bath: np.ndarray, axes,
                      dx: float) -> InducedSeries:
    """Rung 1: subtract the GS bath from every frame."""
    delta = (bath_series - gs_bath[None]).astype(np.float32)
    return InducedSeries(delta=delta, axes=axes, dx=dx, rungs_applied=["GS"])


def remove_mean(ind: InducedSeries) -> InducedSeries:
    """Rung 2 (real-space): remove the time-mean (common rigid component proxy)."""
    ind.delta = (ind.delta - ind.delta.mean(axis=0, keepdims=True)).astype(np.float32)
    ind.rungs_applied.append("mean/rigid")
    return ind


def time_reduce_q(delta: np.ndarray, dx: float, reducer, **kw):
    """Apply a radial-q spectrum reducer per frame, then median over time.

    reducer(frame, dx, **kw) -> (q, amp). Returns (q, robust_amp) where the
    robust amplitude is the per-q median over frames (temporal reduction).
    """
    qs = None
    amps = []
    for k in range(delta.shape[0]):
        q, a = reducer(delta[k], dx, **kw)
        qs = q
        amps.append(a)
    amps = np.asarray(amps)
    return qs, np.median(amps, axis=0)
