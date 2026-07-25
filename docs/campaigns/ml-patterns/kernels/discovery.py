"""Track-B governing-PDE discovery on a real bulk-jellium cell.

Campaign-local (ml-patterns redo, T11/T12). Bridges the run database + VTI series
to the PDE-FIND kernel:
    cell -> bath series -> GS-subtracted induced field dn(r,t)
         -> AXIAL reduction n(z,t) along the beam axis (z)
         -> discover_pde_1d + the two in-kernel walls (forward-predict, bootstrap)

The pinned calibration/held-out CELL split (Wall 1) is enforced by the
orchestrator; this module discovers on ONE cell and reports the per-cell walls.

Axial reduction: the projectile travels along z and the induced wake is dominated
by its axial structure, so we reduce the transverse (x,y) plane by MEAN to get a
robust 1D field n(z,t). (Mean, not sum, keeps the amplitude O(1); a constant
prefactor is absorbed into the discovered coefficients and does not change the
PDE structure.) A 2D (r,z) extension can reuse the same kernel later.

Never np.fft.fftshift a VTI; series load via inqview.load_vti (normaliser).
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass

from . import normaliser as NRM
from . import pdefind as PF


# ---------------------------------------------------------------------------
def axial_reduce(delta: np.ndarray, axes) -> tuple[np.ndarray, float]:
    """(T,nx,ny,nz) induced field -> (T,nz) axial profile + dz.

    Transverse mean over (x,y). z is the last axis (physical order; axis index 2
    in the (x,y,z) convention used across inqkit VTIs).
    """
    u = delta.mean(axis=(1, 2)).astype(np.float64)          # (T, nz)
    z = np.asarray(axes[2])
    dz = float(z[1] - z[0])
    return u, dz


@dataclass
class CellField:
    u: np.ndarray          # (T, Nz) axial induced field
    dz: float
    dt: float              # frame_dt (a.u.)
    n_frames: int
    which: str             # 'cl' or 'wp'
    note: str = ""


def load_cell_axial(cell: dict, which: str, max_frames: int = 200) -> CellField:
    """Load a cell's bath series, GS-subtract, axial-reduce.

    which='cl' -> classical partner; 'wp' -> wavepacket. Uses the bath-only dir
    resolved by celldb (density_system when bath-only, else total). GS from the
    matching density_gs_system.vti.
    """
    bdir = cell[f"{which}_bath_dir"]
    gs = cell[f"{which}_gs"]
    dt = cell.get(f"frame_dt_au_{which}") or cell.get("frame_dt_au_wp") or 0.04
    series, axes, idx = NRM.load_series(bdir, max_frames=max_frames)
    gs_bath = NRM.load_gs(gs)
    # co-grid GS onto the series grid if needed (single-run: usually identical)
    if gs_bath.shape != series.shape[1:]:
        gs_bath = NRM.cogrid(gs_bath, axes, axes)  # identity unless axes differ
    ind = NRM.induced_from_bath(series, gs_bath, axes, dx=float(axes[0][1] - axes[0][0]))
    u, dz = axial_reduce(ind.delta, axes)
    # effective frame spacing after subsampling (idx are original frame indices)
    if len(idx) > 1:
        dt = dt * float(np.mean(np.diff(idx)))
    return CellField(u=u, dz=dz, dt=float(dt), n_frames=u.shape[0], which=which)


# ---------------------------------------------------------------------------
DEFAULT_CFG = dict(order=2, poly=3, deriv_order=3, threshold=0.05,
                   lam=1e-5, smooth_t=1.0, smooth_x=1.0, x_margin=4, t_margin=3,
                   pod_rank=8)


def pod_denoise(u: np.ndarray, rank: int) -> tuple[np.ndarray, float]:
    """SVD-truncate the (T,Nz) field to `rank` spatiotemporal modes (denoise).

    Returns (reconstructed field, retained energy fraction). The time-mean is
    kept out of the SVD so rank counts fluctuation modes. This is the
    'latent-ODE supports field-PDE' denoiser: PDE-FIND on a POD-truncated field
    is standard practice to control derivative noise (Rudy 2017 SI).
    """
    if not rank or rank <= 0:
        return u, 1.0
    m = u.mean(0, keepdims=True)
    U, S, Vt = np.linalg.svd(u - m, full_matrices=False)
    r = min(rank, S.size)
    rec = (U[:, :r] * S[:r]) @ Vt[:r] + m
    ef = float((S[:r] ** 2).sum() / (S ** 2).sum()) if S.size else 1.0
    return rec, ef


def discover_cell(cf: CellField, cfg: dict | None = None,
                  bootstrap: int = 15, fwd_frac: float = 0.5) -> dict:
    """Discover the axial PDE for one cell + its two in-kernel walls.

    Returns a dict with the model equation, active terms, forward-prediction
    rel-L2 (Wall 2), bootstrap term-stability (Wall 3), and post-hoc physics
    interpretation. Skip+log if the field is too short / non-finite.
    """
    c = dict(DEFAULT_CFG)
    if cfg:
        c.update(cfg)
    pod_rank = c.pop("pod_rank", 0)
    u = cf.u
    if u.shape[0] < 8 or not np.isfinite(u).all():
        return {"ok": False, "reason": f"field too short/nonfinite (T={u.shape[0]})",
                "which": cf.which}
    u, energy_frac = pod_denoise(u, pod_rank)
    # nondimensionalise the field to unit std so nonlinear library terms (u^2,
    # u*u_x, ...) compete on the same scale as linear ones (SINDy assumes O(1)
    # variables). Coefficients are reported in NONDIMENSIONAL field units;
    # structure + the classical-vs-WP comparison are scale-consistent.
    field_scale = float(u.std()) or 1.0
    u = u / field_scale
    model = PF.discover_pde_1d(u, cf.dz, cf.dt, **c)
    bs = PF.bootstrap_stability(u, cf.dz, cf.dt, order=c["order"], n_boot=bootstrap,
                                poly=c["poly"], deriv_order=c["deriv_order"],
                                threshold=c["threshold"], lam=c["lam"],
                                x_margin=c["x_margin"], t_margin=c["t_margin"])
    interp = PF.interpret(model)
    # Wall 3: a term is "admitted" if active AND bootstrap-stable (>=0.6 of resamples)
    stable = {k for k, v in bs["frac_active"].items() if v >= 0.6}
    admitted = [t for t in interp if t["term"] in stable]
    # Wall 2: forward-predict the ADMITTED (parsimonious) equation, not the raw
    # STRidge output (whose stiff high-order transient terms blow up explicit
    # integration). Compose the walls: we validate the equation we claim.
    admitted_model = PF.mask_to_terms(model, stable)
    fwd_rel_full, _ = PF.forward_score(model, u, fit_frac=fwd_frac)
    fwd_rel_admit, _ = (PF.forward_score(admitted_model, u, fit_frac=fwd_frac)
                        if stable else (float("nan"), None))
    return {
        "ok": True, "which": cf.which,
        "equation": model.pretty(),
        "admitted_equation": admitted_model.pretty(),
        "order": c["order"], "dz": cf.dz, "dt": cf.dt, "n_frames": cf.n_frames,
        "pod_rank": pod_rank, "pod_energy_frac": energy_frac,
        "field_scale": field_scale,
        "residual_rel": model.residual_rel,
        "forward_rel_l2": fwd_rel_admit,       # the wall verdict (admitted eq)
        "forward_rel_l2_full": fwd_rel_full,   # raw STRidge (diagnostic)
        "active_terms": model.active,
        "bootstrap_active": bs["frac_active"],
        "admitted": admitted,
        "interpretation": interp,
        "config": c,
    }


def coeff_vector(result: dict, term_order: list) -> np.ndarray:
    """Extract admitted-term coefficients aligned to a fixed term ordering
    (for classical-vs-WP comparison, T13)."""
    d = {t["term"]: t["coeff"] for t in result.get("admitted", [])}
    return np.array([d.get(t, 0.0) for t in term_order])
