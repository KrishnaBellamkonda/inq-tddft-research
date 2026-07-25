"""PDE-FIND: sparse identification of a governing PDE from a field u(x,t).

Campaign-local kernel (ml-patterns, bulk-jellium PDE-discovery redo, T9).
FORMULA-BEARING -> pre-gated (formula-validation + code-test + catalogue) before
any headline use.

Method (cited):
  - Sparse PDE identification via STRidge (sequentially thresholded ridge):
    Rudy, Brunton, Proctor, Kutz, "Data-driven discovery of partial differential
    equations", Science Advances 3, e1602614 (2017). Build a candidate library
    Theta(u) of spatial operators and monomials, regress the time-derivative
    b = d^m u/dt^m onto Theta, and enforce sparsity by thresholding small
    coefficients and refitting (STRidge / STLSQ).
  - Base SINDy sparse-regression idea: Brunton, Proctor, Kutz, PNAS 113, 3932
    (2016).

Design for this campaign:
  - BROAD, AGNOSTIC library (minimal physics priors) — powers of u up to `poly`
    times spatial derivatives up to `deriv_order`, with cross products. Physical
    names of surviving terms are assigned POST-HOC, never seeded (ADR 0012).
  - Target order `m` selectable: m=1 -> d u/dt = Theta c ; m=2 -> d^2u/dt^2 =
    Theta c (plasma oscillation is intrinsically 2nd order).
  - THREE VALIDATION WALLS (ADR 0012) live here as: (2) `forward_integrate` +
    `forward_score` (temporal held-out prediction), (3) `bootstrap_stability`
    (term persistence under resampling). Wall (1) — the pinned calib/held-out
    CELL split — is enforced by the orchestrator, not this kernel.

Works on 1D spatial fields u of shape (T, Nx) (e.g. the axial reduction n(z,t)).
A 2D (r,z) extension shares `spatial_derivatives_nd`. Derivatives use central
finite differences on a Gaussian-smoothed field (noise control, Rudy 2017 SI);
interior points are subsampled to avoid edge artefacts.

Never np.fft.fftshift a VTI; fields arrive already in physical order.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from scipy.ndimage import gaussian_filter1d


# ---------------------------------------------------------------------------
# Derivatives
# ---------------------------------------------------------------------------
def smooth_field(u: np.ndarray, sigma_t: float = 0.0, sigma_x: float = 0.0) -> np.ndarray:
    """Gaussian-smooth u(T, Nx) in time and/or space (0 => no smoothing)."""
    out = u.astype(np.float64)
    if sigma_t and sigma_t > 0:
        out = gaussian_filter1d(out, sigma_t, axis=0, mode="nearest")
    if sigma_x and sigma_x > 0:
        out = gaussian_filter1d(out, sigma_x, axis=1, mode="nearest")
    return out


def d_dt(u: np.ndarray, dt: float, order: int = 1) -> np.ndarray:
    """Central finite-difference time derivative of u(T, Nx), given order 1 or 2."""
    if order == 1:
        return np.gradient(u, dt, axis=0, edge_order=2)
    if order == 2:
        # second derivative via central 2nd difference (uniform dt)
        d2 = np.empty_like(u)
        d2[1:-1] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / (dt * dt)
        d2[0] = d2[1]
        d2[-1] = d2[-2]
        return d2
    raise ValueError("order must be 1 or 2")


def d_dx(u: np.ndarray, dx: float, order: int) -> np.ndarray:
    """order-th spatial derivative of u(T, Nx) along axis=1 by repeated gradient."""
    out = u
    for _ in range(order):
        out = np.gradient(out, dx, axis=1, edge_order=2)
    return out


# ---------------------------------------------------------------------------
# Candidate library (broad, agnostic)
# ---------------------------------------------------------------------------
@dataclass
class Library:
    Theta: np.ndarray            # (n_samples, n_terms)
    names: list                  # length n_terms
    spatial: list                # per-term spatial-derivative order (for integration)
    poly: list                   # per-term polynomial power of u


def build_library_1d(u: np.ndarray, dx: float, poly: int = 3,
                     deriv_order: int = 3, smooth_x: float = 0.0,
                     include_const: bool = True):
    """Broad agnostic library for a 1D field u(T, Nx).

    Terms = u^p * d^d u/dx^d for p in 0..poly, d in 0..deriv_order (excluding the
    trivial p=0,d=0 constant unless include_const). Returns a Library plus the
    pre-computed derivative stack so forward-integration can reuse the operators.
    """
    us = smooth_field(u, sigma_x=smooth_x) if smooth_x else u.astype(np.float64)
    # deriv FACTOR: D_0 = 1 (no derivative), D_d = d^d u/dx^d for d >= 1.
    # library term(p, d) = u^p * D_d  (so p=0,d=0 -> the constant 1; p=1,d=0 -> u).
    derivs = {d: d_dx(us, dx, d) for d in range(1, deriv_order + 1)}
    cols, names, spat, pw = [], [], [], []
    for p in range(poly + 1):
        up = us ** p if p > 0 else np.ones_like(us)
        for d in range(deriv_order + 1):
            if p == 0 and d == 0 and not include_const:
                continue
            dfac = derivs[d] if d > 0 else 1.0
            term = up * dfac
            cols.append(term.reshape(-1))
            nm = _term_name(p, d)
            names.append(nm)
            spat.append(d)
            pw.append(p)
    Theta = np.stack(cols, axis=1)
    return Library(Theta=Theta, names=names, spatial=spat, poly=pw), derivs


def _term_name(p: int, d: int) -> str:
    base = "1" if p == 0 else ("u" if p == 1 else f"u^{p}")
    if d == 0:
        return base
    dname = {1: "u_x", 2: "u_xx", 3: "u_xxx"}.get(d, f"u_x{d}")
    if p == 0:
        return dname
    ub = "u" if p == 1 else f"u^{p}"
    return f"{ub}*{dname}"


# ---------------------------------------------------------------------------
# Sparse regression: STRidge (Rudy et al. 2017)
# ---------------------------------------------------------------------------
def _ridge(A, b, lam):
    n = A.shape[1]
    return np.linalg.lstsq(A.T @ A + lam * np.eye(n), A.T @ b, rcond=None)[0]


def stridge(Theta: np.ndarray, b: np.ndarray, threshold: float,
            lam: float = 1e-5, max_iter: int = 20,
            normalize: bool = True):
    """Sequentially thresholded ridge regression (STRidge).

    Returns coefficients in the ORIGINAL (un-normalized) physical scale.

    Support selection is done in a DIMENSIONLESS space: library columns are
    L2-normalized to unit norm AND the target b is normalized by its L2 norm, so
    the magnitude `threshold` is scale-invariant to both the library-term units
    and the overall amplitude of the field (crucial here — induced densities are
    ~1e-5, so an absolute threshold would cull every term). The final unbiased
    least-squares refit and de-normalization recover physical coefficients.
    """
    A = Theta.astype(np.float64)
    b = b.astype(np.float64).reshape(-1)
    if normalize:
        scale = np.linalg.norm(A, axis=0)
        scale[scale == 0] = 1.0
        An = A / scale
        bscale = np.linalg.norm(b) or 1.0
    else:
        scale = np.ones(A.shape[1])
        An = A
        bscale = 1.0
    bn = b / bscale
    w = _ridge(An, bn, lam)
    keep = np.abs(w) >= threshold
    for _ in range(max_iter):
        if keep.sum() == 0:
            break
        w = np.zeros(A.shape[1])
        w[keep] = _ridge(An[:, keep], bn, lam)
        new_keep = np.abs(w) >= threshold
        new_keep &= keep  # only ever remove terms (monotone)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    # final least-squares refit on kept terms (unbiased), then de-normalize
    w = np.zeros(A.shape[1])
    if keep.sum() > 0:
        w[keep] = np.linalg.lstsq(An[:, keep], bn, rcond=None)[0]
    return w * bscale / scale


# ---------------------------------------------------------------------------
# Sample selection (interior, subsampled)
# ---------------------------------------------------------------------------
def interior_mask(T: int, Nx: int, t_margin: int, x_margin: int) -> np.ndarray:
    m = np.zeros((T, Nx), dtype=bool)
    m[t_margin:T - t_margin, x_margin:Nx - x_margin] = True
    return m.reshape(-1)


# ---------------------------------------------------------------------------
# Main discovery entry point
# ---------------------------------------------------------------------------
@dataclass
class PDEModel:
    coeffs: np.ndarray
    names: list
    spatial: list
    poly: list
    order: int                    # time-derivative order (1 or 2)
    dx: float
    dt: float
    residual_rel: float
    active: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def pretty(self, sig=2) -> str:
        lhs = "u_t" if self.order == 1 else "u_tt"
        terms = [f"{c:.{sig}g}*{n}" for c, n in zip(self.coeffs, self.names)
                 if abs(c) > 0]
        return f"{lhs} = " + (" + ".join(terms) if terms else "0")


def discover_pde_1d(u: np.ndarray, dx: float, dt: float, order: int = 2,
                    poly: int = 3, deriv_order: int = 3, threshold: float = 0.02,
                    lam: float = 1e-5, smooth_t: float = 0.0, smooth_x: float = 0.0,
                    t_margin: int = 2, x_margin: int = 3,
                    subsample: int = 20000, seed: int = 0) -> PDEModel:
    """Discover d^order u/dt^order = Theta(u) c for a 1D field u(T, Nx)."""
    u = u.astype(np.float64)
    if smooth_t:
        u = smooth_field(u, sigma_t=smooth_t)
    T, Nx = u.shape
    lib, _ = build_library_1d(u, dx, poly=poly, deriv_order=deriv_order,
                              smooth_x=smooth_x)
    b = d_dt(u, dt, order=order).reshape(-1)
    mask = interior_mask(T, Nx, t_margin, x_margin)
    idx = np.where(mask)[0]
    rng = np.random.default_rng(seed)
    if subsample and idx.size > subsample:
        idx = rng.choice(idx, subsample, replace=False)
    Theta = lib.Theta[idx]
    bb = b[idx]
    w = stridge(Theta, bb, threshold=threshold, lam=lam)
    pred = Theta @ w
    denom = np.linalg.norm(bb) or 1.0
    resid = float(np.linalg.norm(bb - pred) / denom)
    active = [n for n, c in zip(lib.names, w) if abs(c) > 0]
    return PDEModel(coeffs=w, names=lib.names, spatial=lib.spatial, poly=lib.poly,
                    order=order, dx=dx, dt=dt, residual_rel=resid, active=active,
                    meta={"poly": poly, "deriv_order": deriv_order,
                          "threshold": threshold, "lam": lam, "n_samples": int(idx.size)})


# ---------------------------------------------------------------------------
# Wall (2): forward integration + score (temporal held-out prediction)
# ---------------------------------------------------------------------------
def _rhs_1d(u_row: np.ndarray, model: PDEModel) -> np.ndarray:
    """Evaluate Theta(u) c at a single time-row u_row (Nx,).

    Wrapped in errstate('ignore'): a stiff/unstable discovered PDE can overflow
    during forward integration; the resulting inf/nan is DETECTED by
    forward_integrate (isfinite guard) and forward_score (returns 1e3), so the
    forward-prediction wall rejects it cleanly — the warnings are just noise.
    """
    dx = model.dx
    out = np.zeros_like(u_row)
    # cache spatial derivatives up to max order present
    maxd = max(model.spatial) if model.spatial else 0
    derivs = {}
    cur = u_row
    for d in range(1, maxd + 1):
        cur = np.gradient(cur, dx, edge_order=2)
        derivs[d] = cur
    with np.errstate(over="ignore", invalid="ignore"):
        for c, d, p in zip(model.coeffs, model.spatial, model.poly):
            if c == 0:
                continue
            up = u_row ** p if p > 0 else 1.0
            dfac = derivs[d] if d > 0 else 1.0
            out = out + c * up * dfac
    return out


def forward_integrate(model: PDEModel, u0: np.ndarray, nsteps: int,
                      u0_prev: np.ndarray | None = None, nsub: int = 8) -> np.ndarray:
    """Integrate the discovered PDE forward from initial condition(s), returning
    one frame per outer step (nsteps+... frames).

    Each frame interval dt is sub-divided into `nsub` RK4 sub-steps (dt/nsub) so
    the explicit scheme stays inside its CFL limit for wave/diffusion terms —
    without sub-stepping, any dispersive (u_xxx) term blows up and the
    forward-prediction wall would reject good PDEs on purely numerical grounds.

    order=1: du/dt = RHS(u), RK4.
    order=2: recast as the stable first-order system du/dt = v, dv/dt = RHS(u),
    RK4; v0 estimated from (u0 - u0_prev)/dt.
    """
    dt = model.dt
    h = dt / nsub
    if model.order == 1:
        u = u0.copy()
        traj = [u.copy()]
        for _ in range(nsteps):
            for _ in range(nsub):
                k1 = _rhs_1d(u, model)
                k2 = _rhs_1d(u + 0.5 * h * k1, model)
                k3 = _rhs_1d(u + 0.5 * h * k2, model)
                k4 = _rhs_1d(u + h * k3, model)
                u = u + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
                if not np.isfinite(u).all():
                    traj.append(u.copy())
                    return np.stack(traj, axis=0)
            traj.append(u.copy())
        return np.stack(traj, axis=0)
    # order == 2 as (u, v) system, RK4 with sub-stepping
    if u0_prev is None:
        raise ValueError("order=2 forward integration needs u0_prev")
    u = u0.copy()
    v = (u0 - u0_prev) / dt
    traj = [u0_prev.copy(), u.copy()]
    for _ in range(nsteps - 1):
        for _ in range(nsub):
            a1 = _rhs_1d(u, model)
            u2, v2 = u + 0.5 * h * v, v + 0.5 * h * a1
            a2 = _rhs_1d(u2, model)
            u3, v3 = u + 0.5 * h * v2, v + 0.5 * h * a2
            a3 = _rhs_1d(u3, model)
            u4, v4 = u + h * v3, v + h * a3
            a4 = _rhs_1d(u4, model)
            u = u + (h / 6.0) * (v + 2 * v2 + 2 * v3 + v4)
            v = v + (h / 6.0) * (a1 + 2 * a2 + 2 * a3 + a4)
            if not np.isfinite(u).all():
                traj.append(u.copy())
                return np.stack(traj, axis=0)
        traj.append(u.copy())
    return np.stack(traj, axis=0)


def forward_score(model: PDEModel, u_true: np.ndarray, fit_frac: float = 0.5):
    """Fit-window forward-prediction score (Wall 2).

    Integrate from the boundary of the fit window and compare to the true later
    frames. Returns (rel_l2, predicted_stack). rel_l2 is the relative L2 error
    over the prediction window (lower is better; a governing PDE should predict
    dynamics it was not fit to).
    """
    T = u_true.shape[0]
    t0 = max(2, int(fit_frac * T))
    nsteps = T - t0
    if nsteps < 2:
        return float("nan"), None
    if model.order == 1:
        pred = forward_integrate(model, u_true[t0], nsteps)
        true_win = u_true[t0:t0 + pred.shape[0]]
    else:
        pred = forward_integrate(model, u_true[t0], nsteps, u0_prev=u_true[t0 - 1])
        true_win = u_true[t0 - 1:t0 - 1 + pred.shape[0]]
    n = min(pred.shape[0], true_win.shape[0])
    pred, true_win = pred[:n], true_win[:n]
    if not np.isfinite(pred).all():
        # unstable forward integration (stiff/blown-up PDE) -> the wall REJECTS it
        return 1e3, pred
    denom = np.linalg.norm(true_win) or 1.0
    rel = float(np.linalg.norm(pred - true_win) / denom)
    return rel, pred


# ---------------------------------------------------------------------------
# Wall (3): bootstrap coefficient stability
# ---------------------------------------------------------------------------
def bootstrap_stability(u: np.ndarray, dx: float, dt: float, order: int = 2,
                        n_boot: int = 20, frac: float = 0.6, seed: int = 0,
                        **kw) -> dict:
    """Refit on resampled temporal subsets; report per-term active fraction.

    A term is 'stable' if it stays active across resamples. Returns
    {name: fraction_active} plus median coefficient over resamples.
    """
    rng = np.random.default_rng(seed)
    T = u.shape[0]
    counts, coeff_acc = {}, {}
    k = max(4, int(frac * T))
    for bset in range(n_boot):
        rows = np.sort(rng.choice(T, k, replace=False))
        m = discover_pde_1d(u[rows], dx, dt, order=order, seed=bset, **kw)
        for c, n in zip(m.coeffs, m.names):
            coeff_acc.setdefault(n, []).append(c)
            counts[n] = counts.get(n, 0) + (1 if abs(c) > 0 else 0)
    frac_active = {n: counts.get(n, 0) / n_boot for n in coeff_acc}
    med_coeff = {n: float(np.median(v)) for n, v in coeff_acc.items()}
    return {"frac_active": frac_active, "median_coeff": med_coeff, "n_boot": n_boot}


# ---------------------------------------------------------------------------
# Post-hoc physics interpretation (rule-based, deterministic; ADR 0012)
# ---------------------------------------------------------------------------
INTERPRET = {
    "u": "restoring / local response (compare -omega_p^2 n)",
    "u_x": "advection / drift (compare v.grad n)",
    "u_xx": "diffusion / dispersion (compare Bohm-Gross 3 v_th^2 grad^2 n)",
    "u_xxx": "dispersive (3rd-order) transport",
    "u*u_x": "nonlinear advection (Burgers-type)",
    "u^2": "quadratic (Barkas / nonlinear response)",
    "1": "constant source / drift offset",
}


def mask_to_terms(model: PDEModel, keep_names) -> PDEModel:
    """Return a copy of `model` with all coefficients zeroed except `keep_names`.

    Used to build the bootstrap-ADMITTED equation (the one actually claimed) so
    the forward-prediction wall scores the parsimonious governing law, not the
    raw STRidge output whose high-order transient terms can be numerically stiff.
    """
    keep = set(keep_names)
    coeffs = np.array([c if n in keep else 0.0
                       for c, n in zip(model.coeffs, model.names)])
    m = PDEModel(coeffs=coeffs, names=model.names, spatial=model.spatial,
                 poly=model.poly, order=model.order, dx=model.dx, dt=model.dt,
                 residual_rel=model.residual_rel,
                 active=[n for c, n in zip(coeffs, model.names) if c != 0],
                 meta=dict(model.meta, masked_to=list(keep)))
    return m


def interpret(model: PDEModel, min_abs: float = 0.0) -> list:
    """Name surviving terms against known operators (post-hoc, never seeded)."""
    out = []
    for c, n in zip(model.coeffs, model.names):
        if abs(c) <= min_abs:
            continue
        out.append({"term": n, "coeff": float(c),
                    "physics": INTERPRET.get(n, "uninterpreted term")})
    return out
