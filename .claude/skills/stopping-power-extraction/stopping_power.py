"""Electronic stopping-power extraction kernels (skill-local, shippable).

Two run-geometry-specific methods, plus shared diagnostics. All functions take
plain numpy arrays (loading is the caller's job) so the module is portable.

Grounding: Correa 2018, *Comput. Mater. Sci.* 150, 291 — Eq.(10) S=<dE/dt>/v
(equivalently dE/dx); transient/steady-state Fig.8. See
docs/sources/correa-2018-electronic-stopping-power.md.

METHOD A  (continuous-traversal runs: bulk jellium, projectile ploughs through a
           homogeneous medium):
    S = slope of dE_total(x) over a window [x0, xT], FREE-INTERCEPT fit.
        x0 = transient cutoff (fixed-fraction OR agent slope-plateau).
        xT = upper bound (default x_max; reduce to exclude post-exit / periodic
             image re-entry -> the endpoint check).

METHOD B  (localised slab: projectile enters -> deposits -> exits a finite slab):
    S = [E_total(t_final) - E_total(t_0)] / L_z,  L_z = slab thickness.
        REQUIRES a convergence gate: E_total(t_final) must have settled (deposit
        complete), else the slab transit isn't finished and S is an underestimate.

DIAGNOSTICS (run every time; flag large deviations to the user):
    - projectile  dKE  channel (energy conservation cross-check);
    - int F.v cumulative (== dKE when F=m dv/dt -> a profile/consistency check,
      NOT an independent physics channel);
    - N(t) ~ const guard: if a CAP drains the bath, raw E_total is NOT the
      deposited energy and BOTH methods are invalid.

Units: keep E and x in one consistent system (Ha/Bohr, or eV & Bohr -> eV/Bohr).
"""
from __future__ import annotations
import numpy as np

HA_TO_EV = 27.211386245988
HA_PER_BOHR_TO_EV_PER_A = 27.211386245988 / 0.529177210903


# --------------------------------------------------------------------------- fit
def free_fit(x, E, x0, xT):
    """Free-intercept least squares dE = S*x + E0 over [x0, xT].

    Free intercept is mandatory: the transient deposits a fixed energy E0; forcing
    the line through the origin biases S high. Returns None if <5 points.
    """
    x = np.asarray(x, float); E = np.asarray(E, float)
    m = (x >= x0) & (x <= xT); n = int(m.sum())
    if n < 5:
        return None
    A = np.vstack([x[m], np.ones(n)]).T
    (S, c), *_ = np.linalg.lstsq(A, E[m], rcond=None)
    r = E[m] - (S * x[m] + c)
    dof = max(n - 2, 1)
    sxx = np.sum((x[m] - x[m].mean()) ** 2)
    se = np.sqrt(np.sum(r ** 2) / dof / sxx) if sxx > 0 else np.inf
    sst = np.sum((E[m] - E[m].mean()) ** 2)
    r2 = 1.0 - np.sum(r ** 2) / sst if sst > 0 else np.nan
    return dict(S=float(S), E0=float(c), se=float(se), r2=float(r2), n=n,
                x0=float(x0), xT=float(xT), resid=r)


# ------------------------------------------------------------------- METHOD A (1)
def fixed_fraction_window(x, E, frac=0.20, xT=None):
    """Method A, fixed-fraction transient cut in x: x0 = x_min + frac*(x_max-x_min)."""
    x = np.asarray(x, float)
    xT = x.max() if xT is None else xT
    f = free_fit(x, E, x.min() + frac * (x.max() - x.min()), xT)
    if f is None:
        return dict(status="range_too_short")
    f["status"] = "ok"; f["frac"] = frac; f["basis"] = "x"
    return f


def fixed_time_fraction(t, x, E, frac=0.20, xT=None):
    """Method A DEFAULT for continuous/bulk runs: discard the first `frac` of the
    SIMULATION TIME as transient, then free-intercept fit dE(x) over the remainder.

    This is the locked default (user, 2026-06-25) for every NON-slab run: a single,
    deterministic rule — no plateau search, no tuning. `t`, `x` (projectile
    displacement) and `E` (dE_total) are aligned arrays. The cut is on TIME
    (t0 = t_min + frac*(t_max - t_min)); x0 is the displacement reached at t0. For a
    ~constant-velocity bulk run this is ~frac of the path too, but cutting on time is
    the literal rule. S = slope, error bar = regression standard error.

    The stricter `detect_x0_and_stopping_power` (agent slope-plateau + endpoint check)
    remains available as an OPTIONAL diagnostic, not the default.
    """
    t = np.asarray(t, float); x = np.asarray(x, float); E = np.asarray(E, float)
    o = np.argsort(t); t, x, E = t[o], x[o], E[o]
    t_cut = t.min() + frac * (t.max() - t.min())
    x0 = float(np.interp(t_cut, t, x))
    xT = x.max() if xT is None else xT
    f = free_fit(x, E, x0, xT)
    if f is None:
        return dict(status="range_too_short")
    f["status"] = "ok"; f["frac"] = frac; f["t_cut"] = float(t_cut); f["basis"] = "time"
    return f


# ------------------------------------------------------------------- METHOD A (2)
def detect_x0_and_stopping_power(x, E, xT=None, remain_min=0.30, rel_tol=0.02,
                                 k_sigma=2.0, gate=0.40, grid_step=None):
    """Method A, agent slope-plateau detector + endpoint check + 40% gate.

    Sweeps x0; x0 = smallest start where S(x0) enters a tolerance band and stays.
    `endpoint_status='endpoint_contaminated'` when the late-window slope departs
    from the mid-window slope (post-exit flattening / image re-entry) -> lower xT.
    Returns status in {ok, no_plateau, range_too_short, endpoint_contaminated}.
    """
    x = np.asarray(x, float); E = np.asarray(E, float)
    o = np.argsort(x); x, E = x[o], E[o]
    x_min, x_max = x.min(), x.max(); L = x_max - x_min
    if xT is None: xT = x_max
    if grid_step is None: grid_step = max(np.median(np.diff(x)), L / 400)
    ep = "ok"
    end = free_fit(x, E, 0.7 * (xT - x_min) + x_min, xT)
    mid = free_fit(x, E, x_min + 0.4 * (xT - x_min), x_min + 0.7 * (xT - x_min))
    if end and mid and abs(end["S"] - mid["S"]) > rel_tol * abs(mid["S"]) + 3 * mid["se"]:
        ep = "endpoint_contaminated"
    x0cap = x_min + (1 - remain_min) * (xT - x_min)
    grid = np.arange(x_min, x0cap, grid_step)
    fits = [free_fit(x, E, g, xT) for g in grid]
    keep = [f is not None for f in fits]
    grid = grid[keep]; fits = [f for f, k in zip(fits, keep) if k]
    if len(grid) < 3:
        return dict(status="range_too_short", x0=None, endpoint_status=ep, grid=grid)
    S = np.array([f["S"] for f in fits]); serr = np.array([f["se"] for f in fits])
    rem = (xT - grid) / (xT - x_min); ref = (rem >= 0.30) & (rem <= 0.55)
    if ref.sum() < 3: ref = np.ones_like(rem, bool)
    Spl = float(np.median(S[ref]))
    tol = max(rel_tol * abs(Spl), k_sigma * float(np.median(serr[ref])))
    ok = np.abs(S - Spl) <= tol
    x0 = None
    for i in range(len(grid)):
        if ok[i] and ok[i:].all():
            x0 = float(grid[i]); break
    if x0 is None:
        return dict(status="no_plateau", x0=None, Spl=Spl, tol=float(tol),
                    endpoint_status=ep, grid=grid, Sgrid=S, segrid=serr)
    f = free_fit(x, E, x0, xT)
    tf = (x0 - x_min) / L
    f.update(status=(ep if ep != "ok" else ("ok" if tf <= gate else "range_too_short")),
             endpoint_status=ep, x0=x0, Spl=Spl, tol=float(tol),
             transient_fraction=float(tf), grid=grid, Sgrid=S, segrid=serr)
    return f


# --------------------------------------------------------------------- METHOD B
def slab_stopping_power(t, E_total, L_z, *, converge_frac=0.15, converge_tol=0.05):
    """Method B (localised slab): S = [E_total(t_final) - E_total(t_0)] / L_z.

    L_z = slab thickness = traversal length (e.g. 25 Bohr). E_total and L_z in a
    consistent unit -> S in (E-unit)/Bohr.

    CONVERGENCE GATE (mandatory): the deposit must be complete at t_final. We
    require the energy change over the final `converge_frac` of the run to be
    <= `converge_tol` of the total deposit. If not, status='not_converged': the
    projectile has not finished depositing (extend the run); the reported S is a
    LOWER BOUND on the converged value.
    """
    t = np.asarray(t, float); E = np.asarray(E_total, float)
    o = np.argsort(t); t, E = t[o], E[o]
    dE = float(E[-1] - E[0])
    tail = t >= (t[-1] - converge_frac * (t[-1] - t[0]))
    tail_change = float(E[tail][-1] - E[tail][0])
    tail_frac = tail_change / dE if dE != 0 else np.inf
    converged = abs(tail_change) <= converge_tol * abs(dE)
    return dict(S=dE / L_z, dE=dE, L_z=float(L_z), converged=bool(converged),
                tail_change=tail_change, tail_frac_of_total=float(tail_frac),
                converge_frac=converge_frac, converge_tol=converge_tol,
                status="ok" if converged else "not_converged")


# -------------------------------------------------------------------- diagnostics
def conservation_guard(N_t, tol=0.02):
    """N(t) ~ const guard. If a CAP drains the bath, raw E_total is dominated by
    CAP energy and is NOT the deposited energy -> both methods invalid.
    Returns ok=True iff |dN|/N0 <= tol."""
    N = np.asarray(N_t, float)
    drained = abs(N[-1] - N[0]) / abs(N[0]) if N[0] != 0 else np.inf
    return dict(ok=bool(drained <= tol), drained_frac=float(drained), tol=tol)


def kinetic_channel(s_track, ke_track, x_eval, x0, xT):
    """Sanity (a): projectile -dKE/dx over [x0, xT], from the track (independent
    of the electronic E_total channel; their agreement is energy conservation)."""
    dKE_loss = ke_track[0] - np.interp(x_eval, s_track, ke_track)
    return free_fit(x_eval, dKE_loss, x0, xT)


def force_power_channel(t_track, v_track, mass=1.0):
    """Sanity (b): cumulative int(-F.v)dt with F=m*dv/dt, on the track time grid.
    Returns (t, deposited_energy). CAVEAT: equals dKE analytically (F=m dv/dt) ->
    a deposition-profile/discretisation check, NOT a third independent physics
    channel. The caller interpolates deposited_energy onto displacement to fit."""
    t = np.asarray(t_track, float); v = np.asarray(v_track, float)
    a = np.gradient(v, t)
    P = mass * a * v                                  # dKE/dt (<0 while decelerating)
    work = np.concatenate([[0.0], np.cumsum(0.5 * (P[1:] + P[:-1]) * np.diff(t))])
    return t, -work                                   # energy deposited to the medium


# ------------------------------------------------------------------------- test
def _selftest():
    rng = np.linspace(0, 10, 60)
    # free_fit recovers a known slope+intercept
    E = 2.0 * rng + 5.0
    f = free_fit(rng, E, 0, 10)
    assert abs(f["S"] - 2.0) < 1e-9 and abs(f["E0"] - 5.0) < 1e-9, f
    # fixed-fraction window on a clean line
    assert abs(fixed_fraction_window(rng, E, 0.2)["S"] - 2.0) < 1e-9
    # fixed-TIME-fraction (the bulk default): cut first 20% of time, recover slope
    tt = np.linspace(0, 5, 60); xx = 2.0 * tt; EE = 2.0 * xx + 5.0
    ft = fixed_time_fraction(tt, xx, EE, 0.2)
    assert abs(ft["S"] - 2.0) < 1e-9 and ft["basis"] == "time", ft
    assert abs(ft["t_cut"] - 1.0) < 1e-9, ft        # 20% of [0,5] = 1.0
    # slab: ramp (deposit) then perfect plateau -> converged, S = dE/L_z
    t = np.linspace(0, 20, 100)
    Eslab = np.where(t < 10, t, 10.0)          # rises to 10 then flat
    r = slab_stopping_power(t, Eslab, L_z=25.0)
    assert r["converged"] and abs(r["S"] - 10.0 / 25.0) < 1e-9, r
    # slab still rising at the end -> not converged
    r2 = slab_stopping_power(t, t.copy(), L_z=25.0)
    assert r2["status"] == "not_converged", r2
    # conservation guard
    assert conservation_guard([234, 233.8])["ok"]
    assert not conservation_guard([162, 4.6])["ok"]
    print("stopping_power._selftest: all assertions passed")


if __name__ == "__main__":
    _selftest()
