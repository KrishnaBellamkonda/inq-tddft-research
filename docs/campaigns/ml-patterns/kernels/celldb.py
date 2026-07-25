"""Cell resolution from the run database (ml-patterns campaign).

Selects the form-factor cut (E=100, sigma sweep) and wake cut (sigma=5, E sweep)
WP runs and their MATCHED point-classical partners, matching WITHIN a cut on the
jellium background (r_s, box L_z, spacing dx) and velocity — NOT the DB
`best_twin_id` (which matched on energy only and points across densities).

Bath density (project memory reference_canonical_bath_density):
  bath-only field = density_system when it integrates to N_bath electrons
  (the `_wf` convention), else density_total - density_wp.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

REPO = "/local/data/public/skcb2/tddft"
DB_CSV = os.path.join(REPO, "docs/run_database.csv")

# Pinned splits (campaign <resolved_decisions> / ADR 0011)
FF_CALIB_SIGMA = [1.0, 5.0]
FF_HELDOUT_SIGMA = [0.5, 3.0, 8.0]
FF_ENERGY = 100.0
WAKE_SIGMA = 5.0


def load_db() -> pd.DataFrame:
    return pd.read_csv(DB_CSV)


def _results_dir(run_path: str, rel: str) -> str:
    """DB density dirs are relative to <run_path>/results/."""
    if not isinstance(rel, str) or rel == "" or rel == "nan":
        return ""
    return os.path.join(run_path, "results", rel)


def _abs_exists(p: str) -> bool:
    return bool(p) and os.path.isdir(p)


def bath_series_dir(row) -> tuple[str, str]:
    """Return (kind, dir) for the bath-only series.

    kind in {'system_bathonly', 'total_minus_wp', 'total_wpincl', 'none'}.
    Prefers density_system if it is bath-only (we verify by electron count at
    load time); falls back to total-minus-wp if density_wp exists.
    """
    rp = row.run_path
    sys_dir = _results_dir(rp, row.density_system_vti_dir)
    tot_dir = _results_dir(rp, row.density_total_vti_dir)
    wp_dir = _results_dir(rp, row.density_wp_vti_dir)
    if _abs_exists(sys_dir):
        return "system", sys_dir
    if _abs_exists(tot_dir):
        return "total", tot_dir
    return "none", ""


def gs_path(row) -> str:
    rp = row.run_path
    # density_gs_system.vti lives beside the series dirs
    cand = os.path.join(rp, "results/raw/vti/density_gs_system/density_gs_system.vti")
    if os.path.isfile(cand):
        return cand
    # fallback: search
    base = os.path.join(rp, "results/raw/vti/density_gs_system")
    if os.path.isdir(base):
        for f in os.listdir(base):
            if f.endswith(".vti"):
                return os.path.join(base, f)
    return ""


def _has_density(row) -> bool:
    k, d = bath_series_dir(row)
    return _abs_exists(d) and bool(gs_path(row))


def _pick_classical(df, rs, Lz, dx, velocity, vtol=0.02):
    """Best point-classical (coulombic) partner at matched background+velocity."""
    cl = df[(df.system == "jellium") & (df.wp_enabled == False) &
            (df.classical_potential_form == "coulombic")].copy()
    cl = cl[np.isclose(cl.r_s, rs, rtol=0.02) &
            np.isclose(cl.cell_z, Lz, atol=1.0) &
            np.isclose(cl.spacing_bohr, dx, atol=0.01) &
            (np.abs(cl.velocity_au - velocity) <= vtol * max(velocity, 1e-9))]
    cl = cl[cl.apply(_has_density, axis=1)]
    if cl.empty:
        return None
    # prefer most frames
    cl = cl.assign(_nf=cl.density_total_vti_nframes.fillna(0))
    return cl.sort_values("_nf", ascending=False).iloc[0]


def _pick_wp(df, sigma, energy, rs=5.68986719379413, Lz=50.0, dx=0.4,
             prefer_wf=True):
    """WP run at given sigma & energy in the L50 rs=5.69 background, with density."""
    wp = df[(df.system == "jellium") & (df.wp_enabled == True) &
            np.isclose(df.sigma_wp_bohr.fillna(-1), sigma, atol=0.05) &
            np.isclose(df.energy_ev.fillna(-1), energy, rtol=0.02) &
            np.isclose(df.r_s.fillna(-1), rs, rtol=0.02) &
            np.isclose(df.cell_z.fillna(-1), Lz, atol=1.0) &
            np.isclose(df.spacing_bohr.fillna(-1), dx, atol=0.01)].copy()
    wp = wp[wp.apply(_has_density, axis=1)]
    if wp.empty:
        return None
    # prefer runs whose run_name ends in _wf (bath-only density_system) and most frames
    wp = wp.assign(
        _wf=wp.run_name.str.endswith("_wf").astype(int),
        _wpnf=wp.density_wp_vti_nframes.fillna(0),
        _nf=wp.density_system_vti_nframes.fillna(0),
    )
    sort_cols = (["_wf", "_nf"] if prefer_wf else ["_nf"])
    return wp.sort_values(sort_cols, ascending=False).iloc[0]


def resolve_form_factor_cells():
    """Return dict role -> list of resolved cells for the E=100 sigma sweep."""
    df = load_db()
    out = {"calibration": [], "heldout": [], "skipped": []}
    plan = [("calibration", s) for s in FF_CALIB_SIGMA] + \
           [("heldout", s) for s in FF_HELDOUT_SIGMA]
    for role, sigma in plan:
        wp = _pick_wp(df, sigma, FF_ENERGY)
        if wp is None:
            out["skipped"].append({"sigma": sigma, "role": role,
                                   "reason": "no WP run with bath density"})
            continue
        cl = _pick_classical(df, wp.r_s, wp.cell_z, wp.spacing_bohr, wp.velocity_au)
        if cl is None:
            out["skipped"].append({"sigma": sigma, "role": role,
                                   "reason": "no matched coulombic classical"})
            continue
        out[role].append(_cell(wp, cl, sigma=sigma))
    return out


def resolve_wake_cells():
    """Return dict role -> list of resolved cells for the sigma=5 energy sweep.

    Split: sort matched energies ascending by velocity; calibration = even index,
    held-out = odd (deterministic).
    """
    df = load_db()
    wp = df[(df.system == "jellium") & (df.wp_enabled == True) &
            np.isclose(df.sigma_wp_bohr.fillna(-1), WAKE_SIGMA, atol=0.05) &
            np.isclose(df.r_s.fillna(-1), 5.68986719379413, rtol=0.02) &
            np.isclose(df.cell_z.fillna(-1), 50.0, atol=1.0)].copy()
    wp = wp[wp.apply(_has_density, axis=1)]
    # one WP per energy (most frames)
    wp = wp.assign(_nf=wp.density_system_vti_nframes.fillna(0))
    best = {}
    for _, r in wp.iterrows():
        e = round(float(r.energy_ev), 1)
        if e not in best or r._nf > best[e]._nf:
            best[e] = r
    energies = sorted(best.keys(), key=lambda e: best[e].velocity_au)
    out = {"calibration": [], "heldout": [], "skipped": []}
    for i, e in enumerate(energies):
        wpr = best[e]
        cl = _pick_classical(df, wpr.r_s, wpr.cell_z, wpr.spacing_bohr,
                             wpr.velocity_au)
        role = "calibration" if i % 2 == 0 else "heldout"
        if cl is None:
            out["skipped"].append({"energy": e, "v": float(wpr.velocity_au),
                                   "role": role, "reason": "no matched classical"})
            continue
        out[role].append(_cell(wpr, cl, energy=e))
    return out


def _cell(wp, cl, sigma=None, energy=None):
    kbw, dwp = bath_series_dir(wp)
    kbc, dcl = bath_series_dir(cl)
    return {
        "sigma_wp": float(wp.sigma_wp_bohr) if sigma is None else float(sigma),
        "sigma_pot": float(wp.sigma_pot_bohr),
        "energy_ev": float(wp.energy_ev) if energy is None else float(energy),
        "velocity_au": float(wp.velocity_au),
        "r_s": float(wp.r_s),
        "omega_p_ev": float(wp.omega_p_ev),
        "n0": float(wp.n0),
        "kF": float(wp.kF),
        "dx": float(wp.spacing_bohr),
        "frame_dt_au_wp": float(wp.frame_dt_au) if not pd.isna(wp.frame_dt_au) else None,
        "frame_dt_au_cl": float(cl.frame_dt_au) if not pd.isna(cl.frame_dt_au) else None,
        "wp_run": wp.run_id,
        "cl_run": cl.run_id,
        "wp_bath_dir": dwp, "wp_bath_kind": kbw,
        "cl_bath_dir": dcl, "cl_bath_kind": kbc,
        "wp_gs": gs_path(wp), "cl_gs": gs_path(cl),
        "wp_wp_dir": _results_dir(wp.run_path, wp.density_wp_vti_dir),
        "wp_total_dir": _results_dir(wp.run_path, wp.density_total_vti_dir),
    }
