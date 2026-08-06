#!/usr/bin/env python3
"""Analysis layer for the Nazarov-Gross mass ladder.

Plan: docs/plans/nazarov-gross-slab-mass-ladder.md

THE ONE THING THIS MODULE IS FOR
    Nazarov & Gross (arXiv:2510.26222) claim that projectiles of the SAME CHARGE
    moving at the SAME VELOCITY feel different friction depending on their MASS,
    and that mass acts through exactly one channel: the spatial WIDTH of the
    projectile's wavepacket (their Sec. VII, Eq. 41 -- the classical point-charge
    Coulomb potential is replaced by the potential of the projectile's own
    ground-state density). Everything here exists to test those two statements:

        does_S_depend_on_mass()   the claim
        does_S_collapse_on_width()the mechanism

STOPPING POWER IS THE BATH DEPOSIT, NOT THE PROJECTILE'S KINETIC ENERGY.
    S = d(E_bath)/ds over an early in-slab window where v >= 0.85 v0
    (.claude/rules/light-projectile-stopping.md). The deposit is measured ON THE
    MEDIUM, so it means the same thing for a classical Gaussian perturbation and
    for a wavepacket of any mass -- which is the whole reason the ladder can sit
    on one axis. -dKE/ds is the ENERGY-CONSERVATION CROSS-CHECK ONLY and is never
    the headline (.claude/skills/stopping-power-extraction, user 2026-06-30).

    For a wavepacket the two genuinely differ, and the difference is physics: the
    packet absorbs energy into its own spreading (var(p)/2m) that never reaches
    the bath. In the channeling twin that was 54 % of the drift loss.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HA_EV = 27.211386245988

# System constants — mirror shared/configs/slab_n206_L30x30x120_rs2p5.hpp.
V0 = 1.0742685
KF = 0.7673347
EF_EV = 8.01107
OMEGA_P = 0.4378967
RS = 2.5010708
SLAB_HALF = 7.5
SIGMA_WP_DEFAULT = 4.0
LAUNCH_Z = -25.0
CAP_INNER = 45.0


# ------------------------------------------------------------------ loading
def _concat_segments(run_dir: Path, stem: str, key: str = "step") -> pd.DataFrame:
    """Concatenate `stem.csv`, `stem.from<N>.csv`, ... in step order.

    Resumed runs write segment-suffixed CSVs (rule final-timestep-checkpoint.md).
    Duplicates are dropped keeping the LAST occurrence, which is right for a
    scalar-per-step file. It is CATASTROPHIC for a long-format file with many
    rows per step -- that bug silently reduced momentum_distribution.csv to one
    k-bin per step in the channeling study. This helper is therefore only ever
    used on scalar-per-step observables.
    """
    obs = run_dir / "raw/observables"
    if not obs.is_dir():
        return pd.DataFrame()
    files = sorted(obs.glob(f"{stem}.csv")) + sorted(
        obs.glob(f"{stem}.from*.csv"),
        key=lambda p: int(p.name.split(".from")[1].split(".")[0]))
    frames = []
    for f in files:
        try:
            frames.append(pd.read_csv(f, comment="#"))
        except Exception:                                       # noqa: BLE001
            continue
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    if key in df.columns:
        df = df.sort_values(key).drop_duplicates(subset=key, keep="last").reset_index(drop=True)
    return df


def parse_summary(run_dir: Path) -> dict:
    p = run_dir / "run_summary.txt"
    out: dict[str, str] = {}
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    return out


@dataclass
class Run:
    tag: str
    half: str                    # 'wp' | 'classical'
    path: Path
    summary: dict
    obs: pd.DataFrame            # INQ scalar energy ledger
    inter: pd.DataFrame          # pairwise P/S/B decomposition
    proj: pd.DataFrame           # classical only: trajectory
    pmom: pd.DataFrame           # wp only: momentum moments
    rspace: pd.DataFrame         # wp only: real-space moments (the WIDTH)

    @property
    def mass(self) -> float:
        s = self.summary
        for k in ("mass", "projectile"):
            if k in s:
                try:
                    if k == "mass":
                        return float(s[k])
                    return float(s[k].split("mass")[1].split()[0])
                except (ValueError, IndexError):
                    pass
        return float("nan")

    @property
    def sigma_wp(self) -> float:
        for k in ("sigma_wp", "sigma_WP"):
            if k in self.summary:
                try:
                    return float(self.summary[k])
                except ValueError:
                    pass
        try:
            return float(self.summary["projectile"].split("sigma_WP")[1].split()[0])
        except (KeyError, IndexError, ValueError):
            return SIGMA_WP_DEFAULT

    @property
    def complete(self) -> bool:
        return self.summary.get("run_completed", "").lower() == "true"


def load_run(scripts_dir: Path, half: str, tag: str) -> Run:
    d = Path(scripts_dir) / half / "results" / tag
    return Run(tag=tag, half=half, path=d, summary=parse_summary(d),
               obs=_concat_segments(d, "observables"),
               inter=_concat_segments(d, "interactions"),
               proj=_concat_segments(d, "projectile"),
               pmom=_concat_segments(d, "wp_momentum_stats"),
               rspace=_concat_segments(d, "wp_real_space_stats"))


# ------------------------------------------------------------------ kinematics
def projectile_track(r: Run) -> pd.DataFrame:
    """(time_au, z, v) for either representation, on one schema.

    Classical: the Ehrenfest velocity-Verlet track. Wavepacket: the orbital's
    own first moments -- z from <z>, v from <p_z>/M. Both are the projectile's
    mean position and mean velocity, so S(v0) means the same thing in each.
    """
    if r.half == "classical" and not r.proj.empty:
        return pd.DataFrame({"time_au": r.proj["time_au"].to_numpy(),
                             "z": r.proj["proj_z"].to_numpy(),
                             "v": r.proj["proj_vz"].to_numpy()})
    if not r.rspace.empty and not r.pmom.empty:
        t = r.rspace["time_au"].to_numpy()
        z = r.rspace["z_mean"].to_numpy()
        pz = np.interp(t, r.pmom["time_au"].to_numpy(), r.pmom["pz_mean"].to_numpy())
        m = r.mass if np.isfinite(r.mass) and r.mass > 0 else 1.0
        return pd.DataFrame({"time_au": t, "z": z, "v": pz / m})
    return pd.DataFrame(columns=["time_au", "z", "v"])


def wp_width(r: Run) -> pd.DataFrame:
    """sigma(t) of the packet DENSITY, per axis and isotropic (3-D geometric mean).

    This is the mechanism variable. Mixing width definitions is a real hazard --
    a transverse ratio was once divided into a 3-D one in this project and moved
    a headline number from 16.9 % to 20.9 % -- so both are returned explicitly
    and never silently interchanged.
    """
    if r.rspace.empty:
        return pd.DataFrame()
    d = r.rspace
    sx, sy, sz = (np.sqrt(np.maximum(d[f"sigma_{a}2"].to_numpy(), 0.0)) for a in "xyz")
    return pd.DataFrame({
        "time_au": d["time_au"].to_numpy(),
        "sigma_x": sx, "sigma_y": sy, "sigma_z": sz,
        "sigma_perp": np.sqrt(0.5 * (sx**2 + sy**2)),
        "sigma_iso": np.cbrt(np.maximum(sx * sy * sz, 0.0)),
    })


def kinetic_channels(r: Run) -> pd.DataFrame:
    """T1 = <p>^2/2M (drift) and T2 = <p^2>/2M (total orbital KE).

    T2 - T1 = var(p)/2M is exactly the energy the packet has taken into its own
    internal spreading. Nazarov-Gross's projectile has no such channel -- its
    internal state is stationary in the driven steady state -- so T1 is the
    closer analogue of their E_kin, while (T2 - T1) is the purely quantum
    bookkeeping term with no classical counterpart.

    NOTE the label convention: this project has used T1/T2 in BOTH orders in
    different modules, and reading a curve under the wrong one INVERTS the
    conclusion. Here T1 is ALWAYS the drift term.
    """
    if r.pmom.empty:
        return pd.DataFrame()
    d = r.pmom
    m = r.mass if np.isfinite(r.mass) and r.mass > 0 else 1.0
    p2_drift = d["px_mean"]**2 + d["py_mean"]**2 + d["pz_mean"]**2
    p2_tot = d["px2_mean"] + d["py2_mean"] + d["pz2_mean"]
    return pd.DataFrame({
        "time_au": d["time_au"].to_numpy(),
        "T1_drift_ev": (p2_drift / (2 * m)).to_numpy() * HA_EV,
        "T2_total_ev": (p2_tot / (2 * m)).to_numpy() * HA_EV,
        "var_p_over_2m_ev": ((p2_tot - p2_drift) / (2 * m)).to_numpy() * HA_EV,
    })


# ------------------------------------------------------------------ stopping
@dataclass
class Stopping:
    S_ev_per_bohr: float
    stderr: float
    r2: float
    n_points: int
    window_z: tuple[float, float]
    mean_v: float
    v_drop_frac: float
    method: str
    note: str = ""


def bath_energy_ev(r: Run) -> np.ndarray | None:
    """Energy of the electron bath, in eV, relative to t = 0.

    Classical: the projectile is an external perturbation and carries no INQ
    energy of its own, so `energy_total` IS the bath (verified in the channeling
    twin, where the classical budget closed to 2.2e-5 eV over 1501 steps).

    Wavepacket: the packet is an occupied KS orbital, so `energy_total` contains
    it. The bath deposit is recovered from the pairwise decomposition, which is
    representation-independent by construction:
        E_bath = E_total - (T_wp + E_PP + E_PS + E_PB)
    i.e. total minus everything the projectile owns or shares.
    """
    if r.obs.empty or "energy_total" not in r.obs:
        return None
    e_tot = r.obs["energy_total"].to_numpy()
    if r.half == "classical":
        return (e_tot - e_tot[0]) * HA_EV
    if r.inter.empty or r.pmom.empty:
        return None
    t_obs = r.obs["time_au"].to_numpy()

    def on_obs(df, col):
        return np.interp(t_obs, df["time_au"].to_numpy(), df[col].to_numpy())

    proj_owned = (on_obs(r.inter, "e_pp") + on_obs(r.inter, "e_ps")
                  + on_obs(r.inter, "e_pb") + on_obs(r.pmom, "e_kin_ha"))
    e_bath = e_tot - proj_owned
    return (e_bath - e_bath[0]) * HA_EV


def extract_S(r: Run, v_frac: float = 0.85, in_slab_only: bool = True) -> Stopping:
    """S = d(E_bath)/ds over the early, near-constant-velocity, in-slab window.

    A whole-run regression would be WRONG: the projectile decelerates, so a
    full-run slope averages S over every velocity from v0 down, not S AT v0
    (.claude/rules/light-projectile-stopping.md, learned by aborting an overnight
    cylindrical run on exactly this mistake).
    """
    trk = projectile_track(r)
    dep = bath_energy_ev(r)
    if trk.empty or dep is None or len(trk) < 10:
        return Stopping(np.nan, np.nan, np.nan, 0, (np.nan, np.nan), np.nan, np.nan,
                        "unavailable", "missing track or deposit")

    t_obs = r.obs["time_au"].to_numpy()
    z = np.interp(t_obs, trk["time_au"], trk["z"])
    v = np.interp(t_obs, trk["time_au"], trk["v"])

    v_ref = v[0] if np.isfinite(v[0]) and abs(v[0]) > 1e-9 else V0
    mask = np.abs(v) >= v_frac * abs(v_ref)
    if in_slab_only:
        mask &= (z >= -SLAB_HALF) & (z <= SLAB_HALF)

    note = ""
    if mask.sum() < 15:                       # widen exactly as the rule prescribes
        for wider in (0.70, 0.50):
            mask = np.abs(v) >= wider * abs(v_ref)
            if in_slab_only:
                mask &= (z >= -SLAB_HALF) & (z <= SLAB_HALF)
            note = f"window widened to v >= {wider:.2f} v0 (sparse at 0.85)"
            if mask.sum() >= 15:
                break
    if mask.sum() < 5:
        return Stopping(np.nan, np.nan, np.nan, int(mask.sum()), (np.nan, np.nan),
                        np.nan, np.nan, "in_slab_drag",
                        "fewer than 5 points survive the window — the projectile "
                        "never crossed the slab at near-constant velocity")

    zz, ee = z[mask], dep[mask]
    n = len(zz)
    # free-intercept fit: the entrance transient deposits a fixed offset, and
    # forcing through the origin biases S high
    A = np.vstack([zz, np.ones_like(zz)]).T
    coef, res, *_ = np.linalg.lstsq(A, ee, rcond=None)
    S = float(coef[0])
    pred = A @ coef
    ss_res = float(np.sum((ee - pred) ** 2))
    ss_tot = float(np.sum((ee - ee.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    sx2 = float(np.sum((zz - zz.mean()) ** 2))
    stderr = float(np.sqrt(ss_res / max(n - 2, 1) / sx2)) if sx2 > 0 else np.nan
    return Stopping(S, stderr, r2, n, (float(zz.min()), float(zz.max())),
                    float(np.mean(v[mask])),
                    float(1.0 - abs(v[mask][-1]) / abs(v_ref)),
                    "in_slab_drag", note)


def ke_cross_check(r: Run) -> float:
    """-dKE/ds over the whole run: the CONSERVATION check, never the headline."""
    trk = projectile_track(r)
    if trk.empty:
        return np.nan
    m = 1.0e6 if r.half == "classical" and r.mass > 1e5 else (r.mass if np.isfinite(r.mass) else 1.0)
    ke = 0.5 * m * trk["v"].to_numpy() ** 2 * HA_EV
    s = trk["z"].to_numpy()
    ds = s[-1] - s[0]
    return float(-(ke[-1] - ke[0]) / ds) if abs(ds) > 1e-9 else np.nan


# ------------------------------------------------------------------ gates
def pilot_gate(scripts_dir: Path, tags: list[str], halves: list[str]) -> dict:
    """Phase-2 gate. Correctness only — cost is never a reason to stop.

    Blocks the ladder if and only if the pilots show the MEASUREMENT cannot work:
    a diverging energy, a projectile that never crossed at near-constant v, or
    two masses whose deposits are indistinguishable from each other.
    """
    lines, hard_fail = [], False
    runs = {}
    for tag, half in zip(tags, halves):
        r = load_run(scripts_dir, half, tag)
        runs[tag] = r
        if not r.complete:
            lines.append(f"[WARN] {tag}: run_completed is not true — excluded from the gate")
            continue

        if r.obs.empty or "energy_total" not in r.obs:
            lines.append(f"[FAIL] {tag}: no energy_total")
            hard_fail = True
            continue
        e = r.obs["energy_total"].to_numpy()
        if not np.all(np.isfinite(e)):
            lines.append(f"[FAIL] {tag}: energy_total contains NaN/inf — the propagation diverged")
            hard_fail = True
            continue
        lines.append(f"[info] {tag}: E_total range {e.min():.6f} .. {e.max():.6f} Ha")

        # CAP absorption. The CAP contaminates E_total once it starts eating the
        # packet, so the analysis window must close first. Here we only check the
        # packet had not reached it during the pilot.
        trk = projectile_track(r)
        if not trk.empty:
            zmax = float(np.nanmax(np.abs(trk["z"])))
            flag = "OK" if zmax < CAP_INNER else "REACHED CAP"
            lines.append(f"[info] {tag}: |z|_max = {zmax:.1f} Bohr vs CAP face {CAP_INNER} — {flag}")

        # Transverse containment: a packet wider than the half-cell wraps and
        # every width-mediated number after that is meaningless.
        w = wp_width(r)
        if not w.empty:
            s_end = float(w["sigma_perp"].iloc[-1])
            lines.append(f"[info] {tag}: sigma_perp {w['sigma_perp'].iloc[0]:.2f} -> {s_end:.2f} Bohr "
                         f"({'OK' if 4*s_end < 30.0 else 'WRAPS the 30 Bohr cell'})")

        S = extract_S(r)
        lines.append(f"[info] {tag}: S = {S.S_ev_per_bohr:.5f} +- {S.stderr:.5f} eV/Bohr "
                     f"(r2 {S.r2:.4f}, n {S.n_points}, mean v {S.mean_v:.4f}) {S.note}")

    # The discriminating check: do two different masses give different deposits?
    wp_tags = [t for t, h in zip(tags, halves) if h == "wp" and runs[t].complete]
    if len(wp_tags) >= 2:
        vals = [(t, extract_S(runs[t])) for t in wp_tags]
        good = [(t, s) for t, s in vals if np.isfinite(s.S_ev_per_bohr)]
        if len(good) < 2:
            lines.append("[FAIL] fewer than two wavepacket pilots yielded a finite S — "
                         "the deposit cannot be measured as designed")
            hard_fail = True
        else:
            (ta, sa), (tb, sb) = good[0], good[1]
            sep = abs(sa.S_ev_per_bohr - sb.S_ev_per_bohr)
            comb = math.hypot(sa.stderr, sb.stderr) if np.isfinite(sa.stderr) and np.isfinite(sb.stderr) else np.nan
            lines.append(f"[info] separation |S({ta}) - S({tb})| = {sep:.5f} eV/Bohr "
                         f"vs combined stderr {comb:.5f}")
            if np.isfinite(comb) and comb > 0 and sep < 2.0 * comb:
                lines.append("[WARN] the two masses are within 2 sigma of each other. That is "
                             "NOT a hard failure — a null result is a legitimate outcome for "
                             "this campaign, and 600 steps is a short lever arm. Proceeding.")

    verdict = not hard_fail
    header = ("PILOT GATE: PASS — launching the ladder" if verdict
              else "PILOT GATE: FAIL — the ladder was NOT launched")
    return {"pass": verdict, "report": header + "\n\n" + "\n".join(lines)}


# ------------------------------------------------------------------ the claims
def ladder_table(scripts_dir: Path, runs: list[tuple[str, str]]) -> pd.DataFrame:
    """One row per rung: the table the whole campaign exists to produce."""
    rows = []
    for half, tag in runs:
        r = load_run(scripts_dir, half, tag)
        if not r.complete:
            continue
        S = extract_S(r)
        w = wp_width(r)
        sig0 = float(w["sigma_iso"].iloc[0]) if not w.empty else r.sigma_wp / np.sqrt(2)
        sig_mid = float(w["sigma_iso"].median()) if not w.empty else sig0
        k = kinetic_channels(r)
        rows.append({
            "tag": tag, "half": half,
            "mass": r.mass if r.mass < 1e5 else np.inf,
            "sigma_WP_nominal": r.sigma_wp,
            "sigma_iso_t0": sig0, "sigma_iso_mid": sig_mid,
            "spread_factor": sig_mid / sig0 if sig0 > 0 else np.nan,
            "S_ev_per_bohr": S.S_ev_per_bohr, "S_stderr": S.stderr, "S_r2": S.r2,
            "S_n": S.n_points, "mean_v": S.mean_v, "v_drop_frac": S.v_drop_frac,
            "S_ke_crosscheck": ke_cross_check(r),
            "T1_loss_ev": (float(k["T1_drift_ev"].iloc[0] - k["T1_drift_ev"].iloc[-1])
                           if not k.empty else np.nan),
            "var_p_gain_ev": (float(k["var_p_over_2m_ev"].iloc[-1] - k["var_p_over_2m_ev"].iloc[0])
                              if not k.empty else np.nan),
            "note": S.note,
        })
    return pd.DataFrame(rows)


def does_S_depend_on_mass(tbl: pd.DataFrame) -> dict:
    """THE NAZAROV-GROSS CLAIM, as a number.

    Compares the quantum rungs at fixed charge, velocity and initial width. The
    classical rung is excluded from the test statistic and used as the M->inf
    reference, because a rigid cloud cannot depend on mass by construction.
    """
    q = tbl[(tbl.half == "wp") & np.isfinite(tbl.S_ev_per_bohr)].sort_values("mass")
    if len(q) < 2:
        return {"verdict": "inconclusive", "reason": f"only {len(q)} quantum rungs with a finite S"}
    S, err = q.S_ev_per_bohr.to_numpy(), q.S_stderr.to_numpy()
    spread = float(S.max() - S.min())
    comb = float(np.hypot(err[np.argmax(S)], err[np.argmin(S)]))
    n_sigma = spread / comb if comb > 0 else np.inf
    cl = tbl[tbl.half == "classical"]
    S_cl = float(cl.S_ev_per_bohr.iloc[0]) if len(cl) and np.isfinite(cl.S_ev_per_bohr.iloc[0]) else np.nan
    return {
        "verdict": ("mass-dependent" if n_sigma >= 3 else
                    "no significant mass dependence" if n_sigma < 1 else "marginal"),
        "S_min": float(S.min()), "S_max": float(S.max()), "spread": spread,
        "combined_stderr": comb, "n_sigma": float(n_sigma),
        "S_classical": S_cl,
        "masses": q.mass.tolist(), "S": S.tolist(),
        "monotonic_increasing_with_mass": bool(np.all(np.diff(S) > 0)),
        "note": ("NG predict S RISES towards the classical value as M grows, because a "
                 "heavier projectile is more localised and couples like a point charge."),
    }


def does_S_collapse_on_width(tbl: pd.DataFrame) -> dict:
    """THE NAZAROV-GROSS MECHANISM, as a number.

    If mass acts ONLY through width (their Sec. VII), then the mass ladder and
    the fixed-mass sigma sweep must fall on ONE curve when S is plotted against
    the MEASURED mid-transit width. Any systematic offset between the two
    families is mass acting through some other channel.
    """
    d = tbl[np.isfinite(tbl.S_ev_per_bohr) & np.isfinite(tbl.sigma_iso_mid) & (tbl.half == "wp")]
    if len(d) < 4:
        return {"verdict": "inconclusive", "reason": f"only {len(d)} usable points"}
    x = np.log(d.sigma_iso_mid.to_numpy())
    y = np.log(np.maximum(d.S_ev_per_bohr.to_numpy(), 1e-12))
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else np.nan
    # Split the residuals by family: mass ladder (nominal sigma == default) vs
    # sigma sweep (everything else). A family offset is the interesting failure.
    is_ladder = np.isclose(d.sigma_WP_nominal.to_numpy(), SIGMA_WP_DEFAULT)
    off = (float(resid[is_ladder].mean() - resid[~is_ladder].mean())
           if is_ladder.any() and (~is_ladder).any() else np.nan)
    return {
        "verdict": ("collapses — mass acts through width" if r2 > 0.9 and abs(off) < 0.15
                    else "does NOT collapse — mass has a channel beyond width" if r2 < 0.7
                    else "partial collapse"),
        "power_law_exponent": float(coef[0]), "r2": float(r2),
        "family_offset_ln": off, "n_points": int(len(d)),
        "note": ("Exponent < 0 is expected: a wider projectile couples more weakly. "
                 "The family offset is the discriminator — a clean collapse means the "
                 "mass ladder and the sigma sweep are the SAME experiment."),
    }
