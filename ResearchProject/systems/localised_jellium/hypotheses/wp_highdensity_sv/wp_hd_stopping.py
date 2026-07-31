"""
wp_hd_stopping — analysis adapter for the WP twin of the high-density classical
S(v) benchmark (campaign classical-highdensity-sv).

Thin layer over ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/
ks_stopping.py, which already implements exactly the KE/position definitions this
study uses (docs/plans/bulk-jellium-ks-stopping.md section 4):

    T1 = <p^2>/2m   T2 = <p>^2/2m   s3 = circular centroid   s4 = integral <p_z> dt

Two things it adds:

1. LAYOUT. These runs write to  <scripts>/wp_highdensity_sv/wp/results/<name>/raw/
   observables/, whereas ks_stopping.load_wp_run expects <run_dir>/results/raw/
   observables/. load_run() below points the loader at the right directory.

2. CAP BASELINE — the scientifically important part. The cap_check replica showed
   that in an EMPTY box with no bath and no forces, the CAP alone drags the
   surviving packet's <p_z> from 2.00 to 0.61 over 48 a.u., because sigma_WP = 0.5
   spreads fast enough that the packet's LEADING edge is preferentially absorbed.
   A stopping power fitted where the CAP is active would therefore measure the
   CAP. cap_corrected() differences each slab run against its vacuum twin (same
   grid, same k0, same step count, no bath) so what remains is the bath's doing.

Reference classical values (the benchmark being compared against) come from
hypotheses/classical_highdensity_sv/sv_sweep/S_summary.csv.

Plan: docs/plans/wavepacket-highdensity-sv-twin.md
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
_KS = REPO / "ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping"
if str(_KS) not in sys.path:
    sys.path.insert(0, str(_KS))

import ks_stopping as K  # noqa: E402

HA_TO_EV = K.HA_TO_EV

# Geometry, mirrored from shared/configs/slab_n100_L35x35x85.hpp. run_summary.txt
# is the authority if these ever disagree.
LX = LY = 35.0
LZ = 85.0
DX_PRODUCTION = 0.40
SLAB_HALF = 12.5
N_ELECTRONS = 100
R_S = 4.183
SIGMA_WP = 0.5
LAUNCH_Z = -24.0
DT = 0.04
CAP_L = 12.5
CAP_INNER = LZ / 2 - CAP_L          # +30.0 Bohr
CAP_ETA = -1.0

# Fit window: set by two independent limits that agree (plan section 6b).
#   transverse periodic images overlap at 6 sigma_d = L_xy  ->  t = 4.12 a.u.
#   CAP attrition stays below ~0.3 % of norm for t <~ 4 a.u.
# Both are sigma-dependent, so set_campaign() recomputes them; the values here are
# the sigma = 0.5 defaults.
FIT_T0 = 0.5
FIT_T1 = 4.12
T_TRANSVERSE = 4.12

# The four production points (v = 4.0 and 4.5 dropped: momentum aliasing at
# dx = 0.40 reaches +17.9 % and +55.1 % in sigma_pz^2 — plan section 6b item 6).
# That exclusion is a sigma = 0.5 property; see aliasing_bias_pct().
VELOCITIES = (2.0, 2.5, 3.0, 3.5)
N_STEPS = {2.0: 3623, 2.5: 2898, 3.0: 2415, 3.5: 2070}

SCRIPTS = REPO / "ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv"
WP_RESULTS = SCRIPTS / "wp" / "results"
VAC_RESULTS = SCRIPTS / "cap_check" / "results"
CLASSICAL_SUMMARY = (REPO / "ResearchProject/systems/localised_jellium/hypotheses"
                     / "classical_highdensity_sv/sv_sweep/S_summary.csv")

# --------------------------------------------------------------------------
# Sigma campaigns (user instruction 2026-07-31)
# --------------------------------------------------------------------------
# Three campaigns share this module: sigma_WP = 0.5 (the original twin of the
# classical benchmark), 2.0 and 3.0. Only sigma differs — same GS, dx, launch z,
# CAPs, dt, step counts, velocities — so one code path serves all three and the
# traces are comparable point for point.
#
# Run directories are PREFIXED s<sigma>_ for every campaign except 0.5, whose
# names stay bare so the completed runs, their notebooks and wp_S_summary.csv
# keep resolving unchanged.
SIGMAS = (0.5, 2.0, 3.0)
_SIGMA_TAG = ""


def sigma_tag(sigma: float) -> str:
    """'' for the legacy 0.5 campaign, 's2p0_' / 's3p0_' otherwise."""
    return "" if abs(sigma - 0.5) < 1e-9 else "s" + f"{sigma:.1f}".replace(".", "p") + "_"


def transverse_overlap_time(sigma: float, l_xy: float = LX) -> float:
    """
    When the packet's own periodic images start to overlap: 6 sigma_d(t) = L_xy.

    Solving sigma^2/2 + t^2/(2 sigma^2) = (L_xy/6)^2 for t gives
        t = sigma * sqrt(2 (L_xy/6)^2 - sigma^2).
    Reproduces the hand-derived 4.12 a.u. at sigma = 0.5, and shows why the larger
    campaigns are so much better behaved: 16.0 a.u. at sigma = 2, 23.1 at sigma = 3.
    """
    target = 2.0 * (l_xy / 6.0) ** 2
    return float(sigma * np.sqrt(max(target - sigma**2, 0.0)))


def set_campaign(sigma: float) -> None:
    """
    Point the module at one sigma campaign. Mutates module state, so notebooks
    call it once at the top and everything downstream (run paths, fit window,
    labels) follows. Idempotent.
    """
    global SIGMA_WP, _SIGMA_TAG, FIT_T1, T_TRANSVERSE
    SIGMA_WP = float(sigma)
    _SIGMA_TAG = sigma_tag(SIGMA_WP)
    T_TRANSVERSE = transverse_overlap_time(SIGMA_WP)
    FIT_T1 = round(T_TRANSVERSE, 2)


def current_sigma() -> float:
    return SIGMA_WP


def name_for(v: float, sigma: float | None = None) -> str:
    tag = _SIGMA_TAG if sigma is None else sigma_tag(sigma)
    return tag + "v" + f"{v:.1f}".replace(".", "p")


def vac_name_for(v: float, sigma: float | None = None) -> str:
    """Matching CAP-only vacuum control for a production point."""
    return "vac_" + name_for(v, sigma)


def campaign_label(sigma: float | None = None) -> str:
    s = SIGMA_WP if sigma is None else sigma
    return f"sigma_WP = {s:g} Bohr"


def has_campaign(sigma: float) -> bool:
    """True when at least one velocity of this campaign has produced observables."""
    return any((WP_RESULTS / name_for(v, sigma) / "raw" / "observables"
                / "wp_momentum_stats.csv").exists() for v in VELOCITIES)


def sigma_d(t: np.ndarray | float, sigma: float | None = None) -> np.ndarray | float:
    """Free-Gaussian density width. Verified on this engine by the vacuum sweep."""
    s = SIGMA_WP if sigma is None else sigma
    return np.sqrt(s**2 / 2 + np.asarray(t) ** 2 / (2 * s**2))


def aliasing_bias_pct(v: float, sigma: float | None = None,
                      dx: float = DX_PRODUCTION) -> dict[str, float]:
    """
    Momentum-aliasing bias in <p_z> and sigma_pz^2 from folding at k_Nyq = pi/dx.

    The packet's k-distribution is Gaussian at k0 = v with std sigma_p =
    1/(sqrt2 sigma). The FFT grid represents it modulo G = 2 pi/dx, so weight
    beyond k_Nyq wraps back and corrupts the measured moments. This integrates the
    wrapped distribution directly rather than quoting a tail fraction, which badly
    understates the moment error (that mistake was made and corrected 2026-07-30).

    ACCURACY, stated honestly. This is a continuum fold of an ideal Gaussian; the
    engine measures a discrete grid in a finite box. Against the one MEASURED
    aliased point (sigma = 0.5, dx = 0.50, v = 4.5) it gives <p_z> = 3.20 and
    sigma_pz^2 = 10.27 where the run measured 3.44 and 9.05 — i.e. it captures a
    -24 %/+350 % catastrophe to within ~7 %/~13 % of itself. That is ample for its
    only job, sorting "negligible" from "fatal", but it is NOT a calibrated
    correction and must never be used to un-bias a measured moment.

    It reproduces the recorded dx = 0.40 sigma = 0.5 table exactly (0.05, 0.26,
    1.24, 5.06, 17.90, 55.05 %) because that table was generated by this model.

    Returns percent deviations. For sigma >= 2 at dx = 0.40 both are zero to
    machine precision at every velocity up to 4.5 — sigma_p = 1/(sqrt2 sigma) is
    0.354 and 0.236 against k_Nyq = 7.85, so nothing reaches the fold. That is why
    the sigma = 2 and 3 campaigns could also carry v = 4.0/4.5 (not requested; the
    grid is held at the same four points for comparability).
    """
    s = SIGMA_WP if sigma is None else sigma
    sigma_p = 1.0 / (np.sqrt(2.0) * s)
    k_nyq = np.pi / dx
    grid = 2.0 * k_nyq
    k = np.linspace(v - 12 * sigma_p, v + 12 * sigma_p, 200001)
    w = np.exp(-0.5 * ((k - v) / sigma_p) ** 2)
    w /= w.sum()
    k_wrapped = (k + k_nyq) % grid - k_nyq
    p_mean = float((w * k_wrapped).sum())
    p_var = float((w * (k_wrapped - p_mean) ** 2).sum())
    return {
        "p_mean": p_mean,
        "p_mean_err_pct": 100.0 * (p_mean - v) / v,
        "sigma_pz2": p_var,
        "sigma_pz2_err_pct": 100.0 * (p_var - sigma_p**2) / sigma_p**2,
        "k_nyq": float(k_nyq),
        "sigma_p": float(sigma_p),
    }


def _load_at(run_dir: Path, z0: float = LAUNCH_Z) -> K.WPRun:
    """
    Build a ks_stopping.WPRun from <run_dir>/raw/observables.

    ks_stopping.load_wp_run() hard-codes the layout <run_dir>/results/raw/
    observables, but these runs write <scripts>/wp/results/<name>/raw/
    observables — the run NAME sits where "results" would be. So the loader body
    is reproduced here against an explicit directory rather than shimmed; an
    earlier duck-typed shim failed because load_wp_run calls Path(run_dir)
    internally.

    NOTE on the two norm columns: the merge suffixes them _p (momentum-space
    Parseval constant, ~5e7, FFT-prefactor dependent) and _r (real-space
    integral |psi|^2 dV, ~1 at t=0). This uses the REAL-SPACE one — it is the
    physical norm, it is what the CAP correction needs, and it is the one that
    means anything when plotted.
    """
    obs = Path(run_dir) / "raw" / "observables"
    if not obs.is_dir():
        raise FileNotFoundError(f"no observables under {run_dir}")
    mom = K._concat_segments(obs, "wp_momentum_stats")
    pos = K._concat_segments(obs, "wp_real_space_stats")
    df = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))

    t = df["time_au"].to_numpy()
    T1 = df["e_kin_ha"].to_numpy()
    px, py, pz = (df[c].to_numpy() for c in ("px_mean", "py_mean", "pz_mean"))
    T2 = 0.5 * (px**2 + py**2 + pz**2)

    if "z_mean_circ" not in df.columns:
        raise KeyError("z_mean_circ missing — the naive centroid is unusable near "
                       "a cell face and this box wraps.")
    zc = df["z_mean_circ"].to_numpy()
    s3 = K.unwrap_periodic(zc, LZ)
    s3 = s3 - s3[0] + zc[0]        # unwrap fixes increments, not the offset
    s4 = z0 + np.concatenate([[0.0],
                              np.cumsum(0.5 * (pz[1:] + pz[:-1]) * np.diff(t))])

    norm = df["norm_check_r"].to_numpy() if "norm_check_r" in df \
        else df["norm_check"].to_numpy()

    return K.WPRun(run_dir=Path(run_dir), box_length_z=LZ, t=t,
                   step=df["step"].to_numpy(), T1=T1, T2=T2, pz=pz,
                   s3=s3, s3_naive=df["z_mean"].to_numpy(), s4=s4,
                   norm=norm, sigma_z=df["sigma_z_circ"].to_numpy())


def load_run(v: float, results_root: Path = WP_RESULTS) -> K.WPRun:
    """Load one slab (production) velocity point."""
    return _load_at(results_root / name_for(v))


def load_vacuum(v: float) -> K.WPRun | None:
    """Load the CAP-only vacuum control at the same velocity (None if absent)."""
    d = VAC_RESULTS / vac_name_for(v)
    try:
        return _load_at(d)
    except (FileNotFoundError, KeyError):
        return None


@dataclass
class CapCorrected:
    """Slab minus vacuum, on the slab run's time grid."""
    t: np.ndarray
    dT1: np.ndarray          # Ha, (T1_slab - T1_vac)
    dT2: np.ndarray          # Ha
    dpz: np.ndarray          # Bohr^-1
    norm_slab: np.ndarray
    norm_vac: np.ndarray

    @property
    def cap_only_pz_drop(self) -> np.ndarray:
        """How much of the apparent deceleration the CAP alone accounts for."""
        return self.norm_vac  # placeholder kept explicit; see cap_corrected()


def cap_corrected(run: K.WPRun, vac: K.WPRun) -> CapCorrected:
    """
    Difference the slab run against its vacuum twin.

    Both were launched identically and differ only by the presence of the bath,
    so at each time the vacuum run carries the pure CAP-attrition drift and the
    difference is the bath's contribution. Interpolated onto the slab time grid
    in case a resume left the two on slightly different step sets.
    """
    def on(x):
        return np.interp(run.t, vac.t, x)

    return CapCorrected(
        t=run.t,
        dT1=run.T1 - on(vac.T1),
        dT2=run.T2 - on(vac.T2),
        dpz=run.pz - on(vac.pz),
        norm_slab=run.norm,
        norm_vac=on(vac.norm),
    )


def slab_window(v: float, z0: float = LAUNCH_Z) -> tuple[float, float]:
    """
    Time interval over which the WP CENTROID is inside the slab, [-12.5, +12.5].

    This is the interval in which a stopping power can physically be measured:
    the projectile has to be IN the medium. With the campaign-matched launch at
    z0 = -24 there is 11.5 Bohr of vacuum standoff first, so

        v = 2.0 -> t in [ 5.75, 18.25]      v = 3.0 -> t in [3.83, 12.17]
        v = 2.5 -> t in [ 4.60, 14.60]      v = 3.5 -> t in [3.29, 10.43]

    Note this starts AFTER the localised window [FIT_T0, FIT_T1] closes at
    4.12 a.u. The two do not overlap for v <= 3.0. Fitting over the localised
    window alone measures the packet ACCELERATING down the slab's attractive
    gradient while still in vacuum, which is why it returns a negative S.
    Both windows are reported in the notebooks; neither is silently preferred.
    """
    return (-SLAB_HALF - z0) / v, (SLAB_HALF - z0) / v


def fit_all(run: K.WPRun, t0: float = FIT_T0, t1: float = FIT_T1):
    """The four S_ij = -dT_i/ds_j over the given window."""
    return K.fit_all_wp(run, t0, t1)


def fit_all_in_slab(run: K.WPRun, v: float):
    """The four S_ij fitted over the in-slab transit (see slab_window)."""
    t0, t1 = slab_window(v)
    return K.fit_all_wp(run, t0, min(t1, run.t[-1]))


def load_interactions(run_dir: Path) -> pd.DataFrame:
    """Pairwise Coulomb ledger, with the closure residuals attached."""
    df = K._concat_segments(Path(run_dir) / "raw" / "observables", "interactions")
    obs = K._concat_segments(Path(run_dir) / "raw" / "observables", "observables")
    m = pd.merge(df, obs, on="step", suffixes=("", "_obs"))
    # compute_coulomb_wp closure (interaction_energies.hpp):
    #   E_hartree = E_SS + E_PS + E_PP ;  E_external = E_SB + E_PB
    if "energy_hartree" in m:
        m["hartree_residual"] = m["e_hartree_check"] - m["energy_hartree"]
        m["hartree_from_parts"] = m["e_ss"] + m["e_ps"] + m["e_pp"]
    if "energy_external" in m:
        m["external_residual"] = m["e_external_check"] - m["energy_external"]
        m["external_from_parts"] = m["e_sb"] + m["e_pb"]
    return m


def wp_kinetic_norm_correction(run_dir: Path) -> pd.DataFrame:
    """
    Correct INQ's energy ledger for its norm-divided kinetic term.

    INQ reports each orbital's kinetic energy as occ*<psi|T|psi>/<psi|psi>
    (inq/src/hamiltonian/energy.hpp:50-55, used at :83). Every other term is
    density-based and extensive. Under a CAP the WP orbital's norm decays, so its
    kinetic contribution keeps reporting the per-particle MEAN instead of leaving
    the ledger — energy_kinetic, and hence energy_total, are inflated.

    VERIFIED 2026-07-30 that this is still true on this machine: inq-study's
    energy.hpp is BYTE-IDENTICAL to stock inq's (the inq-study fork carries only
    the muon per-state mass work, the CAP complexification in self_consistency,
    and absorbing_monomial.hpp). A real-time column is being added upstream; until
    that lands, this post-processing route is the correction, and afterwards the
    two can be cross-checked against each other.

        E_total_corrected(t) = E_total_reported(t) - occ * <T>(t) * (1/norm(t) - 1)
                             = E_total_reported(t) - occ * T1(t) * (1 - norm(t))

    since T1 (`e_kin_ha` from WPMomentumStats) IS the norm-divided mean <T>/norm.
    occ = 1 (inject_into_last_extra_state with occupation 1.0).

    Everything needed is written EVERY step by these runs, so unlike the original
    `wp_cap_energy_plateau/wp_kinetic_normalization_fix.py` (which reconstructed
    <T> and the norm from ~100 sparse wavefunction VTIs) this is exact at full
    cadence and needs no VTI reconstruction.

    NOTE the two `norm_check` columns are different quantities: use the REAL-SPACE
    one (wp_real_space_stats, ~1 at t=0), not the momentum-space Parseval constant.

    Returns a frame with time_au, norm_wp, T1_ev, correction_ev,
    e_total_raw_ev, e_total_corrected_ev.
    """
    obs_dir = Path(run_dir) / "raw" / "observables"
    mom = K._concat_segments(obs_dir, "wp_momentum_stats")
    pos = K._concat_segments(obs_dir, "wp_real_space_stats")
    obs = K._concat_segments(obs_dir, "observables")

    df = pd.merge(mom[["step", "time_au", "e_kin_ha"]],
                  pos[["step", "norm_check"]], on="step")
    df = pd.merge(df, obs[["step", "energy_total"]], on="step")

    occ = 1.0
    df["norm_wp"] = df["norm_check"]
    df["T1_ev"] = df["e_kin_ha"] * HA_TO_EV
    df["correction_ev"] = occ * df["T1_ev"] * (1.0 - df["norm_wp"])
    df["e_total_raw_ev"] = df["energy_total"] * HA_TO_EV
    df["e_total_corrected_ev"] = df["e_total_raw_ev"] - df["correction_ev"]
    # Bare (extensive) WP kinetic content — the independent cross-check: its
    # decrease must track the corrected total's drift, up to what the bath absorbs.
    df["wp_kinetic_bare_ev"] = df["T1_ev"] * df["norm_wp"]
    return df[["step", "time_au", "norm_wp", "T1_ev", "correction_ev",
               "e_total_raw_ev", "e_total_corrected_ev", "wp_kinetic_bare_ev"]]


# Production ground state (dx = 0.40) — the bath these runs actually load. NOT the
# dx = 0.50 fidelity value 207.18322156141 Ha that the CLASSICAL campaign used;
# the two agree to 9e-6 Ha, but the deposit must be referenced to the GS the run
# was started from.
E_GS_HA_DX040 = 207.18323030158


def deposit_stopping(v: float, e_gs_ha: float = E_GS_HA_DX040,
                     sigma: float | None = None) -> dict[str, float]:
    """
    The CLASSICAL Definition-2 estimator applied to a wavepacket run:

        S = (E_total(t_final) - E_GS) / L_slab      L_slab = 25 Bohr

    Returned both from INQ's raw ledger and from the norm-corrected one
    (see wp_kinetic_norm_correction).

    READ THE CAVEAT. This estimator is clean CLASSICALLY because there the
    projectile is an EXTERNAL perturbation: it is never part of the electronic
    ledger, so plateau - E_GS is purely the slab's energy gain, and the CAP-free
    z-open box lets the projectile leave without taking ledger energy with it.

    Neither holds here. The wavepacket IS part of the system, and the CAP removes
    it — so E_total(t_final) - E_GS is "what is left in the box", i.e. the bath's
    retained excitation plus whatever of the packet has not yet been absorbed,
    with everything the CAP took already subtracted. It is therefore a LOWER
    BOUND on the deposit, not the classical quantity, and the two should not be
    read as the same measurement.

    The raw-vs-corrected contrast is a diagnostic in its own right: on the raw
    ledger this estimator comes out velocity-INDEPENDENT at ~2.44 eV/Bohr, which
    is unphysical (S must fall with v in this regime) and is the signature of the
    norm-divided kinetic term dominating the residual. After the correction it
    falls monotonically with v, as a stopping power must.
    """
    d = wp_kinetic_norm_correction(WP_RESULTS / name_for(v, sigma))
    e_gs_ev = e_gs_ha * HA_TO_EV
    ef_raw = float(d.e_total_raw_ev.iloc[-1]) - e_gs_ev
    ef_cor = float(d.e_total_corrected_ev.iloc[-1]) - e_gs_ev

    # COMPLETENESS GATE. "t_final" is only t_final if the run actually finished.
    # Reading a still-propagating (or killed) run gives a mid-flight snapshot in
    # which the packet has not yet deposited and the CAP has not yet removed it,
    # so E_total - E_GS is nowhere near the deposit — and it looks perfectly
    # plausible. Caught live on 2026-07-31 while the sigma = 2 array was running:
    # 86 of 3623 steps produced S_deposit = 2.35 eV/Bohr with norm_final = 1.000.
    # Callers must check `complete` before quoting the number.
    steps_done = int(d.step.iloc[-1])
    target = N_STEPS.get(v)
    complete = target is not None and steps_done >= target
    return {
        "E_deposit_raw_eV": ef_raw,
        "E_deposit_corrected_eV": ef_cor,
        "S_deposit_raw": ef_raw / (2.0 * SLAB_HALF),
        "S_deposit_corrected": ef_cor / (2.0 * SLAB_HALF),
        "norm_final": float(d.norm_wp.iloc[-1]),
        "t_final_au": float(d.time_au.iloc[-1]),
        "steps_done": steps_done,
        "steps_target": target if target is not None else -1,
        "complete": bool(complete),
    }


def classical_reference() -> pd.DataFrame:
    """The six classical benchmark points this study is compared against."""
    return pd.read_csv(CLASSICAL_SUMMARY)


_KV = re.compile(r"(\w+)\s*=\s*([^=]*?)(?=\s+\w+\s*=|$)")


def run_summary(run_dir: Path) -> dict[str, str]:
    """
    Parse run_summary.txt into a dict.

    MUST use a regex, NOT split-on-first-'='. run_summary.txt packs SEVERAL
    `key = value` pairs onto one line, e.g.

        save_every = 7  wf_every = 21  stats_every = 1  ckpt_every = 414

    so partitioning on the first '=' returns
    "7  wf_every = 21  stats_every = 1  ckpt_every = 414" for save_every, and
    int() on it throws. This is the same trap the classical campaign hit and
    documented ("use regex grab(), not split-on-first-'='"); it is re-fixed here
    rather than re-learned.
    """
    out: dict[str, str] = {}
    p = Path(run_dir) / "run_summary.txt"
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        for k, val in _KV.findall(line):
            out[k.strip()] = val.strip()
    return out
