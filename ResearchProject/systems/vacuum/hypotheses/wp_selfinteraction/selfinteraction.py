"""Quantifying the wavepacket self-interaction in vacuum.

Plan: docs/plans/wp-self-interaction-correction.md
Runs: ResearchProject/systems/vacuum/scripts/wp_selfinteraction/results/
      {noninteracting,hartree,lda,sic_h,sic_pzrun}

THE MEASUREMENT
---------------
A single electron has, exactly, no self-interaction. So one electron alone in a
vacuum box must follow free-particle dispersion, which is known in closed form.
Running the SAME initial Gaussian at three theory levels therefore measures the
self-interaction by DIFFERENCE, with no self-interaction correction implemented:

    noninteracting   no self-interaction         <- the reference
    hartree          + Hartree self-interaction
    lda              + Hartree AND LDA xc self-interaction

    lda - noninteracting      the TOTAL self-interaction error
    hartree - noninteracting  its HARTREE part
    lda - hartree             its XC part

Two further runs propagate under LDA but apply the projected SIC kick
(inqkit::SelfInteractionCorrection), turning the difference MEASUREMENT into an
INTERVENTION test:

    sic_h            LDA + Hartree-only correction   (predicted to OVER-correct)
    sic_pzrun        LDA + full run-consistent SIC   (must match noninteracting)

`sic_pzrun` is gated with the SAME numerics_gate as the reference — the Tier V
acceptance criterion of the plan: with the self-field removed the packet is free
again, so the closed form must re-emerge to the accuracy of the grid itself.

The `noninteracting` run is not merely a control: it carries the SAME grid, the
same propagator, the same time step and the same injected packet as the other
two, so comparing against it (rather than against the analytic formula) cancels
discretisation error to first order. The analytic formula is used separately, to
check that the numerics are converged at all.

THE CLOSED-FORM REFERENCE (atomic units, m_e = 1)
-------------------------------------------------
For psi_0 ~ exp(-r^2 / (2 sigma^2)) the DENSITY |psi|^2 has per-axis std

    sigma_dens(t) = sqrt( sigma^2/2 + t^2/(2 sigma^2) )

and the momentum distribution never changes at all:

    <p_d> = k0_d          var(p_d) = 1/(2 sigma^2)      for all t.

var(p) is the sharpest of the three gates because it is EXACTLY conserved under
free evolution — a growing var(p) cannot be blamed on discretisation of a
spreading packet, because nothing about it is supposed to change.

WHAT E_PP IS HERE
-----------------
In vacuum the wavepacket is the only charge in the box, so the pairwise
decomposition collapses to a single non-zero term and INQ's own `energy_hartree`
must equal our offline `E_PP` exactly for the interacting runs. That identity is
a closure gate, checked by `closure()`. For the `noninteracting` run INQ reports
`energy_hartree = 0` while the offline `E_PP` stays non-zero — there it is a
pure diagnostic of the packet's size, not an energy the packet feels.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]

HA_TO_EV = 27.211386245988

RESULTS = Path(os.environ.get(
    "WPSI_RESULTS",
    REPO / "ResearchProject/systems/vacuum/scripts/wp_selfinteraction/results"))

THEORIES = ("noninteracting", "hartree", "lda", "sic_h", "sic_pzrun")
# runs whose propagation feels an interaction (everything but the reference),
# in the order the plots and tables present them
INTERACTING = ("hartree", "lda", "sic_h", "sic_pzrun")
LABEL = {"noninteracting": "non-interacting (reference)",
         "hartree": "Hartree self-interaction",
         "lda": "Hartree + LDA xc self-interaction",
         "sic_h": "LDA + Hartree-only SIC (over-corrects)",
         "sic_pzrun": "LDA + full SIC (corrected)"}
COLOR = {"noninteracting": "tab:blue", "hartree": "tab:orange", "lda": "tab:red",
         "sic_h": "tab:purple", "sic_pzrun": "tab:green"}


# ---------------------------------------------------------------------------
# Closed-form free-particle reference
# ---------------------------------------------------------------------------

def sigma_dens_free(t, sigma_wp: float) -> np.ndarray:
    """Per-axis density std of a freely dispersing Gaussian, Bohr."""
    t = np.asarray(t, dtype=float)
    return np.sqrt(sigma_wp**2 / 2.0 + t**2 / (2.0 * sigma_wp**2))


def var_p_free(sigma_wp: float) -> float:
    """Per-axis momentum variance. CONSTANT under free evolution."""
    return 1.0 / (2.0 * sigma_wp**2)


def localisation_ev(sigma_wp: float) -> float:
    """(3/2) sum_d var(p_d)/2 = 3/(4 sigma^2), in eV. Constant when free."""
    return 3.0 / (4.0 * sigma_wp**2) * HA_TO_EV


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def summary_of(run_dir: Path) -> dict:
    out: dict[str, str] = {}
    p = Path(run_dir) / "run_summary.txt"
    if not p.is_file():
        return out
    for line in p.read_text().splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    return out


@dataclass
class SIRun:
    """One theory level of the vacuum self-interaction experiment."""
    theory: str
    run_dir: Path
    sigma_wp: float
    t: np.ndarray
    # real space (circular estimators; the packet is centred so they agree with
    # the naive ones until it approaches the box face — their divergence IS the
    # wrap diagnostic, exposed as `wrap_indicator`)
    sigma_x: np.ndarray
    sigma_y: np.ndarray
    sigma_z: np.ndarray
    sigma_naive_z: np.ndarray
    # momentum
    pz: np.ndarray
    var_px: np.ndarray
    var_py: np.ndarray
    var_pz: np.ndarray
    e_kin_ha: np.ndarray
    # energies
    e_total: np.ndarray
    e_hartree: np.ndarray
    e_xc: np.ndarray
    e_kin_inq: np.ndarray
    # pairwise
    e_pp: np.ndarray
    closure_pp: np.ndarray
    norm_wp: np.ndarray
    complete: bool
    # SIC diagnostics from sic.csv — None when no correction was active.
    # e_corrected = E_KS - U_self - Exc_self is the quantity that IS conserved
    # under an active correction (E_total is not, by design — plan section 0/D2).
    u_self: np.ndarray | None = None
    exc_self: np.ndarray | None = None
    e_corrected: np.ndarray | None = None
    cum_norm_removed: np.ndarray | None = None

    @property
    def sigma_iso(self) -> np.ndarray:
        """Geometric-mean density width, (sx sy sz)^(1/3).

        The packet is isotropic when k0 = 0, so this is just a lower-noise
        estimate of the single width; for a moving packet it is the natural
        scalar summary of an anisotropic cloud.
        """
        return (self.sigma_x * self.sigma_y * self.sigma_z) ** (1.0 / 3.0)

    @property
    def var_p3d(self) -> np.ndarray:
        return self.var_px + self.var_py + self.var_pz

    @property
    def var_term_ev(self) -> np.ndarray:
        """var(p)/2m in eV — the internal (non-drift) kinetic energy."""
        return 0.5 * self.var_p3d * HA_TO_EV

    @property
    def e_pp_ev(self) -> np.ndarray:
        return self.e_pp * HA_TO_EV

    @property
    def wrap_indicator(self) -> np.ndarray:
        """|sigma_naive_z - sigma_z| / sigma_z.

        The naive and circular second moments agree for a packet well inside the
        box and diverge once density reaches the face. Non-zero values mean the
        run has started to wrap and later times are not trustworthy.
        """
        return np.abs(self.sigma_naive_z - self.sigma_z) / self.sigma_z

    def sigma_free(self) -> np.ndarray:
        return sigma_dens_free(self.t, self.sigma_wp)

    def excess_spreading(self) -> np.ndarray:
        """Measured width / free-dispersion width. 1.0 means no excess."""
        return self.sigma_iso / self.sigma_free()


def _read(obs: Path, stem: str) -> pd.DataFrame:
    files = sorted(obs.glob(f"{stem}*.csv"))
    if not files:
        raise FileNotFoundError(f"no {stem}*.csv under {obs}")
    df = pd.concat([pd.read_csv(f, comment="#") for f in files], ignore_index=True)
    return df.sort_values("step").drop_duplicates(subset="step",
                                                  keep="last").reset_index(drop=True)


def load(theory: str, name: str | None = None) -> SIRun:
    if theory not in THEORIES:
        raise ValueError(f"theory must be one of {THEORIES}, got {theory!r}")
    run_dir = RESULTS / (name or theory)
    obs = run_dir / "raw" / "observables"
    if not obs.is_dir():
        raise FileNotFoundError(f"no observables under {run_dir}")

    smry = summary_of(run_dir)
    sigma_wp = float(smry.get("wp", "").split("sigma")[1].split()[0]) if "sigma" in smry.get("wp", "") else 4.0

    mom = _read(obs, "wp_momentum_stats")
    pos = _read(obs, "wp_real_space_stats")
    en = _read(obs, "energies")
    ie = _read(obs, "interactions")

    df = mom.merge(pos, on=["step", "time_au"], suffixes=("_p", "_r"))
    df = df.merge(en, on=["step", "time_au"], how="left")
    df = df.merge(ie[["step", "e_pp", "closure_pp_minus_hartree", "norm_wp"]],
                  on="step", how="left")

    # sic.csv exists for every run of the extended binary but carries data rows
    # only when a correction was active (header-only otherwise)
    have_sic = False
    try:
        sic = _read(obs, "sic")
        if len(sic):
            df = df.merge(sic[["step", "u_self_ha", "exc_self_ha",
                               "e_corrected_ha", "cum_norm_removed"]],
                          on="step", how="left")
            have_sic = True
    except FileNotFoundError:
        pass

    return SIRun(
        theory=theory, run_dir=run_dir, sigma_wp=sigma_wp,
        t=df["time_au"].to_numpy(),
        sigma_x=df["sigma_x_circ"].to_numpy(),
        sigma_y=df["sigma_y_circ"].to_numpy(),
        sigma_z=df["sigma_z_circ"].to_numpy(),
        sigma_naive_z=np.sqrt(np.maximum(df["sigma_z2"].to_numpy(), 0.0)),
        pz=df["pz_mean"].to_numpy(),
        var_px=df["sigma_px2"].to_numpy(),
        var_py=df["sigma_py2"].to_numpy(),
        var_pz=df["sigma_pz2"].to_numpy(),
        e_kin_ha=df["e_kin_ha"].to_numpy(),
        e_total=df["total"].to_numpy(),
        e_hartree=df["hartree"].to_numpy(),
        e_xc=df["xc"].to_numpy(),
        e_kin_inq=df["kinetic"].to_numpy(),
        e_pp=df["e_pp"].to_numpy(),
        closure_pp=df["closure_pp_minus_hartree"].to_numpy(),
        norm_wp=df["norm_wp"].to_numpy(),
        complete=smry.get("run_completed", "false").lower() == "true",
        u_self=df["u_self_ha"].to_numpy() if have_sic else None,
        exc_self=df["exc_self_ha"].to_numpy() if have_sic else None,
        e_corrected=df["e_corrected_ha"].to_numpy() if have_sic else None,
        cum_norm_removed=df["cum_norm_removed"].to_numpy() if have_sic else None,
    )


def load_all(suffix: str = "") -> dict[str, SIRun]:
    """Load every theory whose run directory exists.

    Tolerant of absent runs (e.g. the `_smoke` suffix exists only for the three
    original theories) — but the reference is mandatory, since every measurement
    here is a difference against it.
    """
    out = {th: load(th, name=th + suffix)
           for th in THEORIES if (RESULTS / (th + suffix)).is_dir()}
    if "noninteracting" not in out:
        raise FileNotFoundError(f"reference run 'noninteracting{suffix}' "
                                f"missing under {RESULTS}")
    return out


# ---------------------------------------------------------------------------
# Gates
# ---------------------------------------------------------------------------

def numerics_gate(ref: SIRun, tol_sigma: float = 0.005,
                  tol_varp: float = 1e-3) -> dict:
    """Does the NON-INTERACTING run reproduce the analytic free solution?

    This licenses everything else. If the reference run does not disperse at the
    analytic rate, the grid or the time step is not converged and any
    "self-interaction" extracted by difference is contaminated by that instead.
    """
    sig_meas = ref.sigma_iso
    sig_exact = ref.sigma_free()
    rel = np.abs(sig_meas / sig_exact - 1.0)
    vp0 = var_p_free(ref.sigma_wp)
    vdrift = np.abs(ref.var_pz / vp0 - 1.0)
    return {
        "max_rel_sigma_error": float(rel.max()),
        "sigma_ok": bool(rel.max() < tol_sigma),
        "max_var_p_drift": float(vdrift.max()),
        "var_p_ok": bool(vdrift.max() < tol_varp),
        "e_total_drift_ev": float((ref.e_total[-1] - ref.e_total[0]) * HA_TO_EV),
        "max_wrap_indicator": float(ref.wrap_indicator.max()),
        "passed": bool(rel.max() < tol_sigma and vdrift.max() < tol_varp),
    }


def closure(run: SIRun, tol_ha: float = 1e-8) -> dict:
    """Offline E_PP vs INQ's own energy_hartree.

    In vacuum the packet is the ONLY charge, so for an interacting theory these
    must be the same number. For `noninteracting` INQ reports 0 by construction
    and the comparison is not meaningful — reported, not gated.
    """
    resid = np.abs(run.closure_pp)
    gated = run.theory != "noninteracting"
    return {
        "theory": run.theory,
        "max_abs_residual_ha": float(resid.max()),
        "gated": gated,
        "passed": bool((not gated) or resid.max() < tol_ha),
    }


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------

@dataclass
class SIEffect:
    """Self-interaction, extracted by difference against the reference run."""
    theory: str
    t: np.ndarray
    sigma_ratio: np.ndarray        # width relative to the non-interacting run
    sigma_ratio_free: np.ndarray   # width relative to the analytic free solution
    d_var_term_ev: np.ndarray      # excess internal kinetic energy vs reference
    e_pp_ev: np.ndarray            # the self-Hartree energy itself
    d_e_pp_ev: np.ndarray          # its change from t=0

    def at(self, t_query: float) -> dict:
        i = int(np.argmin(np.abs(self.t - t_query)))
        return {"t": float(self.t[i]),
                "sigma_ratio_vs_reference": float(self.sigma_ratio[i]),
                "sigma_ratio_vs_free": float(self.sigma_ratio_free[i]),
                "excess_var_energy_ev": float(self.d_var_term_ev[i]),
                "e_pp_ev": float(self.e_pp_ev[i]),
                "d_e_pp_ev": float(self.d_e_pp_ev[i])}


def effect(run: SIRun, ref: SIRun) -> SIEffect:
    """Isolate the self-interaction of `run` against the non-interacting `ref`.

    Ratios are taken against the REFERENCE RUN rather than the analytic formula
    so that grid and propagator error — identical between the two by
    construction — cancels. The analytic ratio is carried alongside so the two
    can be seen not to disagree.
    """
    n = min(run.t.size, ref.t.size)
    return SIEffect(
        theory=run.theory,
        t=run.t[:n],
        sigma_ratio=run.sigma_iso[:n] / ref.sigma_iso[:n],
        sigma_ratio_free=run.sigma_iso[:n] / run.sigma_free()[:n],
        d_var_term_ev=run.var_term_ev[:n] - ref.var_term_ev[:n],
        e_pp_ev=run.e_pp_ev[:n],
        d_e_pp_ev=run.e_pp_ev[:n] - run.e_pp_ev[0],
    )


def summary_table(runs: dict[str, SIRun], t_query: float | None = None) -> pd.DataFrame:
    """One row per interacting theory: the self-interaction it introduces."""
    ref = runs["noninteracting"]
    t_query = float(ref.t[-1]) if t_query is None else t_query
    rows = []
    for th in INTERACTING:
        if th not in runs:
            continue
        e = effect(runs[th], ref)
        a = e.at(t_query)
        rows.append({
            "theory": th,
            "t": a["t"],
            "sigma_ratio_vs_reference": a["sigma_ratio_vs_reference"],
            "excess_width_pct": 100.0 * (a["sigma_ratio_vs_reference"] - 1.0),
            "E_PP_start_eV": float(e.e_pp_ev[0]),
            "E_PP_end_eV": a["e_pp_ev"],
            "E_PP_released_eV": -a["d_e_pp_ev"],
            "excess_var_energy_eV": a["excess_var_energy_ev"],
        })
    df = pd.DataFrame(rows)
    if {"hartree", "lda"} <= set(df.get("theory", [])):
        h = df[df.theory == "hartree"].iloc[0]
        l = df[df.theory == "lda"].iloc[0]
        df = pd.concat([df, pd.DataFrame([{
            "theory": "xc part (lda - hartree)",
            "t": l["t"],
            "sigma_ratio_vs_reference": float("nan"),
            "excess_width_pct": l["excess_width_pct"] - h["excess_width_pct"],
            "E_PP_start_eV": float("nan"), "E_PP_end_eV": float("nan"),
            "E_PP_released_eV": float("nan"),
            "excess_var_energy_eV": l["excess_var_energy_eV"] - h["excess_var_energy_eV"],
        }])], ignore_index=True)
    return df


# ---------------------------------------------------------------------------
# The link back to the channeling study
# ---------------------------------------------------------------------------

# Measured on the channeling runs at t = 30 a.u., sigma_WP = 4.
#
# DEFINITION DISCIPLINE — this bit was wrong once and the error was invisible.
# All THREE numbers below use the SAME width definition as this vacuum analysis:
# the 3-D geometric-mean circular width sigma_iso = (sx sy sz)^(1/3), divided by
# the analytic free sigma_dens(30) = 6.0104 Bohr. An earlier version compared the
# vacuum 3-D ratio against 1.467, which is the TRANSVERSE <r_perp>/free ratio
# from the 5th handover entry — mixing a 3-D number with a transverse one
# understated the fraction (16.9 % instead of 20.9 %). The channeling packet is
# strongly anisotropic (transverse 1.52, longitudinal 1.13), so the two
# definitions genuinely differ and must never be crossed.
#
# Provenance:
#   scripts/channeling_twin/wp/results/wp        (uncorrected LDA)
#   scripts/channeling_sic/wp/results/wp_sic     (full Perdew-Zunger SIC, job 32615191)
CHANNELING_EXCESS = 1.3778              # uncorrected LDA wavepacket, 3-D iso
CHANNELING_EXCESS_SIC = 1.2990          # same run with the self-interaction REMOVED
CHANNELING_EXCESS_TRANSVERSE = 1.467    # <r_perp>/free — cross-reference only
CHANNELING_T_END = 30.0


def channeling_comparison(runs: dict[str, SIRun]) -> dict:
    """How much of the channeling run's excess spreading is self-interaction?

    Two independent answers, which is the point of reporting both:

    - `*_fraction_of_channeling` — PREDICTED from this vacuum experiment, by
      assuming the vacuum self-interaction carries over unchanged into the bore.
      That assumption is not obviously safe: LDA evaluates v_xc at the TOTAL
      local density, which in the bore is bath + packet, so the ~82 % xc
      cancellation measured in vacuum need not hold there.
    - `sic_measured_fraction` — MEASURED directly, by re-running the channeling
      case with full Perdew-Zunger SIC and differencing. No transfer assumption.

    The two agreeing is what licenses the vacuum number as transferable. They are
    kept separate so that agreement stays a result rather than an input.
    """
    ref = runs["noninteracting"]
    excess = CHANNELING_EXCESS - 1.0
    out = {
        "channeling_excess_vs_free": CHANNELING_EXCESS,
        "channeling_excess_vs_free_SIC": CHANNELING_EXCESS_SIC,
        "sic_measured_fraction": (
            (CHANNELING_EXCESS - CHANNELING_EXCESS_SIC) / excess
            if excess > 0 else float("nan")),
    }
    for th in ("hartree", "lda"):
        if th not in runs:
            continue
        e = effect(runs[th], ref)
        i = int(np.argmin(np.abs(e.t - CHANNELING_T_END)))
        vac = float(e.sigma_ratio[i])
        out[f"{th}_vacuum_excess"] = vac
        out[f"{th}_fraction_of_channeling"] = (
            (vac - 1.0) / excess if excess > 0 else float("nan"))
    return out
