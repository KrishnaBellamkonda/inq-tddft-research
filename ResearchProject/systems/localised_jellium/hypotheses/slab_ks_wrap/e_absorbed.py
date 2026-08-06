"""
E_absorbed stopping power for the localised jellium slab.

THE DEFINITION (user, 2026-08-01 -- the primary measure for this project):

    S  =  E_absorbed / L_slab_z            L_slab_z = 2 * SLAB_HALF_WIDTH = 25 Bohr
    E_absorbed  =  E_total(t_final) - E_GS

E_total is INQ's electronic total energy (observables.csv `energy_total`); E_GS is
the converged ground state the run was launched from. L_slab_z is the slab's
extent along the propagation axis -- the path length over which the projectile
can deposit energy.

WHY THIS MODULE EXISTS. The number was previously computed only inside the
wp_highdensity_sv synthesis notebooks, as one line with no derivation shown. It is
the measure the project compares on, so it gets its own module with the
conventions, the caveats and the closure evidence in one place.

--------------------------------------------------------------------------------
THREE THINGS THAT MUST BE RIGHT, AND ARE EASY TO GET WRONG
--------------------------------------------------------------------------------
1. SEGMENTS. Resumed runs write observables.from<N>.csv alongside observables.csv.
   The final energy is the last row IN STEP ORDER across all segments, not the
   last row of the base file. Segments are concatenated and de-duplicated here.

2. E_GS IS PER DENSITY AND PER GRID. Using the wrong one silently shifts every S
   by a constant. The two slab_ks_wrap ground states differ by 176 Ha:
       n40  (r_s 5.67)   E_GS =  31.5295278631 Ha
       n100 (r_s 4.18)   E_GS = 207.1832303016 Ha
   Both at dx = 0.40, periodicity(2), L = 35x35x85.

3. THE CAP MAKES E_total NON-CONSERVED BY CONSTRUCTION. These runs carry two
   absorbing bands (12.5 Bohr per z face, eta = -1 Ha). E_total(t_final) - E_GS is
   therefore "what is still in the box at t_final", i.e. the medium's retained
   excitation plus any unabsorbed projectile, with everything the CAP removed
   already subtracted. For a WAVEPACKET run the projectile is part of the
   electronic ledger and the CAP eventually removes it, so this is a LOWER BOUND
   on the energy the medium absorbed. For a CLASSICAL run the projectile is an
   external perturbation that was never in the ledger, so the same expression is
   the medium's gain directly. The two are therefore NOT the same estimator even
   though they share a formula -- see `closure` in the returned record.

Additionally, for wavepacket runs INQ divides the orbital kinetic energy by the
orbital norm (inq/src/hamiltonian/energy.hpp), so as the CAP eats the packet the
kinetic term is inflated by 1/norm. `corrected` removes that; `raw` keeps it. On
the raw ledger the estimator comes out velocity-INDEPENDENT (~2.44 eV/Bohr at
sigma = 0.5), which no stopping power can be -- that flatness is the signature of
the artefact, and the reason `corrected` is the number to quote.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HA_TO_EV = 27.211386245988
SLAB_HALF = 12.5
L_SLAB_Z = 2.0 * SLAB_HALF          # 25 Bohr, the deposition path length

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SCRIPTS = REPO / "ResearchProject/systems/localised_jellium/scripts/slab_ks_wrap"

# Ground states, read from the GS run logs (slabks-gs-32484286/7). Per density.
E_GS_HA = {"n40": 31.5295278631, "n100": 207.1832303016}
# n0 = N / (35*35*25); r_s = (3/(4 pi n0))^(1/3)
R_S = {"n40": 5.674, "n100": 4.183}
DENSITIES = ("n100", "n40")
VELOCITIES = (2.0, 2.5, 3.0, 3.5)
HALVES = ("classical", "wp")

_SEG = re.compile(r"\.from(\d+)\.csv$")


def _concat(obs: Path, stem: str) -> pd.DataFrame:
    """Base CSV plus every .from<N> segment, ordered by step, duplicates dropped.

    A resumed run recomputes nothing, but a run that was rewound to a checkpoint
    can legitimately have its base file end after a segment begins; keeping the
    LAST occurrence of each step takes the most recent computation of it.
    """
    parts = []
    base = obs / f"{stem}.csv"
    if base.exists():
        parts.append(pd.read_csv(base, comment="#"))
    for f in sorted(obs.glob(f"{stem}.from*.csv"),
                    key=lambda p: int(_SEG.search(p.name).group(1))):
        parts.append(pd.read_csv(f, comment="#"))
    if not parts:
        raise FileNotFoundError(f"no {stem}.csv under {obs}")
    df = pd.concat(parts, ignore_index=True)
    return (df.drop_duplicates(subset="step", keep="last")
              .sort_values("step").reset_index(drop=True))


@dataclass(frozen=True)
class Absorbed:
    """One run's E_absorbed measurement, with the evidence needed to trust it."""
    run_dir: Path
    half: str
    density: str
    v: float
    r_s: float
    e_gs_ev: float
    e_final_ev: float           # raw ledger
    e_final_corr_ev: float      # norm-corrected (wp only; == raw for classical)
    E_absorbed_eV: float        # corrected - E_GS   <- the headline
    E_absorbed_raw_eV: float
    S_eV_per_Bohr: float        # E_absorbed / L_SLAB_Z
    S_raw_eV_per_Bohr: float
    t_final_au: float
    steps_done: int
    norm_final: float           # wp only; NaN for classical
    plateau_drift_eV: float     # |E(t_f) - E(0.9 t_f)|: has it settled?

    @property
    def settled(self) -> bool:
        """Plateau criterion: the last 10 % of the run moves E by < 1 % of E_absorbed."""
        return abs(self.plateau_drift_eV) < 0.01 * max(abs(self.E_absorbed_eV), 1e-12)


def measure(half: str, density: str, v: float,
            scripts: Path = SCRIPTS, suffix: str = "") -> Absorbed:
    """E_absorbed for one slab_ks_wrap run (suffix="_cap" for the CAP'd twins)."""
    name = f"{density}_v" + f"{v:.1f}".replace(".", "p") + suffix
    run_dir = scripts / half / "results" / name
    return measure_dir(run_dir, E_GS_HA[density], half, density, v, R_S[density])


def measure_dir(run_dir: Path, e_gs_ha: float, half: str,
                density: str, v: float, r_s: float) -> Absorbed:
    """Core measurement for an arbitrary run directory (any sweep with the
    standard raw/observables layout, e.g. wp_highdensity_sv)."""
    obs = run_dir / "raw" / "observables"
    d = _concat(obs, "observables")

    e_ev = d["energy_total"].to_numpy() * HA_TO_EV
    t = d["time_au"].to_numpy()
    e_gs_ev = e_gs_ha * HA_TO_EV

    # Norm correction: undo INQ's division of the WP orbital kinetic energy by its
    # (CAP-decaying) norm. Classical runs have no WP orbital, so raw == corrected.
    norm_final = float("nan")
    e_corr = e_ev.copy()
    if half == "wp":
        mom = _concat(obs, "wp_momentum_stats")
        pos = _concat(obs, "wp_real_space_stats")
        m = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))
        m = m[m.step.isin(d.step)]
        norm = (m["norm_check_r"] if "norm_check_r" in m else m["norm_check"]).to_numpy()
        T1_ev = m["e_kin_ha"].to_numpy() * HA_TO_EV
        n = min(len(e_ev), len(norm))
        # occupation 1 for the injected WP state; E_corr = E_raw - T1*(1 - norm)
        e_corr = e_ev.copy()
        e_corr[:n] = e_ev[:n] - T1_ev[:n] * (1.0 - norm[:n])
        norm_final = float(norm[n - 1])

    i90 = int(0.9 * (len(e_corr) - 1))
    return Absorbed(
        run_dir=run_dir, half=half, density=density, v=v, r_s=r_s,
        e_gs_ev=e_gs_ev,
        e_final_ev=float(e_ev[-1]), e_final_corr_ev=float(e_corr[-1]),
        E_absorbed_eV=float(e_corr[-1]) - e_gs_ev,
        E_absorbed_raw_eV=float(e_ev[-1]) - e_gs_ev,
        S_eV_per_Bohr=(float(e_corr[-1]) - e_gs_ev) / L_SLAB_Z,
        S_raw_eV_per_Bohr=(float(e_ev[-1]) - e_gs_ev) / L_SLAB_Z,
        t_final_au=float(t[-1]), steps_done=int(d["step"].iloc[-1]),
        norm_final=norm_final,
        plateau_drift_eV=float(e_corr[-1] - e_corr[i90]),
    )


def table(scripts: Path = SCRIPTS, suffix: str = "",
          halves: tuple = HALVES) -> pd.DataFrame:
    """Every slab_ks_wrap run: 2 densities x 4 velocities x {classical, wp}.

    suffix="_cap" selects the CAP'd WP twins (halves=("wp",) — no classical
    _cap runs exist; the classical half is CAP-free by design)."""
    rows = []
    for half in halves:
        for dens in DENSITIES:
            for v in VELOCITIES:
                try:
                    a = measure(half, dens, v, scripts, suffix)
                except (FileNotFoundError, KeyError) as e:
                    print(f"  MISSING {half}/{dens}{suffix}/v={v}: {type(e).__name__}")
                    continue
                rows.append({
                    "half": a.half, "density": a.density, "r_s": a.r_s, "v": a.v,
                    "E_absorbed_eV": a.E_absorbed_eV,
                    "S_eV_per_Bohr": a.S_eV_per_Bohr,
                    "S_raw_eV_per_Bohr": a.S_raw_eV_per_Bohr,
                    "t_final_au": a.t_final_au, "steps_done": a.steps_done,
                    "norm_final": a.norm_final,
                    "plateau_drift_eV": a.plateau_drift_eV, "settled": a.settled,
                })
    return pd.DataFrame(rows)


def energy_trace(half: str, density: str, v: float,
                 scripts: Path = SCRIPTS, suffix: str = "") -> pd.DataFrame:
    """t, E_total raw and corrected, in eV relative to E_GS — for plotting."""
    name = f"{density}_v" + f"{v:.1f}".replace(".", "p") + suffix
    obs = scripts / half / "results" / name / "raw" / "observables"
    d = _concat(obs, "observables")
    e_ev = d["energy_total"].to_numpy() * HA_TO_EV
    e_gs_ev = E_GS_HA[density] * HA_TO_EV
    out = pd.DataFrame({"t": d["time_au"].to_numpy(),
                        "step": d["step"].to_numpy(),
                        "dE_raw": e_ev - e_gs_ev})
    if half == "wp":
        mom = _concat(obs, "wp_momentum_stats")
        pos = _concat(obs, "wp_real_space_stats")
        m = pd.merge(mom, pos, on=["step", "time_au"], suffixes=("_p", "_r"))
        m = m[m.step.isin(d.step)]
        norm = (m["norm_check_r"] if "norm_check_r" in m else m["norm_check"]).to_numpy()
        T1_ev = m["e_kin_ha"].to_numpy() * HA_TO_EV
        n = min(len(e_ev), len(norm))
        corr = e_ev - e_gs_ev
        corr[:n] = corr[:n] - T1_ev[:n] * (1.0 - norm[:n])
        out["dE_corr"] = corr
        out["norm"] = np.concatenate([norm[:n], np.full(len(e_ev) - n, np.nan)])
    else:
        out["dE_corr"] = out["dE_raw"]
        out["norm"] = np.nan
    return out
