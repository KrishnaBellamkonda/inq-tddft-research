"""
Stopping power for the slab->bulk L_slab sweep (`lz_bulk_sweep`).

THE MEASUREMENT (plan: docs/plans/jellium-slab-extend-Lz.md):

    S(L) = [E_total(t_final) - E_GS - E_PS(t_final)] / L_slab,
    extrapolated in 1/L_slab:  S(L) = S_bulk + c/L

with L_slab in {15, 25, 35} Bohr at fixed n0 (r_s = 4.183). The L = 25 points
are the EXISTING anchors (sigma56_sv twins for sigma = 5; the wp_highdensity_sv
sigma = 0.5 runs and the classical CAP sweep for sigma = 0.5) — this module
measures only the new L = 15/35 boxes and knows how to READ the anchors.

The arithmetic is NOT reimplemented: `slab_ks_wrap/e_absorbed.py::measure_dir`
does the segment concatenation, the WP norm correction and the plateau check
(validated against the wp_highdensity_sv notebooks to <= 3e-8). This module is
the ADAPTER: box presets, run naming, per-box E_GS, the per-box L_slab divisor
(measure_dir's own S columns assume 25 Bohr and are RESCALED here), the E_PS
monopole-tail correction, and the dispersion geometry.

GEOMETRY IS PER SIGMA FAMILY (user, 2026-08-05): each family's launch standoff
and face->CAP gap replicate its own L = 25 anchor, so arrival width depends on
(sigma, v) only. Cross-family comparisons are bracket-level.

The tables below are the python mirror of shared/bin/lzb-params.sh and
shared/configs/lzb_boxes.hpp — the three MUST agree (self-test: `python
lzb_stopping.py` recomputes the step formula and asserts).
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
LJ = REPO / "ResearchProject/systems/localised_jellium"
SCRIPTS = LJ / "scripts/lz_bulk_sweep"

# Reuse the validated engine rather than reimplementing the estimator.
sys.path.insert(0, str(LJ / "hypotheses/slab_ks_wrap"))
from e_absorbed import measure_dir, _concat, HA_TO_EV   # noqa: E402

DT = 0.04
DX = 0.40
LXY = 35.0
R_S = 4.183
N_BATH = {  # bath electron count per family (the monopole tail is N_e/z)
    0.5: {"L15": 60, "L25": 100, "L35": 140},
}

VELOCITIES = (2.0, 2.5, 3.0, 3.5)
PILOT_V = 3.0

# ---------------------------------------------------------------------------
# Box presets — MIRROR of lzb_boxes.hpp / lzb-params.sh. Literal, not derived.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Box:
    name: str
    LZ: float
    half: float
    n_e: int
    extra_states: int
    launch_z: float
    sigma: float
    gs_tag: str

    @property
    def l_slab(self) -> float: return 2.0 * self.half
    @property
    def standoff(self) -> float: return -self.half - self.launch_z


CFGS: dict[str, Box] = {
    "s0p5_L15": Box("s0p5_L15",  75.0,  7.5,  60, 15, -19.0, 0.5, "slab_n60_L35x35x75_dx0p4_per2"),
    "s0p5_L35": Box("s0p5_L35",  95.0, 17.5, 140, 34, -29.0, 0.5, "slab_n140_L35x35x95_dx0p4_per2"),
    "s5p0_L15": Box("s5p0_L15",  95.0,  7.5,  60, 15, -22.5, 5.0, "slab_n60_L35x35x95_dx0p4_per2"),
    "s5p0_L35": Box("s5p0_L35", 115.0, 17.5, 140, 34, -32.5, 5.0, "slab_n140_L35x35x115_dx0p4_per2"),
}

# N_STEPS = round(4.36 * (|launch_z| + L_z/2) / (v * dt)) — literal (bash mirror).
STEPS_TARGET: dict[str, dict[float, int]] = {
    "s0p5_L15": {2.0: 3079, 2.5: 2463, 3.0: 2053, 3.5: 1760},
    "s0p5_L35": {2.0: 4169, 2.5: 3335, 3.0: 2780, 3.5: 2382},
    "s5p0_L15": {2.0: 3815, 2.5: 3052, 3.0: 2543, 3.5: 2180},
    "s5p0_L35": {2.0: 4905, 2.5: 3924, 3.0: 3270, 3.5: 2803},
}


def steps_formula(cfg: str, v: float) -> int:
    b = CFGS[cfg]
    return int(round(4.36 * (abs(b.launch_z) + b.LZ / 2.0) / (v * DT)))


def v_tag(v: float) -> str:
    return "v" + f"{v:.1f}".replace(".", "p")


def run_name(cfg: str, v: float, half: str = "wp") -> str:
    """wp -> s0p5_L15_v3p0; classical -> cl_...; vac -> vac_... (dispatch names)."""
    base = f"{cfg}_{v_tag(v)}"
    return {"wp": base, "classical": "cl_" + base, "vac": "vac_" + base}[half]


def run_dir(cfg: str, v: float, half: str = "wp") -> Path:
    sub = {"wp": "wp", "classical": "classical", "vac": "vac"}[half]
    return SCRIPTS / sub / "results" / run_name(cfg, v, half)


def summary_kv(path: Path) -> dict[str, str]:
    """run_summary.txt key = value pairs (several may share a line)."""
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for m in re.finditer(r"(\w+)\s*=\s*([^\s].*?)(?=\s+\w+\s*=|\s*$)",
                        path.read_text(), re.M):
        out.setdefault(m.group(1), m.group(2).strip())
    return out


def e_gs_ha(cfg: str) -> float:
    """Per-box E_GS, read from the GS run — never hard-coded, never shared
    across boxes (each thickness is a different electronic system)."""
    summary = SCRIPTS / "gs" / "results" / cfg / "run_summary.txt"
    if not summary.exists():
        raise FileNotFoundError(
            f"no ground-state summary at {summary} — run-lzb-gs.slurm {cfg} first")
    kv = summary_kv(summary)
    if "ground_state_energy_ha" not in kv:
        raise KeyError(f"ground_state_energy_ha not found in {summary}")
    return float(kv["ground_state_energy_ha"])


# ---------------------------------------------------------------------------
# Dispersion geometry (per-family launch)
# ---------------------------------------------------------------------------
def sigma_d(t, sigma: float):
    return np.sqrt(sigma**2 / 2.0 + np.asarray(t, float)**2 / (2.0 * sigma**2))


def transit_window(cfg: str, v: float) -> tuple[float, float]:
    b = CFGS[cfg]
    return ((abs(b.launch_z) - b.half) / v, (abs(b.launch_z) + b.half) / v)


def transverse_overlap_time(sigma: float, l_xy: float = LXY) -> float:
    target = 2.0 * (l_xy / 6.0) ** 2
    return float(sigma * np.sqrt(max(target - sigma**2, 0.0)))


def mean_sigma_d(cfg: str, v: float, n: int = 4001) -> float:
    b = CFGS[cfg]
    ti, to = transit_window(cfg, v)
    t = np.linspace(ti, to, n)
    return float(np.trapezoid(sigma_d(t, b.sigma), t) / (to - ti))


def sigma_eq(cfg: str, v: float) -> float:
    return float(np.sqrt(2.0) * mean_sigma_d(cfg, v))


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Point:
    cfg: str
    sigma_wp: float
    L_slab: float
    inv_L: float
    v: float
    half: str
    run: str
    S_eV_per_Bohr: float           # norm-corrected deposit / L_slab (no E_PS cut)
    S_raw_eV_per_Bohr: float
    E_absorbed_eV: float
    t_final_au: float
    steps_done: int
    steps_target: int
    complete: bool
    settled: bool
    plateau_drift_eV: float
    norm_final: float
    sigma_d_entry: float
    sigma_d_exit: float
    sigma_eq: float
    e_ps_final_eV: float           # projectile-bath interaction still in the ledger
    z_proj_final: float
    S_deposit_eV_per_Bohr: float   # E_PS-corrected — USE THIS ONE


def e_ps_final(cfg: str, v: float, half: str) -> float:
    """E_PS at the last recorded step, in eV. NaN if unreadable (incomplete runs
    are filtered on `complete` downstream)."""
    obs = run_dir(cfg, v, half) / "raw" / "observables"
    try:
        return float(_concat(obs, "interactions")["e_ps"].to_numpy()[-1]) * HA_TO_EV
    except (FileNotFoundError, KeyError, IndexError, ValueError,
            pd.errors.EmptyDataError):
        return float("nan")


def measure(cfg: str, v: float, half: str = "wp") -> Point:
    """One (box, v, half) point of the new runs, evidence attached.

    measure_dir's own S columns divide by the 25-Bohr L_SLAB_Z constant of the
    legacy campaigns; every S here is RE-derived from E_absorbed_eV with the
    box's own L_slab — the whole point of this sweep.
    """
    b = CFGS[cfg]
    d = run_dir(cfg, v, half)
    if not d.exists():
        raise FileNotFoundError(d)
    a = measure_dir(d, e_gs_ha(cfg), half, cfg, v, R_S)
    target = STEPS_TARGET[cfg][v]
    ti, to = transit_window(cfg, v)

    # E_PS monopole-tail correction (sigma56 finding, 2026-08-03): the classical
    # projectile keeps travelling and its N_e/z tail never decays within any
    # affordable run; the WP is annihilated by the CAP so its E_PS(t_f) ~ 0.
    # Subtracting E_PS(t_f) makes both halves measure the same quantity.
    eps_f = e_ps_final(cfg, v, half)
    kv = summary_kv(d / "run_summary.txt")
    z_f = float(kv.get("proj_z_final", "nan")) if half == "classical" else float("nan")

    return Point(
        cfg=cfg, sigma_wp=b.sigma, L_slab=b.l_slab, inv_L=1.0 / b.l_slab,
        v=v, half=half, run=d.name,
        S_eV_per_Bohr=a.E_absorbed_eV / b.l_slab,
        S_raw_eV_per_Bohr=a.E_absorbed_raw_eV / b.l_slab,
        E_absorbed_eV=a.E_absorbed_eV,
        t_final_au=a.t_final_au,
        steps_done=a.steps_done, steps_target=target,
        complete=(a.steps_done >= target),
        settled=a.settled, plateau_drift_eV=a.plateau_drift_eV,
        norm_final=a.norm_final,
        sigma_d_entry=float(sigma_d(ti, b.sigma)),
        sigma_d_exit=float(sigma_d(to, b.sigma)),
        sigma_eq=sigma_eq(cfg, v),
        e_ps_final_eV=eps_f,
        z_proj_final=z_f,
        S_deposit_eV_per_Bohr=(a.E_absorbed_eV - (0.0 if np.isnan(eps_f) else eps_f))
                              / b.l_slab,
    )


def table() -> pd.DataFrame:
    """Every NEW point of the sweep (anchors are separate — `anchors()`).
    Missing/unfinished runs are reported and kept with complete=False."""
    rows = []
    for cfg in CFGS:
        for v in VELOCITIES:
            for half in ("wp", "classical"):
                try:
                    rows.append(asdict(measure(cfg, v, half)))
                except (FileNotFoundError, KeyError) as e:
                    print(f"  MISSING {run_name(cfg, v, half)}: {type(e).__name__}")
    if not rows:
        return pd.DataFrame()
    return (pd.DataFrame(rows)
            .sort_values(["sigma_wp", "half", "L_slab", "v"]).reset_index(drop=True))


# ---------------------------------------------------------------------------
# The L = 25 anchors (existing campaigns, used as-is — user decision 2026-08-05)
# ---------------------------------------------------------------------------
def anchors() -> pd.DataFrame:
    """L_slab = 25 points from the completed campaigns, on the SAME corrected-
    deposit estimator. Columns: sigma_wp, v, half, L_slab, inv_L, S, source.
    Anything unreadable is skipped WITH A PRINT — no silent gaps."""
    rows: list[dict] = []

    # sigma = 5: the sigma56_sv twins (E_PS-corrected S_deposit column).
    s56 = LJ / "hypotheses/sigma56_sv/s56_S_summary.csv"
    try:
        t = pd.read_csv(s56)
        t = t[(t.sigma_wp == 5.0) & t.complete & t.cap]
        for r in t.itertuples():
            rows.append({"sigma_wp": 5.0, "v": r.v, "half": r.half, "L_slab": 25.0,
                         "S": r.S_deposit_eV_per_Bohr, "source": "sigma56_sv"})
    except Exception as e:                                        # noqa: BLE001
        print(f"  ANCHOR SKIPPED sigma=5 ({s56.name}): {type(e).__name__}: {e}")

    # sigma = 0.5 WP: the wp_highdensity_sv deposit sweep (norm-corrected).
    # NOTE: no E_PS correction needed — the CAP annihilates the packet.
    swp = LJ / "hypotheses/wp_highdensity_sv/sigma_sweep_S_deposit.csv"
    try:
        t = pd.read_csv(swp)
        t = t[(t.sigma == 0.5) & t.complete]
        for r in t.itertuples():
            rows.append({"sigma_wp": 0.5, "v": r.v, "half": "wp", "L_slab": 25.0,
                         "S": r.S_deposit_corrected, "source": "wp_highdensity_sv"})
    except Exception as e:                                        # noqa: BLE001
        print(f"  ANCHOR SKIPPED sigma=0.5 wp ({swp.name}): {type(e).__name__}: {e}")

    # sigma = 0.5 classical: the merged CAP sweep. Its S_B_Eabs carries the
    # monopole tail, so re-derive: S = (E_absorbed - N_e/z_final) / 25 with
    # N_e = 100. Verified against the recorded 0.760 at v = 2.0 (and against
    # S_A_keloss to 0.3-5 %, an estimator sharing no machinery).
    scl = (LJ / "hypotheses/classical_highdensity_sv/dyn_direct/S_of_v_cap.csv")
    try:
        t = pd.read_csv(scl)
        for r in t.itertuples():
            monopole_ev = 100.0 / r.z_final * HA_TO_EV
            rows.append({"sigma_wp": 0.5, "v": r.v, "half": "classical",
                         "L_slab": 25.0,
                         "S": (r.E_absorbed_eV - monopole_ev) / 25.0,
                         "source": "cl_hd_sv_cap(EPS-corr)"})
    except Exception as e:                                        # noqa: BLE001
        print(f"  ANCHOR SKIPPED sigma=0.5 classical ({scl.name}): {type(e).__name__}: {e}")

    df = pd.DataFrame(rows)
    if not df.empty:
        df["inv_L"] = 1.0 / df["L_slab"]
    return df


def anchor_S(sigma: float, v: float, half: str) -> float:
    """One anchor value, NaN when absent — for INFO comparisons in the gate."""
    a = anchors()
    if a.empty:
        return float("nan")
    m = a[(a.sigma_wp == sigma) & (a.v == v) & (a.half == half)]
    return float(m["S"].iloc[0]) if len(m) else float("nan")


if __name__ == "__main__":
    print("lz_bulk_sweep — box presets and self-test\n")
    bad = 0
    for cfg, b in CFGS.items():
        n0 = b.n_e / (LXY * LXY * b.l_slab)
        rs = (3.0 / (4.0 * np.pi * n0)) ** (1.0 / 3.0)
        print(f"  {cfg}: L_z={b.LZ} L_slab={b.l_slab} N={b.n_e} launch={b.launch_z} "
              f"standoff={b.standoff:.1f} sigma={b.sigma} r_s={rs:.4f}")
        assert abs(rs - R_S) < 5e-3, f"{cfg}: r_s mismatch"
        for v in VELOCITIES:
            lit, form = STEPS_TARGET[cfg][v], steps_formula(cfg, v)
            if lit != form:
                print(f"    STEP TABLE MISMATCH at v={v}: literal {lit} != formula {form}")
                bad += 1
    print("\ndispersion geometry:")
    print(f"  {'cfg':>9} {'v':>4} {'t_in':>6} {'t_out':>6} {'sd_in':>7} {'sd_out':>7} {'sig_eq':>7} {'t_ov':>6}")
    for cfg, b in CFGS.items():
        for v in VELOCITIES:
            ti, to = transit_window(cfg, v)
            print(f"  {cfg:>9} {v:>4.1f} {ti:>6.2f} {to:>6.2f} {sigma_d(ti, b.sigma):>7.2f} "
                  f"{sigma_d(to, b.sigma):>7.2f} {sigma_eq(cfg, v):>7.2f} "
                  f"{transverse_overlap_time(b.sigma):>6.1f}")
    print("\nanchors:")
    a = anchors()
    print(a.to_string(index=False) if not a.empty else "  (none readable)")
    t = table()
    if not t.empty:
        print("\nnew points:")
        print(t.to_string(index=False))
    raise SystemExit(1 if bad else 0)
