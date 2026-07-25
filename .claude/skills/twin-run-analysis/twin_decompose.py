#!/usr/bin/env python3
"""twin_decompose — deterministic energy-decomposition engine for a twin run PAIR.

A *twin pair* is two localised-jellium runs identical in every physical parameter
EXCEPT the projectile representation:
  - ``wp/``        projectile = wavepacket (a real KS electron of width sigma_WP)
  - ``classical/`` projectile = Gaussian-charge perturbation (or ghost UPF)

This module does ONLY the arithmetic that has a known formula (the "deterministic
Python" half of the twin-run-analysis skill). It parses both runs, asserts config
parity, and at every timestep computes the per-term classical-vs-WP difference,
the residual ``d(E_H+E_ext) - U_proj_bg`` (= the WP self-Hartree), and the known
attributions (localisation kinetic, WP self-Hartree, one-electron SIE). It emits a
structured *findings table*; the physical narrative is the agent's job (see SKILL.md).

Native INQ energies are in Hartree; everything reported here is converted to eV.

Sign / convention notes (see docs/handovers/localised-jellium-energy-book-keeping.md):
  d(X)          = X_wp - X_classical                (WP minus classical)
  residual R    = d(E_H + E_ext) - U_proj_bg        = WP self-Hartree  E_H[WP-WP]
  SIE           = R + dXC                            = LDA one-electron self-interaction error
  dKin(at rest) = 3/(4 sigma_WP^2)                   = WP localisation zero-point   (+ k0^2/2 if k0!=0)

Only numpy + pandas — self-contained and shippable with the skill.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np
import pandas as pd

HA_EV = 27.211386          # 1 Hartree in eV
OBS_REL = "raw/observables/observables.csv"
SUMMARY = "run_summary.txt"

# Fields that MUST agree between the two twins (a real difference => not a valid twin).
PARITY_FIELDS = ("periodicity", "lz", "spacing", "n", "sigma_wp", "launch_z", "gs_dir")
# Absolute tolerance for float parity (spacing/sigma/launch_z are set to few decimals).
PARITY_ATOL = 1e-4


# ----------------------------------------------------------------------------- parsing
def parse_summary(path: str | Path) -> dict:
    """Parse a run_summary.txt into a flat lowercased dict.

    Robust to INQ's free-form summaries: multiple ``key = value`` tokens per line,
    parenthetical asides (``launch_z = -24.5  (r_from_face = 12 Bohr)``), and the
    ``key : value`` variant. Later occurrences win. Values are kept as strings;
    use :func:`_num` to coerce.
    """
    text = Path(path).read_text()
    out: dict[str, str] = {}
    for key, val in re.findall(r"([A-Za-z_]\w*)\s*[=:]\s*(\S+)", text):
        out[key.lower()] = val
    return out


def _num(d: dict, *keys, default=None):
    """First key in *keys* present in *d*, coerced to float; else *default*."""
    for k in keys:
        if k in d:
            try:
                return float(d[k])
            except ValueError:
                return d[k]
    return default


def infer_representation(d: dict) -> str | None:
    """Classify the classical projectile from free-form summary text.

    'perturbation' (Gaussian charge) | 'pseudopotential' (ghost UPF) | 'wavepacket'
    | None. An explicit ``representation = ...`` key wins; else infer from the
    projectile/run/mode strings.
    """
    if "representation" in d:
        return d["representation"].lower()
    blob = " ".join(str(v) for v in d.values()).lower()
    if "wavepacket" in blob or "wp" in (d.get("mode", "").lower()):
        if "perturbation" not in blob and "ghost" not in blob:
            return "wavepacket"
    if "perturbation" in blob:
        return "perturbation"
    if "ghost" in blob or "pseudopotential" in blob or ".upf" in blob or "z_valence" in blob:
        return "pseudopotential"
    return None


@dataclass
class RunConfig:
    """Normalised, engine-facing view of one run's provenance."""
    periodicity: float | None = None
    lz: float | None = None
    spacing: float | None = None
    n: float | None = None
    sigma_wp: float | None = None
    launch_z: float | None = None
    k0: float = 0.0
    gs_dir: str | None = None
    u_proj_bg_ev: float | None = None
    representation: str | None = None
    raw: dict = field(default_factory=dict)

    @classmethod
    def from_summary(cls, path: str | Path) -> "RunConfig":
        d = parse_summary(path)
        u_ev = _num(d, "u_proj_bg_ev")
        if u_ev is None and "u_proj_bg_ha" in d:
            u_ev = float(d["u_proj_bg_ha"]) * HA_EV
        return cls(
            periodicity=_num(d, "periodicity"),
            lz=_num(d, "lz"),
            spacing=_num(d, "spacing"),
            n=_num(d, "n"),
            sigma_wp=_num(d, "sigma_wp"),         # sigma_wp / sigma_WP both lowercased
            launch_z=_num(d, "launch_z"),
            k0=_num(d, "k0", default=0.0) or 0.0,
            gs_dir=d.get("gs_dir"),
            u_proj_bg_ev=u_ev,
            representation=infer_representation(d),
            raw=d,
        )


# Optional per-step auxiliary CSVs a run emits alongside observables.csv
# (kept run-local instead of bloating the shared ObservablesWriter schema).
#  - projectile.csv   : trajectory + proj KE + per-step U_proj_bg (dynamic classical)
#  - wp_centroid.csv   : WP centroid / spread
#  - interactions.csv  : the pairwise P/S/B Coulomb decomposition (both runs)
AUX_CSVS = ("projectile.csv", "wp_centroid.csv", "interactions.csv")

# The pairwise Coulomb terms (Hartree in the CSV). P=projectile, S=slab, B=background.
PAIRWISE_TERMS = ("e_ss", "e_pp", "e_ps", "e_sb", "e_pb", "e_bb")
# Physically-identical between the twins (slab & background unchanged) → Δ must be
# ~0 if there is no inter-run gauge (reference_twin_pairwise_decomposition).
GAUGE_INVARIANT_TERMS = ("e_ss", "e_sb", "e_bb")


def load_run(run_dir: str | Path) -> tuple[pd.DataFrame, RunConfig]:
    """Load one run directory -> (observables DataFrame in eV, RunConfig).

    Merges any auxiliary dynamic CSVs (projectile.csv, wp_centroid.csv) by step.
    Every energy_* column (native Hartree) gets an eV twin `<col>_ev`.
    """
    run_dir = Path(run_dir)
    obs_path = run_dir / OBS_REL
    obs = pd.read_csv(obs_path)
    cfg = RunConfig.from_summary(run_dir / SUMMARY)

    for aux in AUX_CSVS:                                   # merge dynamic trajectory/KE
        p = obs_path.parent / aux
        if not p.exists():
            continue
        adf = pd.read_csv(p)
        new = [c for c in adf.columns if c not in obs.columns and c != "time_au"]
        obs = obs.merge(adf[["step"] + new], on="step", how="left")

    for col in [c for c in obs.columns if c.startswith("energy_")]:
        obs[col + "_ev"] = obs[col] * HA_EV
    return obs, cfg


# ----------------------------------------------------------------------------- parity
@dataclass
class ParityReport:
    ok: bool
    mismatches: list[dict]
    checked: list[str]

    def as_text(self) -> str:
        if self.ok:
            return "config parity: OK (%s)" % ", ".join(self.checked)
        lines = ["config parity: FAIL"]
        for m in self.mismatches:
            lines.append(f"  {m['field']}: wp={m['wp']!r} classical={m['classical']!r}")
        return "\n".join(lines)


def check_parity(wp: RunConfig, cl: RunConfig, fields=PARITY_FIELDS) -> ParityReport:
    """Assert the two twins share every physical parameter (projectile excepted)."""
    mism, checked = [], []
    for f in fields:
        a, b = getattr(wp, f), getattr(cl, f)
        if a is None or b is None:      # not recorded in both summaries -> skip
            continue
        checked.append(f)
        same = (
            math.isclose(a, b, abs_tol=PARITY_ATOL)
            if isinstance(a, (int, float)) and isinstance(b, (int, float))
            else a == b
        )
        if not same:
            mism.append({"field": f, "wp": a, "classical": b})
    return ParityReport(ok=not mism, mismatches=mism, checked=checked)


# ----------------------------------------------------------------------------- physics
def loc_kinetic_ev(sigma_wp: float, k0: float = 0.0) -> float:
    """Expected WP kinetic surplus = localisation zero-point 3/(4 sigma^2) + k0^2/2 (Ha->eV)."""
    return (3.0 / (4.0 * sigma_wp**2) + 0.5 * k0**2) * HA_EV


def self_hartree_freespace_ev(sigma_wp: float) -> float:
    """Free-space Gaussian self-Hartree 1/(sigma_WP*sqrt(2*pi)) (Ha->eV).

    Reference only: the boundary-matched (open-z) value is ~0.9 eV lower; the
    empirically-correct target is the measured residual itself. See the handover.
    """
    return (1.0 / (sigma_wp * math.sqrt(2.0 * math.pi))) * HA_EV


# ----------------------------------------------------------------------------- engine
def _u_proj_bg_series(cl_obs: pd.DataFrame, cl_cfg: RunConfig) -> pd.Series:
    """Per-step U_proj_bg. Prefer a per-step column (dynamic runs); else the
    constant from run_summary broadcast over all steps (at-rest runs)."""
    if "energy_proj_bg_ideal_ev" in cl_obs.columns:
        return cl_obs["energy_proj_bg_ideal_ev"].reset_index(drop=True)
    if cl_cfg.u_proj_bg_ev is None:
        raise ValueError(
            "U_proj_bg unavailable: no energy_proj_bg_ideal column and no "
            "U_proj_bg in the classical run_summary.txt")
    return pd.Series([cl_cfg.u_proj_bg_ev] * len(cl_obs))


def _opt(obs: pd.DataFrame, col: str, n: int):
    """A per-step column if present (eV auto-suffixed for energies), else None."""
    name = col + "_ev" if (col + "_ev") in obs.columns else col
    return obs[name].to_numpy()[:n] if name in obs.columns else None


def decompose(wp_dir: str | Path, classical_dir: str | Path,
              drift_flag_ev: float = 0.5) -> "TwinResult":
    """Full deterministic decomposition of a twin pair (static or dynamic).

    Returns a :class:`TwinResult` with the per-step table, the attribution
    (findings) table, the parity report, and — when the Rung-2 columns are present
    — energy-conservation and trajectory tracking. Defensive to absent dynamic
    columns, so static runs behave exactly as Rung 1. *drift_flag_ev* is the
    per-step change above which a term is flagged as "moving" (onset of dynamics).
    """
    wp_obs, wp_cfg = load_run(wp_dir)
    cl_obs, cl_cfg = load_run(classical_dir)
    parity = check_parity(wp_cfg, cl_cfg)
    representation = cl_cfg.representation or "perturbation"

    n = min(len(wp_obs), len(cl_obs))
    wp_obs, cl_obs = wp_obs.iloc[:n].reset_index(drop=True), cl_obs.iloc[:n].reset_index(drop=True)

    def d(col):     # WP - classical, in eV, per step
        return wp_obs[col + "_ev"].to_numpy() - cl_obs[col + "_ev"].to_numpy()

    dKin, dH, dXC, dExt = d("energy_kinetic"), d("energy_hartree"), d("energy_xc"), d("energy_external")
    d_H_ext = dH + dExt
    u_proj = _u_proj_bg_series(cl_obs, cl_cfg).to_numpy()[:n]
    # U_proj_bg sign is representation-dependent: for a z_valence=0 GHOST UPF, INQ's
    # E_external OMITS the projectile<->background compensation term (the h0 note's
    # "re-add int v_ghost*n_+"), so U_proj_bg must be ADDED, not subtracted.
    u_sign = -1.0 if representation == "pseudopotential" else +1.0
    residual = d_H_ext - u_sign * u_proj      # = WP self-Hartree
    sie = residual + dXC                       # = LDA one-electron SIE (clean only for perturbation)

    cols = {
        "step": wp_obs["step"].to_numpy(), "time_au": wp_obs["time_au"].to_numpy(),
        "dKin": dKin, "dHartree": dH, "dXC": dXC, "dExt": dExt,
        "d_H_ext": d_H_ext, "U_proj_bg": u_proj, "residual": residual, "sie": sie,
    }

    # --- Rung-2 dynamic columns (all optional) --------------------------------
    # Classical projectile KE (perturbation: energy_proj_ke; ghost ion: energy_ion_kinetic).
    proj_ke = _opt(cl_obs, "energy_proj_ke", n)
    if proj_ke is None:
        proj_ke = _opt(cl_obs, "energy_ion_kinetic", n)
    proj_ke = np.zeros(n) if proj_ke is None else proj_ke
    cols["proj_ke_classical"] = proj_ke
    # Motional-matched localisation: subtract the classical projectile KE so only
    # the WP zero-point spread remains (== dKin when at rest / no proj_ke).
    cols["dKin_localisation"] = dKin - proj_ke

    # Trajectory (classical proj_z vs WP centroid) and WP spread.
    proj_z = _opt(cl_obs, "proj_z", n)
    wp_z = _opt(wp_obs, "wp_centroid_z", n)
    wp_sig = _opt(wp_obs, "wp_sigma_z", n)
    if proj_z is not None:
        cols["proj_z"] = proj_z
    if wp_z is not None:
        cols["wp_centroid_z"] = wp_z
    if wp_sig is not None:
        cols["wp_sigma_z"] = wp_sig
    if proj_z is not None and wp_z is not None:
        cols["separation_z"] = wp_z - proj_z

    # Energy conservation: E_electronic (+ classical proj KE + U_proj_bg) must be flat.
    def _total(obs):
        if "energy_total_ev" in obs.columns:
            return obs["energy_total_ev"].to_numpy()[:n]
        parts = [c for c in ("energy_kinetic_ev", "energy_hartree_ev", "energy_xc_ev",
                             "energy_external_ev", "energy_nonlocal_ev", "energy_ion_ev")
                 if c in obs.columns]
        return obs[parts].sum(axis=1).to_numpy()[:n]
    cl_tot, wp_tot = _total(cl_obs), _total(wp_obs)
    cols["E_conserved_classical"] = cl_tot + proj_ke + u_proj
    cols["E_conserved_wp"] = wp_tot
    # Quantum stopping proxy = total electronic energy DEPOSITED (NOT projectile KE!):
    # rise of the WP-run electronic total above its t=0 value.
    cols["E_deposited_wp"] = wp_tot - wp_tot[0]

    steps = pd.DataFrame(cols)
    sig = wp_cfg.sigma_wp or cl_cfg.sigma_wp
    k0 = wp_cfg.k0 or 0.0
    r0 = steps.iloc[0]

    # Representation-aware residual interpretation.
    sh = self_hartree_freespace_ev(sig) if sig else None
    if representation == "pseudopotential":
        res_label = "residual R = d(E_H+E_ext) + U_proj_bg"   # ghost: ADD (omitted comp. term)
        res_note = ("WP self-Hartree; for a GHOST-UPF projectile U_proj_bg is ADDED "
                    "(INQ omits the z_valence=0 background-compensation term). The shortfall "
                    "vs the ~21 eV self-Hartree (and vs the clean perturbation ~20.8 eV) is "
                    "the KNOWN ghost-UPF tail-aliasing (~12-14 eV), NOT missing physics")
        sie_note = ("R + dXC; NOT a clean SIE for the pseudopotential representation — R is "
                    "aliasing-corrupted, so this is not the physical one-electron SIE "
                    "(use the perturbation twin for the clean 4.34 eV)")
    else:
        res_label = "residual R = d(E_H+E_ext) - U_proj_bg"
        res_note = ("WP self-Hartree E_H[WP-WP]; the ~0.9 eV shortfall vs the free-space "
                    "ref is the open-z gauge, not missing physics")
        sie_note = ("LDA one-electron self-interaction error (Perdew-Zunger): "
                    "the part of the WP self-Hartree that XC fails to cancel")

    findings = [
        _finding("dKin_localisation (motional-matched)", r0.dKin_localisation,
                 loc_kinetic_ev(sig, 0.0) if sig else None,
                 "WP localisation zero-point 3/(4*sigma^2); classical proj KE subtracted"),
        _finding("dXC (xc surplus)", r0.dXC, None,
                 "additional XC of the WP alone; local => r-independent"),
        _finding(res_label, r0.residual, sh, res_note),
        _finding("SIE = R + dXC", r0.sie, None, sie_note),
    ]

    # Drift flags + energy-conservation gate.
    drift = {}
    for col in ("dKin", "dHartree", "dXC", "dExt", "residual"):
        series = steps[col].to_numpy()
        dmax = float(np.max(np.abs(np.diff(series)))) if len(series) > 1 else 0.0
        drift[col] = {"max_step_change_ev": dmax, "moving": dmax > drift_flag_ev}
    conservation = {}
    for tag, key in (("classical", "E_conserved_classical"), ("wp", "E_conserved_wp")):
        s = steps[key].to_numpy()
        conservation[tag] = float(np.max(np.abs(s - s[0]))) if len(s) > 1 else 0.0

    # Pairwise Coulomb decomposition (present when BOTH runs emit interactions.csv).
    # Every difference is physically attributable; the gauge-invariant terms (slab &
    # background, physically identical between the twins) must have Δ≈0 → no gauge.
    pairwise, gauge = None, None
    if all(t in wp_obs.columns for t in PAIRWISE_TERMS) and all(t in cl_obs.columns for t in PAIRWISE_TERMS):
        pw = {"step": steps.step.to_numpy(), "time_au": steps.time_au.to_numpy()}
        for t in PAIRWISE_TERMS:
            c = cl_obs[t].to_numpy()[:n] * HA_EV
            w = wp_obs[t].to_numpy()[:n] * HA_EV
            pw[t + "_cl"], pw[t + "_wp"], pw["d_" + t] = c, w, w - c
        pairwise = pd.DataFrame(pw)
        g0 = pairwise.iloc[0]
        gauge = {t: float(g0["d_" + t]) for t in GAUGE_INVARIANT_TERMS}
        gauge["max_invariant_delta_ev"] = max(abs(v) for v in gauge.values())
        gauge["no_gauge"] = gauge["max_invariant_delta_ev"] < 1e-2

    return TwinResult(steps=steps, findings=findings, parity=parity,
                      sigma_wp=sig, k0=k0, drift=drift, conservation=conservation,
                      representation=representation, pairwise=pairwise, gauge=gauge,
                      wp_dir=str(wp_dir), classical_dir=str(classical_dir))


def _finding(name, value, expected, interpretation):
    rem = None if expected is None else float(value - expected)
    return {
        "term": name,
        "value_ev": float(value),
        "expected_ev": None if expected is None else float(expected),
        "unexplained_ev": rem,
        "interpretation": interpretation,
    }


@dataclass
class TwinResult:
    steps: pd.DataFrame
    findings: list[dict]
    parity: ParityReport
    sigma_wp: float | None
    k0: float
    drift: dict
    conservation: dict
    representation: str
    wp_dir: str
    classical_dir: str
    pairwise: "pd.DataFrame | None" = None   # per-step P/S/B pairwise terms (cl/wp/Δ), eV
    gauge: dict | None = None                # gauge test on the physically-identical terms

    def findings_table(self) -> pd.DataFrame:
        return pd.DataFrame(self.findings)

    def pairwise_table(self, step: int = 0) -> "pd.DataFrame | None":
        """classical vs WP vs Δ for each pairwise Coulomb term at a given step (eV)."""
        if self.pairwise is None:
            return None
        r = self.pairwise[self.pairwise.step == step]
        if r.empty:
            r = self.pairwise.iloc[[0]]
        r = r.iloc[0]
        rows = [{"term": t, "classical": round(r[t + "_cl"], 3), "wavepacket": round(r[t + "_wp"], 3),
                 "delta_wp_minus_cl": round(r["d_" + t], 4)} for t in PAIRWISE_TERMS]
        return pd.DataFrame(rows)

    @property
    def is_dynamic(self) -> bool:
        return "separation_z" in self.steps.columns or bool(self.steps["proj_ke_classical"].any())

    def to_dict(self) -> dict:
        return {
            "wp_dir": self.wp_dir, "classical_dir": self.classical_dir,
            "representation": self.representation,
            "sigma_wp": self.sigma_wp, "k0": self.k0,
            "parity_ok": self.parity.ok, "parity_mismatches": self.parity.mismatches,
            "findings": self.findings, "drift": self.drift,
            "conservation": self.conservation, "is_dynamic": self.is_dynamic,
            "gauge": self.gauge,
            "pairwise": None if self.pairwise is None else self.pairwise.to_dict(orient="list"),
            "steps": self.steps.to_dict(orient="list"),
        }

    def report(self) -> str:
        """Human-scannable summary (rounded); the agent narrates from this."""
        lines = [
            f"Twin decomposition  representation={self.representation}  "
            f"sigma_WP={self.sigma_wp}  k0={self.k0}  "
            f"{'DYNAMIC' if self.is_dynamic else 'static'}",
            self.parity.as_text(),
            "",
            "Findings (step 0, eV):",
            f"  {'term':44s} {'value':>8s} {'expected':>9s} {'unexpl.':>8s}",
        ]
        for f in self.findings:
            exp = "-" if f["expected_ev"] is None else f"{f['expected_ev']:9.2f}"
            rem = "-" if f["unexplained_ev"] is None else f"{f['unexplained_ev']:8.2f}"
            lines.append(f"  {f['term']:44s} {f['value_ev']:8.2f} {exp} {rem}")
            lines.append(f"      -> {f['interpretation']}")
        lines.append("")
        lines.append("Per-step drift (max |Δ| between steps, eV):")
        for k, v in self.drift.items():
            tag = "  MOVING" if v["moving"] else ""
            lines.append(f"  {k:12s} {v['max_step_change_ev']:8.4f}{tag}")
        if self.pairwise is not None:
            lines.append("")
            lines.append("Pairwise Coulomb decomposition (step 0, eV):")
            lines.append(f"  {'term':10s} {'classical':>12s} {'wavepacket':>12s} {'Δ(WP-cl)':>10s}")
            g0 = self.pairwise.iloc[0]
            for t in PAIRWISE_TERMS:
                lines.append(f"  {t:10s} {g0[t+'_cl']:12.3f} {g0[t+'_wp']:12.3f} {g0['d_'+t]:10.4f}")
            if self.gauge is not None:
                verdict = "NO GAUGE (zero-points agree)" if self.gauge["no_gauge"] else "GAUGE PRESENT"
                lines.append(f"  gauge test: max|Δ| of slab/bg terms = "
                             f"{self.gauge['max_invariant_delta_ev']:.4f} eV -> {verdict}")
                lines.append("  -> every non-zero Δ is physically attributable (projectile terms only)")
        lines.append("")
        lines.append("Energy conservation (max |E(t)-E(0)|, eV; dynamics correctness gate):")
        for k, v in self.conservation.items():
            lines.append(f"  {k:12s} {v:8.4f}")
        if self.is_dynamic:
            lines.append("\nNote: quantum stopping = total electronic energy deposited "
                         "(E_deposited_wp), NOT projectile KE. Projectile-KE stopping is "
                         "CLASSICAL-only (use stopping-power-extraction on proj_ke_classical).")
        return "\n".join(lines)


# ----------------------------------------------------------------------------- CLI
def _resolve_pair(args) -> tuple[str, str]:
    if args.pair:
        p = Path(args.pair)
        return str(p / "wp"), str(p / "classical")
    if not (args.wp and args.classical):
        raise SystemExit("give either <pair_dir> or --wp DIR --classical DIR")
    return args.wp, args.classical


def main(argv=None):
    ap = argparse.ArgumentParser(description="Twin-run energy decomposition (deterministic engine).")
    ap.add_argument("pair", nargs="?", help="pair dir containing wp/ and classical/")
    ap.add_argument("--wp", help="wavepacket run dir")
    ap.add_argument("--classical", help="classical/perturbation run dir")
    ap.add_argument("--json", help="write full result JSON to this path")
    ap.add_argument("--csv", help="write per-step table to this path")
    args = ap.parse_args(argv)

    wp_dir, cl_dir = _resolve_pair(args)
    res = decompose(wp_dir, cl_dir)
    print(res.report())
    if args.json:
        Path(args.json).write_text(json.dumps(res.to_dict(), indent=2))
        print("\nwrote", args.json)
    if args.csv:
        res.steps.to_csv(args.csv, index=False)
        print("wrote", args.csv)
    return res


if __name__ == "__main__":
    main()
