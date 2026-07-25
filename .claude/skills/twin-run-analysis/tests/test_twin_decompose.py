#!/usr/bin/env python3
"""Known-answer tests for the twin_decompose engine.

Two layers:
  1. A SYNTHETIC fixture built from the documented decomposition table — tests the
     arithmetic in isolation, no run data required (portable, always runs).
  2. The GOLDEN on-disk pair (localised_jellium proj_perturbation, sigma=0.5, r=12)
     — the real end-to-end known answer:
         dKin 81.74 | dXC -16.47 | residual 20.81 | SIE 4.34 eV
     Skipped automatically if that data is not present.

Run:  /local/data/public/skcb2/tddft/venv/bin/python3 -m pytest <this file> -q
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import twin_decompose as td  # noqa: E402

HA = td.HA_EV

# Documented golden numbers (docs/notes/gaussian-pertubation-for-classical-simul,
# perturbation_method study). eV.
GOLD = dict(dKin=81.74, dHartree=-119.02, dXC=-16.47, dExt=274.51,
            d_H_ext=155.49, U_proj_bg=134.69, residual=20.81, sie=4.34)

GOLD_PAIR = Path(
    "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
    "scripts/localised_jellium_dynamics/proj_perturbation")
GOLD_WP = GOLD_PAIR / "stress_scratch/s0p5_r12_lz120_p2/results/wp"
GOLD_CL = GOLD_PAIR / "results/proj_pert_dx0p5"


# --------------------------------------------------------------------- synthetic
def _write_run(dirpath: Path, energies_ha: dict, summary: str):
    """Write a minimal one-step run (observables.csv + run_summary.txt) in Ha."""
    obs = dirpath / "raw/observables"
    obs.mkdir(parents=True, exist_ok=True)
    cols = {"step": 0, "time_au": 0.0}
    cols.update(energies_ha)
    pd.DataFrame([cols]).to_csv(obs / "observables.csv", index=False)
    (dirpath / "run_summary.txt").write_text(summary)


@pytest.fixture
def synthetic_pair(tmp_path):
    """Build a twin pair whose step-0 differences equal the documented table."""
    # classical baseline (arbitrary but fixed), WP = classical + documented delta.
    base = dict(energy_kinetic=74.97, energy_hartree=-2194.99,
                energy_xc=-230.87, energy_external=3854.13)
    delta = dict(energy_kinetic=GOLD["dKin"], energy_hartree=GOLD["dHartree"],
                 energy_xc=GOLD["dXC"], energy_external=GOLD["dExt"])
    cl = {k: v / HA for k, v in base.items()}                      # eV -> Ha
    wp = {k: (base[k] + delta[k]) / HA for k in base}
    common = "periodicity = 2  Lz = 120  spacing = 0.5  launch_z = -24.5\n"
    _write_run(tmp_path / "classical", cl,
               common + f"sigma_wp = 0.5\nU_proj_bg_ev = {GOLD['U_proj_bg']}\ngs_dir = /gs\n")
    _write_run(tmp_path / "wp", wp,
               common + "sigma_WP = 0.5\nk0 = 0\ngs_dir = /gs\n")
    return tmp_path


def test_synthetic_reproduces_documented_table(synthetic_pair):
    res = td.decompose(synthetic_pair / "wp", synthetic_pair / "classical")
    r0 = res.steps.iloc[0]
    for k, want in GOLD.items():
        assert math.isclose(r0[k], want, abs_tol=0.02), f"{k}: {r0[k]:.3f} != {want}"


def test_synthetic_parity_ok(synthetic_pair):
    res = td.decompose(synthetic_pair / "wp", synthetic_pair / "classical")
    assert res.parity.ok, res.parity.as_text()


def test_parity_catches_mismatch(synthetic_pair):
    # Corrupt the WP launch_z -> parity must FAIL (not a valid twin).
    s = (synthetic_pair / "wp/run_summary.txt").read_text().replace("-24.5", "-30.0")
    (synthetic_pair / "wp/run_summary.txt").write_text(s)
    res = td.decompose(synthetic_pair / "wp", synthetic_pair / "classical")
    assert not res.parity.ok
    assert any(m["field"] == "launch_z" for m in res.parity.mismatches)


def test_localisation_formula():
    assert math.isclose(td.loc_kinetic_ev(0.5), 3 / (4 * 0.25) * HA, rel_tol=1e-9)
    assert math.isclose(td.loc_kinetic_ev(0.5), 81.634, abs_tol=0.01)


def test_sie_identity(synthetic_pair):
    res = td.decompose(synthetic_pair / "wp", synthetic_pair / "classical")
    r0 = res.steps.iloc[0]
    assert math.isclose(r0.sie, r0.residual + r0.dXC, abs_tol=1e-9)


# --------------------------------------------------------------------- golden on-disk
golden_present = GOLD_WP.exists() and GOLD_CL.exists()


@pytest.mark.skipif(not golden_present, reason="golden proj_perturbation pair not on disk")
def test_golden_pair_reproduces_known_answer():
    res = td.decompose(GOLD_WP, GOLD_CL)
    assert res.parity.ok, res.parity.as_text()
    r0 = res.steps.iloc[0]
    assert math.isclose(r0.dKin, GOLD["dKin"], abs_tol=0.1)
    assert math.isclose(r0.dXC, GOLD["dXC"], abs_tol=0.1)
    assert math.isclose(r0.residual, GOLD["residual"], abs_tol=0.1)
    assert math.isclose(r0.sie, GOLD["sie"], abs_tol=0.1)


@pytest.mark.skipif(not golden_present, reason="golden proj_perturbation pair not on disk")
def test_golden_at_rest_is_static():
    """At-rest twin: no difference term should drift across the 2-3 steps."""
    res = td.decompose(GOLD_WP, GOLD_CL)
    assert not any(v["moving"] for v in res.drift.values()), res.drift


# --------------------------------------------------------------------- dynamic (Rung 2)
def _write_dynamic_run(dirpath: Path, energies_ha_by_step: list[dict],
                       extra_by_step: list[dict], summary: str):
    obs = dirpath / "raw/observables"
    obs.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, (e, x) in enumerate(zip(energies_ha_by_step, extra_by_step)):
        row = {"step": i, "time_au": 0.01 * i}
        row.update(e); row.update(x)
        rows.append(row)
    pd.DataFrame(rows).to_csv(obs / "observables.csv", index=False)
    (dirpath / "run_summary.txt").write_text(summary)


@pytest.fixture
def dynamic_pair(tmp_path):
    """3-step dynamic twin: classical projectile marches at constant v (linear proj_z),
    energy strictly conserved (E_electronic + proj_ke + U_proj_bg = const), WP centroid
    advances too. Representation = pseudopotential (ghost) to exercise routing."""
    HA = td.HA_EV
    v, dt, m = 0.3, 0.01, 1.0
    steps = 3
    U = 134.69  # constant here for a simple conservation check
    base = dict(kin=74.97, hart=-2194.99, xc=-230.87, ext=3854.13)
    cl_e, cl_x, wp_e, wp_x = [], [], [], []
    for i in range(steps):
        ke = 0.5 * m * v**2 * HA                      # constant v -> constant KE
        # keep classical electronic total fixed so E_conserved is exactly flat
        cl_e.append({k: base[b] / HA for k, b in
                     [("energy_kinetic", "kin"), ("energy_hartree", "hart"),
                      ("energy_xc", "xc"), ("energy_external", "ext")]}
                    | {"energy_total": sum(base.values()) / HA})
        cl_x.append({"energy_proj_ke": ke / HA, "energy_proj_bg_ideal": U / HA,
                     "proj_z": -24.5 + v * dt * i, "proj_vz": v})
        wpb = {"energy_kinetic": (base["kin"] + GOLD["dKin"]) / HA,
               "energy_hartree": (base["hart"] + GOLD["dHartree"]) / HA,
               "energy_xc": (base["xc"] + GOLD["dXC"]) / HA,
               "energy_external": (base["ext"] + GOLD["dExt"]) / HA}
        wp_e.append(wpb | {"energy_total": sum(wpb.values())})
        wp_x.append({"wp_centroid_z": -24.5 + v * dt * i, "wp_sigma_z": 0.5})
    common = "periodicity = 2  Lz = 120  spacing = 0.5  launch_z = -24.5\n"
    _write_dynamic_run(tmp_path / "classical", cl_e, cl_x,
                       common + "sigma_wp = 0.5\nprojectile = ghost UPF z_valence 0\n"
                       "representation = pseudopotential\ngs_dir = /gs\nrun_completed = true\n")
    _write_dynamic_run(tmp_path / "wp", wp_e, wp_x,
                       common + "sigma_WP = 0.5\nk0 = 0.3\nprojectile = wavepacket\n"
                       "representation = wavepacket\ngs_dir = /gs\nrun_completed = true\n")
    return tmp_path


def test_dynamic_flags_and_conservation(dynamic_pair):
    res = td.decompose(dynamic_pair / "wp", dynamic_pair / "classical")
    assert res.is_dynamic
    assert res.representation == "pseudopotential"
    # constant-v, fixed electronic total => conserved energy is flat to machine tol
    assert res.conservation["classical"] < 1e-6, res.conservation
    # trajectory recovered: linear proj_z, separation computed
    s = res.steps
    assert "separation_z" in s.columns
    assert np.allclose(np.diff(s.proj_z), 0.3 * 0.01, atol=1e-9)
    # quantum stopping proxy present and is total deposition, not projectile KE
    assert "E_deposited_wp" in s.columns


def test_pseudopotential_residual_note(dynamic_pair):
    res = td.decompose(dynamic_pair / "wp", dynamic_pair / "classical")
    resid = next(f for f in res.findings if f["term"].startswith("residual"))
    assert "ghost-UPF tail-aliasing" in resid["interpretation"]


def test_motional_matched_localisation(dynamic_pair):
    """dKin_localisation must subtract the classical projectile KE."""
    res = td.decompose(dynamic_pair / "wp", dynamic_pair / "classical")
    r0 = res.steps.iloc[0]
    assert math.isclose(r0.dKin_localisation, r0.dKin - r0.proj_ke_classical, abs_tol=1e-9)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
