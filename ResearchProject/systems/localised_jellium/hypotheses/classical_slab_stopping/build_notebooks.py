#!/usr/bin/env python3
"""Build the per-run run-notebooks for the classical-slab-stopping campaign.

Thin driver over the skill-local ``run_notebook_builder`` (run-notebook skill):
it reconciles the proj_dyn run's on-disk schema to what the shared builder
expects, then calls ``build()``. Nothing run-specific is added to the skill.

Two phases share one geometry (twin of the WP run ``p5_wp_v1p3``):
  P1  p1_ehrenfest_v1p3   free Ehrenfest projectile (decelerates); S = initial drag
  P2  p2_constv_v1p3      prescribed constant-velocity replica; S = deposit / L_slab

Schema bridges (proj_dyn writer -> builder expectations), all non-destructive:
  1. density frames live in  frames/total/density_t*.vti  (physical-order inqkit
     VTIs); the builder/battery read  raw/vti/density_total  (GIF battery) and
     raw/vti/density_system  (z-t carpets + lead GIF).  -> symlink both to
     frames/total.
  2. detect_type() keys "classical" off electron_track.csv; this run writes
     projectile.csv.  -> materialise electron_track.csv (step,time_au,z,vz,
     ke_ion_ha) from projectile.csv so the classical transport battery lights up
     and N(t) is read from the density VTIs.
  3. the builder reads dt_au / cell_bohr; the summary writes dt / Lz.  -> append
     dt_au and cell_bohr aliases (same values) to run_summary.txt if absent.

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_notebooks.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/local/data/public/skcb2/tddft")
SKILL = ROOT / ".claude/skills/run-notebook"
HYP = ROOT / "ResearchProject/systems/localised_jellium/hypotheses/classical_slab_stopping"
RESULTS = ROOT / "ResearchProject/systems/localised_jellium/scripts/classical_slab_stopping/results"
# WP twin (σ_WP=0.5, k0=1.3) — differenced against each classical run to give the
# WP−classical energy-difference bar GIF (user request, slow pace). Same geometry.
WP_TWIN = ROOT / "ResearchProject/systems/localised_jellium/scripts/qsp_phase5/wp/results/p5_wp_v1p3"

sys.path.insert(0, str(SKILL))
import run_notebook_builder as B  # noqa: E402

# shared geometry (both phases) — twin of p5_wp_v1p3
GEOM = dict(rs=5.667, v0=1.3, launch_z=-23.75, proj_sigma=0.3536,
            l_slab=25.0, cap_inner=45.0)  # cap_inner=45 = box half-z (CAP-FREE run);
#                                           marks the box edge, NOT a fake CAP.

# per-phase measured S (eV/Bohr -> Ha/Bohr) overlaid on the analytical Lindhard S(v).
HA_EV = 27.211386
PHASES = [
    dict(name="p1_ehrenfest_v1p3", nb="p1_ehrenfest_v1p3.ipynb",
         measured_s=0.4926 / HA_EV,   # initial-drag -dKE/ds (light-projectile rule)
         measured_v=1.3),
    dict(name="p2_constv_v1p3", nb="p2_constv_v1p3.ipynb",
         measured_s=0.4315 / HA_EV,   # deposit / L_slab (const-v: KE is prescribed)
         measured_v=1.3),
]


def bridge(run_dir: Path):
    """Reconcile proj_dyn on-disk layout to the builder's expectations."""
    # (1) density frames  frames/total  ->  raw/vti/{density_total,density_system}
    src = run_dir / "frames" / "total"
    vti = run_dir / "raw" / "vti"
    vti.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        for link_name in ("density_total", "density_system"):
            link = vti / link_name
            if not link.exists():
                link.symlink_to(os.path.relpath(src, vti))

    # (2) electron_track.csv from projectile.csv (classical detection + transport).
    #     Add fz = mass·dvz/dt (mass=1) — the pipeline `stopping` phase requires an
    #     `fz` column; projectile.csv logs no force, so recover it from the velocity
    #     track by finite difference.
    obs = run_dir / "raw" / "observables"
    proj = pd.read_csv(obs / "projectile.csv").drop_duplicates("step")
    t = proj["time_au"].to_numpy()
    vz = proj["proj_vz"].to_numpy()
    fz = np.gradient(vz, t) if len(t) > 1 else np.zeros_like(vz)   # mass=1
    pd.DataFrame({
        "step": proj["step"], "time_au": t,
        "z": proj["proj_z"], "vz": vz,
        "ke_ion_ha": proj["energy_proj_ke"], "fz": fz,
    }).to_csv(obs / "electron_track.csv", index=False)

    # (3) summary aliases dt_au / cell_bohr (builder key names)
    summ = run_dir / "run_summary.txt"
    txt = summ.read_text()
    add = []
    if "dt_au" not in txt:
        add.append("dt_au = 0.04")
    if "cell_bohr" not in txt:
        add.append("cell_bohr = 50x50x90")
    if add:
        with summ.open("a") as f:
            f.write("\n# --- run-notebook builder aliases (derived, non-canonical) ---\n")
            f.write("  ".join(add) + "\n")


def main():
    HYP.mkdir(parents=True, exist_ok=True)
    for ph in PHASES:
        run_dir = RESULTS / ph["name"]
        assert run_dir.is_dir(), f"missing run dir: {run_dir}"
        print(f"\n===== bridging + building {ph['name']} =====")
        bridge(run_dir)
        out = HYP / ph["nb"]
        B.build(
            str(run_dir), str(out),
            rs=GEOM["rs"], v0=GEOM["v0"], launch_z=GEOM["launch_z"],
            proj_sigma=GEOM["proj_sigma"], l_slab=GEOM["l_slab"],
            cap_inner=GEOM["cap_inner"], lindhard_mode="both",
            measured_s=ph["measured_s"], measured_v=ph["measured_v"],
            gif_seconds=16.0,
            twin_wp=str(WP_TWIN),     # WP−classical energy-diff bar GIF (slow)
            bar_gif_seconds=0.5,      # deliberately slow so the bars can be read
        )
        print(f"[done] {out}")


if __name__ == "__main__":
    main()
