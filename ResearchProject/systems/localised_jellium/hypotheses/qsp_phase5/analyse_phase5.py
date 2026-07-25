#!/usr/bin/env python3
"""qsp_phase5 — per-run quantum stopping power (WP energy method) + S(E) state.

ONE run at a time. Reads a completed WP run's observables, computes the quantum
stopping power S = [E_total(t_f) - E_GS] / L_z (the geometry-correct Method B for
a localised slab, with the WP-specific E_GS anchor — the WP's drift KE lives in
E_total(0), so anchoring to E_total(t0) is wrong), runs the convergence gate +
N-conservation guard from the `stopping-power-extraction` skill, then UPSERTS the
run's row into `se_state.csv` (the cumulative S(E) table the plot/email read).

Validated: reproduces analyse_phase4's energy_method exactly (2.39 eV/Bohr at 54 eV).

Usage:
    python3 analyse_phase5.py <results_dir>           # e.g. .../results/p5_wp_v3p0
    python3 analyse_phase5.py <results_dir> <tag> <v> # explicit overrides
"""
from __future__ import annotations
import json, os, re, sys
import importlib.util as _ilu
import numpy as np
import pandas as pd

HA      = 27.211386
# Bare-slab GS (Ha) — the WP energy-method anchor. Default = the h=0.5 production GS;
# override via P5_EGS for finer-grid reruns (the GS energy is grid-dependent).
E_GS    = float(os.environ.get("P5_EGS", "-70.22568216820937"))
L_SLAB  = 25.0                    # slab thickness / traversal length (Bohr)
SIGMA_WP = 0.5
HERE    = os.path.dirname(os.path.abspath(__file__))
STATE   = os.path.join(HERE, "se_state.csv")
ROOT    = "/local/data/public/skcb2/tddft"
SPK_PATH = os.path.join(ROOT, ".claude/skills/stopping-power-extraction/stopping_power.py")


def _load_spk():
    try:
        spec = _ilu.spec_from_file_location("spk", SPK_PATH)
        m = _ilu.module_from_spec(spec); spec.loader.exec_module(m); return m
    except Exception as exc:  # noqa: BLE001
        print(f"[analyse_phase5] WARN: stopping skill unavailable ({exc}); guard skipped")
        return None


def _read_summary(rdir):
    """Return (k0, E_drift_eV, wall_s) parsed from run_summary.txt (None if absent)."""
    p = os.path.join(rdir, "run_summary.txt")
    k0 = e_ev = wall = None
    if os.path.exists(p):
        txt = open(p).read()
        m = re.search(r"wp_k0\s*=\s*([-\d.eE+]+)", txt);        k0 = float(m.group(1)) if m else None
        m = re.search(r"wp_E_drift_eV\s*=\s*([-\d.eE+]+)", txt); e_ev = float(m.group(1)) if m else None
        m = re.search(r"wall_time_s\s*=\s*([-\d.eE+]+)", txt);   wall = float(m.group(1)) if m else None
    return k0, e_ev, wall


def analyse(rdir, tag=None, v_override=None):
    obs = os.path.join(rdir, "raw", "observables")
    o = pd.read_csv(os.path.join(obs, "observables.csv"))
    t, E = o["time_au"].values.astype(float), o["energy_total"].values.astype(float)
    s = pd.read_csv(os.path.join(obs, "wp_real_space_stats.csv"), comment="#")
    norm = s["norm_check"].values.astype(float)
    k0, e_ev, wall = _read_summary(rdir)
    v = float(v_override if v_override is not None else (k0 if k0 is not None else float("nan")))
    if e_ev is None and np.isfinite(v):
        e_ev = 0.5 * v * v * HA
    tag = tag or os.path.basename(rdir.rstrip("/"))

    # --- Method B (slab), WP-anchored to E_GS ---
    deposited = (float(E[-1]) - E_GS) * HA
    S = deposited / L_SLAB

    # --- convergence gate: WP fully absorbed (norm->0) AND E_total plateaued ---
    norm_f = float(norm[-1])
    m = t >= 0.85 * t.max()
    late_slope = float(np.polyfit(t[m], E[m], 1)[0] * HA) if m.sum() > 1 else float("nan")
    converged = bool(norm_f < 0.02 and abs(late_slope) < 0.2)
    # residual WP carries +energy still in box; CAP removes it -> E_total falls ->
    # converged deposit < current => current is an UPPER bound (late_slope<0 confirms).
    bound = "exact" if converged else ("upper" if late_slope < 0 else "lower")

    # --- N-conservation guard (skill) ---
    drained = float("nan"); guard_ok = None
    enf = os.path.join(obs, "electron_number.csv")
    if os.path.exists(enf):
        try:
            en = pd.read_csv(enf); N = en["N_total"].values.astype(float)
            spk = _load_spk()
            if spk is not None and N.size >= 2:
                g = spk.conservation_guard([float(N[0]), float(N[-1])])
                drained, guard_ok = g["drained_frac"], bool(g["ok"])
        except Exception as exc:  # noqa: BLE001
            print(f"[analyse_phase5] guard read failed: {exc}")

    res = dict(tag=tag, v=v, E_eV=float(e_ev) if e_ev is not None else float("nan"),
               S_eVbohr=float(S), deposited_eV=float(deposited),
               converged=converged, norm_f=norm_f, late_slope_eV_au=late_slope,
               bound=bound, N_drain_frac=drained, N_guard_ok=guard_ok,
               t_final_au=float(t[-1]), wall_s=wall, E_GS_Ha=E_GS, L_slab_bohr=L_SLAB)
    with open(os.path.join(HERE, f"results_{tag}.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    _upsert_state(res)
    print(f"[analyse_phase5] {tag}: v={v:.2f} E={res['E_eV']:.1f} eV  S={S:.3f} eV/Bohr "
          f"[{bound}{'' if converged else ' bound'}] norm_f={norm_f:.3f} slope={late_slope:+.3f} "
          f"Ndrain={drained if drained==drained else float('nan'):.4f}")
    return res


def _upsert_state(res):
    cols = ["tag", "v", "E_eV", "S_eVbohr", "deposited_eV", "converged",
            "norm_f", "late_slope_eV_au", "bound", "N_drain_frac", "wall_s"]
    row = {k: res.get(k) for k in cols}
    if os.path.exists(STATE):
        df = pd.read_csv(STATE)
        df = df[df["tag"] != res["tag"]]                 # drop any prior row for this tag
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    else:
        df = pd.DataFrame([row])
    df = df.sort_values("E_eV").reset_index(drop=True)
    df.to_csv(STATE, index=False)
    print(f"[analyse_phase5] se_state.csv now has {len(df)} point(s): "
          f"{', '.join(f'{e:.0f}eV' for e in df['E_eV'])}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(2)
    rdir = sys.argv[1]
    tag = sys.argv[2] if len(sys.argv) > 2 else None
    v = float(sys.argv[3]) if len(sys.argv) > 3 else None
    analyse(rdir, tag, v)
