#!/usr/bin/env python3
"""Assemble the sigma1_masspair notebooks (plan:
docs/plans/sigma1-masspair-decay-runs.md):

  wp_m2_k4p5_run_notebook.ipynb   per-run, m=2 (E=138 eV)
  wp_m3_k4p5_run_notebook.ipynb   per-run, m=3 (E=92 eV)
  sigma1_masspair_study.ipynb     cross-run: clock + kinetic-ledger attribution

All numbers computed live from the run CSVs; figures written to figs/ by the
notebooks' own code cells.

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_sigma1_masspair_report.py
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from _nbreport import md, code, embed, setup_cell, set_outdir, build

GIF_SYS = {"total": "total density", "bath": "bath (total − WP)", "wp": "wavepacket |ψ|²"}
GIF_VIEW = {"total": "n(x,z,t)", "dfirst": "Δn = n(t)−n(0)", "dprev": "Δn = n(t)−n(t−Δt)"}


def gif_cells(prefix, section="§7"):
    """Density animation GIFs (required in every run notebook — memory
    feedback_run_notebooks_require_density_gifs). Embeds what exists, links
    are relative so the GIFs travel beside the notebook in figs/."""
    cells = [md(f"""## {section} — Density animations (x–z plane, y=0)

{{total, bath = total − WP, wavepacket}} × {{n, Δn vs t=0, Δn per frame}}.
Dashed lines: slab faces ±12.5; dotted: CAP inner faces ±35.""")]
    missing = []
    for key in ("total", "bath", "wp"):
        for view in ("total", "dfirst", "dprev"):
            p = os.path.join(HERE, "figs", f"{prefix}_{key}_{view}.gif")
            if os.path.exists(p):
                cells.append(embed(p, f"{GIF_SYS[key]} · {GIF_VIEW[view]}", width=360))
            else:
                missing.append(os.path.basename(p))
    if missing:
        cells.append(md("*Missing GIFs (regenerate via `_density_views.render_decomposition_views`): "
                        + ", ".join(f"`{m}`" for m in missing) + "*"))
    return cells

WPRES = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
         "scripts/sigma1_masspair/wp/results")
P3RES = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
         "scripts/qsp_phase3/wp/results/p3_wp_m1_rerun")

COMMON = f"""import pandas as pd, numpy as np, matplotlib.pyplot as plt, os
HA = 27.211386
WPRES = {WPRES!r}
P3RES = {P3RES!r}
os.makedirs("figs", exist_ok=True)
def load(run):
    d = {{}}
    d["obs"] = pd.read_csv(f"{{WPRES}}/{{run}}/raw/observables/observables.csv")
    d["ix"]  = pd.read_csv(f"{{WPRES}}/{{run}}/raw/observables/interactions.csv")
    d["rs"]  = pd.read_csv(f"{{WPRES}}/{{run}}/raw/observables/wp_real_space_stats.csv", skiprows=1)
    d["ms"]  = pd.read_csv(f"{{WPRES}}/{{run}}/raw/observables/wp_momentum_stats.csv", skiprows=1)
    d["ne"]  = pd.read_csv(f"{{WPRES}}/{{run}}/raw/observables/electron_number.csv")
    d["se"]  = pd.read_csv(f"{{WPRES}}/{{run}}/raw/observables/state_energies.csv")
    return d"""


def summary_cell(run):
    return code(f"""rows = []
for line in open({WPRES!r} + "/{run}/run_summary.txt"):
    line = line.rstrip("\\n")
    if not line.strip(): continue
    k, _, v = line.partition("=")
    rows.append((k.strip(), v.strip()))
from IPython.display import Markdown
Markdown("| key | value |\\n|---|---|\\n" + "\\n".join(
    f"| {{k}} | {{v.replace('|', chr(92)+'|')}} |" for k, v in rows))""")


def run_nb_cells(run, m, ekin):
    v = 4.5 / m
    cells = [setup_cell()]
    cells.append(md(f"""# sigma1_masspair · `{run}` — run notebook

**Config:** clean qsp_phase3 geometry (50×50×90, dx 0.5, N=82, two-sided CAP
η=−0.7 at ±35..±45, dt 0.04, 2500 steps → τ=100); WP density-width σ_ρ=1.0
Bohr (`sigma(√2)`, house label σ_WP≈1.41 — see plan "σ convention
correction"), k0=4.5, **m={m}** (v={v:.2f}, E={ekin} eV), launch z=−16.5
(4 Bohr from the slab face). Full ledger + per-step pairwise `interactions.csv`
+ checkpoints every 200 steps.

**Goal was a monotonic energy decay (as the σ=0.5/m=1 clean run). Headline:
the drain-then-rise artifact RETURNS in this clean geometry — see §2/§6.**

Plan: `docs/plans/sigma1-masspair-decay-runs.md` · Handover:
`docs/handovers/energy-oscillation-debugging.md`"""))
    cells.append(md("## §1 — Run summary (verbatim `run_summary.txt`)"))
    cells.append(summary_cell(run))
    cells.append(md("""## §2 — Energy ledger

Left: ΔE_total(t). Right: the four ledger components (Δ from t=0). The
post-minimum rise lives entirely in `energy_kinetic`."""))
    cells.append(code(COMMON + f"""
d = load({run!r})
o = d["obs"]; t, E = o.time_au.values, o.energy_total.values
im = int(E.argmin())
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(t, (E-E[0])*HA)
ax[0].axvline(t[im], ls=":", color="0.5", label=f"E_min t={{t[im]:.1f}}")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("dE_total (eV)"); ax[0].legend()
for c in ("energy_kinetic","energy_hartree","energy_xc","energy_external"):
    ax[1].plot(t, (o[c]-o[c].iloc[0])*HA, label=c.replace("energy_",""))
ax[1].axvline(t[im], ls=":", color="0.5")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("d(component) (eV)"); ax[1].legend()
fig.tight_layout(); fig.savefig("figs/{run}_ledger.png", dpi=150); plt.show()
print(f"drain {{(E[0]-E[im])*HA:.1f}} eV   E_min at t={{t[im]:.1f}}   "
      f"rise to tau {{(E[-1]-E[im])*HA:+.1f}} eV")
print(f"d(kinetic) over the rise: {{(o.energy_kinetic.iloc[-1]-o.energy_kinetic.iloc[im])*HA:+.1f}} eV")"""))
    cells.append(md(f"""## §3 — WP transport and spreading

σ_z(t) against free flight with the *measured* birth σ_pz; the slab face is
reached at t≈{4/v:.1f}. Free-flight prediction stops applying inside the
slab."""))
    cells.append(code(f"""rs, ms = d["rs"], d["ms"]
tr = rs.time_au.values
sz = np.sqrt(rs.z2_mean.values - rs.z_mean.values**2)
sp = float(np.sqrt(ms.sigma_pz2.iloc[0])); s0 = sz[0]
pred = np.sqrt(s0**2 + (sp*tr/{m})**2)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(tr, rs.z_mean.values); ax[0].axhline(-12.5, ls="--", color="0.6", label="slab face")
ax[0].axhline(35, ls="--", color="0.4", label="CAP inner face")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("z_mean (Bohr)"); ax[0].legend()
ax[1].plot(tr, sz, label="measured sigma_z")
ax[1].plot(tr, pred, "--", label="free flight (measured sigma_pz)")
ax[1].set_xlim(0, 12); ax[1].set_ylim(0, 6)
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("sigma_z (Bohr)"); ax[1].legend()
fig.tight_layout(); fig.savefig("figs/{run}_transport.png", dpi=150); plt.show()
zi = int(np.argmax(rs.z_mean.values >= -12.5))
print(f"sigma_pz(0)={{sp:.3f}} (min-unc 0.50)   spreading at face: {{(sz[zi]/s0-1)*100:.1f}}%")"""))
    cells.append(md("""## §4 — Absorption: N_total and the WP state norm

t_min coincides with the WP norm crossing ~5 % (the CAP has essentially
finished absorbing the packet); the kinetic rise happens as norm → 1e-4."""))
    cells.append(code(f"""ix, ne = d["ix"], d["ne"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(ne.time_au, ne.N_total); ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("N_total")
ax[1].semilogy(ix.time_au, ix.norm_wp)
ax[1].axvline(t[im], ls=":", color="0.5", label=f"t_min={{t[im]:.1f}}")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("norm_wp (log)"); ax[1].legend()
fig.tight_layout(); fig.savefig("figs/{run}_absorption.png", dpi=150); plt.show()
j = int(np.searchsorted(ix.time_au.values, t[im]))
print(f"absorbed {{ne.N_total.iloc[0]-ne.N_total.iloc[-1]:.3f}} e   "
      f"norm_wp at t_min = {{ix.norm_wp.iloc[j]:.4f}}   at tau = {{ix.norm_wp.iloc[-1]:.2e}}")"""))
    cells.append(md("""## §5 — Pairwise Coulomb decomposition (P/S/B)

Every channel Δ from t=0, plus the closure residual against INQ's ledger
(must be ~1e-9). The Coulomb channels are QUIET through the rise — the
artifact is not electrostatic."""))
    cells.append(code(f"""fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for c in ("e_pp","e_ps","e_ss","e_sb","e_pb"):
    ax[0].plot(ix.time_au, (ix[c]-ix[c].iloc[0])*HA, label=c)
ax[0].axvline(t[im], ls=":", color="0.5")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("d(channel) (eV)"); ax[0].legend(ncol=2)
mrg = o.merge(ix, on="step", suffixes=("", ".ix"))
ax[1].semilogy(mrg.time_au, (mrg.e_hartree_check-mrg.energy_hartree).abs()+1e-16, label="|hartree closure|")
ax[1].semilogy(mrg.time_au, (mrg.e_external_check-mrg.energy_external).abs()+1e-16, label="|external closure|")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("closure residual (Ha)"); ax[1].legend()
fig.tight_layout(); fig.savefig("figs/{run}_pairwise.png", dpi=150); plt.show()
print("max closure:", f"{{(mrg.e_hartree_check-mrg.energy_hartree).abs().max():.1e}}",
      f"{{(mrg.e_external_check-mrg.energy_external).abs().max():.1e}} Ha")"""))
    cells.append(md(f"""## §6 — Verdict

- The energy does **not** decay to a fixed value: it drains, reaches its
  minimum exactly as the CAP finishes absorbing the WP state (norm ≈ 5 %),
  then **rises** — the rise sits entirely in the norm-divided kinetic ledger
  term while all pairwise Coulomb channels stay quiet.
- Mechanically the run is clean: full 2500 steps, closure ≤ 1e-9 Ha,
  checkpoints + segments intact, WP fully absorbed (~1.0 e).
- Cross-run implications (clock scaling, mass correlation, comparison to the
  clean m=1 runs) are in `sigma1_masspair_study.ipynb`."""))
    cells.extend(gif_cells(run, section="§7"))
    return cells


def study_cells():
    cells = [setup_cell()]
    cells.append(md("""# sigma1_masspair study — the artifact returns in the clean geometry

**Design goal** (plan `docs/plans/sigma1-masspair-decay-runs.md`): reproduce
the clean run's monotonic energy decay with a σ_ρ=1, heavier-mass packet.
**Result: refuted** — both runs show the drain-then-rise artifact at full
strength, in the geometry that decays monotonically at σ=0.5/m=1/100 eV.
That failure is highly informative; this notebook extracts the three new
facts.

| run | m | v | E_kin | drain | t_min | rise after min |
|---|---|---|---|---|---|---|
| wp_m2_k4p5 | 2 | 2.25 | 138 eV | 60 eV | 32.6 | **+180 eV** |
| wp_m3_k4p5 | 3 | 1.50 | 92 eV | 44 eV | 50.2 | **+125 eV** |
| p3_wp_m1_rerun (clean ref) | 1 | 2.71 | 100 eV | 126 eV | 97.9 | +0.11 eV |
"""))
    cells.append(md("## §1 — ΔE_total overlays (both runs + the clean m=1 reference)"))
    cells.append(code(COMMON + """
runs = {"wp_m2_k4p5": 2, "wp_m3_k4p5": 3}
D = {r: load(r) for r in runs}
p3 = pd.read_csv(f"{P3RES}/raw/observables/observables.csv")
fig, ax = plt.subplots(figsize=(8, 4.5))
for r in runs:
    o = D[r]["obs"]
    ax.plot(o.time_au, (o.energy_total-o.energy_total.iloc[0])*HA, label=f"{r} (m={runs[r]})")
ax.plot(p3.time_au, (p3.energy_total-p3.energy_total.iloc[0])*HA, "k--", label="p3_wp_m1_rerun (clean, m=1)")
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("dE_total (eV)"); ax.legend()
fig.tight_layout(); fig.savefig("figs/study_overlay.png", dpi=150); plt.show()"""))
    cells.append(md("""## §2 — The clock: t_min tracks the projectile's CAP arrival

t_min ratio between the runs equals the velocity ratio, and in both runs
t_min ≈ 1.44 × the ballistic launch→far-CAP-face time — while the WP-state
norm crosses the same ~5 % at t_min. The clock is the completion of CAP
absorption of the projectile, not slow-spill arrival (which the earlier
campaigns assumed)."""))
    cells.append(code("""rows = []
fig, ax = plt.subplots(figsize=(8, 4))
for r, m in runs.items():
    o, ix = D[r]["obs"], D[r]["ix"]
    t, E = o.time_au.values, o.energy_total.values
    im = int(E.argmin()); v = 4.5/m
    t_cap = (35-(-16.5))/v
    j = int(np.searchsorted(ix.time_au.values, t[im]))
    rows.append((r, m, v, t[im], t_cap, t[im]/t_cap, ix.norm_wp.iloc[j]))
    ax.semilogy(ix.time_au, ix.norm_wp, label=f"{r}")
    ax.axvline(t[im], ls=":", color="0.5")
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("norm_wp (log)"); ax.legend()
fig.tight_layout(); fig.savefig("figs/study_clock.png", dpi=150); plt.show()
print(f"{'run':14s} {'m':>3} {'v':>5} {'t_min':>6} {'t_cap':>6} {'ratio':>6} {'norm@t_min':>10}")
for r, m, v, tm, tc, ra, nw in rows:
    print(f"{r:14s} {m:3d} {v:5.2f} {tm:6.1f} {tc:6.1f} {ra:6.2f} {nw:10.4f}")
print(f"t_min ratio {rows[1][3]/rows[0][3]:.2f} vs velocity ratio {rows[0][2]/rows[1][2]:.2f}")"""))
    cells.append(md("""## §3 — Attribution: the rise is 100 % the kinetic ledger term

Change of each ledger component and each pairwise Coulomb channel over
t_min → τ. Kinetic carries the entire rise; Hartree/external nearly cancel
in their slab parts (e_ss vs e_sb); every projectile channel is ≤ 1.4 eV."""))
    cells.append(code("""labels, kin, coul = [], [], {}
for r, m in runs.items():
    o, ix = D[r]["obs"], D[r]["ix"]
    t, E = o.time_au.values, o.energy_total.values
    im = int(E.argmin()); j = int(np.searchsorted(ix.time_au.values, t[im]))
    labels.append(r)
    kin.append((o.energy_kinetic.iloc[-1]-o.energy_kinetic.iloc[im])*HA)
    for c in ("e_pp","e_ps","e_ss","e_sb","e_pb"):
        coul.setdefault(c, []).append((ix[c].iloc[-1]-ix[c].iloc[j])*HA)
x = np.arange(len(labels)); w = 0.13
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(x-2.5*w, kin, w, label="d kinetic")
for k, (c, vals) in enumerate(coul.items()):
    ax.bar(x+(k-1.5)*w, vals, w, label="d "+c)
ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylabel("change t_min -> tau (eV)")
ax.legend(ncol=3); fig.tight_layout(); fig.savefig("figs/study_attribution.png", dpi=150); plt.show()
for i, r in enumerate(labels):
    print(f"{r}: d(kinetic)={kin[i]:+.1f} eV;  " +
          "  ".join(f"{c}={coul[c][i]:+.2f}" for c in coul))"""))
    cells.append(md("""## §4 — The mechanism at state level

The WP state's raw energy expectation ⟨φ|H|φ⟩ (occupation-weighted) decays to
zero as the CAP absorbs it — but INQ's kinetic ledger divides by the state
norm (`energy.hpp`: occ·⟨φ|T|φ⟩/⟨φ|φ⟩). The *ratio* for the surviving sliver
climbs (≈5.5 → ≈14 Ha for m=2): the CAP preferentially removes the packet's
core and leaves a high-kinetic residue whose normalised energy diverges as
norm → 1e-4. The m=1 clean runs end with a *slow* residue instead (rise
0.11 eV). Why the massive/faster-k0 packet leaves a fast residue (reflected
high-k components near Nyquist?) is the open question."""))
    cells.append(code("""fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for r, m in runs.items():
    se = D[r]["se"]; w = se[se.state_index==60]
    ax[0].plot(w.time_au, w.E_expect_ha, label=f"{r} raw <H>")
    ix = D[r]["ix"]
    nrm = np.interp(w.time_au.values, ix.time_au.values, ix.norm_wp.values)
    ax[1].semilogy(w.time_au, w.E_expect_ha/np.maximum(nrm, 1e-12), label=f"{r} <H>/norm")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("state-60 <H> (Ha, occ-weighted)"); ax[0].legend()
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("<H>/norm (Ha, log)"); ax[1].legend()
fig.tight_layout(); fig.savefig("figs/study_state60.png", dpi=150); plt.show()"""))
    cells.append(md("""## §5 — Implications for the oscillation investigation

1. **H1 (large CAP standoff protects) is refuted as sufficient** — identical
   clean geometry, violent artifact.
2. **The mass/k0 correlation is restored**: every oscillating run to date has
   m_eff>1 (2.10 family; now 2 and 3) and/or k0 well above 2.7; every clean
   run is m=1, k0=2.711. The engine-level ledger audit (m=1 bit-identity)
   does not cover the m>1 dynamics.
3. **New clock law**: t_min = completion of CAP absorption of the projectile
   state (norm ≈ 5 %), scaling as 1/v — reinterprets the cap_fix
   "period-lengthening ladder" as absorption-completion shifts.
4. Sharpest next discriminators: (a) m=1 twin at k0=4.5/σ_ρ=1 in this exact
   geometry (separates mass from k0/σ); (b) recompute E_total excluding the
   WP state's norm-divided term (zero-GPU, from state_energies +
   observables) to confirm the ledger identity; (c) τ-extension of a clean
   run (the checkpoints allow it directly).

Review path: per-run notebooks → this study → handover
`docs/handovers/energy-oscillation-debugging.md`.

**Density animations** (per-run GIFs, §7 of each run notebook — passive links):
[m2 total](figs/wp_m2_k4p5_total_total.gif) ·
[m2 bath Δn](figs/wp_m2_k4p5_bath_dfirst.gif) ·
[m2 wp](figs/wp_m2_k4p5_wp_total.gif) ·
[m3 total](figs/wp_m3_k4p5_total_total.gif) ·
[m3 bath Δn](figs/wp_m3_k4p5_bath_dfirst.gif) ·
[m3 wp](figs/wp_m3_k4p5_wp_total.gif)"""))
    return cells


set_outdir(HERE)
build(run_nb_cells("wp_m2_k4p5", 2, 138), os.path.join(HERE, "wp_m2_k4p5_run_notebook.ipynb"), timeout=900)
build(run_nb_cells("wp_m3_k4p5", 3, 92), os.path.join(HERE, "wp_m3_k4p5_run_notebook.ipynb"), timeout=900)
build(study_cells(), os.path.join(HERE, "sigma1_masspair_study.ipynb"), timeout=900)
