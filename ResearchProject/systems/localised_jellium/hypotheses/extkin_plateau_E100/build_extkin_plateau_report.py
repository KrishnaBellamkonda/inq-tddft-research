#!/usr/bin/env python3
"""Study-notebook builder — extkin_plateau_E100 (first jellium run with the
IN-RUN norm-division fix via OrbitalKineticStats).

Reads scripts/extkin_plateau_E100/wp/results/cap; builds
extkin_plateau_E100_study.ipynb (executed, 0 errors, GIF path-referenced).
Partial-tolerant: builds from whatever files exist.

Run:  PYTHONPATH=$ROOT/inq-stack/python $ROOT/venv/bin/python3 build_extkin_plateau_report.py
"""
import importlib.util
import os
import sys
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
SYS = ROOT / "ResearchProject/systems/localised_jellium"
RUN = SYS / "scripts/extkin_plateau_E100/wp/results/cap"
GS = SYS / "shared_gs/slab_n92_L35x35x120_w0p5_h0p5"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "inq-stack/python"))

spec = importlib.util.spec_from_file_location("_nbreport", SYS / "hypotheses/_nbreport.py")
nb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nb)
nb.set_outdir(HERE)

# _nbreport cells are RETURNED (not appended); collect them here. `anchor` is
# accepted for call-site readability but derived automatically by tag_builder.
CELLS = []


def md(text, anchor=None):
    CELLS.append(nb.md(text))


def code(src, anchor=None):
    CELLS.append(nb.code(src))

# ---------------------------------------------------------------- pre: GIF ---
# Density-matrix GIF (rule notebook-density-gif): xz mid-y slice, total + wp +
# induced, shared log scale for total/wp, per-frame linear for the WP is NOT
# used here (dispersing WP -> wide log per reference_dispersing_wp_gif_scaling).
GIF = HERE / "fig_extkin_plateau_density.gif"


def build_gif():
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, TwoSlopeNorm
    from PIL import Image
    from inqview.visualisation.field_io import load_vti

    vt = sorted((RUN / "raw/vti/density_total").glob("*.vti"))
    vw = sorted((RUN / "raw/vti/density_wp").glob("*.vti"))
    if not vt:
        print("no density frames — GIF skipped")
        return False
    # subsample to <=60 frames
    step = max(1, len(vt) // 60)
    vt, vw = vt[::step], vw[::step]
    f0 = load_vti(str(vt[0]), expect_centered_axis="z")
    d0, x, z = f0.data, f0.x, f0.z
    iy = d0.shape[1] // 2
    tot0 = d0[:, iy, :]
    vmax = max(float(np.nanmax(load_vti(str(p)).data[:, iy, :])) for p in vt[:: max(1, len(vt)//8)])
    vmin = vmax * 1e-7
    frames = []
    for i, p in enumerate(vt):
        ft = load_vti(str(p)).data[:, iy, :]
        fw = load_vti(str(vw[i])).data[:, iy, :] if i < len(vw) else np.zeros_like(ft)
        ind = ft - tot0
        fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
        for a, f, ttl, kind in ((ax[0], ft, "total (log)", "log"),
                                 (ax[1], fw, "wp (log)", "log"),
                                 (ax[2], ind, "induced Δn", "lin")):
            ext = [z[0], z[-1], x[0], x[-1]]
            if kind == "log":
                a.imshow(np.clip(f, vmin, None), origin="lower", aspect="auto",
                         extent=ext, norm=LogNorm(vmin=vmin, vmax=vmax), cmap="inferno")
            else:
                m = max(1e-12, float(np.nanmax(np.abs(ind))))
                a.imshow(f, origin="lower", aspect="auto", extent=ext,
                         norm=TwoSlopeNorm(0, -m, m), cmap="RdBu_r")
            for zz in (-10.07, 10.07):
                a.axvline(zz, color="w" if kind == "log" else "k", ls="--", lw=0.7)
            for zz in (-45, 45):
                a.axvline(zz, color="cyan", ls=":", lw=0.7)
            a.set_title(ttl, fontsize=9)
            a.set_xlabel("z (Bohr)")
        ax[0].set_ylabel("x (Bohr)")
        fig.suptitle(f"n(x,z) mid-y   {p.stem}", fontsize=9)
        fig.tight_layout()
        fig.canvas.draw()
        w, h = fig.canvas.get_width_height()
        frames.append(Image.frombuffer("RGBA", (w, h), fig.canvas.buffer_rgba()).convert("P"))
        plt.close(fig)
    frames[0].save(GIF, save_all=True, append_images=frames[1:], duration=120, loop=0)
    print(f"GIF: {GIF.name} ({len(frames)} frames)")
    return True


have_gif = GIF.exists() or build_gif()

# ------------------------------------------------------------------- cells ---
md("# extkin_plateau_E100 — E_plateau with the norm-division fix IN-RUN",
   anchor="title")
md("""**Question.** How much total electronic energy does a σ=1.5 Bohr, 100 eV
wavepacket deposit in a localised-jellium slab (N=92, r_s=4.0, 20.13 Bohr
thick), measured WITHOUT the per-particle-kinetic artifact — i.e. what is
E_plateau of the corrected extensive energy after the packet is absorbed by
the two-sided CAP?

This is the FIRST jellium run where the fix runs in-run
(`inqkit::observables::OrbitalKineticStats`, all states, every step), replacing
the post-hoc `E_ext = E_reported − e_kin_ha·(1−norm)` route of the replica
campaign. Design + decision log:
`docs/plans/norm-corrected-stopping-power.md` "Run design (2026-07-29)".

⚠ **Recorded caveats (user-approved scope cuts):** first CAP run at dt=0.04
(absorption quality at this dt not separately gated); no no-CAP twin (no
in-system conservation control).""", anchor="question")

md("""## Conventions

Hartree atomic units unless stated; 1 Ha = 27.211 eV.

| symbol | meaning | value |
|---|---|---|
| σ | WP wavefunction width (σ_WP convention) | 1.5 Bohr |
| k₀ = v | drift momentum = velocity (m=1) | 2.711 a.u. |
| η, W | CAP strength / width per side | −1.0 Ha, 15 Bohr |
| E_corr | total − kinetic + kin_bare (extensive total) | per step |
| E_plateau | ⟨E_corr − E_GS⟩ over the post-absorption window | result |
""", anchor="conventions")

code("""import numpy as np, pandas as pd, re
from pathlib import Path
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.visualisation.style import apply_theme
apply_theme()
HA = 27.211386
RUN = Path("%s"); GS = Path("%s")
""" % (RUN, GS), anchor="imports")

md("Dispersion time of a free Gaussian packet: τ = σ₀² (width ratio "
   "R(t) = √(1+(t/τ)²)). At σ=1.5 the packet reaches the slab face "
   "(7.5 Bohr, t=2.77) already ×1.58 wider — the compact-projectile choice, "
   "accepted in the design interview.", anchor="disp-formula")
code("""sigma, v = 1.5, 2.711
tau = sigma**2
R_entry = np.sqrt(1 + (7.5/v/tau)**2)
print(f"tau = {tau:.2f} a.u.;  growth at slab entry x{R_entry:.2f}")""",
     anchor="disp-calc")

md("## Simulation setup (verbatim run_summary.txt)", anchor="setup")
code("""p = RUN/"run_summary.txt"
if p.exists():
    rows = [l.split("=",1) for l in p.read_text().splitlines() if "=" in l]
    display(pd.DataFrame([(a.strip(),b.strip()) for a,b in rows],
                         columns=["key","value"]))
else:
    print("RUN INCOMPLETE — no run_summary.txt yet")""", anchor="setup-table")

md("""## Source files

| role | path |
|---|---|
| run definition | `ResearchProject/systems/localised_jellium/scripts/extkin_plateau_E100/wp/run.cpp` |
| GS definition | `.../scripts/extkin_plateau_E100/gs/run.cpp` |
| dispatcher | `.../scripts/extkin_plateau_E100/run_extkin_plateau.sh` |
| config (SoT) | `.../shared/configs/slab_n92_L35x35x120_w0p5.hpp` |
| observable | `inq-stack/include/inqkit/observables/orbital_kinetic_stats.hpp` |
| this builder | `.../hypotheses/extkin_plateau_E100/build_extkin_plateau_report.py` |
| vacuum validation | `.../systems/vacuum/hypotheses/cap_norm_investigation/extensive_kinetic/extkin_study.ipynb` |
""", anchor="sources")

if have_gif:
    md("## Visual intuition — density evolution (xz, mid-y)\n\n"
       "Dashed white/black: slab faces ±10.07; dotted cyan: CAP inner edges "
       "±45.\n\n![density evolution](fig_extkin_plateau_density.gif)",
       anchor="gif")
else:
    md("## Visual intuition\n\n*(no density frames found at build time — "
       "rebuild after the run completes)*", anchor="gif")

md("""## Result 1 — corrected vs reported total energy

E_corr(t) = total − kinetic + kin_bare_total, all from in-run CSVs. The
reported total pins near its initial value as the CAP removes the WP; the
corrected series must drop by the WP's remaining (undeposited) energy and go
FLAT once absorption completes.""", anchor="r1")
code("""en = pd.read_csv(RUN/"raw/observables/energies.csv")
ek = pd.read_csv(RUN/"raw/observables/orbital_kinetic_stats.csv", comment="#")
m = en.merge(ek, on="step", suffixes=("","_ek"))
t = m["time_au"].to_numpy()
E_rep = m["total"].to_numpy(); K_rep = m["kinetic"].to_numpy()
K_bare = m["kin_bare_total_ha"].to_numpy()
E_corr = E_rep - K_rep + K_bare
ident = (m["kin_normdiv_total_ha"] - m["kinetic"]).abs().max()
print(f"identity max|kin_normdiv - kinetic_INQ| = {ident:.2e} Ha")
gs_txt = (GS/"run_summary.txt").read_text()
E_gs = float(re.search(r"gs_energy_ha\\s*=\\s*([\\-0-9.eE+]+)", gs_txt).group(1))
fig, ax = plt.subplots(figsize=(8,4))
ax.plot(t, (E_rep-E_gs)*HA, label="reported − E_GS")
ax.plot(t, (E_corr-E_gs)*HA, label="corrected (extensive) − E_GS")
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("E − E_GS (eV)"); ax.legend()
fig.savefig("fig_energy_corrected.png", dpi=160)""", anchor="r1-code")

md("""## Result 2 — E_plateau

Plateau window: last 10 a.u. AND WP-orbital norm < 1e-3 (fully absorbed).
The WP orbital is the last state column in orbital_kinetic_stats.""",
   anchor="r2")
code("""ncols = [c for c in ek.columns if c.startswith("norm_")
         and c != "norm_total"]
wp_norm = m[ncols[-1]].to_numpy()

rows = []
for lo, hi in [(40, 50), (50, 55), (55, 60)]:
    w = (t >= lo) & (t < hi)
    v = (E_corr[w] - E_gs) * HA
    rows.append((f"[{lo},{hi})", v.mean(), np.polyfit(t[w], v, 1)[0],
                 wp_norm[w][-1]))
tab = pd.DataFrame(rows, columns=["t window (a.u.)", "mean E−E_GS (eV)",
                                  "drift (eV/a.u.)", "WP norm at end"])
display(tab.round({"mean E−E_GS (eV)": 2, "drift (eV/a.u.)": 3}))

E_pl = rows[-1][1]; drift = rows[-1][2]
bound = wp_norm[-1] * 100.4   # residual norm × injected energy (eV)
converged = abs(drift) < 0.02 and bound < 0.5
print(f"E_plateau = {E_pl:.1f} eV  (final-window drift {drift:+.3f} eV/a.u.; "
      f"residual WP norm {wp_norm[-1]:.1e} bounds further drain at "
      f"~{bound:.1f} eV)")
print("CONVERGED" if converged else
      "PROVISIONAL — drift or residual norm too large; extend via WP_RESUME=1")""",
     anchor="r2-code")

md("## Result 3 — observable cost", anchor="r3")
code("""w = ek["wall_ms"].to_numpy()
summ = (RUN/"run_summary.txt").read_text() if (RUN/"run_summary.txt").exists() else ""
ps = re.search(r"per_step_ms\\s*=\\s*([0-9.eE+-]+)", summ)
print(f"OrbitalKineticStats: {np.mean(w[1:]):.1f} ms/step "
      f"({len(ek)} rows, {len(ncols)} states)")
if ps: print(f"run per-step: {float(ps.group(1)):.0f} ms -> overhead "
             f"{100*np.mean(w[1:])/float(ps.group(1)):.1f}%")""",
     anchor="r3-code")

md(r"""## Orbital-free stopping power — $S = E_\mathrm{plateau} / L_\mathrm{slab}$

The primary (orbital-free) stopping power divides the TOTAL deposited energy by
the slab thickness. This is the identifiability-safe route: it needs no notion
of "the projectile orbital" after injection.

**Analytical expectation (bulk, point charge).** Linear response for a charge
$Z$ at velocity $v$ in a homogeneous electron gas (Lindhard 1954; Echenique,
Flores & Ritchie, Sol. St. Phys. 43, 1990):

$$S=\frac{2Z^2}{\pi v^2}\int\frac{dq}{q}\int_0^{qv} d\omega\,\omega\,
\mathrm{Im}\!\left[\frac{-1}{\varepsilon(q,\omega)}\right]$$

For $v\gg v_F$ (here $v/v_F=5.7$) the loss function is plasmon-pole dominated
and the $q$-integral collapses to a log between $q_\mathrm{min}=\omega_p/v$
(adiabatic cutoff) and $q_\mathrm{max}=2v$ (closest binary collision):

$$S \approx \frac{\omega_p^2}{v^2}\,\ln\!\frac{2v^2}{\omega_p}$$""",
   anchor="sfree")
code("""L_SLAB = 20.134            # Bohr (2*half_width, r_s=4.0 exact)
S_free = E_pl / L_SLAB     # E_pl from Result 2 (eV)
n0 = 3/(4*np.pi*4.0**3); w_p = np.sqrt(4*np.pi*n0); v = 2.711
S_bethe = (w_p**2/v**2)*np.log(2*v**2/w_p)*HA
print(f"S_free   = {E_pl:.2f} eV / {L_SLAB:.2f} Bohr = {S_free:.2f} eV/Bohr")
print(f"S_Bethe  = (w_p^2/v^2)ln(2v^2/w_p) = {S_bethe:.2f} eV/Bohr "
      f"(w_p={w_p*HA:.1f} eV, ln={np.log(2*v**2/w_p):.1f})")
print(f"ratio sim/analytic = {S_free/S_bethe:.2f}")
print("Known suppressions vs bulk point-charge (inference): slab (20 Bohr) << "
      "wake wavelength 2*pi*v/w_p = "
      f"{2*np.pi*v/w_p:.0f} Bohr (collective channel barely develops); "
      "packet form factor cuts q > ~1/sigma_rho ~ 0.5 (close collisions); "
      "2-3 subbands / surface spill-out.")""", anchor="sfree-calc")

md(r"""## Orbital-dependent stopping power

Here we attribute a trajectory and a kinetic energy to the WP's KS orbital
(state 61) and extract $S=-dE_\mathrm{kin}/dz$ along its path — the scheme of
the stopping-power-extraction skill, applied to a quantum orbital.

> **Caveat (standing project decision):** after injection the WP orbital is
> not strictly identifiable as "the projectile" (it hybridises with the bath),
> so this is a SECONDARY estimate; the orbital-free $S$ above is primary.

### Position of the packet — two routes (both available in this run)

1. **Momentum integration**: $z_p(t) = z_0 + \int_0^t \langle p_z\rangle\,dt'$
   ($m=1$), with $\langle p_z\rangle$ every step from `wp_momentum_stats.csv`.
2. **Centroid (center of density)**: directly recorded as `z_mean` every 100
   steps (`wp_real_space_stats.csv`), densified to every 10 steps from the
   `density_wp` VTI frames.

Both are norm-weighted means, so they stay meaningful only while the orbital
norm is ~1 (before CAP contact); fits below are restricted to that window.""",
   anchor="sorb-pos")
code("""mom = pd.read_csv(RUN/"raw/observables/wp_momentum_stats.csv", comment="#")
tm, pz = mom["time_au"].to_numpy(), mom["pz_mean"].to_numpy()
z0 = -17.5
z_p = z0 + np.concatenate([[0.0], np.cumsum(0.5*(pz[1:]+pz[:-1])*np.diff(tm))])

rs = pd.read_csv(RUN/"raw/observables/wp_real_space_stats.csv", comment="#")

from inqview.visualisation.field_io import load_vti
zc_t, zc, zc_mass = [], [], []
for p in sorted((RUN/"raw/vti/density_wp").glob("*.vti")):
    f = load_vti(str(p))
    prof = f.data.sum(axis=(0, 1))
    m0 = prof.sum()
    zc_t.append(int(p.stem.split("_t")[1])*0.04)
    zc.append(float((prof*f.z).sum()/m0)); zc_mass.append(m0)
zc_t, zc, zc_mass = map(np.array, (zc_t, zc, zc_mass))
zc_mass = zc_mass/zc_mass[0]

ok = zc_mass > 0.995                       # centroid valid while norm ~1
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tm, z_p, label=r"$z_0+\\int\\langle p_z\\rangle dt$ (every step)")
ax.plot(zc_t[ok], zc[ok], "o", ms=3, label="VTI centroid (every 10 steps)")
ax.plot(zc_t[~ok], zc[~ok], "o", ms=3, mfc="none", alpha=0.4,
        label="centroid, norm<0.995 (biased by CAP)")
ax.plot(rs["time_au"], rs["z_mean"], "s", ms=5, mfc="none",
        label="z_mean (wp_real_space_stats)")
ax.axhspan(-10.07, 10.07, color="0.9", label="slab")
ax.axhline(45, color="c", ls=":"); ax.axhline(-45, color="c", ls=":")
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("z (Bohr)"); ax.legend(fontsize=8)
fig.savefig("fig_wp_position.png", dpi=160)
zi = np.interp(zc_t[ok], tm, z_p)
print(f"route agreement while norm>0.995: max |z_p - centroid| = "
      f"{np.max(np.abs(zi - zc[ok])):.3f} Bohr")""", anchor="sorb-pos-calc")

md("""### Kinetic energy of the WP orbital vs time

`e_kin_ha` is the norm-weighted mean KE of state 61 (cross-checked against
`tkin_61/norm_61` from the extkin observable). Free dispersion conserves the
mean KE, so the curve should be FLAT in vacuum and drop only in the slab;
entering the attractive background well transiently RAISES the KE
(acceleration by the well), returned on exit.""", anchor="sorb-ke")
code("""ke_mom = np.interp(t, tm, mom["e_kin_ha"].to_numpy())*HA
ke_ext = (m["tkin_61"]/m["norm_61"]).to_numpy()*HA
print(f"cross-check e_kin_ha vs tkin_61/norm_61: max diff "
      f"{np.max(np.abs(ke_mom-ke_ext)):.2e} eV")
wpn61 = m["norm_61"].to_numpy()
val = wpn61 > 0.995
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t[val], ke_ext[val], label="mean KE of state 61 (norm>0.995)")
ax.plot(t[~val], ke_ext[~val], alpha=0.3, label="remnant-mean (norm<0.995)")
tin = np.interp(-10.07, z_p, tm); tout = np.interp(10.07, z_p, tm)
ax.axvspan(tin, tout, color="0.9")
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("E_kin (eV)"); ax.legend(fontsize=8)
fig.savefig("fig_wp_ke_time.png", dpi=160)""", anchor="sorb-ke-calc")

md(r"""### KE vs position, and $S=-dE_\mathrm{kin}/dz$

Two extractions over the valid (norm>0.995) window:
(a) **endpoint**: vacuum KE plateau before vs after the slab, divided by
$L_\mathrm{slab}$ — immune to the in-well KE bump;
(b) **mid-slab gradient**: linear fit of KE($z$) for $z\in[-8,8]$ — rides on
top of the well bump, so treated as the cross-check.""", anchor="sorb-s")
code("""zq = np.interp(t, tm, z_p)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(zq[val], ke_ext[val])
ax.axvspan(-10.07, 10.07, color="0.9")
ax.set_xlabel("z of packet (Bohr)"); ax.set_ylabel("E_kin (eV)")
fig.savefig("fig_wp_ke_position.png", dpi=160)

pre = val & (zq > -15) & (zq < -11)          # vacuum, past injection transient
post = val & (zq > 12) & (zq < 25)           # vacuum, before CAP contact
ke_pre, ke_post = ke_ext[pre].mean(), ke_ext[post].mean()
S_end = (ke_pre - ke_post)/20.134
mid = val & (np.abs(zq) < 8)
S_grad = -np.polyfit(zq[mid], ke_ext[mid], 1)[0]
print(f"KE plateau before {ke_pre:.2f} eV, after {ke_post:.2f} eV "
      f"-> Delta {ke_pre-ke_post:.2f} eV")
print(f"S_orb (endpoint)  = {S_end:.2f} eV/Bohr")
print(f"S_orb (mid-slab gradient) = {S_grad:.2f} eV/Bohr  (well-bump caveat)")
print()
print(f"{'method':<28}{'S (eV/Bohr)':>12}")
for k, s in [("orbital-free E_pl/L", S_free), ("orbital endpoint", S_end),
             ("orbital mid-slab grad", S_grad), ("Bethe bulk point-charge", S_bethe)]:
    print(f"{k:<28}{s:>12.2f}")""", anchor="sorb-s-calc")

md("""## Takeaway

- **E_plateau = 4.1 eV** deposited in the slab (window means 4.57/4.18/4.12 eV
  over t∈[40,50)/[50,55)/[55,60), drift −0.006 eV/a.u. in the final window —
  converged; residual WP norm 9.7e-4 bounds any further drain at ~0.1 eV).
- **The artifact would have said 22 eV**: reported total − E_GS ends at
  0.81 Ha (22 eV), inflated 5× by the 17.8 eV norm-division kinetic artifact
  the in-run fix removes.
- **Identity EXACT in the interacting system**: max |Σocc·T/N − kinetic_INQ|
  = 0.0 Ha across all 1501 steps × 62 states.
- **Absorption**: WP 99.9% absorbed (survival 9.7e-4, vs ~8e-4 predicted for
  η=−1.0, W=15 at v=2.711 — the dt=0.04 CAP behaves as calibrated at dt=0.01,
  retiring most of the ungated-dt caveat empirically).
- **Cost**: OrbitalKineticStats 143 ms/step for 62 states at every-step
  cadence = 0.6% of the 23.7 s step — the fix is effectively free.
- **Stopping power**: orbital-free S = E_pl/L = **0.20 eV/Bohr**; the
  orbital-attributed −dKE/dz gives 0.07 eV/Bohr (both endpoint and mid-slab
  gradient) — the WP orbital's KE loss (1.5 eV) carries only ~1/3 of the
  total deposition, quantitatively confirming that quantum stopping must be
  measured from total energy deposition, not projectile KE. Bulk point-charge
  Bethe expects 0.73 eV/Bohr; the 3.6× suppression is consistent with the
  20-Bohr slab ≪ 79-Bohr wake wavelength + the packet form factor (inference).
- Open: no no-CAP twin (user scope cut) → no in-system conservation control;
  E_plateau is single-source but internally cross-checked.""",
   anchor="takeaway")

nb.build(CELLS, str(HERE / "extkin_plateau_E100_study.ipynb"))
