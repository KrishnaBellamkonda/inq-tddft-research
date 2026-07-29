#!/usr/bin/env python3
"""Build + EXECUTE the master study notebook for the graphene+CAP campaign.

Covers the full campaign journey up to now:
  Stage 1  ground state    -> validated gapless semimetal (24 C, 3x2)
  Stage 2  WP + CAP        -> 4 runs DELIVERED (centroid/channeling x cap/nocap)
  Stage 3  classical ens.  -> DEFERRED (log-grounded post-mortem)

House narrative (notebook-making skill): title+question -> conventions/symbols ->
setup+deviations -> source files -> results -> takeaway. Per the per-run
visual-intuition rule, EVERY significant run carries an xz density GIF + its
energetics BEFORE the aggregate comparison plots; an anomalous/deferred run gets
a post-mortem grounded in the actual log output (not a remembered claim).

Partial-tolerant: reads whatever runs exist. EXECUTED with nbconvert (0 errors);
the committed .ipynb carries its own outputs.

Usage:
    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_report.py
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
SYS = os.path.normpath(os.path.join(HERE, "..", ".."))   # systems/graphene
FIG = os.path.join(HERE, "figs"); os.makedirs(FIG, exist_ok=True)

nb = new_notebook()
C = nb.cells
def md(s): C.append(new_markdown_cell(s))
def code(s): C.append(new_code_cell(s))

# ============================================================================
# 1. TITLE + THE QUESTION
# ============================================================================
md(r"""# Graphene + complex absorbing potential — electron-scattering study notebook

**Feasibility replica of Yao & Schleife**, *wave-packet electron dynamics on
graphene*. A 100 eV electron is fired perpendicularly at a graphene sheet; a
two-sided complex absorbing potential (CAP) on the $z$-ends removes outgoing flux
so periodic images cannot re-enter. We ask:

> Does a CAP reliably absorb the scattered electron flux in a *graphene* TDDFT
> calculation (as it does in our jellium work), and how does a **quantum wave
> packet (WP)** compare with an **ensemble of classical point electrons** of the
> same Gaussian width, fired through an **atom (centroid)** vs a **hexagon hollow
> (channeling)**?

This is **methodology-faithful but NOT the paper's converged numbers** (reduced
cell / ensemble — see the deviations table). **All CAP results are PROVISIONAL**
until the `inq-study` engine regression (Task #7) validates the complexified
scalar potential.

| Where this sits | |
|---|---|
| `hypotheses/cap_scattering/` (this notebook) | what the graphene CAP runs MEAN |
| `scripts/gs/` | ground-state machinery (Stage 1) |
| `scripts/cap/` | WP+CAP build-once binary + dispatchers (Stage 2/3) |
| jellium `cap_*` sweeps | sibling CAP studies where the engine was first exercised |

Plan: `docs/plans/graphene-cap.md` · Handover: `docs/handovers/graphene-cap.md`.""")

# ============================================================================
# 2. CONVENTIONS + SYMBOLS (no formula dump)
# ============================================================================
md(r"""## 1. Conventions & symbols

Hartree atomic units throughout ($\hbar=m_e=e=1$); eV alongside
($1\,\mathrm{Ha}=27.2114\,\mathrm{eV}$). The WP and CAP construction follow the
INQ `perturbations::absorbing` API and our jellium CAP practice; the WP form
follows Yao & Schleife.

| symbol | meaning | value / range |
|---|---|---|
| $E$ | projectile kinetic energy | 100 eV |
| $k_0$ | injection wavevector, $z$ | $\sqrt{2E}=2.711\,a_0^{-1}$ |
| $d,\ \sigma_r$ | WP real-space width ($\sigma_r=d/\sqrt2$) | $d=1.1\,\text{Å}$, $\sigma_r=1.47\,a_0$ |
| $b=(c_x,c_y,z_0)$ | WP launch centre | $z_0=-12.65\,a_0$; $(c_x,c_y)$ per trajectory |
| $W=\eta$ | CAP strength (imaginary) | $-0.5\,\mathrm{Ha}$ |
| $L$ | total CAP width (each end) | $20\,a_0$ |
| $z_{in}$ | inner (free) half-region edge | $20\,a_0$ |
| $\varepsilon(t)$ | survival fraction in $|z|<z_{in}$ | dimensionless |
| $N_0$ | WP norm at $t=0$ | $1$ |

**Numerical requirement:** the CAP makes the Hamiltonian **non-Hermitian**, so
the propagator must be **ETRS** (enforced-time-reversal-symmetry), *not*
Crank–Nicolson (CN renormalises the state each step and would defeat absorption).""")

# ============================================================================
# 3. SETUP (reconstructable) + DEVIATIONS
# ============================================================================
md(r"""## 2. Simulation setup (reconstructable) & deviations from the paper

**Cell / geometry.** $3\times2$ rectangular graphene, **24 C atoms**, box
$13.9462\times16.1037\times60\,a_0$, orthorhombic **periodic**, $\Gamma$-only.
The zigzag multiplicity $n_x=3$ is chosen so the Dirac point at BZ **K folds onto
$\Gamma$** ($n_x$ a multiple of 3 ⇒ metallic; $n_x=4$ would open a spurious gap).
The sheet lies at $z=0$; the free region is $|z|<20\,a_0$, the CAP occupies
$20<|z|<30\,a_0$ at each end.

**Electronic structure.** XC **LDA/ALDA**, ONCV carbon pseudopotential, plane-wave
cutoff **50 Ha**, Fermi–Dirac smearing **0.10 eV** (needed for the gapless
semimetal). 96 electrons, 72 KS states (48 occ + 24 extra).

**Dynamics.** ETRS propagator, $dt=0.02$ a.u., $N_\text{STEPS}=1319$
($\tau\approx26.4$ a.u.), 8 LEED screens at $z=\pm4,\pm8,\pm12,\pm16\,a_0$,
density (total/system/WP) + WP-orbital VTIs at $\sim$63 frames.

**Trajectories.** centroid = atom at $(0,0)$; channeling = hexagon hollow at
$(4.6655,-2.684)\,a_0$.

**Initial state.** WP $\psi_{WP}$ injected into the lowest empty KS state,
orthogonalised against all occupied states (max overlap $\sim10^{-4}$, confirmed).
Classical: Gaussian-smeared $-1$ projectile of width $\sigma_r$, with an mt19937
ensemble (pos $\sigma=1.47$, mom $\sigma_k=0.481$).

| quantity | paper | this replica | why |
|---|---|---|---|
| supercell | 112 C | 24 C ($3\times2$) | one-GPU overnight; $n_x=3$ folds K→Γ |
| $z$-vacuum | 100 $a_0$ | 60 $a_0$ | grid × step-count cost |
| ensemble | $>100$ | 3 / type | overnight budget |
| classical width | 0.5 $a_0$ (spec) | 1.47 $a_0$ | match WP for a fair comparison |
| CAP engine | paper code | `inq-study` complexified | **PROVISIONAL** until Task #7 |""")

# ============================================================================
# 4. SOURCE FILES
# ============================================================================
md(r"""## 3. Source files (every artefact traceable)

| role | path (repo-relative `ResearchProject/systems/graphene/`) |
|---|---|
| GS config | `shared/configs/graphene_gs.hpp` |
| GS geometry generator | `scripts/gs/gen_geometry.py` → `shared/geometry/graphene_3x2.xyz` |
| GS run | `scripts/gs/run.cpp` → checkpoint `shared_gs/gs_3x2_50ha` |
| WP+CAP binary (build-once, env-driven) | `scripts/cap/run.cpp` → frozen `scripts/cap/run_wp` |
| WP dispatcher (4 runs) | `scripts/cap/dispatch_wp.sh` |
| classical binary / dispatcher (deferred) | `scripts/cap/run_classical.cpp`, `scripts/cap/dispatch_cl.sh` |
| per-run quicklook + email | `scripts/cap/post_and_email.py` |
| projectile pseudopotential | `shared/pseudopotentials/electron_gaussian_sigma1p47.upf` |
| VTI loader / slice conventions reused | `inqview.pipeline.density` |
| this builder | `hypotheses/cap_scattering/build_report.py` |
| per-run provenance | `cap_scattering/run_wp_*/results/run_summary.txt`, `.../raw/observables/*.csv`, `.../raw/vti/<cat>/` |""")

# ---- shared setup code cell (helpers, incl. xz density GIF) ----
code(f"""import os, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from IPython.display import Image, display
from inqview.visualisation import style
style.apply_theme()

SYS = {SYS!r}
FIG = os.path.join(SYS, "hypotheses", "cap_scattering", "figs")
HA  = 27.211386245988  # eV per Hartree
Z_IN = 20.0            # free-region half-width (a0); CAP for |z|>Z_IN

def read_summary(p):
    o = {{}}
    if os.path.exists(p):
        for ln in open(p):
            if "=" in ln:
                k, v = ln.split("=", 1); o[k.strip()] = v.strip()
    return o

def read_cols(p):
    if not os.path.exists(p):
        return {{}}
    import csv
    rows = list(csv.DictReader(open(p)))
    if not rows:
        return {{}}
    out = {{}}
    for k in rows[0]:
        try:
            out[k] = np.array([float(r[k]) for r in rows])
        except (ValueError, TypeError):
            pass
    return out

def _load_vti(path):
    # reuse inqview.pipeline.density conventions: VTK x-fastest -> (nx,ny,nz)
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(str(path)); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64)
    cube = flat.reshape((nz, ny, nx)).transpose(2, 1, 0)
    return cube, dict(origin=img.GetOrigin(), spacing=img.GetSpacing(), shape=(nx, ny, nz))

def make_xz_gif(cat_dir, out_gif, title="", stride=3, dpi=80, pct=99.7,
                cmap="inferno", fps=8):
    \"\"\"xz-plane (propagation plane) density GIF, fixed colour scale, CAP edges
    and sheet (z=0) marked. Returns out_gif or None if no frames.\"\"\"
    files = sorted(glob.glob(os.path.join(cat_dir, "*.vti")))[::stride]
    if not files:
        return None
    slabs, meta = [], None
    for f in files:
        cube, meta = _load_vti(f)
        slabs.append(cube.take(cube.shape[1] // 2, axis=1))   # mid-y -> (nx,nz)
    ox, oy, oz = meta["origin"]; dx, dy, dz = meta["spacing"]; nx, ny, nz = meta["shape"]
    flat = np.concatenate([s.ravel() for s in slabs])
    vmax = float(np.percentile(flat, pct)); vmin = max(0.0, float(np.percentile(flat, 100 - pct)))
    if vmax <= vmin:
        vmax = vmin + 1e-12
    extent = [ox, ox + nx * dx, oz, oz + nz * dz]   # x horizontal, z vertical
    tmp = out_gif + ".frames"; os.makedirs(tmp, exist_ok=True)
    pngs = []
    for i, (s, f) in enumerate(zip(slabs, files)):
        fig, ax = plt.subplots(figsize=(4.6, 3.6), dpi=dpi)
        im = ax.imshow(s.T, origin="lower", aspect="auto", cmap=cmap,
                       vmin=vmin, vmax=vmax, extent=extent)
        for zc in (Z_IN, -Z_IN):
            ax.axhline(zc, color="cyan", lw=0.7, ls="--", alpha=0.7)
        ax.axhline(0.0, color="white", lw=0.6, ls=":", alpha=0.6)  # graphene sheet
        step = int(os.path.basename(f).split("_t")[1].split(".")[0])
        ax.set(xlabel="x (bohr)", ylabel="z (bohr)", title=f"{{title}}  step {{step}}")
        fig.colorbar(im, ax=ax, shrink=0.85, label="density")
        fig.tight_layout()
        p = os.path.join(tmp, f"f{{i:04d}}.png"); fig.savefig(p); plt.close(fig)
        pngs.append(p)
    imageio.mimsave(out_gif, [imageio.imread(p) for p in pngs], fps=fps, loop=0)
    for p in pngs:
        os.remove(p)
    os.rmdir(tmp)
    return out_gif

def energetics_fig(obs_csv, title, out_png):
    \"\"\"Per-run energy components vs time (Delta from t=0).\"\"\"
    c = read_cols(obs_csv)
    fig, ax = plt.subplots(figsize=(7, 4))
    comps = [("energy_total", "ΔE_total"), ("energy_kinetic", "ΔE_kin"),
             ("energy_hartree", "ΔE_Hartree"), ("energy_xc", "ΔE_xc")]
    for k, lab in comps:
        if k in c and "time_au" in c:
            ax.plot(c["time_au"], c[k] - c[k][0], label=lab)
    ax.set(xlabel="time (a.u.)", ylabel="ΔE (Ha)", title=title)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out_png, dpi=120); plt.show()

print("paths + helpers ready (xz GIF + energetics); SYS =", SYS)""")

# ============================================================================
# 5. RESULTS
# ============================================================================
md("""## 4. Results

The campaign ran in three stages. Each *significant* run leads with its **xz
density GIF** (visual intuition first) and its **energetics**, before any
aggregate comparison; the deferred classical stage gets a log-grounded
post-mortem.""")

# ---- Stage 1: GS ----
md(r"""### 4.1 Stage 1 — ground state: a gapless Dirac semimetal

The replica is only meaningful if the static graphene is correct. Two checks:

1. **Energy per atom** — size-consistency of the cohesive energy:
$$ E_\text{atom} = \frac{E_\text{total}}{N_C}. $$""")

code(r"""gs = read_summary(os.path.join(SYS, "scripts", "gs", "results", "run_summary.txt"))
E_tot = float(gs["ground_state_energy_ha"])
N_C   = 24
print(f"E_total   = {E_tot:.4f} Ha")
print(f"E/atom    = {E_tot/N_C:.4f} Ha   (size-consistent cohesive energy)")
print(f"electrons = {gs['num_electrons']},  KS states = {gs['num_states']},  "
      f"xc = {gs['xc']},  cutoff = {gs['cutoff_ha']} Ha,  smearing = {gs['temperature_ev']} eV")""")

md(r"""2. **Gaplessness (Dirac point).** Graphene is a semimetal: valence and
conduction bands touch at $E_F$. With $\Gamma$-only sampling and $n_x=3$ folding,
this shows up as **partially-occupied states clustered at $E_F$** rather than a
clean HOMO/LUMO gap. We flag any state with fractional occupation
$0.05<f<1.95$ and measure their energy spread $\Delta E_\text{Dirac}$ — a true
semimetal has $\Delta E_\text{Dirac}\to0$ (a few meV at this cell size).""")

code(r"""eig = pd.read_csv(os.path.join(SYS, "scripts", "gs", "results",
                                "eigenvalues", "eigenvalues.csv"))
partial = eig[(eig.occupation > 0.05) & (eig.occupation < 1.95)]
spread_meV = (partial.eigenvalue_ev.max() - partial.eigenvalue_ev.min()) * 1000.0
print("Partially-occupied (Dirac-point) states:")
print(partial[["state_index", "eigenvalue_ev", "occupation"]].to_string(index=False))
print(f"\nDirac-point energy spread = {spread_meV:.1f} meV  ->  GAPLESS semimetal "
      f"({len(partial)} states at E_F)")
assert len(partial) >= 2 and spread_meV < 100.0, "expected a gapless cluster at E_F"

fig, ax = plt.subplots(figsize=(7.2, 4.3))
sc = ax.scatter(eig.state_index, eig.eigenvalue_ev, c=eig.occupation,
                cmap="viridis", s=22, edgecolor="none")
ax.axhline(partial.eigenvalue_ev.mean(), ls="--", lw=1, color="0.4",
           label=f"$E_F$ (Dirac, {len(partial)} states, {spread_meV:.0f} meV)")
ax.set(xlabel="state index", ylabel="eigenvalue (eV)",
       title="Graphene ground-state spectrum (colour = occupation)")
ax.legend(); ax.grid(alpha=0.3)
fig.colorbar(sc, label="occupation")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_gs_spectrum.png"), dpi=120)
plt.show()""")

md("""**Stage-1 verdict.** $E/\\text{atom}=-6.00$ Ha (size-consistent) and a
$\\sim$10 meV cluster of partially-occupied states at $E_F$ with no gap → the
static graphene is a correct gapless Dirac semimetal. The earlier $n_x=4$
(32-atom) attempt gave a spurious $\\sim$2 eV gap and was discarded.""")

# ---- Stage 2: WP + CAP ----
md(r"""### 4.2 Stage 2 — wave packet + CAP (DELIVERED, 4 runs)

The injected wave packet
$$ \psi_{WP}(\mathbf r) \propto \exp\!\Big(-\tfrac{(\mathbf r-\mathbf b)^2}{2d^2}\Big)\,
   e^{i k_0 z}, \qquad k_0=\sqrt{2E}, $$
is propagated under the two-sided sin² CAP
$$ V_\text{CAP}(z) = -iW\sin^2\!\Big(\tfrac{\pi (z-z_s)}{2\,d_z}\Big)\ \ \text{on each $z$-end}, $$
composed as `absorbing(η,+mid)+absorbing(η,−mid)` with $W=\eta=-0.5$ Ha.

Each of the four runs (centroid/channeling × CAP/no-CAP) is *significant* (a
distinct corner of the comparison), so each gets its own **xz density GIF of the
wave packet** (`density_rt_wp`, propagation plane; dashed cyan = CAP edges at
$z=\pm20$, dotted white = graphene at $z=0$) and its **energetics**, before the
aggregate survival comparison.""")

WP_RUNS = [
    ("centroid_nocap",   "centroid · no CAP"),
    ("centroid_cap",     "centroid · CAP on"),
    ("channeling_nocap", "channeling · no CAP"),
    ("channeling_cap",   "channeling · CAP on"),
]
for tag, nice in WP_RUNS:
    md(f"#### {nice}")
    code(f"""tag, nice = {tag!r}, {nice!r}
d = os.path.join(SYS, "cap_scattering", f"run_wp_{{tag}}")
s = read_summary(os.path.join(d, "results", "run_summary.txt"))
if s.get("run_completed") == "true":
    print(f"{{nice}}:  ε_survival = {{float(s['epsilon_survival']):.3f}}   "
          f"absorbed = {{float(s['absorbed_fraction']):.3f}}   "
          f"CAP = {{s['cap']}}   wall = {{float(s['wall_s'])/60:.0f}} min")
    # (1) visual intuition first: xz density GIF of the wave packet
    gif = make_xz_gif(os.path.join(d, "results/raw/vti/density_rt_wp"),
                      os.path.join(FIG, f"xz_wp_{{tag}}.gif"), title=nice)
    if gif:
        display(Image(filename=gif))
    # (2) energetics of this run
    energetics_fig(os.path.join(d, "results/raw/observables/observables.csv"),
                   f"Energetics — {{nice}}", os.path.join(FIG, f"energy_{{tag}}.png"))
else:
    print(f"{{nice}}: run not present / incomplete — skipped (partial-tolerant)")""")

md(r"""**Aggregate comparison — survival & absorption.** Figure of merit:
$$ \varepsilon(t)=\frac{1}{N_0}\int_{|z|<z_{in}}|\psi_{WP}|^2\,dV,\qquad
   \text{absorbed}(t)=1-\frac{\|\psi_{WP}(t)\|^2}{N_0}. $$
A working CAP drives $\varepsilon$ down and `absorbed` up; the no-CAP baseline
keeps $\|\psi_{WP}\|\approx1$ (flux exits $|z|<z_{in}$ but is never removed).""")

code(r"""runs = sorted(glob.glob(os.path.join(SYS, "cap_scattering", "run_wp_*")))
rows = []
for d in runs:
    s = read_summary(os.path.join(d, "results", "run_summary.txt"))
    if s.get("run_completed") != "true":
        continue
    rows.append(dict(run=os.path.basename(d).replace("run_wp_", ""), CAP=int(s["cap"]),
                     eps_survival=round(float(s["epsilon_survival"]), 4),
                     absorbed=round(float(s["absorbed_fraction"]), 4),
                     wall_min=round(float(s["wall_s"]) / 60.0, 1)))
wp = pd.DataFrame(rows).sort_values(["CAP", "run"]).reset_index(drop=True)
print("WP+CAP run inventory (PROVISIONAL — inq-study Task #7):")
display(wp)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for d in runs:
    s = read_summary(os.path.join(d, "results", "run_summary.txt"))
    if s.get("run_completed") != "true":
        continue
    name = os.path.basename(d).replace("run_wp_", "")
    c = read_cols(os.path.join(d, "results/raw/observables/inner_norm_vs_time.csv"))
    if "time_au" not in c:
        continue
    ls = "-" if int(s["cap"]) else "--"
    ax[0].plot(c["time_au"], c["survival_inner_over_N0"], ls, label=name)
    ax[1].plot(c["time_au"], c["absorbed_fraction"], ls, label=name)
ax[0].set(xlabel="time (a.u.)", ylabel=r"survival $\varepsilon(t)$",
          title=r"WP survival in free region $|z|<20\,a_0$")
ax[1].set(xlabel="time (a.u.)", ylabel="absorbed fraction", title="WP flux absorbed by CAP")
for a in ax:
    a.legend(fontsize=8); a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_survival.png"), dpi=120)
plt.show()""")

md("""**Stage-2 verdict (PROVISIONAL).** The per-run GIFs show the wave packet
launch at $z\\approx-12.65$, traverse the sheet, and — with the CAP on — be
extinguished as it reaches $z=\\pm20$; the no-CAP packets instead persist and wrap.
Quantitatively the CAP drives $\\varepsilon$ from $0.72$ to $0.12$, absorbing
$\\sim$85% of the outgoing flux; baselines absorb $\\sim$0 and conserve energy
(per-run energetics). Centroid $\\approx$ channeling at 100 eV. **The absorbing
boundary works in graphene.** LEED screens (8) and density / WP-wavefunction VTIs
are saved per run for downstream diffraction and current analysis.""")

# ---- Stage 3: classical (deferred) — log-grounded post-mortem ----
md(r"""### 4.3 Stage 3 — classical ensemble: bug DIAGNOSED + FIXED, runs LAUNCHED

The classical comparison (6 runs: centroid×3 + channeling×3) was blocked by an
energy anomaly. As of **2026-06-21 it is root-caused and fixed**; the ensemble is
running. Below we keep the grounded post-mortem of *what was wrong* (it is the
instructive part) and then state the fix and its validation.

**What actually happened, per the logs:**

1. **Frozen-ion bug (first attempt, `build_cl.log`).** INQ freezes ions by
   default. The 2-step smoke reported `e` constant to 10 digits
   ($-36.736273116 \to -36.736273164$) and **`KE_loss=-0.00 eV`** — the
   projectile never moved. Three production runs (`dispatch_cl.out`) finished
   `rc=0` but were physically null and were deleted. Fix: add `.ehrenfest()`.

2. **No "crash at step 24" — correction.** A prior handover note claimed the
   Ehrenfest run *crashed at step 24*. **The on-disk logs do not corroborate
   this.** The 30-step Ehrenfest smoke (`build_cl2.log`) ran to
   `step 30 … real-time propagation ended normally` and wrote
   `run_completed = true`. No error/abort/NaN line exists in any classical log.
   The record is corrected here.

3. **The real anomaly — energy reference + spurious vacuum loss.** Two grounded
   concerns remain, shown below from `results_smoke2`:
   - the classical run starts at $e\approx-36.74$ Ha vs the WP's $-139.8$ Ha
     (same GS, same C pseudo) — the projectile-**as-ion** (bare $-1$ + Gaussian
     UPF in a periodic cell) sits on a different energy reference than the
     projectile-**as-orbital**;
   - $e$ **rises** monotonically and the projectile loses **6.70 eV** of KE in
     just 30 steps while still $\sim$11 $a_0$ from the sheet (in vacuum) — too
     much loss, too far away.""")

code(r"""sm = os.path.join(SYS, "scripts", "cap_cl", "results_smoke2")
ssum = read_summary(os.path.join(sm, "run_summary.txt"))
if ssum:
    print("classical Ehrenfest smoke (results_smoke2):")
    print(f"  run_completed = {ssum.get('run_completed')}   "
          f"propagator = {ssum.get('propagator')}   N_STEPS = {ssum.get('N_STEPS')}")
    print(f"  proj r0 = {ssum.get('proj_r0')}  ->  rf = {ssum.get('proj_rf')}")
    print(f"  KE_final = {ssum.get('KE_final_eV')} eV   KE_loss = {ssum.get('KE_loss_eV')} eV")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
c = read_cols(os.path.join(sm, "raw/observables/observables.csv"))
if "energy_total" in c:
    ax[0].plot(c["time_au"], c["energy_total"], color="#b03030")
    ax[0].set(xlabel="time (a.u.)", ylabel=r"$e$ (Ha)",
              title=f"Electronic energy RISES (ref ≈ {c['energy_total'][0]:.1f} Ha, not −139.8)")
    ax[0].grid(alpha=0.3)
t = read_cols(os.path.join(sm, "raw/observables/electron_track.csv"))
if "vz" in t:
    ke = 0.5 * (t["vx"]**2 + t["vy"]**2 + t["vz"]**2) * HA
    ax[1].plot(t["time_au"], ke, color="#1a4ea0")
    ax[1].set(xlabel="time (a.u.)", ylabel="projectile KE (eV)",
              title="Projectile loses ~6.7 eV in vacuum (spurious)")
    ax[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_classical_postmortem.png"), dpi=120)
plt.show()""")

md(r"""**Visual intuition for the smoke** — xz total-density GIF of the classical
projectile (the bright $-1$ Gaussian moving up from $z\approx-12.65$). Useful to
*see* that the projectile is in vacuum, far from the sheet, while it is already
"losing" energy — supporting the spurious-loss diagnosis.""")

code(r"""sm = os.path.join(SYS, "scripts", "cap_cl", "results_smoke2")
gif = make_xz_gif(os.path.join(sm, "raw/vti/density_rt_total"),
                  os.path.join(FIG, "xz_classical_smoke.gif"),
                  title="classical smoke · total density", stride=1)
if gif:
    display(Image(filename=gif))
else:
    print("no classical VTI frames present")""")

md(r"""**ROOT CAUSE (confirmed at source + by smoke).** The projectile was inserted
as an ion with **`z_valence = 0`** but a $+1/r$ *repulsive* local potential. INQ's
ion–ion Ewald charge array is exactly `valence_charge()`
(`inq-study/src/ionic/interaction.hpp:329`), so a `z_valence=0` projectile is
**invisible to the ion–ion sum**: it felt the carbon **electrons** (repulsion, via
its local pp) but **not the carbon nuclei** (attraction). The two should cancel at
long range for neutral graphene; instead only the repulsion survived → the
projectile climbed a spurious potential and decelerated $\sim$6.7 eV while still
$\sim$11 $a_0$ out in vacuum. The same `z_valence=0` also made the local-potential
$G{=}0$/$\alpha$ term inconsistent → the constant **$+103$ Ha** offset
($-36.74$ vs $-143.94$). **Jellium was immune** (no nuclei; uniform medium → the
overlap is position-independent → zero force), which is why the identical UPF
worked there for months. The WP run is stable for a different reason: its $-1$
electron is a non-back-reacting probe orbital (occ 0), never an ion.

**FIX.** A corrected projectile UPF `electron_gaussian_sigma1p47_zm1.upf` with
**`z_valence = -1`** (now consistent with its $+1/r$ tail → the projectile is a
proper $-1$ charge in the Ewald sum and feels the nuclei), plus
**`.extra_electrons(+1)`** in `run.cpp` so the *quantum* electron count stays 96
(graphene); the cell then carries the physical net $-1$ with a $+1$ uniform
compensating background.

**VALIDATION (30-step smoke, `scripts/cap_cl/results_smoke_zm1/`).**
Number of electrons = 96 ✓; step-0 $e = -143.9047$ Ha (the $+103$ Ha offset is
**gone**) ✓; energy drift over 30 steps $= -1.4\times10^{-4}$ Ha ($-0.004$ eV) vs
the old $+6.7$ eV ✓; `proj_vf.z = 2.71111` vs `v0 = 2.71106`
(KE\_loss $= -0.004$ eV — **ballistic in vacuum**, no spurious drag) ✓. The
6-run ensemble (`run_cl_*`) is running; its KE-stopping curves auto-render below
once the tracks exist.""")

md(r"""**Classical KE stopping (auto-renders as `run_cl_*` complete).** Projectile
kinetic energy $KE(t)=\tfrac12 m(v_x^2+v_y^2+v_z^2)$ from `electron_track.csv`.
Stopping is measured from each member's **own** initial KE (the jittered $v_0$),
not the nominal 100 eV — fixing the earlier `KE_loss=-99` bookkeeping bug.
Centroid (solid) vs channeling (dashed); ensemble mean $\pm$ s.d. in the title.""")

code(r"""HA = 27.211386245988
cls = sorted(glob.glob(os.path.join(SYS, "cap_scattering", "run_cl_*")))
rows = []
fig, ax = plt.subplots(figsize=(7.4, 4.4))
for d in cls:
    tk = read_cols(os.path.join(d, "results/raw/observables/electron_track.csv"))
    if not tk or "vz" not in tk or len(tk["vz"]) < 3:
        continue
    ke = 0.5 * (tk["vx"]**2 + tk["vy"]**2 + tk["vz"]**2) * HA
    name = os.path.basename(d).replace("run_cl_", "")
    is_cent = "centroid" in name
    ax.plot(tk["time_au"], ke, lw=1.3, ls="-" if is_cent else "--",
            color="#1a4ea0" if is_cent else "#b03030", alpha=0.8, label=name)
    rows.append((name, "centroid" if is_cent else "channeling",
                 float(ke[0]), float(ke[-1]), float(ke[0] - ke[-1])))
if rows:
    import statistics as _st
    losses = [r[4] for r in rows]
    mu = _st.fmean(losses); sd = _st.pstdev(losses) if len(losses) > 1 else 0.0
    ax.set(xlabel="time (a.u.)", ylabel="projectile KE (eV)",
           title=f"Classical KE stopping — {len(rows)} run(s); "
                 f"ΔKE = {mu:.2f} ± {sd:.2f} eV")
    ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_classical_stopping.png"), dpi=130)
    plt.show()
    print(f"{'member':<18}{'kind':<12}{'KE0(eV)':>9}{'KEf(eV)':>9}{'ΔKE(eV)':>9}")
    for n, k, k0, kf, dl in rows:
        print(f"{n:<18}{k:<12}{k0:>9.2f}{kf:>9.2f}{dl:>9.2f}")
else:
    plt.close(fig)
    print("no completed run_cl_* tracks yet — re-run this builder after the ensemble finishes")""")

# ============================================================================
# 4.4 FIELD-LEVEL DIAGNOSTICS — planar Δn(z,t) + LEED diffraction
# ============================================================================
md(r"""## 4.4 Field-level diagnostics: planar $\Delta n(z,t)$ and LEED diffraction

Two field observables harvested from data already on disk (no extra runs), using
the validated `inqview.analysis` kernels:

- **Planar-integrated $\Delta n(z,t)$** (the paper's Fig. 1 trace):
  $$\Delta n(z,t)=\iint\!\big[n(x,y,z,t)-n(x,y,z,t_0)\big]\,dx\,dy
   =(dx\,dy)\!\sum_{x,y}\big[n(\cdot,t)-n(\cdot,t_0)\big],$$
  the transverse-integrated density change vs propagation coordinate $z$ and time.
  $\Delta n>0$ = pile-up, $\Delta n<0$ = depletion. Kernel:
  `inqview.analysis.planar_delta_map`.
- **Kinematic LEED pattern** on each of the 8 $z$-screens: the squared 2-D Fourier
  transform (Hann-windowed) of the time-integrated real-space screen density,
  $I(k_x,k_y)=|\,\mathrm{FFT2}[\bar n(x,y)\,w]\,|^2$ — forward ($+z$) screens =
  *transmission*, backward ($-z$) = *reflection*. Kernel:
  `inqview.analysis.diffraction_pattern`. This is a single-scattering replica
  diagnostic (Van Hove–Weinberg–Chan kinematic limit), **not** a dynamical LEED
  calculation.""")

code(r"""from inqview.analysis.planar_density import planar_delta_map
from inqview.analysis.diffraction import diffraction_pattern
from inqview.io.leed import load_leed_pattern

def _load_vti_series(cat_dir, dt=0.02, stride=1):
    # Load a density_rt_* category dir as (cubes, times, z) ordered by step.
    files = sorted(glob.glob(os.path.join(cat_dir, "*.vti")),
                   key=lambda f: int(os.path.basename(f).split("_t")[1].split(".")[0]))[::stride]
    if not files:
        return None
    cubes, times, meta = [], [], None
    for f in files:
        cube, meta = _load_vti(f)
        cubes.append(cube)
        times.append(int(os.path.basename(f).split("_t")[1].split(".")[0]) * dt)
    ox, oy, oz = meta["origin"]; dx, dy, dz = meta["spacing"]; nx, ny, nz = meta["shape"]
    z = oz + np.arange(nz) * dz
    return cubes, np.array(times), z, dx * dy

def planar_dn_fig(run_dir, category, title, out_png, stride=1):
    cat = os.path.join(run_dir, "results/raw/vti", category)
    ser = _load_vti_series(cat, stride=stride)
    if ser is None:
        print(f"no VTI frames in {category} for {os.path.basename(run_dir)}"); return None
    cubes, times, z, area = ser
    dmap = planar_delta_map(cubes, times, z, axis=2, cell_area=area)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    vmax = np.percentile(np.abs(dmap.dn), 99.5) or 1e-12
    im = ax.imshow(dmap.dn, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, extent=dmap.extent)
    ax.axhline(0.0, color="k", lw=0.7, ls=":")                 # graphene sheet
    for zc in (Z_IN, -Z_IN):
        ax.axhline(zc, color="0.4", lw=0.7, ls="--")           # CAP edges
    ax.set(xlabel="time (a.u.)", ylabel="z (bohr)", title=title)
    fig.colorbar(im, ax=ax, label=r"$\Delta n(z,t)$  (electrons / bohr)")
    fig.tight_layout(); fig.savefig(out_png, dpi=130); plt.show(); return out_png

# Prefer a completed classical run (real back-reaction); else the WP CAP run.
cl = sorted(glob.glob(os.path.join(SYS, "cap_scattering", "run_cl_*")))
cl = [d for d in cl if glob.glob(os.path.join(d, "results/raw/vti/density_rt_system/*.vti"))]
if cl:
    _pdrun, _cat, _lbl = cl[0], "density_rt_system", "classical projectile (system back-reaction)"
else:
    _pdrun = os.path.join(SYS, "cap_scattering", "run_wp_centroid_cap")
    _cat, _lbl = "density_rt_total", "WP centroid+CAP (total density)"
planar_dn_fig(_pdrun, _cat,
              f"planar Δn(z,t) — {_lbl}\n(dashed=CAP edges z=±20, dotted=sheet z=0)",
              os.path.join(FIG, "fig_planar_dn.png"))""")

md(r"""**LEED diffraction across the 8 screens.** Each panel is $|\mathrm{FFT2}|^2$
of the time-integrated density crossing that $z$-plane (log colour, DC removed,
Hann-windowed). The top row is *transmission* ($z>0$), the bottom *reflection*
($z<0$); columns are increasing $|z|$. Intensity localised near the centre = mostly
specular/forward; structure away from the centre on the periodic reciprocal-lattice
positions = diffraction off the graphene lattice.""")

code(r"""SCR_Z = [4, 8, 12, 16]
def leed_run(run_dir, title, out_png):
    sdir = os.path.join(run_dir, "results/raw/screens/total")
    if not os.path.isdir(sdir):
        print(f"no screens for {os.path.basename(run_dir)}"); return None
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.4))
    any_panel = False
    for col, zc in enumerate(SCR_Z):
        for row, sgn in enumerate(("+", "-")):       # row0 = +z transmission
            ax = axes[row, col]
            f = os.path.join(sdir, f"scr{sgn}{zc:02d}.dat")
            if not os.path.exists(f):
                ax.axis("off"); continue
            pat = load_leed_pattern(f)
            d = diffraction_pattern(pat.data, pat.dx_bohr, pat.dy_bohr, hann=True)
            I = d.intensity; floor = (I[I > 0].min() if np.any(I > 0) else 1e-12)
            ax.imshow(np.log10(I + floor), origin="lower", aspect="equal",
                      cmap="magma", extent=d.extent)
            ax.set_title(f"{'transmission' if sgn=='+' else 'reflection'} z={sgn}{zc}",
                         fontsize=9)
            ax.set_xlabel(r"$k_x$ (rad/$a_0$)", fontsize=7)
            ax.set_ylabel(r"$k_y$ (rad/$a_0$)", fontsize=7)
            ax.tick_params(labelsize=6); any_panel = True
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(); fig.savefig(out_png, dpi=130)
    if any_panel: plt.show()
    else: plt.close(fig)
    return out_png

leed_run(os.path.join(SYS, "cap_scattering", "run_wp_centroid_cap"),
         "Kinematic LEED — WP centroid + CAP (8 screens, log |FFT2|²)",
         os.path.join(FIG, "fig_leed_centroid_cap.png"))""")

# ============================================================================
# 6. TAKEAWAY
# ============================================================================
md(r"""## 5. Takeaway (PROVISIONAL)

- **Ground state ✓** — 24-C $3\times2$ graphene is a correct **gapless Dirac
  semimetal** ($E/\text{atom}=-6.00$ Ha, $\sim$10 meV cluster at $E_F$, no gap).
  The $n_x=4$ trap (spurious 2 eV gap) was caught and avoided via $n_x=3$ K→Γ folding.
- **CAP works in graphene ✓ (PROVISIONAL)** — two-sided sin² CAP drives the WP
  survival $\varepsilon$ from $0.72$ (no CAP) to $0.12$, absorbing $\sim$85% of the
  outgoing flux; per-run GIFs show the packet being extinguished at $z=\pm20$.
  Centroid $\approx$ channeling at 100 eV.
- **Classical ensemble — FIXED, running ✓** — root cause was a `z_valence=0`
  projectile (invisible to the ion–ion Ewald → felt carbon electrons but not
  nuclei → spurious vacuum drag + $+103$ Ha offset). Fix: `z_valence=-1` UPF +
  `extra_electrons(+1)`. Smoke validated (step-0 $e=-143.90$, ballistic vacuum);
  the 6-run ensemble is running and its KE-stopping curves render in §4.3.
- **Field diagnostics ✓** — planar $\Delta n(z,t)$ and kinematic LEED (8 screens)
  harvested from on-disk data via tested `inqview.analysis` kernels (§4.4).
- **Provisionality gate.** Every CAP number here is **PROVISIONAL until the
  `inq-study` engine regression (Task #7)**. This is a *feasibility replica* —
  reduced cell/ensemble, not the paper's converged values. Do not quote as physical.""")

# ============================================================================
# EXECUTE + WRITE
# ============================================================================
ep = ExecutePreprocessor(timeout=2400, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": HERE}})
out = os.path.join(HERE, "cap_scattering_study.ipynb")
nbf.write(nb, out)
print(f"wrote + executed {out}")
