#!/usr/bin/env python3
"""Build + execute the CAP-in-jellium Baseline-1 study notebook (ADR 0007).

Baseline 1 = the two-sided sin^2 CAP switched on over an equilibrium jellium bath
with NO projectile: the bath-drainage reference for the later B2 (classical) and
B3 (WP) wake runs. This builder reads the run provenance (run_summary.txt,
electron_number.csv) + the region-resolved analysis (cap_b1_region_drainage.csv,
from build_b1_drainage.py) + the depth-sweep (eta_compare) and renders the
house-narrative notebook. Partial-tolerant: builds from whatever runs exist.

    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_cap_baselines_report.py

All absorption numbers PROVISIONAL until the inq-study engine regression (Task #7).
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
           "scripts/cap_baselines/results")

nb = new_notebook()
C = nb.cells


def md(s):
    C.append(new_markdown_cell(s))


def code(s):
    C.append(new_code_cell(s))


# ----------------------------------------------------------------- 1. title + Q
md(r"""# Baseline 1 — does a two-sided CAP leave a usable jellium bath?

**System.** Cubic 50³ Bohr periodic jellium, $N=162$ electrons ($r_s=5.69$),
LDA, with a **two-sided complex absorbing potential (CAP)** switched on at $t=0$.

**The question.** Our jellium stopping-power runs are about to gain CAPs so the
projectile and its wake can *leave* the box instead of wrapping around. But a CAP
is non-Hermitian: it removes density wherever $|\psi|^2\neq 0$ — and in jellium the
equilibrium gas **permanently fills** the absorbing slabs. Before trusting any
wake/stopping result we must answer: **how much, and how fast, does the CAP drain
the equilibrium bath?** Baseline 1 (CAP on, *no projectile*) isolates exactly this,
and becomes the subtraction reference for B2 (classical) and B3 (wavepacket).

> **PROVISIONAL.** This is the first use of the inq-study built-in CAP (a
> complexified scalar potential) in an *interacting* bath. Every absorption number
> here stays provisional until the inq-study engine regression (Task #7) passes.

| Baseline | projectile | role |
|---|---|---|
| B0 | — | the Hermitian ground state (CAP-free) = the $t=0$ reference |
| **B1 (this)** | **none** | **CAP-on-bath drainage reference** |
| B2 | classical σ=0.5 e⁻, 100 eV | wake = B2 − B1 |
| B3 | σ=0.5 Gaussian WP, 100 eV | wake = B3 − B1 |
""")

# --------------------------------------------------------- 2. conventions + syms
md(r"""## Conventions & symbols

**Units: Hartree atomic units** ($\hbar=m_e=e=1$, $4\pi\varepsilon_0=1$); lengths
in Bohr, energies in Hartree (1 Ha = 27.2114 eV), time in a.u.

| symbol | meaning | value / range |
|---|---|---|
| $L$ | cubic box side | 50 Bohr |
| $dx$ | grid spacing | 0.40 Bohr |
| $N$ | bath electrons | 162 (closed shell) |
| $r_s$ | Wigner–Seitz radius | 5.69 Bohr |
| $\eta$ | CAP depth | $-0.5$, $-0.10$ Ha (swept $-0.05\!\to\!-0.5$) |
| $w$ | CAP width per side | 10 Bohr (slabs $|z|\in[15,25]$) |
| free region | absorber-free zone | $|z|<15$ Bohr |
| $dt$ | time step | 0.02 a.u. |
| $\tau$ | total propagated time | 140 a.u. (7000 steps) |

**Method notes.** The CAP is added only in real time (the ground state is the
plain Hermitian jellium — verified CAP-independent). Because the CAP is
non-Hermitian the propagator **must be ETRS, not Crank–Nicolson** (CN renormalises
the orbitals each step and would silently undo the absorption). The two-sided CAP
follows De Giovannini, Larsen & Rubio (2014): $V(z)=i\,\eta\,\sin^2[\,]$ in each
slab.""")

# --- formula: CAP, at point of use
md(r"""### The CAP potential

$$V(z) = i\,\eta\,\sin^2\!\Big(\frac{(z-z_0)\,\pi}{2w}\Big),\quad z\in\text{slab};
\qquad \text{two-sided} = V_{+}\!\oplus V_{-}.$$

In INQ fractional (centred) coordinates the slabs are at $|z|>0.3$ (i.e. $|z|>15$
Bohr): `absorbing(eta, +0.4, 0.2) + absorbing(eta, -0.4, 0.2)`.""")

# --- jellium scales: one quantity per cell, dependency order
md(r"""### Bath scales — set the drainage clock

The CAP drains the slabs fast, then the free-region gas flows in at the Fermi
velocity. So the relevant clock for *free-region* survival is the Fermi transit
time. We build it up from $r_s$.

Density $n=N/L^3$, then $r_s=(3/4\pi n)^{1/3}$.""")
code(r"""import numpy as np
N, L = 162, 50.0
n = N / L**3
rs = (3.0 / (4*np.pi*n))**(1/3)
print(f"n   = {n:.6f} e/Bohr^3")
print(f"r_s = {rs:.3f} Bohr")""")

md(r"""Fermi wavevector $k_F=(3\pi^2 n)^{1/3}$.""")
code(r"""kF = (3*np.pi**2 * n)**(1/3)
print(f"k_F = {kF:.4f} 1/Bohr")""")

md(r"""Fermi velocity $v_F=k_F/m_e=k_F$ (a.u.).""")
code(r"""vF = kF
print(f"v_F = {vF:.4f} a.u.")""")

md(r"""Plasmon energy $\omega_p=\sqrt{4\pi n}$ (the bath's collective scale).""")
code(r"""wp = np.sqrt(4*np.pi*n)
print(f"omega_p = {wp:.4f} Ha = {wp*27.2114:.2f} eV")""")

md(r"""Fermi transit time of the free half-width (15 Bohr): $t_F=15/v_F$.
This is the timescale on which bath inflow refills/depletes the free region.""")
code(r"""tF = 15.0 / vF
print(f"t_F (15 Bohr / v_F) = {tF:.1f} a.u.")
print(f"projectile (v=2.71) crosses the 30-Bohr free region in {30/2.711:.1f} a.u.")""")

# --------------------------------------------------------------- 3. setup table
md(r"""## Simulation setup (fully reconstructable)

- **Cell:** cubic $L=50$ Bohr, periodic; grid $dx=0.40$ ⇒ $125^3$ points,
  cutoff $\pi^2/2dx^2=30.84$ Ha.
- **Electrons:** $N=162$ (closed $|G|^2=6$ shell), LDA, $T=100$ K smearing,
  extra states per `Common_E100_L50_cubic`. Ground state reused from
  `checkpoints/gs_L50_cubic_N162_dx0p40` (pure jellium — no projectile in the GS).
- **Dynamics:** ETRS, $dt=0.02$ a.u., $N_{\rm STEPS}=7000$ ⇒ $\tau=140$ a.u.
  Density VTI every 23 steps (~300 frames); $N(t)$ every step.
- **CAP:** two-sided sin², $w=10$ Bohr/side, slabs $|z|\in[15,25]$, free $|z|<15$.
  Two runs: $\eta=-0.5$ (the chosen production config) and $\eta=-0.10$ (gentler).
- **No projectile** (this is B1).""")

# --------------------------------------------------------------- 4. source files
md(r"""## Source files

| file | role |
|---|---|
| `ResearchProject/systems/jellium/scripts/cap_baselines/run.cpp` | the b1/b2/b3 run (built vs inq-study) |
| `ResearchProject/systems/jellium/scripts/cap_baselines/run_b1_launch.sh` | launcher (emails on completion) |
| `ResearchProject/systems/jellium/hypotheses/cap_baselines/build_b1_drainage.py` | region-N + carpet analysis |
| `ResearchProject/systems/jellium/hypotheses/cap_baselines/precompute_b2.py` | B2 GIFs (density/wake/E_z) + projectile overlay |
| `ResearchProject/systems/jellium/hypotheses/cap_baselines/precompute_b3.py` | B3 GIFs + WP centroid |
| `inq-stack/python/inqview/analysis/efield.py` | E-field kernel (FFT-Poisson) |
| `inq/src/perturbations/absorbing.hpp` (inq-study) | the sin² CAP |
| `…/scripts/cap_baselines/results/b{1_eta0p50,1_eta0p10,2_classical_E100,3_wp_E100}/run_summary.txt` | per-run provenance |
| this builder | `hypotheses/cap_baselines/build_cap_baselines_report.py` |

Method source: **De Giovannini, Larsen & Rubio, *Eur. Phys. J. B* 88, 56 (2014)**
(sin² CAP / absorbing boundaries).""")

# --------------------------------------------------------------- 5. results
md(r"""## Result 1 — drainage vs CAP depth (early transient, 2 a.u.)

A short depth sweep (100 steps = 2.0 a.u.) of the *total* absorbed fraction. Even
the gentlest depth drains a few percent of the whole bath immediately — because
the slabs (40% of the box volume) are full of equilibrium gas. Sub-linear in
$|\eta|$.""")
code(r"""from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from inqview.visualisation import style
style.apply_theme()

RES = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
           "scripts/cap_baselines/results")

def kv(p):
    d = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition('=')
        d[k.strip()] = v.strip()
    return d

etas, abss = [], []
for sub in sorted((RES/'eta_compare').glob('eta_*')):
    f = sub/'run_summary.txt'
    if f.exists():
        d = kv(f); etas.append(abs(float(d['cap_eta_ha']))); abss.append(float(d['absorbed_frac']))
# include the eta=-0.5 100-step pilot if present
pil = RES/'pilot_b1'/'run_summary.txt'
if pil.exists():
    d = kv(pil); etas.append(abs(float(d['cap_eta_ha']))); abss.append(float(d['absorbed_frac']))
order = np.argsort(etas); etas = np.array(etas)[order]; abss = np.array(abss)[order]
fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot(etas, 100*abss, 'o-')
ax.set_xlabel(r'CAP depth $|\eta|$ (Ha)'); ax.set_ylabel('bath absorbed in 2 a.u. (%)')
ax.set_title('Total drainage vs CAP depth (B1, 100 steps)')
for e, a in zip(etas, abss):
    ax.annotate(f'{100*a:.1f}%', (e, 100*a), textcoords='offset points', xytext=(4, 4), fontsize=8)
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b1_drainage_vs_eta.png', dpi=150); plt.show()
print(dict(zip([f'-{e:.2f}' for e in etas], [f'{100*a:.1f}%' for a in abss])))""")

md(r"""## Result 2 — total electron number $N(t)$ over the full run

$N(t)=\int n(\mathbf r,t)\,d^3r$ (every step). Both depths drive the *total* bath
to near-zero by $\tau=140$ a.u. — a relentless, near-exponential collapse with no
plateau. (This is the alarming-looking number; Result 3 shows it is mostly the
slabs + late time.)""")
code(r"""import csv
fig, ax = plt.subplots(figsize=(7.0, 4.2))
summ = {}
for sub, lab, c in [('b1_eta0p50', r'$\eta=-0.5$', 'C3'), ('b1_eta0p10', r'$\eta=-0.10$', 'C0')]:
    f = RES/sub/'raw/observables/electron_number.csv'
    if not f.exists(): continue
    t, Nt = [], []
    with open(f) as fh:
        r = csv.DictReader(fh)
        for row in r:
            t.append(float(row['time_au'])); Nt.append(float(row['N_total']))
    ax.plot(t, Nt, color=c, lw=2, label=lab)
    s = kv(RES/sub/'run_summary.txt'); summ[sub] = s
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('total bath electrons N(t)')
ax.set_title('B1: total bath collapse under the CAP'); ax.legend()
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b1_total_N.png', dpi=150); plt.show()
for k, v in summ.items():
    print(f"{k}: N0={v.get('N0')} N_final={v.get('N_final')} absorbed={v.get('absorbed_frac')} nan={v.get('nan_seen')}")""")

md(r"""## Result 3 — region-resolved survival: the FREE region holds through the wake window

Split $N(t)$ into the **free region** $|z|<15$ (where the wake is measured) and the
**CAP slabs** $|z|\ge15$ (the absorber). The free region holds **97.2 e⁻** at
$t=0$; the table below is its survival. The projectile crosses the 30-Bohr free
region in ~10 a.u., over which the free bath is ~90% intact at $\eta=-0.5$ and
~94% at $\eta=-0.10$. *That* is the usable window.""")
code(r"""import csv
rows = {}
csvf = HERE/'cap_b1_region_drainage.csv'
with open(csvf) as fh:
    for row in csv.DictReader(fh):
        rows.setdefault(row['run'], []).append((float(row['t_au']), float(row['N_free']),
                                                 float(row['N_slab']), float(row['N_total'])))
fig, ax = plt.subplots(figsize=(7.0, 4.2))
for sub, lab, c in [('b1_eta0p50', r'$\eta=-0.5$', 'C3'), ('b1_eta0p10', r'$\eta=-0.10$', 'C0')]:
    if sub not in rows: continue
    a = np.array(rows[sub])
    ax.plot(a[:,0], a[:,1], color=c, lw=2, label=f'{lab}: free |z|<15')
    ax.plot(a[:,0], a[:,2], color=c, lw=1, ls='--', label=f'{lab}: slab |z|>=15')
ax.axvspan(0, 12, color='0.6', alpha=0.15)
ax.set_xlim(0, 60)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('bath electrons in region')
ax.set_title('B1 region-resolved: free region survives the wake window (shaded)')
ax.legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b1_region_drainage_nb.png', dpi=150); plt.show()

print('free-region survival (% of t=0):')
for sub in ('b1_eta0p50', 'b1_eta0p10'):
    if sub not in rows: continue
    a = np.array(rows[sub]); N0f = a[0,1]
    msg = '  '.join(f't={tt}: {100*a[np.argmin(np.abs(a[:,0]-tt)),1]/N0f:.0f}%' for tt in (5,10,15,20,30))
    print(f'  {sub}: {msg}')""")

md(r"""## Result 4 — density carpet $n(z,t)$ ($\eta=-0.5$)

Transverse-integrated linear density vs $(z,t)$: the bath visibly collapses inward
from the CAP edges (dashed). Generated by `build_b1_drainage.py`.""")
code(r"""from IPython.display import Image, display
carpet = HERE/'fig_b1_density_carpet.png'
display(Image(filename=str(carpet))) if carpet.exists() else print('carpet not built yet')""")

md(r"""## Result 5 — total electronic density GIF (xz slice, $\eta=-0.5$)

The animated **xz slice through the cell mid-plane** (mid-$y$): $x$ horizontal,
$z$ vertical, colour = total electronic density $n(x,z;y{=}0,t)$, frame = time.
With no projectile this *is* the bath; you watch it drain from the CAP slabs
(dashed cyan, $|z|=15$) inward. Frames strided to ≤120.""")
code(r"""import re
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.animation as animation

def _vol(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(str(path)); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    o, s = img.GetOrigin(), img.GetSpacing()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return o, s, a

gif_out = HERE/'fig_b1_density_xz.gif'
vdir = RES/'b1_eta0p50'/'raw/vti/density_system'
frames = sorted(vdir.glob('density_t*.vti'),
                key=lambda p: int(re.search(r'_t(\d+)', p.name).group(1)))
if len(frames) > 120:
    frames = frames[::len(frames)//120 + 1]
if frames:
    o, s, a0 = _vol(frames[0]); nz, ny, nx = a0.shape; midy = ny//2
    x = o[0] + s[0]*np.arange(nx); z = o[2] + s[2]*np.arange(nz)
    slabs = [_vol(p)[2][:, midy, :] for p in frames]          # each (nz, nx)
    steps = [int(re.search(r'_t(\d+)', p.name).group(1)) for p in frames]
    # SHARED density colour scale (same as the B3 density GIFs) — fixed clim.
    import json
    _clf = HERE/'cap_b3_clim.json'
    vmin, vmax = (json.load(open(_clf))['density'] if _clf.exists()
                  else (0.0, max(float(sl.max()) for sl in slabs)))
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    im = ax.imshow(slabs[0], origin='lower', aspect='auto',
                   extent=[x[0], x[-1], z[0], z[-1]], vmin=vmin, vmax=vmax, cmap='inferno')
    for ze in (-15, 15):
        ax.axhline(ze, color='cyan', ls='--', lw=1)
    ax.set_xlabel('x (Bohr)'); ax.set_ylabel('z (Bohr)')
    cb = fig.colorbar(im, ax=ax, label=r'$n$  (e/Bohr$^3$)')
    ttl = ax.set_title('')
    def upd(i):
        im.set_data(slabs[i]); ttl.set_text(f'B1 $\\eta=-0.5$   t = {steps[i]*0.02:6.1f} a.u.')
        return im, ttl
    animation.FuncAnimation(fig, upd, frames=len(slabs), interval=120, blit=False
                            ).save(gif_out, writer='pillow', dpi=90)
    plt.close(fig)
    print(f'wrote {gif_out}  ({len(slabs)} frames)')
else:
    print('no density_system frames found')""")
md(r"""![B1 density xz](fig_b1_density_xz.gif)""")

# ===================== Baseline 2 — the classical projectile =================
md(r"""# Baseline 2 — the classical σ=0.5 electron projectile

A **classical** $\sigma=0.5$ Bohr electron (a soft Gaussian pseudo-ion, mass
$=m_e$, charge $-1$), launched at $z_0=-13$ with $v_0=2.711$ a.u. (100 eV) moving
$+z$ under **Ehrenfest** dynamics, with the **same CAP ($\eta=-0.5$) and window as
B1** so B1 is the exact subtraction reference. This is the classical-projectile
counterpart of B3 (the quantum wavepacket).

> **Key subtlety — the CAP does *not* absorb the classical projectile.** The CAP is
> a potential on the *electron density*; the projectile lives in `ions`, so it is
> never absorbed. It coasts through the **periodic** box and decelerates against the
> bath: its final position is $z=346$ Bohr — about **7 box traversals**. So "the
> projectile exits through the far CAP" (true for the B3 wavepacket, which *is*
> density) does **not** hold here. The clean stopping measurement is therefore the
> **first traversal** — the transit window $t\le t^\star$ — exactly as for B3.

**Reference time** $t^\star=(15-(-13))/v_0$: the projectile reaches the far edge of
the free region ($z=+15$) at its initial speed. Drawn as a dashed line throughout.""")
code(r"""import csv, numpy as np
# shared loader + reference time (defined here; B3 reuses it)
def load_obs(sub):
    t, E, kin, har, xc, jz = [], [], [], [], [], []
    with open(RES/sub/'raw/observables/observables.csv') as fh:
        for r in csv.DictReader(fh):
            t.append(float(r['time_au'])); E.append(float(r['energy_total']))
            kin.append(float(r['energy_kinetic'])); har.append(float(r['energy_hartree']))
            xc.append(float(r['energy_xc'])); jz.append(float(r['current_z']))
    return dict(t=np.array(t), E=np.array(E), kin=np.array(kin),
               har=np.array(har), xc=np.array(xc), jz=np.array(jz))
B1 = load_obs('b1_eta0p50')
k0 = 2.711063; z0, zedge = -13.0, 15.0
t_star = (zedge - z0)/k0

B2 = load_obs('b2_classical_E100')
# classical projectile track (dedup the duplicate step-0 row)
tk, zk, vzk, vk = [], [], [], []
seen = set()
with open(RES/'b2_classical_E100/raw/observables/electron_track.csv') as fh:
    for r in csv.DictReader(fh):
        st = int(r['step'])
        if st in seen: continue
        seen.add(st)
        tk.append(float(r['time_au'])); zk.append(float(r['z']))
        vx, vy, vz = float(r['vx']), float(r['vy']), float(r['vz'])
        vzk.append(vz); vk.append((vx*vx + vy*vy + vz*vz)**0.5)
tk = np.array(tk); zk = np.array(zk); vzk = np.array(vzk); vk = np.array(vk)
KE = 0.5*vk**2            # projectile KE (Ha); m_e = 1 a.u.
HA = 27.2114
print(f't_star = {t_star:.2f} a.u.   (projectile reaches z=+15)')
print(f'projectile: z0={zk[0]:.1f} -> z_final={zk[-1]:.1f} Bohr '
      f'({zk[-1]/50:.1f} box traversals)')
print(f'v_z: {vzk[0]:.4f} -> {vzk[-1]:.4f} a.u.   '
      f'KE: {KE[0]:.4f} -> {KE[-1]:.4f} Ha  (lost {(KE[0]-KE[-1])*HA:.1f} eV total)')""")

md(r"""## Result B2.1 — projectile stopping: $v_z(t)$ and $KE(t)$ *(the headline)*

The projectile **decelerates** as it ploughs through the bath. $KE=\tfrac12 m_e v^2$
(a.u., $m_e=1$). The drop up to $t^\star$ is the energy lost over the clean first
traversal (28 Bohr); after that it keeps losing more slowly as it re-crosses the
draining periodic bath.""")
code(r"""fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8))
axs[0].plot(tk, vzk, lw=2, color='C3'); axs[0].axvline(t_star, ls='--', color='0.4')
axs[0].set_ylabel(r'$v_z$ (a.u.)'); axs[0].set_title('projectile speed')
axs[1].plot(tk, KE*HA, lw=2, color='C1'); axs[1].axvline(t_star, ls='--', color='0.4')
axs[1].set_ylabel('KE (eV)'); axs[1].set_title('projectile kinetic energy')
for a in axs:
    a.set_xlabel('time (a.u.)')
fig.suptitle('B2 classical projectile stopping (dashed = t*)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_stopping.png', dpi=150); plt.show()
# stopping over the clean transit (z: -13 -> +15, 28 Bohr)
i_star = int(np.argmin(np.abs(tk - t_star)))
dKE_tr = (KE[0] - KE[i_star]) * HA
dz_tr = zk[i_star] - zk[0]
print(f'transit: KE {KE[0]*HA:.2f} -> {KE[i_star]*HA:.2f} eV over dz={dz_tr:.1f} Bohr')
print(f'  => mean stopping power S = dKE/dz = {dKE_tr/dz_tr:.3f} eV/Bohr  '
      f'(lost {dKE_tr:.1f} eV in the first traversal)')""")

md(r"""## Result B2.2 — projectile trajectory $z(t)$ (absolute and folded)

Left: the **absolute** $z(t)$ climbing past the box to $z\approx346$ Bohr — the
periodic re-crossings that the CAP cannot stop. Right: $z(t)$ **folded** into the
box $[-25,25)$, marking the free-region edges ($\pm15$) and $t^\star$.""")
code(r"""def fold(z, L=50.0): return ((z + L/2) % L) - L/2
fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8))
axs[0].plot(tk, zk, lw=2, color='C0'); axs[0].axvline(t_star, ls='--', color='0.4')
axs[0].axhline(15, ls=':', color='0.6')
axs[0].set_ylabel('z absolute (Bohr)'); axs[0].set_title('absolute (periodic re-crossings)')
axs[1].plot(tk, fold(zk), '.', ms=2, color='C0'); axs[1].axvline(t_star, ls='--', color='0.4')
for ze in (-15, 15): axs[1].axhline(ze, ls='--', color='0.5')
axs[1].set_ylabel('z folded (Bohr)'); axs[1].set_title('folded into box [-25,25)')
for a in axs: a.set_xlabel('time (a.u.)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_trajectory.png', dpi=150); plt.show()""")

md(r"""## Result B2.3 — energy difference $\Delta E(t)=E_{B1}-E_{B2}$

$$\Delta E(t) = E_{\rm B1}(t) - E_{\rm B2}(t).$$

This is the raw electronic-energy difference the user asked for. **Read it with
care:** both runs are losing energy fast as the CAP drains tens of bath electrons,
and the projectile slightly *rearranges* the bath, changing how much density sits in
the absorbing slabs — so $\Delta E$ mixes the projectile's perturbation with a
**difference in drainage trajectory** between the two runs. It is therefore **not**
a clean "energy deposited in the bath" number (the drainage difference is tens of
Ha, dwarfing the projectile's eV-scale loss). The **robust** projectile energy loss
is its own KE drop in Result B2.1 (6.6 eV over the transit) — that number does not
depend on the bath drainage at all.""")
code(r"""dE2 = B1['E'] - B2['E']
fig, ax = plt.subplots(figsize=(7.0, 4.2))
ax.plot(B2['t'], dE2, lw=2, color='C2'); ax.axvline(t_star, ls='--', color='0.4')
ax.annotate(f'projectile reaches z=+15\n(t* = {t_star:.1f} a.u.)',
            xy=(t_star, dE2[np.argmin(np.abs(B2['t']-t_star))]),
            xytext=(t_star+8, dE2.min()+0.4*(dE2.max()-dE2.min())),
            fontsize=8, color='0.3', arrowprops=dict(arrowstyle='->', color='0.5'))
ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\Delta E = E_{B1}-E_{B2}$ (Ha)')
ax.set_title('Energy difference: baseline minus classical-projectile run')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_energy_diff.png', dpi=150); plt.show()
print(f'dE(0)={dE2[0]:.3f}  dE(t*)~{dE2[np.argmin(np.abs(B2["t"]-t_star))]:.3f}  '
      f'dE(end)={dE2[-1]:.3f} Ha')""")

md(r"""## Result B2.4 — B2 energetics (components vs time)
Total / kinetic / Hartree / xc electronic energies of B2. As in B1 the collapse
toward zero is the bath draining into the CAP; the projectile's perturbation rides
on top.""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.2))
for key, lab in [('E','total'), ('kin','kinetic'), ('har','Hartree'), ('xc','xc')]:
    ax.plot(B2['t'], B2[key], lw=1.6, label=lab)
ax.axvline(t_star, ls='--', color='0.4')
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('energy (Ha)')
ax.set_title('B2 electronic energy components'); ax.legend()
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_energetics.png', dpi=150); plt.show()""")

md(r"""## Result B2.5 — integrated bath current $\int J_z\,dV(t)$
The box-integrated $z$-current is now the **bath's induced** current (the projectile
itself carries no electronic current — it is a classical ion). It tracks the drag
the moving projectile exerts on the gas.

> The spatially-resolved current-density *field* remains a deferred observable
> (not written in this run).""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.2))
ax.plot(B2['t'], B2['jz'], lw=2, color='C0')
ax.axvline(t_star, ls='--', color='0.4'); ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\int J_z\,dV$ (a.u.)')
ax.set_title('B2 integrated bath z-current (drag of the moving projectile)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_current.png', dpi=150); plt.show()""")

md(r"""## Result B2.6 — total electronic (bath) density GIF (xz slice)
Mid-$y$ slice of the bath density; the **green tick** marks the *exact* classical
projectile $z$ (folded into the box), CAP edges dashed. Same fixed colour scale as
B1/B3. With no electronic projectile this is the bath responding to (and draining
around) the moving classical charge.""")
code(r"""from IPython.display import Image, display
g = HERE/'fig_b2_density_xz.gif'
display(Image(filename=str(g))) if g.exists() else print('run precompute_b2.py first')""")

md(r"""## Result B2.7 — induced bath density (wake) GIF, $\delta n=$ B2 $-$ B1
Subtracting the B1 bath reference isolates the **pure bath wake** the classical
projectile drives — cleaner than B3 (no wavepacket spike). Red = pile-up ahead /
around, blue = depletion behind. Green tick = projectile position.""")
code(r"""g = HERE/'fig_b2_wake_xz.gif'
display(Image(filename=str(g))) if g.exists() else print('run precompute_b2.py first')""")

md(r"""## Result B2.8 — electric field $E_z$ GIF (FFT-Poisson on the bath density)
$E_z$ of the **bath** density (same kernel/scale as B3). Note this is the field of
the *electrons only*; the projectile's own bare Coulomb field is the driver, not
shown here.""")
code(r"""g = HERE/'fig_b2_efield_xz.gif'
display(Image(filename=str(g))) if g.exists() else print('run precompute_b2.py first')""")

# ----------------- Baseline 2 — TRANSIT WINDOW (t ≤ t*) ---------------------
md(r"""# Baseline 2 — transit window only ($t\le t^\star\approx10.3$ a.u.)

The **clean first traversal**: launch to the far free-region edge, before the
projectile re-enters the periodic box and the bath collapses. This is the window
the stopping power is read from.""")
code(r"""mask2 = B2['t'] <= t_star + 0.30
maskk = tk <= t_star + 0.30
print(f'transit: {mask2.sum()} obs samples, {maskk.sum()} track samples, '
      f't in [0,{B2["t"][mask2][-1]:.2f}]')""")

md(r"""## TB2.1 — stopping $v_z(t)$, $KE(t)$ (transit) + stopping power""")
code(r"""fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8))
axs[0].plot(tk[maskk], vzk[maskk], 'o-', ms=3, color='C3')
axs[0].axvline(t_star, ls='--', color='0.4')
axs[0].set_ylabel(r'$v_z$ (a.u.)'); axs[0].set_title('projectile speed (transit)')
axs[1].plot(tk[maskk], KE[maskk]*HA, 'o-', ms=3, color='C1')
axs[1].axvline(t_star, ls='--', color='0.4')
axs[1].set_ylabel('KE (eV)'); axs[1].set_title('projectile KE (transit)')
for a in axs: a.set_xlabel('time (a.u.)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_stopping_transit.png', dpi=150); plt.show()
ii = int(np.argmin(np.abs(tk - t_star)))
dKE = (KE[0]-KE[ii])*HA; dz = zk[ii]-zk[0]
print(f'clean transit stopping: lost {dKE:.2f} eV over {dz:.1f} Bohr  '
      f'=> S = {dKE/dz:.3f} eV/Bohr')""")

md(r"""## TB2.2 — energy difference $\Delta E=E_{B1}-E_{B2}$ (transit)""")
code(r"""dE2t = (B1['E']-B2['E'])[mask2]; tt2 = B2['t'][mask2]
fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(tt2, dE2t, 'o-', ms=3, color='C2'); ax.axvline(t_star, ls='--', color='0.4')
ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\Delta E$ (Ha)')
ax.set_title('Energy difference B1-B2 (transit)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_energy_diff_transit.png', dpi=150); plt.show()
print(f'raw dE change over transit = {(dE2t[-1]-dE2t[0]):+.3f} Ha '
      f'(drainage-trajectory difference, NOT the projectile energy loss --')
print(f' the clean projectile loss is the KE drop in B2.1: ~6.6 eV)')""")

md(r"""## TB2.3 — energetics (transit)""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.0))
for key, lab in [('E','total'), ('kin','kinetic'), ('har','Hartree'), ('xc','xc')]:
    ax.plot(B2['t'][mask2], B2[key][mask2], 'o-', ms=2, lw=1.4, label=lab)
ax.axvline(t_star, ls='--', color='0.4')
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('energy (Ha)')
ax.set_title('B2 energy components (transit)'); ax.legend()
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_energetics_transit.png', dpi=150); plt.show()""")

md(r"""## TB2.4 — integrated bath current (transit)""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(B2['t'][mask2], B2['jz'][mask2], 'o-', ms=3, color='C0')
ax.axvline(t_star, ls='--', color='0.4'); ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\int J_z\,dV$ (a.u.)')
ax.set_title('B2 integrated bath z-current (transit)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2_current_transit.png', dpi=150); plt.show()""")

md(r"""## TB2.5 — density / wake / $E_z$ GIFs (transit)
The clean first traversal, same fixed colour scales as the full run and as B1/B3.""")
code(r"""for g in ('fig_b2_density_xz_transit.gif', 'fig_b2_wake_xz_transit.gif',
          'fig_b2_efield_xz_transit.gif'):
    p = HERE/g
    display(Image(filename=str(p))) if p.exists() else print(f'missing {g}')""")

# ========================== Baseline 3 — the σ=0.5 wavepacket ================
md(r"""# Baseline 3 — the σ=0.5 wavepacket projectile

Now a $\sigma=0.5$ Bohr Gaussian wavepacket (WP), $k_0=2.711$ Bohr$^{-1}$ (100 eV),
launched at $z_0=-13$ moving $+z$, **same CAP ($\eta=-0.5$) and window as B1** so
B1 is the exact subtraction reference. The WP is injected as one extra electron
($N_0=163$). It crosses the 30-Bohr free region and exits through the far CAP.

**Reference time** — at its initial speed the WP reaches the far edge of the free
zone ($z=+15$) at $t^\star=(15-(-13))/k_0$. This marks where the "clean" transit
ends and the WP enters the absorber; it is drawn as a dashed line throughout.""")
code(r"""import csv
def load_obs(sub):
    t, E, kin, har, xc, jz = [], [], [], [], [], []
    with open(RES/sub/'raw/observables/observables.csv') as fh:
        for r in csv.DictReader(fh):
            t.append(float(r['time_au'])); E.append(float(r['energy_total']))
            kin.append(float(r['energy_kinetic'])); har.append(float(r['energy_hartree']))
            xc.append(float(r['energy_xc'])); jz.append(float(r['current_z']))
    return dict(t=np.array(t), E=np.array(E), kin=np.array(kin),
               har=np.array(har), xc=np.array(xc), jz=np.array(jz))
B3 = load_obs('b3_wp_E100'); B1 = load_obs('b1_eta0p50')
k0 = 2.711063; z0, zedge = -13.0, 15.0
t_star = (zedge - z0)/k0
print(f't_star (WP reaches free-zone edge z=+15) = {t_star:.2f} a.u.')""")

md(r"""## Result 6 — energy difference $\Delta E(t)=E_\mathrm{baseline}-E_\mathrm{WP}$

$$\Delta E(t) = E_{\rm B1}(t) - E_{\rm B3}(t).$$

B1 and B3 share the same draining bath and CAP. At $t=0$, $\Delta E\approx-E_{\rm
WP}$ (the WP is one extra electron, so B3 is *higher* in energy). As time advances,
$\Delta E$ also picks up the **difference in drainage** between a 163- and a
162-electron bath (the CAP removes the WP itself once it reaches the far slab), so
the later evolution is not a pure "energy deposited" curve — interpret the early
transit, and cross-check against the WP's momentum/current decay.""")
code(r"""dE = B1['E'] - B3['E']           # aligned (same time grid)
fig, ax = plt.subplots(figsize=(7.0, 4.2))
ax.plot(B3['t'], dE, lw=2, color='C2')
ax.axvline(t_star, ls='--', color='0.4')
ax.annotate(f'WP reaches free-zone edge\n(t* = {t_star:.1f} a.u.)',
            xy=(t_star, ax.get_ylim()[0]), xytext=(t_star+8, dE.min()+0.4*(dE.max()-dE.min())),
            fontsize=8, color='0.3', arrowprops=dict(arrowstyle='->', color='0.5'))
ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\Delta E = E_{B1}-E_{B3}$ (Ha)')
ax.set_title('Energy difference: baseline minus wavepacket')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_energy_diff.png', dpi=150); plt.show()
print(f'dE(0)={dE[0]:.3f}  dE(t*)~{dE[np.argmin(np.abs(B3["t"]-t_star))]:.3f}  dE(end)={dE[-1]:.3f} Ha')""")

md(r"""## Result 7 — B3 energetics (components vs time)

Total / kinetic / Hartree / xc energies of the full B3 run. The collapse toward
zero is the bath draining into the CAP (as in B1); the WP rides on top early.""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.2))
for key, lab in [('E','total'), ('kin','kinetic'), ('har','Hartree'), ('xc','xc')]:
    ax.plot(B3['t'], B3[key], lw=1.6, label=lab)
ax.axvline(t_star, ls='--', color='0.4')
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('energy (Ha)')
ax.set_title('B3 energy components'); ax.legend()
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_energetics.png', dpi=150); plt.show()""")

md(r"""## Result 8 — integrated probability current $\int J_z\,dV(t)$

The box-integrated $z$-current — the wavepacket's net current. It starts at the WP
momentum ($|k_0|=2.711$) and decays to ~0 as the WP is absorbed in the far CAP.

> **Note.** This is the *box-integrated* current. The spatially-resolved
> **current-density field** $\mathbf J(\mathbf r,t)$ is a deferred observable (needs
> the `current_density` VTI writer — see handover); it is **not** in this run, so
> the field-view of the current is not shown here.""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.2))
ax.plot(B3['t'], B3['jz'], lw=2, color='C0')
ax.axvline(t_star, ls='--', color='0.4'); ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\int J_z\,dV$ (a.u.)')
ax.set_title('B3 integrated z-current (WP momentum, decaying under absorption)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_current.png', dpi=150); plt.show()""")

md(r"""## Result 9 — momentum distribution $n(k)$, before vs after

$n(k)=\sum_i f_i|\tilde\psi_i(k)|^2$ binned in $|k|$: **total** (bath shells + WP)
and **WP-only**. Comparing $t=0$ to $t=\tau$ shows the WP peak at $k_0=2.711$
broaden/vanish (absorption + scattering) and the bath shells redistribute — the
background e–e scattering fingerprint.""")
code(r"""import csv
def mom_at(sub, want_step):
    k, ntot, nwp = [], [], []
    with open(RES/sub/'raw/observables/momentum_distribution.csv') as fh:
        for ln in fh:
            if ln.startswith('#') or ln.startswith('step'): continue
            p = ln.split(',')
            if int(p[0]) == want_step:
                k.append(float(p[2])); ntot.append(float(p[3])); nwp.append(float(p[4]))
    return np.array(k), np.array(ntot), np.array(nwp)
steps_m = []
with open(RES/'b3_wp_E100/raw/observables/momentum_distribution.csv') as fh:
    for ln in fh:
        if ln[0].isdigit(): steps_m.append(int(ln.split(',')[0]))
s0, sN = min(steps_m), max(steps_m)
k0v, t0_tot, t0_wp = mom_at('b3_wp_E100', s0)
kNv, tN_tot, tN_wp = mom_at('b3_wp_E100', sN)
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.6))
axs[0].plot(k0v, t0_tot, label='t=0'); axs[0].plot(kNv, tN_tot, label='t=end')
axs[0].axvline(k0, ls=':', color='0.5'); axs[0].set_title('total n(k)')
axs[1].plot(k0v, t0_wp, label='t=0'); axs[1].plot(kNv, tN_wp, label='t=end')
axs[1].axvline(k0, ls=':', color='0.5'); axs[1].set_title('WP-only n(k)')
for a in axs:
    a.set_xlabel(r'$|k|$ (1/Bohr)'); a.set_ylabel('n(k)'); a.legend(fontsize=8); a.set_xlim(0, 5)
fig.suptitle('Momentum distribution before/after (k0=2.711 dotted)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_momentum.png', dpi=150); plt.show()""")

md(r"""## Result 10 — wavepacket centroid track $z(t)$

WP position via the centroid of the *positive* WP-induced density (B3 − B1). Valid
while the WP is in the box; it crosses $z_0=-13\to+15$ and is absorbed near $t^\star$.""")
code(r"""import csv
t_c, zc = [], []
with open(HERE/'cap_b3_wp_centroid.csv') as fh:
    for r in csv.DictReader(fh):
        t_c.append(float(r['t_au'])); zc.append(float(r['wp_centroid_z']))
fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(t_c, zc, 'o-', ms=3, color='C3')
ax.axhline(15, ls='--', color='0.5'); ax.axhline(-15, ls='--', color='0.5')
ax.axvline(t_star, ls='--', color='0.4')
ax.text(1, 15.6, 'free-zone edge z=+15', fontsize=8, color='0.4')
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('WP centroid z (Bohr)')
ax.set_title('B3 wavepacket centroid (from positive wake density)')
ax.set_xlim(0, 40)
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_centroid.png', dpi=150); plt.show()""")

md(r"""## Result 11 — total electronic density GIF (xz slice)
Mid-$y$ slice, total density (bath + WP), time-animated; CAP edges dashed.""")
code(r"""from IPython.display import Image, display
g = HERE/'fig_b3_density_xz.gif'
display(Image(filename=str(g))) if g.exists() else print('run precompute_b3.py first')""")

md(r"""## Result 12 — WP-induced density (wake) GIF, $\delta n = $ B3 $-$ B1

Subtracting the B1 bath reference isolates the wavepacket and the density it
displaces — the wake. Red = enhancement (the WP), blue = depletion.""")
code(r"""g = HERE/'fig_b3_wake_xz.gif'
display(Image(filename=str(g))) if g.exists() else print('run precompute_b3.py first')""")

md(r"""## Result 13 — electric field $E_z$ GIF (FFT-Poisson on the total density)

$\mathbf E=-\nabla\phi$ with $\nabla^2\phi=-4\pi\rho$, $\rho=-n$, $G{=}0$ removed
(neutralizing background) — `inqview.analysis.electric_field`. Mid-$y$ slice of
$E_z$: where the electrostatic field concentrates around the WP and the draining
edges.""")
code(r"""g = HERE/'fig_b3_efield_xz.gif'
display(Image(filename=str(g))) if g.exists() else print('run precompute_b3.py first')""")

# ================= Baseline 3 — TRANSIT WINDOW (t ≤ t*) ======================
md(r"""# Baseline 3 — transit window only ($t \le t^\star \approx 10.3$ a.u.)

Everything above recomputed on the **clean transit window**: from launch to when
the wavepacket (at its initial speed) reaches the far edge of the free region
($z=+15$). This is the physically meaningful window — afterwards the WP is in the
absorber and the bath is collapsing. The B1 free-region survival here is ≳ the
$t=10$ value (~89% at $\eta=-0.5$), so the subtraction reference is still clean.""")
code(r"""mask = B3['t'] <= t_star + 0.30     # include the frame at/just past t*
print(f'transit window: {mask.sum()} samples, t in [0, {B3["t"][mask][-1]:.2f}] a.u.')""")

md(r"""## T6 — energy difference $\Delta E=E_{B1}-E_{B3}$ (transit)""")
code(r"""dE = (B1['E'] - B3['E'])[mask]; tt = B3['t'][mask]
fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(tt, dE, 'o-', ms=3, color='C2')
ax.axvline(t_star, ls='--', color='0.4')
ax.annotate(f't* = {t_star:.1f} a.u.', xy=(t_star, dE.max()),
            xytext=(t_star-3.5, dE.max()), fontsize=8, color='0.3')
ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\Delta E$ (Ha)')
ax.set_title('Energy difference (transit window)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_energy_diff_transit.png', dpi=150); plt.show()
print(f'dE(0)={dE[0]:.3f}  dE(t*)={dE[-1]:.3f} Ha   -> deposited ~{dE[0]-dE[-1]:.3f} Ha')""")

md(r"""## T7 — energetics (transit)""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.0))
for key, lab in [('E','total'), ('kin','kinetic'), ('har','Hartree'), ('xc','xc')]:
    ax.plot(B3['t'][mask], B3[key][mask], 'o-', ms=2, lw=1.4, label=lab)
ax.axvline(t_star, ls='--', color='0.4')
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('energy (Ha)')
ax.set_title('B3 energy components (transit)'); ax.legend()
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_energetics_transit.png', dpi=150); plt.show()""")

md(r"""## T8 — integrated current $\int J_z\,dV$ (transit)""")
code(r"""fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(B3['t'][mask], B3['jz'][mask], 'o-', ms=3, color='C0')
ax.axvline(t_star, ls='--', color='0.4'); ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\int J_z\,dV$ (a.u.)')
ax.set_title('B3 integrated z-current (transit)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_current_transit.png', dpi=150); plt.show()""")

md(r"""## T9 — momentum distribution: $t=0$ vs $t\approx t^\star$

"After" is now taken at the transit edge (not $\tau$): the WP peak at $k_0$ before
it enters the absorber, so this is the in-flight scattering of the bath, with the
WP still present.""")
code(r"""st_star = min(steps_m, key=lambda s: abs(s*0.02 - t_star))
kSv, tS_tot, tS_wp = mom_at('b3_wp_E100', st_star)
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.6))
axs[0].plot(k0v, t0_tot, label='t=0'); axs[0].plot(kSv, tS_tot, label=f't={st_star*0.02:.1f}')
axs[0].axvline(k0, ls=':', color='0.5'); axs[0].set_title('total n(k)')
axs[1].plot(k0v, t0_wp, label='t=0'); axs[1].plot(kSv, tS_wp, label=f't={st_star*0.02:.1f}')
axs[1].axvline(k0, ls=':', color='0.5'); axs[1].set_title('WP-only n(k)')
for a in axs:
    a.set_xlabel(r'$|k|$ (1/Bohr)'); a.set_ylabel('n(k)'); a.legend(fontsize=8); a.set_xlim(0, 5)
fig.suptitle('Momentum distribution: t=0 vs transit edge')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_momentum_transit.png', dpi=150); plt.show()""")

md(r"""## T10 — wavepacket centroid $z(t)$ (transit)""")
code(r"""import csv
t_c, zc = [], []
with open(HERE/'cap_b3_wp_centroid_transit.csv') as fh:
    for r in csv.DictReader(fh):
        t_c.append(float(r['t_au'])); zc.append(float(r['wp_centroid_z']))
fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(t_c, zc, 'o-', ms=4, color='C3')
ax.axhline(15, ls='--', color='0.5'); ax.axhline(-13, ls=':', color='0.6')
ax.axvline(t_star, ls='--', color='0.4')
ax.text(0.5, 15.4, 'free-zone edge z=+15', fontsize=8, color='0.4')
ax.text(0.5, -12.4, 'launch z=-13', fontsize=8, color='0.5')
ax.set_xlabel('time (a.u.)'); ax.set_ylabel('WP centroid z (Bohr)')
ax.set_title('B3 wavepacket centroid (transit)')
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b3_centroid_transit.png', dpi=150); plt.show()
# crude transit speed from a linear fit
import numpy as np
zc = np.array(zc); t_c = np.array(t_c); ok = np.isfinite(zc)
if ok.sum() > 2:
    v = np.polyfit(t_c[ok], zc[ok], 1)[0]
    print(f'centroid speed ~ {v:.2f} a.u.  (launch k0 = {k0:.2f})')""")

md(r"""## T11 — total density GIF (transit)""")
code(r"""from IPython.display import Image, display
g = HERE/'fig_b3_density_xz_transit.gif'
display(Image(filename=str(g))) if g.exists() else print('run precompute_b3.py T_MAX=10.6 SUFFIX=_transit')""")

md(r"""## T12 — WP-induced density (wake) GIF (transit)""")
code(r"""g = HERE/'fig_b3_wake_xz_transit.gif'
display(Image(filename=str(g))) if g.exists() else print('missing transit wake gif')""")

md(r"""## T13 — electric field $E_z$ GIF (transit)""")
code(r"""g = HERE/'fig_b3_efield_xz_transit.gif'
display(Image(filename=str(g))) if g.exists() else print('missing transit efield gif')""")

# ===================== cross-baseline comparison (B2 vs B3) ==================
md(r"""# B2 vs B3 — classical vs quantum projectile

> **⚠ Width caveat (σ-convention, documented 2026-06-21).** B2 and B3 were both
> labelled "σ=0.5", but under the unified wavepacket convention they are **not the
> same physical width**: B2 (classical) used charge std 0.5 (= unified σ_wp **0.707**),
> while B3 (WP) used σ_wp 0.5 (density std 0.354). The two clouds differ by √2, so the
> comparison below conflates a *width* difference with the *classical-vs-quantum*
> difference — read it qualitatively only. The dedicated quantum-vs-classical study
> uses width-matched projectiles (unified σ; see CONTEXT.md "σ-convention unification").

Both projectiles (100 eV, $z_0=-13$, same CAP/window). The energy the bath absorbs,
$\Delta E=E_{B1}-E_{\rm run}$, is broadly comparable. The wake GIFs (B2$-$B1 vs
B3$-$B1) share one colour scale, so the bath response to a *classical* charge vs a
*quantum wavepacket* can be read side by side (modulo the √2 width caveat above).""")
code(r"""fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(B2['t'], B1['E']-B2['E'], lw=2, color='C0', label=r'B2 classical: $E_{B1}-E_{B2}$')
ax.plot(B3['t'], B1['E']-B3['E'], lw=2, color='C3', label=r'B3 wavepacket: $E_{B1}-E_{B3}$')
ax.axvline(t_star, ls='--', color='0.4'); ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlim(0, 40)
ax.set_xlabel('time (a.u.)'); ax.set_ylabel(r'$\Delta E$ (Ha)')
ax.set_title('Energy difference vs baseline: classical (B2) vs wavepacket (B3)')
ax.legend()
fig.tight_layout(); fig.savefig(HERE_FIG/'fig_b2b3_energy_diff.png', dpi=150); plt.show()
print('Note: B3 dE starts near -KE_WP (the WP energy is IN the electron total);')
print('      B2 dE is purely the bath response (classical KE is not in E_total).')""")

# --------------------------------------------------------------- 6. takeaway
md(r"""## Takeaway

- **The inq-study CAP works in an interacting jellium bath** — stable ETRS
  propagation, energy real and finite, no NaN, to $\tau=140$ a.u. (engine gate
  passed), for all four baselines. *PROVISIONAL until Task #7.*
- **A CAP on the whole bath drains it almost completely given time** (95–97% by
  140 a.u. at every depth/run) — the vacuum-optimal $\eta=-0.5$/20 Bohr is *not* a
  "leave the bath alone" absorber, because the slabs are permanently full of gas.
- **But the free region survives the measurement window**: ~89% intact at $t=10$
  ($\eta=-0.5$), ~94% ($\eta=-0.10$), while the projectile crosses in ~10 a.u. So
  **B2/B3 wake runs are viable with B1 as the exact subtraction reference**, read
  off the early transit — not the late collapse.
- **Baseline 2 (classical) — stopping is resolved.** The classical projectile
  **decelerates** $v_z: 2.711\!\to\!2.540$ (KE $3.675\!\to\!3.226$ Ha), losing
  **12.2 eV total** and **~6.6 eV over the clean first traversal** (28 Bohr)
  ⇒ mean stopping power **$S\approx0.24$ eV/Bohr** at 100 eV. **Caveat:** the CAP
  does *not* absorb a classical ion, so it flies through the periodic box ~7 times
  (final $z=346$ Bohr); only the **first traversal** (transit window) is a clean
  single-pass measurement.
- **Baseline 3 (WP) is done** ($N_0=163$, absorbed 97% by $\tau$, no NaN). The WP
  reaches the free-zone edge at $t^\star\approx10.3$ a.u. and *is* absorbed in the
  far CAP (it is electron density); its integrated current and momentum peak at
  $k_0=2.711$ decay to ~0, and the wake ($\delta n=$ B3$-$B1) + $E_z$ field show
  the packet crossing and draining.
- **$\eta=-0.10$ preserves the bath noticeably better** while still absorbing;
  a candidate production depth to compare against $-0.5$.
- **Open:** the spatially-resolved current-density *field* + flux screens (deferred
  observables) and total-system $n(k)$ field remain to add; a clean single-pass
  classical stopping number wants either a non-periodic launch or stopping the
  projectile at the first CAP crossing.""")

# ----------------------------------------------------------------- execute+write
# inject HERE_FIG into the namespace via a preamble cell
C.insert(0, new_code_cell(f"from pathlib import Path\nHERE = HERE_FIG = Path(r'{HERE}')"))

ep = ExecutePreprocessor(timeout=900, kernel_name='python3')
ep.preprocess(nb, {'metadata': {'path': str(HERE)}})
out = HERE / 'cap_baselines_study.ipynb'
with open(out, 'w') as fh:
    nbf.write(nb, fh)
print(f'wrote {out}  ({len(nb.cells)} cells)')
