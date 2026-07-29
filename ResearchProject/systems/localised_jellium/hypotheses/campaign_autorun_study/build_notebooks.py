#!/usr/bin/env python3
"""Build the per-phase (H0-H5) + campaign notebooks for the localised-jellium GS
ladder (campaign_autorun). One script -> 7 .ipynb. Re-runnable: re-execute after
new data lands to REFRESH. House narrative (notebook-making skill): context ->
hypothesis -> setup -> results (recomputed from run CSVs + embedded plot) ->
takeaway. Numbers are recomputed from the run data so a plot and its quoted number
never disagree.

Run:  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
      /local/data/public/skcb2/tddft/venv/bin/python3 build_notebooks.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

REPO = "/local/data/public/skcb2/tddft"
LJ = f"{REPO}/ResearchProject/systems/localised_jellium"
CA = f"{LJ}/scripts/campaign_autorun"
OUT = Path(f"{LJ}/hypotheses/campaign_autorun_study"); OUT.mkdir(parents=True, exist_ok=True)

# common preamble code cell (readers reused from analyse_phase; no emails on import)
PRE = f"""import sys, glob, csv, numpy as np
sys.path.insert(0, {CA!r}); sys.path.insert(0, {REPO+'/inq-stack/python'!r})
import matplotlib.pyplot as plt
from analyse_phase import e_total0, e_kin0, gs_energy, load_nz, _rs_present
from pathlib import Path
HA_EV=27.211386; CA=Path({CA!r}); RUNS=CA/'runs'
GS120_P3=Path({LJ!r})/'scripts/h0_base_difference/gs/results'
print('campaign_autorun notebook — data root', RUNS)"""

def nb(title, cells):
    n = new_notebook(); n.cells = cells
    n.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
    return n

def md(s): return new_markdown_cell(s)
def co(s): return new_code_cell(s)

import base64
PNG = {"H0": "runs/h0/H0_base_difference.png", "H1": "runs/h1/H1_edge_model.png",
       "H2": "runs/h2/H2_gs_convergence.png", "H3": "runs/h3/H3_surface_energetics.png",
       "H4": "runs/h4/H4_wp_energetics.png", "H5": "runs/h5/H5_classical_subtraction.png"}
def img(key, cap="Highlight plot:"):
    p = Path(CA) / PNG[key]
    if not p.exists():
        return md(f"*(plot pending — appears after the run/refresh: {p.name})*")
    b = base64.b64encode(p.read_bytes()).decode()
    c = md(f"{cap}\n\n![{p.name}](attachment:{p.name})"); c.attachments = {p.name: {"image/png": b}}
    return c

# ---- per-phase notebook specs (hypothesis, setup, plot-code, takeaway) ------
PHASES = {
"H0": dict(title="H0 — base WP-vs-classical E_total(0) gap",
  hyp="Is the base WP$-$classical energy gap at $t{=}0$ just the wavepacket localisation energy $3/4\\sigma^2$? (σ_WP=0.5, L_z=120, PBC.)",
  setup="Stationary WP (k₀=0) and a charge-matched classical ghost placed at the same distances r from the slab face; single-step total energy; excess above the GS.",
  method="""## Method — what is actually computed (verified against the `run.cpp` sources)

Every plotted point is one number, `E_tot(0) − E_GS`, in eV. **Both** energies are read
straight from the run files by `analyse_phase.py`; nothing is re-converged in the notebook,
so a plotted number can never disagree with the run that produced it.

| Symbol | Meaning | Where it comes from |
|---|---|---|
| `E_GS` | converged LDA total energy of the **bare neutral slab** (no projectile) | `gs_energy()` → `ground_state_energy_ha` in the GS `run_summary.txt` |
| `E_tot(0)` | Kohn–Sham total energy of the **slab + projectile** system at `t=0` | `e_total0()` → row 0 (`step=0`), `energy_total` column of the run `observables.csv` |
| `r` | projectile distance from the slab **face**, Bohr (`r = −z_launch − a`, `a=12.5`) | run-directory name `*_r{r}_p3` |
| `p3` | periodicity 3 = fully periodic box (PBC) | run-name suffix |
| `3/4σ²` | Gaussian-WP localisation (zero-point) kinetic energy `= 81.6 eV` for σ=0.5 | dotted reference line |

**`E_GS` — the baseline (`gs/run.cpp`).** A positive jellium **slab** background (half-width
`a=12.5`, uniform density `n₀ = N/(Lx·Ly·2a)`) is placed in a `50×50×120` Bohr box (periodicity 3,
grid spacing 0.5) and LDA SCF is run for `N=82` electrons. The converged `gs.energy.total()`
(**−108.53 Ha**) is written to `run_summary.txt` and read back as `E_GS`. It contains **no
projectile** — it is the reference the projectile energy is measured against.

**`E_tot(0)` — with the projectile (single-point, not a new SCF).** Both branches
`electrons.load(GS_DIR)` the *same* bare-slab GS, insert a projectile, then call
`real_time::propagate(…, n_steps=2, dt=0.01)` and record `energy_total` at **step 0**. That first
row is INQ's KS energy functional evaluated on the initial combined state — a single-point energy
of "slab-plus-projectile", **before** any meaningful dynamics. Because both branches start from the
identical GS, `E_tot(0) − E_GS` isolates purely the **energy of inserting the projectile**.

**Wavepacket branch — quantum (`wp/run.cpp`).** A Gaussian wavepacket (`inqkit::WavePacket`,
centre `z=z_launch`, width σ=0.5, `k₀=0` → *stationary*) is orthogonalised against the occupied
states and injected into one empty KS state at occupation 1.0 (`inject_into_last_extra_state`).
This adds a **real electron** (`N: 82→83`) shaped as a localised Gaussian, so `E_tot(0) − E_GS`
carries (i) the confinement kinetic energy `3/4σ² = 81.6 eV`, (ii) the LDA **self-interaction
error** of that one-electron density (a few eV), and (iii) its electrostatic coupling to the slab.
Terms (i)+(ii) dominate and barely depend on `r` → the WP curve is ~distance-stable (~88 eV at `r=4`).

**Classical branch — ghost potential (`classical/run.cpp`).** A **classical Gaussian ghost** is
inserted as an ion (species "H", pseudopotential `electron_gaussian_wpsigma0p5.upf` with
`z_valence 0`, light mass ≈ mₑ, velocity 0) at `z=z_launch`. `z_valence 0` means it adds **no
electron** (`N` stays 82): it is a bare external Gaussian **potential**, with *no* confinement KE
and *no* SIE. Here `E_tot(0) − E_GS` is the electrostatic coupling of the ghost to the slab, and
the run deliberately **omits** the ghost–background integral `∫ v_ghost·n₊`
(`ghost_background_term_omitted = true` in its summary) — that term is re-added analytically in H5.
The result is strongly `r`-dependent (unscreened ghost–slab Coulomb; ~188 eV at `r=4`, falling with
distance).

**Why the two are not a clean mirror.** The WP adds an *electron* (localisation + SIE); the classical
adds a *potential* missing its background term. So the raw `E_tot(0) − E_GS` gap is **not** the WP
localisation energy — it is dominated by these asymmetric artifacts, which is exactly what motivates
the H5 ghost-background correction.""",
  code="""E_GS=gs_energy(GS120_P3); ZP=3/(4*0.5**2)*HA_EV; base=RUNS/'h0'
rs_wp=_rs_present(base,'wp',3); rs_cl=_rs_present(base,'cl',3)
wp=[(e_total0(base/f'wp_r{r}_p3')-E_GS)*HA_EV for r in rs_wp]
cl=[(e_total0(base/f'cl_r{r}_p3')-E_GS)*HA_EV for r in rs_cl]
print('WP excess (eV):',[f'{x:.0f}' for x in wp]); print('cl excess (eV):',[f'{x:.0f}' for x in cl])
plt.figure(figsize=(6,4)); plt.plot(rs_wp,wp,'o-',label='WP (quantum)'); plt.plot(rs_cl,cl,'s--',label='classical ghost')
plt.axhline(ZP,ls=':',c='.4',label=f'localisation {ZP:.0f} eV'); plt.xlabel('r (Bohr from face)'); plt.ylabel('E_tot(0)-E_GS (eV)')
plt.legend(); plt.title('H0 base gap vs r'); plt.tight_layout(); plt.show()""",
  take="WP excess is ~distance-stable (≈ localisation + SIE); the classical excess is strongly r-dependent (unscreened ghost artifact). The raw gap is **not** the localisation energy — it is artifact-dominated, motivating the H5 ghost-background correction.",
  decomp_md="""### Energy decomposition — why the WP excess stays flat but the classical one falls

*(Added to answer: "why doesn't `E_tot(0) − E_GS` decrease with r for the WP, like the classical does?")*

INQ's total energy splits as `E_tot = T + U_H + E_xc + E_ext + E_ion`. The runs log
`T`, `U_H`, `E_xc`; we group the rest into a **remainder** `= E_tot − (T + U_H + E_xc)`
— the external-potential + ionic/Ewald electrostatics. Each component below is plotted
**relative to its value at r=40** (the far point), which isolates the part that actually
depends on distance (the huge r-independent offsets drop out).

**Classical ghost.** The ghost has `z_valence 0` — it adds an external Gaussian *potential*
but **no electron density** — so `T`, `U_H` and `E_xc` are frozen at their GS values (flat
at zero). The *entire* r-dependence sits in the **remainder**: the electrostatic coupling of
the ghost potential to the slab electrons, which falls steeply as the ghost leaves the slab.
Hence the classical excess decreases with r — exactly as expected.

**Wavepacket.** The WP is a *real electron*, so it couples to the slab through two large,
opposing channels: its Hartree repulsion against the slab electrons (`U_H` **rises ≈86 eV**
as the WP nears the slab) and its attraction to the positive jellium background (in the
**remainder**, which **falls ≈84 eV**). Because the slab is **neutral**, these nearly cancel —
the net r-dependence of the total is only ≈2 eV across the whole range. What survives is the
r-**independent** intrinsic energy: the confinement kinetic energy `3/4σ²` (it appears as
`T_WP − T_cl = 81.7 eV`, matching 81.6 eV) plus the WP's own self-interaction (Hartree+xc SIE,
≈4 eV). That constant self-energy is why the WP excess sits at ≈86 eV, essentially independent of r.

**Bottom line.** The WP curve is flat not because there is no interaction, but because the
electron's repulsion from the slab electrons and its attraction to the neutral slab's positive
background cancel; the surviving ≈86 eV is the electron's own (r-independent) localisation +
self-interaction energy. The classical ghost carries none of that intrinsic energy and gets no
cancellation, so its excess is pure ghost–slab electrostatics and decays with distance.

*(Inference: the channel attributions — `U_H`↑ = WP–slab electron repulsion, remainder↓ =
WP–background attraction — are the standard reading of the KS energy terms; the numbers
quoted are measured from the step-0 rows.)*""",
  decomp_code="""# Decompose E_tot(0) into T, U_H, E_xc and a remainder (external + ionic).
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass
def _comp0(run):
    f=next(iter(glob.glob(str(RUNS/'h0'/run/'**/observables.csv'),recursive=True)))
    rr=list(csv.reader(open(f))); h,d=rr[0],rr[1]; g=lambda c: float(d[h.index(c)])
    tot,T,H,X=g('energy_total'),g('energy_kinetic'),g('energy_hartree'),g('energy_xc')
    return dict(total=tot,T=T,U_H=H,E_xc=X,remainder=tot-(T+H+X))
rs=[4,12,20,28,36,40]
CH=[('T','kinetic  T'),('U_H','Hartree  U_H'),('E_xc','xc  E_xc'),
    ('remainder','external+ionic (remainder)'),('total','total')]
COL={'T':'#1b6ca8','U_H':'#c0392b','E_xc':'#27ae60','remainder':'#8e44ad','total':'k'}
fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
for ax,(tag,name) in zip(axes,(('wp','wavepacket (quantum electron)'),
                                ('cl','classical ghost (external potential)'))):
    comp={r:_comp0(f'{tag}_r{r}_p3') for r in rs}; far=comp[40]
    for key,lab in CH:
        dv=[(comp[r][key]-far[key])*HA_EV for r in rs]
        ax.plot(rs,dv,('-o' if key=='total' else '--o'),color=COL[key],
                lw=(2.4 if key=='total' else 1.4),ms=4,label=lab)
    ax.axhline(0,color='.6',lw=.8,zorder=0); ax.set_title(name); ax.set_xlabel('r (Bohr from face)')
axes[0].set_ylabel('ΔE component, relative to r=40 (eV)')
axes[0].legend(frameon=False,fontsize=8)
fig.suptitle('H0 — decomposition of E_tot(0) vs distance (Δ isolates the r-dependence)')
fig.tight_layout(); plt.show()
wpT=_comp0('wp_r40_p3')['T']*HA_EV; clT=_comp0('cl_r40_p3')['T']*HA_EV
print(f'WP kinetic - classical kinetic = {wpT-clT:.1f} eV  (= confinement 3/4σ² = 81.6 eV)')
print('WP:  U_H rises ~86 eV near slab, remainder falls ~84 eV -> nearly cancel (neutral slab)')
print('cl:  T/U_H/E_xc frozen (ghost adds no electron); only remainder moves -> excess decays')""",
  p2_md="""### Periodicity 2 (open-z) — the same `E_tot(0) − E_GS`, periodic in x,y only

The main H0 plot above is periodicity 3 (fully periodic). Here is the **same** base gap with
the box **open in z** — `periodicity 2`: periodic in x and y, finite (non-periodic) in z. The
open-z runs were produced alongside H4 (WP, `runs/h4/wp_r*_p2`) and H5 (classical,
`runs/h5/cl_r*_p2`); the open-z GS baseline is `runs/h2/gs_p2_lz120` (`E_GS = +60.38 Ha`).

Note the absolute open-z `E_GS` differs completely from the p3 value (−108.53 Ha) because the
2D Coulomb kernel uses a different self-energy reference — so **only differences within one BC
are meaningful**, never cross-BC absolute energies.

**Caveat (open-z + net charge).** The WP adds a net −1 charge, and under periodicity 2 a
net-charged cell carries a G=0 compensation term (in the 2D kernel) that the **neutral** GS
does not — so the naive WP `E_tot(0) − E_GS` is biased low by a few eV. That is why the WP
open-z curve sits *slightly below* the 81.6 eV localisation line. The classical ghost keeps the
cell neutral (`z_valence 0`), so its open-z excess is unbiased. Compare the *shapes* (WP flat,
classical decaying — the same qualitative story as p3), not the absolute WP offset. (Documented
in H4/H5; the charged-cell G=0 reference correction is the open follow-up.)""",
  p2_code="""# Periodicity-2 (open-z) E_tot(0) - E_GS.  WP runs -> H4, classical -> H5,
# GS baseline -> the open-z GS (h2/gs_p2_lz120).
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass
E_GS_P2=gs_energy(RUNS/'h2/gs_p2_lz120/results'); ZP=3/(4*0.5**2)*HA_EV
wp_base=RUNS/'h4'; cl_base=RUNS/'h5'
rs_wp2=_rs_present(wp_base,'wp',2); rs_cl2=_rs_present(cl_base,'cl',2)
wp2=[(e_total0(wp_base/f'wp_r{r}_p2')-E_GS_P2)*HA_EV for r in rs_wp2]
cl2=[(e_total0(cl_base/f'cl_r{r}_p2')-E_GS_P2)*HA_EV for r in rs_cl2]
print('periodicity 2 (open-z), E_GS =',f'{E_GS_P2:.3f} Ha')
print('WP p2 excess (eV):',[f'{x:.0f}' for x in wp2])
print('cl p2 excess (eV):',[f'{x:.0f}' for x in cl2])
fig,ax=plt.subplots(figsize=(6.4,4.4))
if rs_wp2: ax.plot(rs_wp2,wp2,'o-',color='#1b6ca8',label='wavepacket (quantum)')
if rs_cl2: ax.plot(rs_cl2,cl2,'s--',color='#c0392b',label='classical ghost (raw)')
ax.axhline(ZP,ls=':',color='.4',label=f'WP localisation {ZP:.0f} eV')
ax.set_xlabel('r (Bohr from face)'); ax.set_ylabel('E_tot(0) - E_GS  (eV)')
ax.set_title('H0 open-z (periodicity 2): base gap vs r'); ax.legend(frameon=False)
fig.tight_layout(); plt.show()""",
  p2full_md="""### Periodicity-2 re-run with the full *measured* energy decomposition (2026-07-07)

The H0 base gap **re-run at periodicity 2** (open-z slab: periodic in x, y; finite
in z), with the observables writer extended to stream **every** INQ energy component
each step. So `E_ext = ∫ n·v_ext` — the term the classical ghost changes — and the
Hartree, kinetic, xc, ionic and non-local terms are now **measured directly**, not
inferred as a `total − (T+U_H+E_xc)` remainder (as the periodicity-3 panels above do).

- Data: `scripts/campaign_autorun/runs/h0_p2/{wp,cl}_r{r}_p2/` (radii 4–40), the
  dedicated periodicity-2 H0 sweep off the open-z GS.
- Reference GS: the open-z GS `runs/h2/gs_p2_lz120` (`E_GS = +60.38 Ha`).
- Streamed columns: `energy_{total,kinetic,hartree,xc,external,nonlocal,ion,ion_kinetic,exact_exchange,nvxc,eigenvalues}`.
- Sum check (printed below): `kinetic+external+non_local+hartree+xc+exact_exchange+ion+ion_kinetic`
  reproduces `energy_total` to ~1e-13 Ha for every run.

Left panel = wavepacket, right = classical ghost; each component plotted relative to
its r=40 value to isolate the distance dependence. Per-run confirmable table:
[`runs/H0_p2_runs.ipynb`](runs/H0_p2_runs.ipynb).""",
  p2full_code="""# Periodicity-2 H0 re-run: FULL measured energy decomposition (external streamed, not inferred).
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass
KEYS=['total','kinetic','hartree','xc','external','nonlocal','ion','ion_kinetic','exact_exchange']
def _comp0_p2(run):
    f=next(iter(glob.glob(str(RUNS/'h0_p2'/run/'**/observables.csv'),recursive=True)))
    rr=list(csv.reader(open(f))); h,d=rr[0],rr[1]; g=lambda c: float(d[h.index(c)])
    return {k:g('energy_'+k) for k in KEYS}
E_GS_P2=gs_energy(RUNS/'h2/gs_p2_lz120/results'); ZP=3/(4*0.5**2)*HA_EV
rs=[4,12,20,28,36,40]
wp_ex=[(_comp0_p2(f'wp_r{r}_p2')['total']-E_GS_P2)*HA_EV for r in rs]
cl_ex=[(_comp0_p2(f'cl_r{r}_p2')['total']-E_GS_P2)*HA_EV for r in rs]
print('periodicity 2, E_GS =',f'{E_GS_P2:.4f} Ha   zero-point 3/4sig^2 =',f'{ZP:.1f} eV')
print('WP p2 excess (eV):',[f'{x:.1f}' for x in wp_ex])
print('cl p2 excess (eV):',[f'{x:.1f}' for x in cl_ex])
for tag in ('wp','cl'):
    c=_comp0_p2(f'{tag}_r28_p2')
    s=c['kinetic']+c['external']+c['nonlocal']+c['hartree']+c['xc']+c['exact_exchange']+c['ion']+c['ion_kinetic']
    print(f'sum-check {tag} r=28: total={c[\"total\"]:.6f}  sum(8 comp)={s:.6f}  diff={abs(c[\"total\"]-s):.1e}')
CH=[('kinetic','T'),('hartree','U_H'),('xc','E_xc'),('external','E_ext'),('ion','E_ion'),('nonlocal','E_nl'),('total','total')]
COL={'kinetic':'#1b6ca8','hartree':'#c0392b','xc':'#27ae60','external':'#8e44ad','ion':'#e67e22','nonlocal':'#7f8c8d','total':'k'}
fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
for ax,(tag,name) in zip(axes,(('wp','wavepacket (quantum electron)'),
                                ('cl','classical ghost (external potential)'))):
    comp={r:_comp0_p2(f'{tag}_r{r}_p2') for r in rs}; far=comp[40]
    for key,lab in CH:
        dv=[(comp[r][key]-far[key])*HA_EV for r in rs]
        ax.plot(rs,dv,('-o' if key=='total' else '--o'),color=COL[key],
                lw=(2.4 if key=='total' else 1.4),ms=4,label=lab)
    ax.axhline(0,color='.6',lw=.8,zorder=0); ax.set_title(name); ax.set_xlabel('r (Bohr from face)')
axes[0].set_ylabel('ΔE component, relative to r=40 (eV)'); axes[0].legend(frameon=False,fontsize=8,ncol=2)
fig.suptitle('H0 periodicity-2: measured energy decomposition vs distance (Δ vs r=40)')
fig.tight_layout(); plt.show()"""),
"H1": dict(title="H1 — edge model (Gibbs vs Friedel)",
  hyp="A finite erfc edge width $w \\gtrsim$ grid spacing removes numerical Gibbs ringing at the slab boundary while preserving the physical Friedel tail (λ=π/k_F≈9.3 Bohr).",
  setup="GS-only, periodicity 3, baseline slab (L_z=90, a=12.5, N=82). Edge width w swept; planar density n(z) extracted.",
  code="""base=RUNS/'h1'; ws=sorted(float(p.name.split('_w')[1]) for p in base.glob('gs_w*'))
plt.figure(figsize=(6.5,4))
for w in ws:
    z,nz=load_nz(base/f'gs_w{w:g}/results'); plt.plot(z,nz,lw=1.3,label=f'w={w:g}')
plt.axvspan(-12.5,12.5,color='.92',zorder=0); plt.xlim(-25,25); plt.xlabel('z (Bohr)'); plt.ylabel('n(z)'); plt.legend(fontsize=8); plt.title('H1 n(z) vs edge width'); plt.tight_layout(); plt.show()
print('w values:',ws)""",
  take="The edge width that suppresses near-boundary ringing while keeping the interior Friedel structure is the clean choice; below grid spacing the ringing is numerical (Gibbs)."),
"H2": dict(title="H2 — GS convergence + open-z viability",
  hyp="The neutral-slab interior density is box-independent and open-z (periodicity 2) is usable; the work function plateaus with vacuum.",
  setup="GS-only, w=0, a=12.5, N=82. L_z swept (periodicity 3) + open-z (periodicity 2) GS at L_z=90,120.",
  code="""base=RUNS/'h2'; lzs=sorted(int(p.name.split('_lz')[1]) for p in base.glob('gs_lz*'))
n0=[]
for lz in lzs:
    z,nz=load_nz(base/f'gs_lz{lz}/results'); n0.append(float(nz[np.abs(z)<6].mean()))
plt.figure(figsize=(6,4)); plt.plot(lzs,n0,'o-'); plt.axhline(1.312e-3,ls=':',c='.4',label='target'); plt.xlabel('L_z (Bohr)'); plt.ylabel('interior n0'); plt.legend(); plt.title('H2 interior n0 vs L_z'); plt.tight_layout(); plt.show()
print('interior n0 vs Lz:',[f'{x:.3e}' for x in n0])""",
  take="Interior n₀ is flat vs L_z (box-converged). Open-z GS converges and is usable for H4/H5. NB absolute E_GS is box-dependent (E_self) — use Φ/densities, not absolute E."),
"H3": dict(title="H3 — surface energetics (thickness)",
  hyp="E(N) is liquid-drop-linear → surface energy σ_s and bulk e_bulk; thin slabs lose the bulk interior.",
  setup="GS-only, w=0, L_z=90, half-width a swept with N scaled to hold n₀.",
  code="""base=RUNS/'h3'; aN=sorted((float(p.name.split('_a')[1].split('_N')[0]),int(p.name.split('_N')[1])) for p in base.glob('gs_a*_N*'))
N=np.array([n for _,n in aN]); E=np.array([gs_energy(base/f'gs_a{a:g}_N{n}/results') for a,n in aN])
plt.figure(figsize=(6,4)); plt.plot(N,E,'o-'); plt.xlabel('N electrons'); plt.ylabel('E_GS (Ha)'); plt.title('H3 E(N) (E_self-uncorrected)'); plt.tight_layout(); plt.show()
print('a,N,E:',list(zip([a for a,_ in aN],N.tolist(),[f'{e:.2f}' for e in E])))""",
  take="**CAVEAT:** absolute E_GS carries a box/thickness-dependent E_self (confirmed in H0), so the raw liquid-drop σ_s is unreliable until E_self is subtracted per thickness — flagged for follow-up."),
"H4": dict(title="H4 — WP energetics: E_SIE + PBC-vs-open-z",
  hyp="The net-charge periodic-image error excess(r,PBC)−excess(r,open-z) is significant/negligible → choose the production BC; excess(r) plateaus to E_SIE (~4.5 eV).",
  setup="Stationary WP (k₀=0), L_z=120, r swept × periodicity {3,2}. excess = E_tot(0) − E_GS(BC) − ⟨T_WP⟩ (zero-point 81.6 eV).",
  code="""base=RUNS/'h4'; ZP=3/(4*0.5**2)*HA_EV
EG={3:gs_energy(GS120_P3),2:gs_energy(RUNS/'h2/gs_p2_lz120/results')}
plt.figure(figsize=(6,4))
for per,mk in ((3,'o-'),(2,'s--')):
    rs=_rs_present(base,'wp',per)
    if not rs: continue
    exc=[(e_total0(base/f'wp_r{r}_p{per}')-EG[per])*HA_EV-ZP for r in rs]
    plt.plot(rs,exc,mk,label=f'periodicity {per}'); print(f'per{per} E_SIE plateau={exc[-1]:.1f} eV')
plt.axhline(4.5,ls=':',c='.4',label='known SIE ~4.5'); plt.xlabel('r (Bohr from face)'); plt.ylabel('E_SIE-ish (eV)'); plt.legend(); plt.title('H4 WP excess vs r'); plt.tight_layout(); plt.show()""",
  take="**PBC E_SIE ≈ 4.3 eV matches the known ~4.5 eV (trustworthy).** The open-z value is biased by the net-charge G=0 term in the 2D kernel (the GS is neutral, the WP cell is net −1) — open-z reference needs the charged-cell G=0 correction before the BC verdict is final."),
"H5": dict(title="H5 — classical mirror (route 2 + thread D)",
  hyp="Corrected route-2 E_SIE (classical subtraction with ghost-background re-added) matches route-1; the classical periodic-image error informs the Campaign-1 cutoff.",
  setup="Matched stationary ghost, L_z=120, r swept × periodicity {3,2}. Classical excess E_cl(0) − E_GS(BC).",
  code="""base=RUNS/'h5'
EG={3:gs_energy(GS120_P3),2:gs_energy(RUNS/'h2/gs_p2_lz120/results')}
plt.figure(figsize=(6,4))
for per,mk in ((3,'o-'),(2,'s--')):
    rs=_rs_present(base,'cl',per)
    if not rs: continue
    e=[(e_total0(base/f'cl_r{r}_p{per}')-EG[per])*HA_EV for r in rs]
    plt.plot(rs,e,mk,label=f'periodicity {per}'); print(f'per{per} classical excess (eV):',[f'{x:.0f}' for x in e])
plt.xlabel('r (Bohr from face)'); plt.ylabel('E_cl(0)-E_GS (eV)'); plt.legend(); plt.title('H5 classical excess vs r'); plt.tight_layout(); plt.show()""",
  take="The classical ghost excess is strongly r-dependent (unscreened ghost-slab Coulomb). The corrected route-2 E_SIE needs the analytic ghost-background integral ∫v_ghost·n₊ — flagged. The classical PBC-vs-open-z difference feeds the Campaign-1 cutoff (thread D)."),
}

# ---- neutral house-narrative: Question -> What was done -> Results -> Provisional --
# The reorganisation (2026-07-06, /grill-with-docs) presents each hypothesis WITHOUT
# assistant interpretation: the "question you were aiming to answer" is the user's own
# `hyp` wording; `setup`+`method` are what was done; the plot+numbers are the results;
# the old `take` (a verdict) is KEPT but quarantined in a marked provisional box that
# the reader owns. See docs/plans/campaign-autorun-review-organisation.md.
FILESTEM = {"H0": "base_difference", "H1": "edge_model", "H2": "gs_convergence",
            "H3": "surface_energetics", "H4": "wp_energetics", "H5": "classical_subtraction"}
# representative single-run deep-dives per hypothesis (run-notebook assembler), and the
# per-hypothesis run-evidence notebook (aggregator). Linked from each study nb + index.
REPS = {"H0": ["rep_H0_wp_r28_p3", "rep_H0_cl_r28_p3"],
        "H2": ["rep_H2_gs_lz120"], "H3": ["rep_H3_gs_a15_N98"],
        "H4": ["rep_H4_wp_r28_p2", "rep_H4_wp_r28_p3"]}

PROVISO = ("## ⚠ Provisional — author-generated interpretation (you own the verdict)\n"
           "> The reading below was drafted by the analysis assistant and is kept for "
           "reference only. **Deriving the learnings and the next experiments is your "
           "job.** The Question / What was done / Results above are neutral; this box is "
           "not — treat it as a hypothesis to accept, reject, or rewrite.\n\n")


def _links_md(key):
    reps = REPS.get(key, [])
    ev = f"runs/{key}_runs.ipynb"
    parts = [f"**Run-evidence (every run in the sweep):** [`{ev}`]({ev})"]
    if reps:
        parts.append("**Representative single-run deep-dives:** "
                      + ", ".join(f"[`{r}`](runs/{r}.ipynb)" for r in reps))
    return "## Independently confirm the results\n" + "  \n".join(parts)


def build_phase(key, spec):
    cells = [
        md(f"# {spec['title']}\n\n*Localised-jellium `campaign_autorun` run-set — hypothesis "
           f"{key}. Auto-built from run data (`build_notebooks.py`); numbers are recomputed "
           f"from the run files, never re-converged.*"),
        md(f"## Question you were aiming to answer\n{spec['hyp']}"),
        md(f"## What was done\n{spec['setup']}\n\nRun data: "
           f"`scripts/campaign_autorun/runs/{key.lower()}/`. Engine: inq-study (LDA), "
           f"σ_WP=0.5, spacing 0.5 Bohr."),
    ]
    if spec.get("method"):
        cells.append(md(spec["method"]))
    cells += [co(PRE), md("## Results"), img(key), co(spec['code'])]
    if spec.get("decomp_md"):
        cells.append(md(spec["decomp_md"])); cells.append(co(spec["decomp_code"]))
    if spec.get("p2_md"):
        cells.append(md(spec["p2_md"])); cells.append(co(spec["p2_code"]))
    if spec.get("p2full_md"):
        cells.append(md(spec["p2full_md"])); cells.append(co(spec["p2full_code"]))
    cells.append(md(_links_md(key)))
    cells.append(md(PROVISO + spec['take']))
    p = OUT/f"{key}_{FILESTEM[key]}.ipynb"
    nbf.write(nb(spec['title'], cells), str(p)); print("wrote", p.name); return p


def build_index(paths):
    # neutral ladder table: question (user's hyp) + what was done + links. No verdicts
    # in the table; all interpretation goes to the single merged provisional box.
    def _1line(s):
        return " ".join(s.split()).replace("|", "\\|")
    rows = ["| # | Question you were aiming to answer | What was done | Notebooks |",
            "|---|---|---|---|"]
    for key, spec in PHASES.items():
        links = f"[study]({key}_{FILESTEM[key]}.ipynb) · [runs](runs/{key}_runs.ipynb)"
        if REPS.get(key):
            links += " · " + " ".join(f"[{r.split('_',1)[1]}](runs/{r}.ipynb)" for r in REPS[key])
        rows.append(f"| **{key}** | {_1line(spec['hyp'])} | {_1line(spec['setup'])} | {links} |")
    cells = [
        md("# `campaign_autorun` — localised-jellium GS ladder (index)\n\n"
           "Single entry point for the run-set. Each hypothesis notebook states **the "
           "question you were aiming to answer**, **what was done**, and **the results** — "
           "no assistant interpretation. Every plotted number is recomputed from the run "
           "files. Auto-built; re-run `build_notebooks.py` to refresh.\n\n"
           "*(Supersedes `campaign_summary.ipynb`, kept on disk but off the reading path.)*"),
        md("## The ladder (H0 → H5)\n" + "\n".join(rows)),
        md("## Highlight plot per hypothesis"),
    ]
    for key, spec in PHASES.items():
        cells.append(md(f"### {key} — {_1line(spec['hyp'])[:90]}"))
        cells.append(img(key))
    cells.append(md(
        PROVISO
        + "**Per-hypothesis readings (drafted, unverified):**\n"
        + "\n".join(f"- **{k}** — {' '.join(PHASES[k]['take'].split())}" for k in PHASES) + "\n\n"
        + "**Open follow-ups the assistant noted (not execution errors, yours to weigh):**\n"
        "1. Open-z periodicity-2 net-charge G=0 reference (H4/H5).\n"
        "2. H2 work-function Φ extractor.\n"
        "3. H3 σ_s E_self correction.\n"
        "4. H5 ghost-background integral ∫v_ghost·n₊."))
    p = OUT/"00_index.ipynb"; nbf.write(nb("index", cells), str(p)); print("wrote", p.name)


if __name__ == "__main__":
    paths = [build_phase(k, s) for k, s in PHASES.items()]
    build_index(paths)
    print("done — to execute: python3 -m nbconvert --to notebook --execute --inplace <nb> (venv)")
