#!/usr/bin/env python3
"""Build + execute the two-sided CAP-vs-mask study notebook (ADR 0007; built per
the notebook-making skill). Plan: docs/plans/twosided-cap-vs-mask.md.

Parses twosided_cap_vs_mask/run_{cap,mask}_*/results/{epsilon.txt,inner_norm_vs_time.csv},
writes the combined CSV and the 7 locked results, then the executed notebook.
Partial-tolerant: builds from whatever runs exist.

    PYTHONPATH=.../inq-stack/python /local/.../venv/bin/python3 build_twosided_report.py

ε PROVISIONAL until the inq-study engine regression (Task #7).
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path
import csv

HERE = Path(__file__).resolve().parent
VAC = HERE.parent.parent
SWEEP = VAC / "twosided_cap_vs_mask"
ANCHOR_E, ETA_STAR, ANCHOR_L = 10, -0.50, 20

def parse(p):
    out = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(' ')
        try: out[k] = float(v)
        except ValueError: out[k] = v
    return out

recs = []
for mode in ("cap", "mask"):
    for d in sorted(SWEEP.glob(f"run_{mode}_*")):
        f = d / "results/epsilon.txt"
        if f.exists():
            try:
                r = parse(f); r["name"] = d.name; recs.append(r)
            except Exception: pass
cols = ["name", "mode", "E_eV", "k0", "L_total", "Lhalf", "eta_Ha", "epsilon",
        "absorbed_fraction", "z_in", "Lcell_z", "tau", "N_STEPS", "dx", "sigma"]
with open(HERE / "twosided_combined.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in recs: w.writerow([r.get(c, "") for c in cols])
print(f"{len(recs)} runs -> twosided_combined.csv")

nb = new_notebook(); C = nb.cells

C.append(new_markdown_cell(r"""# Two-sided absorbers — sin² CAP vs mask, comfortable-region search

**Question.** For a free electron wavepacket leaving a finite box, which absorber —
the built-in **sin² complex absorbing potential (CAP)** or the **mask-function
absorber (MFA)** — and at what **width L** and (for the CAP) **depth η**, holds the
boundary **reflection error ε** low across the energy band? Here both absorbers span
**both** z-boundaries (width L split L/2 at each end). The output is the ε(E,L) maps
from which the **comfortable region** (L*, η*) for production is read off (by the
user — not auto-declared).

## Absorbers

**Two-sided sin² CAP** — two built-in `perturbations::absorbing` slabs, one per end,
composed with `perturbations::sum` (no new engine code; `inq/` and `inq-study`
untouched):
$$ V_\mathrm{CAP}(z)=i\,\eta\,\sin^2\!\Big(\tfrac{\pi(|z|-z_\mathrm{in})}{2\,L_\mathrm{half}}\Big),\quad z_\mathrm{in}<|z|<z_\mathrm{in}+L_\mathrm{half} $$

**Two-sided mask (MFA)** — `inqkit::absorbers::TwoSidedMaskAbsorber`, a real envelope
multiplied onto the orbital after each ETRS step (paper Eq. 12–13):
$$ M(z)=1-\sin^2\!\Big(\tfrac{\pi(|z|-z_\mathrm{in})}{2\,L_\mathrm{half}}\Big),\quad z_\mathrm{in}<|z|<z_\mathrm{in}+L_\mathrm{half};\ \ M=1\ \text{inner},\ M=0\ \text{wall} $$

| symbol | meaning |
|---|---|
| $\eta$ | CAP depth (Ha), $\eta<0$ absorbs. Mask has **no depth knob** — its only knob is L |
| $L$ / $L_\mathrm{half}=L/2$ | total / per-end absorber width (Bohr) |
| $z_\mathrm{in}=6\sigma$ | inner-region half-width; absorber occupies $z_\mathrm{in}<|z|<z_\mathrm{in}+L_\mathrm{half}$ |
| ETRS | required: absorption is non-Hermitian/non-unitary; Crank–Nicolson would cancel it |

## Reflection error ε (figure of merit) and absorption time
$$ \varepsilon=\frac{\int_{|z|<z_\mathrm{in}}|\psi_\mathrm{wp}(\tau)|^2\,dV}{N_0},\qquad
t_\mathrm{abs}=\min\{t:\ \varepsilon(t)<0.01\} $$
ε = surviving un-absorbed fraction (lower = better; symmetric inner region $|z|<z_\mathrm{in}$).
$t_\mathrm{abs}$ = time for the inner norm to fall below 1 % of $N_0$ (how *fast*).

> **Provisional.** All CAP ε are PROVISIONAL until the inq-study engine regression
> (Task #7). Mask ε do not depend on it. Source: De Giovannini, Larsen & Rubio,
> *Eur. Phys. J. B* **88**, 56 (2014), §IV & Eq. 12–13; Riss & Meyer 1993.

## Simulation setup (reconstructable)

Energy-scaled quasi-monochromatic packet `σ=4√2/k0` (≈12 % momentum spread), so the
box scales with energy and the packet stays a clean beam. Box `Lcell_z=12σ+L`,
thin transverse cell; WP launched at `z0=−z_in+4σ` (≥4σ from the near absorber ⇒
negligible density in its CAP region), moving +z toward the far absorber. ETRS,
`dt=0.01`, `dx=clamp(0.75/k0,0.18,0.30)`, `τ=2(z_in+L)/k0`. Grids: E ∈ {1…1000} eV;
L ∈ {10,16,20,26,30}; η ∈ {−0.3,−0.5,−0.7,−1.0} (CAP); anchor E=10 eV.

## Source files
| role | path |
|---|---|
| run binary (cap+mask modes) | `ResearchProject/systems/vacuum/scripts/twosided_cap_vs_mask/run.cpp` |
| dispatcher (3 phases, emailed) | `ResearchProject/systems/vacuum/scripts/twosided_cap_vs_mask/dispatch.py` |
| two-sided mask + ε | `inq-stack/include/inqkit/absorbers/mask_absorber.hpp` (`TwoSidedMaskAbsorber`, `inner_region_norm_twosided`) |
| mask shape + pure test | `inq-stack/include/inqkit/absorbers/mask_shape.hpp`, `inq-stack/tests/include/inqkit/absorbers/test_mask_shape.cpp` |
| this builder | `ResearchProject/systems/vacuum/hypotheses/twosided_cap_vs_mask/build_twosided_report.py` |
"""))

C.append(new_code_cell(r"""import csv, re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from inqview.visualisation import style
style.apply_theme()
HERE=Path(r"%s"); SWEEP=Path(r"%s"); ANCHOR_E=%d; ETA_STAR=%g
rows=list(csv.DictReader(open(HERE/'twosided_combined.csv')))
for r in rows:
    for k in ('E_eV','k0','L_total','Lhalf','eta_Ha','epsilon','absorbed_fraction','z_in','sigma'):
        try: r[k]=float(r[k])
        except: r[k]=np.nan
cap=[r for r in rows if r['mode']=='cap']; mask=[r for r in rows if r['mode']=='mask']
print(f"{len(cap)} cap + {len(mask)} mask runs")
def parsef(p):
    o={}
    for ln in Path(p).read_text().splitlines():
        k,_,v=ln.partition(' ')
        try:o[k]=float(v)
        except:o[k]=v
    return o
""" % (HERE.as_posix(), SWEEP.as_posix(), ANCHOR_E, ETA_STAR)))

# 1. eps(E) by eta (CAP, L=20) — log | linear(0-1) side by side, 1% guide line
C.append(new_markdown_cell(r"""## 1. ε(E) by CAP depth η (L=20 Bohr) — *how plotted:* each curve fixes η, plots ε vs energy at the reference width L=20 Bohr. **Left = log scale** (resolves the small-ε differences); **right = linear scale 0–1** (reads off "how close to a perfect absorber"). Dashed line = the **1 % reflection** level (ε=0.01)."""))
C.append(new_code_cell(r"""sub=[r for r in cap if abs(r['L_total']-20)<0.5]
fig,(axL,axR)=plt.subplots(1,2,figsize=(11,4.6))
for eta in sorted({round(r['eta_Ha'],3) for r in sub}):
    pts=sorted([r for r in sub if abs(r['eta_Ha']-eta)<1e-6],key=lambda r:r['E_eV'])
    if pts:
        E=[p['E_eV'] for p in pts];eps=[p['epsilon'] for p in pts]
        axL.plot(E,eps,'o-',label=f'η={eta}');axR.plot(E,eps,'o-',label=f'η={eta}')
for ax in (axL,axR):
    ax.set_xscale('log');ax.set_xlabel('E (eV)');ax.set_ylabel('ε (reflectivity)')
    ax.axhline(0.01,ls='--',lw=1,color='0.4')
    ax.text(0.01,0.012,'1%',transform=ax.get_yaxis_transform(),fontsize=7,color='0.4',va='bottom')
axL.set_yscale('log');axL.set_title('log scale')
axR.set_ylim(0,1);axR.set_title('linear scale (0–1)')
axL.legend(fontsize=8)
fig.suptitle('Two-sided CAP: ε(E) by depth η  (L=20 Bohr)')
fig.tight_layout();fig.savefig(HERE/'fig_eps_vs_E_by_eta.png',dpi=140)
"""))

# 1b. eps vs eta at L=20, one curve per energy (orthogonal cut of the same data)
C.append(new_markdown_cell(r"""## 1b. ε vs CAP depth η at L=20 Bohr — *how plotted:* the orthogonal cut of §1 — **η on the x-axis**, one curve per energy. Shows whether deepening the absorber keeps helping or reverses (over-absorption / reflection off too-steep a potential). The η sweep exists **only at L=20**; energies with all four η present are 2, 10, 32, 100, 300 eV. **Left = log scale**, **right = linear scale 0–1**; dashed line = **1 %** (ε=0.01)."""))
C.append(new_code_cell(r"""sub20=[r for r in cap if abs(r['L_total']-20)<0.5]
fig,(axL,axR)=plt.subplots(1,2,figsize=(11,4.6))
for E in sorted({round(r['E_eV'],3) for r in sub20}):
    pts=sorted([r for r in sub20 if abs(r['E_eV']-E)<0.6],key=lambda r:r['eta_Ha'])
    if len({round(p['eta_Ha'],3) for p in pts})>1:
        et=[p['eta_Ha'] for p in pts];eps=[p['epsilon'] for p in pts]
        axL.plot(et,eps,'o-',label=f'{E:g} eV');axR.plot(et,eps,'o-',label=f'{E:g} eV')
for ax in (axL,axR):
    ax.set_xlabel('CAP depth η (Ha)');ax.set_ylabel('ε (reflectivity)')
    ax.axhline(0.01,ls='--',lw=1,color='0.4')
    ax.text(0.01,0.012,'1%',transform=ax.get_yaxis_transform(),fontsize=7,color='0.4',va='bottom')
axL.set_yscale('log');axL.set_title('log scale')
axR.set_ylim(0,1);axR.set_title('linear scale (0–1)')
axL.legend(fontsize=8,title='energy')
fig.suptitle('Two-sided CAP: ε vs depth η at L=20 Bohr (one curve per energy)')
fig.tight_layout();fig.savefig(HERE/'fig_eps_vs_eta_L20.png',dpi=140)
"""))

# 2. eps(E) by L, CAP vs mask — 2x2: top row log, bottom row linear(0-1), 1% guide line
C.append(new_markdown_cell(r"""## 2 & 3. ε(E) by width L — CAP (η=−0.5) vs mask — *how plotted:* one curve per total width L; **columns** = CAP at η=−0.5 (left) vs mask (no depth knob, right); **top row = log scale**, **bottom row = linear scale 0–1**. Dashed line = the **1 % reflection** level (ε=0.01)."""))
C.append(new_code_cell(r"""capL=[r for r in cap if abs(r['eta_Ha']-ETA_STAR)<1e-6]
fig,axes=plt.subplots(2,2,figsize=(11,8.4),sharex=True)
panels=[(capL,'CAP (η=-0.5)'),(mask,'mask (no depth knob)')]
for col,(data,ttl) in enumerate(panels):
    for L in sorted({round(r['L_total'],1) for r in data}):
        pts=sorted([r for r in data if abs(r['L_total']-L)<0.5],key=lambda r:r['E_eV'])
        if pts:
            E=[p['E_eV'] for p in pts];eps=[p['epsilon'] for p in pts]
            axes[0,col].plot(E,eps,'o-',label=f'L={int(L)}');axes[1,col].plot(E,eps,'o-',label=f'L={int(L)}')
    axes[0,col].set_title(f'{ttl} — log');axes[1,col].set_title(f'{ttl} — linear (0–1)')
for ax in axes.flat:
    ax.set_xscale('log');ax.axhline(0.01,ls='--',lw=1,color='0.4')
    ax.text(0.01,0.012,'1%',transform=ax.get_yaxis_transform(),fontsize=7,color='0.4',va='bottom')
for ax in axes[0,:]: ax.set_yscale('log')
for ax in axes[1,:]: ax.set_ylim(0,1);ax.set_xlabel('E (eV)')
for ax in axes[:,0]: ax.set_ylabel('ε (reflectivity)')
axes[0,0].legend(fontsize=8)
fig.suptitle('Two-sided absorbers: ε(E) by width L  (top: log, bottom: linear 0–1)')
fig.tight_layout()
fig.savefig(HERE/'fig_eps_vs_E_by_L_sidebyside.png',dpi=140)
"""))

# 4. eps vs L at 10 eV — log | linear(0-1) side by side, 1% guide line
C.append(new_markdown_cell(r"""## 4. ε vs total width L at 10 eV — *how plotted:* the "how thin can we go" cut at 10 eV; CAP(η=−0.5) and mask. **Left = log scale**, **right = linear scale 0–1**. Dashed line = the **1 % reflection** level (ε=0.01)."""))
C.append(new_code_cell(r"""def at_anchor(data): return sorted([r for r in data if abs(r['E_eV']-ANCHOR_E)<0.6],key=lambda r:r['L_total'])
fig,(axL,axR)=plt.subplots(1,2,figsize=(11,4.6))
for data,lab in [(capL,'CAP η=-0.5'),(mask,'mask')]:
    pts=at_anchor(data)
    if pts:
        L=[p['L_total'] for p in pts];eps=[p['epsilon'] for p in pts]
        axL.plot(L,eps,'s-',label=lab);axR.plot(L,eps,'s-',label=lab)
for ax in (axL,axR):
    ax.set_xlabel('total width L (Bohr)');ax.set_ylabel('ε (reflectivity)')
    ax.axhline(0.01,ls='--',lw=1,color='0.4')
    ax.text(0.01,0.012,'1%',transform=ax.get_yaxis_transform(),fontsize=7,color='0.4',va='bottom')
axL.set_yscale('log');axL.set_title('log scale')
axR.set_ylim(0,1);axR.set_title('linear scale (0–1)')
axL.legend(fontsize=9)
fig.suptitle(f'ε vs total width L at {ANCHOR_E} eV')
fig.tight_layout();fig.savefig(HERE/'fig_eps_vs_L_anchor.png',dpi=140)
"""))

# 5. heatmaps — 2x2: top row log10(eps), bottom row linear eps, CAP | mask
C.append(new_markdown_cell(r"""## 5. ε(E,L) heatmaps — *how plotted:* ε over the E×L grid, CAP(η=−0.5) and mask. **Top row = log₁₀ ε** (resolves the small-ε comfortable corner), **bottom row = linear ε** (absolute reflectivity). Colormap **reversed** so **bright/yellow = low reflectivity (comfortable)**, dark = high ε."""))
C.append(new_code_cell(r"""def grid(data):
    Es=sorted({round(r['E_eV'],3) for r in data});Ls=sorted({round(r['L_total'],1) for r in data})
    M=np.full((len(Ls),len(Es)),np.nan)
    for r in data:
        try:
            i=Ls.index(round(r['L_total'],1));j=Es.index(round(r['E_eV'],3));M[i,j]=r['epsilon']
        except: pass
    return Es,Ls,M
allv=[r['epsilon'] for r in capL+mask if r['epsilon']>0]
lvmin,lvmax=(np.log10(min(allv)),np.log10(max(allv))) if allv else (-3,0)
linmax=max([r['epsilon'] for r in capL+mask] or [1.0])
def draw(ax,Es,Ls,A,vmin,vmax):
    im=ax.imshow(A,aspect='auto',origin='lower',vmin=vmin,vmax=vmax,cmap='viridis_r',
                 extent=[0,len(Es),0,len(Ls)])
    ax.set_xticks(np.arange(len(Es))+0.5);ax.set_xticklabels([f'{e:g}' for e in Es],rotation=45,fontsize=7)
    ax.set_yticks(np.arange(len(Ls))+0.5);ax.set_yticklabels([f'{int(l)}' for l in Ls],fontsize=8)
    ax.set_xlabel('E (eV)');ax.set_ylabel('L total (Bohr)');return im
fig,axes=plt.subplots(2,2,figsize=(11,8.6))
for col,(data,ttl) in enumerate([(capL,'CAP η=-0.5'),(mask,'mask')]):
    Es,Ls,M=grid(data)
    imL=draw(axes[0,col],Es,Ls,np.log10(M),lvmin,lvmax);axes[0,col].set_title(f'{ttl} — log₁₀ ε')
    imN=draw(axes[1,col],Es,Ls,M,0,linmax);axes[1,col].set_title(f'{ttl} — linear ε')
fig.colorbar(imL,ax=axes[0,:],label='log₁₀ ε',shrink=0.85)
fig.colorbar(imN,ax=axes[1,:],label='ε (reflectivity)',shrink=0.85)
fig.suptitle('ε(E,L) — bright/yellow = comfortable, low ε (top: log₁₀, bottom: linear)')
fig.savefig(HERE/'fig_eps_heatmap.png',dpi=140)
"""))

# 6. t_absorb
C.append(new_markdown_cell(r"""## 6. Time-to-absorb — *how plotted:* from each run's inner-norm time series, the first time ε(t)<1%·N₀; shown vs energy for the widest CAP and mask (fastest, cleanest)."""))
C.append(new_code_cell(r"""def t_abs(name):
    f=SWEEP/name/'results/raw/observables/inner_norm_vs_time.csv'
    if not f.exists(): return np.nan
    import csv as _c
    rr=list(_c.DictReader(open(f)))
    for row in rr:
        if float(row['inner_norm_over_N0'])<0.01: return float(row['time_au'])
    return np.nan
fig,ax=plt.subplots(figsize=(7,4.6))
for data,lab in [(capL,'CAP η=-0.5'),(mask,'mask')]:
    Lmax=max({r['L_total'] for r in data}) if data else None
    pts=sorted([r for r in data if abs(r['L_total']-Lmax)<0.5],key=lambda r:r['E_eV']) if Lmax else []
    xs=[p['E_eV'] for p in pts];ys=[t_abs(p['name']) for p in pts]
    ax.plot(xs,ys,'o-',label=f'{lab} (L={int(Lmax) if Lmax else 0})')
ax.set_xscale('log');ax.set_xlabel('E (eV)');ax.set_ylabel('t_absorb (a.u.)')
ax.set_title('Time for inner norm to fall below 1% of N₀');ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(HERE/'fig_t_absorb.png',dpi=140)
"""))

# 7. density GIFs of best CAP and best mask
C.append(new_markdown_cell(r"""## 7. Density GIFs — best CAP and best mask (lowest ε at 10 eV)."""))
C.append(new_code_cell(r"""import matplotlib.animation as animation
try:
    import vtk; from vtk.util.numpy_support import vtk_to_numpy
    def zprof(p):
        rd=vtk.vtkXMLImageDataReader();rd.SetFileName(str(p));rd.Update();img=rd.GetOutput()
        dm=img.GetDimensions();a=vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(dm[2],dm[1],dm[0])
        sp=img.GetSpacing();org=img.GetOrigin();z=org[2]+sp[2]*np.arange(dm[2]);return z,a.sum(axis=(1,2))
    def best(data):
        pool=[r for r in data if abs(r['E_eV']-ANCHOR_E)<0.6] or data
        return min(pool,key=lambda r:r['epsilon']) if pool else None
    for data,out,lab in [(capL,'fig_best_cap_density.gif','CAP'),(mask,'fig_best_mask_density.gif','mask')]:
        sc=best(data)
        if not sc: continue
        d=SWEEP/sc['name'];vti=d/'results/raw/vti/density_wp'
        frames=sorted(vti.glob('density_wp_t*.vti'),key=lambda p:int(re.search(r'_t(\d+)',p.name).group(1)))
        if not frames: print(f'{lab}: no frames');continue
        z0,_=zprof(frames[0]);profs=[zprof(p)[1] for p in frames];ymax=max(p.max() for p in profs)*1.05
        r=parsef(d/'results/epsilon.txt');zin=r['z_in'];Lh=r['Lhalf']
        fig,ax=plt.subplots(figsize=(6.4,3.2));(line,)=ax.plot(z0,profs[0])
        for sgn in (+1,-1):
            ax.axvspan(sgn*zin,sgn*(zin+Lh),color='C3',alpha=0.15)
        ax.set_xlabel('z (Bohr)');ax.set_ylabel(r'$n_{\rm WP}$');ax.set_ylim(0,ymax)
        ax.set_title(f'best {lab}: L={int(r["L_total"])} ε={r["epsilon"]:.2e}')
        ttl=ax.text(0.02,0.92,'',transform=ax.transAxes,fontsize=8)
        def upd(k,line=line,profs=profs,ttl=ttl): line.set_ydata(profs[k]);ttl.set_text(f'{k+1}/{len(profs)}');return line,ttl
        animation.FuncAnimation(fig,upd,frames=len(profs),interval=120,blit=False).save(HERE/out,writer='pillow',dpi=90);plt.close(fig)
        print('wrote',out)
except Exception as ex:
    print('GIF step skipped:',ex)
"""))

C.append(new_markdown_cell(r"""![best CAP](fig_best_cap_density.gif)
![best mask](fig_best_mask_density.gif)

## 8. Density carpets (z–t) — CAP η=−0.5, L=30 Bohr, at 2 / 10 / 100 eV — *how plotted:* a static 2D map of the ⟂-summed WP density: **x = z (Bohr), y = time (a.u.), colour = density**. One panel per energy (lowest=2 eV, 10 eV, closest-to-100=100 eV). Dashed cyan = inner-region edge $|z|=z_\mathrm{in}$ (absorber onset); dotted cyan = box wall $|z|=z_\mathrm{in}+L/2$. The beam streak travels +z, reaches the absorber, and fades — a clean diagonal that vanishes at the wall means low ε."""))
C.append(new_code_cell(r"""targets=[('run_cap_E2_L30_eta0.50','2 eV (lowest)'),
         ('run_cap_E10_L30_eta0.50','10 eV'),
         ('run_cap_E100_L30_eta0.50','100 eV (≈ band top)')]
try:
    import vtk; from vtk.util.numpy_support import vtk_to_numpy
    def zprof(p):
        rd=vtk.vtkXMLImageDataReader();rd.SetFileName(str(p));rd.Update();img=rd.GetOutput()
        dm=img.GetDimensions();a=vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(dm[2],dm[1],dm[0])
        sp=img.GetSpacing();org=img.GetOrigin();z=org[2]+sp[2]*np.arange(dm[2]);return z,a.sum(axis=(1,2))
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.6))
    for ax,(name,lab) in zip(axes,targets):
        d=SWEEP/name;vti=d/'results/raw/vti/density_wp'
        frames=sorted(vti.glob('density_wp_t*.vti'),key=lambda p:int(re.search(r'_t(\d+)',p.name).group(1)))
        if not frames: ax.set_title(f'{lab}: no frames');ax.axis('off');continue
        steps=np.array([int(re.search(r'_t(\d+)',p.name).group(1)) for p in frames])
        z,_=zprof(frames[0]);carpet=np.array([zprof(p)[1] for p in frames]);t=steps*0.01  # dt=0.01 locked
        r=parsef(d/'results/epsilon.txt');zin=r['z_in'];Lh=r['Lhalf']
        im=ax.imshow(carpet,aspect='auto',origin='lower',extent=[z[0],z[-1],t[0],t[-1]],cmap='inferno')
        for sgn in (+1,-1):
            ax.axvline(sgn*zin,color='c',lw=0.8,ls='--');ax.axvline(sgn*(zin+Lh),color='c',lw=0.8,ls=':')
        ax.set_xlabel('z (Bohr)');ax.set_title(f'{lab}  (ε={r["epsilon"]:.2e})')
        fig.colorbar(im,ax=ax,shrink=0.85,label=r'$n_{\rm WP}$')
    axes[0].set_ylabel('time (a.u.)')
    fig.suptitle('Density carpets (z–t): CAP η=−0.5, L=30 Bohr — dashed=inner edge, dotted=box wall')
    fig.tight_layout();fig.savefig(HERE/'fig_density_carpet_L30.png',dpi=140);print('wrote fig_density_carpet_L30.png')
except Exception as ex:
    print('carpet step skipped:',ex)
"""))
C.append(new_markdown_cell(r"""![density carpets](fig_density_carpet_L30.png)

## 9. xz-plane density GIFs — CAP η=−0.5, L=30 Bohr, at 2 / 10 / 100 eV — *how plotted:* the animated 2D **xz slice** through the cell mid-plane (mid-y): **x (Bohr) horizontal, z (Bohr) vertical, colour = density**, frame = time. The beam moves up (+z) into the far absorber. Dashed cyan = inner edge $|z|=z_\mathrm{in}$; dotted cyan = box wall. Complements the static carpet of §8 — here you watch the packet enter the absorber and disappear."""))
C.append(new_code_cell(r"""xz_targets=[('run_cap_E2_L30_eta0.50','2 eV','fig_xz_density_E2.gif'),
            ('run_cap_E10_L30_eta0.50','10 eV','fig_xz_density_E10.gif'),
            ('run_cap_E100_L30_eta0.50','100 eV','fig_xz_density_E100.gif')]
try:
    import vtk; from vtk.util.numpy_support import vtk_to_numpy
    import matplotlib.animation as animation
    def volume(p):
        rd=vtk.vtkXMLImageDataReader();rd.SetFileName(str(p));rd.Update();img=rd.GetOutput()
        dm=img.GetDimensions();A=vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(dm[2],dm[1],dm[0])
        sp=img.GetSpacing();org=img.GetOrigin()
        x=org[0]+sp[0]*np.arange(dm[0]);z=org[2]+sp[2]*np.arange(dm[2]);return x,z,A
    for name,lab,out in xz_targets:
        d=SWEEP/name;vtidir=d/'results/raw/vti/density_wp'
        frames=sorted(vtidir.glob('density_wp_t*.vti'),key=lambda p:int(re.search(r'_t(\d+)',p.name).group(1)))
        if not frames: print(f'{lab}: no frames');continue
        if len(frames)>120: frames=frames[::len(frames)//120+1]
        x,z,A0=volume(frames[0]);midy=A0.shape[1]//2
        slabs=[volume(p)[2][:,midy,:] for p in frames]   # each (nz,nx): z rows, x cols
        vmax=max(s.max() for s in slabs)*1.02 or 1.0
        r=parsef(d/'results/epsilon.txt');zin=r['z_in'];Lh=r['Lhalf']
        fig,ax=plt.subplots(figsize=(3.0,5.4))
        im=ax.imshow(slabs[0],origin='lower',aspect='auto',extent=[x[0],x[-1],z[0],z[-1]],
                     cmap='inferno',vmin=0,vmax=vmax)
        for sgn in (+1,-1):
            ax.axhline(sgn*zin,color='c',lw=0.7,ls='--');ax.axhline(sgn*(zin+Lh),color='c',lw=0.7,ls=':')
        ax.set_xlabel('x (Bohr)');ax.set_ylabel('z (Bohr)')
        ax.set_title(f'xz density — {lab}\nL=30, η=-0.5, ε={r["epsilon"]:.2e}',fontsize=9)
        tt=ax.text(0.04,0.97,'',transform=ax.transAxes,color='w',fontsize=8,va='top')
        fig.colorbar(im,ax=ax,shrink=0.7,label='density');fig.tight_layout()
        def upd(k,im=im,slabs=slabs,tt=tt): im.set_data(slabs[k]);tt.set_text(f'{k+1}/{len(slabs)}');return im,tt
        animation.FuncAnimation(fig,upd,frames=len(slabs),interval=110,blit=False).save(HERE/out,writer='pillow',dpi=90);plt.close(fig)
        print('wrote',out)
except Exception as ex:
    print('xz gif step skipped:',ex)
"""))
C.append(new_markdown_cell(r"""<table><tr>
<td>![xz 2 eV](fig_xz_density_E2.gif)</td>
<td>![xz 10 eV](fig_xz_density_E10.gif)</td>
<td>![xz 100 eV](fig_xz_density_E100.gif)</td>
</tr></table>

## Takeaway (read from the maps)
- The **comfortable region** (L*, η*) for the CAP and L* for the mask — the user reads
  it off §4–5: the smallest total width holding ε below the chosen threshold across the
  band, at the best depth η.
- CAP vs mask: which reaches lower ε at equal width (§2–3), and which absorbs faster (§6).
- Recall the two-sided split: a total L absorbs the forward beam with only L/2, so the
  comfortable L here is ~2× a single-sided width.

*All CAP ε PROVISIONAL until Task #7. The locked (L*, η*) feeds the σ=0.5 production
baseline task — `docs/campaigns/absorbing_boundary/sigma0p5_baselines_with_locked_params.md`.*"""))

ep = ExecutePreprocessor(timeout=3600, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
out = HERE / "twosided_cap_vs_mask_study.ipynb"
with open(out, "w") as fh: nbf.write(nb, fh)
print("wrote", out)
