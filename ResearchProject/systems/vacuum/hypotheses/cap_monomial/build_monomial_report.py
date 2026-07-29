#!/usr/bin/env python3
"""Build + execute the monomial-CAP benchmark notebook (ADR 0007 location).

Compares the inq-study MONOMIAL CAP (V=iη·s^n) against the built-in sin² CAP at
fixed thin L=5 Bohr, E=10 eV: does the wall-peaked ramp beat sin²'s ~0.20 floor,
and how does it compare to the width sin² needs (cap_Lopt_E10)? Writes the combined
CSV, ε(η) curves per order, a density GIF, and cap_monomial_study.ipynb.

    PYTHONPATH=.../inq-stack/python /local/.../venv/bin/python3 build_monomial_report.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path
import csv

HERE = Path(__file__).resolve().parent
VAC = HERE.parent.parent
MONO = VAC / "cap_monomial"

def parse_kv(p):
    out = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(' ')
        try: out[k] = float(v)
        except ValueError: out[k] = v
    return out

recs = []
for d in sorted(MONO.glob('run_mono_*')):
    f = d / 'results/epsilon.txt'
    if f.exists():
        r = parse_kv(f); r['name'] = d.name; recs.append(r)
recs.sort(key=lambda r: (int(r.get('order', 0)), r.get('eta_Ha', 0)))
cols = ['name', 'order', 'eta_Ha', 'L_abs', 'E_eV', 'epsilon', 'absorbed_fraction']
with open(HERE / 'cap_monomial_combined.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in recs: w.writerow([r.get(c, '') for c in cols])
print(f'{len(recs)} monomial runs -> cap_monomial_combined.csv')

nb = new_notebook(); C = nb.cells
C.append(new_markdown_cell(r"""# Monomial CAP vs built-in sin² — thin absorber (L=5 Bohr), E=10 eV

Can a **monomial** complex absorbing potential (the inq-study ramp
$V=i\eta\,s^n$, $s\in[0,1]$ across the slab, peak absorption at the box wall) beat
the built-in **sin² hump** at a *thin* L=5 Bohr absorber — i.e. reach low reflection
without widening the box?

## Method
- `perturbations::absorbing_monomial` — NEW inq-study perturbation (inq/ immutable,
  untouched). Validated by `tests/monomial_shape_check`: it absorbs, and
  $\varepsilon(n{=}1)<\varepsilon(n{=}4)$ (the order signature of $s^n$).
- Benchmark: order $n\in\{1,2,3,4\}$ × depth $\eta\in\{-0.1,-0.2,-0.3,-0.5\}$ at
  fixed L=5 Bohr, E=10 eV (k₀≈0.857), ETRS, full free-WP observable set.
- Baselines: built-in sin² at L=5 (cap_thin_L5) and the sin² **width** it takes to
  reach the same ε (cap_Lopt_E10: L8→2.7%, L10→1.1%, L12→0.49%).

> **Provisional** until the inq-study engine regression (Task #7). Source:
> De Giovannini–Larsen–Rubio 2014 §IV; Riss & Meyer 1996; transmission-free CAP
> (Manolopoulos 2002) is the stretch target for ε→0 at short L.
"""))

C.append(new_code_cell(r"""import numpy as np, csv, re
from pathlib import Path
import matplotlib.pyplot as plt
from inqview.visualisation import style
style.apply_theme()
HERE=Path(r"%s"); VAC=Path(r"%s"); MONO=Path(r"%s")
def parse(p):
    o={}
    for ln in Path(p).read_text().splitlines():
        k,_,v=ln.partition(' ')
        try:o[k]=float(v)
        except ValueError:o[k]=v
    return o
mono=list(csv.DictReader(open(HERE/'cap_monomial_combined.csv')))
for r in mono:
    for k in ('order','eta_Ha','epsilon','absorbed_fraction'): r[k]=float(r[k])
orders=sorted({int(r['order']) for r in mono})
# sin2 baseline at L=5,E=10 (cap_thin_L5) and widen-L (cap_Lopt_E10)
def eps_of(globpat, key):
    out={}
    for d in (VAC).glob(globpat):
        f=d/'results/epsilon.txt'
        if f.exists():
            r=parse(f); out[round(r[key],2)]=r['epsilon']
    return out
sin2_L5=eps_of('cap_thin_L5/run_cap_k0.86_L5_eta*','eta_Ha')   # E=10 sin2 at L5 vs eta
widenL =eps_of('cap_Lopt_E10/run_cap_*','L_abs')               # best-ish sin2 vs L (mixed eta; take min per L below)
print('orders:',orders,'| sin2 L5 (eta->eps):',{k:round(v,3) for k,v in sin2_L5.items()})
"""%(str(HERE),str(VAC),str(MONO))))

C.append(new_markdown_cell(r"""## 1. Reflection ε vs depth, per monomial order — and the sin² floor
Each curve is one order n; the dashed line is the best built-in sin² at L=5. Does
any (n, η) drop below it?"""))
C.append(new_code_cell(r"""fig,ax=plt.subplots(figsize=(6.6,4.4))
for n in orders:
    c=sorted([r for r in mono if int(r['order'])==n],key=lambda r:abs(r['eta_Ha']))
    x=[abs(r['eta_Ha']) for r in c]; y=[r['epsilon'] for r in c]
    ax.semilogy(x,y,'o-',label=f'monomial n={n}')
if sin2_L5:
    best_sin2=min(sin2_L5.values())
    ax.axhline(best_sin2,ls='--',color='0.4',lw=1.4,label=f'sin² L=5 best ({best_sin2*100:.1f}%)')
ax.set_xlabel(r'depth $|\eta|$ (Ha)'); ax.set_ylabel(r'reflection $\varepsilon$')
ax.set_title('Monomial CAP at L=5, E=10 eV vs sin² floor'); ax.legend(fontsize=8); ax.grid(True,which='both',alpha=0.3)
fig.tight_layout(); fig.savefig('fig_monomial_vs_eta.png',dpi=140); plt.show()
best=min(mono,key=lambda r:r['epsilon'])
print(f"best monomial: eps={best['epsilon']:.3e} at n={int(best['order'])}, eta={best['eta_Ha']:+.2f}")
if sin2_L5: print(f"best sin2 L5 : eps={min(sin2_L5.values()):.3e}")
"""))

C.append(new_markdown_cell(r"""## 2. The width trade-off — monomial L=5 vs the L sin² needs
Where does the best monomial(L=5) land on the sin² width ladder? If best-monomial
≈ sin²(L=10), the ramp buys ~5 Bohr of effective width at the same physical L."""))
C.append(new_code_cell(r"""# sin2 min eps per L (across the eta we ran in cap_Lopt_E10)
perL={}
for d in (VAC).glob('cap_Lopt_E10/run_cap_*'):
    f=d/'results/epsilon.txt'
    if f.exists():
        r=parse(f); L=round(r['L_abs']); perL.setdefault(L,[]).append(r['epsilon'])
perL={L:min(v) for L,v in perL.items()}
fig,ax=plt.subplots(figsize=(6.2,4.0))
if perL:
    Ls=sorted(perL); ax.semilogy(Ls,[perL[L] for L in Ls],'s-',color='C2',label='sin² (best η per L)')
bm=min(mono,key=lambda r:r['epsilon'])['epsilon']
ax.axhline(bm,ls='--',color='C0',lw=1.5,label=f'best monomial @ L=5 ({bm*100:.2f}%)')
ax.set_xlabel('sin² absorber width L (Bohr)'); ax.set_ylabel(r'reflection $\varepsilon$')
ax.set_title('Best monomial(L=5) on the sin² width ladder (E=10 eV)'); ax.legend(fontsize=8); ax.grid(True,which='both',alpha=0.3)
fig.tight_layout(); fig.savefig('fig_monomial_width_equiv.png',dpi=140); plt.show()
print('sin2 best-eps per L:',{L:round(perL[L],4) for L in sorted(perL)})
"""))

C.append(new_markdown_cell(r"""## 3. Density GIF — best monomial run (WP meeting the ramp CAP)"""))
C.append(new_code_cell(r"""import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.animation as animation
sc=min(mono,key=lambda r:r['epsilon']); d=MONO/sc['name']
vti=d/'results/raw/vti/density_wp'
frames=sorted(vti.glob('density_wp_t*.vti'),key=lambda p:int(re.search(r'_t(\d+)',p.name).group(1)))
print(f"{len(frames)} frames from {sc['name']} (n={int(sc['order'])}, eta={sc['eta_Ha']:+.2f}, eps={sc['epsilon']:.3e})")
def zprof(p):
    rd=vtk.vtkXMLImageDataReader(); rd.SetFileName(str(p)); rd.Update(); img=rd.GetOutput()
    dm=img.GetDimensions(); a=vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(dm[2],dm[1],dm[0])
    sp=img.GetSpacing(); org=img.GetOrigin(); z=org[2]+sp[2]*np.arange(dm[2]); return z,a.sum(axis=(1,2))
if frames:
    z0,_=zprof(frames[0]); profs=[zprof(p)[1] for p in frames]; ymax=max(p.max() for p in profs)*1.05
    r=parse(d/'results/epsilon.txt')
    fig,ax=plt.subplots(figsize=(6,3.2)); (line,)=ax.plot(z0,profs[0])
    ax.axvspan(r['z_abs0'],r['z_abs0']+r['L_abs'],color='C3',alpha=0.15,label='monomial CAP (last 5 a₀)')
    ax.set_xlabel('z (Bohr)'); ax.set_ylabel('WP density (transverse sum)'); ax.set_ylim(0,ymax); ax.legend(loc='upper left',fontsize=8)
    ttl=ax.set_title('')
    def upd(k): line.set_ydata(profs[k]); ttl.set_text(f'frame {k+1}/{len(profs)}'); return line,ttl
    animation.FuncAnimation(fig,upd,frames=len(profs),interval=120,blit=False).save('fig_monomial_density.gif',writer='pillow',dpi=90); plt.close(fig)
    print('wrote fig_monomial_density.gif')
"""))

C.append(new_markdown_cell(r"""![gif](fig_monomial_density.gif)

## 4. Verdict (filled from the curves)
- Does any monomial (n, η) at L=5 beat sin²'s L=5 floor (~0.20)? By how much?
- What sin² width does the best monomial(L=5) equal? (effective-width gain)
- Every run carries the full free-WP minimum observable set (tiers 1–3 pass);
  energy-drift / norm-band = absorption signature.
- **Provisional** until Task #7. If the monomial gain is small, the
  transmission-free CAP (Manolopoulos 2002) is the next inq-study target.
"""))
C.append(new_code_cell(r"""for r in sorted(mono,key=lambda r:r['epsilon'])[:6]:
    print(f"n={int(r['order'])} eta={r['eta_Ha']:+.2f}  eps={r['epsilon']:.3e}  absorbed={r['absorbed_fraction']:.4f}  {r['name']}")
"""))

if __name__ == "__main__":
    ep = ExecutePreprocessor(timeout=2400, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': str(HERE)}})
    out = HERE / 'cap_monomial_study.ipynb'
    nbf.write(nb, out)
    print(f'wrote {out}')
