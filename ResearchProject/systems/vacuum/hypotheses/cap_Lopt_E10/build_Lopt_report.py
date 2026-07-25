#!/usr/bin/env python3
"""Build + execute the widen-L sin² CAP baseline notebook (ADR 0007 location).

The "how thin can the built-in sin² CAP go at E=10 eV?" baseline: fix E=10 eV
(k0≈0.857) and 2D-sweep absorber width L × depth η with perturbations::absorbing,
to map the reflection error ε(L, η) and locate the smallest L (and best η) that
reaches a low ε. This is the run-SET that carries the "if you relax thinness you
reach ~1%" half of the thin-absorber decision (vs cap_thin_L5 sin² @ L=5 and
cap_monomial @ L=5). Writes the combined CSV, ε(L) curves per η, ε(η) curves per
L, and cap_Lopt_E10_study.ipynb.

    PYTHONPATH=.../inq-stack/python /local/.../venv/bin/python3 build_Lopt_report.py

ε PROVISIONAL until the inq-study engine regression (Task #7).
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path
import csv

HERE = Path(__file__).resolve().parent
VAC = HERE.parent.parent
LOPT = VAC / "cap_Lopt_E10"


def parse_kv(p):
    out = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(' ')
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out


recs = []
for d in sorted(LOPT.glob('run_cap_*')):
    f = d / 'results/epsilon.txt'
    if f.exists():
        r = parse_kv(f)
        r['name'] = d.name
        recs.append(r)
recs.sort(key=lambda r: (r.get('L_abs', 0), r.get('eta_Ha', 0)))
cols = ['name', 'L_abs', 'eta_Ha', 'E_eV', 'k0', 'epsilon', 'absorbed_fraction',
        'N0', 'tau', 'N_STEPS', 'dt', 'dx', 'Lcell_z', 'width_frac', 'mid_frac']
with open(HERE / 'cap_Lopt_combined.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(cols)
    for r in recs:
        w.writerow([r.get(c, '') for c in cols])
print(f'{len(recs)} widen-L runs -> cap_Lopt_combined.csv')

nb = new_notebook()
C = nb.cells

C.append(new_markdown_cell(r"""# Built-in sin² CAP — how thin can we go at E = 10 eV?  (widen-L baseline)

**Question.** For a single-electron Gaussian wavepacket leaving the box at
E = 10 eV, what absorber **width L** does the *built-in* INQ complex absorbing
potential (CAP) need to suppress the boundary **reflection error** ε to the
low-percent level — and what depth η is optimal? This is the reference baseline
against which the *thin* (L = 5 Bohr) absorbers are judged:

| sweep | absorber | thin? | best ε @ 10 eV |
|---|---|---|---|
| `cap_thin_L5` | sin² hump, L = 5 | ✓ | ~20.9 % (floor) |
| `cap_monomial` | inq-study ramp $i\eta s^n$, L = 5 | ✓ | 8.3 % (n=1, η=−0.5) |
| **`cap_Lopt_E10`** *(this study)* | **sin² hump, L = 6…15** | ✗ | **see below** |

So this notebook answers: *if we are allowed to spend box width, how cheaply does
the stock absorber reach ε ≈ 1 %?*

## The absorber — built-in `perturbations::absorbing`

INQ's only built-in CAP adds an imaginary (non-Hermitian) potential over a slab of
width $L$ centred at fractional position `mid_pos`, hard-coded to a **sin² hump**:

$$ V_\mathrm{cap}(z) \;=\; i\,\eta\,\sin^2\!\Big(\pi\,\frac{z-z_\mathrm{abs0}}{L}\Big),
\qquad z\in[z_\mathrm{abs0},\,z_\mathrm{abs0}+L] $$

| symbol | meaning |
|---|---|
| $\eta$ | absorber **depth** (`amplitude`, in Ha); $\eta<0$ ⇒ absorbing. Swept ∈ {−0.3, −0.5, −1.0} |
| $L$ | absorber **width** (Bohr), the slab thickness. Swept ∈ {6, 8, 10, 12, 15} |
| $z_\mathrm{abs0}$ | inner edge of the absorber slab; the **inner region** is $z<z_\mathrm{abs0}$ |
| sin² | the (only) hard-coded shape: zero at *both* slab edges, peak at the centre |

Absorption is intrinsically non-Hermitian, so it **requires the ETRS
propagator** — Crank–Nicolson renormalises the wavefunction each step and silently
cancels the absorption.

## Reflection error ε — the figure of merit

$$ \varepsilon \;=\; \frac{\displaystyle\int_{z<z_\mathrm{abs0}}|\psi(z,\tau)|^2\,dz}{N_0}
\;=\; \frac{\text{WP norm left in the inner region at }\tau}{\text{initial norm }N_0} $$

| symbol | meaning |
|---|---|
| $\tau$ | stop time; chosen so a *perfectly* absorbed packet has fully exited the inner region |
| $N_0$ | initial wavepacket norm (≈ 0.041 here, the single electron's projection on the slab) |
| $\varepsilon$ | fraction **not** absorbed — what reflects/lingers. **Lower = better.** ε→0 is the goal |

ε is a direct reflectivity proxy: ε = 0.21 means 21 % of the packet was reflected
back into the physical region by the absorber.

> **Provisional.** All ε here are PROVISIONAL until the inq-study engine regression
> (Task #7) confirms the scalar-potential complexification. Source for the CAP
> method and the sin² form: De Giovannini, Larsen & Rubio, *Eur. Phys. J. B* **88**,
> 56 (2014), §IV; Riss & Meyer, *J. Phys. B* **26**, 4503 (1993). The
> transmission-free CAP (Manolopoulos, *J. Chem. Phys.* **117**, 9552 (2002)) is the
> stretch target for ε→0 at *short* L.

## Files that produced this study

| role | path |
|---|---|
| run binary (built-in sin² CAP, env-parameterised) | `ResearchProject/systems/vacuum/scripts/cap_sweep/run.cpp` |
| dispatcher (15 runs, 5 L × 3 η, GPUs 0–1) | `ResearchProject/systems/vacuum/scripts/cap_Lopt_E10/dispatch.py` |
| this builder | `ResearchProject/systems/vacuum/hypotheses/cap_Lopt_E10/build_Lopt_report.py` |
| per-run provenance | `cap_Lopt_E10/run_cap_*/results/epsilon.txt` |
"""))

C.append(new_code_cell(r"""import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from inqview.visualisation import style
style.apply_theme()

HERE = Path(r"%s"); VAC = Path(r"%s"); LOPT = Path(r"%s")
rows = list(csv.DictReader(open(HERE / 'cap_Lopt_combined.csv')))
for r in rows:
    for k in ('L_abs', 'eta_Ha', 'epsilon', 'absorbed_fraction', 'E_eV', 'k0'):
        r[k] = float(r[k])
Ls   = sorted({r['L_abs']  for r in rows})
etas = sorted({r['eta_Ha'] for r in rows}, reverse=True)   # -0.3, -0.5, -1.0
def eps(L, eta):
    for r in rows:
        if r['L_abs'] == L and r['eta_Ha'] == eta:
            return r['epsilon']
    return np.nan
print(f"{len(rows)} runs | L = {Ls} Bohr | eta = {etas} Ha | E = {rows[0]['E_eV']} eV")
""" % (HERE.as_posix(), VAC.as_posix(), LOPT.as_posix())))

C.append(new_markdown_cell(r"""## Result 1 — ε vs absorber width L (one curve per depth η)

The headline relationship: how fast the reflection error falls as the absorber is
widened, for each depth. The dashed lines mark the two **thin (L = 5)** competitors."""))

C.append(new_code_cell(r"""fig, ax = plt.subplots(figsize=(7.0, 4.6))
for eta in etas:
    y = [eps(L, eta) for L in Ls]
    ax.plot(Ls, [v*100 for v in y], 'o-', label=fr'$\eta={eta:.1f}$ Ha')
# thin L=5 competitors
ax.axhline(20.9, ls='--', lw=1.2, color='0.45', label='sin² L=5 floor (20.9%)')
ax.axhline(8.3,  ls=':',  lw=1.4, color='C3',   label='monomial n=1 L=5 (8.3%)')
ax.set_yscale('log')
ax.set_xlabel('absorber width  L  (Bohr)')
ax.set_ylabel(r'reflection error  $\varepsilon$  (%)')
ax.set_title('Built-in sin² CAP at E = 10 eV: ε falls with width')
ax.legend(fontsize=8, framealpha=0.9)
fig.tight_layout(); fig.savefig(HERE / 'fig_Lopt_eps_vs_L.png', dpi=140)
print('best (lowest ε):',
      min(rows, key=lambda r: r['epsilon'])['name'],
      f"= {min(r['epsilon'] for r in rows)*100:.2f}%")
"""))

C.append(new_markdown_cell(r"""## Result 2 — the depth optimum (ε vs η at each L)

Deeper is **not** always better: too large an |η| makes the absorber onset itself
reflect. The sweep brackets an optimum near η ≈ −0.5."""))

C.append(new_code_cell(r"""fig, ax = plt.subplots(figsize=(7.0, 4.6))
for L in Ls:
    y = [eps(L, eta) for eta in etas]
    ax.plot([abs(e) for e in etas], [v*100 for v in y], 's-', label=fr'$L={int(L)}$')
ax.set_yscale('log')
ax.set_xlabel(r'absorber depth  $|\eta|$  (Ha)')
ax.set_ylabel(r'reflection error  $\varepsilon$  (%)')
ax.set_title('Depth optimum: η ≈ −0.5 beats η = −1.0 at every width')
ax.legend(fontsize=8, ncol=2, title='width L (Bohr)')
fig.tight_layout(); fig.savefig(HERE / 'fig_Lopt_eps_vs_eta.png', dpi=140)
"""))

C.append(new_markdown_cell(r"""## Summary table — ε(L, η) in percent"""))

C.append(new_code_cell(r"""hdr = 'L\\η  ' + '  '.join(f'{eta:+.1f}' for eta in etas)
print(hdr); print('-'*len(hdr))
for L in Ls:
    print(f'{int(L):<4} ' + '  '.join(f'{eps(L,eta)*100:5.2f}' for eta in etas))
print('\nunits: percent reflection error ε; E = 10 eV; built-in sin² CAP (ETRS).')
"""))

C.append(new_markdown_cell(r"""## Takeaway

- The built-in sin² CAP reaches the **low-percent** regime only by spending box
  width: at the optimal depth **η = −0.5**, ε ≈ **2.7 % (L=8) → 1.1 % (L=10) →
  0.49 % (L=12) → 0.15 % (L=15)**.
- **η = −0.5 is optimal**, beating both the weaker η = −0.3 (under-absorbs) and the
  deeper η = −1.0 (onset reflection) at every width.
- Against the **thin (L=5)** competitors, no width here is L=5 — the stock sin²
  floors at ~21 % there, and the inq-study monomial ramp gets to 8.3 %. Reaching
  ε ≈ 1 % currently **costs ≥ L≈10 Bohr** of absorber. Whether thinness can be
  bought back at ε ≈ 1 % is the open question the transmission-free CAP targets.

*All ε PROVISIONAL until Task #7.*"""))

ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
ep.preprocess(nb, {'metadata': {'path': str(HERE)}})
out = HERE / 'cap_Lopt_E10_study.ipynb'
with open(out, 'w') as fh:
    nbf.write(nb, fh)
print('wrote', out)
