#!/usr/bin/env python3
"""build_run_evidence.py — per-hypothesis RUN-EVIDENCE notebooks for campaign_autorun.

One notebook per hypothesis (H0_runs … H5_runs) listing EVERY run in that sweep with
its raw evidence — run_summary config, converged/step-0 energy, interior density — so
each data point behind a study-notebook plot can be independently confirmed. The
single-run run-notebook assembler cannot span a sweep; this aggregator fills that gap.

Neutral by construction: it TABULATES what each run produced (numbers read from the run
files, never re-converged) and points at each run_summary.txt. No interpretation.

Run:  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
      /local/data/public/skcb2/tddft/venv/bin/python3 build_run_evidence.py
"""
from __future__ import annotations
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

REPO = "/local/data/public/skcb2/tddft"
LJ = f"{REPO}/ResearchProject/systems/localised_jellium"
CA = f"{LJ}/scripts/campaign_autorun"
OUT = Path(f"{LJ}/hypotheses/campaign_autorun_study/runs"); OUT.mkdir(parents=True, exist_ok=True)

# shared preamble: readers from analyse_phase (no emails on import) + pandas
PRE = f"""import sys, glob, csv, numpy as np, pandas as pd
sys.path.insert(0, {CA!r}); sys.path.insert(0, {REPO+'/inq-stack/python'!r})
from analyse_phase import e_total0, gs_energy, load_nz, _rs_present
from pathlib import Path
HA_EV=27.211386; CA=Path({CA!r}); RUNS=CA/'runs'
GS_P3=Path({LJ!r})/'scripts/h0_base_difference/gs/results'   # p3 bare-slab GS
def _egs(per): return gs_energy(GS_P3) if per==3 else gs_energy(RUNS/'h2/gs_p2_lz120/results')
def _summ(run_results):
    p=Path(run_results)/'run_summary.txt'
    return p.read_text() if p.exists() else '(missing)'
pd.set_option('display.max_rows', 200); pd.set_option('display.width', 160)
print('run-evidence — data root', RUNS)"""

# per-hypothesis evidence code. Each builds and DISPLAYS a DataFrame `ev` (one row per
# run) plus prints where every run_summary.txt lives.
EV = {
"H0": ("H0 — base WP-vs-classical E_total(0) gap",
       "Every WP and classical single-point run in the radius sweep (periodicity 3): "
       "its step-0 total energy and the excess over the bare-slab GS.",
"""base=RUNS/'h0'; EGS=_egs(3); rows=[]
for tag in ('wp','cl'):
    for r in _rs_present(base,tag,3):
        rd=base/f'{tag}_r{r}_p3'/'results'/f'{tag}_r{r}_p3'
        Et=e_total0(base/f'{tag}_r{r}_p3')
        rows.append(dict(run=f'{tag}_r{r}_p3', kind=tag, r_bohr=r,
                         E_tot0_Ha=round(Et,4), excess_eV=round((Et-EGS)*HA_EV,2),
                         summary=str(rd/'run_summary.txt')))
ev=pd.DataFrame(rows).sort_values(['kind','r_bohr']); print('E_GS(p3)=',round(EGS,4),'Ha'); ev"""),

"H0_p2": ("H0 (periodicity 2) — full measured energy decomposition",
       "Every WP and classical single-point run in the periodicity-2 radius sweep "
       "(re-run 2026-07-07), with the full streamed energy decomposition — total, "
       "external (E_ext=∫n·v_ext), Hartree, kinetic, xc — and the excess over the "
       "open-z GS. Each row's total is verified equal to the sum of its eight "
       "components (sum_minus_total column, ~1e-13 Ha).",
"""base=RUNS/'h0_p2'; EGS=_egs(2); rows=[]
import glob as _glob
KEYS=['total','kinetic','hartree','xc','external','nonlocal','ion','ion_kinetic','exact_exchange']
def _comp(rd):
    f=next(iter(_glob.glob(str(rd)+'/**/observables.csv',recursive=True)))
    rr=list(csv.reader(open(f))); h,d=rr[0],rr[1]; g=lambda c: float(d[h.index(c)])
    return {k:g('energy_'+k) for k in KEYS}
for tag in ('wp','cl'):
    for r in (4,12,20,28,36,40):
        rd=base/f'{tag}_r{r}_p2'
        if not list(rd.glob('**/observables.csv')): continue
        c=_comp(rd); s=sum(c[k] for k in ['kinetic','external','nonlocal','hartree','xc','exact_exchange','ion','ion_kinetic'])
        rows.append(dict(run=f'{tag}_r{r}_p2', kind=tag, r_bohr=r,
                         E_tot0_Ha=round(c['total'],4), excess_eV=round((c['total']-EGS)*HA_EV,1),
                         E_ext_Ha=round(c['external'],3), U_H_Ha=round(c['hartree'],3),
                         T_Ha=round(c['kinetic'],3), E_xc_Ha=round(c['xc'],3),
                         sum_minus_total=f'{abs(s-c["total"]):.1e}',
                         summary=str(next(iter(rd.glob('**/run_summary.txt'))))))
ev=pd.DataFrame(rows).sort_values(['kind','r_bohr']); print('E_GS(p2)=',round(EGS,4),'Ha'); ev"""),

"H1": ("H1 — edge model (Gibbs vs Friedel)",
       "Every ground-state run in the edge-width sweep: converged energy and the "
       "interior planar density n0 (mean of n(z) for |z|<6).",
"""base=RUNS/'h1'; rows=[]
for p in sorted(base.glob('gs_w*')):
    w=float(p.name.split('_w')[1]); res=p/'results'
    z,nz=load_nz(res); n0=float(nz[np.abs(z)<6].mean())
    rows.append(dict(run=p.name, edge_w=w, E_GS_Ha=round(gs_energy(res),4),
                     interior_n0=f'{n0:.3e}', summary=str(res/'run_summary.txt')))
ev=pd.DataFrame(rows).sort_values('edge_w'); ev"""),

"H2": ("H2 — GS convergence + open-z viability",
       "Every ground-state run in the box-length sweep (periodicity 3) plus the two "
       "open-z (periodicity 2) runs: converged energy and interior n0.",
"""base=RUNS/'h2'; rows=[]
for p in sorted(base.glob('gs_lz*'))+sorted(base.glob('gs_p2_lz*')):
    res=p/'results'; z,nz=load_nz(res); n0=float(nz[np.abs(z)<6].mean())
    per=2 if 'p2' in p.name else 3
    lz=int(p.name.split('_lz')[1])
    rows.append(dict(run=p.name, periodicity=per, Lz_bohr=lz,
                     E_GS_Ha=round(gs_energy(res),4), interior_n0=f'{n0:.3e}',
                     summary=str(res/'run_summary.txt')))
ev=pd.DataFrame(rows).sort_values(['periodicity','Lz_bohr']); ev"""),

"H3": ("H3 — surface energetics (thickness)",
       "Every ground-state run in the slab-thickness sweep (half-width a, electron "
       "count N scaled to hold n0): converged energy and interior n0.",
"""base=RUNS/'h3'; rows=[]
for p in sorted(base.glob('gs_a*_N*')):
    a=float(p.name.split('_a')[1].split('_N')[0]); N=int(p.name.split('_N')[1]); res=p/'results'
    z,nz=load_nz(res); n0=float(nz[np.abs(z)<a-2].mean()) if a>3 else float('nan')
    rows.append(dict(run=p.name, half_width_a=a, N_elec=N,
                     E_GS_Ha=round(gs_energy(res),4), interior_n0=f'{n0:.3e}',
                     summary=str(res/'run_summary.txt')))
ev=pd.DataFrame(rows).sort_values('half_width_a'); ev"""),

"H4": ("H4 — WP energetics: E_SIE + PBC-vs-open-z",
       "Every stationary-WP single-point run in the radius sweep at periodicity 2 and "
       "3: step-0 energy, excess over the GS (per BC), and excess minus the 81.6 eV "
       "zero-point (the E_SIE-ish residual).",
"""base=RUNS/'h4'; ZP=3/(4*0.5**2)*HA_EV; rows=[]
for per in (2,3):
    EGS=_egs(per)
    for r in _rs_present(base,'wp',per):
        rd=base/f'wp_r{r}_p{per}'; Et=e_total0(rd)
        exc=(Et-EGS)*HA_EV
        rows.append(dict(run=f'wp_r{r}_p{per}', periodicity=per, r_bohr=r,
                         E_tot0_Ha=round(Et,4), excess_eV=round(exc,2),
                         excess_minus_ZP_eV=round(exc-ZP,2),
                         summary=str(rd/'results'/f'wp_r{r}_p{per}'/'run_summary.txt')))
ev=pd.DataFrame(rows).sort_values(['periodicity','r_bohr']); print('ZP=3/4σ²=',round(ZP,1),'eV'); ev"""),

"H5": ("H5 — classical mirror (route 2 + thread D)",
       "Every classical-ghost single-point run in the radius sweep at periodicity 2 "
       "and 3: step-0 energy and excess over the GS (per BC).",
"""base=RUNS/'h5'; rows=[]
for per in (2,3):
    EGS=_egs(per)
    for r in _rs_present(base,'cl',per):
        rd=base/f'cl_r{r}_p{per}'; Et=e_total0(rd)
        rows.append(dict(run=f'cl_r{r}_p{per}', periodicity=per, r_bohr=r,
                         E_tot0_Ha=round(Et,4), excess_eV=round((Et-EGS)*HA_EV,2),
                         summary=str(rd/'results'/f'cl_r{r}_p{per}'/'run_summary.txt')))
ev=pd.DataFrame(rows).sort_values(['periodicity','r_bohr']); ev"""),
}


def nb(cells):
    n = new_notebook(); n.cells = cells
    n.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
    return n


def build(key):
    title, what, code = EV[key]
    cells = [
        new_markdown_cell(
            f"# {title} — run evidence\n\n*Per-run raw evidence for every run behind the "
            f"{key} study plot. Numbers read from the run files (never re-converged); open "
            f"the linked `run_summary.txt` to confirm any row. Auto-built by "
            f"`build_run_evidence.py`.*"),
        new_markdown_cell(f"## What this table shows\n{what}\n\nRun data: "
                          f"`scripts/campaign_autorun/runs/{key.lower()}/` — one row per run, "
                          f"with the absolute path to each run's `run_summary.txt` in the last "
                          f"column."),
        new_code_cell(PRE),
        new_markdown_cell("## Evidence table (one row per run)"),
        new_code_cell(code),
        new_markdown_cell("## Confirm any data point\nEach row's `summary` column is the "
                          "absolute path to that run's `run_summary.txt` (full config: cell, "
                          "grid, electrons, background, propagator). The energy columns are "
                          "read from that run's `results/**/observables.csv` (step 0) or the "
                          "GS `run_summary.txt` — nothing here is recomputed by a new SCF."),
    ]
    p = OUT / f"{key}_runs.ipynb"
    nbf.write(nb(cells), str(p)); print("wrote", p.relative_to(OUT.parent)); return p


if __name__ == "__main__":
    for k in ("H0", "H0_p2", "H1", "H2", "H3", "H4", "H5"):
        build(k)
    print("done — execute: python3 -m nbconvert --to notebook --execute --inplace runs/H*_runs.ipynb (venv)")
