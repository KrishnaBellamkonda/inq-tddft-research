#!/usr/bin/env python3
"""twin_notebook_builder — assemble the executed twin-run ANALYSIS notebook.

Builds a self-contained `.ipynb` that runs the deterministic engine
(`twin_decompose`) on a twin pair and lays out: provenance, the findings table,
the per-step decomposition plot, and the residual/SIE — with markdown narrative
cells the agent fills using the interpretation rules in SKILL.md.

The generated notebook lives in the run-set's `hypotheses/<sweep>/` folder
(ADR 0007), never inside the skill. Figures are PNG (project rule).

Usage:
  twin_notebook_builder.py --pair DIR         --out study.ipynb
  twin_notebook_builder.py --wp DIR --classical DIR --out study.ipynb [--title "..."]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

SKILL_DIR = str(Path(__file__).resolve().parent)


def _cells(wp_dir: str, cl_dir: str, title: str) -> list:
    return [
        new_markdown_cell(
            f"# {title}\n\n"
            "**Twin-run energy decomposition** — a matched classical (Gaussian-charge "
            "perturbation) vs wavepacket pair, identical except the projectile "
            "representation. Every difference `d(·) = WP − classical` is a *quantum "
            "effect* of the wavepacket treatment.\n\n"
            "Deterministic numbers come from `twin_decompose.py`; the narrative "
            "(marked _Narrative_) is filled from the interpretation rules in the "
            "`twin-run-analysis` skill."),
        new_code_cell(
            f"import sys; sys.path.insert(0, {SKILL_DIR!r})\n"
            "import numpy as np, pandas as pd, matplotlib.pyplot as plt\n"
            "import twin_decompose as td\n"
            "try:\n"
            "    from inqview.visualisation import style as _s; _s.apply()\n"
            "except Exception:\n"
            "    plt.rcParams.update({'figure.dpi': 120, 'font.size': 10})\n"
            f"WP_DIR = {wp_dir!r}\nCL_DIR = {cl_dir!r}\n"
            "res = td.decompose(WP_DIR, CL_DIR)\n"
            "print(res.report())"),
        new_markdown_cell(
            "## 0. Visual intuition — the density matrix (animated xz slices)\n\n"
            "The full 2D electron density in the propagation **x–z plane** (mid-y slice), animated "
            "over time — the *density matrix* the run produces. A 3×3 grid:\n\n"
            "- **columns** — `density` n(x,z,t), `induced` Δn=n(t)−n(0), `instantaneous` "
            "Δn=n(t)−n(t−Δt);\n"
            "- **rows** — `classical`, `wavepacket`, and `WP−classical`.\n\n"
            "Each tile is a LINEAR | LOG (or ±symlog) panel pair; dashed lines mark the slab faces "
            "(|z|=12.5). Rendered by `inqview.visualisation.make_twin_density_matrix` — mid-y xz "
            "slice, physical-order VTIs (never fftshift'd, per the vti-coordinate-mapping rule). "
            "Classical and WP rows share one colour scale per column; the WP−classical row owns "
            "its own symmetric scale."),
        new_code_cell(
            "from inqview.visualisation import make_twin_density_matrix\n"
            "from IPython.display import HTML, display, Markdown\n"
            "import os\n"
            "DT = float(td.parse_summary(CL_DIR + '/run_summary.txt').get('dt', 0.05))\n"
            "mgifs = make_twin_density_matrix(CL_DIR, WP_DIR, '.', dt=DT, slab_face=12.5,\n"
            "                                 cap_inner=None, frames_max=30, fps=8)\n"
            "if mgifs:\n"
            "    by = {(r, c): os.path.basename(p) for r, c, p, _ in mgifs}\n"
            "    cols = ['density', 'induced', 'instantaneous']\n"
            "    rows = ['classical', 'wp', 'wp_minus_cl']\n"
            "    rl = {'classical': 'Classical', 'wp': 'Wavepacket', 'wp_minus_cl': 'WP − classical'}\n"
            "    cl_ = {'density': 'density  n(x,z,t)', 'induced': 'induced  \\u0394n=n(t)\\u2212n(0)',\n"
            "           'instantaneous': 'instantaneous  \\u0394n=n(t)\\u2212n(t\\u2212\\u0394t)'}\n"
            "    h = \"<table><tr><td></td>\" + ''.join(f\"<td align='center'><b>{cl_[c]}</b></td>\" for c in cols) + '</tr>'\n"
            "    for r in rows:\n"
            "        h += f\"<tr><td valign='middle'><b>{rl[r]}</b></td>\"\n"
            "        for c in cols:\n"
            "            fn = by.get((r, c))\n"
            "            h += f\"<td><img src='{fn}' width='340'></td>\" if fn else '<td>\\u2014</td>'\n"
            "        h += '</tr>'\n"
            "    h += '</table>'\n"
            "    display(HTML(h))\n"
            "    print('wrote', len(mgifs), 'density-matrix GIFs:', sorted(by.values()))\n"
            "else:\n"
            "    display(Markdown('_No density frames saved for this pair \\u2014 density matrix unavailable._'))"),
        new_markdown_cell(
            "_Narrative — first read._ Read the matrix by column and row. **density** row: does the "
            "WP track the classical charge, or spread / reflect / split / tunnel? **induced** "
            "(Δn vs t=0): where does charge pile up or deplete — the wake, the reflected lobe, the "
            "captured cloud? **instantaneous** (frame-to-frame): the live flux, sharpest at the "
            "surface crossing. The **WP−classical** row isolates the pure quantum difference in "
            "each; the decomposition (§1b onward) attributes it term-by-term."),
        new_markdown_cell(
            "### Pairwise-energy evolution (animated)\n\n"
            "The Δ(WP−classical) pairwise Coulomb bars over time — the energy companion to the "
            "density matrix."),
        new_code_cell(
            "import imageio, io\n"
            "if res.pairwise is not None and len(res.pairwise) > 2:\n"
            "    P = res.pairwise; terms = ['e_ss', 'e_pp', 'e_ps', 'e_sb', 'e_pb']\n"
            "    lo = min(P['d_'+t].min() for t in terms); hi = max(P['d_'+t].max() for t in terms)\n"
            "    idx = P.step.values[::max(1, len(P)//40)]; pf = []\n"
            "    for st in idx:\n"
            "        r = P[P.step == st].iloc[0]\n"
            "        fig, ax = plt.subplots(figsize=(4.6, 3.2))\n"
            "        ax.bar(terms, [r['d_'+t] for t in terms], color='C0'); ax.axhline(0, color='k', lw=.5)\n"
            "        ax.set_ylim((lo-0.05*abs(lo)-1e-6)*1.05, (hi+0.05*abs(hi)+1e-6)*1.05)\n"
            "        ax.set_title(f'step {int(st)}   \\u0394(WP\\u2212cl)  eV'); fig.tight_layout()\n"
            "        buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=72); plt.close(fig); buf.seek(0)\n"
            "        pf.append(imageio.imread(buf))\n"
            "    imageio.mimsave('pairwise_evolution.gif', pf, duration=0.15)\n"
            "    display(HTML(\"<img src='pairwise_evolution.gif' width='430'>\"))\n"
            "else:\n"
            "    display(Markdown('_Too few steps for a pairwise GIF._'))"),
        new_markdown_cell("## 1. Provenance & parity"),
        new_code_cell(
            "print('sigma_WP =', res.sigma_wp, ' k0 =', res.k0)\n"
            "print(res.parity.as_text())\n"
            "assert res.parity.ok, 'NOT a valid twin pair — fix generation before analysing'"),
        new_markdown_cell(
            "## 1b. Energy accounting at the ground state (t=0)\n\n"
            "Before any dynamics: confirm the twins are consistent at step 0. At t=0 the "
            "classical Gaussian *charge* equals the WP *density* by construction, so **all six "
            "pairwise Coulomb terms must be ~equal** (the slab/background terms are the gauge "
            "test — exactly 0). INQ's total energy is charged-cell-convention-dependent in "
            "absolute value, but the **twin difference must close exactly**: "
            "`ΔE_total = ΔE_kin + ΔE_hartree + ΔE_external + ΔE_xc`, with `ΔE_kin` equal to the "
            "analytic WP kinetic surplus `½k₀² + 3/(4σ²)`. Finally the residual "
            "(= WP self-Hartree) and the one-electron **SIE** are read off. The scorecard "
            "summarises all of this in one table."),
        new_code_cell(
            "clo,_ = td.load_run(CL_DIR); wpo,_ = td.load_run(WP_DIR)\n"
            "c0, w0 = clo.iloc[0], wpo.iloc[0]\n"
            "HA = td.HA_EV\n"
            "def comp(col): return float(c0[col]), float(w0[col]), float(w0[col]-c0[col])\n"
            "# (a) pairwise Coulomb terms at t=0 (all 6; slab/bg are the gauge test)\n"
            "pw0 = res.pairwise_table(step=0)\n"
            "# (b) INQ E_total component accounting (Hartree)\n"
            "acct = pd.DataFrame([dict(quantity=lbl, classical=round(comp(col)[0],4),\n"
            "        WP=round(comp(col)[1],4), delta_WP_minus_cl=round(comp(col)[2],4))\n"
            "    for lbl,col in [('E_total (Ha)','energy_total'),('E_kinetic (Ha)','energy_kinetic'),\n"
            "        ('E_hartree (Ha)','energy_hartree'),('E_external (Ha)','energy_external'),\n"
            "        ('E_xc (Ha)','energy_xc')]])\n"
            "# checks\n"
            "dtot = comp('energy_total')[2]\n"
            "sumd = sum(comp(c)[2] for c in ['energy_kinetic','energy_hartree','energy_external','energy_xc'])\n"
            "dKin = comp('energy_kinetic')[2]; kin_an = 0.5*res.k0**2 + 0.75/res.sigma_wp**2\n"
            "R   = [f for f in res.findings if f['term'].startswith('residual')][0]['value_ev']\n"
            "SIE = [f for f in res.findings if f['term'].startswith('SIE')][0]['value_ev']\n"
            "# closure of E_H,E_ext from the pairwise terms (Poisson linearity)\n"
            "p0 = res.pairwise.iloc[0]\n"
            "cl_H=p0.e_ss_cl/HA; cl_x=(p0.e_sb_cl+p0.e_ps_cl)/HA\n"
            "wp_H=(p0.e_ss_wp+p0.e_ps_wp+p0.e_pp_wp)/HA; wp_x=(p0.e_sb_wp+p0.e_pb_wp)/HA\n"
            "clo_err=max(abs(cl_H-c0.energy_hartree),abs(cl_x-c0.energy_external))\n"
            "wpo_err=max(abs(wp_H-w0.energy_hartree),abs(wp_x-w0.energy_external))\n"
            "gmax=res.gauge['max_invariant_delta_ev']; pmax=max(abs(pw0.delta_wp_minus_cl))\n"
            "scorecard = pd.DataFrame([\n"
            "  dict(check='pairwise slab/bg  Δ≈0  (no inter-run gauge)', value=f'{gmax:.2e} eV', verdict='PASS' if res.gauge['no_gauge'] else 'FAIL'),\n"
            "  dict(check='pairwise projectile terms  Δ  (t=0 density match)', value=f'max {pmax:.3f} eV', verdict='PASS' if pmax<1.0 else 'LARGE'),\n"
            "  dict(check='ΔE_total = ΣΔ(components)  (energy closure)', value=f'{dtot-sumd:.1e} Ha', verdict='PASS' if abs(dtot-sumd)<1e-6 else 'FAIL'),\n"
            "  dict(check='ΔE_kin = ½k₀²+3/4σ²  (WP kinetic surplus)', value=f'{dKin:.3f} vs {kin_an:.3f} Ha', verdict='PASS' if abs(dKin-kin_an)<1e-2 else 'CHECK'),\n"
            "  dict(check='E_H,E_ext reconstruct from pairwise (Poisson)', value=f'err {max(clo_err,wpo_err):.1e} Ha', verdict='PASS' if max(clo_err,wpo_err)<1e-6 else 'FAIL'),\n"
            "  dict(check='residual R = Δ(E_H+E_ext) − U_proj_bg = WP self-Hartree', value=f'{R:.2f} eV', verdict='—'),\n"
            "  dict(check='SIE = R + ΔE_xc  (LDA one-electron self-interaction)', value=f'{SIE:.2f} eV', verdict='—'),\n"
            "])\n"
            "print('=== pairwise Coulomb terms at t=0 (eV) ==='); display(pw0)\n"
            "print('=== INQ E_total component accounting (Hartree) ==='); display(acct)\n"
            "print('=== ground-state energy-accounting scorecard ==='); display(scorecard)"),
        new_markdown_cell(
            "_Narrative — accounting._ Read the scorecard top-to-bottom: the slab/background "
            "terms are gauge-free (Δ=0), the total-energy difference closes into its components, "
            "the WP kinetic surplus matches `½k₀²+3/(4σ²)` to the analytic value, and the pairwise "
            "terms reconstruct E_hartree/E_external exactly — so the residual is a clean WP "
            "self-Hartree and the SIE is the physical one-electron self-interaction. Any FAIL here "
            "invalidates the downstream physics and must be fixed before interpreting dynamics."),
        new_markdown_cell(
            "## 2. Findings — the decomposition\n\n"
            "`value` = measured `d(·)`; `expected` = the known attribution; "
            "`unexplained` = the remainder to interpret."),
        new_code_cell(
            "ft = res.findings_table()\n"
            "ft_disp = ft.copy()\n"
            "for c in ['value_ev','expected_ev','unexplained_ev']:\n"
            "    ft_disp[c] = ft_disp[c].astype(float).round(2)\n"
            "ft_disp[['term','value_ev','expected_ev','unexplained_ev','interpretation']]"),
        new_markdown_cell(
            "_Narrative — the physics._ Fill from the interpretation rules:\n"
            "- **dKin** → WP localisation zero-point `3/(4σ²)` (+ `k0²/2`).\n"
            "- **residual R** → WP self-Hartree `E_H[WP–WP]`; the `unexplained` gap "
            "vs the free-space reference is the open-z gauge (~0.9 eV), not missing physics.\n"
            "- **SIE = R + dXC** → LDA one-electron self-interaction error (the irreducible residue).\n"
            "- **Gauge caveat**: interpret only `d(E_H+E_ext)`, never `dHartree`/`dExt` alone."),
        new_markdown_cell(
            "## 3. Pairwise Coulomb decomposition — every Δ physically attributable\n\n"
            "The lumped `E_hartree`/`E_external` resolved into the P/S/B pairwise terms; "
            "classical vs WP with Δ. Physically-identical terms (E_SS, E_SB, E_BB) must "
            "have Δ≈0 → the gauge test."),
        new_code_cell(
            "pt = res.pairwise_table(step=1) if res.pairwise is not None else None\n"
            "if pt is not None:\n"
            "    print('gauge test:', res.gauge)\n"
            "    if res.gauge and not res.gauge['no_gauge']:\n"
            "        print('WARNING: gauge present — gauge-correct before interpreting Δ')\n"
            "    p = res.pairwise\n"
            "    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3), constrained_layout=True)\n"
            "    for t in ['e_ss','e_sb','e_bb']:\n"
            "        ax[0].plot(p.step, p['d_'+t], marker='.', label='Δ'+t)\n"
            "    ax[0].set_title('gauge-invariant Δ (must stay ~0)'); ax[0].axhline(0,color='k',lw=.5)\n"
            "    ax[0].set_xlabel('step'); ax[0].set_ylabel('Δ (eV)'); ax[0].legend(fontsize=7)\n"
            "    for t in ['e_pp','e_ps','e_pb']:\n"
            "        ax[1].plot(p.step, p['d_'+t], marker='.', label='Δ'+t)\n"
            "    ax[1].set_title('projectile Δ (the physical differences)')\n"
            "    ax[1].set_xlabel('step'); ax[1].set_ylabel('Δ (eV)'); ax[1].legend(fontsize=7)\n"
            "    fig.savefig('twin_pairwise.png', dpi=150)\n"
            "    display(pt)\n"
            "else:\n"
            "    print('no interactions.csv — pairwise decomposition unavailable')"),
        new_markdown_cell(
            "_Narrative._ With the gauge test passing (ΔE_SS=ΔE_SB=ΔE_BB≈0), every "
            "non-zero Δ is physical — attribute each: ΔE_PP → WP self-Hartree vs "
            "dispersion (∝1/σ); ΔE_PS → projectile-slab shape/polarisation; ΔE_PB → "
            "projectile-background. The kinetic surplus (ΔKE) is WP localisation + motion."),
        new_markdown_cell("## 4. Per-step decomposition (twin differences)"),
        new_code_cell(
            "s = res.steps\n"
            "fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)\n"
            "for col, lab in [('dKin','dKin'),('dXC','dXC'),('d_H_ext','d(E_H+E_ext)'),\n"
            "                 ('residual','residual R'),('sie','SIE')]:\n"
            "    ax[0].plot(s.step, s[col], marker='o', ms=3, label=lab)\n"
            "ax[0].set_xlabel('step'); ax[0].set_ylabel('energy diff (eV)')\n"
            "ax[0].legend(fontsize=7); ax[0].set_title('difference terms vs step')\n"
            "ax[1].plot(s.step, s.residual, 'o-', label='residual R')\n"
            "ax[1].plot(s.step, s.sie, 's-', label='SIE')\n"
            "ax[1].set_xlabel('step'); ax[1].set_ylabel('eV'); ax[1].legend(fontsize=8)\n"
            "ax[1].set_title('WP self-Hartree & SIE')\n"
            "fig.savefig('twin_decomposition.png', dpi=150)\n"
            "print('drift:', {k: round(v['max_step_change_ev'],4) for k,v in res.drift.items()})"),
        new_markdown_cell(
            "_Narrative — dynamics._ For an at-rest pair every term is flat "
            "(drift ≈ 0). For a **dynamic** run: interpret the first few steps "
            "explicitly, then the general trend. Sections 4–6 below make this concrete."),
        new_markdown_cell("## 4. Energy conservation (dynamics correctness gate)"),
        new_code_cell(
            "print('max |E(t)-E(0)| (eV):', {k: round(v,4) for k,v in res.conservation.items()})\n"
            "if res.is_dynamic:\n"
            "    fig, axc = plt.subplots(figsize=(5,3), constrained_layout=True)\n"
            "    for key,lab in [('E_conserved_classical','classical: E_elec+E_proj_KE+U_proj_bg'),\n"
            "                    ('E_conserved_wp','wp: E_elec')]:\n"
            "        axc.plot(s.step, s[key]-s[key].iloc[0], marker='.', label=lab)\n"
            "    axc.set_xlabel('step'); axc.set_ylabel('E(t)-E(0) (eV)'); axc.legend(fontsize=7)\n"
            "    axc.set_title('conserved-total drift'); fig.savefig('twin_conservation.png', dpi=150)\n"
            "else:\n"
            "    print('static pair — conservation trivially flat')"),
        new_markdown_cell(
            "_Narrative — gate._ A flat conserved total (classical: `E_elec + E_proj_KE + "
            "U_proj_bg`) validates the Ehrenfest force; a drift means an integrator/force bug, "
            "not physics."),
        new_markdown_cell(
            "## 5. Dynamics — trajectory, stopping, residual collapse\n\n_Dynamic pairs only._"),
        new_code_cell(
            "if res.is_dynamic and 'proj_z' in s.columns:\n"
            "    fig, ax = plt.subplots(1, 3, figsize=(12, 3.3), constrained_layout=True)\n"
            "    ax[0].plot(s.step, s.proj_z, 'o-', ms=3); ax[0].set_title('classical projectile z(t)')\n"
            "    ax[0].set_xlabel('step'); ax[0].set_ylabel('z (Bohr)')\n"
            "    ax[1].plot(s.step, s.proj_ke_classical, 'o-', ms=3, label='proj KE (classical)')\n"
            "    ax[1].plot(s.step, s.U_proj_bg, 's-', ms=3, label='U_proj_bg')\n"
            "    ax[1].set_xlabel('step'); ax[1].set_ylabel('eV'); ax[1].legend(fontsize=7)\n"
            "    ax[1].set_title('classical stopping / proj-bg')\n"
            "    ax[2].plot(s.step, s.residual, 'o-', ms=3, color='crimson')\n"
            "    ax[2].set_xlabel('step'); ax[2].set_ylabel('residual R (eV)')\n"
            "    ax[2].set_title('WP self-Hartree residual collapse')\n"
            "    fig.savefig('twin_dynamics.png', dpi=150)\n"
            "    dke = s.proj_ke_classical.iloc[0]-s.proj_ke_classical.iloc[-1]\n"
            "    dz  = abs(s.proj_z.iloc[0]-s.proj_z.iloc[-1])\n"
            "    print(f'classical stopping: dKE={dke:.3f} eV over {dz:.2f} Bohr -> S~{dke/dz:.4f} eV/Bohr')\n"
            "    print(f'residual R: {s.residual.iloc[0]:.2f} -> {s.residual.iloc[-1]:.2f} eV')\n"
            "    print(f'quantum stopping proxy E_deposited_wp(final): {s.E_deposited_wp.iloc[-1]:.3f} eV')\n"
            "else:\n"
            "    print('static pair — no dynamics section')"),
        new_markdown_cell(
            "_Narrative — the quantum effect._ The classical projectile decelerates "
            "(proj KE → classical stopping via `stopping-power-extraction`). The residual "
            "collapses as the WP disperses (§6). Quantum stopping is `E_deposited_wp` (total "
            "electronic energy), NOT the WP kinetic — the WP orbital is not the projectile."),
        new_markdown_cell(
            "## 6. WP dispersion — the mechanism of the residual collapse\n\n"
            "σ_z(t) from the saved WP frames vs analytic free-packet dispersion "
            "σ(t)=σ₀·√(1+(t/2mσ₀²)²)."),
        new_code_cell(
            "import glob\n"
            "frames = sorted(glob.glob(WP_DIR + '/frames/wp/density_t*.vti'))\n"
            "if res.is_dynamic and frames:\n"
            "    from inqview import load_vti\n"
            "    sig0 = (res.sigma_wp or 0.5)/np.sqrt(2)   # density std at t=0\n"
            "    rows=[]\n"
            "    for f in frames:\n"
            "        d=load_vti(f); n=d.data; xx,yy,zz=d.x,d.y,d.z\n"
            "        dxx=xx[1]-xx[0]; dyy=yy[1]-yy[0]; dzz=zz[1]-zz[0]\n"
            "        nz=n.sum(axis=(0,1))*dxx*dyy; nn=nz.sum()*dzz\n"
            "        zc=(zz*nz).sum()*dzz/nn; sz=np.sqrt(((zz-zc)**2*nz).sum()*dzz/nn)\n"
            "        st=int(f.split('_t')[-1].split('.')[0]); rows.append((st,zc,sz))\n"
            "    fr=pd.DataFrame(rows,columns=['step','wp_centroid_z','wp_sigma_z'])\n"
            "    dt=(s.time_au.iloc[1]-s.time_au.iloc[0]) if len(s)>1 else 0.05\n"
            "    m=1.0; t=fr.step*dt; sig_an=sig0*np.sqrt(1+(t/(2*m*sig0**2))**2)\n"
            "    fig, ax = plt.subplots(1,2, figsize=(8,3.2), constrained_layout=True)\n"
            "    ax[0].plot(fr.step, fr.wp_sigma_z, 'o-', label='WP σ_z (measured)')\n"
            "    ax[0].plot(fr.step, sig_an, 'k--', label='analytic free dispersion')\n"
            "    ax[0].set_xlabel('step'); ax[0].set_ylabel('σ_z (Bohr)'); ax[0].legend(fontsize=7)\n"
            "    ax[0].set_title('WP dispersive spreading')\n"
            "    ax[1].plot(fr.step, fr.wp_centroid_z, 'o-', label='WP centroid')\n"
            "    if 'proj_z' in s.columns:\n"
            "        ax[1].plot(s.step, s.proj_z, '-', label='classical proj_z')\n"
            "    ax[1].set_xlabel('step'); ax[1].set_ylabel('z (Bohr)'); ax[1].legend(fontsize=7)\n"
            "    ax[1].set_title('WP centroid vs classical (tracking?)')\n"
            "    fig.savefig('wp_spreading.png', dpi=150)\n"
            "    print(fr.round(3).to_string(index=False))\n"
            "else:\n"
            "    print('no WP frames — σ(t) overlay skipped')"),
        new_markdown_cell(
            "_Narrative._ If σ_z tracks the analytic curve, the spreading is correct QM (not "
            "numerics). R ∝ 1/σ, so dispersion crushes the WP self-Hartree — the dominant "
            "classical-vs-quantum difference when the centroids track."),
        new_markdown_cell(
            "## 7. Density n(z,t) — classical vs WP, and Δn\n\n"
            "z-lineouts (∫dx dy) of each run's density over time, and the difference "
            "Δn = n_WP − n_classical. Slab faces at z=±12.5 (white lines)."),
        new_code_cell(
            "import glob\n"
            "from inqview import load_vti\n"
            "def zlineout(vti):\n"
            "    d=load_vti(vti); n=d.data; xx,yy,zz=d.x,d.y,d.z\n"
            "    return zz, n.sum(axis=(0,1))*(xx[1]-xx[0])*(yy[1]-yy[0])\n"
            "clf=sorted(glob.glob(CL_DIR+'/frames/total/density_t*.vti'))\n"
            "wpf=sorted(glob.glob(WP_DIR+'/frames/total/density_t*.vti'))\n"
            "if clf and wpf:\n"
            "    m=min(len(clf),len(wpf)); clf,wpf=clf[:m],wpf[:m]\n"
            "    steps=[int(f.split('_t')[-1].split('.')[0]) for f in wpf]\n"
            "    z,_=zlineout(wpf[0])\n"
            "    NC=np.array([zlineout(f)[1] for f in clf]); NW=np.array([zlineout(f)[1] for f in wpf])\n"
            "    fig,ax=plt.subplots(1,3,figsize=(12,3.4),constrained_layout=True)\n"
            "    ext=[z[0],z[-1],steps[-1],steps[0]]\n"
            "    ax[0].imshow(NC,aspect='auto',extent=ext,cmap='viridis'); ax[0].set_title('classical n(z,t)')\n"
            "    ax[1].imshow(NW,aspect='auto',extent=ext,cmap='viridis'); ax[1].set_title('WP n(z,t)')\n"
            "    im=ax[2].imshow(NW-NC,aspect='auto',extent=ext,cmap='RdBu_r'); ax[2].set_title('Δn = WP − classical')\n"
            "    for a in ax: a.axvline(-12.5,color='w',lw=.6); a.axvline(12.5,color='w',lw=.6); a.set_xlabel('z (Bohr)')\n"
            "    ax[0].set_ylabel('step'); fig.colorbar(im,ax=ax[2]); fig.savefig('density_carpets.png',dpi=150)\n"
            "else:\n"
            "    print('density frames missing in one/both runs')"),
        new_markdown_cell(
            "_Narrative._ The classical carpet is a rigid stripe moving in z; the WP carpet "
            "spreads (and may split/reflect). Δn localises WHERE the quantum treatment departs "
            "from classical — the reflected lobe, the spread, the polarisation of the slab."),
        new_markdown_cell(
            "## 8. WP − classical energy budget (bar plot)\n\n"
            "Every component's Δ at the final step — where the quantum energy difference lives.\n\n"
            "_(The animated density and pairwise GIFs are at the top — §0 — for visual intuition.)_"),
        new_code_cell(
            "sf=res.steps.iloc[-1]; comps={'dKin':float(sf.dKin),'dXC':float(sf.dXC)}\n"
            "if res.pairwise is not None:\n"
            "    pf=res.pairwise.iloc[-1]\n"
            "    comps.update({'ΔE_SS':float(pf.d_e_ss),'ΔE_PP':float(pf.d_e_pp),'ΔE_PS':float(pf.d_e_ps),\n"
            "                  'ΔE_SB':float(pf.d_e_sb),'ΔE_PB':float(pf.d_e_pb)})\n"
            "fig,ax=plt.subplots(figsize=(7,3.6),constrained_layout=True)\n"
            "ks=list(comps); vs=[comps[k] for k in ks]\n"
            "ax.bar(ks,vs,color=['C0' if v>=0 else 'C3' for v in vs]); ax.axhline(0,color='k',lw=.5)\n"
            "ax.set_ylabel('WP − classical (eV)'); ax.set_title(f'energy-budget difference @ step {int(sf.step)}')\n"
            "plt.setp(ax.get_xticklabels(),rotation=45,ha='right'); fig.savefig('wp_minus_cl_bars.png',dpi=150)"),
        new_markdown_cell(
            "## 9. Summary\n\n_Narrative._ Static ledger: `dKin` (localisation) + `residual` "
            "(WP self-Hartree) + `dXC` → irreducible LDA one-electron **SIE**. Dynamic: the "
            "wavepacket disperses (σ grows), its self-Hartree collapses, reflects/tunnels at the "
            "slab, and cannot be held rigid like the classical projectile — the quantum effect, "
            "now decomposed term-by-term (gauge-clean) and visualised in n(z,t), the bar budget, "
            "and the pairwise GIF."),
    ]


def build(wp_dir: str, cl_dir: str, out_path: str, title: str, execute: bool = True):
    # Absolute paths: the notebook executes with cwd = its own directory.
    wp_dir, cl_dir = str(Path(wp_dir).resolve()), str(Path(cl_dir).resolve())
    nb = new_notebook(cells=_cells(wp_dir, cl_dir, title))
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    if execute:
        from nbclient import NotebookClient
        NotebookClient(nb, timeout=600, kernel_name="python3",
                       resources={"metadata": {"path": str(Path(out_path).parent or ".")}}).execute()
    nbformat.write(nb, out_path)
    return out_path


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build the executed twin-run analysis notebook.")
    ap.add_argument("--pair"); ap.add_argument("--wp"); ap.add_argument("--classical")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="Twin-run energy decomposition")
    ap.add_argument("--no-exec", action="store_true", help="assemble without executing")
    args = ap.parse_args(argv)
    if args.pair:
        wp, cl = str(Path(args.pair) / "wp"), str(Path(args.pair) / "classical")
    elif args.wp and args.classical:
        wp, cl = args.wp, args.classical
    else:
        ap.error("give --pair DIR or --wp DIR --classical DIR")
    path = build(wp, cl, args.out, args.title, execute=not args.no_exec)
    print("wrote", path)


if __name__ == "__main__":
    main()
