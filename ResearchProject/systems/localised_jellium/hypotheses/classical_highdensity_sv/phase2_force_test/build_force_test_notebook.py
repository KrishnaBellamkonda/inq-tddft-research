#!/usr/bin/env python3
"""Build + execute the Phase-2 Test-A notebook: INQ projectile_force_z vs the
closed-form two-Gaussian Coulomb force. Run with venv python once force_test.csv exists.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent
NB = HERE / "phase2_force_test.ipynb"
RES = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/classical_highdensity_sv/force_test/results")

cells = []
cells.append(new_markdown_cell(
    "# Phase 2 · Test A — force validation vs the analytic two-Gaussian force\n\n"
    "Campaign `classical-highdensity-sv`. Independent, closed-form check of the "
    "Ehrenfest force operator `inqkit::dynamics::projectile_force_z` (a symmetric "
    "finite difference of the on-grid Poisson interaction integral).\n\n"
    "**Setup:** a fixed unit Gaussian source (σ_s) at z=0; the projectile Gaussian "
    "(σ_pot) swept along z in a **finite (periodicity 0)** box so the Poisson "
    "potential is free-space and the closed form applies.\n\n"
    "**Analytic anchor** (two Gaussians, separation d, σ12=√(σ_pot²+σ_s²), a=1/(√2·σ12)):\n"
    "- energy  E(d) = erf(a·d)/d\n"
    "- force   F(d) = −dE/dd = erf(a·d)/d² − (2a/√π)·e^(−a²d²)/d"))

cells.append(new_code_cell(
    "import numpy as np, csv, matplotlib.pyplot as plt\n"
    "from scipy.special import erf\n"
    f"RES = {RES!r}\n"
    "import re\n"
    "txt = open(RES+'/run_summary.txt').read()\n"
    "def grab(key, default):\n"
    "    m = re.search(key+r'\\s*=\\s*([-0-9.eE]+)', txt); return float(m.group(1)) if m else default\n"
    "sig_p = grab('sigma_pot', 0.35355); sig_s = grab('sigma_s', 0.5); src_norm = grab('source_norm', float('nan'))\n"
    "print('sigma_pot=%.5f  sigma_s=%.5f  source_norm=%.4f (should be 1)' % (sig_p, sig_s, src_norm))\n"
    "rows=list(csv.DictReader(open(RES+'/force_test.csv')))\n"
    "d=np.array([float(r['d']) for r in rows])\n"
    "E_num=np.array([float(r['E_num']) for r in rows])\n"
    "F_num=np.array([float(r['F_num']) for r in rows])\n"
    "s12=np.sqrt(sig_p**2+sig_s**2); a=1/(np.sqrt(2)*s12)\n"
    "E_ana=erf(a*d)/d\n"
    "F_ana=erf(a*d)/d**2 - (2*a/np.sqrt(np.pi))*np.exp(-a**2*d**2)/d\n"))

cells.append(new_markdown_cell(
    "## Energy and force: INQ (points) vs analytic (line)"))
cells.append(new_code_cell(
    "fig,ax=plt.subplots(1,2,figsize=(12,4.2))\n"
    "ax[0].plot(d,E_ana,'-',label='analytic erf(ad)/d'); ax[0].plot(d,E_num,'o',ms=4,label='INQ drag_energy')\n"
    "ax[0].set_xlabel('separation d (Bohr)'); ax[0].set_ylabel('E(d) (Ha)'); ax[0].set_title('interaction energy'); ax[0].legend()\n"
    "ax[1].plot(d,F_ana,'-',label='analytic −dE/dd'); ax[1].plot(d,F_num,'o',ms=4,label='INQ projectile_force_z')\n"
    "ax[1].set_xlabel('separation d (Bohr)'); ax[1].set_ylabel('F_z(d) (Ha/Bohr)'); ax[1].set_title('Hellmann-Feynman force'); ax[1].legend()\n"
    "fig.tight_layout(); plt.show()"))

cells.append(new_markdown_cell(
    "## Residuals + normalisation\n"
    "The ratio F_num/F_ana should be a constant ≈ 1 (a Poisson-convention factor "
    "would show as a constant ≠ 1; a *varying* ratio would mean the force LAW is "
    "wrong). We report the median ratio and the max shape deviation."))
cells.append(new_code_cell(
    "m = d > 1.5   # avoid the very-near field where the finite grid under-resolves\n"
    "kE = np.median(E_num[m]/E_ana[m]); kF = np.median(F_num[m]/F_ana[m])\n"
    "shape_err = np.max(np.abs((F_num[m]/F_ana[m])/kF - 1))\n"
    "print('median E_num/E_ana = %.4f' % kE)\n"
    "print('median F_num/F_ana = %.4f  (1.00 => atomic-unit Poisson, correct force law)' % kF)\n"
    "print('max shape deviation of F_num/F_ana about its median = %.2e' % shape_err)\n"
    "fig,ax=plt.subplots(figsize=(7,3.6))\n"
    "ax.axhline(1,ls=':',c='k'); ax.plot(d,F_num/F_ana,'o-',ms=4)\n"
    "ax.set_xlabel('separation d (Bohr)'); ax.set_ylabel('F_num / F_analytic')\n"
    "ax.set_title('force ratio (flat = correct force law)'); fig.tight_layout(); plt.show()\n"
    "verdict = 'PASS' if (abs(kF-1)<0.03 and shape_err<0.03) else 'CHECK'\n"
    "print('\\nVERDICT:', verdict)"))

cells.append(new_markdown_cell(
    "## Verdict\n\n"
    "If F_num tracks the closed-form two-Gaussian force (flat ratio ≈ 1), the "
    "Ehrenfest force operator is validated independently — the anchor for Test B "
    "(energy-conserving trajectory) and Test C (perturbation vs pseudopotential). "
    "**Manual gate: yours to accept.**"))

nb = new_notebook(cells=cells, metadata={"kernelspec": {
    "display_name": "Python 3", "language": "python", "name": "python3"}})
ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
nbf.write(nb, NB)
print("WROTE + EXECUTED", NB)
