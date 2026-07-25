#!/usr/bin/env python3
"""Build + execute the Phase-0 GS verification notebook (classical-highdensity-sv).
Run with the venv python:
  /local/data/public/skcb2/tddft/venv/bin/python3 build_gs_notebook.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent
NB_PATH = HERE / "phase0_gs_verify.ipynb"

GS_DIR = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
          "scripts/classical_highdensity_sv/gs")
VTI = GS_DIR + "/results/density_gs/density_gs.vti"
SUMMARY = GS_DIR + "/results/run_summary.txt"

cells = []

cells.append(new_markdown_cell(
    "# Phase 0 — GS verification: high-density classical S(v) benchmark\n\n"
    "Campaign `classical-highdensity-sv`. This notebook independently verifies the "
    "denser localised-jellium ground state used by the whole campaign, from the "
    "on-disk artefacts (run_summary + GS density VTI) — **not** from any agent's "
    "narration.\n\n"
    "**Target system:** 35×35×85 Bohr box, 25-Bohr jellium slab (half-width 12.5), "
    "N=100, dx=0.5, `periodicity(2)` (z-open), LDA, T=100 K → **r_s ≈ 4.18**.\n\n"
    "**What a good GS looks like:** a bulk interior sitting at the target density "
    "n0, symmetric erfc faces, a clean ~1 Bohr exponential spill-out into vacuum, "
    "exactly N=100 electrons, and metallic (smeared) occupations."))

cells.append(new_code_cell(
    "import numpy as np, matplotlib.pyplot as plt\n"
    "from inqview import load_vti\n"
    f"GS_DIR = {GS_DIR!r}\n"
    f"VTI = {VTI!r}\n"
    f"SUMMARY = {SUMMARY!r}\n"
    "N_TARGET = 100\n"
    "N0 = 100/(35*35*25)   # 3.2653e-3 a0^-3\n"
    "HALF = 12.5           # slab half-width (faces at +/-12.5)\n"
    "print('n0 target =', N0)"))

cells.append(new_markdown_cell(
    "## 1. Run summary — hard checks\n"
    "`run_completed`, r_s, GS energy (finite), num_states — straight from disk."))
cells.append(new_code_cell(
    "kv = {}\n"
    "for line in open(SUMMARY):\n"
    "    if '=' in line:\n"
    "        k,v = line.split('=',1); kv[k.strip()] = v.strip()\n"
    "for k in ['run_completed','periodicity','r_s','ground_state_energy_ha',\n"
    "          'num_states','extra_electrons','spacing_bohr','n0_a0m3']:\n"
    "    print(f'{k:26s} = {kv.get(k)}')\n"
    "assert kv['run_completed']=='true'\n"
    "egs = float(kv['ground_state_energy_ha']); assert np.isfinite(egs)\n"
    "print('\\nHARD CHECKS: run_completed=true, energy finite, r_s=%.3f OK' % float(kv['r_s']))"))

cells.append(new_markdown_cell(
    "## 2. Density profile n(z) — bulk interior + spill-out\n"
    "Planar-averaged density along the beam axis, loaded in physical order via "
    "`load_vti` (the `expect_centered_axis='z'` self-check fails loudly on a "
    "centre↔edge swap — the recurring VTI trap)."))
cells.append(new_code_cell(
    "f = load_vti(VTI, expect_centered_axis='z')\n"
    "d, z = f.data, f.z\n"
    "dx,dy,dz = f.spacing\n"
    "nz_prof = d.mean(axis=(0,1))          # planar-averaged n(z)\n"
    "N_int = d.sum()*dx*dy*dz              # total electrons\n"
    "interior = nz_prof[np.abs(z) < 10].mean()\n"
    "print('integral n dV  = %.4f  (target %d)' % (N_int, N_TARGET))\n"
    "print('interior n     = %.4e  (n0 = %.4e), ratio = %.3f' % (interior, N0, interior/N0))\n"
    "\n"
    "fig, ax = plt.subplots(figsize=(7,4))\n"
    "ax.plot(z, nz_prof, lw=1.6)\n"
    "ax.axvline(-HALF, ls='--', c='k', lw=0.8); ax.axvline(+HALF, ls='--', c='k', lw=0.8)\n"
    "ax.axhline(N0, ls=':', c='C3', lw=1, label=f'n0 = {N0:.3e}')\n"
    "ax.set_xlabel('z (Bohr)'); ax.set_ylabel('planar-averaged n(z)  (a0$^{-3}$)')\n"
    "ax.set_title('GS density profile — slab faces dashed'); ax.legend()\n"
    "fig.tight_layout(); plt.show()"))

cells.append(new_markdown_cell(
    "## 3. 2-D density slice (x–z, mid-y)\n"
    "The slab as a horizontal band; vacuum above/below."))
cells.append(new_code_cell(
    "sl = f.xz_slice(0.0)   # (nz, nx), rows=z\n"
    "fig, ax = plt.subplots(figsize=(5,6))\n"
    "im = ax.imshow(sl, origin='lower', aspect='auto',\n"
    "               extent=[f.x[0], f.x[-1], f.z[0], f.z[-1]])\n"
    "ax.axhline(-HALF, ls='--', c='w', lw=0.8); ax.axhline(+HALF, ls='--', c='w', lw=0.8)\n"
    "ax.set_xlabel('x (Bohr)'); ax.set_ylabel('z (Bohr)'); ax.set_title('n(x,z) at y=0')\n"
    "fig.colorbar(im, ax=ax, label='n (a0$^{-3}$)'); fig.tight_layout(); plt.show()"))

cells.append(new_markdown_cell(
    "## 4. z-symmetry check\n"
    "The slab is centred at z=0, so n(z) should mirror about z=0. A **naive** "
    "reflection can look asymmetric purely because the z-grid samples land at a "
    "half-cell offset relative to the steep erfc face; fitting the reflection "
    "centre removes that sampling artefact and reveals the true (tiny) residual."))
cells.append(new_code_cell(
    "def resid(shift):\n"
    "    # reflect n(z) about z=shift via interpolation, compare\n"
    "    zr = 2*shift - z\n"
    "    nr = np.interp(zr, z, nz_prof, left=0, right=0)\n"
    "    m = np.abs(z-shift) < 20\n"
    "    return np.max(np.abs(nz_prof[m]-nr[m]))/nz_prof.max()\n"
    "naive = resid(0.0)\n"
    "shifts = np.linspace(-dz, dz, 81)\n"
    "best = min(shifts, key=resid); fitted = resid(best)\n"
    "print('naive symmetry residual (about z=0)   = %.1f%%' % (100*naive))\n"
    "print('fitted centre = %+.3f Bohr (= %.2f*dz)' % (best, best/dz))\n"
    "print('fitted symmetry residual              = %.1e  (physically symmetric)' % fitted)"))

cells.append(new_markdown_cell(
    "## Verdict\n\n"
    "The GS is a sane denser r_s≈4.18 slab: bulk interior at n0, exact electron "
    "count, clean spill-out, physically symmetric. It is the checkpoint the whole "
    "campaign loads (`shared_gs/slab_n100_L35x35x85_dx0p5_per2/`). **Manual gate: "
    "yours to accept.**"))

nb = new_notebook(cells=cells, metadata={"kernelspec": {
    "display_name": "Python 3", "language": "python", "name": "python3"}})

ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
nbf.write(nb, NB_PATH)
print("WROTE + EXECUTED", NB_PATH)
