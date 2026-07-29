#!/usr/bin/env python3
"""WP spreading tradeoff plot for choosing (sigma_WP, E, spread-threshold).

Decision aid for the quantum-vs-classical stopping runs: in the wavepacket
convention (sigma_WP, density std = sigma_WP/sqrt(2)), a free Gaussian disperses as
    sigma_WP(t) = sigma_WP * sqrt(1 + (t/sigma_WP^2)^2)        (a.u., m_e=1)
so the fractional spread over the transit (length X, time t* = X/sqrt(2E)) is
    f(sigma_WP, E) = sqrt(1 + (t*/sigma_WP^2)^2) - 1.
Inverting at fixed f gives the minimum energy for a width:
    E_min(sigma_WP, f) = X^2 / (2 sigma_WP^4 ((1+f)^2 - 1)).

Two panels:
  LEFT  spread% vs sigma_WP at E = 100/300/600 eV, with candidate % thresholds.
  RIGHT E_min vs sigma_WP for several spread thresholds, target band 100-600 eV.

Saved as fig_spread_tradeoff.png for the cap_baselines study notebook.
Transit X is the current geometry's launch->far-edge (z=-13 -> +15 = 28 Bohr);
change X if the CAP geometry changes for the new (larger-sigma) runs.

    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 make_spread_tradeoff.py
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from inqview.visualisation import style

style.apply_theme()
HERE = Path(__file__).resolve().parent
HA = 27.2114
X = 28.0  # transit length (Bohr)


def spread_frac(sw, E_eV):
    E = E_eV / HA
    tstar = X / np.sqrt(2 * E)
    return np.sqrt(1 + (tstar / sw**2) ** 2) - 1.0


def E_min_eV(sw, f):
    return X**2 / (2 * sw**4 * ((1 + f) ** 2 - 1)) * HA


sw = np.linspace(0.4, 3.5, 400)
fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.6))

# ---- LEFT: spread% vs sigma_WP at fixed energies ------------------------------
for E_eV, c in [(100, 'C3'), (300, 'C1'), (600, 'C0')]:
    axL.plot(sw, 100 * spread_frac(sw, E_eV), color=c, lw=2, label=f'E = {E_eV} eV')
for thr, lab in [(20, '20%'), (50, '50%'), (100, '100%'), (200, '200%')]:
    axL.axhline(thr, color='0.6', ls=':', lw=1)
    axL.text(3.46, thr, lab, fontsize=7, color='0.4', va='center', ha='left')
axL.axvline(0.5, color='0.3', ls='--', lw=1)
axL.annotate('sigma_WP=0.5:\n1590-4030% spread\n(packet balloons 16-40x)',
             xy=(0.5, 1590), xytext=(0.85, 800), fontsize=8, color='0.25',
             arrowprops=dict(arrowstyle='->', color='0.5'))
axL.set_yscale('log')
axL.set_xlabel(r'$\sigma_{\rm WP}$ (Bohr)')
axL.set_ylabel('spread over transit (%)')
axL.set_title('How much the WP spreads (E fixed)')
axL.legend(loc='upper right', fontsize=8)
axL.set_xlim(0.4, 3.5)

# ---- RIGHT: E_min vs sigma_WP for spread thresholds, target band --------------
for f, lab, c in [(0.20, '20%', 'C0'), (0.50, '50%', 'C2'),
                  (1.00, '100%', 'C1'), (2.00, '200%', 'C4'), (4.00, '400%', 'C3')]:
    axR.plot(sw, E_min_eV(sw, f), color=c, lw=2, label=f'spread = {lab}')
axR.axhspan(100, 600, color='0.5', alpha=0.18, zorder=0)
axR.text(3.45, 245, 'target\n100-600 eV', fontsize=8, color='0.35', ha='right')
axR.axvline(0.5, color='0.3', ls='--', lw=1)
axR.set_yscale('log')
axR.set_ylim(50, 5e4)
axR.set_xlim(0.4, 3.5)
axR.set_xlabel(r'$\sigma_{\rm WP}$ (Bohr)')
axR.set_ylabel('energy needed (eV)')
axR.set_title('Energy to hit a spread threshold')
axR.legend(loc='upper right', fontsize=8)

fig.suptitle('Wavepacket spreading tradeoff — pick (sigma_WP, E, spread%) '
             f'(transit X={X:.0f} Bohr)', fontsize=11)
fig.tight_layout()
out = HERE / 'fig_spread_tradeoff.png'
fig.savefig(out, dpi=150)
print(f'wrote {out}')

# ---- a small decision table printed alongside ---------------------------------
print('\nspread% over transit:')
print(f"{'sigma_WP':>8} {'dens std':>9} | E=100  E=300  E=600 eV")
for s in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
    print(f'{s:>8.2f} {s/np.sqrt(2):>9.2f} | '
          f'{100*spread_frac(s,100):>5.0f}  {100*spread_frac(s,300):>5.0f}  '
          f'{100*spread_frac(s,600):>5.0f}')
