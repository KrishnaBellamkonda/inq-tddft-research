#!/usr/bin/env python3
"""
plot_convergence.py — Visualise the two jellium convergence tests.

Reads results/convergence_results.csv (written by run_convergence.cpp).
Run from the 02_ground_state_convergence/ directory:

    python3 plot_convergence.py

Produces two figures:

  Figure 1 — Grid-spacing (E_cut) convergence
      E_total and T_s vs E_cut for N=40, L=40 bohr.
      Both axes in Ha. Convergence reference = finest E_cut run.
      Secondary panel shows ΔE_total = E(h) − E(h_min).

  Figure 2 — Shell-closure convergence
      T_s/N vs N for closed-shell magic numbers at fixed r_s.
      Overlaid dashed line: Thomas-Fermi bulk limit (3/5)E_F.
      Oscillations show the finite-size shell effect.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({
    'font.family':   'serif',
    'font.size':     11,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'legend.framealpha': 0.95,
})

# ── Parse output from run_convergence ────────────────────────────────────────

def parse_convergence(lines):
    """Parse the CSV output written by run_convergence.cpp."""
    data_A, data_B = [], []
    T_TF_per_N = None

    for line in lines:
        line = line.strip()
        if not line or not line.startswith('#'):
            continue

        line = line.lstrip('#').strip()

        if line.startswith('TEST_A '):
            try:
                vals = line[len('TEST_A '):].split(',')
                data_A.append([float(v) for v in vals])
            except ValueError:
                pass  # skip column-header line

        elif line.startswith('TEST_B '):
            try:
                vals = line[len('TEST_B '):].split(',')
                row = [float(vals[i]) if i != 0 else int(vals[i]) for i in range(len(vals))]
                data_B.append(row)
            except ValueError:
                pass  # skip column-header line

        elif 'T_TF_per_N' in line and '=' in line:
            T_TF_per_N = float(line.split('=')[1].split()[0])

    return (np.array(data_A) if data_A else None,
            np.array(data_B) if data_B else None,
            T_TF_per_N)


# ── Load data ─────────────────────────────────────────────────────────────────

CSV_PATH = 'results/convergence_results.csv'

if not __import__('os').path.exists(CSV_PATH):
    print(f'ERROR: {CSV_PATH} not found. Run inq-run in this directory first.')
    raise SystemExit(1)

with open(CSV_PATH) as f:
    lines = f.readlines()

A, B, T_TF_per_N = parse_convergence(lines)

if A is None and B is None:
    print("No data found. Run:  ./run_convergence > convergence_results.csv")
    sys.exit(1)

# ── Figure 1: E_cut convergence ────────────────────────────────────────────

if A is not None:
    # Columns: spacing, E_cut, E_total, T_s, E_xc, n_iter
    h        = A[:, 0]
    E_cut    = A[:, 1]
    E_total  = A[:, 2]
    T_s      = A[:, 3]
    E_xc     = A[:, 4]

    # Reference: finest grid
    E_ref = E_total[np.argmax(E_cut)]

    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                                     gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.06})

    # Top: E_total and T_s
    ax1.plot(E_cut, E_total, 'o-', color='#2166ac', lw=1.8, ms=6,
             label=r'$E_{\rm total}$ (INQ)')
    ax1.plot(E_cut, T_s,     's--', color='#d73027', lw=1.4, ms=5,
             label=r'$T_s$ (kinetic)')

    ax1.set_ylabel('Energy  (Ha)', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(ls=':', alpha=0.4)
    ax1.set_title(r'Grid-spacing convergence  ($N=40$, $L=40\,a_0$, $r_s=7.26\,a_0$)',
                  fontsize=12, pad=8)

    # Bottom: ΔE_total
    dE = np.abs(E_total - E_ref) * 1000   # mHa
    ax2.semilogy(E_cut, np.where(dE > 0, dE, 1e-10), 'o-', color='#2166ac',
                 lw=1.8, ms=6)
    ax2.axhline(10,   color='#e08000', ls='--', lw=1.2, label='10 mHa')
    ax2.axhline(1,    color='#777777', ls=':', lw=1.2,  label='1 mHa')
    ax2.set_xlabel(r'$E_{\rm cut} = \pi^2 / (2h^2)$  (Ha)', fontsize=12)
    ax2.set_ylabel(r'$|\Delta E_{\rm total}|$  (mHa)', fontsize=12)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(ls=':', alpha=0.4)

    # Annotate spacing on top axis
    ax_top = ax1.twiny()
    ax_top.set_xlim(ax1.get_xlim())
    ticks = E_cut
    ax_top.set_xticks(ticks)
    ax_top.set_xticklabels([f'{hi:.2f}' for hi in h], fontsize=9)
    ax_top.set_xlabel(r'Grid spacing $h$  (bohr)', fontsize=11, labelpad=8)

    plt.tight_layout()
    fig1.savefig('convergence_Ecut.png', bbox_inches='tight', dpi=150)
    print('Saved convergence_Ecut.png')

# ── Figure 2: Shell-closure convergence ───────────────────────────────────────

if B is not None and T_TF_per_N is not None:
    # Columns: N, L, k0, Ts, Ts_per_N, T_TF_per_N, n_iter
    N_vals     = B[:, 0].astype(int)
    Ts_per_N   = B[:, 4]
    T_TF       = T_TF_per_N   # scalar from header comment

    # Analytical shell-sum values (computed in jellium_utils, reproduced here)
    # These are the T=0 discrete shell sums at each closed N.
    Ts_per_N_analytical = Ts_per_N   # INQ values ARE the reference

    fig2, ax = plt.subplots(figsize=(7, 4.5))

    ax.plot(N_vals, Ts_per_N * 1000, 'o-', color='#2166ac',
            lw=1.8, ms=8, zorder=3, label=r'$T_s/N$  (INQ, closed shells)')

    ax.axhline(T_TF * 1000, color='#d73027', ls='--', lw=1.6,
               label=r'$(3/5)\,E_F$  (Thomas-Fermi bulk limit)')

    # Add vertical dashed lines at each N and label them
    for N, Ts_N in zip(N_vals, Ts_per_N):
        ax.annotate(f'$N={N}$',
                    xy=(N, Ts_N * 1000),
                    xytext=(N + 1, Ts_N * 1000 + 0.2),
                    fontsize=9, color='#2166ac',
                    arrowprops=dict(arrowstyle='->', color='#aaaaaa', lw=0.8))

    ax.set_xlabel(r'Electron count $N$  (closed-shell magic numbers)', fontsize=12)
    ax.set_ylabel(r'$T_s/N$  (mHa per electron)', fontsize=12)
    ax.set_title(
        r'Shell-closure finite-size effect  ($r_s = 7.26\,a_0$ fixed)',
        fontsize=12, pad=8
    )
    ax.legend(fontsize=10)
    ax.grid(ls=':', alpha=0.4)
    ax.set_xlim(-5, N_vals[-1] + 10)

    # Inset text explaining physics
    ax.text(0.62, 0.15,
            'Oscillations arise from discrete\n'
            r'$\Gamma$-point shell structure.' '\n'
            r'Amplitude $\propto 1/N$ as $N\to\infty$.',
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle='round', fc='white', alpha=0.8, ec='#cccccc'))

    plt.tight_layout()
    fig2.savefig('convergence_shells.png', bbox_inches='tight', dpi=150)
    print('Saved convergence_shells.png')

plt.show()
