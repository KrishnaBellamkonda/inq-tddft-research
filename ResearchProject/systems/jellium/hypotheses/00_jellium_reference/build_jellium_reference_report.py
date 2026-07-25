#!/usr/bin/env python3
"""Build + execute the jellium electron-gas analytical reference notebook.

A pedagogical, run-independent reference for the homogeneous electron gas (HEG)
that the jellium bath realises: N=162 electrons in a cubic periodic box, density
set by a SINGLE knob `RS` (the Wigner-Seitz radius r_s; r_s=5.69 reproduces the
L=50 Bohr box). Every analytical quantity is derived FROM SCRATCH and cross-checked
against `inqview.analysis.lindhard_elf` and the canonical shell table in
`inqkit::jellium::shells` (shells.hpp) via in-notebook assertions.

STYLE (notebook-making skill, formula-placement rule): each formula is restated as
display LaTeX immediately above the single cell that computes that quantity, and
quantities are derived ONE AT A TIME in dependency order (k_F, then v_F, then
E_F, ...). No up-front formula dump; no batch calculation cell.

This notebook analyses NO single run-set, so it has no dispatcher auto-build tail;
regenerate it by re-running this builder:

    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_jellium_reference_report.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent

nb = new_notebook()
C = nb.cells
def md(s): C.append(new_markdown_cell(s))
def code(s): C.append(new_code_cell(s))

# ============================================================================
# 1. TITLE + THE QUESTION
# ============================================================================
md(r"""# The homogeneous electron gas behind the jellium bath — an analytical reference

**System.** The jellium runs propagate a charged projectile through a *bath* of
$N=162$ electrons in a cubic periodic box, neutralised by a uniform positive
background. That bath is a finite realisation of the **homogeneous electron gas
(HEG)** — the model whose only parameter is its density.

**The question this notebook answers.** *Before* looking at any time-dependent
result, what do we already know analytically about this electron gas? Given one
number — the density — what are its Fermi scales, its discrete single-particle
energy levels, how it screens a charge, its energy budget, and where its plasmon
lives in $(q,\omega)$? These are the quantities to *know by heart* to reason about
stopping power, plasmon excitation, and finite-size effects.

**One knob.** Everything below is driven by a single control parameter `RS`
(the Wigner-Seitz radius $r_s$). $N=162$ is held fixed because it is a *closed
shell* (a jellium magic number), and the box side $L$ is **derived** from $r_s$.
Change `RS` and every table, plot and number recomputes.

**How to read this notebook.** Each quantity is derived *one at a time*: the
formula is shown immediately above the cell that evaluates it, in dependency order.

| Where this sits | |
|---|---|
| `00_jellium_reference/` (this notebook) | analytical HEG reference — no run data |
| `06_sigma_convergence/` | numerical S(v) → point-charge Lindhard limit |
| `inqview.analysis.lindhard_elf` | the production RPA loss-function / stopping code we cross-check against |
| `inqkit::jellium::shells` (`shells.hpp`) | the production shell/degeneracy table we cross-check against |
""")

# ============================================================================
# 2. CONVENTIONS + SYMBOLS (no formula dump)
# ============================================================================
md(r"""## 1. Conventions & symbols

All quantities are in **Hartree atomic units** ($\hbar=m_e=e=4\pi\varepsilon_0=1$),
with eV given alongside ($1\,\mathrm{Ha}=27.2114\,\mathrm{eV}$). HEG relations
follow Giuliani & Vignale, *Quantum Theory of the Electron Liquid* (2005, ch. 1 &
4) and Ashcroft & Mermin, *Solid State Physics* (1976). The formulas themselves
appear at their point of use below.

| Symbol | Meaning | Units |
|---|---|---|
| $r_s$ (`RS`) | Wigner-Seitz radius — **the knob** | Bohr ($a_0$) |
| $N$ | electron count (fixed, closed shell) | — |
| $L$ | cubic box side (derived) | Bohr |
| $n$ | electron number density | $a_0^{-3}$ |
| $k_F,\,v_F$ | Fermi wavevector / velocity ($v_F=k_F$ in a.u.) | $a_0^{-1}$ / a.u. |
| $E_F,\,T_F,\,\lambda_F$ | Fermi energy / temperature / wavelength | Ha / K / $a_0$ |
| $|G|^2$ | $m_x^2+m_y^2+m_z^2$, shell label | — |
| $\omega_p$ | plasma frequency | Ha |
| $k_{TF},\,\lambda_{TF}$ | Thomas–Fermi wavevector / screening length | $a_0^{-1}$ / $a_0$ |
| $q,\omega$ | momentum / energy transfer | $a_0^{-1}$ / Ha |
| $\eta$ | small Lindhard broadening (numerics) | Ha |
| $q_m=2\pi m/L$ | discrete wavevector the box supports | $a_0^{-1}$ |
| $k_B$ | Boltzmann constant $=3.166811\times10^{-6}$ | Ha/K |
""")

code(r"""import numpy as np

# universal constants / conversions (a.u.)
HA_TO_EV = 27.211386245988      # 1 Hartree in eV
KB_HA_PER_K = 3.166811563e-6    # Boltzmann constant in Ha/K (CODATA)
PI = np.pi
print("conventions loaded (Hartree a.u.)")
""")

# ============================================================================
# 3. THE KNOB -> n -> L  (each a formula + its own cell)
# ============================================================================
md(r"""## 2. Set the knob, then derive the box

`RS = 5.69` is the density of the production jellium bath. **Change `RS` here and
re-run the whole notebook to study a different electron gas.** $N=162$ is fixed.""")
code(r"""# ----------------------------------------------------------------------------
# THE KNOB.  Everything in this notebook is a function of RS alone (N fixed).
# ----------------------------------------------------------------------------
RS = 5.69          # Wigner-Seitz radius r_s  [Bohr]  <-- change me
N  = 162           # electrons; a closed shell (jellium magic number), held fixed
print(f"RS = {RS} Bohr   (the knob)")
print(f"N  = {N}        (fixed, closed shell)")
""")

md(r"""### 2.1 Number density
One electron occupies a sphere of radius $r_s$, so
$$ n=\frac{3}{4\pi r_s^{3}} . $$""")
code(r"""n = 3.0 / (4.0 * PI * RS**3)        # electrons per Bohr^3
print(f"n = {n:.6e} 1/Bohr^3")
""")

md(r"""### 2.2 The derived box side
Holding $N$ fixed, the cubic side that reproduces this density is
$$ L=\Big(\frac{N}{n}\Big)^{1/3}=\Big(\tfrac{4}{3}\pi r_s^{3}N\Big)^{1/3} . $$""")
code(r"""L = (N / n)**(1.0/3.0)              # so that N / L^3 = n exactly
print(f"L = {L:.4f} Bohr   (= 50.0 at RS=5.69 -> the simulated box)")
""")

# ============================================================================
# 4. SOURCE FILES
# ============================================================================
md(r"""## 3. Source files (cross-checked & cited)

| File (repo-relative) | Role |
|---|---|
| `ResearchProject/systems/jellium/hypotheses/00_jellium_reference/build_jellium_reference_report.py` | this builder (generates + executes this notebook) |
| `inq-stack/python/inqview/analysis/lindhard_elf.py` | production RPA $\chi_0$ / $\varepsilon$ / loss / stopping — **cross-checked** below |
| `inq-stack/include/inqkit/jellium/shells.hpp` | production shell/degeneracy table — **cross-checked** below |
| `inq-stack/python/inqview/analysis/plasmon_spectrum.py` | rt-TDDFT plasmon *peak-locator* at the discrete $q_m$ (distinct from the loss function) |
| `ResearchProject/systems/jellium/shared/configs/sv_ladder_L50_sigma0p5.hpp` | the run config naming this $r_s=5.69$ bath |

**Sources.** Lindhard (1954); Giuliani & Vignale (2005); Ashcroft & Mermin (1976);
Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981) (correlation); Legendre's
three-square theorem for the absent $|G|^2$ shells.
""")

# ============================================================================
# 5. FERMI & DENSITY SCALES (bundle A) — ONE QUANTITY PER CELL
# ============================================================================
md(r"""## 4. Fermi & density scales (bundle A)

The foundational scales, derived one at a time. Each builds on the one before.""")

md(r"""### 4.1 Fermi wavevector $k_F$
Invert the spin-2 free-fermion filling $n=k_F^3/3\pi^2$:
$$ k_F=(3\pi^2 n)^{1/3}=\frac{(9\pi/4)^{1/3}}{r_s} . $$""")
code(r"""kF = (3.0 * PI**2 * n)**(1.0/3.0)
print(f"k_F = {kF:.5f} 1/Bohr")
""")

md(r"""### 4.2 Fermi velocity $v_F$
In atomic units $v_F=\hbar k_F/m_e=k_F$ (since $\hbar=m_e=1$):
$$ v_F=k_F . $$""")
code(r"""vF = kF
print(f"v_F = {vF:.5f} a.u.")
""")

md(r"""### 4.3 Fermi energy $E_F$
The kinetic energy of the highest occupied state,
$$ E_F=\frac{\hbar^2 k_F^2}{2m_e}=\tfrac12 k_F^2 . $$""")
code(r"""EF = 0.5 * kF**2
print(f"E_F = {EF:.6f} Ha = {EF*HA_TO_EV:.4f} eV")
""")

md(r"""### 4.4 Fermi wavelength $\lambda_F$
$$ \lambda_F=\frac{2\pi}{k_F} . $$""")
code(r"""lamF = 2.0 * PI / kF
print(f"lambda_F = {lamF:.4f} Bohr")
""")

md(r"""### 4.5 Fermi temperature $T_F$
The temperature at which $k_B T$ reaches $E_F$,
$$ T_F=\frac{E_F}{k_B} . $$""")
code(r"""TF = EF / KB_HA_PER_K
print(f"T_F = {TF:.4g} K   (room temperature is far below -> the gas is degenerate)")
""")

md(r"""### 4.6 Friedel wavevector $2k_F$
The momentum scale of Friedel oscillations / the Kohn anomaly,
$$ q_{\rm Friedel}=2k_F . $$""")
code(r"""two_kF = 2.0 * kF
print(f"2k_F = {two_kF:.5f} 1/Bohr")
""")

md(r"""### 4.7 Mean inter-electron spacing
$$ d=n^{-1/3} . $$""")
code(r"""mean_spacing = n**(-1.0/3.0)
print(f"mean spacing d = {mean_spacing:.4f} Bohr")
""")

md(r"""### 4.8 Recap (bundle A)""")
code(r"""for k, v, u in [("k_F", kF, "1/Bohr"), ("v_F", vF, "a.u."),
                ("E_F", EF, "Ha"), ("E_F", EF*HA_TO_EV, "eV"),
                ("lambda_F", lamF, "Bohr"), ("T_F", TF, "K"),
                ("2k_F", two_kF, "1/Bohr"), ("d (spacing)", mean_spacing, "Bohr"),
                ("r_s (knob)", RS, "Bohr")]:
    print(f"  {k:14s} = {v:12.5g} {u}")
""")

# ============================================================================
# 6. EIGEN-ENERGIES / SHELLS
# ============================================================================
md(r"""## 5. Single-particle eigen-energies, lumped by degeneracy

Plane waves in an $L^3$ periodic box are labelled by integer triples
$\mathbf m=(m_x,m_y,m_z)$ with energy
$$ E(\mathbf m)=\frac12\Big(\frac{2\pi}{L}\Big)^2\big(m_x^2+m_y^2+m_z^2\big)
   \equiv \frac12\Big(\frac{2\pi}{L}\Big)^2|G|^2 . $$
States sharing $|G|^2$ are degenerate → a **shell**. We enumerate every triple,
group by $|G|^2$ (degeneracy = number of triples), sort by energy, and fill $N$
electrons two-per-orbital; cumulative counts at closures are the jellium **magic
numbers**. (Legendre: $|G|^2=4^k(8m+7)$ — 7,15,23,28,… — is never a sum of three
squares, so those shells are absent.)""")
code(r"""from collections import Counter

# enumerate triples and count the degeneracy of each |G|^2 -------------------
M = 6                                          # covers all filled shells (|G|^2<=6)
deg = Counter()
for mx in range(-M, M+1):
    for my in range(-M, M+1):
        for mz in range(-M, M+1):
            deg[mx*mx + my*my + mz*mz] += 1     # one plane-wave orbital per triple

# order shells by energy (= by |G|^2) and fill N electrons -------------------
E_scale = 0.5 * (2.0*PI / L)**2                 # E = E_scale * |G|^2
shells, magic = [], []
cum_orb = cum_e = 0
print(f"{'|G|^2':>5} {'deg':>4} {'E [Ha]':>10} {'E [eV]':>9} "
      f"{'cum.orb':>8} {'cum.e-':>7} {'filled?':>8}")
for gsq in sorted(deg):
    g = deg[gsq]; E = E_scale * gsq
    cum_orb += g; cum_e += 2*g                  # 2 electrons per spatial orbital
    filled = cum_e <= N
    if filled: magic.append(cum_e)
    shells.append((gsq, g, E))
    if gsq <= 9:
        print(f"{gsq:5d} {g:4d} {E:10.5f} {E*HA_TO_EV:9.4f} "
              f"{cum_orb:8d} {cum_e:7d} {'full' if filled else '--':>8}")
print(f"\nJellium magic numbers (cumulative e- at shell closures): {magic}")
top_gsq = max(gsq for gsq, g, E in shells
              if 2*sum(deg[gg] for gg in sorted(deg) if gg <= gsq) <= N)
print(f"Highest filled |G|^2 = {top_gsq}, E_top = {E_scale*top_gsq:.5f} Ha "
      f"vs bulk E_F = {EF:.5f} Ha (finite-size gap is expected).")
""")

md(r"""### 5.1 The discrete level diagram
Each horizontal line is a degenerate shell (degeneracy annotated); filled shells
(below $E_F$) are solid, empty ones dashed.""")
code(r"""import matplotlib.pyplot as plt
from inqview.visualisation import style
style.apply_theme()

fig, ax = plt.subplots(figsize=(6.4, 4.6))
cum = 0
for gsq, g, E in shells:
    if gsq > 9: continue
    after = cum + 2*g; filled = after <= N
    ax.hlines(E*HA_TO_EV, 0, 1, color=('C0' if filled else '0.6'),
              lw=2.2 if filled else 1.0, ls=('-' if filled else '--'))
    ax.text(1.02, E*HA_TO_EV, f"$|G|^2={gsq}$, deg={g}" + ("  (filled)" if filled else ""),
            va='center', fontsize=8, color=('C0' if filled else '0.5'))
    cum = after
ax.axhline(EF*HA_TO_EV, color='C3', ls=':', lw=1.5)
ax.text(0.0, EF*HA_TO_EV, ' bulk $E_F$', color='C3', va='bottom', fontsize=8)
ax.set_xlim(0, 1.8); ax.set_xticks([])
ax.set_ylabel('single-particle energy (eV)')
ax.set_title(f'Jellium box levels  (N={N}, L={L:.1f} Bohr, $r_s$={RS})')
fig.tight_layout(); fig.savefig('fig_box_levels.png', dpi=140); plt.show()
""")

# ============================================================================
# 7. LINEAR RESPONSE / LOSS FUNCTION
# ============================================================================
md(r"""## 6. Linear response: the RPA loss function $L(q,\omega)$

Built one piece at a time: susceptibility $\chi_0$ → dielectric $\varepsilon$ →
loss $L$, then the plasmon frequency and dispersions.""")

md(r"""### 6.1 Lindhard susceptibility $\chi_0(q,\omega)$
With $z=q/2k_F$ and the **full complex** reduced frequency $u=(\omega+i\eta)/(qv_F)$,
$$ \chi_0=-\frac{k_F}{\pi^2}\Big[\tfrac12+\tfrac{1}{8z}\big((1-(z\!-\!u)^2)\ln\tfrac{z-u+1}{z-u-1}+(1-(z\!+\!u)^2)\ln\tfrac{z+u+1}{z+u-1}\big)\Big]. $$
Using only $\mathrm{Re}\,u$ would kill the undamped plasmon pole and break the
f-sum rule — the complex $u$ is essential (`lindhard_elf` docstring).""")
code(r"""def chi0_scratch(q, omega, kF, eta=1e-3):
    '''Lindhard chi_0(q, omega), full complex argument (Giuliani-Vignale sign,
    chi_0 < 0). q, omega broadcast together; returns complex array.'''
    q = np.asarray(q, float); omega = np.asarray(omega, float)
    q, omega = np.broadcast_arrays(q, omega)
    vF = kF
    z = q / (2.0*kF); z = np.where(z == 0.0, 1e-12, z)   # guard q=0
    u = (omega + 1j*eta) / (q*vF)                         # COMPLEX reduced freq
    a, b = z - u, z + u
    clog = lambda x: np.log((x + 1.0)/(x - 1.0))          # principal complex log
    F = 0.5 + (1.0/(8.0*z))*((1.0 - a*a)*clog(a) + (1.0 - b*b)*clog(b))
    return -(kF/PI**2) * F
print("chi0_scratch defined")
""")

md(r"""### 6.2 RPA dielectric function $\varepsilon_{\rm RPA}(q,\omega)$
$$ \varepsilon_{\rm RPA}=1-\frac{4\pi}{q^2}\,\chi_0(q,\omega). $$""")
code(r"""def eps_rpa_scratch(q, omega, kF, eta=1e-3):
    q = np.asarray(q, float)
    return 1.0 - (4.0*PI/np.where(q == 0.0, np.inf, q*q)) * chi0_scratch(q, omega, kF, eta)
print("eps_rpa_scratch defined")
""")

md(r"""### 6.3 Energy-loss function $L(q,\omega)$
$$ L(q,\omega)=\operatorname{Im}\!\Big[\frac{-1}{\varepsilon_{\rm RPA}(q,\omega)}\Big]. $$""")
code(r"""def loss_scratch(q, omega, kF, eta=1e-3):
    return -(1.0/eps_rpa_scratch(q, omega, kF, eta)).imag
print("loss_scratch defined")
""")

md(r"""### 6.4 Plasma frequency $\omega_p$
The $q\to0$ collective oscillation,
$$ \omega_p=\sqrt{4\pi n}. $$""")
code(r"""omega_p = np.sqrt(4.0*PI*n)
print(f"omega_p = {omega_p:.5f} Ha = {omega_p*HA_TO_EV:.3f} eV")
""")

md(r"""### 6.5 Plasmon dispersion (Bohm–Gross)
Small-$q$ correction to the plasmon energy,
$$ \omega_{\rm pl}(q)=\sqrt{\omega_p^2+\tfrac35 v_F^2 q^2}. $$""")
code(r"""def bohm_gross(q):
    return np.sqrt(omega_p**2 + 0.6*vF**2*np.asarray(q, float)**2)
print(f"omega_pl(q=0.1) = {bohm_gross(0.1):.5f} Ha")
""")

md(r"""### 6.6 Electron–hole continuum edges
The band of single particle–hole excitations,
$$ \omega_\pm(q)=\tfrac12 q^2\pm q v_F. $$""")
code(r"""def eh_upper(q): return 0.5*np.asarray(q,float)**2 + np.asarray(q,float)*vF
def eh_lower(q): return np.abs(0.5*np.asarray(q,float)**2 - np.asarray(q,float)*vF)
print("e-h continuum edge functions defined")
""")

md(r"""### 6.7 The $(q,\omega)$ loss map with discrete box modes overlaid
Colour is $\log_{10}L(q,\omega)$. The bright ridge is the plasmon; the white
curves are the e-h continuum edges $\omega_\pm$; the dashed curve is Bohm–Gross.
Vertical ticks mark the discrete $q_m=2\pi m/L$ the finite box supports — the only
$q$ your rt-TDDFT can resolve.""")
code(r"""qg = np.linspace(0.02, 2.0*kF + 0.6, 400)
wg = np.linspace(1e-3, 0.45, 400)
QQ, WW = np.meshgrid(qg, wg)
Lmap = loss_scratch(QQ, WW, kF, eta=1e-2)

fig, ax = plt.subplots(figsize=(6.8, 4.6))
pcm = ax.pcolormesh(qg, wg*HA_TO_EV, np.log10(np.clip(Lmap, 1e-6, None)),
                    cmap=style.cmap_for('sequential'), shading='auto')
fig.colorbar(pcm, ax=ax, label=r'$\log_{10} L(q,\omega)$')
ax.plot(qg, eh_upper(qg)*HA_TO_EV, color='w', lw=1.0, alpha=0.7)
ax.plot(qg, eh_lower(qg)*HA_TO_EV, color='w', lw=1.0, alpha=0.7)
ax.plot(qg, bohm_gross(qg)*HA_TO_EV, 'w--', lw=1.4, label='Bohm-Gross')
for m in range(1, 7):
    q = 2*PI*m/L
    if q <= qg[-1]:
        ax.axvline(q, color='C1', lw=0.8, alpha=0.6)
        ax.text(q, 0.44*HA_TO_EV, f"$q_{{{m}}}$", color='C1', fontsize=7, ha='center')
ax.set_xlabel(r'$q$  (1/Bohr)'); ax.set_ylabel(r'$\omega$  (eV)')
ax.set_title('RPA loss function $L(q,\\omega)$ with box modes $q_m$')
ax.legend(loc='lower right', fontsize=8)
fig.tight_layout(); fig.savefig('fig_loss_map.png', dpi=140); plt.show()
""")

md(r"""### 6.8 Landau-damping cutoff $q_c$
The plasmon stays a sharp excitation until its dispersion meets the top of the
e-h continuum; $q_c$ solves
$$ \omega_{\rm pl}(q_c)=\omega_+(q_c)=\tfrac12 q_c^2+q_c v_F. $$""")
code(r"""qfine = np.linspace(1e-3, 2.0*kF + 1.0, 5000)
qc = qfine[np.argmin(np.abs(bohm_gross(qfine) - eh_upper(qfine)))]
print(f"q_c ~ {qc:.4f} 1/Bohr  (beyond this the plasmon Landau-damps)")
""")

md(r"""### 6.9 Plasmon energy at each discrete box mode $q_m$
For $q_m=2\pi m/L$ we list the Bohm–Gross energy and the measured loss peak
$\arg\max_\omega L(q_m,\omega)$. They agree at small $q$ and diverge once the mode
passes $q_c$ (Landau damping).""")
code(r"""wfine = np.linspace(1e-3, 0.45, 4000)
print(f"{'m':>2} {'q_m[1/Bohr]':>12} {'BohmGross[Ha]':>14} {'[eV]':>8} "
      f"{'peakL[Ha]':>11} {'[eV]':>8} {'damped?':>8}")
for m in range(1, 7):
    q = 2*PI*m/L; wbg = bohm_gross(q)
    Lw = loss_scratch(np.full_like(wfine, q), wfine, kF, eta=5e-3)
    wpk = wfine[int(np.argmax(Lw))]
    print(f"{m:2d} {q:12.4f} {wbg:14.5f} {wbg*HA_TO_EV:8.3f} "
          f"{wpk:11.5f} {wpk*HA_TO_EV:8.3f} {'yes' if q>qc else 'no':>8}")
""")

# ============================================================================
# 6.10 MINIMUM PROPAGATION TIME TO RESOLVE THE SPECTRUM (RAYLEIGH)
# ============================================================================
md(r"""### 6.10 Minimum propagation time to resolve the spectrum (Rayleigh)

To separate two spectral lines by Fourier analysis of a signal of total duration
$T$, their spacing must exceed one FFT bin: the bin width is $\Delta f_{\rm
bin}=1/T$, so two frequencies $f_i,f_j$ are resolved only if $|f_i-f_j|\ge 1/T$
(**Rayleigh frequency-resolution criterion**). Hence the minimum propagation time
to resolve a whole set of frequencies $F$ is
$$ T_{\min}=\frac{1}{\min_{i\ne j}|f_i-f_j|}=\frac{1}{\min(\Delta f)} . $$

**Which frequencies?** In an rt-TDDFT run the induced density is Fourier-analysed
mode by mode at the box wavevectors $q_m=2\pi m/L$ (see `plasmon_spectrum.py`). The
features that appear in those spectra — all analytic functions already built above —
are the **plasmon line** $\omega_{\rm pl}(q_m)$ (§6.5) and the **electron-hole band
edges** $\omega_\pm(q_m)$ (§6.6). (We use the band edges, not the continuum
interior: the loss function is a *band* between $\omega_\pm$, not a set of separate
lines.)

**Angular vs ordinary frequency.** With $\hbar=1$, an energy $\omega$ in Ha *is* an
angular frequency; the ordinary frequency is $f=\omega/2\pi$. Therefore
$$ T_{\min}=\frac{1}{\min(\Delta f)}=\frac{2\pi}{\min(\Delta\omega)} . $$

We pool $\omega_\pm(q_m)$ and $\omega_{\rm pl}(q_m)$ over the resolvable modes into a
combined set and take the global smallest gap — a **conservative** target: a $T$
that resolves the closest pair *anywhere* resolves every individual mode's features
(the true per-mode requirement is $\le T_{\min}$). The pair that *produces*
$\min(\Delta f)$ is reported explicitly — the hardest-to-resolve pair, which may be
two nearby continuum onsets $\omega_-(q_m)$ (these bunch up because $\omega_-$ has a
zero at $q=2k_F$) or the plasmon meeting the upper e-h edge near $q_c$.""")

md(r"""#### Electron-hole excitation set $\{\omega_\pm(q_m)\}$
Two continuum edges per box mode, $\omega_\pm(q_m)=\tfrac12 q_m^2\pm q_m v_F$.""")
code(r"""# number of box modes to include (q_m = 2 pi m / L, m = 1..M_MODES)
M_MODES = 6

F_eh = []   # list of (label, omega_in_Ha)
for m in range(1, M_MODES+1):
    q = 2*PI*m/L
    F_eh.append((f"eh-(q{m})", float(eh_lower(q))))
    F_eh.append((f"eh+(q{m})", float(eh_upper(q))))
print("Electron-hole excitation frequencies (continuum edges omega_+-(q_m)):")
for lab, w in F_eh:
    print(f"  {lab:10s} omega = {w:.5f} Ha = {w*HA_TO_EV:7.3f} eV")
""")

md(r"""#### Plasmon set $\{\omega_{\rm pl}(q_m)\}$
One Bohm–Gross line per box mode, $\omega_{\rm pl}(q_m)=\sqrt{\omega_p^2+\tfrac35
v_F^2 q_m^2}$.""")
code(r"""F_pl = []
for m in range(1, M_MODES+1):
    q = 2*PI*m/L
    F_pl.append((f"pl(q{m})", float(bohm_gross(q))))
print("Plasmon frequencies (Bohm-Gross omega_pl(q_m)):")
for lab, w in F_pl:
    print(f"  {lab:10s} omega = {w:.5f} Ha = {w*HA_TO_EV:7.3f} eV")
""")

md(r"""#### Combined set, smallest $\Delta f$, and $T_{\min}$
Pool both sets, sort, take the smallest ordinary-frequency gap, and report the pair
that produces it.""")
code(r"""AU_TIME_FS = 0.02418884254       # 1 atomic time unit in femtoseconds

F = sorted(F_eh + F_pl, key=lambda t: t[1])     # combined set, sorted by omega
print(f"Combined set F ({len(F)} frequencies), sorted by omega:")
for lab, w in F:
    print(f"  {lab:10s} omega = {w:.5f} Ha   f = omega/2pi = {w/(2*PI):.6f} /a.u.time")

# smallest spacing in ORDINARY frequency f = omega / 2pi --------------------
gaps = [(abs(w2 - w1)/(2*PI), l1, l2, w1, w2)
        for (l1, w1), (l2, w2) in zip(F[:-1], F[1:])]
df_min, la, lb, wa, wb = min(gaps, key=lambda g: g[0])
T_min = 1.0/df_min                              # = 2 pi / min(delta omega)

print(f"\nsmallest delta f = {df_min:.6e} /a.u.time")
print(f"  PRODUCED BY the pair: {la} (omega={wa:.5f} Ha) and {lb} (omega={wb:.5f} Ha)")
print(f"  energy gap = {abs(wb-wa):.5f} Ha = {abs(wb-wa)*HA_TO_EV:.4f} eV")
print(f"\nT_min = 1/min(delta f) = {T_min:.1f} a.u. of time = {T_min*AU_TIME_FS:.3f} fs")
print(f"  (check: 2*pi/min(delta omega) = {2*PI/abs(wb-wa):.1f} a.u.)")
print(f"\nIllustrative: for a timestep dt=0.05 a.u., N_steps = T_min/dt = {T_min/0.05:.0f}")
print("  (swap 0.05 for your propagation dt to get the step count)")
""")

md(r"""#### The combined frequency comb
e-h band edges (grey) and plasmon lines (orange) on the energy axis; the closest
pair — which sets $T_{\min}$ — is dotted in red.""")
code(r"""fig, ax = plt.subplots(figsize=(6.8, 2.8))
for lab, w in F_eh: ax.axvline(w*HA_TO_EV, color='0.55', lw=1.0)
for lab, w in F_pl: ax.axvline(w*HA_TO_EV, color='C1', lw=1.7)
for w in (wa, wb):  ax.axvline(w*HA_TO_EV, color='C3', lw=1.0, ls=':')
ax.plot([], [], color='0.55', lw=1.0, label=r'e-h band edges $\omega_\pm(q_m)$')
ax.plot([], [], color='C1', lw=1.7, label=r'plasmon $\omega_{pl}(q_m)$')
ax.plot([], [], color='C3', lw=1.0, ls=':',
        label=f'closest pair ($T_{{min}}$={T_min*AU_TIME_FS:.2f} fs)')
ax.set_yticks([]); ax.set_xlabel(r'$\omega$ (eV)')
ax.set_title('Combined frequency set — Rayleigh resolution target')
ax.legend(loc='upper right', fontsize=8)
fig.tight_layout(); fig.savefig('fig_frequency_comb.png', dpi=140); plt.show()
""")

# ============================================================================
# 8. SCREENING (bundle B)
# ============================================================================
md(r"""## 7. Screening (bundle B)

How the gas screens the projectile's charge, derived one quantity at a time.""")

md(r"""### 7.1 Thomas–Fermi wavevector $k_{TF}$
$$ k_{TF}=\sqrt{4k_F/\pi}. $$""")
code(r"""k_TF = np.sqrt(4.0*kF/PI)
print(f"k_TF = {k_TF:.5f} 1/Bohr")
""")

md(r"""### 7.2 Screening length $\lambda_{TF}$
$$ \lambda_{TF}=1/k_{TF}. $$""")
code(r"""lam_TF = 1.0/k_TF
print(f"lambda_TF = {lam_TF:.5f} Bohr  (charge screened over ~this distance)")
""")

md(r"""### 7.3 Static dielectric $\varepsilon(q,0)$ — Thomas–Fermi vs full RPA
The TF model is the long-wavelength limit of the static RPA:
$$ \varepsilon^{TF}(q,0)=1+\frac{k_{TF}^2}{q^2},\qquad
   \varepsilon^{RPA}(q,0)=1-\frac{4\pi}{q^2}\chi_0(q,0). $$""")
code(r"""qs = np.linspace(0.05, 2.0*kF, 200)
eps_TF = 1.0 + k_TF**2/qs**2
eps_RPA0 = eps_rpa_scratch(qs, np.full_like(qs, 1e-5), kF, eta=1e-4).real

fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(qs, eps_TF, label=r'Thomas-Fermi $1+k_{TF}^2/q^2$')
ax.plot(qs, eps_RPA0, '--', label=r'full RPA $\varepsilon(q,0)$')
ax.axvline(2*kF, color='0.6', ls=':'); ax.text(2*kF, ax.get_ylim()[1]*0.9, r' $2k_F$',
           color='0.5', fontsize=8)
ax.set_xlabel(r'$q$ (1/Bohr)'); ax.set_ylabel(r'static $\varepsilon(q,0)$')
ax.set_title('Static screening: Thomas-Fermi vs RPA'); ax.legend()
fig.tight_layout(); fig.savefig('fig_screening.png', dpi=140); plt.show()
""")

# ============================================================================
# 9. ENERGY PER ELECTRON (bundle C)
# ============================================================================
md(r"""## 8. Energy per electron (bundle C)

The LDA energy budget, term by term (Ha per electron).""")

md(r"""### 8.1 Kinetic energy
Average over the filled Fermi sphere,
$$ t=\tfrac35 E_F. $$""")
code(r"""t_kin = 0.6*EF
print(f"t   = {t_kin:+.5f} Ha = {t_kin*HA_TO_EV:+.3f} eV")
""")

md(r"""### 8.2 Exchange energy
Exact HEG exchange,
$$ \epsilon_x=-\frac{3}{4\pi}k_F=-\frac{0.458166}{r_s}. $$""")
code(r"""e_x = -3.0/(4.0*PI)*kF
print(f"e_x = {e_x:+.5f} Ha = {e_x*HA_TO_EV:+.3f} eV   (check e_x*r_s = {e_x*RS:.6f}, want -0.458166)")
""")

md(r"""### 8.3 Correlation energy (Perdew–Zunger 1981)
QMC parametrisation; for $r_s\ge1$ (ours is 5.69),
$$ \epsilon_c=\frac{\gamma}{1+\beta_1\sqrt{r_s}+\beta_2 r_s},\quad
   \gamma=-0.1423,\ \beta_1=1.0529,\ \beta_2=0.3334. $$""")
code(r"""def e_correlation_PZ81(rs):
    '''Perdew-Zunger 1981 correlation energy per electron [Ha].'''
    if rs >= 1.0:
        g, b1, b2 = -0.1423, 1.0529, 0.3334
        return g / (1.0 + b1*np.sqrt(rs) + b2*rs)
    A, B, Cc, D = 0.0311, -0.048, 0.0020, -0.0116
    return A*np.log(rs) + B + Cc*rs*np.log(rs) + D*rs
e_c = e_correlation_PZ81(RS)
print(f"e_c = {e_c:+.5f} Ha = {e_c*HA_TO_EV:+.3f} eV")
""")

md(r"""### 8.4 Total energy per electron
$$ \epsilon_{\rm tot}=t+\epsilon_x+\epsilon_c. $$""")
code(r"""e_tot = t_kin + e_x + e_c
print(f"e_tot = {e_tot:+.5f} Ha = {e_tot*HA_TO_EV:+.3f} eV  "
      f"(exchange dominates correlation at this density)")
""")

# ============================================================================
# 10. STOPPING & DYNAMICS (bundle D)
# ============================================================================
md(r"""## 9. Stopping & dynamics (bundle D)

Kinematic scales for projectile energy loss. The e-h continuum (§6.6) and the
Landau cutoff $q_c$ (§6.8) are reused here.""")

md(r"""### 9.1 Velocity & time scales
The Bohr velocity $v_0=1$ a.u. is the reference; the regime is set by $v/v_F$.
The plasmon period is $T_p=2\pi/\omega_p$.""")
code(r"""print(f"Bohr velocity v0       = 1.0 a.u. (reference)")
print(f"Fermi velocity v_F     = {vF:.5f} a.u.")
print(f"plasmon energy hw_p    = {omega_p:.5f} Ha = {omega_p*HA_TO_EV:.3f} eV")
print(f"Landau cutoff q_c       = {qc:.4f} 1/Bohr   (sec 6.8)")
print(f"plasmon period 2pi/w_p = {2*PI/omega_p:.3f} a.u. of time")
""")

md(r"""### 9.2 f-sum rule
A conservation law tying the loss function's total weight to $\omega_p$ alone:
$$ \int_0^\infty \omega\,L(q,\omega)\,d\omega=\frac{\pi}{2}\,\omega_p^2 \quad (\forall q). $$""")
code(r"""q_test = 0.3
wgr = np.linspace(1e-4, 3.0, 60000)
Lwr = loss_scratch(np.full_like(wgr, q_test), wgr, kF, eta=1e-2)
fsum = np.trapezoid(wgr*Lwr, wgr)
print(f"f-sum at q={q_test}: integral = {fsum:.5f},  (pi/2) w_p^2 = {0.5*PI*omega_p**2:.5f},"
      f"  ratio = {fsum/(0.5*PI*omega_p**2):.4f}")
""")

# ============================================================================
# 11. VERIFICATION
# ============================================================================
md(r"""## 10. Verification — cross-checks against the production code

Every from-scratch result is asserted against the canonical implementations. A
green build (`nbconvert`, 0 errors) *is* the proof.

1. **Lindhard $\chi_0$, $\varepsilon$, loss** vs `inqview.analysis.lindhard_elf`.
2. **Fermi / plasmon / screening helpers** ($k_F$, $\omega_p$, $k_{TF}$) vs the same.
3. **Shell magic numbers** vs the `shells.hpp` ladder $[2,14,38,54,66,114,162]$.
4. **f-sum rule** holds to $<10^{-2}$ across several $q$.""")
code(r"""from inqview.analysis import lindhard_elf as LE

# (1) chi0 / eps / loss agree pointwise with the package
qC = np.array([0.1, 0.3, 0.5, 0.7, 1.0]); wC = np.array([0.05, 0.12, 0.18, 0.25, 0.30])
QC, WC = np.meshgrid(qC, wC)
assert np.allclose(chi0_scratch(QC, WC, kF, 1e-3), LE.chi0(QC, WC, kF, eta=1e-3), rtol=1e-10, atol=1e-12)
assert np.allclose(eps_rpa_scratch(QC, WC, kF, 1e-3), LE.epsilon_rpa(QC, WC, kF, eta=1e-3), rtol=1e-10, atol=1e-12)
assert np.allclose(loss_scratch(QC, WC, kF, 1e-3), LE.loss_function(QC, WC, kF, eta=1e-3), rtol=1e-10, atol=1e-12)
print("[1] chi0 / eps / loss match inqview.analysis.lindhard_elf  ... OK")

# (2) scalar helpers
assert np.isclose(kF, LE.kF_from_rs(RS)) and np.isclose(n, LE.density_from_kF(kF))
assert np.isclose(omega_p, LE.omega_p(kF)) and np.isclose(k_TF, LE.k_TF(kF))
print("[2] kF, n, omega_p, k_TF match the package helpers          ... OK")

# (3) magic-number ladder + degeneracies vs shells.hpp
assert magic == [2, 14, 38, 54, 66, 114, 162], magic
for gsq, g_exp in {0:1, 1:6, 2:12, 3:8, 4:6, 5:24, 6:24}.items():
    assert deg[gsq] == g_exp, (gsq, deg[gsq])
for gsq in (7, 15, 23, 28):
    assert deg[gsq] == 0, gsq                    # Legendre exclusions
print("[3] shell degeneracies + magic numbers match shells.hpp     ... OK")
print(f"    magic numbers: {magic}")

# (4) f-sum rule across several q
for q in (0.2, 0.4, 0.6, 0.8):
    wgr = np.linspace(1e-4, 3.0, 60000)
    Lwr = loss_scratch(np.full_like(wgr, q), wgr, kF, eta=1e-2)
    assert abs(np.trapezoid(wgr*Lwr, wgr)/(0.5*PI*omega_p**2) - 1.0) < 1e-2, q
print("[4] f-sum rule holds to <1e-2 across q in {0.2..0.8}        ... OK")
print("\nALL CROSS-CHECKS PASSED.")
""")

# ============================================================================
# 12. TAKEAWAY
# ============================================================================
md(r"""## 11. Takeaway

For the production bath ($r_s=5.69$, $N=162$, $L=50\,a_0$):

- **Fermi scales:** $k_F\approx0.337\,a_0^{-1}$, $E_F\approx0.057\,\mathrm{Ha}
  =1.55\,\mathrm{eV}$, $v_F=k_F$ — a *low-density* gas (lower than Na, $r_s\approx4$),
  so collective effects sit at modest energies.
- **Closed shell:** 162 closes the $|G|^2=6$ shell exactly — why the ground state
  converges cleanly. The degeneracy ladder depends on $N$ only; changing `RS`
  rescales energies but never reopens the shell.
- **Plasmon:** $\hbar\omega_p\approx0.128\,\mathrm{Ha}=3.47\,\mathrm{eV}$, where
  rt-TDDFT spectra should peak as $q\to0$; it Landau-damps beyond $q_c$ (§6.8). The
  box can only excite the discrete $q_m=2\pi m/L$.
- **Screening:** $\lambda_{TF}\approx1.5\,a_0$ — much shorter than the box.
- **Energetics:** exchange dominates correlation at this density.
- **Trust:** every formula reproduces `inqview.analysis.lindhard_elf` and
  `shells.hpp` to machine / $10^{-2}$ precision (§10), and the f-sum rule holds.

*Change `RS` in §2 and re-run to get the full analytical picture of any other
electron-gas density.*
""")

# ============================================================================
# EXECUTE + WRITE
# ============================================================================
if __name__ == "__main__":
    ep = ExecutePreprocessor(timeout=1800, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
    out = HERE / "jellium_reference.ipynb"
    nbf.write(nb, out)
    print(f"wrote {out}")
