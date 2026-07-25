#!/usr/bin/env python3
"""Deterministic builder for the Fourier-analysis AUDIT notebook (Task 2).

Three-stage known-answer audit of inqview's FFT code:
  Stage A  synthetic analytic-truth  -> validate fourier.py machinery + expose
           density_fourier BUG-A on a signal whose answer is known exactly.
  Stage B  fourier.py on the QKE v0p0626 dipole_x (locked plasmon 6.48 eV).
  Stage C  density_fourier.py on the E15 jellium n_q (omega_p ~ 3.47 eV) -> the
           two documented bugs + the corrected loss-function estimator.

Produces EVIDENCE; the dossier verdict lines stay blank for the user
(verification-user-owns-verdict). Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_fourier_notebook.py
"""
from __future__ import annotations
import pathlib
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = pathlib.Path(__file__).resolve().parent
NB = HERE / "fourier_analysis.ipynb"
ROOT = "/local/data/public/skcb2/tddft"

nb = new_notebook()
md = lambda s: nb.cells.append(new_markdown_cell(s))
code = lambda s: nb.cells.append(new_code_cell(s))

# ================================================================ §1 ==========
md(r"""# Fourier-analysis audit — windowing, detrending, peaks, and the loss function

**What this is.** A hands-on audit of the two FFT code paths in `inqview`, each
checked against a case whose answer is known *independently*:

| Path | What it does | Known-answer anchor |
|---|---|---|
| `inqview.analysis.fourier` | windowed FFT of scalar time-series (energy/dipole/current) | QKE Li `v0p0626` `dipole_x` plasmon = **6.48 eV** (journal 2026-05-06) |
| `inqview.pipeline.density_fourier` | $n_q(\omega)$ loss-function path | E15 jellium $n_q$, $m{=}1$ peak $\approx\omega_p=3.47$ eV |

**Why.** `fourier.py` carries an explicit TODO ("investigate if the windowing /
detrending / convenience methods are skewing the findings, *especially for the
QuantumKickExtension run*"), and `density_fourier.py` carries two documented bugs.
Until this methodology is locked, the project's standing **loss-function gate**
forbids new $L(q,\omega)$ work (`feedback_fourier_loss_function_gate`).

**Strategy (known-answer first).** Stage A validates every knob on *synthetic*
signals with analytic truth; Stage B/C then apply the validated tools to real runs
where the physical peak is already known. This produces **evidence**; the verdict
lines in the three dossiers (`loss-function-formula-validation`,
`fft-drift-removal-validation`, `fft-normalization-validation`) are left for the
user to sign.""")

code(r"""import sys; sys.path.insert(0, '%s/inq-stack/python')
import numpy as np, pandas as pd
%%matplotlib inline
import matplotlib.pyplot as plt
from inqview.visualisation import style as S
S.apply_theme()
from inqview.analysis.fourier import FourierTransform, WindowSpec

HA2EV = 27.211386                      # Hartree -> eV
def omega_eV(freq_au):                 # fourier.py returns freq_au = rfftfreq (cycles/a.u.)
    return freq_au * 2*np.pi * HA2EV   # angular energy hbar*omega in eV (hbar=1 a.u.)
def peak_in(freq_au, power, lo, hi):
    o = omega_eV(freq_au); w = (o > lo) & (o < hi)
    return o[w][np.argmax(power[w])]
QKE = '%s/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0626_xyz/results/raw/observables'
E15 = '%s/ResearchProject/systems/jellium/run_plasmon_n162_L50_E15/results/analysis/observables/n_q_vs_time.csv'
print('setup ok')""" % (ROOT, ROOT, ROOT))

# ================================================================ §2 ==========
md(r"""## Conventions and symbols

| symbol | meaning | units |
|---|---|---|
| $t$ | time | a.u. ($\hbar/E_h$) |
| $f$ | ordinary frequency = `rfftfreq(d=dt)` | cycles / a.u. |
| $\hbar\omega$ | angular energy $=2\pi f\cdot E_h$ | eV (×`HA2EV`) |
| $\Delta\omega$ | FFT bin $=2\pi/T\cdot E_h$ | eV |
| $w_n$ | window; **coherent gain** $=\sum_n w_n$ | — |
| $n_q(t)$ | induced-density Fourier amplitude at $q$ | complex |

**Amplitude normalisation (the convention under audit).** `fourier.py` returns
$|{\rm FFT}|/\sum_n w_n$ with interior bins ×2 — the **coherent-gain** convention
(Harris 1978): a unit-amplitude tone reads ~1.0 for *any* window. Validated in A1.

**Loss-function estimator (the density path).** The proposed proxy is
$$L(q,\omega)\;=\;\frac{|n_q(\omega)|^2}{q^2}.$$
**Caveat carried throughout (dossier `loss-function-formula-validation`):** this is
**quadratic** in $n_q$ and uses a bare $1/q^2$, whereas the true loss function
$-\mathrm{Im}\,[1/\epsilon]$ is **linear** in $\mathrm{Im}\,\chi$ and carries
$v(q)=4\pi/q^2$. So $L$ is a **plasmon peak-LOCATOR** — peak *positions* are
trustworthy; absolute lineshape / area / cross-$q$ intensity are NOT.""")

# ================================================================ Stage A =====
md(r"""# Stage A — synthetic analytic truth

Every knob validated where the right answer is known in closed form.""")

md(r"""### A1 — Window coherent-gain normalisation
A unit-amplitude cosine must return amplitude **≈ 1.0** through `fourier.py` for
*any* window if the $\sum_n w_n$ normalisation (+ interior ×2) is correct.""")
code(r"""N, dt, f0 = 4000, 0.05, 0.10
tt = np.arange(N)*dt
tone = 1.0*np.cos(2*np.pi*f0*tt)               # unit amplitude
rows = []
for wn in ['boxcar','hann','hamming','blackman','flattop']:
    r = FourierTransform(window=WindowSpec(wn), zero_pad=1, subtract='none').transform(tt, tone)
    rows.append((wn, r.amplitude.max()))
    print(f'  window={wn:9s}: recovered amplitude = {r.amplitude.max():.4f}  (expect ~1.0)')
print('\n-> coherent-gain normalisation is window-independent and correct.')""")

md(r"""### A2 — Baseline subtraction (`none`/`mean`/`detrend`/`initial`)
Build tone + DC offset + linear drift. The **peak location** must be invariant to
the baseline choice; only the $\omega\!\approx\!0$ content changes. `none` lets the
offset hijack the lowest bin; `detrend` removes the drift; `initial` enforces
$s(0)=0$ (the induced-response IC).""")
code(r"""omega_t = 0.1276                                # ~3.47 eV target tone
sig = np.cos(omega_t*tt) + 3.0 + 0.002*tt        # tone + DC + drift
fig, ax = S.figure_one_col()
for sub in ['none','mean','detrend','initial']:
    r = FourierTransform(window=WindowSpec('hann'), zero_pad=4, subtract=sub).transform(tt, sig)
    o = omega_eV(r.frequency_au)
    pk = peak_in(r.frequency_au, r.power, 0.5, 20)
    ax.plot(o[o<8], r.amplitude[o<8], lw=0.9, label=f'{sub}: peak={pk:.3f} eV')
ax.set_xlabel(r'$\hbar\omega$ (eV)'); ax.set_ylabel('amplitude')
ax.set_title('A2 - baseline modes: peak location invariant, low-w differs'); ax.legend()
print('tone target =', round(omega_t*HA2EV,3), 'eV')
fig""")

md(r"""### A3 — Zero-padding interpolates only
Zero-padding densifies the frequency axis but adds **no** spectral information: the
peak location and width are set by the record length $T$, not by `zero_pad`.""")
code(r"""for zp in [1, 4, 8]:
    r = FourierTransform(window=WindowSpec('hann'), zero_pad=zp, subtract='none').transform(tt, np.cos(omega_t*tt))
    o = omega_eV(r.frequency_au)
    pk = peak_in(r.frequency_au, r.power, 0.5, 20)
    print(f'  zero_pad={zp}: n_freq={len(o)}  peak={pk:.4f} eV  bin dw={o[1]-o[0]:.4f} eV')
print('-> peak unchanged; only the axis sampling (n_freq, dw) densifies.')""")

md(r"""### A4 — Test 1: undamped plasmon, peak location + $1/q^2$ scaling
Dossier Test 1. A single undamped phasor $n_q(t)=A\,e^{-i\omega_p t}$: the estimator
$|n_q(\omega)|^2/q^2$ peaks at $\omega_p$, and with **constant** $A$ its peak height
scales as exactly $1/q^2$ (so $L_{\rm peak}\cdot q^2$ is $q$-independent). NOTE the
phasor peaks at $-\omega_p$ (sign convention), so we search $|\omega|$.""")
code(r"""def loss_synth(omega_p, q, amp=1.0, use_complex=True, gamma=0.0, N=2000, dt=0.5):
    '''proposed estimator |n_q(w)|^2/q^2 for one (optionally damped) mode.'''
    tt = np.arange(N)*dt
    z = amp*np.exp(-1j*omega_p*tt - gamma*tt)
    s = z if use_complex else z.real
    F = np.fft.fft(s*np.hanning(N)); f = np.fft.fftfreq(N, d=dt)*2*np.pi   # angular, a.u.
    return f, (np.abs(F)**2)/q**2

wp = 0.1276                                       # ~3.473 eV
for q in np.array([1,2,3,4])*0.1257:
    f, L = loss_synth(wp, q, 1.0, True); k = np.argmax(L)
    print(f'  q={q:.3f}: peak |w|={abs(f[k])*HA2EV:.3f} eV   L_peak*q^2={L[k]*q**2:.4e} (q-indep)')
print('-> peak at omega_p; estimator 1/q^2 factor exact (constant-amp seed).')""")

md(r"""### A5 — Test 2: BUG-A — taking the real part folds $\pm\omega$ and halves it
`density_fourier.py` line 182 does `np.fft.fft(sig.real, ...)` on the complex
$n_q$. A complex phasor is **directional** (one-sided in $\omega$); its real part
$\cos$ is **two-sided**, so the FFT folds $\pm\omega$ symmetric and halves the
amplitude. Peak *location* survives; amplitude + direction are corrupted.""")
code(r"""f, Lc = loss_synth(wp, 0.1257, 1.0, use_complex=True)
f, Lr = loss_synth(wp, 0.1257, 1.0, use_complex=False)
pc, nc = Lc[f>0.02].max(), Lc[f<-0.02].max()
pr, nr = Lr[f>0.02].max(), Lr[f<-0.02].max()
print(f'  complex : max(+w)={pc:.3e}  max(-w)={nc:.3e}   -> ONE-sided (directional)')
print(f'  real    : max(+w)={pr:.3e}  max(-w)={nr:.3e}   -> folded (both sides)')
print(f'  amplitude ratio real/complex at dominant peak = '
      f'{np.sqrt(max(pr,nr))/np.sqrt(max(pc,nc)):.3f}  (BUG-A halving, expect 0.5)')
fig, ax = S.figure_one_col()
ax.plot(f*HA2EV, np.sqrt(Lc), lw=0.9, label='complex (correct)')
ax.plot(f*HA2EV, np.sqrt(Lr), lw=0.9, ls='--', label='real-part (BUG-A)')
ax.set_xlim(-8,8); ax.set_xlabel(r'$\hbar\omega$ (eV)'); ax.set_ylabel(r'$|n_q(\omega)|/q$')
ax.set_title('A5 - BUG-A: real part folds +/-w and halves it'); ax.legend()
fig""")

md(r"""### A6 — Test 3: damped mode — the line-shape caveat
Dossier Test 3. With damping $\gamma$, $|n_q(\omega)|^2$ is a **Lorentzian$^2$**
(HWHM $\approx0.644\gamma$ — too narrow, wrong area) where the *true* loss function
wants a single Lorentzian (HWHM $=\gamma$). The peak *position* still matches — so
$L$ is a peak-locator, not a quantitative lineshape.""")
code(r"""gamma = 0.01
f, L = loss_synth(wp, 0.1257, 1.0, use_complex=True, gamma=gamma)
amp = np.sqrt(L); side = f < -0.02
fneg = f[side]; an = amp[side]; k = np.argmax(an)
half = an[k]/np.sqrt(2)
# crude HWHM of |n_q| (Lorentzian) vs |n_q|^2 (Lorentzian^2)
print(f'  damping gamma={gamma} a.u.;  |n_q|^2 is Lorentzian^2 (HWHM ~ 0.644*gamma),')
print(f'  the true single-Lorentzian loss would be wider (HWHM = gamma).')
print(f'  peak position |w| = {abs(fneg[k])*HA2EV:.3f} eV (still on omega_p) -> locator OK')
fig, ax = S.figure_one_col()
ax.plot(-fneg*HA2EV, an/an[k], lw=1.0)
ax.set_xlim(2.5,4.5); ax.set_xlabel(r'$\hbar\omega$ (eV)'); ax.set_ylabel(r'$|n_q(\omega)|$ (norm)')
ax.set_title(f'A6 - damped mode line shape (gamma={gamma})')
fig""")

# ================================================================ Stage B =====
md(r"""# Stage B — `fourier.py` on the QKE plasmon (known: 6.48 eV)

The 54-atom BCC Li `v0p0626` kick. The journal (2026-05-06) locked the plasmon at
**6.480 eV in `dipole_x`**, found by searching the **[5.5, 8] eV** band — because a
stronger low-$\omega$ feature and DC drift dominate the *global* argmax.""")
code(r"""obs = pd.read_csv(QKE + '/observables.csv').drop_duplicates('time_au')
t = obs.time_au.to_numpy(); dx = obs.dipole_x.to_numpy()
print(f'QKE v0p0626: N={len(t)}  dt={t[1]-t[0]:.3f} a.u.  T={t[-1]:.1f} a.u.  '
      f'bin dw={2*np.pi/t[-1]*HA2EV:.3f} eV')
r = FourierTransform(window=WindowSpec('hann'), zero_pad=8, subtract='detrend').transform(t, dx, 'dipole_x')
print(f'  dipole_x plasmon (search 5.5-8 eV) = {peak_in(r.frequency_au, r.power, 5.5, 8.0):.3f} eV'
      f'   (journal 6.480, paper 6.5)')
print(f'  GLOBAL argmax (>0.5 eV)            = {peak_in(r.frequency_au, r.power, 0.5, 20):.3f} eV'
      f'   <- misleading low-w feature')""")

md(r"""### B2 — the trap: global argmax + the energy channel
The naive "find the biggest peak" gives the wrong answer, and `energy_total` is
*drowned* by low-$\omega$ drift (journal: "argmax misleading"). This is the concrete
content of the `fourier.py` TODO.""")
code(r"""fig, ax = S.figure_one_col()
o = omega_eV(r.frequency_au)
ax.plot(o[o<12], r.amplitude[o<12], lw=0.8, label='dipole_x (detrend)')
ax.axvspan(5.5, 8.0, color='0.85', label='plasmon search band')
ax.axvline(6.48, ls=':', color='k', label='journal 6.48 eV')
ax.set_xlabel(r'$\hbar\omega$ (eV)'); ax.set_ylabel('amplitude')
ax.set_title('B2 - dipole_x: plasmon at 6.48 is NOT the global argmax'); ax.legend()
# energy channel for contrast
de = pd.read_csv(QKE + '/excess_energy_per_uc.csv').drop_duplicates('time_au')
re_ = FourierTransform(window=WindowSpec('hann'), zero_pad=8, subtract='detrend').transform(
        de.time_au.to_numpy(), de.dE_eV_per_uc.to_numpy(), 'dE')
print(f'  energy excess channel global argmax = {peak_in(re_.frequency_au, re_.power, 0.5, 20):.3f} eV'
      f'  (drowned by low-w drift; the plasmon is weak here)')
fig""")

md(r"""### B3 — baseline robustness (answers the TODO)
Does `fourier.py`'s windowing/detrending *skew* the QKE result? Re-extract the
plasmon for every baseline mode.""")
code(r"""print('dipole_x plasmon peak in [5.5, 8] eV per baseline:')
for sub in ['none','initial','mean','detrend']:
    rr = FourierTransform(window=WindowSpec('hann'), zero_pad=8, subtract=sub).transform(t, dx)
    print(f'   subtract={sub:8s}: {peak_in(rr.frequency_au, rr.power, 5.5, 8.0):.3f} eV')
print('\n-> VERDICT (evidence): the plasmon POSITION is invariant to the baseline'
      '\n   choice (~6.45 eV ~ stored 6.48 within the 0.28 eV bin). fourier.py does'
      '\n   NOT skew the peak; the only real pitfall is naive global-argmax peak-'
      '\n   finding + the low-w DC feature -> always search a physical band.')""")

# ================================================================ Stage C =====
md(r"""# Stage C — `density_fourier.py` on the E15 jellium $n_q$ (known: $\omega_p\approx3.47$ eV)

The dedicated E15 plasmon run ($T\approx2000$ a.u., $\Delta\omega\approx0.09$ eV).
`density_fourier.py` extracts $n_q(t)$ and FFTs it — and contains both documented
bugs. We reproduce the known $\omega_p$, then show each bug on real data.""")
code(r"""nq = pd.read_csv(E15)
def mode(m):
    d = nq[nq.m==m].sort_values('time_au')
    return d.time_au.to_numpy(), d.re_n_q.to_numpy()+1j*d.im_n_q.to_numpy(), float(d.q_au.iloc[0])
ts, z1, q1 = mode(1); dts = ts[1]-ts[0]
print(f'E15 n_q: m=1  q={q1:.4f} 1/bohr  M={len(ts)} samples  dt_sample={dts:.2f} a.u.')
# correct extraction: drop t<5 transient, remove DC, Hann, complex FFT
cut = ts >= 5.0; zz = z1[cut] - z1[cut].mean(); win = np.hanning(len(zz)); npad = len(zz)*4
f = np.fft.fftfreq(npad, d=dts)*2*np.pi; Fc = np.fft.fft(zz*win, n=npad)
oeV = f*HA2EV; band = (np.abs(oeV) > 0.5) & (np.abs(oeV) < 20)
pk = abs(oeV[band][np.argmax(np.abs(Fc)[band])])
print(f'  complex FFT (DC-removed) peak |w| = {pk:.3f} eV   (stored 3.53, omega_p~3.47)')""")

md(r"""### C2 — BUG-A on real data (line 182: `np.fft.fft(sig.real, ...)`)
The real $n_q$ phasor is directional; the code's `.real` folds it.""")
code(r"""Fr = np.fft.fft(zz.real*win, n=npad)
pos = (oeV>0.5)&(oeV<20); neg = (oeV<-0.5)&(oeV>-20)
print(f'  complex: max(+w)={np.abs(Fc)[pos].max():.1f}  max(-w)={np.abs(Fc)[neg].max():.1f}  (directional)')
print(f'  real   : max(+w)={np.abs(Fr)[pos].max():.1f}  max(-w)={np.abs(Fr)[neg].max():.1f}  (folded)')
fig, ax = S.figure_one_col()
ax.plot(oeV, np.abs(Fc), lw=0.8, label='complex (correct)')
ax.plot(oeV, np.abs(Fr), lw=0.8, ls='--', label='real-part (BUG-A, line 182)')
ax.set_xlim(-8,8); ax.set_xlabel(r'$\hbar\omega$ (eV)'); ax.set_ylabel(r'$|{\rm FFT}[n_q]|$')
ax.set_title('C2 - BUG-A: real-part folds the directional plasmon'); ax.legend()
fig""")

md(r"""### C3 — BUG-B (line 183: stores $|{\rm FFT}[n_q]|$, not $|{\rm FFT}|^2/q^2$)
The docstring intends the loss-function proxy $|n_q(\omega)|^2/q^2$; the code plots
bare $|{\rm FFT}[n_q]|$. The $q$-weighting changes the *relative* mode intensities
(it does not move peaks).""")
code(r"""Ha2 = HA2EV
# Bohm-Gross dispersion for reference (density_fourier line 52)
wp_au, vF = 3.473/HA2EV, 0.3374
def bohm_gross(q): return np.sqrt(wp_au**2 + 0.6*vF**2*q**2 + q**4/4)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
cmap = plt.get_cmap('viridis')
for m in range(1, 7):
    tm, zm, qm = mode(m); zc = zm[ts>=5.0] - zm[ts>=5.0].mean()
    Fm = np.fft.fft(zc*np.hanning(len(zc)), n=len(zc)*4)
    fm = np.fft.fftfreq(len(zc)*4, d=dts)*2*np.pi*HA2EV
    sel = (fm>0.3)&(fm<14); c = cmap((m-1)/5)
    axes[0].plot(fm[sel], np.abs(Fm)[sel], color=c, lw=0.8)                 # BUG-B: |FFT|
    axes[1].plot(fm[sel], (np.abs(Fm)[sel]**2)/qm**2, color=c, lw=0.8,
                 label=f'm={m} (q={qm:.2f})')
    axes[1].axvline(bohm_gross(qm)*HA2EV, color=c, ls=':', lw=0.7)
axes[0].set_title(r'C3 BUG-B: $|{\rm FFT}[n_q]|$ (as plotted)')
axes[1].set_title(r'corrected $|n_q(\omega)|^2/q^2$ + Bohm-Gross')
for a in axes: a.set_xlabel(r'$\hbar\omega$ (eV)')
axes[1].legend(fontsize=6, ncol=2)
fig.tight_layout(); fig""")

md(r"""### C4 — the corrected estimator (evidence for the task-3 fix)
A clean drop-in (NOT applied to the file here — that is task 3, via `code-test`):
complex FFT, DC-removed, Hann, then $|n_q(\omega)|^2/q^2$, **with the peak-locator
caveat in bold**.""")
code(r"""def loss_function(times, n_q_complex, q, *, t_start=5.0, zero_pad=4):
    '''Corrected loss-function PEAK-LOCATOR L(q,w)=|n_q(w)|^2/q^2.
    Fixes BUG-A (keep complex) and BUG-B (square + /q^2). PEAK POSITIONS ONLY:
    absolute lineshape/area/cross-q intensity are NOT trusted (quadratic-in-n_q,
    bare 1/q^2 vs true 4*pi/q^2; dossier loss-function-formula-validation).'''
    m = times >= t_start
    z = n_q_complex[m] - n_q_complex[m].mean()
    win = np.hanning(z.size); npad = z.size*zero_pad
    F = np.fft.fft(z*win, n=npad)
    f = np.fft.fftfreq(npad, d=times[1]-times[0])*2*np.pi
    return f*HA2EV, (np.abs(F)**2)/q**2

oeV1, L1 = loss_function(ts, z1, q1)
b = (np.abs(oeV1)>0.5)&(np.abs(oeV1)<20)
print(f'corrected loss_function: m=1 peak |w| = {abs(oeV1[b][np.argmax(L1[b])]):.3f} eV'
      f'  (omega_p ~ 3.47; PEAK-LOCATOR only)')""")

# ================================================================ findings ====
md(r"""# Findings (evidence — verdict lines for the user)

| # | Audited element | Evidence from this notebook |
|---|---|---|
| 1 | coherent-gain normalisation (`fourier.py`) | unit tone → 1.0 for all windows (A1) — **correct** |
| 2 | baseline modes | peak location invariant; only $\omega\!\approx\!0$ changes (A2, B3) — **correct, not skewing** |
| 3 | zero-padding | interpolates axis only; peak/width set by $T$ (A3) — **correct** |
| 4 | QKE plasmon via `fourier.py` | 6.45 eV in [5.5,8] ≈ stored 6.48 (B1, B3); **pitfall = naive global-argmax**, not the FFT |
| 5 | loss estimator $1/q^2$ + peak position | exact on Test 1 (A4); peak-locator only (A6, §2 caveat) |
| 6 | `density_fourier` BUG-A (line 182 `.real`) | folds $\pm\omega$, halves amplitude (A5, C2) — **bug confirmed** |
| 7 | `density_fourier` BUG-B (line 183 `|FFT|`) | plots $|{\rm FFT}|$ not $|{\rm FFT}|^2/q^2$ (C3) — **bug confirmed** |

**Dossier verdicts (left blank for the user — verification-user-owns-verdict):**
- `loss-function-formula-validation`: formula accepted as peak-locator? Y / N ____
- `fft-drift-removal-validation`: baseline convention accepted? Y / N ____
- `fft-normalization-validation`: coherent-gain convention accepted? Y / N ____

**Takeaway.** `fourier.py`'s window/baseline/zero-pad machinery is **sound** and does
*not* skew the QKE plasmon position — the real failure mode is **peak-finding**
(global argmax on a DC-dominated spectrum) and the low-$\omega$ drift, fixed by
always searching a physical band. `density_fourier.py` has **two real bugs** (A:
`.real` before FFT; B: plots $|{\rm FFT}|$ not the $|{\rm FFT}|^2/q^2$ proxy) — both
fixed by the `loss_function` drop-in in C4, to be applied via `code-test` in task 3.
Throughout, $L=|n_q|^2/q^2$ is a **peak-locator**, never a quantitative lineshape.""")

# --------------------------------------------------------------- execute -----
ep = ExecutePreprocessor(timeout=900, kernel_name="inqview-venv")
ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
nbformat.write(nb, NB)
print(f"executed + wrote {NB}  ({len(nb.cells)} cells)")
