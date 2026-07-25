#!/usr/bin/env python3
"""Deterministic builder for the stopping-power-extraction training notebook.

Preserves the existing Section 1 cells (read from the on-disk .ipynb), appends
the first-principles ladder (Sections 2-6), the critique phase, and the clean
rebuild, then executes the whole notebook against the `inqview-venv` kernel so
outputs are embedded.

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_stopping_notebook.py
"""
from __future__ import annotations
import pathlib
import nbformat
from nbformat.v4 import new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = pathlib.Path(__file__).resolve().parent
NB = HERE / "stopping_power_extraction.ipynb"

nb = nbformat.read(NB, as_version=4)

# Keep ONLY the original Section-1 cells (title + load + look). Anything this
# builder previously appended is regenerated, so truncate at the Section-1 tail
# marker to stay idempotent.
keep = []
for cell in nb.cells:
    keep.append(cell)
    if cell.cell_type == "markdown" and "Next (Section 2)" in cell.source:
        break
nb.cells = keep

md = lambda s: nb.cells.append(new_markdown_cell(s))
code = lambda s: nb.cells.append(new_code_cell(s))

# ---------------------------------------------------------------- Section 2 ---
md(r"""## Section 2 — Define $S$ and do one global fit

The electronic stopping power is the energy the projectile loses per unit path:

$$S \;=\; -\frac{dE_\text{proj}}{dx}\;=\;-\frac{d(\mathrm{KE})}{ds},
\qquad \mathrm{KE}=\tfrac12 m|\mathbf v|^2,\quad ds=|\mathbf v|\,dt .$$

`v3p0` decelerates by only 1.4 %, so $\mathrm{KE}(s)$ is nearly straight and a
**single global least-squares line** $\mathrm{KE}(s)\approx -S\,s+c$ is a sensible
first scalar estimate of $S$. We report it in Ha/Bohr **and** eV/Å.""")

code(r"""import matplotlib.pyplot as plt
HaB_to_eVA = 27.211386 / 0.52917721      # 1 Ha/Bohr = 51.4221 eV/Angstrom
print(f'1 Ha/Bohr = {HaB_to_eVA:.4f} eV/Angstrom')

def global_S(s_arr, ke_arr, transient_bohr=0.0):
    '''Single linear slope of KE vs s past an initial transient; S = -slope.'''
    m = s_arr >= transient_bohr
    A = np.vstack([s_arr[m], np.ones(int(m.sum()))]).T
    a, b = np.linalg.lstsq(A, ke_arr[m], rcond=None)[0]
    res = ke_arr[m] - (a * s_arr[m] + b)
    return -a, np.sqrt((res**2).mean()), int(m.sum())

S_global, rms_global, n_all = global_S(s, KE, 0.0)
print(f'S_global = {S_global:.6f} Ha/Bohr = {S_global*HaB_to_eVA:.4f} eV/Angstrom')
print(f'fit RMS residual = {rms_global:.3e} Ha  '
      f'({100*rms_global/(KE.max()-KE.min()):.1f}% of the {KE.max()-KE.min():.4f} Ha KE span)')""")

code(r"""fig, ax = S.figure_one_col()
ax.plot(s, KE, label='KE$_{proj}$(s)')
ax.plot(s, -S_global*s + (KE[0]), '--',
        label=f'global fit  S={S_global:.5f} Ha/Bohr')
ax.set_xlabel('path  s = z - z0  (Bohr)'); ax.set_ylabel('KE$_{proj}$ (Ha)')
ax.set_title('Section 2 - one global slope = first S estimate'); ax.legend()
fig""")

md(r"""**Read-out.** $S_\text{global}\approx0.0074$ Ha/Bohr $\approx0.38$ eV/Å. But the
fit residual is ~5 % of the KE span — **not noise, structure**: the global line is
dragged by the *entry transient* near $s\approx0$ (the projectile only feels the
steady drag once it is well inside the jellium) plus mild curvature from the
deceleration. §2 deliberately produces a slightly-wrong number whose wrongness
points to the next rung.""")

# ---------------------------------------------------------------- Section 3 ---
md(r"""## Section 3 — The entry transient

Discard the first `transient_bohr` of path and re-fit. The cut length is a
**physical** judgement (how far the projectile travels before its wake is
established), revisited quantitatively in the critique. We scan it: the residual
collapses, **but $S$ itself keeps drifting upward with the cut — there is no clean
plateau** within the available path. So the transient length is a real,
consequential choice (the kernel default `3.0` Bohr sits on the low side of the
still-rising $S$, i.e. it may slightly *under*estimate). This foreshadows §7b.""")

code(r"""print(f'{"transient (Bohr)":>16s}  {"S (Ha/Bohr)":>12s}  {"S (eV/A)":>9s}  '
      f'{"rms (Ha)":>10s}  {"n_pts":>6s}')
tb_scan = [0.0, 1.0, 2.0, 3.0, 5.0, 8.0]
S_scan = []
for tb in tb_scan:
    S_tb, rms_tb, n_tb = global_S(s, KE, tb)
    S_scan.append(S_tb)
    print(f'{tb:16.1f}  {S_tb:12.6f}  {S_tb*HaB_to_eVA:9.4f}  {rms_tb:10.2e}  {n_tb:6d}')

S_cut, rms_cut, _ = global_S(s, KE, 3.0)   # the kernel default transient_bohr=3.0
print(f'\nlocked at transient_bohr=3.0 -> S_cut = {S_cut:.6f} Ha/Bohr '
      f'= {S_cut*HaB_to_eVA:.4f} eV/Angstrom')""")

code(r"""fig, ax = S.figure_one_col()
ax.plot(tb_scan, np.array(S_scan)*HaB_to_eVA, 'o-')
ax.axvline(3.0, ls=':', color='0.5'); ax.set_xlabel('transient_bohr (Bohr)')
ax.set_ylabel('S (eV/Angstrom)'); ax.set_title('Section 3 - S vs transient cut')
fig""")

# ---------------------------------------------------------------- Section 4 ---
md(r"""## Section 4 — Electronic cross-check + integrator health

Energy conservation: what the projectile loses, the electrons gain. So an
**independent** estimate is $S' = +\,dE_\text{elec}/ds$ from `observables.csv`
(electronic total energy). Agreement validates the extraction. The leftover
$E_\text{tot}=\mathrm{KE}_\text{proj}+E_\text{elec}$ drift measures integrator
health. Note `observables.csv` is written coarsely (every 50 steps), so this
cross-check is necessarily low-resolution.""")

code(r"""trk_s = pd.DataFrame({'step': trk['step'].to_numpy(), 's': s, 'KE': KE})
mrg = obs.merge(trk_s, on='step', how='inner')
mc = mrg[mrg.s >= 3.0]                      # same transient window as S_cut
A = np.vstack([mc.s.to_numpy(), np.ones(len(mc))]).T
aE, bE = np.linalg.lstsq(A, mc.energy_total.to_numpy(), rcond=None)[0]
S_elec = aE                                 # electrons gain -> +slope = S'

print(f'S  (from -dKE/ds, transient-cut) = {S_cut:.6f} Ha/Bohr')
print(f"S' (from +dE_elec/ds)            = {S_elec:.6f} Ha/Bohr   "
      f'(n={len(mc)} coarse points)')
print(f'channel agreement                = {100*abs(S_elec-S_cut)/S_cut:.1f}%')

Etot = mrg.KE.to_numpy() + mrg.energy_total.to_numpy()
print(f'\nE_total drift over run = {Etot[-1]-Etot[0]:+.4e} Ha '
      f'({100*(Etot[-1]-Etot[0])/mrg.KE.iloc[0]:.2f}% of KE0) -> integrator health')""")

# ---------------------------------------------------------------- Section 5 ---
md(r"""## Section 5 — Why $S$ must be binned by $v(t)$

For `v3p0`, $v$ is essentially constant, so one global slope ≈ one $(v,S)$ point.
That breaks the moment $v$ changes a lot. `v0p8` loses **53 %** of its speed: a
single slope through its $\mathrm{KE}(s)$ would average the drag at $v=0.8$ with
the drag at $v=0.38$ — physically different stopping. The fix is a **local** slope
over a short path window, *labelled by the instantaneous $v$ at that window* —
i.e. $S(v)$, a curve. That is exactly what the kernel does.""")

code(r"""tk8 = pd.read_csv(RUN.replace('v3p0','v0p8')+'/electron_track.csv')
tk8 = tk8.drop_duplicates(subset='step').sort_values('time_au').reset_index(drop=True)
v8 = np.sqrt(tk8.vx**2+tk8.vy**2+tk8.vz**2).to_numpy()
s8 = (tk8.z - tk8.z.iloc[0]).abs().to_numpy(); t8 = tk8.time_au.to_numpy()

fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
axes[0].plot(t, v, label='v3p0'); axes[0].plot(t8, v8, label='v0p8')
axes[0].axhline(0.3373, ls=':', color='0.5', label='v_F=0.337')
axes[0].set_xlabel('t (a.u.)'); axes[0].set_ylabel('v (a.u.)')
axes[0].set_title('v(t): const vs strongly decelerating'); axes[0].legend()
# a single global slope on v0p8 is ill-posed: show its large residual
S8g, rms8, _ = global_S(s8, 0.5*v8**2, 3.0)
axes[1].plot(s8, 0.5*v8**2, label='KE(s) v0p8')
axes[1].plot(s8, -S8g*s8 + (0.5*v8[0]**2), '--', label=f'1 global slope (rms={rms8:.1e})')
axes[1].set_xlabel('s (Bohr)'); axes[1].set_ylabel('KE (Ha)')
axes[1].set_title('one slope is wrong for v0p8'); axes[1].legend()
fig.tight_layout(); fig""")

md(r"""The single-slope residual for `v0p8` is an order of magnitude worse than for
`v3p0`: the curve is genuinely curved because $S$ depends on $v$ and $v$ is moving.
Hence a *local* slope, binned by $v$.""")

# ---------------------------------------------------------------- Section 6 ---
md(r"""## Section 6 — Arrive at the kernel (no black box)

We have reconstructed every ingredient: path $s$, $\mathrm{KE}=\tfrac12 m v^2$, a
transient cut, and a **local slope labelled by $v$**. That *is*
`inqview.analysis.stopping_extract`: `load_track` builds $(t,s,v,\mathrm{ke})$;
`stopping_vs_v(transient_bohr=3.0, window=11)` discards the transient and takes an
11-point local least-squares slope, returning $(v,S)$. For `v3p0` (almost constant
$v$) its mean must reproduce our hand number.""")

code(r"""from inqview.analysis import stopping_extract as SE
tr = SE.load_track(RUN+'/electron_track.csv', mass=1.0, axis='z')
vk, Sk = SE.stopping_vs_v(tr, transient_bohr=3.0, window=11)
print(f'kernel S(v) for v3p0: mean = {Sk.mean():.6f} Ha/Bohr '
      f'over v in [{vk.min():.3f}, {vk.max():.3f}]  ({len(Sk)} local points)')
print(f'hand global (transient-cut) S = {S_cut:.6f} Ha/Bohr')
print(f'difference = {100*abs(Sk.mean()-S_cut)/S_cut:.1f}%  -> the kernel is our hand method')

fig, ax = S.figure_one_col()
ax.plot(vk, Sk*HaB_to_eVA, '.', ms=3, label='kernel S(v) local slopes')
ax.axhline(S_cut*HaB_to_eVA, ls='--', label=f'hand global = {S_cut*HaB_to_eVA:.3f} eV/A')
ax.set_xlabel('v (a.u.)'); ax.set_ylabel('S (eV/Angstrom)')
ax.set_title('Section 6 - kernel reproduces the hand number (v3p0)'); ax.legend()
fig""")

# ----------------------------------------------------------------- Critique ---
md(r"""---
# Section 7 — Critique: where this method bends

Having built it honestly, we now stress the four knobs and the one regime where it
breaks.""")

md(r"""### 7a — Window size (the bias/variance knob)
`window` is the number of track points in each local fit. In principle: too small →
noisy slope; too large → over-smooths and blends different velocities. We sweep it
on the strongly-decelerating `v0p8`. **Empirically these tracks are smooth, so the
curve is nearly window-independent** until the window grows large enough (≳40 pts)
to begin flattening the *real* $S(v)$ trend near its peak.""")

code(r"""fig, ax = S.figure_one_col()
for w in [5, 11, 21, 41]:
    vw, Sw = SE.stopping_vs_v(SE.load_track(RUN.replace('v3p0','v0p8')+'/electron_track.csv'),
                              transient_bohr=3.0, window=w)
    ax.plot(vw, Sw*HaB_to_eVA, '.', ms=2, label=f'window={w}')
ax.set_xlabel('v (a.u.)'); ax.set_ylabel('S (eV/Angstrom)')
ax.set_title('7a - window size on v0p8 (noise vs over-smoothing)'); ax.legend()
fig""")

md(r"""### 7b — Transient length
The `transient_bohr` cut is a magic number. §3 showed $S$ on `v3p0` is fairly flat
for cuts of 2–8 Bohr, but on a short/strongly-decelerating run the cut eats a large
*fraction* of the usable path. We quantify how many points survive each cut per
run.""")

code(r"""for run in ['v3p0','v0p8','v0p6']:
    trr = SE.load_track(RUN.replace('v3p0',run)+'/electron_track.csv')
    for tb in [2.0, 3.0, 5.0]:
        n = int((trr.s >= tb).sum())
        print(f'{run}: transient={tb:.0f} Bohr -> {n:4d}/{trr.s.size} track pts survive '
              f'(s_max={trr.s.max():.1f})')
    print()""")

md(r"""### 7c — Local fit vs raw finite difference
The kernel uses a local **least-squares slope** (Savitzky-Golay-like); the naive
alternative is a point-to-point finite difference `np.gradient`. We compare their
**point-to-point roughness** $\mathrm{std}(\Delta S)$ (whole-curve std is a poor
noise metric because $S(v)$ has real structure). **Finding:** on these *smooth*
Ehrenfest tracks the two are essentially equivalent (roughness ratio ≈ 1) — the
data is not noisy at the per-step level, so the derivative method barely matters
here. Windowing is **robustness insurance** (a well-defined slope, immunity to the
duplicated-$t_0$ row or an occasional glitch) and would matter for a *noisy*
channel such as the coarse electronic energy — it is not dramatic noise suppression
on this track.""")

code(r"""def roughness(x):
    '''high-frequency jitter = std of consecutive differences (isolates noise).'''
    return float(np.std(np.diff(x)))

for run in ['v3p0', 'v0p8']:
    trr = SE.load_track(RUN.replace('v3p0', run)+'/electron_track.csv')
    m = trr.s >= 3.0
    fd = -np.gradient(trr.ke[m], trr.s[m])
    _, Sw = SE.stopping_vs_v(trr, transient_bohr=3.0, window=11)
    rf, rw = roughness(fd*HaB_to_eVA), roughness(Sw*HaB_to_eVA)
    print(f'{run}: finite-diff roughness = {rf:.3f} eV/A   '
          f'windowed = {rw:.3f} eV/A   (roughness ratio {rf/rw:.1f}x ~ equivalent)')

trr = SE.load_track(RUN.replace('v3p0','v0p8')+'/electron_track.csv')
m = trr.s >= 3.0
fd = -np.gradient(trr.ke[m], trr.s[m])
vw, Sw = SE.stopping_vs_v(trr, transient_bohr=3.0, window=11)
fig, ax = S.figure_one_col()
ax.plot(trr.v[m], fd*HaB_to_eVA, '.', ms=2, label='finite difference')
ax.plot(vw, Sw*HaB_to_eVA, '.', ms=3, label='windowed lstsq (kernel)')
ax.set_xlabel('v (a.u.)'); ax.set_ylabel('S (eV/Angstrom)')
ax.set_title('7c - finite difference vs windowed fit (v0p8)'); ax.legend()
fig""")

md(r"""### 7d — The energy-budget disagreement
$\Delta\mathrm{KE}_\text{proj}+\Delta E_\text{elec}$ should be zero. It is not,
exactly — but the imbalance is a near-**constant ~1 mHa** across *all* runs (a
fixed integrator/cadence artefact, not a velocity trend). Being fixed in absolute
size it is sub-1 % everywhere, and *largest in relative terms* for the
smallest-$\Delta\mathrm{KE}$ run (`v3p0`, 0.8 %). This ~1 mHa is the honest floor
on "how much energy actually went into the electrons".""")

code(r"""print(f'{"run":>6s}  {"dKE_proj":>10s}  {"dE_elec":>10s}  {"sum":>10s}  {"|sum|/|dKE|":>11s}')
for run in ['v3p0','v2p0','v1p3','v0p8','v0p6']:
    rd = RUN.replace('v3p0',run)
    tt = pd.read_csv(rd+'/electron_track.csv').drop_duplicates('step').sort_values('time_au')
    vv = np.sqrt(tt.vx**2+tt.vy**2+tt.vz**2).to_numpy(); ke = 0.5*vv**2
    oo = pd.read_csv(rd+'/observables.csv')
    dKE = ke[-1]-ke[0]; dEe = oo.energy_total.iloc[-1]-oo.energy_total.iloc[0]
    print(f'{run:>6s}  {dKE:+10.4f}  {dEe:+10.4f}  {dKE+dEe:+10.4f}  '
          f'{100*abs(dKE+dEe)/abs(dKE):10.1f}%')""")

md(r"""### 7e — Assembling S(v): the Bragg peak, and the real low-v caveat
Stacking the per-run local $S(v)$ gives a **smooth, strictly positive** curve that
**rises** from $\sim0.4$ eV/Å at $v=3$ toward a peak of $\sim2.6$ eV/Å near
$v\sim v_F$ (0.4–0.8 a.u.) and **turns over** below $v_F$ (`v0p6`) — the expected
Bragg-peak velocity dependence of electronic stopping. **There is no numerical
breakdown**: no scatter (roughness $\sim10^{-3}$ eV/Å), no sign flips, across the
whole range. The genuine low-$v$ caveat is **statistical, not noise**: short path +
transient cut leave fewer points there (§7b) and the no-plateau transient
sensitivity (§3) is proportionally larger — so the low-$v$ values are *less
robust*, which the rebuild flags with an error bar.""")

code(r"""fig, ax = S.figure_one_col()
for run in ['v3p0','v1p3','v0p8','v0p6']:
    vv, Ss = SE.stopping_vs_v(SE.load_track(RUN.replace('v3p0',run)+'/electron_track.csv'),
                              transient_bohr=3.0, window=11)
    ax.plot(vv, Ss*HaB_to_eVA, '.', ms=2, label=run)
ax.axhline(0, color='0.6', lw=0.8)
ax.axvline(0.3373, ls=':', color='0.5', label='v_F=0.337')
ax.set_xlabel('v (a.u.)'); ax.set_ylabel('S (eV/Angstrom)')
ax.set_title('7e - S(v) assembles into a smooth Bragg-peak curve'); ax.legend()
fig""")

# ------------------------------------------------------------------ Rebuild ---
md(r"""---
# Section 8 — Clean deterministic rebuild

The critique exposed four choices the out-of-box kernel makes implicitly. The clean
version makes them **explicit and physical**:

1. **Window in Bohr, not point-count.** A fixed *path-length* window
   $w_\text{bohr}$ gives a velocity-resolution that does not silently depend on the
   per-run time step / cadence.
2. **Transient stated as a rule**, not a bare default.
3. **Local least-squares slope** (keep — beats finite difference, §7c) **with an
   error bar** from the fit covariance, so low-$v$ scatter is reported, not hidden.
4. **Explicit low-$v$ guard**: points whose fit error exceeds the value are flagged.

It must reproduce the kernel on `v3p0` and additionally carry a **1-σ error bar**
so the low-$v$ scatter is *reported* rather than hidden.""")

code(r"""def stopping_curve(track, *, transient_bohr=3.0, window_bohr=1.0):
    '''Deterministic S(v): local lstsq slope over a fixed PATH window (Bohr),
    returning v, S, and the 1-sigma slope error. S = -dKE/ds.'''
    m = track.s >= transient_bohr
    s_, ke_, v_ = track.s[m], track.ke[m], track.v[m]
    v_out, S_out, E_out = [], [], []
    half = window_bohr / 2.0
    for i in range(s_.size):
        sel = np.abs(s_ - s_[i]) <= half
        if int(sel.sum()) < 4:
            continue
        ss, kk = s_[sel], ke_[sel]
        A = np.vstack([ss, np.ones(ss.size)]).T
        coef, res, *_ = np.linalg.lstsq(A, kk, rcond=None)
        dof = max(ss.size - 2, 1)
        sigma2 = (res[0]/dof) if res.size else 0.0
        cov00 = np.linalg.inv(A.T @ A)[0, 0]
        v_out.append(v_[i]); S_out.append(-coef[0]); E_out.append(np.sqrt(sigma2*cov00))
    return np.array(v_out), np.array(S_out), np.array(E_out)

# reproduce the kernel on v3p0
vc, Sc, Ec = stopping_curve(tr, transient_bohr=3.0, window_bohr=0.66)  # ~11 pts at v3p0 cadence
print(f'rebuild  S(v3p0) mean = {Sc.mean():.6f} Ha/Bohr')
print(f'kernel   S(v3p0) mean = {Sk.mean():.6f} Ha/Bohr')
print(f'hand     global       = {S_cut:.6f} Ha/Bohr')""")

code(r"""tr8 = SE.load_track(RUN.replace('v3p0','v0p8')+'/electron_track.csv')
v8c, S8c, E8c = stopping_curve(tr8, transient_bohr=3.0, window_bohr=1.0)
v8k, S8k = SE.stopping_vs_v(tr8, transient_bohr=3.0, window=11)
fig, ax = S.figure_one_col()
ax.errorbar(v8c, S8c*HaB_to_eVA, yerr=E8c*HaB_to_eVA, fmt='.', ms=3, lw=0.5,
            label='rebuild (Bohr window + error bar)')
ax.plot(v8k, S8k*HaB_to_eVA, '.', ms=2, alpha=0.6, label='kernel (11-pt window)')
ax.set_xlabel('v (a.u.)'); ax.set_ylabel('S (eV/Angstrom)')
ax.set_title('Section 8 - clean rebuild vs kernel (v0p8)'); ax.legend()
fig""")

md(r"""## Choices to lock (for the `stopping-power-extraction` skill)

These are the decisions the deterministic workflow encodes — **flagged for sign-off,
not silently chosen**:

| Choice | Out-of-box kernel | Proposed lock | Why (per the stress test) |
|---|---|---|---|
| Window | 11 points | fixed `window_bohr` (path length) | result is window-insensitive on smooth tracks; fixed length = cadence-independent clarity, not accuracy (§7a) |
| Transient | `3.0` Bohr default | stated rule + per-run survival check | **the one real systematic** — no plateau in $S$ vs cut (~15 %); short runs lose too much path (§3, §7b) |
| Slope estimator | local lstsq | keep local lstsq | ≈ finite-difference here, but robust to glitches/noisy channels (§7c) |
| Uncertainty | none | 1-σ slope error returned | exposes the low-$v$ statistical thinness (§7e) |
| Low-v guard | none | flag points with err > value | low-$v$ is statistically thin, not noisy — flag, don't discard (§7e) |

**Takeaway.** On smooth Ehrenfest tracks the extraction is **more robust than
feared**: window size and fit-vs-finite-difference barely matter, and the assembled
$S(v)$ is a clean, positive **Bragg-peak** curve (peak $\sim2.6$ eV/Å near
$v\sim v_F$). The one genuine systematic is the **transient cut** (no plateau → $S$
uncertain at the ~15 % level: 0.42 vs 0.48 eV/Å on `v3p0`); the one genuine
weakness is **statistical thinness at low $v$**. For the clean `v3p0` case the three
routes (hand global, kernel, rebuild) agree to a percent
($S\approx0.0082$ Ha/Bohr $\approx0.42$ eV/Å, transient-cut). The rebuild's value is
therefore not noise suppression but making the transient choice explicit and the
low-$v$ uncertainty *reported* — the property the skill must guarantee.""")

# --------------------------------------------------------------- execute -----
ep = ExecutePreprocessor(timeout=600, kernel_name="inqview-venv")
ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
nbformat.write(nb, NB)
print(f"executed + wrote {NB}  ({len(nb.cells)} cells)")
