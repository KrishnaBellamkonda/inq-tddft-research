# Loss-function formula validation (independent)

Independent fresh-context re-derivation and sanity-check of the proposed
real-time-TDDFT jellium energy-loss-function estimator. This dossier produces
evidence neutrally; the **Verdict** section is left blank for the user to fill.

Conventions: atomic units (ℏ = m_e = e = 1, 4πε₀ = 1) throughout, matching the
repo's `inq-stack/python/inqview/postprocess/lindhard.py`. Coulomb in Fourier
space is `v(q) = 4π/q²`.

Labels: statements attributed to a source are cited inline; my own reasoning is
prefixed **Inference:**.

---

## Formula under review

A real-time TDDFT jellium run propagates `n(r,t)`. The proposed loss-function
estimator is

    L(q, ω) = |n_q(ω)|² / q²

with
- `δn(r,t) = n(r,t) − n(r,0)` (induced density),
- `n_q(t)` = 3D spatial FT of `δn(r,t)` at wavevector **q**,
- `n_q(ω)` = temporal FT of `n_q(t)`,
- `q = |q|`.

Two **q**-sampling modes: (a) "axial", **q** = (0,0,2πm/L_z); (b) "3d_binned",
all FFT modes binned by |q|.

As-implemented reference (the code this is meant to describe):
`/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/density_fourier.py`.
The analytic ground-truth dielectric module is
`/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/lindhard.py`
(`loss_function(q,ω,kF) = Im[−1/ε_RPA]`, `vq(q) = 4π/q²`, `epsilon_rpa = 1 − v(q)χ⁰`).

---

## Derivation: |n_q(ω)|² vs −Im[1/ε(q,ω)]

### Step 1 — Linear response of the density

For a homogeneous electron gas perturbed by a weak external scalar potential
`V_ext(r,t)`, linear response gives, mode by mode,

    δn(q, ω) = χ(q, ω) · V_ext(q, ω).                                  (1)

Here χ is the **full** (interacting, density-density) response function. This is
the standard Kubo linear-response result (Giuliani & Vignale, *Quantum Theory of
the Electron Liquid*, §3.2–3.3; Pines & Nozières, *Theory of Quantum Liquids*
Vol. I, Ch. 2). The repo's own note
`docs/sources/dipole_as_q0_density_projection.md` states the same Eq. (1).

### Step 2 — χ, ε, and the loss function

The retarded density response and the inverse dielectric function are tied by

    1/ε(q, ω) = 1 + v(q) χ(q, ω),     v(q) = 4π/q².                     (2)

(Giuliani-Vignale Eq. 5.32 / Mahan §5; this is the *definition* of the test-
charge dielectric function in terms of the full proper response.) Taking the
imaginary part,

    Im[1/ε(q, ω)] = v(q) · Im χ(q, ω).                                  (3)

The energy-loss function is

    L_true(q, ω) ≡ −Im[1/ε(q, ω)] = −v(q) Im χ(q, ω) = −(4π/q²) Im χ.   (4)

So the loss function is **linear** in Im χ, and the Coulomb factor that appears
is exactly `v(q) = 4π/q²`. **Inference:** the `1/q²` in the proposed estimator is
plausibly an attempt to reproduce this `v(q)` factor — but `v(q) = 4π/q²` carries
a `4π`, so a bare `1/q²` differs from the true Coulomb factor by the constant
`4π` (a normalization-only discrepancy, harmless for *locating* peaks or checking
*relative* q-scaling, but wrong for absolute magnitude).

### Step 3 — From Im χ to the power spectrum |n_q(ω)|²

The induced-density power spectrum is, from Eq. (1),

    |δn(q, ω)|² = |χ(q, ω)|² · |V_ext(q, ω)|².                          (5)

This is **quadratic** in the response (|χ|²), whereas the loss function (Eq. 4)
is **linear** (Im χ). They are not the same object in general. Two things must
hold to connect them:

**(i) A flat (white) drive spectrum.** If the external perturbation is an
impulsive temporal delta-kick, `V_ext(r,t) = v₀(r) δ(t)`, then
`V_ext(q,ω) = ṽ₀(q)` is **ω-independent**. Then

    |δn(q, ω)|² = |ṽ₀(q)|² · |χ(q, ω)|².                                (6)

The kick makes the *drive spectrum* flat, so the ω-structure of |δn(q,ω)|² is the
ω-structure of |χ(q,ω)|². This is the Yabana–Bertsch real-time-TDDFT impulse
recipe (Yabana & Bertsch, *Phys. Rev. B* 54, 4484 (1996)); the repo references it
in `docs/handovers/projectile_loss_function_interference.md` and
`docs/sources/dipole_as_q0_density_projection.md`.

**(ii) |χ|² vs Im χ.** Even with a flat drive, |δn(q,ω)|² ∝ |χ|², NOT Im χ. The
correct linear-response spectral quantity is the dynamic structure factor
S(q,ω), tied to Im χ by the **fluctuation–dissipation theorem** (T = 0):

    S(q, ω) = −(1/π) Im χ(q, ω),   ω > 0.                               (7)

(Giuliani-Vignale Eq. 3.110; Pines & Nozières Vol. I Ch. 2; confirmed in standard
electron-gas literature — see Sources.) Combining (4) and (7):

    L_true(q, ω) = −v(q) Im χ = π v(q) S(q, ω) = (4π²/q²) S(q, ω).      (8)

Equivalently the relation quoted in the structure-factor literature,
`S(q,ω) = −(1/[n v(q)]) Im[1/ε]`, rearranges to the same `L = Im[−1/ε] =
n v(q) S`-type proportionality (the `n` vs `π` prefactor difference is a
convention/normalization of how S is defined; both are linear in Im χ).

**Inference (the crux).** `|χ|²` (a Lorentzian-*squared*, peak ∝ 1/Γ², width set
by Γ/√(√2−1) ≈ 1.55Γ) is NOT `Im χ` (a single Lorentzian, peak ∝ 1/Γ, FWHM = 2Γ)
except in the **undamped / sharp-resonance limit** where both collapse to a
δ-function at the same pole ω = ω_res. So:

- **Peak *location*** of `|n_q(ω)|²` and of `L_true` coincide at every pole of χ
  (plasmon, e-h edges), because a pole of χ is a pole of both |χ|² and Im χ. For
  *plasmon detection* (locating ω_p(q)) the estimator is therefore defensible.
- **Peak *shape / width / weight*** do NOT coincide once damping is finite: a
  damped mode gives |χ|² ∝ Lorentzian² (too narrow, wrong area) where L_true wants
  Lorentzian¹. So `|n_q(ω)|²/q²` is **not** a quantitatively correct loss
  function off the undamped limit; it is a peak-*locator*, not the spectral
  density `Im[−1/ε]`.

### Summary of the |n_q(ω)|²/q² claim

| Aspect | Verdict-relevant finding |
|---|---|
| `1/q²` = Coulomb factor? | The *true* factor is `v(q)=4π/q²`. Bare `1/q²` is off by the constant `4π` (normalization only). |
| Power `|·|²` correct? | **No** in general. The loss function is linear in Im χ (∝ S(q,ω)), not quadratic in χ. `|χ|²` only matches at sharp/undamped poles. |
| When is `|n_q(ω)|² ∝ Im χ` valid? | Only when the drive is a flat (delta-kick) spectrum AND the mode is (near-)undamped, so |χ|² and Im χ share the same δ-peak. For finite damping the line shapes differ (Lorentzian² vs Lorentzian). |
| Fit for purpose? | **Inference:** good as a *plasmon-peak locator* (right pole positions, right qualitative q-trend); **not** a quantitatively faithful `−Im[1/ε]` (wrong line shape, wrong area, wrong absolute normalization). |

---

## Longitudinal vs transverse / axial vs 3d_binned

**Density response is intrinsically longitudinal.** The density-density response
χ(q,ω) couples to a *scalar* (longitudinal) external potential and produces a
*charge-density fluctuation* (a compression). The loss function `−Im[1/ε_L]`
built from it is the **longitudinal** loss function — it is the inelastic channel
probed by EELS / IXS / a moving charge. This is standard (Giuliani-Vignale §4–5;
Pines & Nozières Vol. I): plasmons are *longitudinal* collective charge
oscillations. **Confirmed:** there is no "transverse density" loss channel — a
transverse current response χ_T exists and has its own ε_T(q,ω), but it couples to
the *transverse vector potential* (light/magnetic probes), does not modulate the
scalar density, and is irrelevant to a moving-charge / kick stopping-power
calculation. **Inference:** for jellium stopping power and plasmon detection from
δn, only the longitudinal `−Im[1/ε_L]` is in play; the estimator correctly uses
the density and so is automatically the longitudinal channel.

**Axial vs 3d_binned is sampling, not physics.** Jellium is isotropic, so the
*exact* loss function depends only on |q|: L(q,ω) = L(|q|,ω). Therefore:

- "axial" (**q** along z) and "3d_binned" (all directions binned by |q|) sample
  the **same** physical function L(|q|,ω). The distinction is purely a matter of
  **|q|-sampling and statistics**, not anisotropy.
- **3d_binned** averages all FFT modes on a |q|-shell → more samples per |q| bin,
  better signal-to-noise, denser |q| coverage (it reaches |q| values off the
  z-axis grid lines). For an isotropic system this directional average is exact
  and strictly improves statistics.
- **axial** restricts to **q** = (0,0,2πm/L_z): fewer samples, coarse |q| grid
  (only the discrete axial harmonics), but it is the natural sampling when the
  *perturbation itself* is along z. **Inference:** the only physics caveat is that
  a *finite, non-isotropic perturbation* (a z-directed WP) breaks isotropy of the
  *excited state* even though χ is isotropic; off-axis modes then carry a
  different *excitation amplitude*, so 3d_binned mixes shells that were driven
  with different strengths. For pure χ (kick-limit, weak) this is irrelevant; for
  a strong directional WP the two modes can legitimately differ in *amplitude*
  (not in pole position). The code's own TODO at
  `density_fourier.py:28–34` flags exactly this open question.

---

## Proposed reduced analytic test system (numpy pseudocode + expected values + tolerance)

Goal: a **portable, machine-independent, pure-numpy** unit test with a
**closed-form** ground truth, no INQ engine. I evaluate the three candidates,
then specify the chosen one fully.

**Candidate evaluation.**
- (i) **Single undamped plasmon** `n_q(t) ∝ e^{−iω_p t}`. Closed form: |n_q(ω)|²
  is a δ-peak (→ a sinc²-broadened peak on a finite grid) at exactly ω_p, with an
  exact `1/q²` scaling if amplitudes are seeded ∝ q^{−1}. Trivial, exact,
  machine-independent, and it *directly* exercises (a) peak location, (b) the
  `1/q²` scaling, (c) the BUG-A real-part folding. **Chosen.**
- (ii) **Analytic Lindhard at a few (q,ω).** Already covered by the repo's
  `postprocess/test_lindhard.py`; it validates `lindhard.py` (the ground-truth
  dielectric), not the *FFT estimator*. Keep as a *separate* check of the analytic
  reference, not of the time-FFT pipeline. (Not chosen for this test.)
- (iii) **Damped oscillator → Lorentzian.** Best for exposing the |χ|²-vs-Im χ
  line-shape discrepancy (BUG-relevant), because a damped mode makes |n_q(ω)|²
  a Lorentzian² while −Im[1/ε] wants a Lorentzian. **Chosen as a second, optional
  test** specifically to lock the line-shape caveat.

### Test 1 (primary): undamped plasmon — peak location, 1/q² scaling, BUG-A

```python
import numpy as np

def synth_n_q_t(omega_res, amp, t):
    # single undamped complex phasor (induced density of one collective mode)
    return amp * np.exp(-1j * omega_res * t)

def loss_estimator(n_q_t, q, dt, use_complex=True):
    # mirrors the proposed estimator L = |n_q(omega)|^2 / q^2
    sig = n_q_t if use_complex else n_q_t.real      # BUG-A toggle
    F = np.fft.fft(sig)
    freq = np.fft.fftfreq(len(sig), d=dt)
    omega = 2*np.pi*freq
    return omega, (np.abs(F)**2) / q**2

# --- fixed inputs (machine-independent) ---
omega_p = 0.1276          # a.u.  (L=50, N=162 plasma frequency; matches code default)
T  = 4096                 # samples
dt = 1.0                  # a.u. per sample  -> Nyquist omega = pi >> omega_p
t  = np.arange(T) * dt
q_list  = np.array([0.10, 0.20, 0.40])     # three |q| values (a.u.)
amp_list = 1.0 / q_list                    # seed amplitude ~ 1/q  (see note)

omega_peaks, peak_heights = [], []
for q, amp in zip(q_list, amp_list):
    nqt = synth_n_q_t(omega_p, amp, t)     # NB: same omega_p for all q (no dispersion in toy)
    om, L = loss_estimator(nqt, q, dt, use_complex=True)
    pos = (om >= 0)                         # positive-frequency half
    k = np.argmax(L[pos])
    omega_peaks.append(om[pos][k])
    peak_heights.append(L[pos][k])
```

**Exact expected features** (finite-grid, undamped):

1. **Peak location.** `omega_peaks[i] ≈ omega_p` for every q. With `dt=1`,
   `T=4096`, the FFT bin spacing is `Δω = 2π/(T·dt) = 2π/4096 ≈ 1.534e-3` a.u.
   The nearest bin to ω_p = 0.1276 is within ½·Δω ≈ 7.7e-4.
   - Assert: `np.allclose(omega_peaks, omega_p, atol=2e-3)` (≈ one bin; loose
     because ω_p need not land exactly on a bin). Machine-independent: a complex
     exponential FFT puts essentially all weight in the nearest 1–2 bins
     regardless of platform.

2. **1/q² scaling of the loss estimator.** With `amp = 1/q`, the time-FFT peak
   height of `|n_q(ω)|²` scales as `|amp|² = 1/q²`; dividing by `q²` in the
   estimator gives total `1/q⁴`. **Inference / design choice:** to test the
   *estimator's own `1/q²` factor* cleanly, seed **constant** amplitude
   (`amp_list = np.ones_like(q_list)`) instead; then `|n_q(ω)|²` peak is
   q-independent and the estimator output peak scales as exactly `1/q²`:
   - Assert (constant-amp variant): for peaks `h_i` at `q_i`,
     `np.allclose(h_i * q_i**2, h_0 * q_0**2, rtol=1e-6)` — i.e. `h·q²` is
     constant across q. This is exact up to floating-point (same waveform, only
     the `1/q²` prefactor changes), so a tight `rtol=1e-6` is justified and still
     machine-independent (it is a ratio of identical FFTs).

3. **Peak height magnitude (sanity).** For a pure phasor of amplitude `A` over
   `N` samples, `|FFT|` at the resonant bin ≈ `A·N` (Parseval/DFT of a complex
   exponential on-grid; off-grid it leaks but stays O(A·N)). So
   `peak_height ≈ (A·N)² / q²`. Use only as an order-of-magnitude assert
   (`rtol=0.2`) because ω_p is generally off-grid (spectral leakage).

### Test 2 (BUG-A trigger): real-part-only folds ±ω

Reuse Test 1 but compare `use_complex=True` vs `use_complex=False`:

```python
om, L_cplx = loss_estimator(nqt, q, dt, use_complex=True)
_,  L_real = loss_estimator(nqt, q, dt, use_complex=False)
```

**Expected:**
- `L_cplx` has a **single** peak at `+ω_p`, and is ≈0 at `−ω_p`.
- `L_real = |FFT(Re n_q)|²/q²` has **two equal peaks**, at `+ω_p` *and* `−ω_p`,
  each ≈ ¼ the complex peak's height (since `Re[e^{−iω t}] = ½(e^{−iω t}+e^{+iω t})`).
  - Assert: `L_real[at +ω_p] ≈ L_real[at −ω_p]` (`rtol=1e-6`, symmetric by
    construction) and `L_real[+ω_p] ≈ 0.25 * L_cplx[+ω_p]` (`rtol=1e-3`,
    leakage-limited). This is the exact, machine-independent signature of BUG A.

### Test 3 (optional, line-shape caveat): damped oscillator

```python
gamma = 0.01            # a.u. damping, gamma << omega_p
nqt = amp*np.exp(-1j*omega_p*t)*np.exp(-gamma*t)   # decaying phasor
```
**Expected:** the *complex* FFT magnitude `|n_q(ω)|` is a **Lorentzian** of
HWHM = γ centred at ω_p; the estimator's `|n_q(ω)|²` is a **Lorentzian²**
(HWHM = γ·√(√2−1) ≈ 0.644γ — narrower). The *true* loss function for one damped
mode is a single Lorentzian (HWHM γ). Assert the estimator peak is at ω_p
(`atol=2·Δω`) but its FWHM is ≈ 0.64× the Im-χ Lorentzian FWHM — documenting that
`|n_q(ω)|²` is **not** the correct line shape for `−Im[1/ε]`. Tolerance on the
width ratio: `rtol=0.1` (finite-grid Lorentzian fit).

**Tolerance philosophy.** All asserts are either (a) ratios of identical FFTs
(tight `rtol≈1e-6`, exact up to IEEE-754, platform-independent) or (b) peak
positions on a known grid (`atol ≈ one bin = 2π/(T·dt)`). No bit-exact whole-array
comparison is used, so the test is portable across numpy/BLAS builds.

---

## Bug audit (BUG A, BUG B)

Audited against `density_fourier.py` as it stands on disk.

### BUG A — time-FFT taken on the real part only

`density_fourier.py:177–183`:

```python
sig = n_q[mask, m - 1] * win
# ... comment claims: "complex FFT gives the right spectrum."
full = np.fft.fft(sig.real, n=n_pad)     # <-- .real, NOT the complex signal
```

**Physics finding: this is a bug** (the code's *own* comment two lines above says
the complex FFT is the right thing, then takes `sig.real`). `n_q(t)` is a complex
phasor: a single mode behaves as `n_q(t) ∝ e^{−iω_res t}`. Taking the real part,
`Re[e^{−iω_res t}] = ½(e^{−iω_res t}+e^{+iω_res t})`, **injects a spurious +ω_res
component**, so `FFT(Re n_q)` is forced to be Hermitian-symmetric and folds the
+ω and −ω halves together. Consequences:

- The true single-peak spectrum at the physical sign of ω is split into a
  symmetric ±ω pair, each at half amplitude (quarter power). **Inference:** for a
  *symmetric* drive where physical content already sits at both ±ω this is
  cosmetically harmless after taking |·| on the positive half — but it (i) halves
  the recovered amplitude, (ii) discards the sign/direction information in
  `n_q(t)` (which distinguishes a forward-propagating density wave from a
  backward one), and (iii) can alias a genuine −ω feature onto a +ω feature.
- **When real-part-only is acceptable:** only if `n_q(t)` is already known to be
  real (e.g. you symmetrised, or you only ever inspect `|n_q(ω)|` on ω≥0 and accept
  the factor-of-2 amplitude loss and loss of propagation direction). It is NOT
  acceptable when the sign of ω (forward vs backward density wave) matters, nor for
  any quantitative amplitude/area. **Recommended fix:** `np.fft.fft(sig, n=n_pad)`
  on the complex `sig` and take the appropriate ω-half.

### BUG B — plots |n_q(ω)|, not |n_q(ω)|²/q²

`density_fourier.py:183` and the plot at `:202–204`:

```python
spectra[:, m - 1] = np.abs(full[:n_pad // 2 + 1])     # amplitude |n_q(omega)|
...
ax.plot(omega_eV[...], spectra[...])                  # plots |FFT[n_q]|
ax.set_ylabel(r"$|\mathrm{FFT}[n_{q_m}(t)]|$")
```

**Physics finding: confirmed — what is plotted is the amplitude `|n_q(ω)|`, not
the loss-function estimator `|n_q(ω)|²/q²`.** Dimensional/physical difference:

- `|n_q(ω)|` is the **amplitude** of the density Fourier component: linear in the
  perturbation, units of (density × volume × time) = (e · time). It is monotone in
  but not equal to the spectral weight.
- `|n_q(ω)|²/q²` is **quadratic** in the perturbation and carries the `1/q²`
  Coulomb-like weighting; it is the proposed loss-function proxy (∝ structure
  factor / Im χ at a sharp pole, per Eq. 8).

So the plotted quantity differs from the *proposed* estimator by (i) the square,
(ii) the `1/q²` factor. **Inference:** as a **peak locator** the difference is
immaterial — `|n_q(ω)|`, `|n_q(ω)|²`, and `|n_q(ω)|²/q²` all peak at the same
ω(q) for each fixed q, so the plasmon-dispersion read-off is unaffected. But the
plot is **mislabelled relative to its docstring/intent** ("loss function"): it is
neither `−Im[1/ε]` nor even the stated `|n_q|²/q²`. For any cross-q intensity
comparison (relative weight of m=1 vs m=4), the missing square and `1/q²` change
the result. **Recommended fix:** compute `spectra = (np.abs(full[...])**2) /
q_vals[m-1]**2` if the loss-function proxy is intended, and relabel; or relabel
the axis honestly as "amplitude |n_q(ω)|" if only peak location is sought.

---

## Sources cited

Authoritative:
- G. F. Giuliani & G. Vignale, *Quantum Theory of the Electron Liquid*, Cambridge
  (2005): §3.2–3.3 (linear response, Eq. 1), Eq. 3.110 (fluctuation–dissipation,
  Eq. 7), §4–5 + Eq. 5.32 (χ–ε relation, Eq. 2; longitudinal loss function,
  plasmon as Re ε=0). The repo's `lindhard.py` cites G-V Eq. 4.21 for χ⁰.
- D. Pines & P. Nozières, *The Theory of Quantum Liquids* Vol. I, Benjamin (1966):
  Ch. 2 (density fluctuations, structure factor, longitudinal nature of plasmons).
- J. Lindhard, *Mat. Fys. Medd. Dan. Vid. Selsk.* 28, no. 8 (1954): dielectric
  formalism, loss function in stopping power. (Repo: `docs/sources/stopping-power-formulae.md`.)
- K. Yabana & G. F. Bertsch, *Phys. Rev. B* 54, 4484 (1996): real-time-TDDFT
  impulse (delta-kick → flat drive spectrum), Eq. 6 reasoning.
- G. D. Mahan, *Many-Particle Physics*, 3rd ed., Ch. 5: Lindhard/RPA dielectric,
  inverse-dielectric loss function (cross-check of Eqs. 2–4).

Repo-internal (consistency cross-checks, not primary physics authority):
- `inq-stack/python/inqview/postprocess/lindhard.py` — analytic `L=Im[−1/ε]`,
  `v(q)=4π/q²`, RPA ε; the ground-truth dielectric module.
- `inq-stack/python/inqview/postprocess/density_fourier.py` — the audited estimator.
- `docs/sources/dipole_as_q0_density_projection.md` — Eq. (1), kick→χ identity.
- `docs/handovers/projectile_loss_function_interference.md` — Yabana-Bertsch
  Δω=2π/T convention; jellium ω_p/period numbers.
- `docs/sources/stopping-power-formulae.md` — Lindhard `S(v) = (2/πv²)∫(dq/q)∫ ω
  Im[−1/ε] dω`.

Confirmatory web search (standard result, not used as a primary citation):
the structure-factor / inverse-dielectric relation `S(q,ω) = −(1/[n v(q)]) Im[1/ε]`
is reproduced across the electron-gas literature (e.g. finite-T electron-liquid
and EELS reviews on arXiv), matching Eq. (8) here.
- https://arxiv.org/pdf/2603.17699 (finite-T electron liquid, S–ε relation)
- https://arxiv.org/pdf/physics/9905060 (loss function L=−Im ε⁻¹ in plasmas)

---

## Verdict (user — 2026-06-25)

- **formula accepted? Y** — accepted **as a plasmon-peak LOCATOR** with the bold
  caveat (it is NOT a quantitatively faithful `−Im[1/ε]`: wrong line shape, area,
  and absolute `4π` normalisation off the undamped limit). Use the **complex**
  `n_q(t)` (BUG-A) and compute `|n_q(ω)|²/q²` (BUG-B).
- **reduced-test accepted? Y.**
- Applied 2026-06-25 in `pipeline/density_fourier.py::loss_locator` via `code-test`
  (`test_density_fourier_loss.py`, 3/3; catalogue row added).
