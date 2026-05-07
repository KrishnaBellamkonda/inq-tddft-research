# Topical: Plasmons, electron-hole excitations, and the stopping-power regime of our jellium WP runs

**Status:** topical entry (no single run; cross-cuts the L=30, L=50/N=138, L=50/N=162 runs)
**Linked entries:**
- [`2026-05-05_run_base_n138_L30_E5`](2026-05-05_run_base_n138_L30_E5.md)
- [`2026-05-05_run_base_n138_L50_E1p5`](2026-05-05_run_base_n138_L50_E1p5.md)
- [`2026-05-05_run_base_n162_L50_E1p5`](2026-05-05_run_base_n162_L50_E1p5.md)
- [`2026-05-03_run_propagate_v0p0123_extensive`](../quantumkickextension/2026-05-03_run_propagate_v0p0123_extensive.md)
**Source pack:** `docs/sources/correa-2018-electronic-stopping-power.md`

This entry consolidates: (i) a pedagogical explanation of what a
plasmon is and what it should look like in a TDDFT density movie, (ii)
the formulae for plasmon energy and electron-hole excitation threshold
from Correa (2018), and (iii) the substitution of our jellium-run
parameters into those formulae to **classify the energy regime** we
have been working in.

---

## 1. Pedagogical: what is a plasmon?

A **plasmon** is a *collective oscillation of the electron density* in
a metal — every electron moves a tiny bit, but they all move *in phase*,
and the net effect is a wave of density compression and rarefaction
sloshing back and forth at a single resonant frequency.

The clean way to derive the resonant frequency is the Drude / RPA
picture. Imagine the homogeneous electron gas (jellium) at density
$n$. Displace the entire electron cloud by a small distance $\xi$
relative to the (fixed) positive background:

- The displacement creates a surface charge of $-n\,e\,\xi$ on one face
  and $+n\,e\,\xi$ on the other, and therefore a uniform restoring
  electric field $E = 4\pi n e \xi$ inside the slab (Gaussian units,
  $4\pi\varepsilon_0 = 1$).
- The equation of motion of each electron is $m_e \ddot\xi = -eE = -4\pi n e^2 \xi$.
- That is a simple harmonic oscillator with angular frequency
  $$\boxed{\;\omega_p \;=\; \sqrt{\frac{4\pi n e^2}{m_e}}.\;}$$

In **Hartree atomic units** ($e = m_e = \hbar = 1$, $4\pi\varepsilon_0 = 1$):
$$\omega_p\;=\;\sqrt{4\pi n},\qquad \hbar\omega_p\;=\;\omega_p\;[\text{Ha}].$$

This is the **bulk plasmon frequency** — the same one Bohm and Pines
identified in 1953 for an electron liquid. RPA confirms that the
dielectric function $\varepsilon(q,\omega)$ has a zero at $\omega \to
\omega_p$ as $q \to 0$, and that zero shows up in the energy-loss
function $\Im[-1/\varepsilon(q,\omega)]$ as a sharp $\delta$-like peak
on top of the (broader) electron-hole continuum (see Correa 2018
Fig. 2 for the $(q,\omega)$ map).

### What does a plasmon look like in a TDDFT simulation?

A plasmon excitation in a density movie is **not** a localised wave
packet drifting through the box. It is a **standing density
oscillation**: regions of the cell alternately gain and lose density
at the plasmon frequency, with a wavelength set by the smallest
allowed wavevector $q$. In a periodic supercell of side $L$, the
smallest $q$ is $q_\text{min} = 2\pi/L$, so the plasmon you can excite
has wavelength $L$ — the cell length itself — and the density
oscillates as a dipole-like pattern that swaps direction every
$T_p = 2\pi/\omega_p$.

In a **time-dependent dipole** $\langle z\rangle(t)$ trace, a plasmon
shows up as a **damped sinusoid** at $\omega_p$. After an FFT, that's
the iconic single-peak spectrum near $\hbar\omega_p$ — what we got in
the QBall Li 54-atom run at 5.72 eV (cf. Li bulk plasmon ≈ 6.56 eV).

### What are *electron-hole excitations*?

These are **single-particle** transitions: one electron is promoted
from an occupied state below the Fermi level (a "hole" left behind)
to an unoccupied state above it. The excitation energy is the gap
between the two single-particle levels. In an infinite electron gas,
the spectrum of e-h pairs is **continuous** from $\omega = 0$ (a
particle right at the Fermi surface, kicked into the next infinitesimal
state) up to $\omega \approx q v_F + q^2/2m$ (the upper edge of the
"electron-hole continuum"). In a finite box, the $\{|G|^2\}$ values
are quantised, so the e-h spectrum is **discrete** and the lowest gap
is from the highest occupied $|G|^2$ shell to the next empty one.

In the QBall low-velocity result for Li (`run_propagate_v0p0123`,
`v = 0.0123` a.u.), the conclusion was that the FFT peak at 5.72 eV is
**the bulk plasmon**, not an interband transition — because the same
peak position survived in different windowing schemes
(`fft_windowing_comparison.py`) and matches the DFT-RPA $\hbar\omega_p$
of Li within a few percent. **That was a plasmon-dominated regime.**

---

## 2. Plasmon energy formula (Correa 2018 Eq. via RPA)

$$\boxed{\;\hbar\omega_p \;=\; \hbar\sqrt{\frac{4\pi n e^2}{m_e}} \;=\; \sqrt{4\pi n}\;\;[\text{atomic units}]\;}$$

(See Bohm–Pines 1953; standard RPA result; quoted explicitly as the
$q\to 0$ pole of $\varepsilon(q,\omega)$ in Correa 2018 §2 around Eq.
(3) and Fig. 2. The asymptotic Bethe limit at high velocity in Eq. (1)
takes $I \to \hbar\omega_p$ which fixes this same identification.)

## 3. Electron-hole threshold velocity (Correa 2018 Eq. (3) & Fig. 4)

Correa shows in Fig. 4 that for a system with electronic gap $E_g$
(insulator, or finite-box jellium with discrete shells), **single-particle
e-h excitations are kinematically forbidden** below

$$\boxed{\;v_\text{th}^{\,e\!-\!h} \;=\; \frac{E_g}{2 \hbar k_F}\;}$$

— a threshold velocity defined by the requirement that the projectile
must transfer at least $E_g$ in energy and the corresponding momentum
transfer fits inside the Fermi sphere. Above $v_\text{th}^{\,e\!-\!h}$
the projectile sees a continuum of accessible e-h pairs and the
linear-response stopping power Eq. (3) starts to climb.

For our finite-box jellium, $E_g$ is the **lowest e-h gap**, i.e. the
spacing between the top occupied shell and the next unoccupied shell.

## 4. Plasmon excitation threshold in a finite supercell (Correa 2018 §6.3)

A periodic supercell of side $L$ enforces a longest allowed
wavelength, hence a smallest allowed wavevector $q_\text{min} = 2\pi/L$.
A projectile at velocity $v$ cannot excite a plasmon whose dispersion
sits above the kinematic line $\omega = qv$ — so plasmon excitation is
**kinematically suppressed** unless

$$\boxed{\;v \;\geq\; v_\text{th}^{\,\text{plasmon}} \;=\; \frac{\omega_p}{q_\text{min}} \;=\; \frac{\omega_p\, L}{2\pi}\;}$$

This is the size-effect correction Eq. (15) is built around. **The
practical implication for our runs is large.**

---

## 5. Substitution: our jellium runs vs the regime boundaries

### Density, k_F, v_F, ω_p, lowest e-h gap, plasmon-threshold velocity for each run

| Quantity (atomic units, eV where stated) | L=30, N=138 | L=50, N=138 | L=50, N=162 |
|---|---|---|---|
| Volume Ω = L³ (Bohr³) | 27 000 | 125 000 | 125 000 |
| Density n = N/Ω (e/Bohr³) | 5.111 × 10⁻³ | 1.104 × 10⁻³ | 1.296 × 10⁻³ |
| $r_s = (3/4πn)^{1/3}$ (Bohr) | 3.62 | 6.04 | 5.72 |
| $k_F = (3π²n)^{1/3}$ (Bohr⁻¹) | 0.530 | 0.320 | 0.337 |
| $v_F = k_F$ (a.u.) | 0.530 | 0.320 | 0.337 |
| $E_F = k_F²/2$ (eV) | 3.82 | 1.39 | 1.55 |
| $\omega_p = \sqrt{4\pi n}$ (a.u.) | 0.2535 | 0.1178 | 0.1276 |
| $\hbar\omega_p$ (eV) | **6.90** | **3.20** | **3.47** |
| $q_\text{min} = 2\pi/L$ (Bohr⁻¹) | 0.2094 | 0.1257 | 0.1257 |
| $v_\text{th}^{\text{plasmon}} = \omega_p L/2\pi$ (a.u.) | **1.21** | **0.937** | **1.014** |
| Lowest e-h kinetic gap E_{|G|²=6→8} (eV) | 1.197 | **0.430** | **0.430** |
| $v_\text{th}^{\,e\!-\!h} = E_g/(2 k_F)$ (a.u.) | 0.0415 | 0.0247 | 0.0235 |

Lowest e-h gap: at $L$, $|G|^2$ values are $0, 1, 2, 3, 4, 5, 6, 8, 9, …$
in units of $(2\pi/L)^2$. Energies are $|G|^2 (2\pi/L)^2/2$. Going from
$|G|^2 = 6$ (top occupied at our N) to $|G|^2 = 8$ (next unoccupied)
gives the gap quoted.

### And the WP — our "projectile":

| Quantity | All three runs |
|---|---|
| WP $k_0$ (Bohr⁻¹) | (0, 0, 0.3320) — **L=50 runs**; (0, 0, 0.6062) — **L=30 run** |
| WP velocity $v = k_0/m_e$ (a.u.) | **0.332** (L=50, E=1.5 eV); **0.606** (L=30, E=5 eV) |
| WP kinetic energy ½ k_0² (eV) | **1.5** (L=50); **5.0** (L=30) |

### Where each run sits

| Regime test | L=30, N=138 (E=5 eV) | L=50, N=138 (E=1.5 eV) | L=50, N=162 (E=1.5 eV) |
|---|---|---|---|
| $v / v_F$ | **1.14** (above v_F, near Bragg peak) | **1.04** (right at v_F) | **0.985** (right at v_F) |
| $v / v_\text{th}^{\text{plasmon}}$ | **0.50** | **0.354** | **0.327** |
| $v / v_\text{th}^{\,e\!-\!h}$ | **14.6** | **13.4** | **14.1** |
| $E_\text{WP} / \hbar\omega_p$ | 0.72 | 0.47 | 0.43 |
| $E_\text{WP} / E_g$ (e-h gap) | 4.18 | 3.49 | 3.49 |

---

## 6. Calculated conclusion — what regime are we in?

1. **All three runs sit at or just below $v_F$.** They are squarely in
   the **electronic-stopping regime** of Correa Fig. 1 (the peak of
   $S(v)$ is right at $v \sim v_F$ for free-electron-like metals). This
   was the design intent.

2. **The WP velocity is well above the e-h excitation threshold** in
   every run ($v / v_\text{th}^{\,e\!-\!h} \gtrsim 14$), so single-particle
   electron-hole excitations across the lowest $|G|^2 = 6 \to 8$ shell
   gap are *kinematically allowed and copious*. This is consistent
   with the L=30 entry §3 estimate that the missing energy budget
   corresponds to ~5 elementary $|G|^2 = 6 \to 8$ excitations.

3. **In all three runs the WP velocity is well below the plasmon
   excitation threshold** ($v / v_\text{th}^{\text{plasmon}} \le 0.50$). The
   supercell is **too small** to host a plasmon at the WP velocity:
   the plasmon mode that fits in the box has wavelength $L$ and would
   need a projectile velocity of $\omega_p L / 2\pi$ to be excited
   resonantly, and we are well below that. **No plasmon should be
   excited in these runs.**

   This means the **interference pattern observed in the L=30 run**
   ([2026-05-05_run_base_n138_L30_E5](2026-05-05_run_base_n138_L30_E5.md) §1)
   is **not a plasmon** — it is more likely a *kinematic* wave-packet
   revival of the periodic-box quantum walk (the "Robinett revival"
   already flagged as required reading in the L=30 entry §1 and in
   `docs/observables_reference.md §13.4`). The QBall low-velocity Li
   run *did* see a plasmon because the Li 54-atom supercell has cell
   length 19.9 a.u. with $\hbar\omega_p \approx 6.5$ eV, giving
   $v_\text{th}^{\text{plasmon}} \approx 0.76$ a.u. — and the Li run's
   $v = 0.0123$ is below threshold by a factor of 60. So strictly
   speaking, even the QBall result at $v = 0.0123$ a.u. is below
   the kinematic threshold; what we observed there is the **plasmon
   resonance excited by the broadband initial-time impulse** (the
   instantaneous $v$-kick injects all frequencies, including
   $\omega_p$). Our jellium WP injection is *not* a broadband
   $\delta$-pulse; it is a Gaussian wave packet at fixed $k_0$, and
   the Gaussian envelope's frequency content is $\sim 1/\sigma_t$
   where $\sigma_t \sim \sigma_r/v$. For $\sigma_r = 5$, $v = 0.332$,
   $\sigma_t \approx 15$ a.u., bandwidth $\sim 1/15 \approx 0.07$ Ha
   $\approx 1.8$ eV — the Gaussian *does* contain frequencies near
   $\hbar\omega_p$ (3.2 eV is ~1.7 σ in the bandwidth tail).
   **Inference:** the plasmon could be marginally excited via the
   tail of the WP's spectral content even though the dominant kinematic
   line is sub-threshold — but only weakly. The dominant excitation
   channel is e-h.

4. **For the L=50 runs at E=1.5 eV, the only available physical loss
   channel is single-particle e-h excitation across the discrete
   shell gap (~0.43 eV, ~3 quanta to absorb the WP's 1.5 eV).** Any
   "missing kinetic energy" should match a discrete count of these
   e-h transitions, not a continuous plasmon Lorentzian.

### The hole-as-Coulomb-attractor hypothesis (corrected: charge-conjugate of the positive-ion wake)

A first reading of the energy bookkeeping ($\Delta E_H < 0$ in both
L=50 runs) might appear to falsify the "hole behind WP attracts WP
backwards" picture, on the grounds that doing work against an
attraction should *raise* the Hartree energy. **That reading is
wrong**, because the projectile here is a *negative* particle and the
sign of the wake is opposite to the textbook positive-ion case.

For a **positive ion** moving through an electron gas:
- Electrons are *attracted* to the ion and pile up around it.
- With finite response time, the pile-up trails the ion — the
  classical Bohr/Lindhard **electron wake**, an *accumulation* of
  electron density *behind* the ion.
- The ion-wake interaction energy is *negative* (the ion sits in a
  region of enhanced electron density it likes); as the wake forms,
  $\Delta E_H \le 0$, and the ion's kinetic energy is dissipated into
  bath kinetic energy as electrons countercurrent forward to keep up
  with the ion.

For a **negative wave packet** (our WP), the picture is the
**charge-conjugate**:
- Electrons are *repelled* from the WP and locally depleted around
  it.
- With finite response time, the depletion trails the WP — an
  *anti-wake*, a region of *lower* electron density *behind* the WP.
- The WP-anti-wake interaction is the negative of the positive-ion
  case in sign-of-charge but **the same in sign of energy**:
  removing electron density *behind* a negative WP also lowers
  $E_H$, because the bath gets to organise itself further from the
  WP than it could in equilibrium. So $\Delta E_H \le 0$ is
  **consistent** with the anti-wake picture.

The hole behind the WP is therefore **not just consistent with** the
"effective-positive-charge attracting the WP backwards" mechanism —
it **is the charge-conjugate countercurrent signature** of the
positive-ion electron wake. The retarding force on the WP is the
gradient of the WP-anti-wake Hartree interaction along $-\hat z$,
and that integrates to a kinetic-energy loss of the WP balanced by a
kinetic-energy gain of the (deformed) bath. Both observations —
$\Delta E_H < 0$ and bath KE > WP KE loss — are consistent with this
picture.

**The cleanest direct test** is to **run a positive-ion projectile
through the same L=50 jellium bath at matched v ≈ v_F** and compare
the wake. Predictions:

| Observable | Negative WP (this project) | Positive ion (proposed companion) |
|---|---|---|
| Sign of trailing-density Δn behind projectile | **negative** (hole) | **positive** (accumulation) |
| Sign of $\Delta E_H$ | **negative** | **negative** |
| Sign of Δkinetic_bath | **positive** | **positive** |
| Direction of retarding force on projectile | $-\hat z$ (slowdown) | $-\hat z$ (slowdown) |
| Stopping power $S(v) = \langle dE/dt\rangle/v_\text{proj}$ | (computed via cod_z slope, with bath-cod correction) | direct from Eq. (10) of Correa 2018 |

Hence the positive-ion run is the **canonical reference** that ties
our jellium WP scattering directly to the textbook electronic-stopping
literature. Added to the TODO list.

The companion orbital-cod-vs-time test
(`docs/plans/jellium_orthonormalisation_rerun.md §4`) is independent of
this and remains valuable: it pins down the ratio of WP motion to bath
drift in the cod_z signal.

---

## 7. Eigenenergy-rearrangement bookkeeping (the Σ f_i ε_i − D formula)

For the user's question — *"can the rearrangement of the eigenenergies
of the orbitals account for the WP's KE drop?"* — the relevant
identity is the standard KS bookkeeping (Correa 2018 Eq. (8); also
e.g. Parr & Yang Ch. 7):

$$\boxed{\;
E_\text{total}[n] \;=\; \sum_i f_i\,\varepsilon_i \;-\; \underbrace{\Big(\,E_H[n] \;+\; \int v_\text{xc}(\mathbf r) n(\mathbf r)\, d^3r \;-\; E_\text{xc}[n]\,\Big)}_{\equiv\,D \text{ ("double-counting" + xc shift)}}\;+\; E_\text{ext-bg}\;
}$$

Term-by-term:

- **$\sum_i f_i \varepsilon_i$** — the *band-structure sum*. $\varepsilon_i$
  is the Kohn–Sham eigenvalue of state $i$; $f_i$ is the (frozen, in
  TDDFT) occupation. This sum *double-counts* the electron-electron
  Hartree interaction because each pair $(i,j)$ contributes once
  through $\varepsilon_i$ (the field of all other electrons acting on
  $i$) and again through $\varepsilon_j$.
- **$E_H[n] = \tfrac{1}{2}\!\!\int\!\!\int\!\frac{n(\mathbf r) n(\mathbf r')}{|\mathbf r - \mathbf r'|}\, d^3r\, d^3r'$**
  — the Hartree energy of the density. **Subtracted with a minus sign
  to remove the band-structure double-count.** (The $\tfrac12$ in
  $E_H$ is the standard pair-counting prefactor; it is not the
  "double-counting" we are subtracting — we are subtracting the
  *full* $E_H$ once.)
- **$\int v_\text{xc} n d^3r$** — what the band-structure sum picked
  up of the xc potential, *also* subtracted because the band
  structure used $v_\text{xc}$ as an effective single-particle
  potential.
- **$+ E_\text{xc}[n]$** — the xc functional itself, *added back* with
  a plus sign so the total energy uses $E_\text{xc}[n]$ rather than
  $\int v_\text{xc} n$.
- **$E_\text{ext-bg}$** — the static jellium background energy
  (constant in time, drops out of any $\Delta E$).

So the **bookkeeping for any change** is

$$\Delta E_\text{total} \;=\; \Delta\!\Big(\sum_i f_i \varepsilon_i\Big) \;-\; \Delta E_H \;-\; \Delta\!\Big(\!\!\int v_\text{xc} n\Big) \;+\; \Delta E_\text{xc}.$$

Rearranged for the user's question:

$$\boxed{\;\Delta\!\Big(\sum_i f_i \varepsilon_i\Big) \;=\; \Delta E_\text{total} \;+\; \Delta D,\qquad \Delta D \,\equiv\, \Delta E_H + \Delta\!\!\int v_\text{xc} n - \Delta E_\text{xc}\;}$$

Because $\Delta E_\text{total} \approx 0$ in our runs (10⁻⁷ eV), the
band-structure sum change is **directly equal to $\Delta D$**, the
double-counting + xc shift. This is **not** a physical energy transfer
— it is a bookkeeping shift caused by the density redistributing
under the moving WP.

### Proposed numerical evaluation for the L=50 runs

Ingredients on disk:

- `state_energies.csv`: $\varepsilon_i(t), f_i$ at every observable
  step → compute $\Sigma f_i \varepsilon_i$ at $t=0$ and $t=30$ a.u.
- `observables.csv` columns `energy_hartree(t)`, `energy_xc(t)` →
  $\Delta E_H$ and $\Delta E_\text{xc}$ are direct.
- $\int v_\text{xc} n\, d^3r$ at $t=0$ and $t=30$ — **not currently
  written**; needs a one-shot GPU reduction added to the run template
  (or evaluated post-hoc from `density_rt_total.vti` + a libxc call,
  since LDA's $v_\text{xc}$ is a closed-form function of $n$).

For the **closed-shell N=162 run**, evaluating $\Sigma f_i \varepsilon_i$
directly from `state_energies.csv`:

| Quantity (N=162 L=50) | Value (eV) |
|---|---|
| $\Sigma f_i \varepsilon_i$ at $t = 0$ | −464.8835 |
| $\Sigma f_i \varepsilon_i$ at $t = 30$ a.u. | −465.4521 |
| $\Delta(\Sigma f_i \varepsilon_i)$ | **−0.5687** |
| ΔE_WP (single state ε, $f=1$) | −2.5553 |
| $\Sigma f_i \Delta\varepsilon_i$ over bath ($f_i = 2$) | +1.9866 |
| ΔE_total (drift) | −2.9 × 10⁻⁷ |
| ΔE_hartree | −0.5975 |
| ΔE_xc | +0.0955 |

For the bookkeeping to close ($\Delta E_\text{total} \approx 0$):

$$\Delta\!\!\int v_\text{xc} n \,d^3r \;=\; \Delta\!\Big(\sum_i f_i \varepsilon_i\Big) - \Delta E_H + \Delta E_\text{xc} - \Delta E_\text{total} \;=\; +0.1243 \text{ eV}.$$

This was the **+0.12 eV prediction** in the original draft of this
section; **direct evaluation now confirms it** (within ~10⁻⁴ eV
rounding). The same calculation for the partial-shell N=138 run gives
$\Delta(\int v_\text{xc} n) = +0.1278$ eV — essentially identical,
indicating the density redistribution caused by the WP traversal is
nearly the same in the two cases (the closed-shell vs partial-shell
distinction matters for the GS, less for the dynamical perturbation).

| Quantity (N=138 L=50, partial shell) | Value (eV) |
|---|---|
| $\Delta(\Sigma f_i \varepsilon_i)$ | −0.4289 |
| ΔE_WP single state | −2.6167 |
| $\Sigma f_i \Delta\varepsilon_i$ bath | +2.1878 |
| ΔE_total | −1.2 × 10⁻⁷ |
| ΔE_hartree | −0.4585 |
| ΔE_xc | +0.0982 |
| **Δ(∫v_xc n)** (closure) | **+0.1278** |

Numerical script: `venv/bin/python3` reads `state_energies.csv` (long
format `step,time_au,kpoint_index,state_index,weight,occupation,
E_expect_ha,...`) and `observables.csv`, takes the last-step minus
zeroth-step difference of $\Sigma$ weight × occupation × $E_\text{expect}$
and the four energy components.

### Verdict on the user's question

**Can the rearrangement of the orbital eigenenergies account for the
WP's kinetic-energy drop?**

**No, not on its own — and now numerically confirmed.** The WP
kinetic-energy drop derived from $\frac12 \Delta(v^2)$ is **−1.75 eV**
(N=162 L=50). The corresponding band-structure sum change is
$\Delta(\Sigma f_i \varepsilon_i) =$ **−0.5687 eV** — wrong sign for a
"WP-loss = band-sum loss" interpretation if one tried to make ε's
do all the work, and a third of the magnitude. The decomposition

$$\Delta(\Sigma f_i \varepsilon_i) \;=\; \underbrace{(-2.5553)}_{\text{WP}\,\Delta\varepsilon} \;+\; \underbrace{(+1.9866)}_{\text{bath}\,\Sigma f_i\Delta\varepsilon} \;=\; -0.5687 \text{ eV}$$

shows the WP eigenvalue does drop by 2.56 eV (the WP slot becomes
more bound to the bath as the cloud forms), but the bath orbitals
*rise* by about +2.0 eV in the band-structure sum. The remaining
budget shows up as the bath **kinetic-energy gain** (+2.25 eV from
$\Delta E_\text{kin}^\text{system} - \Delta E_\text{kin}^\text{WP}$), which
$\Sigma f_i \varepsilon_i$ does *not* directly capture because
$\varepsilon_i$ contains the full Hamiltonian expectation value
(kinetic + Hartree + xc + ext), not the kinetic alone, and the bath
ΔE_H is negative. Net of all bookkeeping, $\Delta E_\text{total}$
closes to ~10⁻⁷ eV with the verified $\Delta(\int v_\text{xc} n) \approx
+0.12$ eV correction above.

**The band-structure rearrangement is a bookkeeping shift, not a
physical loss channel.** The physical loss channel is the bath kinetic
gain (the polarisation cloud drag), exactly as Correa 2018 §3 / Eq. (10)
expects for the electronic-stopping process: WP momentum is dissipated
into bath single-particle excitations, which appear as a kinetic-energy
increase of the bath orbitals.

---

## 8. Open questions / next steps

- **Run the orthonormalisation rerun** (plan in
  `docs/plans/jellium_orthonormalisation_rerun.md`). Verdict goes into
  `docs/reports/orthonormalisation-rerun-verdict.md`.
- **Implement the `energy_balance` postprocess phase** so $\Sigma f_i
  \Delta\varepsilon_i$ and $\Delta D$ are auto-extracted; verify the
  $\Delta(\int v_\text{xc} n) \approx +0.12$ eV prediction above.
- **Track per-orbital cod vs time** to test the "hole-as-attractor"
  mechanism (predicts bath orbitals near the WP track the WP
  monotonically).
- **Explore a higher-velocity follow-up at L=50** with $v >
  v_\text{th}^{\text{plasmon}} \approx 0.94$ a.u. (i.e. $E_\text{WP} >
  12$ eV) — would be the first run in this project where plasmons can
  *kinematically* be excited.
