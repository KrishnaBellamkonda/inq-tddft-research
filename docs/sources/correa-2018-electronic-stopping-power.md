# Source: Correa 2018 — Calculating electronic stopping power in materials from first principles

## Full citation

Alfredo A. Correa, *Calculating electronic stopping power in materials
from first principles*, **Computational Materials Science** **150**,
291–303 (2018). DOI: 10.1016/j.commatsci.2018.03.064.
LLNL preprint, special issue on Radiation Effects.

PDF on disk: `ResearchProject/literature/resources/calculating-electronic-stopping-power-in-materials-from-first-principles.pdf`.

## Relevance to this project

Authoritative review of the TDDFT methodology used for electronic
stopping-power calculations in jellium and beyond. The author is the
co-developer of the QBall first-principles code we use as the
reference for the QuantumKickExtension benchmark. Sets the language
(transient vs. steady-state, channelling, off-channelling, plasmon vs
single-particle excitations, supercell long-wavelength size effect)
that all our jellium WP scattering work inherits.

## Key claims used

### Equations (with page/section)

- **Eq. (1)** (p. 292) — **Bethe formula** (high velocity):
  $$S(v) = \frac{4\pi n}{m_e v^2}\, Z^2 k_e^2 e^4 \log\!\Big(\frac{2 m_e v^2}{I}\Big),$$
  with $I$ the mean excitation energy. Diverges as $v \to 0$, breaks
  down for $v \lesssim \sqrt{I/2m_e}$.
- **Eq. (2)** (p. 292) — **Fermi–Teller formula** (low velocity, linear
  in $v$):
  $$S(v) = Z^2\, \frac{2v}{3\pi}\, m_e\, \frac{k_e^2 e^4}{\hbar^3}\, \log(\pi \hbar v_F\,k_e\,e^2).$$
  First formula to take into account degeneracy of the electrons.
- **Eq. (3)** (p. 292) — **Lindhard stopping power** (homogeneous
  electron gas, RPA dielectric function $\varepsilon(k,\omega)$):
  $$S(v) = \frac{2 Z^2 e^2}{\pi v^2}\int_0^{\infty}\!\!\frac{dk}{k}\!\int_0^{kv}\!\!\omega\,d\omega\,\Im\!\Big(\frac{-1}{\varepsilon(k,\omega)}\Big).$$
  Asymptotically equal to Bethe at high $v$ if $I = \hbar\omega_p$. The
  integration domain is the **triangular double integral** in
  $(\omega, k)$-space depicted in Fig. 2.
- **Eq. (4)** (p. 294) — TDDFT Kohn–Sham equation. Confirms the
  one-electron orbitals $\psi_i$ form an "effective orthonormal set"
  with **time-independent occupations** $f \int \psi^*_i \psi_j d^3r =
  \delta_{ij}$ (no $f$ prefix in the orthonormality), i.e. the
  occupations decorate the determinant statically.
- **Eq. (8)** (p. 294) — **TDDFT total-energy decomposition**:
  $$E(t) = \sum_i \int dr\,\psi_i^* \Big(-\frac{\hbar^2 \nabla^2}{2 m_e} + V_\text{ext}\Big) \psi_i \;+\; E_{HXC}[n] \;+\; V_\text{ion-ion} + \sum_J \tfrac{1}{2} M_J \dot{R}_J^2.$$
  Closed-system $E$ is conserved; the $E_{HXC}$ piece breaks into
  Hartree + xc as in the standard ALDA bookkeeping. **Used directly
  in our journal entries' energy bookkeeping.**
- **Eq. (10)** (p. 295) — Stopping power from energy slope:
  $S = \langle dE(t)/dt\rangle / v_\text{proj}$. **Used as the model
  for our cod_z-slope-vs-component-energy diagnostic.**

### Figures used as physical anchors

- **Fig. 1** (p. 292) — qualitative S(v) curve, electronic vs nuclear
  stopping. Establishes that electronic stopping dominates at $v
  \gtrsim$ thermal-ion-velocity scale and peaks near $v \approx v_F$.
- **Fig. 2** (p. 293) — integration domain of Eq. (3) on
  $(k, \omega)$-axes. The **plasmon dispersion line** (red curve) and
  the **electron-hole continuum** (blue/colored region) are shown
  superimposed; the projectile velocity $v$ defines the triangular
  integration cone $\omega < kv$.
- **Fig. 3 + Fig. 4** (p. 293) — for an insulator with electronic gap
  $E_g$ between valence and conduction, the energy-loss function is
  zero for $v < E_g/(2\hbar k_F)$, defining a **threshold velocity**
  $v_\text{th} = E_g / (2 \hbar k_F)$.
- **Fig. 8** (p. 297) — **transient vs. steady-state** distinction.
  "This 'kick' or shake creates a transient in the electronic system…
  this transient state of the system eventually disappears as the
  system enters into a steady state… the process involved in the
  transient is not representative of the asymptotic electronic
  stopping power process, and these points should be discarded from
  the analysis." **Direct justification for our rule that spectra
  must be computed with the transient excluded
  (`docs/observables_reference.md §13.6`).**
- **Section 6.3** (p. 299, "Plasmons and size effects") — periodic
  supercell of size $L$ enforces a longest allowed wavelength; charge
  oscillations at velocities below $\hbar \omega_p L / 2\pi$ would
  decay faster than in an infinite system, breaking the Bethe limit
  in size-limited supercells. **Direct constraint on our L=30 vs
  L=50 vs L=60 comparisons.**
- **Eq. (15)** (p. 300) — long-wavelength size correction
  $\Delta S(v) = (2 Z^2 e^2/\pi v^2)\int_0^{2\pi/T}\!dk/k \int_0^{kv}\!\omega\,d\omega\,\Im(-1/\varepsilon(k,\omega))$
  — the supercell "missing piece" we need to add at low $v$.

## Limitations / uncertainties

- The review is for **uniform / quasi-uniform electron gas + ion
  channelling** geometry; our WP-jellium scattering is geometrically
  closer to the channelling case (the WP plays the role of the
  projectile) but with the projectile being a *delocalised electron
  packet* rather than a point ion. The Coulomb-singularity-induced
  effects (Bethe logarithm, core electrons) are absent in our setup.
- Bethe formula numerator $4\pi n$ vs. $4\pi n_\text{val}$ ambiguity
  (Sec. 6.2, p. 299): for jellium, $n$ is unambiguous (homogeneous
  density); for our WP-on-jellium, the relevant $n$ is the bath
  density.
- The plasmon energy referenced (Li bulk $\hbar\omega_p \approx 6.56$
  eV from DFT-RPA) is **species-specific**; for our jellium at
  $r_s \approx 5.97$ Bohr (N=138 at L=50) the corresponding $\hbar\omega_p$
  is much smaller — see the topical journal entry
  `plasmons-and-stopping-power.md` for the substitution.
- Our jellium runs use **ALDA** xc (adiabatic LDA); Eq. (4) framework
  in this paper allows for any TDDFT approximation — the
  physical conclusions about transient/steady-state apply
  irrespective of the xc kernel.

## Cross-references

- Plans:
  - `docs/plans/jellium_orthonormalisation_rerun.md`
- Handovers:
  - `docs/handovers/jellium_l50_n162_observables.md`
- Journal entries:
  - `docs/journals/researchproject/2026-05-05_run_base_n138_L50_E1p5.md`
  - `docs/journals/researchproject/2026-05-05_run_base_n162_L50_E1p5.md`
  - `docs/journals/researchproject/plasmons-and-stopping-power.md` (topical)
- Reports:
  - `docs/reports/qball-spectra-comparison.md`
- Code (logic adapted): see comments in
  `inq-stack/python/inqview/postprocess/observables.py` (transient
  exclusion, once implemented per `observables_reference.md §13.6`).
