# Interpretable ML on induced electron-density fields — research scoping for the classical-vs-quantum projectile campaign

Deep-research synthesis. Author: research agent (Opus). Date: 2026-06-30.
Scope: scoping the `ml-patterns` campaign — using interpretable ML as a
discovery/representation tool on the induced electron-density fields n(r,t)
produced by rt-TDDFT runs of a light projectile (classical Gaussian-charge
electron vs quantum wavepacket) in jellium, to find and explain spatial &
dynamical differences "the scalar S(v) cannot see."

Method note: every substantive external claim below carries an author-year +
venue/arXiv tag. Papers I fetched and read in full are marked **[fetched]**;
claims taken only from search-result abstracts are marked **[abstract-only]** and
should be upgraded by reading the source before they enter a manuscript. My own
reasoning is tagged **Inference:**. Trust ordering follows the project rule
(peer-reviewed > textbook > official docs > lecture notes; forums/blogs excluded).

---

## Executive summary — the 5 most important findings

1. **The exact campaign (ML to EXPLAIN the induced-density structure of a quantum
   vs classical projectile) has not been done — but every adjacent half exists.**
   ML *surrogates* for the **scalar** stopping force from TDDFT exist and are mature
   (Ward et al., npj Comput. Mater. 2024 — 10⁷× speedup) **[fetched]**, and the
   quantum-vs-classical-projectile *physics* is being actively worked out in
   rt-TDDFT (Kononov et al., arXiv:2511.00759, 2025; Nazarov & Gross,
   arXiv:2510.26222, 2025) **[fetched / project note]**. But these reduce the
   density to a number (S, or an effective charge). **No located work uses
   interpretable ML on the 3D(+t) induced-density field itself to characterise the
   classical↔quantum difference.** That is the campaign's genuine gap.

2. **The single strongest piece of precedent — and the template to beat — is Ward,
   Blaiszik, Lee, Martin, Foster & Schleife, npj Comput. Mater. 10, 155 (2024)**
   (arXiv:2311.00787) **[fetched]**: they feed the *ground-state density* +
   projectile position/velocity into an ML model and predict the scalar stopping
   force, 10⁷× faster than TDDFT. It proves density-field features carry the
   stopping physics, but it throws away exactly the object our campaign studies (the
   *induced, time-dependent* field). Our campaign is the interpretive inverse of
   theirs: keep the field, discover the structure.

3. **Method-per-rung recommendation (see Q2):**
   - **Rung 1 (static spatial structure):** **POD/PCA on the matched
     difference-field Δn = n_quantum − n_classical**, plus **persistent-homology /
     topological descriptors** for fringe/hole counting. PCA is the honest, fully
     interpretable baseline; it is the linear backbone of every method below and
     must be run first.
   - **Rung 2 (spatiotemporal dynamics):** **DMD / Koopman** is the best-fit primary
     tool — it returns spatial modes each tagged with a complex frequency ω + decay
     rate, which maps directly onto plasmon/wake physics (λ = 2πv/ω_p). DMD is
     already proven on first-principles real-time electronic dynamics (Yin et al.,
     npj Comput. Mater. 2024, electron–phonon Boltzmann) **[abstract-only]**.
     **SINDy** is a high-payoff secondary for discovering a reduced governing
     equation in the POD/DMD latent space.

4. **Of the four target signatures, two are physically robust and likely
   observable; two are at serious risk of being dominated by artifacts.** Robust:
   (iii) the **collective wake** (phase/wavelength/decay differences are first-order
   and large) and (iv) **projectile form-factor softening** exp(−q²σ²/2) (a clean,
   analytic, q-space signature). At risk: (i) the **exchange/Fermi-hole** difference
   and (ii) **quantum diffraction fringes** — both are subtle and live on the same
   length scale as, and can be mimicked by, the **WP self-interaction error (~7 eV)
   and dispersion**. These two MUST be analysed only on the SIE-subtracted field
   (vacuum-WP control), or the "discovery" is an artifact (Q3, Q5).

5. **The dominant methodological risk is spurious discovery from un-normalised
   fields and a coordinate-convention trap.** Any ML run on raw n(r,t) will "find"
   (a) the ground state, (b) the trivial rigid translation of the projectile, and
   (c) the SIE/dispersion of the WP — none of which is the physics of interest.
   The required hygiene is a **subtraction ladder**: subtract GS → subtract the
   common rigid-projectile motion → subtract linear-response (Lindhard) → subtract
   vacuum-WP (SIE). And the project's standing rule holds with full force: **these
   VTIs are in PHYSICAL order — never `np.fft.fftshift` a loaded field** (a
   centre↔edge swap produces plausible-but-wrong "modes").

**Bottom line on the premise:** the campaign is *not* already-solved and *not*
weak — but its headline claim ("differences the scalar S cannot see") is only
defensible for signatures (iii) and (iv) unless the SIE control is done rigorously;
signatures (i) and (ii) are where the campaign could most easily fool itself.

---

## Q1 — Precedent: has this been done?

I searched four sub-areas. Summary verdict: **the pieces exist; the synthesis does
not.**

### (a) ML / reduced-order methods on rt-TDDFT or TDDFT electron-density dynamics
- **Ward, Blaiszik, Lee, Martin, Foster, Schleife, npj Comput. Mater. 10, 155
  (2024)** (arXiv:2311.00787) **[fetched]** — ML surrogate for electronic stopping.
  Inputs = projectile position + velocity + **ground-state electron density**;
  output = **scalar** stopping force; 10⁷× fewer core-hours than TDDFT; demonstrated
  on proton-in-Al, predicts Bragg-peak depth vs incidence angle. **Does not analyse
  the induced/time-dependent density** — it is a forward surrogate, not a discovery
  tool. *This is the closest and strongest precedent.*
- **Chiang, Choi, Osei-Kuffuor, arXiv:2509.00169 (2025)** **[fetched]** —
  "Generative Latent Space Dynamics of Electron Density": 3D convolutional
  autoencoder + latent diffusion model to roll out electron-density *trajectories*.
  **Important caveat:** trained on **ground-state AIMD** densities (liquid Li, 800 K),
  **not rt-TDDFT** induced response. Proves the AE-on-3D-density-trajectory machinery
  works and is stable over long rollouts; says nothing about interpretability of the
  latent space. Directly relevant as a Rung-2 *architecture* precedent.
- **ML time-propagators for the TDDFT density** — Gupta et al. / "Accelerating
  Electron Dynamics Simulations through Machine Learned Time Propagators"
  (arXiv:2407.09628; also IOP *Mach. Learn.: Sci. Technol.*, 2025) and the
  moment-propagation-theory ML work (arXiv:2412.05260, *JCTC* 2025)
  **[abstract-only]** — autoregressive neural operators that propagate n(r,t).
  These are *acceleration* tools, not interpretive; but they confirm the density
  field is a learnable dynamical object.
- **DMD on first-principles real-time dynamics** — Yin et al., "Dynamic mode
  decomposition of nonequilibrium electron–phonon dynamics," npj Comput. Mater.
  (2024, s41524-024-01308-4) **[abstract-only]**: DMD on ~10% of a real-time
  Boltzmann trajectory extrapolates to steady state and **identifies dominant
  electronic processes**. Plus DMD/Koopman for density-matrix dynamics
  extrapolation (e.g. arXiv:2203.14892 on non-equilibrium Green's functions)
  **[abstract-only]**. These establish DMD as a *legitimate, interpretable* tool for
  electronic real-time data — the cornerstone of the Rung-2 recommendation.

### (b) ML for electronic stopping power (surrogates, feature discovery)
- Ward 2024 (above) — interpolation surrogate.
- Stacking-ensemble ML for stopping powers (arXiv:2208.00227) **[abstract-only]** —
  empirical/tabular regression (SRIM-like data), not field-based.
- **Gap:** all stopping-power ML I found predicts the *number*. None does *feature
  discovery on the induced field* to explain *why* S differs.

### (c) Projectile/ion WAKES in an electron gas / plasma — data-driven analysis
- Wake *theory* is mature and is the physics anchor, not an ML precedent:
  Echenique, Ritchie & Brandt, Phys. Rev. B 20, 2567 (1979); Echenique & Ritchie,
  Phys. Rev. B 27, 4117 (1983) (RPA wake) — see project note
  `docs/sources/echenique-ritchie-image-wake-stopping.md`.
- **DMD/POD are heavily used on plasma and fluid WAKES** (e.g. ExB-plasma DMD,
  Kawashima/Boeuf-style; PIC space-charge ROM via DMD, arXiv:2303.16286; magnetised-
  plasma DMD, *Phys. Plasmas* 27, 032108 (2020)) **[abstract-only]**. This is strong
  *methodological* precedent that modal decomposition extracts wake structures with
  physical frequencies — but applied to classical plasma PIC/fluid data, never (that
  I found) to a quantum induced-density wake.
- **Inference:** the ExB/PIC-wake DMD literature is the best *transfer* evidence that
  DMD will cleanly pull the plasmon-wake mode (ω ≈ ω_p, wavelength 2πv/ω_p) out of
  our n(r,t). It is a method import, not a duplicated result.

### (d) Quantum-vs-classical projectile response / stopping
- **Kononov, Hentschel, Hansen, Baczewski, arXiv:2511.00759 (2025)** **[fetched]** —
  nonlinear (charge-sign / Z) effects in light-ion stopping in rt-TDDFT for warm
  dense matter; nonlinear processes shift S by ~10% near/below the Bragg peak;
  uses *fractionally-charged* projectiles to isolate the linear-response limit; and
  states plainly that **"induced density analysis has been previously used as a tool
  to understand the molecular-level details of the stopping process."** → precedent
  that induced-density *inspection* is standard practice, but it is done **by eye /
  by hand, not by ML**.
- **Nazarov & Gross, arXiv:2510.26222 (2025)** (project note
  `docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`) — quantum
  (wavepacket) vs classical point-charge stopping in an electron liquid; friction is
  mass/width dependent (a purely quantum effect); classical M→∞ limit recovers
  Lindhard. *This is the theoretical backbone of the campaign's quantum-vs-classical
  premise* — but it is analytic/response-theory, not ML, and not a real-space
  field-structure study.

### The GAP this campaign fills (one sentence)
> Existing work either (i) compresses the TDDFT density to a scalar via ML, or
> (ii) inspects the induced density by hand to explain quantum-vs-classical stopping
> physics — **no one has used interpretable ML (POD/DMD/SINDy + contrastive
> representation) on the full 3D(+t) induced-density field to systematically
> discover and explain the spatial and dynamical differences between a quantum
> wavepacket and a matched classical projectile.** The campaign occupies the empty
> intersection of (a)+(d).

---

## Q2 — Interpretable ML families for 3D(+t) density fields

For each: assumption / data need / interpretability / failure mode / rung.

### POD / PCA (a.k.a. Karhunen–Loève, EOF)
- **What it is / assumes:** SVD of the snapshot matrix; assumes the dominant
  variance lives in a low-dim *linear* subspace; modes are orthogonal and ranked by
  captured variance ("energy") (Lumley 1967; review arXiv:2111.04829)
  **[abstract-only]**.
- **Data need:** modest — works on tens–hundreds of snapshots; the field is
  flattened to a vector per time/per run. Resolution only needs to resolve the
  feature (here ≲ k_F⁻¹ and ≲ σ).
- **Interpretability:** highest of all methods — each mode is a real spatial field;
  the spectrum tells you the effective dimensionality. For *static* comparison, run
  PCA on the **stacked set of matched difference-fields** Δn = n_q − n_cl across
  (v, σ) to see what few spatial patterns the quantum-vs-classical difference lives
  in.
- **Failure modes:** purely linear/variance-optimal → smears travelling waves across
  many modes (a moving wake is NOT compactly represented by POD — this is the
  classic POD weakness that motivates DMD); sensitive to un-subtracted offsets (will
  put the ground state / mean in mode 1); mode mixing when two phenomena have
  comparable energy.
- **Rung:** **Rung 1 (primary baseline).** Also the pre-compression front-end for
  DMD/SINDy in Rung 2.

### DMD / Koopman
- **What it is / assumes:** fits a best-fit linear operator A advancing snapshots
  x_{k+1} ≈ A x_k; eigen-decomposition gives **spatial modes each with a single
  complex eigenvalue → oscillation frequency ω and growth/decay rate** (Schmid 2010;
  Tu et al. 2014; review *Annu. Rev. Fluid Mech.* 2025) **[abstract-only]**.
  Koopman theory justifies it as a linear representation of nonlinear dynamics in an
  observable space.
- **Data need:** uniform time sampling (our 300-frame VTI cadence is well-suited);
  needs the sampling interval dt small enough to resolve ω_p (Nyquist: dt < π/ω_p).
  Works on a single trajectory.
- **Interpretability:** very high *for oscillatory/wave physics* — a mode at
  Re(ω)≈ω_p with the spatial wavelength 2πv/ω_p IS the plasmon wake; the decay rate
  is the wake damping. Directly testable against analytic wake theory.
- **Failure modes:** assumes (approximately) linear-time-invariant dynamics — a
  *decelerating* light projectile (see project rule `light-projectile-stopping.md`)
  is NON-stationary, so standard DMD over the whole run will blur frequencies;
  mitigate with **windowed / multi-resolution DMD** over the near-constant-velocity
  early window, or **DMD with control**. Sensitive to noise (use exact/optimised DMD,
  TLS-DMD, or rank-truncated DMD); spurious modes if rank over-chosen.
- **Rung:** **Rung 2 (primary).**

### SINDy (sparse governing-equation discovery)
- **What it is / assumes:** sparse regression of time-derivatives onto a library of
  candidate terms; assumes the dynamics is *sparse* in a chosen basis (Brunton,
  Proctor, Kutz, PNAS 2016) **[abstract-only]**.
- **Data need:** clean time-derivatives (noise-amplifying) → needs either dense
  sampling or weak-form/integral SINDy; best applied in a **low-dim latent
  coordinate** (POD or autoencoder coefficients), not on the raw 10⁶-voxel field.
- **Interpretability:** the *highest* in principle — outputs an actual equation
  (e.g. ȧ = −γa + ω b for a wake amplitude). Strong adjacency to physics: can test
  whether the WP-coefficient dynamics carry an extra (SIE/dispersion) term absent in
  the classical case.
- **Failure modes:** library must contain the true terms; fails on strongly
  non-sparse or chaotic dynamics; derivative noise; coordinate-dependent (garbage
  latent → garbage equation). PDE-find variants (on the field) need very clean data.
- **Rung:** **Rung 2 (secondary, high-payoff).**

### Autoencoders / VAEs (incl. β-VAE, latent diffusion)
- **What it is / assumes:** nonlinear compression to a latent code; β-VAE encourages
  *disentangled* latent axes (Higgins et al. 2017; reviews 2025) **[abstract-only]**;
  latent-diffusion rollout proven on 3D density trajectories (Chiang 2025)
  **[fetched]**.
- **Data need:** **large** — hundreds–thousands of fields to avoid overfitting;
  3D conv AEs are data- and compute-hungry. Our run count (tens of matched pairs ×
  300 frames) is borderline; augmentation (translations, the matched pairs
  themselves) helps.
- **Interpretability:** **low by default** — latent axes are not physical unless
  forced (β-VAE disentanglement, or contrastive/physics-informed priors). Use only
  *after* PCA/DMD as a nonlinear stress-test, with disentanglement regularisation and
  post-hoc probing.
- **Failure modes:** posterior collapse, entangled latents, hallucinated structure
  (generative models can invent plausible-but-wrong density), drift in rollout
  (mitigated by JS-divergence regularisation per Chiang 2025). High spurious-
  discovery risk → must be cross-checked against the linear methods.
- **Rung:** Rung 1 (disentangled VAE on static fields) and Rung 2 (latent dynamics),
  both **tertiary** — confirmatory, never the sole evidence.

### Contrastive learning on matched pairs
- **What it is / assumes:** learn a representation where matched/similar inputs are
  close and dissimilar far; can "inject physics knowledge … to represent the same
  phenomenon similarly despite outward differences" (HEP applications; review
  *RASTI* 2023) **[abstract-only]**.
- **Fit to this campaign:** *natural* — the campaign IS a matched-pair design
  (quantum vs classical at identical v, σ). A contrastive objective that pulls
  together same-(v,σ) and pushes apart different-(v,σ), then asks *what axis
  separates quantum from classical*, is a principled discovery framing.
- **Data need / interpretability / failure:** needs enough pairs; latent is not
  inherently interpretable (couple with disentanglement / linear probes); risk of
  learning the trivial confounder (SIE, dispersion, mean shift) as the "difference"
  axis unless those are subtracted first.
- **Rung:** Rung 1 (primary *framing*), feeding a probe back to PCA modes.

### Neural ODEs
- **What it is:** learn dx/dt = f_θ(x) as a continuous-time NN; flexible latent
  dynamics. **Interpretability low**, data-hungry, training-unstable on stiff/
  oscillatory systems. **Inference:** weaker fit than DMD/SINDy here because it
  trades the very interpretability the campaign needs for flexibility we don't.
  **Rung:** optional Rung-2 tertiary.

### Topological / structural descriptors (persistent homology)
- **What it is / assumes:** sub-level-set filtration of the scalar field → barcodes/
  persistence diagrams counting connected components, loops, voids; **stable under
  perturbation, machine-readable AND physically interpretable** (review
  arXiv:2411.14390; cosmology arXiv:2412.15405) **[abstract-only]**.
- **Fit:** excellent for *counting and quantifying* the discrete structures the
  campaign predicts — wake oscillation lobes, interference fringes, the depth/extent
  of the exchange hole. A persistence diagram of Δn gives an artifact-robust fringe/
  lobe count that a human can defend.
- **Data need:** low; per-field. **Interpretability:** high for structure-counting,
  low for dynamics. **Failure:** insensitive to *phase*/sign unless signed
  filtration used; threshold choices.
- **Rung:** **Rung 1 (primary descriptor alongside PCA).**

### Recommended mapping
| Rung | Primary | Secondary | Confirmatory |
|---|---|---|---|
| **1 — static spatial** | POD/PCA on Δn; persistent homology | contrastive matched-pair representation | β-VAE (disentangled) |
| **2 — spatiotemporal** | DMD/Koopman (windowed) | SINDy in latent coords | latent-diffusion AE; neural-ODE |

---

## Q3 — Physics grounding of the four target signatures

For each: is the classical↔quantum difference physically expected, and is it
expected to be observable/dominant *vs the SIE artifact*?

### (i) Exchange / Fermi-xc hole (WP indistinguishable from bath vs distinct classical particle)
- **Physics:** the classical Gaussian-charge electron is a distinguishable external
  potential; the bath screens it (Friedel/exchange-correlation screening cloud). The
  quantum WP is *one orbital of the same KS determinant* — it is antisymmetrised with
  the bath, so it carries its own Fermi/exchange hole and is (partly)
  indistinguishable. The xc hole of the homogeneous electron gas is ~few-Å scale and
  strongly depresses pair density at short range (RPA/plasmon-pole; arXiv:2301.05590
  "time-dependent xc hole of the electron gas") **[abstract-only]**.
- **Expected difference:** YES in principle — antisymmetrisation changes the
  short-range induced-density structure around the projectile.
- **Observable / dominant vs SIE?** **HIGH RISK.** This signature lives at exactly
  the length/energy scale of the WP **self-interaction error (~7 eV, the WP feeling
  its own Hartree)** and **dispersion**. In a (semi-)local XC functional the WP's
  "exchange hole with the bath" is itself approximated, and the SIE *is* a
  mistreatment of the projectile's self-exchange. **Inference:** without the
  vacuum-WP subtraction this signature is essentially *indistinguishable from the
  SIE artifact*; it should be reported only on the SIE-subtracted field and even then
  flagged as functional-dependent. This is the campaign's weakest claim.

### (ii) Quantum diffraction / interference fringes
- **Physics:** a quantum WP scattering off the density inhomogeneity can split and
  self-interfere, producing fringes in |ψ|² (and hence in n_wp) with spacing set by
  the de Broglie wavelength / local Fermi wavelength (self-interference in electron
  scattering, arXiv:1710.02583; fringe spacing ≈ ½ λ_F near constrictions,
  arXiv:1206.1371) **[abstract-only]**. The classical projectile has no such
  internal interference.
- **Expected difference:** YES — fringes are a genuinely quantum, classical-absent
  feature, and they are a *single-particle* (WP) effect, so they appear in n_wp.
- **Observable / dominant vs SIE?** **MODERATE RISK.** Fringes are real but *subtle*
  and on the same scale as **WP dispersion ripples and grid/Gibbs artifacts**. The
  discriminator: true diffraction fringes scale with v and σ in a de-Broglie-
  predictable way and survive vacuum-WP subtraction (vacuum-WP has dispersion but no
  bath-scattering interference); spurious ripples do not. Persistent homology (fringe
  count) + the σ/v scaling test is the right adversarial check.

### (iii) Collective WAKE phase / wavelength / decay (λ = 2πv/ω_p; linear vs nonlinear, Barkas Z³)
- **Physics:** a swift charge drives an oscillating induced-density wake trailing it,
  with characteristic wavelength **λ = 2πv/ω_p** (Echenique–Ritchie–Brandt, Phys.
  Rev. B 20, 2567 (1979); RPA wake, Phys. Rev. B 27, 4117 (1983)) — see
  `docs/sources/echenique-ritchie-image-wake-stopping.md`. At leading order the wake
  is **charge-even** (∝ Z²); the **charge-odd Barkas Z³ correction** distinguishes +
  from − projectiles (Lindhard, Nucl. Instrum. Methods 132, 1 (1976), project note
  `lindhard-1976-barkas-effect.md`). rt-TDDFT now resolves these nonlinear shifts
  (~10% near the Bragg peak; Kononov et al. 2511.00759, 2025 **[fetched]**).
- **Expected difference (quantum vs classical):** the *wake wavelength* is set by v
  and ω_p and should match between matched runs; the **differences** are expected in
  (a) the **near-field amplitude/phase** (the WP's distributed, spreading charge
  softens the source — see (iv)), and (b) **nonlinear/Barkas-type terms** that
  depend on the projectile's self-consistent potential rather than a point Coulomb
  (Nazarov & Gross 2025 — the quantum version replaces the Coulomb potential with the
  KS potential of the distributed charge).
- **Observable / dominant vs SIE?** **LOW RISK / most robust.** The wake is a
  first-order, large-amplitude, *bath* feature (lives in n_bath = n_total − n_wp),
  spatially separated from the projectile self-interaction. Subtracting the
  linear-response (Lindhard) wake isolates the nonlinear residual cleanly. **This is
  the signature most likely to yield a defensible ML "discovery."** DMD is purpose-
  built to extract it (mode at ω_p, wavelength 2πv/ω_p, plus its decay rate).

### (iv) Projectile spreading / form-factor softening exp(−q²σ²/2)
- **Physics:** a Gaussian charge of width σ has form factor exp(−q²σ²/2), so it
  couples to the electron gas only up to q ≲ 1/σ — high-q response is suppressed
  relative to a point charge. A quantum WP additionally **spreads in time** (free
  Gaussian disperses: σ(t)² = σ₀² + (ħt/2mσ₀)²), so its *effective* form factor
  softens further as the run proceeds; the classical Gaussian charge keeps fixed σ.
- **Expected difference:** YES and **analytically clean** — the time-growing high-q
  suppression is a sharp, predictable quantum signature. Nazarov & Gross (2025) make
  width/mass dependence the central quantum effect.
- **Observable / dominant vs SIE?** **LOW RISK / clean.** This is a **q-space**
  signature: Fourier-transform n_wp(r,t) and track the q-dependence and its
  *time-broadening* — directly comparable to exp(−q²σ(t)²/2). It is *separable* from
  SIE because SIE shifts energy/phase but the form-factor softening is a width
  measurement. **Caveat:** WP dispersion and the (artifactual) SIE-driven spreading
  both broaden σ(t); the vacuum-WP control measures the *intrinsic* dispersion+SIE
  spreading, and subtracting it isolates the *bath-induced* part.

### Signature scorecard
| Signature | Difference expected? | Robust vs SIE? | Best home (field) | Verdict |
|---|---|---|---|---|
| (i) exchange/Fermi hole | Yes | **No — high risk** | n_total short-range | weakest; SIE-confounded |
| (ii) diffraction fringes | Yes | Moderate | n_wp | promising w/ σ,v scaling test |
| (iii) collective wake | Yes (near-field + Barkas) | **Yes — robust** | n_bath | **strongest** |
| (iv) form-factor softening | Yes | **Yes — clean (q-space)** | n_wp (FT) | **strong, analytic** |

---

## Q4 — Adjacent questions this dataset could answer (ranked by payoff/feasibility)

1. **[High / High] Map the wake's nonlinear (beyond-Lindhard) residual vs v and σ.**
   Subtract the linear-response Lindhard wake from n_bath; DMD/POD the residual. This
   is the cleanest "new physics" extraction and reuses the strongest signature (iii).
   Directly connects to Kononov 2025's ~10% nonlinear shifts but in real space.
2. **[High / High] Velocity-dependence of the wake wavelength as a self-consistency
   check λ(v) = 2πv/ω_p.** DMD wavelength vs v across runs validates the whole
   pipeline against analytic theory before any "discovery" claim — a built-in
   falsification test.
3. **[High / Med] Effective-charge / form-factor extraction from the induced field.**
   Fit the induced density to an effective point-charge + form factor; track
   Z_eff(t) (partial neutralisation à la Kononov 2025) and σ_eff(t) (dispersion).
   Turns the field into two interpretable scalars the scalar S cannot give.
4. **[Med / High] Energy-flux / wake-damping map.** From the time-resolved induced
   density and its currents, localise *where* energy is deposited (near-field vs
   radiated plasmon), complementing the global S with a spatial energy-loss map.
5. **[Med / Med] σ-sweep toward the point-charge (Lindhard) limit.** As σ→0 the
   quantum and classical fields should converge to the point-charge response; the
   *rate* and *manner* of convergence is itself a result (and a consistency check on
   the σ-matching convention).
6. **[Med / Low] Onset of nonlinearity / response saturation.** Vary projectile
   charge (fractional charges, as Kononov 2025 do) to find where linear response
   breaks down in the *field structure*, not just in S.

---

## Q5 — Pitfalls, resolution & normalisation requirements

### The subtraction ladder (mandatory normalisation)
Run every ML method on the *most-subtracted* field appropriate to the signature:
1. **Subtract the ground state:** δn = n(r,t) − n_GS(r). (Removes the dominant,
   physics-irrelevant background — otherwise PCA mode 1 is just n_GS.)
2. **Remove the common rigid-projectile motion:** co-move to the projectile frame or
   subtract the rigid-translation component shared by both runs — otherwise the
   leading "mode" is the trivial translation, identical for quantum and classical.
3. **Subtract linear response (Lindhard):** for the wake (iii), the *nonlinear*
   residual is the interesting object; the Lindhard wake is known analytically.
4. **Subtract the vacuum-WP (SIE + dispersion control):** for any n_wp-based or
   short-range signature (i, ii, iv), subtract the no-bath WP run at matched (v, σ)
   to remove the **~7 eV self-interaction and the intrinsic dispersion**, which have
   no classical counterpart. **This is the single most important control in the
   campaign** — without it, signatures (i) and (ii) are artifacts.

### Method-specific pitfalls
- **POD/PCA:** travelling waves smear across many modes → do not over-interpret mode
  count for the wake; offsets dominate mode 1 → always subtract GS first.
- **DMD/Koopman:** the **decelerating light projectile is non-stationary** (project
  rule `light-projectile-stopping.md`) → standard DMD blurs ω; use **windowed DMD on
  the early near-constant-velocity window** (vz ≥ 0.85 v₀) or DMD-with-control;
  ensure **dt < π/ω_p** (Nyquist) or the plasmon mode aliases.
- **SINDy:** derivative noise → use weak/integral form; only in low-dim latent coords;
  verify the recovered equation predicts held-out time.
- **AE/VAE/diffusion:** generative hallucination + posterior collapse → never the sole
  evidence; cross-check every latent "axis" against a linear PCA mode; report
  reconstruction error.
- **Contrastive:** will latch onto the largest confounder (SIE, dispersion, mean) as
  the "quantum axis" → only run on the subtracted field; validate the discovered axis
  scales with a physical parameter (v or σ).
- **Persistent homology:** sign/phase-blind unless signed filtration; threshold
  choices → report barcodes across thresholds.

### Resolution requirements
- **Spatial:** must resolve the smallest target scale — the exchange hole and
  fringes need grid spacing ≲ k_F⁻¹ and ≲ σ; the wake needs the box long enough to
  hold ≳ 1 wavelength 2πv/ω_p downstream. The existing 4σ/1σ boundary rule and dx
  choices (`shared/configs/boundary_rule.hpp`) already encode this.
- **Temporal:** DMD needs dt < π/ω_p AND enough frames to span ≳ a few plasmon
  periods *within the constant-velocity window*; the 300-frame VTI cadence target is
  appropriate. For light projectiles the usable window is the early drag window, not
  the whole run.

### Coordinate-convention trap (MUST heed)
- **inqkit VTIs are written in PHYSICAL order** (origin −L/2, fft_shift applied at
  write). **Never `np.fft.fftshift` a loaded VTI/field** — it swaps centre↔edge and
  produces plausible-but-wrong "modes" that ML will happily fit (project rule
  `vti-coordinate-mapping.md`). Load every field through `inqview.load_vti`; use its
  returned axes for coordinates. The ONLY field needing fftshift is LEED `.dat`
  (not relevant here). A modal-decomposition pipeline is *exactly* the kind of
  loader that re-derives the wrong "INQ is FFT-natural" assumption — guard against it.
- **Bath density definition:** n_bath = n_total − n_wp is the canonical,
  run-independent bath (project note `reference_canonical_bath_density`); confirm
  which `density_system` convention each run used (WP-included vs bath-only) before
  differencing.

### Data-hygiene checklist (anti-spurious-discovery)
- Always run the **linear baseline (PCA) first**; a nonlinear method must beat it to
  earn its complexity.
- Every "discovery" must **scale predictably with a physical knob** (v, σ, or charge)
  across runs — a feature that appears in one run only is suspect.
- **Validate the pipeline on a known case first:** λ(v) = 2πv/ω_p from DMD must match
  theory before trusting any novel mode (this is the falsification gate, Q4 item 2).
- Keep full-precision data in code; round reported numbers to 2 s.f. (project rule).

---

## Source notes written (verified, fetched papers)
- `docs/sources/ward-2024-ml-stopping-power-surrogate.md` — Ward et al., npj Comput.
  Mater. 2024 (strongest precedent; scalar surrogate).
- `docs/sources/chiang-2025-latent-density-dynamics.md` — Chiang, Choi, Osei-Kuffuor,
  arXiv:2509.00169 (Rung-2 architecture precedent; AE + latent diffusion on density).
- `docs/sources/kononov-2025-nonlinear-stopping-rt-tddft.md` — Kononov et al.,
  arXiv:2511.00759 (quantum/nonlinear stopping physics + induced-density-analysis
  precedent and the gap).

Existing project notes reused: `echenique-ritchie-image-wake-stopping.md`,
`lindhard-1976-barkas-effect.md`, `nazarov-gross-2025-quantum-projectile-stopping.md`.

## Key uncertainties / honesty flags
- Items marked **[abstract-only]** (DMD-on-electron-phonon, ML time-propagators,
  POD/DMD reviews, persistent-homology reviews, contrastive-physics, β-VAE) were read
  only at abstract/search-snippet level. They are *methodological* citations whose
  general claims are well-established, but exact numbers/scope must be verified before
  any go into a manuscript.
- I did NOT independently re-derive the WP SIE = ~7 eV figure or the dispersion
  formula here — they are taken from the campaign brief / project memory and should
  be cited to the project's own validation, not to an external source.
- The strongest premise risk (restated): signatures (i) exchange hole and (ii)
  diffraction fringes are at real risk of being SIE/dispersion artifacts; the
  campaign's defensible core is (iii) wake and (iv) form-factor in q-space.
