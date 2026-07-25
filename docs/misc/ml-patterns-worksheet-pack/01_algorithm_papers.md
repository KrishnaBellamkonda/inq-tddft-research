# 01 — Algorithm papers (seminal references for the methods used)

Scope: the **data-driven / numerical algorithms actually applied** in the
`ml-patterns` campaign. Physics-context papers (Lindhard linear response, wake
theory, ML-stopping surrogates) are deliberately **out of scope** here — that was
the user's explicit choice ("algorithm papers only"). They are catalogued in the
campaign research doc
`docs/campaigns/ml-patterns/research/ml_induced_density_research.md` if needed.

Every citation below is transcribed from the campaign kernel docstrings (which
carry them inline) — none are invented. Each entry lists: the seminal reference,
what the method is, the campaign kernel that implements it, and which task used it.

---

## 1. POD / PCA / Karhunen–Loève (via truncated SVD)

- **Lumley, J. L. (1967)** — "The structure of inhomogeneous turbulent flows,"
  in *Atmospheric Turbulence and Radio Wave Propagation* — origin of Proper
  Orthogonal Decomposition.
- **Eckart, C. & Young, G. (1936)**, *Psychometrika* 1, 211 — low-rank optimality
  of the truncated SVD (the rank-k SVD minimises the Frobenius reconstruction
  error; this is *why* POD modes are the optimal energy-ranked basis).
- **Halko, N., Martinsson, P.-G. & Tropp, J. (2011)**, *SIAM Review* 53, 217 —
  randomized SVD, used for the memory-bounded truncation on large voxel fields.
- Textbook/review: **Brunton, S. & Kutz, J. N.**, *Data-Driven Science and
  Engineering*, Ch. 1; review "Modal Analysis of Fluid Flows" (arXiv:2111.04829).

What it does: given a snapshot matrix `X` (rows = spatial voxels, columns =
time snapshots), the economy SVD `X = U S Vᵀ` yields orthonormal **spatial modes**
(columns of `U`), a **mode-energy spectrum** `E_i = S_i²` and its fraction
`f_i = S_i²/Σ S_j²`, and **temporal coefficients** `a = S Vᵀ`. The number of modes
needed to reach 90 % cumulative energy ("POD rank") and the leading-mode energy
fraction are the descriptors reported.

- **Kernel:** `docs/campaigns/ml-patterns/kernels/pod.py`
- **Used in:** T1 (PCA on the matched difference field Δn); the **POD/DMD
  bath-structure sweep** (`pod_rank90`, `lead_energy` per run).

---

## 2. Dynamic Mode Decomposition (exact DMD) / Koopman

- **Schmid, P. J. (2010)**, *J. Fluid Mech.* 656, 5 — Dynamic Mode Decomposition.
- **Tu, J., Rowley, C., Luchtenburg, D., Brunton, S. & Kutz, J. N. (2014)**,
  *J. Comput. Dynamics* 1, 391 — "On Dynamic Mode Decomposition" (**exact DMD**,
  the variant implemented).
- Textbook: **Brunton & Kutz**, *Data-Driven Science and Engineering*, Ch. 7.
  (Koopman-operator theory, Koopman 1931, is the underlying justification for DMD
  as a linear representation of nonlinear dynamics.)

What it does: from time-shifted snapshot pairs `X, X'` it fits the best-fit linear
operator `X' ≈ A X`, reduces via SVD to `Ã`, and eigendecomposes. Each **DMD mode**
is tagged with a complex eigenvalue → a **continuous frequency** `ω_i =
ln(λ_i)/dt` giving an oscillation frequency (reported in eV) and a growth/decay
rate. Applied over a chosen snapshot window (the near-constant-velocity early
stretch) because the decelerating light projectile makes the dynamics
non-stationary.

- **Kernel:** `docs/campaigns/ml-patterns/kernels/dmd.py`
- **Used in:** the **POD/DMD bath-structure sweep** (`dmd_omega_ev`,
  `dmd_growth` per run); T3 wake analysis.

---

## 3. SINDy — Sparse Identification of Nonlinear Dynamics

- **Brunton, S., Proctor, J. & Kutz, J. N. (2016)**, *PNAS* 113, 3932 — SINDy:
  sparse regression of a time-derivative onto a candidate library of nonlinear
  terms, assuming the governing dynamics is sparse in that basis.

What it does: builds a library `Θ` of candidate terms (monomials, products) and
solves `d/dt(state) = Θ c` with a sparsity-promoting regression so that only a few
`c` survive — yielding a compact ODE. In the campaign it was applied in a **low-dim
latent coordinate** (POD coefficients) as a reduced-order-model cross-check, and it
is the conceptual base of the PDE-FIND kernel (entry 4).

- **Used in:** T5 (2-mode latent ODE); base idea behind Track-B PDE discovery.

---

## 4. PDE-FIND / STRidge — data-driven discovery of a governing PDE

- **Rudy, S., Brunton, S., Proctor, J. & Kutz, J. N. (2017)**, *Science Advances*
  3, e1602614 — "Data-driven discovery of partial differential equations." Builds
  a candidate library `Θ(u)` of spatial operators × monomials, regresses the
  time-derivative `b = ∂ₜᵐ u` onto it, and enforces sparsity by **STRidge**
  (sequentially thresholded ridge regression / STLSQ).

What it does here: a **broad, agnostic** library (powers of `u` up to `poly` ×
spatial derivatives up to `deriv_order`, with cross-products), target order
`m ∈ {1, 2}` (m = 2 for intrinsically second-order plasma oscillation), and two
in-kernel validation "walls" — **forward-integration prediction** on a held-out
temporal window and **bootstrap coefficient stability**. Physical names of
surviving terms are assigned **post-hoc**, never seeded (ADR 0012).

- **Kernel:** `docs/campaigns/ml-patterns/kernels/pdefind.py`
- **Used in:** T11 (discover `PDE_classical`), T12 (discover `PDE_WP`),
  and the synthetic known-PDE **recovery demo** notebook.

---

## 5. Form-factor / linear-response residual method

The two form-factor kernels are built on a single standard analytic identity plus
a campaign-designed test procedure (there is no external "seminal paper" for the
residual test itself — it was designed by the scientific panel, 2026-07-06).

- **Analytic form factor:** the Fourier transform of a normalised 3-D Gaussian
  charge cloud of std `σ` is `F(q) = exp(−q²σ²/2)` — a standard result (e.g.
  **Jackson, *Classical Electrodynamics*, Ch. 3**, form factor of a charge
  distribution).
- **Method (`formfactor_residual.py`):** in linear response the induced density
  factorises as `n_ind(q,ω) = χ(q,ω)·V_ext(q)` with `χ` a property of the medium
  (identical for both projectiles); the WP is then a **low-pass-filtered** point
  charge, `n_WP(q,t) = F(q)·n_cl(q,t)`. Taking the **time-domain ratio**
  `|R(q,t)| = |n_WP(q,t)| / |n_cl(q,t)|` cancels `χ` without any frequency
  binning, and fitting `|R(q)| ~ exp(−a q²)` recovers the filter width.

- **Kernels:** `docs/campaigns/ml-patterns/kernels/formfactor.py`,
  `docs/campaigns/ml-patterns/kernels/formfactor_residual.py`
- **Used in:** T2 (q-space form-factor ratio); the **linear-response residual
  test** (2026-07-06).

> Note: the *physics* null this test challenges (linear-response χ, the Lindhard
> function) is grounded in physics papers held out of this list per the
> "algorithm papers only" scope. See the research doc if the worksheet needs them.

---

## Supporting numerical results (used implicitly)

- **Eckart–Young–Mirsky theorem** — optimal low-rank approximation (underpins POD
  truncation).
- **Randomized SVD** — Halko, Martinsson & Tropp (2011), as above.
- **STRidge / STLSQ** thresholded-regression sparsity (Rudy 2017; Brunton 2016).

## One-line method → task map

| Algorithm | Seminal paper | Kernel | Tasks |
|---|---|---|---|
| POD/PCA (SVD) | Lumley 1967 | `pod.py` | T1, POD/DMD sweep |
| DMD (exact) | Schmid 2010; Tu 2014 | `dmd.py` | T3, POD/DMD sweep |
| SINDy | Brunton–Proctor–Kutz 2016 | (in `discovery`/T5) | T5 latent ODE |
| PDE-FIND/STRidge | Rudy et al. 2017 | `pdefind.py` | T11, T12, recovery demo |
| Form-factor / residual | Jackson Ch. 3 (+ panel design) | `formfactor*.py` | T2, linear-response residual |
