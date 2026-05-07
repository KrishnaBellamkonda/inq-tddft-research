# Source note: WP-orthonormalisation algorithm (professor PDF)

## Source

`ResearchProject/literature/misc/orthonormalization-by-professor.pdf`
(2 pages, dated 2025-06-19, authored by the supervising professor;
no other attribution on file).

## What the document specifies

A single-step Gram–Schmidt projection that produces a wave-packet trial
function $\chi(\mathbf r)$ orthogonal to a given orthonormal set of
Kohn–Sham (KS) orbitals $\{\psi_n^{KS}\}$, **without** modifying the KS
orbitals.

### (i) Test function

$$\phi(\mathbf r) = e^{i \mathbf k_0 \cdot \mathbf r}\, g(\mathbf r),
\qquad
g(\mathbf r) = e^{-|\mathbf r - \mathbf r_0|^2 / (2\sigma^2)}.$$

A Gaussian envelope of width $\sigma$ centred at $\mathbf r_0$, modulated
by a plane wave with wavevector $\mathbf k_0$.

### (ii) Projection of $\phi$ onto $\{\psi_n^{KS}\}$ in momentum space

Plane-wave expansions on the supercell-permitted reciprocal grid $\{\mathbf G\}$:

$$\psi_n^{KS}(\mathbf r) = \frac{1}{\sqrt\Omega}\sum_{\mathbf G} C_n(\mathbf G)\, e^{i\mathbf G \cdot \mathbf r},
\qquad
\phi(\mathbf r) = \frac{1}{\sqrt\Omega}\sum_{\mathbf G} F(\mathbf G)\, e^{i(\mathbf k_0 + \mathbf G)\cdot \mathbf r},$$

with the explicit Gaussian Fourier amplitude

$$F(\mathbf G) = \frac{1}{\sqrt\Omega}\int_\Omega g(\mathbf r)\, e^{-i\mathbf G \cdot \mathbf r}\, d^3 r.$$

The closed-form projection coefficient is then

$$\boxed{\langle \psi_n^{KS} | \phi \rangle
= \sum_{\mathbf G} C_n^{*}(\mathbf G)\, F(\mathbf G - \mathbf k_0)}.$$

The shift by $\mathbf k_0$ on the Gaussian Fourier transform encodes the
plane-wave envelope analytically.

### (iii) Orthogonalised function

$$|\chi\rangle = |\phi\rangle - \sum_n \langle \psi_n^{KS} | \phi \rangle\, |\psi_n^{KS}\rangle,$$

then **renormalise** $\chi$ before injection.

This is **classical Gram–Schmidt** with a **single projection sweep** —
each overlap $\langle \psi_n^{KS} | \phi \rangle$ is computed against the
*original* $\phi$, not against the running residual.

## What the inq-stack code does (current behaviour)

`inq-stack/include/inqkit/wavepacket/wavepacket.hpp::inject_into_last_extra_state`
(method body around lines 200–270) implements:

1. Loop $i = 0, \ldots, n_\text{occ} - 1$.
2. For each $i$, compute the overlap by **real-space** GPU reduction:
   $$\langle \psi_i | \psi_\text{wp} \rangle = \int \psi_i^{*}(\mathbf r)\, \psi_\text{wp}(\mathbf r)\, d^3 r,$$
   discretised on the INQ real-space grid (real and imaginary parts as
   two separate reductions, lines 207–228).
3. Subtract the projection in place:
   $\psi_\text{wp} \leftarrow \psi_\text{wp} - \langle \psi_i | \psi_\text{wp} \rangle\, \psi_i$
   (lines 232–246).
4. After the loop, renormalise $\psi_\text{wp}$ (lines 252–271).

## Concrete differences

| Aspect | Professor PDF | inq-stack code |
|---|---|---|
| Domain of inner product | momentum space (sum over $\mathbf G$) | real space (integral on grid) |
| Gaussian Fourier amplitude | analytic $F(\mathbf G)$ closed-form | implicit (lives in the gridded $\psi_\text{wp}$) |
| Gram–Schmidt variant | **classical** (single sweep against fixed $\phi$) | **modified** (overlaps computed against the running residual $\psi_\text{wp}$, which is updated in place between iterations) |
| Re-orthogonalisation | not specified (single pass) | not implemented (single pass) |
| Normalisation step | "needs normalizing" (explicit) | implemented (final renormalisation) |

## Are the two algorithms equivalent?

For a perfectly orthonormal $\{\psi_n^{KS}\}$, classical and modified
Gram–Schmidt produce **algebraically identical** results: every $\psi_i$
is orthogonal to every other, so subtracting the overlap with $\psi_i$
does not affect future overlaps with $\psi_j$ ($j \ne i$). The
projection is well-defined and order-independent.

Numerically the two differ in floating-point rounding: **modified GS is
generally more stable** when the basis is not exactly orthonormal (e.g.
when the KS orbitals carry small numerical orthogonality defects from
the SCF). Because of this we should not consider the inq-stack code
"wrong" relative to the professor's PDF — it is a numerically robust
re-implementation of the same projection in a different basis.

The choice of *real-space integration on the INQ grid* vs the
*plane-wave momentum-space sum* produces the same value up to grid
quadrature error, which on jellium's uniform grid is identical to
truncation of the plane-wave sum at the same cutoff (both are
finite-Fourier representations of the same operator).

## Inferences (clearly labelled)

- **Inference**: the professor's exposition is best read as a
  *derivation* (how the analytic projection is obtained from the test
  function), not as a *prescription* for how to implement it. The
  algorithm written down is the classical projector $1 - \sum_n
  |\psi_n\rangle\langle\psi_n|$, applied once to $\phi$, then
  renormalised. Our code computes the same projector, in real space,
  using the modified Gram–Schmidt order. No discrepancy of substance.
- **Inference**: a momentum-space implementation could be slightly more
  accurate at large $|\mathbf k_0|$ where the real-space Gaussian is
  rapidly oscillating; in our regime ($|\mathbf k_0| \sim 3.83$ Bohr$^{-1}$
  on a $\Delta x = 0.5$ Bohr grid: $k_0 \Delta x \sim 1.92$ rad, well
  below the Nyquist $\pi$) we expect grid quadrature to be acceptable.
  This is a candidate for a future cross-check rather than an immediate
  fix.

## Recommendation

No code change is currently required. If we later see large
post-injection overlaps (`max_overlap` field of `injection_report`)
that exceed a few × $10^{-5}$, two improvements are available:

1. **Re-orthogonalise twice** (CGS2 / "twice is enough"): repeat the
   projection sweep a second time after the first renormalisation. This
   eliminates almost all residual non-orthogonality without restructuring
   the algorithm.
2. **Switch to momentum-space projection** as in the PDF, exploiting the
   analytic $F(\mathbf G)$. This is more invasive but tracks the
   reference algorithm exactly.

Both should be deferred until the base run reports actual injection
defects; with $\sigma = 1.00$ Bohr and the closed-shell N=128 KS basis
we expect `max_overlap` $\lesssim 10^{-4}$ from the current code.
