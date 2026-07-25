# Localised Jellium — Theory Worksheet

> **Purpose.** Everything you need to *independently validate* a localised-jellium
> implementation in INQ, built up from first principles. Each block ends with a
> **Validation Checkpoint (VC)** — a number or shape you can compute by hand and
> then assert against the code's output. Read it as a scaffold: the ▢ EXERCISE
> lines are deliberately left for you to fill in; the "Answer / check" lines give
> you the target so you know when you're right.
>
> **Scope note.** This worksheet is the *theory*. The design decisions (geometry,
> edge softening, r_s, projectile model, boundary conditions) are still being
> grilled out in the session and will land in `docs/plans/` + `docs/adr/`. Where a
> decision is still open I say so.
>
> **Honesty flags.** Anything marked **⚠ VERIFY** is a number I am *not* certain of
> and that you (or a literature check) should pin down before trusting it. The
> surface-energy target `86.4 erg/cm²` is the main one.

---

## Part 0 — Units, notation, and the one conversion that bites

Work in **Hartree atomic units** (INQ's internal units): ℏ = mₑ = e = 4πε₀ = 1.

| Quantity | a.u. | SI / cgs value of 1 a.u. |
|---|---|---|
| length | Bohr (a₀) | 0.529177 Å = 5.29177×10⁻⁹ cm |
| energy | Hartree (Ha) | 27.21138 eV = 4.35974×10⁻¹¹ erg |
| density | a₀⁻³ | — |
| 1 Ry | ½ Ha | 13.6057 eV |

Symbols used throughout:

- **n₀** — the target uniform interior density of the background, `n₀ = 3/(4π r_s³)`.
- **n₊(r)** — the positive background charge density (the jellium "ions", smeared).
- **n(r)** or **n₋(r)** — the electron number density (always ≥ 0 in DFT convention).
- **N** — number of electrons = ∫ n₊ d³r (neutrality).
- **r_s** — Wigner–Seitz / density parameter (Bohr).
- **R_cl** — the confinement radius (sphere) or half-width (slab/box).
- **k_F, E_F** — Fermi wavevector and energy of the *interior* HEG.
- **v_bg(r)** — the electrostatic potential energy an electron feels from n₊ (the
  confining well). This is the object the new perturbation adds to the KS potential.

**The conversion that bites** (surface energy):
```
σ[erg/cm²] = σ[Ha/Bohr²] × 1.55690×10⁶        (and 1 erg/cm² = 1 mJ/m²)
⇒ 86.4 erg/cm²  =  5.55×10⁻⁵ Ha/Bohr²
```
Jellium surface energies are *tiny* numbers in a.u. — getting this factor wrong by
10× is the classic way to "fail" a benchmark you actually passed.

---

## Part 1 — The homogeneous electron gas (the reference your cluster must approach)

**Learning objective:** know, cold, the four numbers `n₀, k_F, E_F, ε_HEG(r_s)` so you
can predict the interior of any jellium from r_s alone.

### 1.1 Why jellium
Jellium = electrons in a *uniform* positive background. It is the reference system
of solid-state DFT: the LDA exchange–correlation functional *is* the energy per
electron of this system, tabulated vs density. Your localised jellium is a finite
chunk of it, so deep in the interior it must reproduce HEG numbers.

### 1.2 Density parameter
Each electron occupies a sphere of radius r_s:
```
(4/3) π r_s³ = 1/n₀     ⇒     n₀ = 3 / (4π r_s³)
```
▢ **EXERCISE 1.2** — For r_s = 4 a₀, compute n₀ in a₀⁻³.
*Answer / check:* n₀ = 3/(4π·64) = **3.73×10⁻³ a₀⁻³**.

### 1.3 Fermi wavevector and energy
For a spin-degenerate (paramagnetic) gas:
```
k_F = (3π² n₀)^{1/3} = (9π/4)^{1/3} / r_s = 1.91916 / r_s          [a₀⁻¹]
E_F = k_F² / 2                                                      [Ha]
```
▢ **EXERCISE 1.3** — For r_s = 4: k_F and E_F (Ha and eV).
*Answer / check:* k_F = 0.4798 a₀⁻¹; E_F = 0.1151 Ha = **3.13 eV**.

### 1.4 Energy per electron — the HEG limit
Per electron (paramagnetic), in Hartree:
```
kinetic:     t(r_s)  = (3/5) E_F = (3/10) k_F²  = 1.10495 / r_s²
exchange:    ε_x(r_s) = −(3/4)(3/π)^{1/3} n₀^{1/3} = −0.458165 / r_s
correlation: ε_c(r_s) = (Perdew–Zunger 1981 or Perdew–Wang 1992 parametrisation)
ε_HEG(r_s)  = t + ε_x + ε_c
```
(`t + ε_x` are the famous `2.21/r_s² − 0.916/r_s` Ry Gell-Mann–Brueckner terms.)

▢ **EXERCISE 1.4** — For r_s = 4, compute t and ε_x in Ha. Look up ε_c (PZ81 ≈
−0.0450 Ha at r_s=4) and form ε_HEG.
*Answer / check:* t = 0.0691, ε_x = −0.1145, ε_c ≈ −0.0450 ⇒ **ε_HEG ≈ −0.090 Ha/e**.

> **VC-1 (interior energy).** A large-N localised jellium at r_s should have an
> energy-per-electron in the deep interior tending to ε_HEG(r_s). For a *finite*
> cluster the approach is via the liquid-drop expansion (Part 4.3) — do not expect
> the small-N cluster to sit *on* ε_HEG; expect it to *extrapolate* to it.

---

## Part 2 — Electrostatics of the background (where INQ's G=0 trick lives)

**Learning objective:** understand why today's jellium is delocalised "for free",
and exactly what electrostatic object the new perturbation must add.

### 2.1 Charge neutrality is not optional
A periodic box with net charge has infinite Coulomb energy (the G=0 / monopole
term diverges). Plane-wave codes enforce neutrality by **dropping the G=0 Fourier
component** of the potential. For electrons-only (`extra_electrons(N)`, no ions),
this is *identical* to adding a uniform positive background of density N/V over the
**whole cell** — i.e. delocalised jellium. **This is the status quo in your repo.**

### 2.2 The localised model
Replace the whole-cell uniform background with a confined one:
```
sphere:  n₊(r) = n₀ · Θ(R_cl − |r − r₀|)
slab:    n₊(r) = n₀ · Θ(R_cl − |z − z₀|)      (infinite/periodic in x,y)
box:     n₊(r) = n₀ · Π(region)
```
Neutrality fixes the electron count to the background's charge:
```
sphere:  N = n₀ · (4/3)π R_cl³
slab:    N = n₀ · (2 R_cl) · A_xy
```
▢ **EXERCISE 2.2** — You want N = 40 electrons at r_s = 4. What R_cl (sphere)?
*Answer / check:* R_cl = r_s · N^{1/3} = 4 · 40^{1/3} = **13.68 a₀**. (General result:
`R_cl = r_s N^{1/3}` for a sphere — memorise this, it's your sizing formula.)

### 2.3 Potential of a uniform charged sphere (the confining well) — DERIVE THIS
This is the single most useful analytic check you have. A uniform sphere of
*positive* charge Q = N, radius R_cl, produces a potential energy for an electron
(charge −1):
```
            ⎧ −(N / 2R_cl)·(3 − r²/R_cl²)     r ≤ R_cl   (parabolic well inside)
v_bg(r) =   ⎨
            ⎩ −N / r                          r ≥ R_cl   (point-charge outside)
```
(Units: Hartree. Sign: electrons are *attracted*, so v_bg < 0.)

▢ **EXERCISE 2.3** — Derive v_bg(r) from Gauss's law. Then evaluate the **well depth
at the centre** and the **value at the edge** for N=40, R_cl=13.68 a₀.
*Answer / check:* centre `v_bg(0) = −3N/(2R_cl) = −4.39 Ha`; edge
`v_bg(R_cl) = −N/R_cl = −2.92 Ha`. (Centre is 1.5× deeper than the edge.)

> **VC-2 (well shape).** Have the implementation dump `v_bg(r)` along a line through
> r₀. Inside R_cl it must be a **downward parabola**; outside it must fall off as
> **−N/r**; the centre/edge ratio must be exactly **3/2**. If the interior isn't
> parabolic, the Θ-profile or the Poisson sign is wrong. *(With a softened edge —
> Part 8 — expect this only away from the edge transition region.)*

### 2.4 The key identity (why neutrality saves the G=0 drop)
INQ builds the Hartree potential as `v_H = poisson(n)` with G=0 dropped. We add
`v_bg = −poisson(n₊)` (also G=0 dropped). The electron then feels
```
v_es(r) = v_H + v_bg = poisson(n − n₊)
```
Because `∫n = ∫n₊ = N`, the difference `n − n₊` has **zero net charge**, so its G=0
component is zero *anyway* — dropping it is exact, not an approximation. This is the
crux of correctness and it is mathematically the GPAW jellium recipe.

▢ **EXERCISE 2.4** — Convince yourself that if `∫n ≠ ∫n₊`, the dropped G=0 term would
*not* cancel and the absolute energy would be gauge-dependent. (This is why a bug in
the electron count vs background charge shows up as a constant energy offset.)

### 2.5 Background self-energy E_self
The background interacts with *itself*. The classical self-energy of a uniform
sphere is
```
E_self = (3/5) N² / R_cl          [Ha]
```
Whether this appears in the reported total energy depends on the bookkeeping
(Part 3.4). You must know it because the cluster-energy → HEG-limit benchmark only
works once E_self is handled consistently.

▢ **EXERCISE 2.5** — E_self for N=40, R_cl=13.68 a₀.
*Answer / check:* (3/5)·1600/13.68 = **70.2 Ha**. (Large! — hence it must be tracked
explicitly, never left implicit.)

---

## Part 3 — Kohn–Sham DFT for jellium

**Learning objective:** know exactly which term v_bg is, and how the total energy
decomposes, so you can audit the code's energy print-out line by line.

### 3.1 The KS equations
```
[ −½∇² + v_KS(r) ] ψ_i = ε_i ψ_i ,     n(r) = Σ_i f_i |ψ_i|²
v_KS(r) = v_ext(r) + v_H(r) + v_xc(r)
```
For *our* system: **v_ext(r) = v_bg(r)** — the background well IS the external
potential. There are no nuclei. (Contrast a molecule, where v_ext is the ionic
pseudopotential.)

### 3.2 Where v_bg enters in INQ
- `self_consistency::update_hamiltonian` assembles `vscalar = vion + perturbation +
  v_H + v_xc`.
- We have no ions ⇒ `vion = 0`. The perturbation's `.potential(t, vscalar)` injects
  `v_bg`. Because the same `self_consistency` object is used by `ground_state` *and*
  `real_time`, v_bg is present in the SCF (electrons localise) and during the
  projectile flight (background stays put). **Time-independent ⇒ static well.**

### 3.3 Exchange–correlation: use LDA/LSDA
For jellium the LDA is not an approximation of convenience — the LDA energy density
is *defined* as ε_HEG of the local density. Use LDA (or LSDA if you want spin).
This is the literature-standard choice for jellium (Lang–Kohn used the LDA precursor;
Parr & Yang Ch. 7). PBE/GGA is fine too but Lang–Kohn benchmarks are LDA.

### 3.4 Total energy decomposition and double-counting
```
E_tot = T_s[n] + ∫ v_ext n  + E_H[n] + E_xc[n]   (+ E_self of background, see below)
      = Σ f_i ε_i − E_H[n] + (E_xc − ∫ v_xc n)   (+ background terms)
```
The electrostatic energy of the *neutral* system is
```
E_es = ½ ∫∫ (n−n₊)(n−n₊)/|r−r'|  =  E_H[n] + (∫ v_bg n + E_self)
                                       └ electron ┘ └ e–bg ┘ └ bg–bg ┘
```
Note the three pieces: electron–electron repulsion, electron–background attraction,
background–background self-repulsion (E_self, Part 2.5). **A correct implementation
must account for E_self** or the absolute total energy is meaningless (only energy
*differences* at fixed geometry would survive). This is exactly the "E_self" your
spec mentions for the HEG-limit benchmark.

### 3.5 Sign-convention audit (the thing most likely to be wrong in code)
Three signs must be mutually consistent:
1. n₊ is a **positive** charge ⇒ v_bg is an **attractive (negative)** well for
   electrons (Part 2.3).
2. In INQ, `poisson::solve(ρ)` returns the potential of a positive density ρ as a
   *positive* hump; electrons (charge −1) see `+v_H` (repulsion). So the background
   contributes `v_bg = −poisson(n₊)` (attraction). **Check the minus sign.**
3. The electron–background energy `∫ v_bg n` is **negative** (attraction).

> **VC-3 (sign sanity).** With the background on, the SCF must (a) converge, (b)
> bind the electrons *inside* R_cl (electron density peaks where n₊ is non-zero), and
> (c) give a negative electron–background energy term. If electrons flee to the box
> edges, the sign of v_bg is flipped.

---

## Part 4 — Benchmark A: spherical jellium cluster

**Learning objective:** predict the magic numbers and the energy-vs-N trend so you
can confirm the *finite* cluster is physically correct, not just numerically stable.

### 4.1 Spherical shells ≠ your current box shells ⚠ READ THIS
Your existing `inqkit/jellium/shells.hpp` enumerates **plane-wave shells of a cubic
periodic box** (magic numbers 2, 14, 38, 54, 66, 114, **162**). Those come from the
degeneracy of `(2π/L)²(nₓ²+nᵧ²+n_z²)` and are an artefact of the *delocalised*
periodic jellium. **A localised spherical cluster has DIFFERENT magic numbers**,
set by the spherical well's angular-momentum shells (1s, 1p, 1d, 2s, 1f, …):
```
spherical jellium magic N:  2, 8, 18, 20, 34, 40, 58, 92, 138, ...
```
(The famous Na-cluster series — Knight et al., PRL 52, 2141 (1984); de Heer,
Rev. Mod. Phys. 65, 611 (1993).) Confusing the two tables is the single most likely
way to "validate" against the wrong target.

▢ **EXERCISE 4.1** — Fill the spherical shell-filling table: order the levels
(1s,1p,1d,2s,1f,2p,…), give each its degeneracy 2(2ℓ+1), and accumulate to recover
the magic numbers above.
*Answer / check:* 1s(2)→2, 1p(6)→8, 1d(10)→18, 2s(2)→20, 1f(14)→34, 2p(6)→40, …

### 4.2 Self-consistent shell structure
At closed-shell N (2, 8, 20, 40, …) the cluster is especially stable: a clear HOMO–
LUMO gap, spherical density. Off magic numbers the system would Jahn–Teller distort
(not captured by a fixed spherical background). So **run validation clusters at
magic N** to get clean closed shells.

> **VC-4a (shells).** For a magic-N spherical cluster, the KS eigenvalue spectrum
> should show degenerate groups with sizes 2, 6, 10, 2, 14, … and a visible gap above
> the last filled shell. The radial electron density should be smooth, spherical,
> with the right number of radial nodes.

### 4.3 Energy → HEG limit (liquid-drop model)
A finite jellium drop obeys the liquid-drop expansion:
```
E(N)/N = ε_v + ε_s N^{-1/3} + ε_c' N^{-2/3} + ...
         └ volume ┘ └ surface ┘ └ curvature ┘
ε_v = ε_HEG(r_s)        (Part 1.4)
```
So a plot of E/N vs N^{-1/3} extrapolates (N→∞) to **ε_HEG(r_s)**, and its slope is
the surface contribution. *This* is the quantitative HEG-limit check — and it only
works after E_self (Part 2.5) is consistently included/excluded.

> **VC-4b (HEG limit).** Run a ladder of magic-N clusters (e.g. N = 8, 20, 40, 92,
> 138 at fixed r_s), plot E/N vs N^{-1/3}, linear-fit, and check the intercept equals
> ε_HEG(r_s) from Exercise 1.4 within a few %. This is the headline "ground-state
> energy → HEG limit with E_self" benchmark from your spec.

### 4.4 Interior density (the "few-% / increase R_cl" rule)
Deep inside a large cluster, n(r) → n₀. Near the surface it overshoots/undershoots
with **Friedel oscillations** (wavelength π/k_F) and spills out past R_cl. Your spec's
rule — *"if interior density deviates from n₀ by more than a few percent away from the
surface, increase R_cl"* — is literally asking: is the cluster big enough to have a
bulk-like interior? Small clusters are *all surface* and never show a flat n₀ region.

> **VC-4c (interior).** Plot n(r)/n₀ vs r. Away from the edge (say r < R_cl − 2π/k_F)
> it must be flat to within a few %. If not, R_cl is too small — increase N (and hence
> R_cl = r_s N^{1/3}) and repeat. Friedel wavelength `π/k_F` is your yardstick for
> "away from the surface".

---

## Part 5 — Benchmark B: Lang–Kohn slab

**Learning objective:** reproduce the canonical semi-infinite-jellium surface numbers
(density profile, work function, surface energy) — the most stringent test because it
fixes absolute energetics.

### 5.1 Geometry
A **slab**: n₊ = n₀ for |z| < R_cl, zero outside; uniform/periodic in x,y. This models
a semi-infinite metal surface (two surfaces, at z = ±R_cl). Reference: Lang & Kohn,
Phys. Rev. B **1**, 4555 (1970) — density profile & surface energy; Phys. Rev. B **3**,
1215 (1971) — work function.

### 5.2 Surface density profile
Across the surface the electron density:
- **spills out** beyond the background edge (electrons leak into vacuum),
- shows **Friedel oscillations** decaying into the bulk (wavelength π/k_F),
- creates a **dipole layer** (electrons outside, depleted positive just inside).

> **VC-5a (profile).** n(z) should be flat = n₀ in the slab centre, decay smoothly
> through ~1/k_F into vacuum, and ripple with period π/k_F on the metal side. Compare
> shape against Lang–Kohn Fig. 1 for your r_s.

### 5.3 Work function
```
Φ = v_vac(∞) − μ          (μ = chemical potential = Fermi level; v_vac far outside)
```
The dipole layer (5.2) raises the vacuum level relative to the interior; Φ is the
energy to remove an electron. Lang–Kohn 1971 tabulate Φ(r_s) (e.g. a few eV; rises
as r_s decreases). 

> **VC-5b (work function).** From the converged slab, read μ (Fermi energy) and the
> electrostatic potential plateau in vacuum; Φ = v_vac − μ. Compare to Lang–Kohn
> Φ(r_s). ⚠ VERIFY the exact tabulated value for your chosen r_s.

### 5.4 Surface energy and the 86.4 erg/cm² target
The surface energy σ is the energy cost per unit area of creating the surface:
```
σ = [ E_slab − N · ε_HEG(r_s) ] / (2 A_xy)      (factor 2: two surfaces)
```
- In a.u. σ is in **Ha/Bohr²**; convert with Part 0 (×1.5569×10⁶ → erg/cm²).
- Jellium surface energies are notoriously small and even **go negative for high
  density** (small r_s, e.g. Al r_s≈2.07) — the famous Lang–Kohn deficiency that
  motivated stabilised jellium. They are positive and ~tens of erg/cm² for low
  density (large r_s).

⚠ **VERIFY — which r_s gives 86.4 erg/cm²?** I do **not** have the Lang–Kohn / GPAW
table memorised reliably. 86.4 erg/cm² (= 5.55×10⁻⁵ Ha/Bohr²) is a *positive,
modest* value, so it corresponds to a **low-density (largish r_s)** jellium — plausibly
the r_s the GPAW jellium tutorial uses, but **confirm against the actual GPAW run that
produced 86.4** before adopting it as the gate. This is a worksheet to-do, not a
settled fact.

> **VC-5c (surface energy).** Compute σ from the slab via the formula above, convert
> to erg/cm², and compare to the GPAW-reproduced **86.4 erg/cm²** *at the matching
> r_s*. Getting the sign/magnitude right validates the whole electrostatic + XC chain.

---

## Part 6 — The interior-density convergence test (your "change the plan" gate)

This is the cheap, do-it-first benchmark and the one your spec calls out explicitly.

**Procedure:**
1. Pick r_s ⇒ n₀ (Part 1.2) and a target N ⇒ R_cl = r_s N^{1/3} (Part 2.2).
2. SCF the cluster/slab with the background on.
3. Plot n(r)/n₀ vs distance from centre.
4. In the window `r < R_cl − (a few × π/k_F)`, require |n/n₀ − 1| ≲ "a few %".
5. **If it fails, increase R_cl (i.e. increase N to the next magic number) and repeat.**

> **VC-6.** This is the master gate: no surface/energy benchmark is meaningful until
> the interior is bulk-like. A common failure is choosing N too small — then the
> cluster is all surface and Part 4.3 / Part 5 will look "wrong" when really the model
> is just too small.

---

## Part 7 — The projectile (dynamic phase) — brief, machinery already exists

Once the static localised jellium is validated, the projectile is fired at it from
far away. Your repo already has both projectile models:
- **wave-packet electron** (Gaussian pseudopotential / injected WP), and
- **classical ion** (e.g. antiproton ONCV).

The physics you'll read off:
- **energy loss / stopping**: electronic energy transferred to the jellium as the
  projectile traverses it (energy bookkeeping you already compute).
- **edge effects**: unlike infinite jellium, a localised target has *entrance* and
  *exit* surfaces — the projectile sees a finite interaction length ≈ 2R_cl.

This worksheet doesn't re-derive stopping power (the existing S(v) campaign covers
it); the *new* element is simply that the target is finite and localised, so the
projectile starts in true vacuum, enters, traverses ~2R_cl of bulk-like jellium, and
exits — a cleaner analogue of a real scattering experiment than the periodic gas.

(Decision still open in the session: projectile model, launch distance, and boundary
treatment — those land in the plan/ADR.)

---

## Part 8 — Numerical validation knobs (so you can tell physics bugs from grid bugs)

1. **Edge softening.** A sharp Θ-function on a finite grid causes Gibbs ringing in
   both n₊ and v_bg (the FFT-based Poisson solve amplifies it). The GPAW jellium
   template softens the edge (e.g. an erf/Fermi profile of width w over a couple of
   grid points). **Open decision**, but know that a too-sharp edge produces spurious
   density wiggles that masquerade as Friedel oscillations.
   - *Discriminator:* true Friedel oscillations have wavelength **π/k_F** and live on
     the metal side; Gibbs ringing tracks the **grid spacing** and sits right at the
     edge. If your "oscillations" change wavelength when you refine the grid, they're
     numerical.

2. **Grid spacing / cutoff.** Must resolve k_F and the edge width. Converge the
   interior density and total energy w.r.t. spacing before trusting any benchmark.

3. **Poisson / G=0.** Confirm the implementation drops G=0 consistently for both n
   and n₊ (Part 2.4). A residual net charge shows up as a constant energy drift with
   box size.

4. **Box size & periodic images.** A *localised neutral* target (electrons + finite
   background) has no monopole, so periodic images interact only via higher
   multipoles — much weaker than delocalised jellium. Still, leave vacuum ≳ a few a₀
   around the target, and (for the projectile) enough run-up + run-out distance. The
   slab benchmark needs enough vacuum in z that v reaches a flat plateau (for Φ).

> **VC-8.** Refine the grid by ~2× and confirm interior n₀, total energy, and any
> "oscillations" are stable. Anything that moves is numerical, not physical.

---

## Reference list (pin these down as you work)

- **Lang & Kohn, Phys. Rev. B 1, 4555 (1970)** — jellium surface density profile &
  surface energy (Benchmark B core).
- **Lang & Kohn, Phys. Rev. B 3, 1215 (1971)** — work function of jellium surfaces.
- **Knight et al., Phys. Rev. Lett. 52, 2141 (1984)** — electronic shell structure /
  magic numbers of metal clusters (Benchmark A).
- **Ekardt, Phys. Rev. B 29, 1558 (1984)** — self-consistent spherical jellium
  clusters.
- **Brack, Rev. Mod. Phys. 65, 677 (1993)** / **de Heer, Rev. Mod. Phys. 65, 611
  (1993)** — jellium cluster reviews (liquid-drop, shells).
- **Perdew & Zunger, Phys. Rev. B 23, 5048 (1981)** — LDA correlation parametrisation
  (ε_c for Part 1.4).
- **Parr & Yang, *Density-Functional Theory of Atoms and Molecules*** — KS DFT, LDA,
  energy decomposition (Part 3).
- **GPAW jellium documentation / `jellium.py`** — the "mask template" your spec
  references; the source of the 86.4 erg/cm² target ⚠ (confirm r_s there).
- **In-repo:** `inq-stack/include/inqkit/jellium/shells.hpp` — but note it is the
  *periodic-box* shell table (Part 4.1), NOT the spherical-cluster one.

---

## Worksheet self-test (do these before trusting any run)

| # | Check | Source |
|---|---|---|
| VC-1 | E/N extrapolates to ε_HEG(r_s) | Part 1.4 / 4.3 |
| VC-2 | v_bg interior parabolic, exterior −N/r, centre/edge = 3/2 | Part 2.3 |
| VC-3 | electrons bind inside R_cl; e–bg energy negative | Part 3.5 |
| VC-4a | spherical magic numbers 2,8,18,20,34,40,… (NOT 2,14,38…) | Part 4.1 |
| VC-4b | E/N vs N^{-1/3} intercept = ε_HEG | Part 4.3 |
| VC-4c | interior n/n₀ flat to a few % | Part 4.4 / 6 |
| VC-5a | slab profile: flat centre, spillout, Friedel π/k_F | Part 5.2 |
| VC-5b | Φ = v_vac − μ matches Lang–Kohn(r_s) ⚠ | Part 5.3 |
| VC-5c | σ = 86.4 erg/cm² at matching r_s ⚠ | Part 5.4 |
| VC-6 | interior gate passes before any other benchmark | Part 6 |
| VC-8 | refine grid 2×; physics stable | Part 8 |
