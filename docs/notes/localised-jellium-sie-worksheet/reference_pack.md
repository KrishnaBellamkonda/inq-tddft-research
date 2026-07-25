---
title: "Localised Jellium — SIE / Coulomb-vs-Distance / Long-Range Cutoff: Reference Pack"
subtitle: "Source material for the B+C+D worksheet (threads B, C, D of the GS parameter study campaign)"
author: "Localised-jellium GS study — campaign task 1 reference material"
date: "2026-06-27"
---

# How to use this document

This is the **authored reference pack** (the "material") for a derive-it-yourself
worksheet covering three threads of the *localised-jellium ground-state parameter
study* campaign:

- **Thread B — the SIE decomposition.** What extra energy a wavepacket (WP)
  projectile adds at $t=0$, and how to split it cleanly into a distance-independent
  self-interaction error $E_{\mathrm{SIE}}$ and a distance-dependent WP–jellium
  cross term $E_{\mathrm{cross}}(r)$.
- **Thread C — Coulomb-vs-distance & the classical subtraction.** The
  electrostatics of a charge-matched Gaussian projectile near the slab, and the
  subtle reason a *naïve* WP$-$classical subtraction fails (the "+798 eV" puzzle).
- **Thread D — the long-range cutoff.** A defensible cutoff for the classical
  projectile's radial Coulomb potential to avoid periodic-image / loop-around
  self-interaction.

It states **answers, methods, and the in-repo numbers** that any derivation must
reproduce. The companion `worksheet_plan.md` turns the boxed results below into
*derive-it-yourself* problems with validation checks and a self-test. Two
open-access source PDFs accompany this pack in `resources/`
(arXiv:2307.03213, arXiv:1805.01377).

Every numerical anchor is from an existing repo artefact (cited inline). Nothing
here is invented; where a statement is an inference rather than a sourced fact it
is labelled **Inference:**.

---

# Notation and units

Atomic units throughout ($\hbar = m_e = e = 4\pi\varepsilon_0 = 1$): length in Bohr
$a_0$, energy in Hartree (Ha). Conversion $1\ \mathrm{Ha} = 27.2114\ \mathrm{eV}$.
Report human-facing energies at **2–3 significant figures**.

| Symbol | Meaning |
|---|---|
| $\sigma_{\mathrm{WP}}$ | wavepacket width — the *wavefunction* std, $\psi \propto e^{-r^2/2\sigma_{\mathrm{WP}}^2}$. **All labels use this.** |
| $s \equiv \sigma_{\mathrm{pot}}$ | charge-density std $= \sigma_{\mathrm{WP}}/\sqrt2$ (so the WP density and the classical charge present the same cloud) |
| $n_{\mathrm{WP}}(\mathbf r)$ | WP electron density, $\propto e^{-r^2/\sigma_{\mathrm{WP}}^2}$, std $s=\sigma_{\mathrm{WP}}/\sqrt2$ |
| $\rho_{\mathrm{cl}}(\mathbf r)$ | classical Gaussian charge, charge $-1$, **matched** to $n_{\mathrm{WP}}$ |
| $E_{\mathrm{GS}}$ | total energy of the converged jellium slab (no projectile) |
| $E_{\mathrm{tot}}(0)$ | total energy of slab $+$ projectile at $t=0$ (single, non-self-consistent energy evaluation) |
| $\langle T_{\mathrm{WP}}\rangle$ | total kinetic energy of the wavepacket |
| $E_{\mathrm{SIE}}$ | one-electron self-interaction error of the WP (distance-independent) |
| $E_{\mathrm{cross}}(r)$ | WP–jellium electrostatic + XC cross term (distance-dependent, $\to 0$ as $r\to\infty$) |
| $r$ | WP-centroid → near slab-face distance |
| $r_s$ | density parameter; baseline slab $r_s \approx 5.67$ |

**Baseline slab** (source of truth `shared/configs/slab_n82_L50x50x90.hpp`): orthorhombic
$50\times50\times90$ Bohr, slab $50\times50$ face $\times 25$ Bohr thick (half-width
$d=12.5$, centred $z=0$), $N=82$, $n_0=1.31\times10^{-3}\,a_0^{-3}$, $r_s\approx5.67$,
spacing $0.50$ Bohr, $k_F\approx 0.34\,a_0^{-1}$. (A second slab, $N=234$, appears in the
Phase-3 stopping notebooks; the energy anchors in this pack come from that
Phase-3 σ=0.5, $E=100$ eV WP/classical pair —
`hypotheses/03_cap_stopping/qa_jellium_slab_baselines.ipynb` and
`hypotheses/qsp_phase3/results.json`.)

---

# Part B — The SIE decomposition

## B.1 The setup: a single $t=0$ energy evaluation

A wavepacket projectile is injected as **added electron density** $n_{\mathrm{WP}}$ onto
the converged jellium GS density $n_{\mathrm{GS}}$. The relevant quantity is the
**instantaneous** total energy of the combined density $n_{\mathrm{GS}}+n_{\mathrm{WP}}$
evaluated with the KS energy functional — **one Hamiltonian build, no SCF, no
propagation**. This is what makes the whole study cheap: each datum is a single-point
energy, not a trajectory.

> Charge bookkeeping (`qa_jellium_slab_baselines.ipynb §2`): the WP run has $N{+}1$
> electrons (bath $+$ one real WP electron), net charge $-1$, $G{=}0$-compensated by
> INQ's uniform background. The classical run keeps $N$ electrons and adds a
> **chargeless ghost** ($z_{\mathrm{valence}}{:}\,1\to0$) — a pure moving Gaussian
> *potential*, no extra electron, no contribution to neutrality, no SCF seed.

## B.2 The energy difference, term by term

With the KS total-energy functional
$$
E[n] = T_s[n] + \int v_{\mathrm{bg}}\,n + E_H[n] + E_{xc}[n] + E_{\mathrm{self}}^{\mathrm{bg}},
$$
the energy added by injecting the WP is
$$
E_{\mathrm{tot}}(0) - E_{\mathrm{GS}}
= \underbrace{\big(T_s[n_{\mathrm{GS}}{+}n_{\mathrm{WP}}] - T_s[n_{\mathrm{GS}}]\big)}_{=\ \langle T_{\mathrm{WP}}\rangle}
+ \underbrace{\int v_{\mathrm{bg}}\,n_{\mathrm{WP}}}_{\text{WP–background}}
+ \underbrace{\big(E_H[n_{\mathrm{GS}}{+}n_{\mathrm{WP}}] - E_H[n_{\mathrm{GS}}]\big)}_{\text{cross-Hartree} + \text{self-Hartree}}
+ \underbrace{\big(E_{xc}[n_{\mathrm{GS}}{+}n_{\mathrm{WP}}] - E_{xc}[n_{\mathrm{GS}}]\big)}_{\text{cross-XC} + \text{self-XC}} .
$$
Because the KS kinetic functional is additive over (orthogonal) orbitals, the first
bracket is exactly the kinetic energy of the added WP orbital, $\langle T_{\mathrm{WP}}\rangle$.
Expanding the Hartree difference,
$$
E_H[n_{\mathrm{GS}}{+}n_{\mathrm{WP}}] - E_H[n_{\mathrm{GS}}]
= \underbrace{\iint \frac{n_{\mathrm{GS}}(\mathbf r)\,n_{\mathrm{WP}}(\mathbf r')}{|\mathbf r-\mathbf r'|}}_{\text{cross-Hartree (r-dependent)}}
+ \underbrace{\tfrac12\iint \frac{n_{\mathrm{WP}}(\mathbf r)\,n_{\mathrm{WP}}(\mathbf r')}{|\mathbf r-\mathbf r'|}}_{\text{self-Hartree (r-independent)}} .
$$

## B.3 The wavepacket kinetic energy $\langle T_{\mathrm{WP}}\rangle$

For $\psi \propto \exp(-r^2/2\sigma_{\mathrm{WP}}^2)\,e^{i k_0 z}$ the total kinetic energy
splits into a **drift** part and a **zero-point** (localisation) part:
$$
\boxed{\ \langle T_{\mathrm{WP}}\rangle = \tfrac12 k_0^2 \;+\; \frac{3}{4\sigma_{\mathrm{WP}}^2}\ }\qquad(\text{Ha}).
$$
*Derivation of the zero-point term.* $\langle T\rangle = \tfrac12\langle|\nabla\psi|^2\rangle/\langle|\psi|^2\rangle$;
with $\nabla\psi=-(\mathbf r/\sigma_{\mathrm{WP}}^2)\psi$ (drift aside),
$\langle T\rangle_{\mathrm{zp}}=\tfrac12\langle r^2\rangle/\sigma_{\mathrm{WP}}^4$, and for the density
$|\psi|^2\propto e^{-r^2/\sigma_{\mathrm{WP}}^2}$ (variance $\sigma_{\mathrm{WP}}^2/2$ per axis)
$\langle r^2\rangle = 3\sigma_{\mathrm{WP}}^2/2$, giving $3/(4\sigma_{\mathrm{WP}}^2)$. $\square$

**Numbers (σ_WP=0.5, $E_{\mathrm{drift}}=100$ eV):** drift $\tfrac12 k_0^2 = 100$ eV,
zero-point $3/(4\cdot0.25)=3$ Ha $=81.6$ eV, total $\langle T_{\mathrm{WP}}\rangle = 181.6$ eV.
Repo anchor: `qsp_phase3/results.json` — analytic $181.6$ eV, **run-measured $180.8$ eV**
(`T_WP_run_Ha = 6.644`). The two agree to $\sim0.4$%.

> **Critical pitfall (repo-confirmed).** Using "+100 eV" as $\langle T_{\mathrm{WP}}\rangle$
> omits the $\sim$82 eV zero-point term and **overcounts the SIE by $\sim$82 eV**
> (you would report $\sim$85 eV instead of $\sim$4.5 eV). Always subtract the
> **measured** $\langle p^2\rangle/2$, not the nominal drift energy.

## B.4 The decomposition and the definition of $E_{\mathrm{SIE}}$

Grouping the $r$-independent pieces (intrinsic to the isolated WP) versus the
$r$-dependent pieces (interaction with the slab):
$$
\boxed{\ E_{\mathrm{tot}}(0) - E_{\mathrm{GS}} - \langle T_{\mathrm{WP}}\rangle
= E_{\mathrm{SIE}} + E_{\mathrm{cross}}(r)\ }
$$
with
$$
E_{\mathrm{SIE}} \equiv \underbrace{\tfrac12\!\iint\!\frac{n_{\mathrm{WP}}n_{\mathrm{WP}}}{|\mathbf r-\mathbf r'|}}_{E_H[n_{\mathrm{WP}}]}
+ \;E_{xc}[n_{\mathrm{WP}}]
\;\;(\text{one-electron SIE; Perdew–Zunger 1981}),
$$
$$
E_{\mathrm{cross}}(r) = \int v_{\mathrm{bg}}\,n_{\mathrm{WP}}
+ \iint\frac{n_{\mathrm{GS}}\,n_{\mathrm{WP}}}{|\mathbf r-\mathbf r'|}
+ \big(E_{xc}\ \text{cross term}\big)\ \xrightarrow{r\to\infty}\ 0 .
$$

**Why $E_{\mathrm{SIE}}$ is the self-interaction error.** For a true one-electron density
the exact functional gives $E_H[n_1]+E_{xc}[n_1]=0$. LDA breaks this; the non-zero
residual *is* the SIE. It depends only on the WP **shape** ($\sigma_{\mathrm{WP}}$),
not on where the WP sits — hence distance-independent.

**Why $E_{\mathrm{cross}}(r)\to0$.** Far from a **neutral** slab, $n_{\mathrm{WP}}$ sits in
near-zero jellium density, so the XC cross term vanishes; and a neutral slab has no
monopole, so the electrostatic cross term ($\int v_{\mathrm{bg}}n_{\mathrm{WP}}$, attraction,
plus $\iint n_{\mathrm{GS}}n_{\mathrm{WP}}/r$, repulsion) cancels to leading order. **Inference:**
the surviving far-field is the slow image/multipole tail, $O(0.5\text{ eV})$ at $r\sim10$
Bohr (consistent with the measured $0.47$ eV change below).

> **The bug thread B fixes.** The *old* estimate took the entire LHS
> $E_{\mathrm{tot}}(0)-E_{\mathrm{GS}}-\langle T_{\mathrm{WP}}\rangle$ and called it "SIE." At
> finite $r$ that quantity still carries $E_{\mathrm{cross}}(r)$. The corrected SIE is the
> $r\to\infty$ limit (route 1) or the Hartree-matched classical subtraction (route 2).

## B.5 Route 1 — far-launch / vacuum control (the plateau)

Push the WP far from the slab (or remove the slab entirely — the "vacuum-WP
control"); then $E_{\mathrm{cross}}\to0$ and the LHS **plateaus** at $E_{\mathrm{SIE}}$.
This is already partly measured (`qa_jellium_slab_baselines.ipynb §5.1`):

| run | launch $z$ | dist → face | $E_{\mathrm{tot}}(0)$ (Ha) | $\langle T_{\mathrm{WP}}\rangle$ (eV) | excess $=E_{\mathrm{tot}}(0)-E_{\mathrm{GS}}-\langle T_{\mathrm{WP}}\rangle$ |
|---|---|---|---|---|---|
| p3_wp (far) | $-23.0$ | 10.5 Bohr | $-154.180$ | 180.8 | $+0.167$ Ha $= \mathbf{+4.55}$ **eV** |
| p5_wp (near) | $-15.5$ | 3.0 Bohr | $-154.162$ | 180.8 | $+0.185$ Ha $= +5.02$ eV |

**Reading:** moving the WP $7.5$ Bohr farther changes the excess by only $0.47$ eV, so
$E_{\mathrm{cross}}$ is small and the excess is dominated by genuine SIE:
$\boxed{E_{\mathrm{SIE}}\approx 4.5\ \mathrm{eV}}$ (σ_WP=0.5). A true vacuum-WP control (no
slab) would remove the residual $\sim$0.5 eV image/Hartree tail for a Hartree-free
number. The far-launch test is cited to **Nazarov & Gross 2025** and the SIE
literature (`docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`).

## B.6 $\sigma$-scaling of the SIE (why it must be re-measured per $\sigma$)

The **bare** Gaussian self-Hartree (Part C.1) is $E_H[n_{\mathrm{WP}}]=q^2/(2s\sqrt\pi)$ with
$s=\sigma_{\mathrm{WP}}/\sqrt2$, i.e. $E_H \propto 1/\sigma_{\mathrm{WP}}$. For $\sigma_{\mathrm{WP}}=0.5$,
$s=0.354$, $E_H[n_{\mathrm{WP}}]\approx 0.80$ Ha $\approx \mathbf{22\ eV}$. XC cancels most of
it, leaving the net SIE $\approx4.5$ eV. **The 22 eV (bare self-Hartree) and the
4.5 eV (net SIE after XC) are different quantities** — do not conflate them; the
"$\sim$7–21 eV" figures in older notes refer to the bare/partly-cancelled range.

**Consequence:** a larger, near-rigid WP ($\sigma_{\mathrm{WP}}\approx3$) has
$E_H\propto1/\sigma$ about $6\times$ smaller, so SIE $<1$ eV, *and* its zero-point KE
$3/(4\cdot9)\approx2$ eV is tiny — which is the campaign's motivation for the
wide-WP direction. **The SIE is $\sigma$-dependent and must be re-measured for each
$\sigma$ used** (`qa_jellium_slab_baselines.ipynb §5.1`, §7 #1).

---

# Part C — Coulomb-vs-distance & the classical subtraction

## C.1 The Gaussian charge: density, potential, self-energy

A normalised 3-D Gaussian charge of total charge $q$ and density std $s$:
$$
\rho(\mathbf r)=\frac{q}{(2\pi s^2)^{3/2}}\,e^{-r^2/2s^2},\qquad
V(r)=\frac{q}{r}\,\mathrm{erf}\!\Big(\frac{r}{s\sqrt2}\Big),\qquad
E_{\mathrm{self}}=\frac{q^2}{2 s\sqrt\pi}.
$$
*Self-energy derivation (Fourier).* $\rho(\mathbf k)=q\,e^{-k^2 s^2/2}$, so
$E_{\mathrm{self}}=\tfrac12\!\int\!\frac{d^3k}{(2\pi)^3}\frac{4\pi}{k^2}|\rho(\mathbf k)|^2
=\tfrac12(4\pi)q^2\frac{1}{4\pi^{3/2}s}=\frac{q^2}{2s\sqrt\pi}$. $\square$
The $V(r)$ form is exactly the INQ classical-projectile pseudopotential
(`inqview.io.gaussian_psp`, $V=\mathrm{erf}(r/(\sigma_{\mathrm{pot}}\sqrt2))/r$,
$\sigma_{\mathrm{pot}}=s$).

## C.2 Charge matching: $\sigma_{\mathrm{pot}}=\sigma_{\mathrm{WP}}/\sqrt2$

To present the *same cloud* to the bath, the classical charge std must equal the WP
**density** std:
$$
s = \sigma_{\mathrm{pot}} = \sigma_{\mathrm{WP}}/\sqrt2,\qquad
\rho_{\mathrm{cl}}\propto e^{-r^2/\sigma_{\mathrm{WP}}^2}=n_{\mathrm{WP}}.
$$
`generate_gaussian_psp(sigma_wp)` builds the erf charge at std $\sigma_{\mathrm{WP}}/\sqrt2$
internally, so a WP and a classical run at the **same labelled $\sigma_{\mathrm{WP}}$** are
genuinely matched (`docs` memory `reference_sigma_matching_convention`; CONTEXT.md
"σ-convention unification"). **Caution:** legacy `electron_gaussian_sigmaXpY.upf` use
the *old* convention (charge std $=$ label); multiply their label by $\sqrt2$ to get
$\sigma_{\mathrm{WP}}$.

## C.3 The jellium slab potential (the charged-plate skeleton)

For an infinite (in $x,y$) slab of uniform positive background $n_+$ on $|z|<d$, 1-D
Poisson $V''(z)=-4\pi\rho$ gives a **parabolic well inside, flat outside**:
$$
v_{\mathrm{bg}}(z) = \begin{cases}
-2\pi n_+\,(d^2 - z^2) + C, & |z|<d,\\[2pt]
\text{const (flat)}, & |z|>d,
\end{cases}
$$
(electrons, charge $-1$, see this as an attractive well). The *neutral* slab
(electrons $+$ background) has **no far field**: $V\to$ const as $|z|\to\infty$. The
GS density's surface profile (Friedel tail + spill-out) furnishes the alternating
"$+/-$ plate" picture of thread A — the projectile interacts with this neutral plate
stack, which is why $E_{\mathrm{cross}}$ is small. (Benchmarks: Lang–Kohn slab,
worksheet `docs/notes/localised-jellium-theory.md` Part 5.)

## C.4 $E_{\mathrm{cross}}(r)$: a charge interacting with a *neutral* slab

The WP electron feels the **net** (electron $+$ background) slab potential, which is
screened to $\sim0$ outside and nearly flat inside, so
$E_{\mathrm{cross}}(r)=\int n_{\mathrm{WP}}(\mathbf r - r\hat z)\,V_{\mathrm{net}}(\mathbf r)\,d^3r$ is
small and falls off with $r$ (measured: $0.47$ eV over $7.5$ Bohr, B.5). Far away the
leading survivor is the **image** attraction (C.7).

## C.5 The classical-ghost asymmetry — the "+798 eV" puzzle

Here is the crux of thread C. Compare the measured $t=0$ ledger
(`qa_jellium_slab_baselines.ipynb §5`):

| state | $E_{\mathrm{tot}}$ (Ha) | above slab GS |
|---|---|---|
| GS jellium slab | $-160.992$ | 0 (ref) |
| WP run, $t=0$ | $-154.162$ | $+185.9$ eV |
| **classical run, $t=0$** | $-131.654$ | $\mathbf{+798.3}$ **eV** |

The WP's $+185.9$ eV $\approx 100$ (drift) $+82$ (zero-point) $+4$ (SIE/cross) — all
physical. The classical's $+798$ eV is **not** physical and **not** $\sim$100 eV.
**Why:** the chargeless ghost ($z_{\mathrm{valence}}=0$) contributes its erf potential to
$v_{\mathrm{ext}}$ and so adds $\int v_{\mathrm{ghost}}\,n_{\mathrm{GS}}$ — the ghost interacting
with the slab's **bare, unscreened electrons** — but it carries **no charge for
neutrality**, so the compensating **ghost–background** term $\int v_{\mathrm{ghost}}\,n_+$
is *omitted*. The bare electron repulsion is large and uncancelled; applied to the
**unrelaxed** loaded GS density it is a sudden, non-variational jump.

$$
E_{\mathrm{classical}}(0) - E_{\mathrm{GS}} \;=\; \underbrace{\int v_{\mathrm{ghost}}\,n_{\mathrm{GS}}}_{\text{unscreened, large}}
\;\;\neq\;\; E_{\mathrm{cross}}^{\mathrm{WP}}(r)\;(\text{screened, small}).
$$

This is the user's own diagnosis (`§5` reader observation): *"the starting energies …
are not comparable … the Hartree term must be made almost the same."* It is also,
read precisely, thread B's statement that the current SIE estimate "omits the
classical [ghost–jellium] repulsion."

## C.6 The corrected classical subtraction (route 2)

To make WP and classical comparable so that
$E_{\mathrm{SIE}} = E_{\mathrm{WP}}(0) - E_{\mathrm{classical, corrected}}(0) - \langle T_{\mathrm{WP}}\rangle$,
**one correction is mandatory — re-adding the ghost–background term:**

$$
E_{\mathrm{classical, corrected}}(0) = E_{\mathrm{classical}}(0) + \int v_{\mathrm{ghost}}\,n_+ ,
$$

a Gaussian charge $-1$ against the uniform slab background (numeric grid integral;
closed-form erf cross-check). This restores the **neutral**-slab interaction the WP
actually feels (the WP, a real electron, gets both electron-repulsion and
background-attraction; the chargeless ghost gets only the bare electron-repulsion,
because $z_{\mathrm{valence}}=0$ drops the $\int v_{\mathrm{ghost}}n_+$ attraction).

> **Launch-far alone does NOT fix route 2** (unlike route 1). The bare ghost–electron
> Coulomb decays only as $\sim N/r$: at $r=40$ Bohr it is still
> $82/40\approx2.0$ Ha $\approx 56$ eV above GS — far from the WP's $\sim$0.5 eV. The
> WP's far-launch works *because* it sees the neutral slab; the ghost does not until
> the $\int v_{\mathrm{ghost}}n_+$ term is restored. So the ghost–background correction
> is **required at every $r$**, large or small (it is large: ~56 eV at $r{=}40$,
> hundreds near contact).

**Cross-check:** corrected route 2 must agree with route 1's plateau ($\approx4.5$ eV,
σ_WP=0.5) at every $r$, and $E_{\mathrm{cross}}^{\mathrm{WP}}(r)$ must equal the corrected
$E_{\mathrm{cross}}^{\mathrm{cl}}(r)$ (proving the WP–slab interaction is purely classical).
**If they disagree, an assumption is wrong** (unmatched ghost–background term, the
ghost self-energy bookkeeping, or an unconverged plateau). This mutual check is the
campaign's falsifiable validation of the energy reference handed to Campaign 1.

## C.7 Image-charge response (the dynamic regime)

Once electrons **relax** (during propagation, or in a self-consistent treatment), a
charge $q$ at distance $d$ outside the slab feels the classical **image** attraction
$$
W_{\mathrm{im}}(d) = -\frac{q^2}{4(d-d_0)},
$$
with the image plane $d_0$ a fraction of a Bohr outside the jellium edge (Lang–Kohn).
This is the leading long-range survivor of $E_{\mathrm{cross}}$ and the physical
content of the slow approach; it is **absent** in the frozen $t=0$ evaluation (no
response yet) but matters for the dynamic stopping (Campaign 1) and sets the scale
the cutoff (Part D) must not corrupt. (Jackson, *Classical Electrodynamics*; Lang &
Kohn 1973 image-plane position.)

---

# Part D — The long-range cutoff

## D.1 The periodic-image / loop-around problem

The classical projectile's potential $V(r)=-\,\mathrm{erf}(r/s\sqrt2)/r$ has a $-1/r$
tail reaching all space. In a periodic cell it interacts with (i) **its own periodic
images** and (ii) the **images of the slab**, and a propagating projectile that exits
one face **re-enters** the opposite face ("loop-around"), corrupting the energy
dynamics (`classical_projectile_fix.md §description`). INQ's $G{=}0$ drop neutralises
the net cell charge with a uniform background but **does not** remove these
image/wrap interactions.

## D.2 Minimum-image constraint

A truncation/cutoff radius $R_c$ must satisfy
$$
R_c < L_{\min}/2
$$
($L_{\min}$ the smallest box dimension the projectile can approach an image along) so
the projectile never sees a periodic copy. For the baseline ($L_x=L_y=50$, $L_z=90$):
$R_c < 25$ Bohr is set by the in-plane size, not $L_z$.

## D.3 Physical interaction range (screening sets the floor)

$R_c$ must also be **large enough** to keep the real physics. Linear screening in the
electron gas has the Thomas–Fermi length
$$
k_{\mathrm{TF}} = \sqrt{4 k_F/\pi},\qquad \lambda_{\mathrm{TF}} = 1/k_{\mathrm{TF}} .
$$
For $r_s\approx5.67$: $k_F=(9\pi/4)^{1/3}/r_s\approx0.34\,a_0^{-1}$,
$k_{\mathrm{TF}}\approx0.66\,a_0^{-1}$, $\lambda_{\mathrm{TF}}\approx1.5$ Bohr. So the screened
interaction decays over $\sim$1.5 Bohr inside the metal; the **image tail** outside
decays as $1/d$ and is the binding constraint. A defensible window is
$$
\boxed{\ \lambda_{\mathrm{TF}} \ll R_c < L_{\min}/2,\quad\text{e.g. }R_c\sim10\text{–}20\ \text{Bohr for the baseline.}\ }
$$

## D.4 Cutoff prescriptions (options to choose among)

1. **Geometric (no kernel change):** elongate $L_z$ and/or **stop the analysis before
   the projectile re-crosses** a face — the cleanest for a transit-only run; pairs
   with the "transit-only / elongated box" fix in `classical_projectile_fix.md`.
2. **Truncated Coulomb kernel:** replace $1/r$ by a cutoff kernel for $r>R_c$
   (spherical truncation / Wolf-type). Removes images directly; watch the energy
   discontinuity at $R_c$ (use a smoothly switched form).
3. **Martyna–Tuckerman (MT) nonperiodic Poisson** (J. Chem. Phys. **110**, 2810,
   1999): solves open boundary conditions for the cluster + projectile, eliminating
   images without a hard cutoff. Standard image-artefact removal for charged/aperiodic
   systems.

**Inference (recommended default):** for a *transit-only, t=0-anchored* localised-slab
run, prescription 1 (elongated box $+$ truncation of the analysis window) is the
lowest-risk and needs no engine change (`inq/` is immutable); reserve MT / truncated
kernel for cases where the projectile must remain in-box through equilibration.

## D.5 Finite-size error — what to converge and report

The finite-size literature quantifies the residual. The accompanying
**arXiv:2307.03213** (npj Comput. Mater. 2023) reports a **plasmon-cutoff finite-size
error $\sim$8 %** from periodic-image re-crossing into excited density, and recommends
trajectory sampling + Coulomb-cutoff / MT corrections. **arXiv:1805.01377** examines
RT-TDDFT stopping convergence (box, $k$-points, projectile representation). The
prescription handed downstream must state, with numbers, the chosen $R_c$ (or box
$L_z$ and stop-time), and a convergence check (vary $L_z$ / $R_c$, confirm the energy
reference and $S$ change by less than a stated tolerance).

## D.6 The prescription for Campaign 1

Deliverable to `classical_projectile_fix.md`: (a) the chosen scheme (1/2/3) **with
values** ($R_c$ or $L_z$+stop-time); (b) the minimum-image and screening checks
($\lambda_{\mathrm{TF}}\ll R_c<L_{\min}/2$); (c) a finite-size convergence statement with a
numeric tolerance; (d) a one-line note that the static $t=0$ energy reference
(Part B/C) is image-insensitive to leading order (constant offset cancels in
fixed-geometry differences), while the **dynamic** stopping requires the cutoff to
avoid loop-around.

---

# Worked-number sanity table (repo anchors — any derivation must reproduce these)

| Quantity | Value | Source |
|---|---|---|
| $\langle T_{\mathrm{WP}}\rangle$ (σ_WP=0.5, 100 eV), analytic | 181.6 eV | `qsp_phase3/results.json` |
| $\langle T_{\mathrm{WP}}\rangle$, run-measured | 180.8 eV | `qsp_phase3/results.json` |
| zero-point KE $3/(4\sigma_{\mathrm{WP}}^2)$ (σ_WP=0.5) | 81.6 eV (3 Ha) | derived; `qa_…baselines §5` |
| bare self-Hartree $E_H[n_{\mathrm{WP}}]=q^2/(2s\sqrt\pi)$ (σ_WP=0.5) | $\approx$22 eV (0.80 Ha) | derived (C.1) |
| **net SIE** (far-launch, σ_WP=0.5) | **4.5 eV** | `qa_…baselines §5.1` |
| LHS $E_{\mathrm{tot}}(0)-E_{\mathrm{GS}}-\langle T_{\mathrm{WP}}\rangle$ @ production $r$ | 3.93 eV | `qsp_phase3/results.json` |
| $E_{\mathrm{cross}}$ @ production $r$ (LHS $-$ SIE) | $\approx -0.5$ eV ($\approx$0) | `qsp_phase3_study.ipynb` |
| classical ghost $t=0$ excess (unmatched) | +798 eV | `qa_…baselines §5` |
| WP $t=0$ excess | +185.9 eV | `qa_…baselines §5` |
| Thomas–Fermi $\lambda_{\mathrm{TF}}$ ($r_s$=5.67) | $\approx$1.5 Bohr | derived (D.3) |
| finite-size plasmon-cutoff error | $\sim$8 % | arXiv:2307.03213 |

---

# References

**In-repo source notes** (`docs/sources/`):

- `nazarov-gross-2025-quantum-projectile-stopping.md` — quantum projectile stopping;
  one-electron SIE in TDDFT stopping.
- `correa-2018-electronic-stopping-power.md` — RT-TDDFT electronic stopping (review).
- `quijada-2007-cluster-bulk-stopping.md` — cluster vs bulk stopping.
- `stopping-power-formulae.md`, `stopping-power-jellium-anchors.md` — Lindhard/RPA
  anchors and formulae.

**In-repo worksheet:** `docs/notes/localised-jellium-theory.md` — Parts 2 (background
electrostatics, $G{=}0$/neutrality/self-energy), 3 (KS energy decomposition + sign
audit), 5 (Lang–Kohn slab), 8 (numerical knobs).

**External:**

- J. P. Perdew & A. Zunger, *Self-interaction correction to density-functional
  approximations for many-electron systems*, Phys. Rev. B **23**, 5048 (1981) — the
  one-electron SIE.
- N. D. Lang & W. Kohn, Phys. Rev. B **1**, 4555 (1970); **3**, 1215 (1971); **7**,
  3541 (1973) — jellium surface, work function, surface energy, image-plane position.
- G. J. Martyna & M. E. Tuckerman, J. Chem. Phys. **110**, 2810 (1999) — reciprocal-space
  Poisson solver for nonperiodic boundary conditions (image-artefact removal).
- J. D. Jackson, *Classical Electrodynamics* — image charge near a conductor, slab
  electrostatics.
- *Trajectory sampling and finite-size effects in first-principles stopping power
  calculations*, npj Comput. Mater. (2023), **arXiv:2307.03213** (in `resources/`).
- *Examining real-time time-dependent DFT for stopping power*, **arXiv:1805.01377**
  (in `resources/`).
