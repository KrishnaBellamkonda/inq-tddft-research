# Worksheet plan — Localised jellium: SIE / Coulomb-vs-distance / long-range cutoff (B + C + D)

> This is a **plan/spec**, not the worksheet itself. A worksheet-building agent
> consumes this file plus the resources in `resources/` and produces a polished,
> derive-it-yourself worksheet in the project house style. Tags below are
> structural; prose inside them is the instruction to the builder.

<identity>
You are building a graduate-level, **derive-it-yourself physics worksheet** for a
first-principles (TD)DFT researcher. The worksheet teaches the energetics of a
Gaussian projectile interacting with a localised jellium slab, so the reader can
(i) define and measure the wavepacket self-interaction error, (ii) do the
WP$-$classical "Coulomb subtraction" correctly, and (iii) choose a defensible
long-range cutoff. Match the existing house style exactly.
</identity>

<context>
- Campaign: *Localised jellium GS — parameter study + analytical mental models*
  (`docs/campaigns/localised_jellium/ground_state_parameter_study.md`), threads
  **B** (SIE decomposition), **C** (Coulomb-vs-distance), **D** (long-range cutoff).
- This worksheet is on the **critical path**: it produces the corrected energy
  reference `E_jellium(0) = E_tot(0) − ⟨T_WP⟩ − E_SIE` and the cutoff prescription
  that **unblock Campaign 1** (`classical_projectile_fix.md`).
- The physics is already partly computed in the repo; the worksheet **formalises and
  derives** it (it must reproduce the in-repo numbers, not invent new ones).
</context>

<inputs>
- **`resources/localised_jellium_sie_reference_pack.docx`** — the authored reference
  pack: every boxed result, derivation method, and in-repo numeric anchor. Each
  worksheet problem below cites the reference-pack section that holds its answer.
- **`resources/arxiv_2307.03213_finite_size_stopping.pdf`** — finite-size / cutoff
  errors in first-principles stopping (npj Comput. Mater. 2023). Source for Part D.5.
- **`resources/arxiv_1805.01377_rttddft_stopping.pdf`** — RT-TDDFT stopping
  convergence. Supporting source for Part D.
- **`docs/notes/localised-jellium-theory.md`** — the existing 472-line worksheet
  (Parts 0–8). **This new worksheet is its continuation**; reuse its notation, its
  "VC-N" validation-check callouts, and its self-test format. Do NOT restate Parts
  0–8; cross-reference them (e.g. background electrostatics → Part 2; KS energy
  decomposition → Part 3; Lang–Kohn slab → Part 5; numerical knobs → Part 8).
</inputs>

<house_style>
Mirror `localised-jellium-theory.md`:
- Each Part opens with **"Learning objective:"** (one sentence).
- Derivations the reader must do are flagged **"DERIVE THIS"** with the *method* hinted
  and the *boxed answer* deferred to a collapsible/footnote "answer" the reader checks
  against (answers come from the reference pack).
- Validation checks are blockquoted **"> VC-N (name)."** callouts stating an
  observable consequence (e.g. "the plateau must be flat to within X").
- Atomic units; state the `1 Ha = 27.2114 eV` conversion once. Report human-facing
  numbers at **2–3 significant figures**.
- $\sigma$ always means $\sigma_{\mathrm{WP}}$ (wavefunction std); surface
  $\sigma_{\mathrm{pot}}=\sigma_{\mathrm{WP}}/\sqrt2$ only in a methods aside.
- Label any non-sourced statement **"Inference:"**.
- End with a **Worksheet self-test** (answers known) and a **Reference list**.
</house_style>

<worksheet>

<part id="B" title="The wavepacket self-interaction error (SIE)">
<learning_objective>
Define the one-electron SIE, derive the clean decomposition
$E_{\mathrm{tot}}(0)-E_{\mathrm{GS}}-\langle T_{\mathrm{WP}}\rangle = E_{\mathrm{SIE}} + E_{\mathrm{cross}}(r)$,
and learn the two ways to measure $E_{\mathrm{SIE}}$.
</learning_objective>

<problem id="B1" type="setup">
Explain why $E_{\mathrm{tot}}(0)$ is a **single non-SCF energy evaluation** of
$n_{\mathrm{GS}}+n_{\mathrm{WP}}$, and lay out the WP-vs-classical charge bookkeeping
(WP $=$ real extra electron, net $-1$; classical $=$ chargeless ghost potential).
[ref pack B.1]
</problem>

<problem id="B2" type="derive">
**DERIVE THIS:** expand $E_{\mathrm{tot}}(0)-E_{\mathrm{GS}}$ term by term ($T_s$,
$\int v_{\mathrm{bg}}n$, $E_H$ split into cross- and self-Hartree, $E_{xc}$ cross/self).
Identify which terms are $r$-independent. [ref pack B.2]
</problem>

<problem id="B3" type="derive">
**DERIVE THIS:** $\langle T_{\mathrm{WP}}\rangle = \tfrac12 k_0^2 + 3/(4\sigma_{\mathrm{WP}}^2)$.
Show the zero-point term from $\langle T\rangle=\tfrac12\langle r^2\rangle/\sigma_{\mathrm{WP}}^4$
with $\langle r^2\rangle=3\sigma_{\mathrm{WP}}^2/2$. Evaluate for $\sigma_{\mathrm{WP}}=0.5$,
$E_{\mathrm{drift}}=100$ eV → 181.6 eV. [ref pack B.3]
<vc>VC: the analytic 181.6 eV must match the run-measured 180.8 eV (≈0.4%).
Using "+100 eV" instead overcounts SIE by ~82 eV — make the reader compute that error.</vc>
</problem>

<problem id="B4" type="derive">
**DERIVE THIS:** group into $E_{\mathrm{SIE}}=E_H[n_{\mathrm{WP}}]+E_{xc}[n_{\mathrm{WP}}]$
(Perdew–Zunger one-electron SIE, $r$-independent) and $E_{\mathrm{cross}}(r)\to0$.
Argue *why* the cross term vanishes for a **neutral** slab (no monopole; screened
interior). State explicitly the bug this fixes: the old estimate called the whole LHS
"SIE". [ref pack B.4]
</problem>

<problem id="B5" type="measure">
**Route 1 — far-launch / vacuum plateau.** Show that pushing the WP far makes the LHS
plateau at $E_{\mathrm{SIE}}$. Use the repo table (far 10.5 Bohr → +4.55 eV; near 3.0
Bohr → +5.02 eV) to conclude $E_{\mathrm{SIE}}\approx4.5$ eV (σ_WP=0.5). [ref pack B.5]
<vc>VC (plateau flatness): a 7.5-Bohr move changes the excess by only 0.47 eV ⇒ the
residual cross term is small and the plateau is real. A true vacuum control (no slab)
removes the last ~0.5 eV image tail.</vc>
</problem>

<problem id="B6" type="scaling">
**DERIVE THIS:** $E_H[n_{\mathrm{WP}}]=q^2/(2s\sqrt\pi)\propto1/\sigma_{\mathrm{WP}}$
($s=\sigma_{\mathrm{WP}}/\sqrt2$). Evaluate ≈22 eV for σ_WP=0.5 and contrast with the
net SIE ≈4.5 eV (XC cancels most of the bare self-Hartree). Conclude the SIE is
σ-dependent and must be re-measured per σ; a σ_WP≈3 WP has SIE <1 eV. [ref pack B.6]
</problem>
</part>

<part id="C" title="Coulomb-vs-distance and the classical subtraction">
<learning_objective>
Compute the electrostatics of a charge-matched Gaussian projectile near the slab, and
understand why a naïve WP$-$classical subtraction fails — then fix it.
</learning_objective>

<problem id="C1" type="derive">
**DERIVE THIS (Fourier):** for $\rho=q(2\pi s^2)^{-3/2}e^{-r^2/2s^2}$, show
$V(r)=(q/r)\,\mathrm{erf}(r/s\sqrt2)$ and $E_{\mathrm{self}}=q^2/(2s\sqrt\pi)$. Connect
$V(r)$ to the INQ pseudopotential form. [ref pack C.1]
</problem>

<problem id="C2" type="convention">
State the matching $\sigma_{\mathrm{pot}}=\sigma_{\mathrm{WP}}/\sqrt2$ and show both clouds
become $\propto e^{-r^2/\sigma_{\mathrm{WP}}^2}$. Warn about legacy UPF convention (×√2).
[ref pack C.2]
</problem>

<problem id="C3" type="derive">
**DERIVE THIS:** the uniform-slab potential $v_{\mathrm{bg}}(z)$ (parabolic well inside,
flat outside) from 1-D Poisson. Connect to the alternating "+/− plate" picture (thread
A) and to why a neutral slab has no far field. Cross-ref existing Part 2/5. [ref pack C.3]
</problem>

<problem id="C4" type="reason">
Argue that $E_{\mathrm{cross}}(r)$ is small because the WP feels the **net neutral** slab
potential (screened); quote the measured 0.47-eV change over 7.5 Bohr. [ref pack C.4]
</problem>

<problem id="C5" type="derive">
**THE KEY PROBLEM.** Explain the measured **+798 eV** (classical) vs **+185.9 eV** (WP)
$t=0$ excess. Show $E_{\mathrm{classical}}(0)-E_{\mathrm{GS}}=\int v_{\mathrm{ghost}}n_{\mathrm{GS}}$
(bare, unscreened electrons) because the chargeless ghost's compensating
$\int v_{\mathrm{ghost}}n_+$ term is omitted — so it is **not** comparable to the WP's
screened cross term. [ref pack C.5]
<vc>VC (comparability): before any WP−classical subtraction, the two runs' Hartree/ghost
interaction with the slab must be matched; otherwise the difference is dominated by the
unphysical +798 eV jump, not by physics.</vc>
</problem>

<problem id="C6" type="measure">
**Route 2 — corrected classical subtraction.** Show that the ghost–background term
$\int v_{\mathrm{ghost}}n_+$ is **mandatory** (launch-far alone fails: the bare
ghost–electron Coulomb decays only as $\sim N/r$, ~56 eV even at r=40). Re-add it:
$E_{\mathrm{SIE}}=E_{\mathrm{WP}}(0)-[E_{\mathrm{classical}}(0)+\int v_{\mathrm{ghost}}n_+]-\langle T_{\mathrm{WP}}\rangle$.
[ref pack C.6]
<vc>VC (cross-check, the campaign's falsifiable test): corrected route 2 must agree with
route 1's plateau (≈4.5 eV, σ_WP=0.5) at every r, AND E_cross^WP(r) == corrected
E_cross^cl(r). Disagreement ⇒ unmatched ghost–background/self-energy term or an
unconverged plateau.</vc>
</problem>

<problem id="C7" type="concept">
Introduce the image potential $W_{\mathrm{im}}(d)=-q^2/4(d-d_0)$ (electrons relaxed):
the leading long-range survivor, absent at frozen $t=0$ but central to the dynamic
stopping and to sizing the Part-D cutoff. Cite Lang–Kohn image plane. [ref pack C.7]
</problem>
</part>

<part id="D" title="The long-range cutoff for the classical projectile">
<learning_objective>
Choose a defensible cutoff that removes periodic-image / loop-around self-interaction
while preserving the real screened + image physics — the prescription handed to
Campaign 1.
</learning_objective>

<problem id="D1" type="concept">
State the problem: the projectile's $-1/r$ tail interacts with its own and the slab's
periodic images, and a propagating projectile re-enters the opposite face. $G{=}0$ drop
does not fix this. [ref pack D.1; arXiv:2307.03213]
</problem>

<problem id="D2" type="constraint">
**DERIVE/STATE:** minimum-image bound $R_c<L_{\min}/2$ (=25 Bohr for the 50-Bohr in-plane
baseline, not set by $L_z$). [ref pack D.2]
</problem>

<problem id="D3" type="derive">
**DERIVE THIS:** Thomas–Fermi $k_{\mathrm{TF}}=\sqrt{4k_F/\pi}$, $\lambda_{\mathrm{TF}}=1/k_{\mathrm{TF}}$;
evaluate ≈1.5 Bohr for $r_s$=5.67. Conclude the cutoff window
$\lambda_{\mathrm{TF}}\ll R_c<L_{\min}/2$, e.g. $R_c\sim10$–20 Bohr. [ref pack D.3]
</problem>

<problem id="D4" type="compare">
Compare three prescriptions — (1) geometric (elongate $L_z$ + stop before re-crossing),
(2) truncated Coulomb kernel, (3) Martyna–Tuckerman nonperiodic Poisson — with
pros/cons. Note `inq/` is immutable ⇒ prescription 1 is lowest-risk for transit-only.
[ref pack D.4]
</problem>

<problem id="D5" type="ground">
From **arXiv:2307.03213**: the ~8% plasmon-cutoff finite-size error and the recommended
corrections; define the convergence test (vary $L_z$/$R_c$, energy reference and $S$
stable to a stated tolerance). [ref pack D.5; arXiv:1805.01377]
</problem>

<problem id="D6" type="deliverable">
Write the prescription handed to Campaign 1: scheme + values ($R_c$ or $L_z$+stop-time),
the two checks, a numeric finite-size tolerance, and the note that the static $t=0$
reference is image-insensitive to leading order while the dynamic stopping needs the
cutoff. [ref pack D.6]
</problem>
</part>

</worksheet>

<self_test>
Build a short self-test (answers from the reference-pack sanity table); the reader must
reproduce, to 2–3 s.f.:
1. $\langle T_{\mathrm{WP}}\rangle$ for σ_WP=0.5, 100 eV → 181.6 eV (and the zero-point 81.6 eV).
2. Bare self-Hartree $q^2/(2s\sqrt\pi)$ for σ_WP=0.5 → ≈22 eV; net SIE → ≈4.5 eV; explain the gap.
3. Why the classical $t=0$ excess is +798 eV, not ~100 eV.
4. The corrected energy reference $E_{\mathrm{jellium}}(0)=E_{\mathrm{tot}}(0)-\langle T_{\mathrm{WP}}\rangle-E_{\mathrm{SIE}}$.
5. The cutoff window for the baseline: $\lambda_{\mathrm{TF}}\approx1.5\ll R_c<25$ Bohr.
6. Predict $E_{\mathrm{SIE}}$ trend for σ_WP=3 (∝1/σ ⇒ <1 eV).
</self_test>

<references>
Pull the full list from the reference pack's References section (Perdew–Zunger 1981;
Lang–Kohn 1970/71/73; Martyna–Tuckerman 1999; Jackson; arXiv:2307.03213;
arXiv:1805.01377; in-repo `docs/sources/` notes and `localised-jellium-theory.md`).
</references>

<constraints>
- Reproduce in-repo numbers exactly; invent nothing. Label every non-sourced claim
  "Inference:".
- Atomic units; 2–3 s.f. for human-facing numbers; $\sigma\equiv\sigma_{\mathrm{WP}}$.
- This worksheet is the *continuation* of `localised-jellium-theory.md` — cross-ref, do
  not duplicate, Parts 0–8.
- `inq/` is immutable: any cutoff implementation lives in `inqkit`/`inqview` or is
  geometric.
</constraints>

<output_contract>
Deliver a single markdown (or the project's worksheet format) document titled
"Localised Jellium — Theory Worksheet, Parts B/C/D (SIE, Coulomb-vs-distance, cutoff)".
Suggested placement: append as new Parts to `docs/notes/localised-jellium-theory.md`
(continuing its Part numbering), OR a sibling `localised-jellium-theory-BCD.md` that
cross-links it. Must contain: per-Part learning objectives, the DERIVE-THIS problems
above with worked answers, the VC callouts, the self-test, and the reference list.
</output_contract>
