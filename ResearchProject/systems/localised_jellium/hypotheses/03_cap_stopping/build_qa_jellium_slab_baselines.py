#!/usr/bin/env python3
"""Assemble the jellium-slab baselines REVIEW notebook (the full Q&A record).

Path-referenced (markdown + relative image refs). Covers everything established in the
2026-06-23/24 grill-with-docs Q&A on the p5_wp / p5_classical CAP runs: geometry (Q2),
classical projectile construction (Q4), WP injection (Q5), absorption/transmission/
reflection (Q1/Q6), the stopping energy ledger (Q3/Q8), and the loss-function deferral
(Q7) — then the consolidated requirements for the NEXT simulation. Figures are produced
by qa_i..qa_v (run those first). Output: qa_jellium_slab_baselines.ipynb.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell

HERE = os.path.dirname(os.path.abspath(__file__))
cells = []


def md(s):
    cells.append(new_markdown_cell(s))


def fig(png, cap):
    cells.append(new_markdown_cell(f"*{cap}*\n\n![{cap}]({png})"))


# ===================== 0. Overview =====================
md(r"""# Jellium-slab baselines — full review
### localised jellium slab · 100 eV projectile · two-sided sin² CAP · WP vs classical

**Purpose.** A single review surface for the two proof-of-concept CAP runs, gathering
everything established in the Q&A so an informed choice can be made for the **next
campaign**. Every number is computed from run data — none invented. PROVISIONAL until
the inq-study engine regression (Task #7).

**The two runs** (`hypotheses/03_cap_stopping/`, scripts in `scripts/fullsuite_{wp,classical}/`):
- `p5_wp` — quantum **wavepacket** electron, σ_WP = 0.5 (**density std 0.354**), 100 eV.
- `p5_classical` — **classical** Gaussian test-charge ion, **charge std 0.350**, 100 eV.

**σ-convention (validated 2026-06-24, independent agent).** "Width" = wavepacket σ.
A WP of σ_WP has electron **density std σ_WP/√2**; the matched classical charge std is the
same. These runs are a **~1% mismatch** (WP 0.354 vs classical 0.350 — the legacy
`sigma0p35` UPF). The exact-matched UPF `electron_gaussian_wpsigma0p5.upf` (charge std
0.35355) is generated for future runs. See `CONTEXT.md` → "σ-convention".

**What this notebook covers:** §1 geometry (Q2) · §2 classical projectile (Q4) ·
§3 WP injection (Q5) · §4 absorption/transmission/reflection (Q1, Q6) · §5 stopping
energy ledger (Q3, Q8) · §6 loss function (Q7) · §7 next-simulation requirements.
""")

# ===================== 1. Geometry (Q2) =====================
md(r"""## §1 — Geometry: 2D-periodic slab, not 3D bulk (Q2)

**As built** (verified from `localised_background.hpp` + `run.cpp`): a 50³ Bohr cubic
**periodic** cell; an explicit positive jellium background confined to the slab
`|z| < 12.5` (r_s ≈ 4, N = 234); two-sided sin² CAP (η = −0.5) in `[17.5, 25]` and
`[−25, −17.5]`.

- **Periodic & seamless in x, y** — the background depends on **z only** (full 50×50
  face), so periodic images tile with no gap: a genuine **infinite 2D-periodic slab** in
  the xy-plane.
- **z is formally periodic too**, but the background is localised with vacuum padding
  (slab-in-a-supercell). So this is **not** the repo's old whole-cell **periodic-3D bulk**
  jellium — the slab has **two real surfaces** (spill-out, Friedel oscillations, an image
  potential) that bulk does not.

**What generalises / caveats.** Interior stopping should approach bulk if the slab is
thick enough; **surface** physics lives near the faces. Finite-size watch-outs: slab
**z-images are only 50 Bohr apart** (dipole/image coupling through the vacuum), and the
slab is only **25 Bohr thick**. A true bulk comparison would need the periodic-3D run or a
thicker slab + dipole correction.
""")

# ===================== 2. Classical projectile (Q4) =====================
md(r"""## §2 — How the classical projectile is built (Q4)

`run.cpp:80`: `ionic::species("H").pseudo_file(electron_gaussian_sigma0p35.upf).mass(1 mₑ)`,
inserted at z = −15.5, velocity (0,0,2.71), Ehrenfest.

**The two "changes to the ion":**
1. **UPF `z_valence` 1.00 → 0.00** (the file carries this comment). The projectile is a
   **chargeless ghost**: pure moving Gaussian potential (erf-smeared −1 test charge,
   repulsive to electrons), adding **no electron, no nuclear charge to neutrality, no SCF
   seed density**.
2. **Mass = 1 electron mass** (not heavy/fictitious) → free Ehrenfest, so it **decelerates
   measurably** (2.71 → 2.35) and `ΔKE_ion` is a real stopping signal.

`extra_electrons(234)` is the **bath, identical in both runs** (it is large only because
the *slab* is N=234, not because of the projectile). Charge bookkeeping:

| | electrons | neutrality | projectile |
|---|---|---|---|
| **classical** | 234 (bath) | neutral | Z=0 ghost **potential** (test charge) |
| **WP** | 235 (234 bath + 1 WP) | **net −1** (G=0-compensated) | a **real electron** |

**Impact.** The classical projectile is a **test-charge potential** — it cannot be
absorbed (it transmits; §4), and it is not in the electron density (so it is
**reconstructed** for analysis). The WP is a **real electron** — it feels **Pauli
exchange** and its own **self-interaction (SIE, ~7–21 eV)**, neither of which the classical
has. So part of any WP−classical difference is **SIE artifact, not quantum physics** — the
single biggest thing to control next.
""")

# ===================== 3. WP injection (Q5) =====================
md(r"""## §3 — Was the WP absorbed at injection? No (Q5)

`N_wp(0) = 1.000` exactly (density_wp integrates to unity at t=0) and `wp_norm_after = 1`.

- **Not captured by the CAP at t=0:** launch z=−15.5 is **2 Bohr (5.6σ) from the −z CAP
  edge** (density std 0.354) → Gaussian tail ~10⁻⁷, negligible.
- **Orthonormalisation:** full norm retained after `orthogonalise_against_occupied`.

So the WP is **fully present at injection** — no loss to the CAP or orthonormalisation.
Minor caveat: launch is only **3 Bohr from the slab face (−12.5)**, so it begins
interacting with the slab spill-out early (fine here; tune for a cleaner pre-slab baseline
next time). The norm overview below shows N_wp(0)=1 and the drain that follows.
""")
fig("qa_iv_norms.png", "Total / bath / WP norms vs time (WP norm on the right axis); N_wp(0)=1")

# ===================== 4. Absorption / transmission / reflection (Q1, Q6) =====================
md(r"""## §4 — Absorption, transmission & reflection (Q1, Q6)

**Regions (z, Bohr):** −z CAP [−25,−17.5] · left-free [−17.5,−12.5] · slab [−12.5,12.5] ·
right-free [12.5,17.5] · +z CAP [17.5,25]. WP density from `density_wp`; bath = total − wp;
WP current from the complex `wavefunction_wp` (**spectral** z-derivative — finite-difference
undercounts the k₀=2.71 current ~28% at dx=0.5).
""")

md(r"""### 4.1 How many electrons did the CAP absorb? — robust total-density count first

**Decomposition-free measure (the safe one).** The number of electrons the CAP removed is
simply the drop in the *total* electron norm, `N_total(0) − N_total(t)` with
`N_total = ∫_box n_total dV` — **no slab/WP separation required**. The baselines are
integer-clean, so this is unambiguous:
- **WP run:** `235.000 → 234.169` ⟹ **0.831 absorbed**.
- **classical run:** `234.000 → 233.784` ⟹ **0.216 absorbed** — and this *is* the bath
  overflow, because the classical projectile is a ghost **potential**, not in the density.

**Decomposing the WP run — with the caveat.** ⚠️ The WP Kohn–Sham orbital is **not a clean
physical separator**: it hybridises with the bath and the KS orbitals stop being good
descriptors, so **"0.62 of the wavepacket was captured" is _not_ a physical statement.** What
*is* well-defined is the **WP KS-orbital norm** `∫_box |ψ_WP|² dV`, which falls `1.000 → 0.378`
(**−0.622**), while the **bath** overflows **−0.209**. These add back to the robust total:
`0.622 + 0.209 = 0.831` (matches the total-density count to machine precision). Both baths
overflow **~identically** (WP-run bath 0.209 vs classical 0.216) — a **CAP property, not a
projectile effect**. So the one clean WP-vs-classical contrast is that the WP run loses an
**extra 0.622 of WP KS-orbital norm** that the classical (ghost, transmitting) run does not.
**Not equilibrated** at t=18 (KS-orbital only 62% drained).""")
fig("qa_iii_absorbed_norm.png", "Cumulative absorbed norm — WP KS-orbital 0.62 / bath 0.21 / classical bath 0.22")

md(r"""### 4.2 Where the projectile is — region occupancy
WP drains slab → right-free → +z CAP (forward) with a small back-trickle to −z bands;
classical marches ballistically and periodic-wraps to z=−20.""")
fig("qa_i_region_densities.png", "Projectile charge per z-band vs time — WP vs classical")

md(r"""### 4.3 Bath redistribution / drainage (Δ from t=0)
The t=0 baselines are large (~9.7 e bath spill-out between slab and CAPs — the
boundary-width density). Plotted as Δ-from-t0 to expose the dynamics.""")
fig("qa_iv_bath_bands.png", "Bath density change per z-band (Δ from t=0) — WP run")

md(r"""### 4.4 Reflection diagnostic — un-absorbed density between slab & CAPs  ⚠️
After left-free empties (min 0.015) it **refills to ~0.092** (peak t≈9) — density heading
back toward the −z CAP. **This is consistent with reflection BUT ALSO with the trailing
edge of the ~23×-spread WP**; 1D density cannot separate "moving −z (reflection)" from
"spread tail." Beyond t≈14.9 transmitted WP **periodic-wraps**. **Reflection is therefore
NOT cleanly resolved** — do not read it as either negligible or large.""")
fig("qa_iv_reflection_freeregion.png", "Un-absorbed WP between slab & CAPs (full + zoom)")

md(r"""### 4.5 Per-CAP flux — transmission vs reflection (validated method, contaminated split)
T+R = 0.698 vs the independent bookkeeping 0.700 → **method validated**. Forward
(transmission) flux dominates (~0.743); the −z-edge flux comes out net slightly forward —
transmitted WP **wrapping around** the periodic box — so genuine reflection is **masked**.""")
fig("qa_ii_per_cap_flux.png", "Per-CAP cumulative flux: transmission (+z) vs reflection (−z)")

md(r"""### 4.6 The WP centroid stalls — the packet is non-rigid  ⚠️
The WP does **not** translate as a rigid packet, and this **invalidates the naïve
"centroid ± 3σ crosses each face" timing scheme**. From `wp_real_space_stats`:

- Its z-centroid `⟨z⟩(t)` (survival-weighted) rises from the launch point **−15.5**,
  decelerates, and **stalls at ⟨z⟩ ≈ +5.1** (max +5.6) from t≈9 on — it **never reaches the
  far slab face (+12.5)**. The centroid velocity `d⟨z⟩/dt` collapses **2.64 → −1.28** (it even
  drifts *backward* late).
- Meanwhile `σ_z` balloons **0.37 → 11.7 (≈32×)**, so the **±3σ envelope covers the whole
  50-Bohr box** and the "trailing edge" `⟨z⟩−3σ_z` actually runs *backward* to ≈ −30. Only the
  **leading edge** `⟨z⟩+3σ_z` cleanly crosses the faces (near −12.5 at t≈0.4, far +12.5 at
  t≈4.0); the centroid- and trailing-edge crossings of the **far face are undefined**.

**Why this matters physically.** The stall is a **survival-weighted absorption signature**,
*not* the packet physically stopping: the +z CAP removes the fast forward components
(`pz_mean` 2.63→1.71, `e_kin` 6.65→4.31 Ha), so the surviving WP is the slow / back-scattered
remainder and `⟨z⟩` is dragged back. The classical projectile (overlaid) marches
**ballistically** through the slab to the +z box edge, decelerating only **2.71→2.35** — *that*
deceleration is the real, interpretable stopping signal (the classical curve is mainly a
reference here, as expected). This non-rigidity is the same σ=0.5 zero-point/spreading problem
that makes the energy-ledger stopping unmeasurable (§5), and it is what motivates the larger-σ,
near-rigid WP for the next campaign (§7 #1).""")
fig("qa_viii_centroid_trajectory.png", "z-centroid vs time — WP stalls at +5 (non-rigid, survival-weighted); classical marches ballistically")

md(r"""### 4.7 Cumulative current across every internal plane — reflection is negligible *where observable*
Generalising §4.5 to **all four internal planes**: signed WP current `J_z` (spectral derivative),
integrated over each plane and cumulated in time, with the §4.6 kinematic events annotated.
Sign: **+ = +z (forward), − = −z (backward)**. Transport budget (cumulative, at wrap onset
t≈14.9 | final t=18):

| plane | cum ∫J_z dt @t≈14.9 | final | reading |
|---|---|---|---|
| near face (−12.5) | **+1.01** | +1.03 | ~1 e of WP passes **forward** into the slab (backward peak only −0.008) |
| far face (+12.5) | **+0.72** | +0.81 | forward transmission out of the slab |
| +z CAP edge (+17.5) | **+0.63** | +0.74 | forward absorption (matches §4.5 transmission 0.743) |
| **−z CAP edge (−17.5)** | **+0.033** | +0.045 | **≈ zero, and net *forward*** — not backward |

**Reading (reflection).** In the wrap-free window the backward-side planes show **no significant
backward current**: the −z CAP edge integrates to **+0.033 e forward** (not backward) and the
largest instantaneous backward excursion anywhere is **−0.008**. So the left-free density refill
(§4.4) is **more consistent with the spreading tail than with a large reflected wave** — genuine
**near-face** reflection is **negligible** here.

⚠️ **Caveat — far-face reflection is unmeasured, not absent.** A wave reflected at the **far face
(+12.5)** would only reach the −z CAP at **t ≈ 22 a.u.**, *after* this 18-a.u. run ends. So
"reflection negligible" holds **only for what is observable in this run**; the far-face channel is
simply not yet probed. This is precisely why the next campaign needs a **longer run / box** (§7
#2, #3) before any reflection verdict — and a **sign-resolved n(k_z)** (§7 #4) to see a backward
lobe directly. (Frame cadence Δt = 0.2 a.u.; a sub-0.2-a.u. fast-component pulse could be
undersampled — flagged, not hidden.)""")
fig("qa_ix_cumulative_current_regions.png", "Cumulative WP current across all internal planes (+forward/−backward) + backward-side flux (wrap-free)")

# ===================== 5. Stopping energy ledger (Q3, Q8) =====================
md(r"""## §5 — Stopping power from total energy, and why σ=0.5 fails it (Q3, Q8)

**The energies (todo 1).**

| state | E_total (Ha) | E_total (eV) | above slab GS (eV) |
|---|---|---|---|
| GS jellium slab (no projectile) | −160.9921 | −4380.8 | 0 (reference) |
| WP run, t=0 | −154.1618 | −4194.3 | **+185.9** |
| classical run, t=0 | −131.6538 | −3582.1 | **+798.3** |

The WP run sits +185.9 eV above the bath GS (≈ 100 drift + ~82 zero-point + ~4 net
self/interaction). The classical t=0 sits **+798 eV** above the bath GS — the ghost
potential applied to the *unrelaxed* loaded GS density (a sudden, non-variational jump),
**not** a physical 100 eV. This is why the naive Q8 check (WP−classical ≈ 100 eV) fails.

**Method (agreed):** work with **total energy only** (WP and bath cannot be separated
energetically). Deposited energy two ways, S = deposited / 25 Bohr:

| baseline | deposited | S |
|---|---|---|
| **Formula 1** — `E_total(final) − E_total(0) + 100 eV` (subtract the drift) | **32.2 eV** | 1.29 eV/Bohr |
| **Formula 2** — `E_total(final) − E_GS_slab` (−160.99 Ha; exact at full absorption) | **118.1 eV** | 4.72 eV/Bohr |

**Why Formula 2 is "outrageously large" — corrected (todo 7).** Formula 2 is **not** the
deposited stopping at t=18: `E_total(final) − E_GS = (energy deposited in the bath) +
(energy of the residual un-absorbed WP still in the box)`. The WP **KS orbital** is only
**62% absorbed**, so the residual ~0.38 e of WP-orbital norm carries roughly **~70 eV** of WP
energy that *has not yet left* —
and **that** is what inflates Formula 2 to 118 eV. Formula 2 equals the deposited energy
**only at full absorption** (residual → 0). Formula 1 (32 eV) fails differently: it assumes
the WP brought exactly 100 eV, but the table shows it brought **185.9 eV** — the extra ~86 eV
is the WP's zero-point KE (`3/(4σ²)=81.6 eV` for σ=0.5) plus net self-interaction. So **both
formulas are contaminated** (Formula 1 by the wrong baseline, Formula 2 by incomplete
absorption) and **the σ=0.5 WP stopping is not measurable from this run.** The **classical**
number stands: **S = 0.706 eV/Bohr ≈ point-charge Lindhard 0.719**.

**Reading the figure (todo 9).** `E_total(t) − E_total(0)` ends at **−68 eV** — that is the
energy the CAP *removed* (with the absorbed WP), **not** the WP's energy. Adding back the
100 eV drift (dashed curve) gives the **slab-remaining** energy change = Formula 1, ending
at +32 eV.
""")

md(r"""**Reader observation (verbatim) — the initial energies are not comparable.**
> - Here, the starting energies of the classical and the wave packet runs are not comparable. Meaning, the hartree term might be playing an impact.  the wavepacket and the projectile must be initialised at a distance such that the Hartree term is almost the same, and E_wavepacket_total is approximately 100 eV higher. Now, I should expect self interaction to make the energy difference more than 100 eV, we should clearly do this for the next set of simulations to ensure the results are comparable.

**Response.** Agreed, and consistent with the table above. The classical t=0 sits **+798 eV**
above the bath GS because the ghost potential is applied to the *unrelaxed* loaded density (a
sudden, non-variational jump), so `E_WP(0) − E_classical(0)` is **not** a meaningful ~100 eV —
exactly your point. The fix you state is the right one: launch the WP **and** the classical
projectile **far enough from the slab that the electron–slab Hartree term is essentially
identical** in both, so the only intended t=0 difference is the WP's drift + zero-point. As you
note, the WP's **self-interaction** then makes `E_WP(0) − E_classical(0)` exceed 100 eV — which
is precisely why a **vacuum-WP SIE control** is needed to quantify the SIE before any "quantum
component" is reported. Folded into the next-sim requirements (matched-Hartree launch + SIE
control, §7 #6).""")

fig("qa_v_stopping_energy.png", "Total-energy change vs time — WP, WP+100 eV (slab-remaining), classical")

md(r"""### 5.1 Quantifying the wavepacket self-interaction (SIE) — with the *corrected* reference
**Question.** "Just by *having* the wavepacket, how much extra energy is added to the system —
the spurious one-electron self-interaction (SIE)?" In LDA a single electron feels a spurious
Hartree self-repulsion that XC does **not** fully cancel (one-electron SIE; Nazarov & Gross 2025
and the SIE literature). To isolate it, launch the WP **far from the slab** (WP–slab Hartree
cross-term → 0) and compute

`SIE = E_total(0)_WP[far] − E_GS_slab − KE_WP`

⚠️ **Correction to the proposed test.** `KE_WP` must be the WP's **measured total kinetic
energy** `⟨p²⟩/2` (from `wp_momentum_stats`), **not "100 eV"**. The 100 eV is only the *drift*
KE; the WP also carries ~81 eV of **zero-point + transverse** KE. Using "+100 eV" would omit
that 81 eV and **overcount the SIE by ~81 eV** (reporting ~85 eV instead of ~4.5 eV). This is
the single most important correction here.

**Result — no new run needed.** The Phase-3 baseline **`p3_wp` IS the far-launch run**
(`LJ_CAP=0` defaults `launch_z = −23`, so the WP sits 10.5 Bohr from the slab face, no CAP):

| run | launch z | dist → face | E_total(0) (Ha) | KE_WP (Ha / eV) | excess = E_tot(0)−E_GS−KE_WP |
|---|---|---|---|---|---|
| **p3_wp (far)** | −23.0 | 10.5 Bohr | −154.1797 | 6.6451 / 180.8 | **+0.167 Ha = +4.55 eV** |
| p5_wp (near) | −15.5 | 3.0 Bohr | −154.1618 | 6.6457 / 180.8 | +0.185 Ha = +5.02 eV |

- **SIE ≈ 4.5 eV** (single-electron LDA self-interaction of the σ=0.5 WP). The WP–slab
  **Hartree term is small**: moving 7.5 Bohr farther changed the excess by only **0.47 eV**, so
  the excess is dominated by genuine SIE, not Hartree.
- **KE_WP = 180.8 eV** (≈93 drift + ~88 zero-point+transverse) — confirms the "+100 eV"
  reference is wrong by ~81 eV.
- ⚠️ **σ-dependent:** the Gaussian Hartree self-energy scales ~1/s (s = density std), so a
  larger-σ **rigid** WP has a much smaller SIE (σ_WP=3 → ~6× smaller, < 1 eV). The SIE must be
  **re-measured per σ** in the campaign. A true **vacuum-WP control** (no slab) would remove the
  residual ~0.5 eV image/Hartree term for a Hartree-free number (§7 #6).""")

# ===================== 6. Loss function (Q7) =====================
md(r"""## §6 — Loss function L(q,ω): deferred (Q7)

Two distinct objects: the **analytical** RPA `L(q,ω) = Im[−1/ε]` (ungated theory) and the
**run-extracted** `|n_q(ω)|²/q²` (Fourier-from-run-data, **hard-gated** behind Fourier
training). **Both are deferred** (user decision) until that training. Note the resolution
wall regardless: at τ ≈ 18 a.u. the FFT bin is **Δω ≈ 9 eV**, coarser than the ~6 eV
r_s≈4 plasmon — a clean run-extracted L(q,ω) needs a **~10³ a.u.** run.
""")

# ===================== 7. Next-simulation requirements =====================
md(r"""## §7 — Requirements for the next simulation

**Robust conclusions to build on.** The CAP-absorbed electron count is the total-density drop
`N_total(0)−N_total(t)` = **0.831** (WP) / **0.216** (classical); bath overflow is
run-independent (~0.21 both); the WP run loses an extra ~0.62 of **WP KS-orbital norm**;
transmission dominates (~0.74); the classical stopping (0.706 eV/Bohr)
matches point-charge Lindhard; **the σ=0.5 WP zero-point energy (82 eV) makes its E_total
stopping unmeasurable**.

**Open / unresolved.** Reflection extent (reflection vs spreading indistinguishable in 1D);
the WP stopping number; the SIE contamination of WP−classical.

**The σ–energy spreading tradeoff (the central difficulty).** A free Gaussian disperses as
`σ_WP(t)=σ_WP√(1+(t/σ_WP²)²)`, so the fractional spread over the transit (length X, time
t*=X/√(2E)) is `f(σ,E)=√(1+(t*/σ²)²)−1`; inverting at fixed f gives the minimum energy for a
width, `E_min(σ,f)=X²/(2σ⁴((1+f)²−1))` (X≈28 Bohr, launch→far edge). **LEFT** spread% vs σ_WP at
E=100/300/600 eV; **RIGHT** E_min vs σ_WP for several spread thresholds. This is why σ=0.5 fails
and why (σ, E) must be chosen **jointly**: at fixed E spreading falls steeply with σ; at fixed σ
higher E (shorter transit) reduces spreading. The 20%-spread contour is the practical boundary.

![σ–energy spreading tradeoff — the central difficulty](ref_spread_tradeoff.png)

**Design inputs for the next campaign (to settle these):**
1. **Larger σ (≈3)** — zero-point KE drops to ~2 eV (`3/(4·9)`), the WP becomes near-rigid
   (no 23× spread), and the E_total stopping becomes clean. The CONTEXT "matched-pair"
   direction (σ-scan extrapolated to s₀→0).
2. **Run to full absorption** (≥30 a.u.) so Formula 2 (`E_total − E_GS`) converges and the
   residual-WP contamination vanishes.
3. **Longer z-box / stop before t≈15 a.u.** — no periodic wrap, so reflection and the
   per-CAP split are clean.
4. **Momentum-resolved observable (k_z sign / 2D k_z–k_⊥)** — separate reflection (k_z<0)
   from forward spreading.
5. **Fine-cadence flux/current observable** (spectral-derivative current or a CAP-edge flux
   screen) — exact per-CAP split, incl. the classical bath.
6. **Exact-matched UPF** (`electron_gaussian_wpsigma0p5.upf`-style) + a **matched-Hartree
   launch** (WP and projectile started far enough from the slab that the electron–slab Hartree
   term is identical, so `E_WP(0)−E_classical(0)` is only drift+zero-point+SIE) + a **vacuum-WP
   control** to bound the SIE before any "quantum component" is reported.
7. **CAP L/η tuning** from the two-sided ε(E,L) reflection maps, informed by 3–5.

**Momentum-snapshot anchor scheme (run-independent — for the next runs' k_z–k_⊥ observable).**
Because the next runs change the **jellium density and the WP energy**, fixed snapshot *times*
do not transfer — but the **physical events** do. Snapshot the full 3-D momentum distribution
`ψ(r)→ψ(k)` (FFT of the complex `wavefunction_wp`) at events computed *per run* from `⟨z⟩(t)`
and `σ_z(t)`:

| anchor | event | note |
|---|---|---|
| **A0** | launch (t=0) | reference forward lobe at k_z = k₀ |
| **A1** | leading edge `⟨z⟩+3σ_z` reaches near face | first slab contact |
| **A2** | centroid `⟨z⟩` reaches near face | slab entry |
| **A3** | packet fully inside slab | early scattering (convenient mid-point) |
| **A4** | leading edge reaches far face | first transmission opportunity |
| **A5** | centroid at max forward excursion (`d⟨z⟩/dt→0`) | turnaround / stall onset |
| **A6** | `σ_z` saturates / centroid stalls | scattering "complete" |
| **A7** | quasi-steady, **before** any periodic wrap (`t<t_wrap`) | clean (no wrap contamination) |
| **A8** | final state | — |

For *these* baseline runs the anchors land at t ≈ 0, 0.4, 1.2, 2.4, 4.0, 6.0, 9.0, 13.0, 18.0 —
but it is the **scheme, not the times**, that carries forward (every anchor shifts with density
and energy). Diagnostics at each anchor: **sign-resolved n(k_z)** (reflection = a k_z<0 lobe),
**n(k_⊥)** (transverse heating), the **2-D (k_z, k_⊥) map**, and the **survival ratio
n_final(k_z)/n_init(k_z)** — which directly tests whether **high-|k_z| components are absorbed
preferentially** (the "fast components captured / reflected first" hypothesis). Caveat from §4.6:
for a spreading σ=0.5 WP the centroid/trailing anchors (A2, A5, A6) are survival-weighted and
the trailing-edge ones can be undefined — they sharpen once the next campaign uses a larger-σ,
near-rigid WP (#1).

*All figures path-referenced beside this notebook; regenerate via `qa_*.py` then
`build_qa_jellium_slab_baselines.py`.*
""")

# ===================== 8. Round-2 reader observations (verbatim) =====================
md(r"""## §8 — Round-2 reader observations & validations (verbatim)

Observations and TODOs raised on review, **quoted verbatim**, each with a response.
Items 1, 7, 9 are addressed inline in §5 (energies table, Formula-2 correction, +100 eV
curve).""")

md(r"""### Todo 2 — how was the WP's "62% absorbed" defined?
> One important observation I have is that the bath overflow for the wavepacket run and classical run are the same. But visually, the wavepacket run looks like its pushing the density from a much bigger range, however, classsical pushes the density right ahead of it. How was it calculated that teh WP orbtial's 62% density was absorbed. Meaning, was this from the KS orbital? (remind me of how this was defined - was this defined by integrating over the bath region of the simulation cell at final time).

**Response.** It is from the **WP Kohn–Sham orbital's own density**, `density_wp = |ψ_WP|²`,
integrated over the **entire** simulation cell (NOT the bath region):
`N_wp(t) = ∫_box |ψ_WP(r,t)|² dV`. It went **1.000 (t=0) → 0.378 (t=18)**, so
`1 − 0.378 = 0.622 = 62%` is the WP orbital norm the CAP removed from the whole box. It is
**not** an integral over the bath region — it is the WP orbital's total surviving norm.
**Important caveat:** because the WP and bath are not cleanly separable (the KS orbitals stop
being good descriptors as the WP hybridises with the bath), this 62% is a statement about the
**WP KS-orbital norm**, *not* "62% of the physical wavepacket." The **safe, decomposition-free**
number of absorbed electrons is instead the **total-density** drop `N_total(0) − N_total(t)`
(§4.1): **0.831** (WP) / **0.216** (classical) — that one needs no slab/WP separation and is the
figure to quote.
Your visual observation is **correct and is the spreading difference**: the WP density std
grows ~23× (0.354 → tens of Bohr), so it pushes density over a wide range; the classical
Gaussian stays compact (std 0.35, no spreading), pushing right ahead. The final-time xz map
(Todo 10 below) shows this directly. **Validated.**""")

md(r"""### Todo 3 — total norm of the simulation, easy to compare (+ classical projectile)
> The classical overlays in the projectile charge in band make sense. Now, thinking about the total density of electrons in and out, make a plot of total norm of the simulation of the classical and wavepacket runs and make it easy to compare. Also, you need to add in 1 density for the classical run. You also need to validate the statements I want you to validate. Write my words verbatim in markdown cells of the ipynb.

**Response.** Below: both runs plotted starting at **235 charge units** — the classical run's
reconstructed **1-unit projectile is added** to its electron count. WP loses **0.831**
(WP **KS-orbital** 0.62 + bath 0.21); classical loses only **0.216** (bath only; the projectile
is conserved at 1.0 and transmits). **Statements validated:** the bath overflow is the same
in both runs (0.209 vs 0.216), and the entire difference is the WP **KS-orbital norm** being
absorbed. (This 0.831 is exactly the robust `N_total(0)−N_total(t)` total-density count of §4.1
— it needs no slab/WP separation.)""")
fig("qa_vi_total_norm_compare.png", "Total norm of the simulation — WP vs classical (both start at 235)")

md(r"""### Todo 4 — what is t0 in `n(t) − n(0)`?
> what is t0 in the bath density change n(t) - n(0) run?

**Response.** `t0` is the **first RT frame, t = 0 a.u.** — the state at the very start of the
propagation (GS loaded + background + projectile injected, step 0), before any dynamics. For
the bath specifically, `n_bath(0) = n_total(0) − n_wp(0)`. So `n(t) − n(0)` is the change
relative to that initial frame, which **already contains** the injected projectile and the
large spill-out baseline (≈9.7 e between slab and CAPs).""")

md(r"""### Todo 5 — can the inter-region density picture genuine reflection?
> The unabsorbed density between CAP and slab is an important measure. However, I do not know how I can turn this to check if there is genuine reflection. I am currently thinking while typing. If there is a reflection, then there would be essentially being moving towards the simulation box. Can I use this somehow to picture reflection? Because now it might be too vague to answer this question.

**Response.** Your intuition is exactly the right definition: **reflection = density (or
current) moving in the −z direction**, back toward the box edge it came from. The clean way
to *picture* it is a **sign-resolved momentum distribution** `n(k_z)` — reflection appears as
a **k_z < 0 lobe** (backward-moving), distinct from forward spreading (k_z > 0 spread around
k₀). These runs save only **1D |k|** (no sign), and the −z-edge probability current (§4.5) is
contaminated by periodic wrap-around — so with *this* data it stays vague, as you say. This is
exactly **next-sim requirement #4** (k_z-resolved momentum / directional current) **+ #3**
(longer box, no wrap): with those, your reflection picture becomes one clean plot.""")

md(r"""### Todo 6 — what does "WP transmission" mean? whole box or the WP orbital?
> In figure 6 and other figures in general, what does WP transmission mean? Does it mean that the transmission was calculated using the entire simulation box or only the wavepacket KS orbital.

**Response.** **Only the wavepacket KS orbital.** "WP transmission" (§4.5 / `qa_ii`) is the
probability current `J_z = Im(ψ_WP* ∂_z ψ_WP)` of the *complex `wavefunction_wp`* orbital,
integrated across the +17.5 plane — **not** the whole-box (bath + WP) current. The bath's own
flux is not separately available from these runs (it needs the current-density field, a
next-sim observable).""")

md(r"""### Todo 10 — final-time xz log density with region demarcations
> At the last time instant, make a log plot of total density of the box using a xz plot, clearly using dashed lined to demarkate different regions, CAP, bw CAP and jellium-slab, and jellium-slab.

**Response.** Below: final-frame (t=18) total electron density, mid-y xz slice, **log** scale.
**Cyan dashed = jellium-slab faces (±12.5); lime dashed = CAP inner edges (±17.5); white
dotted = box edges (±25).** The spread WP (wide, low-density across many bands) vs the compact
classical projectile is directly visible — confirming the Todo-2 observation.""")
fig("qa_vii_xz_logdensity_final.png", "Final-time total density (xz, log) — cyan=slab faces, lime=CAP edges")

nb = new_notebook()
nb.cells = cells
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = os.path.join(HERE, "qa_jellium_slab_baselines.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
