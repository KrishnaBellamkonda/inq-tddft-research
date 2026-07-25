In this document, I am going to write down my understanding of the muon implementation in INQ. 

To really understand the plan, I need to understand how the INQ calls the scripts in a chain. 
1. The first call comes from the run.cpp file where the actual system is defined. Here, the etrs propagator is called. 

> 💬 **FEEDBACK (mostly right, one precision):** run.cpp calls `real_time::propagate(...)`, NOT `etrs` directly. `propagate()` owns the time-loop and, *each step*, dispatches to EITHER `etrs` OR `crank_nicolson` depending on `opts.propagator()`. So the real chain is: `run.cpp → real_time::propagate → (per step) etrs / crank_nicolson`. This matters for us because `propagate()` is where the Hamiltonian is *built* (and where we inject `inverse_mass`) — see next point.

2. In etrs propagator, because of the structure of this propagator, it calls the exponential function. IN the exponent, there is the KS hamiltonain that needs to be built. 

> 💬 **FEEDBACK (important correction):** the KS Hamiltonian is **built ONCE in `propagate()` BEFORE the loop** (`ks_hamiltonian ham(...)` at propagate.hpp:79), then *passed into* `etrs`/`crank_nicolson`, which pass it to `operations::exponential_*`. The exponential does NOT build the Hamiltonian — it **applies** `ham(phi)` repeatedly to construct the Taylor series of `exp(-iH·dt)`. Each step the potential is refreshed (`sc.update_hamiltonian`), but the object (and its `inverse_mass_`) is constructed once. 👉 This is exactly *why* the mass must live on the Hamiltonian and be keyed by state-index: the exponential re-applies H many times per step, and the muon mass must be identical on every one of those applications.
3. The KS hamiltonian is built in the ks_hamiltonian.hpp file. KS hamiltonian has the kinetic, coulombic and the xc terms. The kinetic term is evaluated in fourier space using the laplacian function. The coulombic repulsion (i presummed is also worked out in fourier space, and then brought back to real space) is applied in real space. Finally, the density of the system the defined the E_xc[n]. So, the total KS Hamiltonian is built here. 

> 💬 **FEEDBACK (correct, with one clarification about *where*):**
> - **Kinetic in Fourier via `laplacian`** — ✅ correct.
> - **Hartree "worked out in Fourier, applied in real"** — ✅ correct and nicely put: the Hartree potential is solved via the Poisson equation `V_H(G) = 4π n(G)/G²` (diagonal in G-space), then FFT'd to real space and applied as `V_H(r)·ψ(r)`.
> - **Clarification on *where* this happens:** these potentials are NOT re-derived inside `ks_hamiltonian::operator()`. They are pre-assembled ONCE per step into a single combined real-space field `scalar_potential_` (= local-pseudopotential + Hartree + XC) by `hamiltonian::self_consistency::update_hamiltonian`. `operator()` then just does ONE real-space multiply (`scalar_potential_add`) plus the kinetic (`laplacian`) plus the non-local projectors. So "the Hamiltonian is built here" is true for the *kinetic + apply* machinery; the *potential assembly* (Hartree/XC/Poisson) lives in `self_consistency`.
> - **`E_xc[n]` from the density** — ✅ correct; V_xc is folded into that same `scalar_potential_`.
> 
> ⭐ **Takeaway:** the muon fork touches **only** the kinetic (Fourier/`laplacian`) branch. Everything in the real-space `scalar_potential_` (Hartree + XC + local) is untouched by mass — which is exactly your next paragraph's insight.

Now, going to the suggested plan. In the suggested plan, the idea is to create add an attribute to the electrons class that stores all the states in it called as inverse_mass. I believe this would be an array (need to design this carefully to ensure that each kpin and each ist in each kpin gets an inverse mass correctly). This arrays provides a mass to each of the states, and hence we control the dial this way. To do this, we change the electrons.hpp file. Then, we change the laplacian file (where the kinetic operator is applied, as this is the most affected in due to the mass). The coulombic repulsion energy I would assume is the same for a muon and an electron with the target system. The E_xc would be different for a muon, but for a second, we are going to ignore this. Other files such as the ks_hamiltonian and the propagate files are changed so that the right arguments are passed to the right functions. 

> 💬 **FEEDBACK (this paragraph is essentially correct — good):**
> - **`inverse_mass` as a `[kpin][ist]` array, mirroring `occupations_`** — ✅ exactly right, and your flag "design carefully so each kpin/ist gets its inverse mass correctly" is the single most important correctness point. The safe rule: the array is stored **per-rank local** (like `occupations_`), and in the GPU kernel it is indexed by the **same local `ist`** the field-set uses — so it stays aligned with the state partition even under MPI. (This alignment is also why we store it on the Hamiltonian, not on the transient `orbital_set` temporaries.)
> - **Change `electrons.hpp` + `laplacian.hpp`** — ✅ the two core files (data home + compute).
> - **"Coulomb repulsion is the same for muon and electron"** — ✅✅ **exactly right, and this is the deep reason the fork is small.** Hartree depends only on **charge** (μ⁻ and e⁻ both have −e) and density, NOT mass. So `V_H` is byte-for-byte unchanged; only the kinetic prefactor changes.
> - **"E_xc would be different for a muon; ignore for now"** — ⚠️ **right to defer, but here's the precise nuance so you defer it knowingly:**
>   - For **muon-in-ELECTRON-jellium** (heavy projectile, distinguishable): the muon genuinely should NOT be inside the electrons' XC (it's a different species → no cross-species exchange, no XC self-energy). Lumping it in is the SIE/cross-species error → we *bound* it with a vacuum control (option a) or give the muon its own spin channel (option b). "Ignore for now" = "accept + bound later." Fine.
>   - For **all-MUON jellium at fixed r_s**: E_xc is actually **NOT different** — a muon HEG at the same dimensionless r_s IS the electron HEG rescaled, so electron-LDA is *exact* there (see engine notes §3 rescaling fact). So "ignore" is not even an approximation in that case.
> - **`ks_hamiltonian.hpp` + `propagate.hpp` for argument passing** — ✅ correct.

Quetions:
1. Would ground state be affected? If so, how does this workflow ensure that the GS is converged correctly given a muon jellium system instead of an electron jellium? 

> 💬 **FEEDBACK / ANSWER (excellent question — it splits the whole design in two):**
> **It depends on WHICH system, and this is the key realisation.**
>
> **Case A — muon (heavy) projectile through an ELECTRON jellium.** The GS is a pure *electron* jellium; the muon does not exist yet at GS time — it is *injected as a wavepacket at t=0 of the real-time run*. So the **GS is NOT affected**: it converges exactly as today with all masses = 1. This is why my original design put the mass opt-in *only* in `propagate.hpp` and left the ground-state constructors (`ground_state/calculator.hpp`, `initial_guess.hpp`) untouched — the GS is provably identical to the unforked engine.
>
> **Case B — all-MUON jellium (the bath itself is muon-mass), and the band-structure campaign.** Here the **GS IS affected** — the *bath particles* carry muon mass, so the SCF ground state must be solved with `-1/(2 m_μ)∇²`. That means `inverse_mass` must ALSO be passed into the ground-state construction sites (`ground_state/calculator.hpp` + `initial_guess.hpp`), not just `propagate.hpp`. **This is the design update from our discussion:** thread `electrons.inverse_mass()` through *all* `ks_hamiltonian` construction sites, defaulting to all-ones (so Case A and every existing run stay bit-for-bit identical).
>
> **How convergence is ensured in Case B:** nothing special is needed — the SCF loop (density mixing, diagonalisation, tolerances) is mass-agnostic; changing the kinetic prefactor just shifts the eigenvalues/orbitals, and the self-consistent cycle converges to the muon-jellium GS the same way it does for electrons. Two things to *check* (validation gates, not blockers): (i) the `m=1` case reproduces the current electron GS **bit-for-bit** (proves the fork is inert when off); (ii) at fixed dimensionless r_s the muon GS energy equals the electron GS energy × `m` in effective units (proves the fork scales correctly). See engine notes §2 "Validation gates."
>
> ⭐ **One-line rule to remember:** *mass affects the GS iff the thing whose mass you changed is present in the GS.* Projectile-injected-at-RT → GS untouched. Bath-is-muon (or band-structure) → GS must carry the mass.

---

## 💬 Added by feedback pass (2026-07-06) — two more things worth writing into your understanding

**2. The energy ledger is the subtle third file-touch.** Besides *applying* kinetic (in `operator()`), the same kinetic operator is used to *measure* kinetic energy for `E_total` (via `kinetic_expectation_value` → `laplacian_expectation_value`, summed in `energy.hpp`). If you scale the apply-kernel by mass but forget the expectation-value kernel, the muon's kinetic energy in the ledger is silently wrong. Both use the SAME per-state factor → the "clean energy ledger" the muon route was chosen for.

**3. The rescaling caveat on "muon dynamics in jellium" (the open Q4 debate).** For an all-muon jellium at fixed r_s, *everything rescales from the electron problem* — including the quantum-vs-classical stopping difference (`S_muon = m²·S_electron` in real units, identical in effective units). New physics appears only when **projectile mass ≠ bath mass** (muon in electron jellium) OR at a strongly-correlated density regime (muon jellium at fixed *physical* density → r_s ~ hundreds). This is still open between us — see engine notes §3 "CRITICAL rescaling fact."