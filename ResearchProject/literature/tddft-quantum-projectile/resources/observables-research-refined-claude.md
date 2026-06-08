# Observables for Quantum Wavepacket Projectiles in Jellium under RT-TDDFT: A Beginner-Oriented Synthesis

## 0. Orientation: What the Problem Really Is

You are running RT-TDDFT on a jellium box, with one extra "thing" — either an additional Kohn–Sham (KS) electronic packet you have grafted onto the Fermi sea (the *indistinguishable* case) or a quantum proton/ion with its own coordinate (the *distinguishable* case). You want to track how energy and momentum flow out of that thing into the host. The fundamental difficulty is that, in TDDFT, the KS orbitals are **auxiliary** mathematical objects: their individual energies ε_i(t) are not, in general, expectation values of any physical observable (they are eigenvalues of an effective one-body Hamiltonian whose interpretation is only fully justified in the ground state, via Janak's theorem and the HOMO Koopmans' theorem). What is physical are densities, currents, density matrices, response functions, and quantities derived from them.

Your prior survey covered the standard "field-like" toolbox (excess γ, control volumes, S(q,ω), detector fluxes, Wigner/Husimi, local energy maps, NTOs). This report goes *beyond* by (i) building the underlying many-body machinery from scratch, (ii) importing concepts from nuclear TDHF, plasma physics, polaron physics, quantum impurity physics, and quantum thermodynamics that are mathematically equivalent to your problem, and (iii) giving you concrete handles for new diagnostics not covered before — Anderson orthogonality, Loschmidt-echo / work statistics, GW spectral-weight transfer, Balian–Vénéroni variance, t-SURFF for matter-wave packets, and quantum-information measures.

A crucial recent paper to ground your project: **Nazarov & Gross, "Stopping power of electron liquid for slow quantum projectiles," arXiv:2510.26222 (2025)** — this is, to our knowledge, the first systematic stopping-power formulation that treats the projectile fully quantum-mechanically using the Exact Factorization, exactly the framework you need for the *distinguishable* case. They derive that, *with the same charge and velocity*, projectiles of different mass have different SPs — a uniquely quantum effect arising because lighter projectiles delocalize more strongly. This is the theoretical home of your problem.

---

## 1. Foundational Theory (built up from scratch)

### 1.1 The one-body reduced density matrix γ(r,r′,t)

For an N-electron pure many-body state |Ψ(t)⟩, the **one-body reduced density matrix (1-RDM)** is

  γ(r,r′,t) = N ∫ Ψ*(r′,r₂,…,r_N,t) Ψ(r,r₂,…,r_N,t) dr₂…dr_N.

It encodes everything needed to evaluate any one-body observable Ô = Σᵢ ô(rᵢ): ⟨Ô⟩ = Tr[ô γ]. The diagonal γ(r,r,t) is the density n(r,t); off-diagonal elements encode coherence and momentum content (its Fourier transform γ(k,k′) gives the momentum distribution n(k) on the diagonal).

In Kohn–Sham TDDFT, the *non-interacting* 1-RDM is exactly γ_KS(r,r′,t) = Σᵢ fᵢ φᵢ*(r′,t) φᵢ(r,t). By Runge–Gross and the van Leeuwen theorem, γ_KS reproduces the *true* density n(r,t) but is otherwise different from the interacting γ. This is critical: **only quantities expressible through n(r,t) (and currents, via the gauge-extended theorems) are formally exact in TDDFT**; observables that depend on the off-diagonal of γ are KS approximations. Nevertheless, in practice the KS 1-RDM is used as a proxy for the interacting one when computing momentum distributions and spectral diagnostics.

The "excess" 1-RDM δγ(t) = γ(t) − γ_eq is exactly the right object because (a) it is gauge-invariant under U(1) phases of orbitals, (b) all conserved quantities (N, P, E in the Hartree+xc-decomposable sense) are linear functionals of γ, and (c) δγ admits a natural splitting into a "packet" part and a "host disturbance" part by a projector decomposition.

**Key reference for the 1-RDM as the central object:** Coleman & Yukalov, *Reduced Density Matrices*, Lecture Notes in Chemistry (2000); Mazziotti, *Reduced-Density-Matrix Mechanics*, Adv. Chem. Phys. 134 (2007).

### 1.2 The dynamic structure factor S(q,ω) and density-density correlations

S(q,ω) is the Fourier transform of the **density–density correlation function**

  S(q,ω) = (1/2πN) ∫ dt e^{iωt} ⟨δn̂(q,t) δn̂(−q,0)⟩,

where δn̂(q) = n̂(q) − ⟨n̂(q)⟩. Physically, S(q,ω) is the probability density (per ω) of exciting the many-body system into a state with momentum transfer q and energy transfer ω — it is exactly what an inelastic scattering experiment measures.

The **fluctuation–dissipation theorem (FDT)** ties S(q,ω) to the imaginary part of the (retarded) density response function χ(q,ω):

  S(q,ω) = −(ℏ/π) [1 + n_B(ω)] Im χ(q,ω)   (zero T:  S(q,ω) = −(ℏ/π) Im χ(q,ω) for ω>0).

This is the bridge between **dissipation** (Im χ — energy absorbed from a perturbing field) and **fluctuation** (S — equilibrium density fluctuations). For a stopping-power problem, **the Bethe sum rule** ∫ω S(q,ω) dω = N q²/2m is conserved exactly, and the f-sum rule gives a direct check of any S(q,ω) you compute from RT-TDDFT.

The classical Bethe–Bloch stopping power can be written as

  S_e(v) = (2Z²/πv²) ∫₀^∞ ω dω ∫_{ω/v}^∞ (dq/q) Im[−1/ε(q,ω)],

so the *whole* low/intermediate-velocity stopping problem is determined by the loss function Im[−1/ε(q,ω)] = (4π/q²)(−Im χ)/|ε|² of the host. This is the reason **S(q,ω) is the single most informative scalar diagnostic** in your simulation: if you can extract it cleanly from RT-TDDFT, you have direct access to the canonical stopping integrand.

### 1.3 Linear response (Kubo) and the Lindhard function

Kubo's formalism says that, to first order in a perturbation V̂_ext(t) = ∫dr v(r,t) n̂(r), the induced density change is

  δn(r,t) = ∫ dt′ ∫ dr′ χ(r,r′,t−t′) v(r′,t′),

with the retarded response function

  χ(r,r′,t) = −(i/ℏ) θ(t) ⟨[n̂(r,t), n̂(r′,0)]⟩.

For non-interacting jellium electrons at T=0, evaluating this commutator on the Fermi-sea ground state gives the **Lindhard function** χ⁰(q,ω):

  χ⁰(q,ω) = 2 ∫ d³k/(2π)³ [f(ε_k) − f(ε_{k+q})] / [ℏω + ε_k − ε_{k+q} + iη].

Inside the **electron–hole continuum** (|q²/2m − vF q| < ω < q²/2m + vF q), Im χ⁰ ≠ 0 and the system can absorb (this is single-particle Landau damping, see §2.4). Outside, only collective modes (plasmons) carry weight. The interacting RPA dielectric is ε^RPA(q,ω) = 1 − v_q χ⁰(q,ω), and the full density response is χ^RPA = χ⁰/[1 − v_q χ⁰].

The pedagogical paper **Mihaila, "Lindhard function of a d-dimensional Fermi gas," arXiv:1111.5337 (2011)** derives all this in closed form — start there. The standard textbook source is **Giuliani & Vignale, *Quantum Theory of the Electron Liquid* (Cambridge, 2005)**, especially Chapters 4, 5 and 8.

### 1.4 Exact factorization (EF) for the distinguishable case

Abedi, Maitra & Gross [PRL 105, 123002 (2010); J. Chem. Phys. 137, 22A530 (2012)] showed that any electron–nuclear wavefunction Ψ(r,R,t) of the full TDSE can be written *exactly* (not approximately, not in BO) as

  Ψ(r,R,t) = χ(R,t) Φ_R(r,t),  with the partial-normalization condition ∫|Φ_R(r,t)|² dr = 1 ∀ R,t.

The marginal probability |χ(R,t)|² is the **exact** nuclear density. Inserting this ansatz into the TDSE produces:
  - a nuclear TDSE i ∂_t χ = [(−iℏ∇+A(R,t))²/2M + ε(R,t)] χ, with a **time-dependent potential energy surface (TDPES)** ε(R,t) = ⟨Φ_R|Ĥ_BO + Û_coup − iℏ∂_t|Φ_R⟩_r and a **time-dependent vector potential** A(R,t) = ⟨Φ_R| − iℏ∇_R |Φ_R⟩_r;
  - a coupled equation for Φ_R(r,t) involving R-derivatives of χ.

For the projectile-in-jellium problem, **R is the projectile coordinate** and r is the electronic many-body coordinate. The electronic stopping force is then *literally* −∇_R ε(R,t), and the friction tensor is encoded in A(R,t). This is exactly the framework the recent **Nazarov–Gross (2025)** paper builds. Crucially:

- The **TDPES exhibits dynamical "steps"** (Abedi, Agostini, Suzuki, Gross, PRL 110, 263001 (2013)) when nuclear wavepackets bifurcate — the analog in your problem is whenever your projectile branches into transmitted/reflected/captured channels.
- A non-trivial geometric phase / vector potential in A(R,t) signals **non-adiabatic coupling that classical Ehrenfest cannot capture**.
- Mean-field/Ehrenfest stopping is recovered as a *specific* approximation (uncorrelated factorization) of EF, which is why standard RT-TDDFT-with-classical-projectile is a special case.

### 1.5 Non-equilibrium Green's functions and Kadanoff–Baym

The full two-time, contour-ordered Green's function G(1,2) = −i⟨T_C ψ̂(1) ψ̂†(2)⟩ on the Keldysh/Schwinger–Kadanoff contour contains *all* one-body information including spectral content. Its components G^>, G^<, G^R, G^A satisfy the **Kadanoff–Baym equations (KBE)**:

  [iℏ∂_{t₁} − h_HF(1)] G^≷(1,2) = ∫_C dt̄ Σ(1,t̄) G^≷(t̄,2),

with self-energy Σ accounting for correlations (e.g. GW, second Born, T-matrix). **For your problem the KBE perspective is valuable because**:
1. Energy *exchange* is not generally conservative at the TDDFT mean-field level when one uses simple xc functionals; KBE with conserving Φ-derivable Σ (Baym 1962) gives strict global energy/momentum/particle conservation by construction.
2. The spectral function A(k,ω) = −(1/π) Im G^R(k,ω) tells you how a sharp KS quasiparticle line broadens and develops satellites (plasmon polaron, plasmaron) when the projectile is dressed.
3. The **generalized Kadanoff–Baym ansatz (GKBA)** of Lipavský–Špička–Velický (1986) reduces two-time KBE to one-time for the 1-RDM, making it numerically comparable to RT-TDDFT while retaining genuine many-body memory; see Hermanns, Balzer & Bonitz, J. Phys. Conf. Ser. 427, 012008 (2013), arXiv:1205.4427.
4. **Stochastic interpretation**: Greiner & Leupold (hep-ph/9809296) showed the KBE for a closed system is exactly the ensemble-averaged Langevin equation with FDT-related noise — i.e. dissipation in your projectile's reduced description **must** come with a corresponding noise term, a fact you can use as a consistency check on any reduced-density-matrix dynamics you extract.

For practical use, see Dahlen & van Leeuwen, PRL 98, 153004 (2007); Stefanucci & van Leeuwen, *Nonequilibrium Many-Body Theory of Quantum Systems* (Cambridge 2013) — the canonical textbook.

### 1.6 Bethe–Bloch and its quantum descendants

Classically, Bohr (1913) and then Bethe (1930) derived

  −dE/dx = (4π Z² e⁴ n_e / m_e v²) ln(2 m_e v² / I)

for a fast charged projectile. Bloch (1933) interpolated between Born (Bethe) and classical (Bohr) regimes; Lindhard & Scharff (1953/61) extended it for slow ions in the electron gas, giving the **friction coefficient** Q(v) = v Σ_l (l+1) sin²(δ_l − δ_{l+1}) [phase shifts of the screened scattering on jellium]. Echenique, Nieminen & Ritchie [Solid State Commun. 37, 779 (1981)] showed how to compute Q with DFT scattering. This *non-linear* DFT-scattering result is the gold standard against which any RT-TDDFT slow-velocity SP must agree, and you should use it for benchmarking.

The full quantum unification of Bethe–Bloch and Lindhard–Scharff regimes, as well as the *new* mass-dependent corrections of the wavepacket projectile, sits inside the Nazarov–Gross EF framework (§1.4).

---

## 2. Cross-Disciplinary Analogues — Where Equivalent Problems Have Been Solved

### 2.1 Nuclear physics: TDHF, TDDFT for nuclei, and the Balian–Vénéroni variance

Nuclear TDHF/TDHFB (Hartree–Fock–Bogoliubov) has wrestled with **exactly your problem** for forty years: a quantum projectile (incident nucleus) collides with a target nucleus, you propagate Slater determinants in real time, and you ask "how much energy was dissipated, how much was transferred, what fragments emerged?" Despite very different physics, the *observables* are mathematically identical to yours. Key transferable ideas:

- **Particle-number projection** to define fragments unambiguously: Simenel, PRL 105, 192701 (2010) — given a Slater determinant Φ(t), the probability that exactly N nucleons are in spatial region V is P_V(N) = ⟨Φ|P̂_V^N|Φ⟩, computed via a determinant of the spatial overlap submatrix. Apply directly: define V as the wake region behind your projectile and compute the probability of N excess electrons there.
- **One-body dissipation (Loebl, Bertsch, Norenberg, Koonin)**: in TDHF, dissipation is mediated by single-particle scattering against a self-consistent moving mean field. Diagnostics include the **momentum-space quadrupole tensor** Q_k^{ij}(r,t) measuring local anisotropy of the momentum distribution — this is your local "thermalization" diagnostic and a natural companion to the Wigner-function quadrupole. See Umar et al., J. Phys. G 37, 064037 (2010); Loebl et al., Phys. Rev. C 84, 034608 (2011).
- **Stopping power and transparency in TDHF**: Reinhard, Maruhn et al., Phys. Rev. C 85, 014614 (2012) (arXiv:1208.5805) introduced precisely the local momentum-space quadrupole as a phase-space relaxation probe in heavy-ion stopping. This is *literally the analog of your jellium stopping problem with a quantum projectile*, only in the nuclear context. The recipe is: compute m-th moments of the local momentum distribution, build the deviatoric tensor, integrate over a control volume, and watch its decay.
- **Balian–Vénéroni (BV) variational principle**: Simenel, PRL 106, 112502 (2011); Williams, Simenel, Phys. Rev. C 88, 064601 (2013). TDHF underestimates **fluctuations** of one-body observables (variances, dispersions); BV gives a *systematically improved* estimator that requires propagating an adjoint state with a slightly perturbed initial condition. Translate to your case: TDDFT-based ⟨E_pkt⟩ may be reasonable, but Var(E_pkt) is generally underestimated; the BV correction would be a quantitatively new diagnostic capturing energy-loss fluctuations beyond mean field.
- **Stochastic mean-field (SMF), Ayik, Lacroix, Yilmaz**: PRC 85, 034616 (2012); PRC 102, 064619 (2020). Initialize an *ensemble* of Slater determinants with quantum-fluctuating initial 1-RDMs; their classical evolution under TDHF reproduces transport coefficients (friction, diffusion) consistent with FDT. A directly portable recipe for adding fluctuations on top of your deterministic RT-TDDFT.

The deepest takeaway: **the nuclear community has already solved the bookkeeping for "quantum projectile + many-fermion target" — borrow it wholesale.**

### 2.2 Plasma physics: Vlasov, Landau damping, phase mixing, and plasma echoes

The Vlasov-Poisson equation is the classical mean-field limit of the electronic TDDFT problem, and **Landau damping** (Landau 1946) is the textbook example of *collisionless, time-reversible* energy transfer from a wave to particles via phase mixing of resonant trajectories. The mathematics carries over to your problem essentially unchanged at the linear level. Useful imports:

- **Phase mixing as the mechanism of dissipation**: Free streaming f(x,v,t) = f_in(x−vt, v) means the Fourier transform of f decays in k-space, even though the entropy is exactly conserved. This is **identical** to the decay of off-diagonal γ elements after a wavepacket is launched — δγ "smears out" in k–k′ space. The Husimi-IPR you already had is a smoothed version of this.
- **Plasma echoes** (Gould, O'Neil, Malmberg 1967): two perturbations at times t₁ < t₂ with momenta k₁, k₂ produce a coherent density spike at t_echo = k₂ t₂/(k₂ − k₁). In RT-TDDFT, **a quantum projectile that scatters twice (e.g. through periodic boundaries) should produce an analogous echo signature** — a non-trivial diagnostic you have not had on your list.
- **Bedrossian, Mouhot & Villani** [arXiv:1712.08498; Acta Math. 207, 29 (2011)] proved nonlinear Landau damping rigorously; their estimates of Gevrey-class regularity decay are exactly the rate at which any "free-streaming" component of δγ should phase-mix away in your simulation. Use them as *a priori* tests on your numerics: if your δγ phase mixes faster than this, your boundaries are absorbing; slower, you are under-resolved.

### 2.3 Polaron physics: a paradigm for the dressed projectile

The Fröhlich polaron — an electron coupled to a phonon bath — is the closest atomic analog of your quantum projectile coupled to the jellium *plasmon* bath (you can think of "plasmaron" rather than polaron). Tracked observables in real time:

- **Polaron formation time** τ_form = time at which the dressed wavepacket separates from a "wake" of radiated bosons. Mishchenko, Prokof'ev et al., Phys. Rev. B 71, 035105 (2005); arXiv:cond-mat/0310226. They define it **as the time required for the polaron to physically separate from the radiated phonons** — exactly the time-scale you should report for your packet to separate from its plasmon wake.
- **Spectral weight transfer**: the bare-electron quasiparticle peak loses weight Z = 1 − ∂Σ/∂ω|_{ω=0}; the "missing" 1−Z spectral weight goes into incoherent (satellite) structure. Compute this via post-RT GW or via the cumulant expansion (Aryasetiawan, Hedin, Karlsson, PRL 77, 2268 (1996)). This is a powerful complement to your packet-survival N_pkt: where your N_pkt is a *projector overlap*, the GW Z is a *spectral* weight.
- **First-principles non-equilibrium polaron formation**: Garcia-Herrero, Emeis, Caruso et al., arXiv:2601.21810 (2026) develop a quantum-kinetic theory for **real-time polaron formation in MgO** under pump–probe excitation. The structure of their observables (electron–phonon population pumping, wavepacket localization metric, dynamical fingerprint) maps cleanly onto your problem with phonons → plasmons.

### 2.4 Quantum impurity physics: Anderson orthogonality, X-ray edge, Kondo

If you treat a *static* quantum impurity (or your slowed, captured projectile) as suddenly switched on, the Fermi sea responds by generating an infinity of low-energy particle-hole pairs. This produces three closely related universal phenomena:

- **Anderson Orthogonality Catastrophe (AOC)** [Anderson, PRL 18, 1049 (1967)]: the overlap |⟨FS_initial|FS_final⟩| ~ L^{−γ/2} → 0 in the thermodynamic limit, where γ = Σ_l (2l+1)(δ_l/π)² and δ_l are scattering phase shifts. **In RT-TDDFT this is computable directly as the Slater-determinant overlap |det S(t)|, where S_{ij} = ⟨φᵢ_eq | φⱼ(t)⟩**. A power-law decay of this overlap with system size at long times is a smoking-gun signature of irreversible (in the thermodynamic-limit sense) projectile capture.
- **X-ray edge / Mahan-Nozières-de Dominicis singularity**: gives a power-law in the absorption spectrum at the threshold, with exponent determined by the same δ_l. Here it tells you the *spectral shape* of the energy-loss distribution at low ω — a sharp non-Lorentzian feature you can measure by Fourier transforming the density–density response.
- **Generalized Loschmidt echo** L(t) = |⟨Ψ(0) | e^{i Ĥ_eq t/ℏ} e^{−i Ĥ(t) t/ℏ} | Ψ(0)⟩|² and its many-body avatar (Dóra, Pollmann, Zaránd, PRL 111, 046402 (2013); Münder, Weichselbaum, Goldstein, Gefen, von Delft, PRB 85, 235104 (2012); Zangara, Pastawski et al., PRA 86, 012322 (2012)). **For Luttinger liquids, the Loschmidt echo decays as a power law with universal exponent** twice the AOC exponent. In your jellium box, computing this overlap (you have all the ingredients in RT-TDDFT) gives an exquisitely sensitive measure of how strongly the projectile reorganizes the host.

This connects directly to **quantum work statistics** and the two-time-measurement (TTM) protocol: the characteristic function of work G(u) = ⟨Ψ_eq | e^{iuĤ(0)/ℏ} e^{−iuĤ_pert/ℏ} | Ψ_eq⟩ is *literally the Loschmidt amplitude*, so the work probability distribution P(W) — i.e. the probability that the host absorbs energy W from the projectile — is the Fourier transform of the Loschmidt echo. See **Goold, Plastina, Gambassi & Silva, "The role of quantum work statistics in many-body physics," in *Thermodynamics in the Quantum Regime* (Springer 2018), arXiv:1804.02805**. This is, in our view, **the single most underused diagnostic for your problem**: it gives a probabilistic energy-loss distribution rather than just a mean, and connects naturally to fluctuation theorems (Jarzynski, Crooks).

### 2.5 Quantum transport and Landauer–Büttiker

For an injected packet entering a "sample" region of jellium and exiting (transmission/reflection), the Landauer formula T(E) gives transmission coefficients per energy channel. Energy-resolved transmission is computed by binning the asymptotic momentum distribution of the packet, exactly analogous to t-SURFF (§3.2). Büttiker's **virtual probe** trick — a fictitious reservoir whose chemical potential and temperature are adjusted to give zero net particle/energy current — is a clever way to define an effective local temperature of the jellium near the projectile track from purely one-body quantities. See Büttiker, Phys. Rev. B 33, 3020 (1986); Engquist & Anderson, Phys. Rev. B 24, 1151 (1981).

### 2.6 Open quantum systems: Caldeira–Leggett mapping

If we coarse-grain over the jellium degrees of freedom, the projectile is a quantum particle in a bosonic bath of plasmons + electron–hole pairs. **Caldeira & Leggett, Physica A 121, 587 (1983)** mapped this to a master equation with a friction kernel γ(t) and a noise correlator ν(t) tied by the FDT. The bath spectral density J(ω) is *exactly* J(ω) = Σ_q (q-projection of dipole coupling)² Im[−1/ε(q,ω)] / ω — i.e. computable from your own RT-TDDFT loss function!

Useful observables that transfer:
- **Purity** P(t) = Tr[ρ_proj(t)²] — decoherence rate of projectile.
- **Decoherence functional** Γ(t) — see Hu, Paz, Zhang, PRD 45, 2843 (1992).
- **Non-Markovianity measures** (BLP, Breuer–Laine–Piilo, PRL 103, 210401 (2009); RHP, Rivas–Huelga–Plenio).
- **Relative entropy of system to instantaneous Gibbs state** as entropy production indicator.

Relating Caldeira–Leggett to Landau damping: Hagstrom & Morrison [arXiv:1008.5190 (2011)] proved CL is *exactly* a quantum analog of the Vlasov–Poisson Landau damping, closing the loop between §2.2 and §2.6.

### 2.7 Quark–gluon plasma jet quenching

Although the energy scales differ by ~25 orders of magnitude, the structure is uncannily similar: a hard parton (your "projectile") propagates through a hot bath (your "jellium"); the **transport coefficient q̂** = mean transverse momentum² gained per unit length plays the role of momentum diffusion. The most-used theoretical objects (BDMPS-Z, GLV, AMY, Higher-Twist) all start from a *medium-induced gluon emission spectrum* — the QCD analog of the Im[−1/ε(q,ω)] loss function. See Qin & Wang, Int. J. Mod. Phys. E 24, 1530014 (2015), arXiv:1511.00790; Cao & Wang, Rep. Prog. Phys. 84, 024301 (2021), arXiv:2002.04028. For *your* purposes the take-home is the **q̂ formalism**: define q̂(t) = d⟨P_⊥²⟩/dt for the projectile, separate it into elastic and radiative parts, and use this scalar as a robust complement to dE/dt.

### 2.8 Ultrafast/attosecond physics: t-SURFF and time-frequency analysis

In strong-field photoionization, the electron leaves the atom and acquires asymptotic momentum that is what is measured. **Time-dependent surface-flux (t-SURFF)** of Tao & Scrinzi, New J. Phys. 14, 013021 (2012); Scrinzi, NJP 14, 085008 (2012) is the clean way to extract momentum-resolved spectra from a finite-box TDDFT/TDSE simulation: place a closed surface S, monitor the outgoing flux of the absorbed wave-packet, and integrate against a Volkov (free) state to obtain the momentum-resolved amplitude b(k):

  b(k) = (i/(2π)^{3/2}) ∫₀^∞ dt e^{ik²t/2} ∮_S dS · [ψ(r,t) ∇χ_k*(r,t) − χ_k*(r,t) ∇ψ(r,t)] e^{−ik·r}.

For your problem, **t-SURFF is exactly the tool to build your "scattered projectile spectrum"** — i.e. the differential cross-section dσ/dE_kdΩ_k of the scattered packet. This decomposes the projectile's energy loss into *channel-resolved* contributions in the same way that ARPES decomposes a photoemission signal. Implementation: see also Mosert & Bauer, Comput. Phys. Commun. 207, 452 (2016) for the QPROP code.

Time–frequency methods you should also have on hand:
- **Gabor transform** for time-resolved S(q,ω) (windowed Fourier);
- **Continuous wavelet transform (Morlet)** for chirped projectile responses;
- **Empirical mode decomposition / Hilbert–Huang** for non-stationary intrinsic modes (these have been used in HHG analysis; see Chini, Wang, Cheng, et al., Sci. Rep. 3, 2941 (2013)).

---

## 3. Practical Implementation in Plane-Wave RT-TDDFT

### 3.1 Computing δγ, free-packet projector, and projector decomposition

In a plane-wave basis with N_b bands and N_k k-points, γ(r,r′,t) = Σ_{nk} f_{nk} φ_{nk}*(r′,t) φ_{nk}(r,t). In practice, store and update the band coefficients C_{Gnk}(t) = ⟨G | φ_{nk}(t)⟩. Then:

- **Excess 1-RDM**: δγ_{GG′}(t) = Σ_{nk} f_{nk} [C*_{G′nk}(t)C_{Gnk}(t) − C*_{G′nk}(0)C_{Gnk}(0)].
- **Free-packet projector**: build a reference Hamiltonian Ĥ_free that contains only the kinetic + (optionally) the static jellium ionic background, **without** the dynamical xc response. Time-evolve a "phantom" packet under Ĥ_free starting from the same initial wavepacket; call its instantaneous 1-RDM γ_free(t). The packet projector is then P_free(t) = γ_free(t) / Tr[γ_free(t)²] (or, more cleanly, the orthogonal projector onto the column span of γ_free's natural orbitals with non-negligible occupations).
- **Survival, energy, momentum**: N_pkt(t) = Tr[P_free δγ]; E_pkt(t) = Tr[P_free Ĥ_KS δγ]; P_pkt(t) = Tr[P_free p̂ δγ]. Note: Ĥ_KS is gauge/xc-dependent — for the **most defensible** energy bookkeeping use the *interacting* total-energy functional E[n] (Hartree + xc + KE + ext) and compute its derivatives, not the sum of orbital eigenvalues.

### 3.2 t-SURFF in RT-TDDFT

Place a smooth absorbing layer (CAP or mask function W(r) = sin²(π(r−r_0)/(2L))) outside a flux surface S inside the simulation box. The outgoing-wave amplitude is

  b(k,t) = (1/(2π)^{3/2}) ∫₀^t dt′ e^{i(k²/2)(t−t′)} ∮_S dS · J_k(r,t′),

where J_k = (1/2i)(ψ ∇e^{−ik·r} − e^{−ik·r}∇ψ) projected onto the surface normal. This works for both the indistinguishable-electron and the distinguishable-projectile case (just track the projectile's wavefunction χ(R,t) for the latter). See Scrinzi (2012) for the multi-electron generalization (channel resolution); for solids see Yamada & Yabana, PRB 99, 245103 (2019) on RT-TDDFT photocurrents.

### 3.3 Boundary conditions and box-size considerations

For a quantum projectile you need:
1. A box larger than 2 × (projectile coherence length + jellium screening length λ_TF ≈ 1/k_TF);
2. A complex absorbing potential or transparent boundary far from the surface S, so that backscattered packet does not pollute t-SURFF;
3. A time step δt ≤ 0.05 / E_max where E_max is the largest plane-wave kinetic energy retained;
4. For Fourier processing of S(q,ω): apply a window (Hann, Blackman) before FFT; the frequency resolution is Δω = 2π/T_total, the maximum unambiguous ω is π/δt;
5. Convergence in N_e (jellium electron count, equivalently r_s box dimensions) must be checked because **AOC, Loschmidt and δγ off-diagonals all scale with system size** — do at least three sizes.

### 3.4 Code support landscape (as of mid-2026)

| Code | Basis | Stopping power | Excess γ accessible? | t-SURFF / flux | Wannier MLWF |
|---|---|---|---|---|---|
| **Octopus** (Marques et al.) | real-space grid | yes (classical projectile, Andrade et al.) | partially via density-matrix module | yes (post-processing module) | via ports |
| **Qb@ll/Qbox** (Schleife, Draeger et al.) | plane wave PW | yes — large-scale (e.g. solvated DNA, INCITE-class runs) | yes (orbital/MLWF output) | needs custom post | yes |
| **SALMON** (Yabana, Noda et al.) | real-space grid | yes (light–matter focus, Sato 2023) | yes via TDKS orbital files | yes | external |
| **exciting** (Kazempour et al., 2021) | LAPW+lo | mainly LR/RT-TDDFT for excitations | yes (full all-electron) | external | LAPW Wannier |
| **eQE / Quantum ESPRESSO RT-TDDFT** | PW PP | yes (linear-scaling extensions) | yes | external | via wannier90 |
| **PETRA / KSSOLV / ChronusQ / NWChem** | mixed | partial | yes (RDM output) | rare | code-dependent |

For the **plane-wave route specifically**, **Qb@ll** is the most battle-tested for stopping-power applications (Schleife et al., J. Chem. Phys. 137, 22A546 (2012); Yost & Kanai et al., J. Chem. Phys. 155, 100901 (2021)). For all-electron precision use exciting LAPW+lo (Pela & Draxl, arXiv:2102.02630). For exact factorization with quantum nuclei you will likely have to hack — exciting and Octopus both expose nuclear wavefunctions in restricted form. The Nazarov–Gross 2025 work uses a custom code.

### 3.5 Post-processing toolkits

- **WannierTools / Wannier90** for MLWF construction along the trajectory (the Schleife-group recipe for projector-resolved energy transfer).
- **TRIQS/cthyb-style cumulant code** for post-RT spectral function reconstruction.
- **Husimi/Wigner**: the Octopus utility `oct-wigner`; standalone Python (QuTiP, qutip-qip) provides Wehrl entropy and Husimi IPR.
- **Scientific Python**: scipy.signal.stft for Gabor transforms; PyWavelets for CWT.
- **NESSi** (Schüler, Eckstein, Werner et al., Comput. Phys. Commun. 257, 107484 (2020)) for parallel Kadanoff–Baym; lets you compute true two-time G(t,t′) and post-process A(k,ω) along the dynamics — a powerful complement to RT-TDDFT for *correlation* diagnostics in jellium.

---

## 4. New Ideas and Unexplored Angles

### 4.1 Quantum-information measures

(a) **System–bath entanglement entropy** S_E = −Tr[ρ_proj log ρ_proj] in the distinguishable case, where ρ_proj = Tr_jellium |Ψ⟩⟨Ψ|. For Slater-determinant electronic states this reduces to a determinant computable from γ alone via the Peschel formula S_E = −Σ_α [n_α ln n_α + (1−n_α) ln(1−n_α)], where n_α are eigenvalues of the spatially restricted γ. **This gives an O(N³) computable, basis-independent measure of how strongly your projectile is dressed by the bath.** See Peschel, J. Phys. A 36, L205 (2003); Latorre & Riera, J. Phys. A 42, 504002 (2009) for free-fermion entanglement.

(b) **Mutual information** I(A:B) between two disjoint regions (e.g. a sphere around the projectile and a far-field region) carries information about how the disturbance propagates non-locally — particularly relevant for plasmon-mediated long-range transfer. See Calabrese & Cardy, J. Phys. A 42, 504005 (2009).

(c) **Quantum Fisher information** F_Q(λ) of any control parameter λ characterizes the *maximum* statistical resolution of energy/momentum readouts: F_Q ≥ Σ_n ω_n² P_n(W) − ⟨ω⟩² (related to work fluctuations). Used in quantum metrology (Tóth & Apellaniz, J. Phys. A 47, 424006 (2014)); a natural figure-of-merit for the precision of any extracted stopping force.

### 4.2 Loschmidt-echo and full work statistics

As argued in §2.4, the Loschmidt amplitude L(t) = ⟨Ψ_eq | Û_eq^†(t) Û(t) | Ψ_eq⟩ is the characteristic function of work, so Fourier-transforming it gives the energy-loss probability distribution P(W,t). For a Slater-determinant initial state and a quadratic (KS) Hamiltonian, this reduces to a determinant of N×N matrices and is directly computable along your RT-TDDFT trajectory. **Three new diagnostics emerge for free:**

- **Mean work** ⟨W⟩(t) = ∫W P(W,t) dW — equal to ΔE_total absorbed; you already have this.
- **Variance Var(W)(t)** — "spread of energy-loss outcomes," missing from prior tooling.
- **Higher cumulants and entropy production** Σ(t) = ⟨W⟩/T − ΔF (Crooks / Jarzynski).

### 4.3 GW spectral-weight transfer as a diagnostic

After the wavepacket equilibrates (or at any fixed time), perform a one-shot G⁰W⁰ on the *instantaneous* KS state. Compute Z = (1 − ∂ReΣ/∂ω)^{−1} of the projectile-localized natural orbital. The deviation 1 − Z(t) is the spectral weight transferred from the coherent quasiparticle peak into incoherent satellites (plasmaron, two-plasmon shake-up). This is the *spectral* analog of the *projector-based* packet-survival N_pkt(t). For jellium, see Caruso et al., PRB 86, 081102 (2012); Hedin's original cumulant treatment (J. Phys.: Condens. Matter 11, R489 (1999)). Vlcek et al., PRB 98, 075107 (2018) (arXiv:1708.03848) discuss QP-satellite split in nano-systems.

### 4.4 Quantum thermodynamics: heat versus work

Decomposition (Alipour, Benatti, Floreanini, Mehboudi, Rezakhani; PRR 4, 023034 (2022)): for a system under non-trivial time-dependent driving by an external degree of freedom,

  dU = dW + dQ,  with dW = Tr[ρ Ḣ] dt,  dQ = Tr[ρ̇ H] dt.

In your simulation, treating the projectile motion as the "driving," dW is the work done by the projectile's mean trajectory on the host (the "mean-field" energy transfer) and dQ is the genuine *heating* of the host (the part one would identify with thermal fluctuations). The split is gauge-dependent on the projectile/host partition but, once partitioned, gives a clean thermodynamic interpretation. Coupled with the FDT-fluctuation tracker, this provides a route to **define an effective electronic temperature near the projectile track** independent of any equilibrium fit.

### 4.5 Machine learning for channel identification and surrogate models

Two recent threads are directly applicable:
- **Surrogate stopping-power models**: Ward et al., arXiv:2311.00787 (2023) train an ML model on RT-TDDFT outputs to predict stopping force from {projectile position, velocity, ground-state density}. This factor-10⁷ acceleration immediately ports to quantum-packet projectiles when augmented with packet width and momentum spread inputs.
- **Channel-clustering with unsupervised learning**: PCA/UMAP on the time-series of natural-orbital occupations of δγ can autonomously identify "transmission/reflection/excitation" channels — a quantitative replacement for visually inspecting NTOs.
- **Symbolic regression** on extracted dE/dx vs (v, mass) data to discover the wave-packet stopping law. Given the small recent data set from Nazarov–Gross 2025, this is an *open* problem.

### 4.6 Kadanoff–Baym ansatz–based reduced descriptions

A highly promising direction not present in the prior list: run Kadanoff–Baym for the *reduced* projectile system in the GKBA + 2nd-Born (or T-matrix) approximation, using the jellium χ⁰ from your RT-TDDFT as the bath-spectral input. This gives a *fully memory-resolved*, *strictly-conserving* reduced description of the projectile that, by construction, satisfies the Bethe sum rule and the FDT. Cf. Hopjan, Karlsson, Ydman, Verdozzi, von Friesen, PRL 116, 236402 (2016); Karlsson, Hopjan & Verdozzi, PRB 97, 125151 (2018). Tools: NESSi (Schüler et al., 2020).

### 4.7 Balian–Vénéroni variance for stopping-power fluctuations

Adapting the Simenel BV recipe (PRL 106, 112502 (2011)) to your problem: starting from the perturbed Slater determinant at time T, propagate *backwards* under a slightly modified Hamiltonian (with a one-body operator Â = position or momentum of the projectile inserted as a δ-perturbation in the boundary condition) and compute the resulting overlap. The departure of this variance from the naive ⟨Â²⟩ − ⟨Â⟩² gives a fluctuation-dissipation-consistent estimate of energy-loss variance beyond mean-field. This is a *new* diagnostic for the TDDFT stopping community.

### 4.8 Topological/geometric: Berry curvature of the EF vector potential

The exact-factorization vector potential A(R,t) carries a curvature ∇×A(R,t) that is, in general, nonzero whenever the electron-projectile coupling is non-adiabatic. In a quantum-stopping context, this curvature is the **electronic Lorentz/anomalous force** on the projectile, missing from any classical-projectile RT-TDDFT. Diagnostic: when classical Ehrenfest agrees with EF, ∇×A ≈ 0; when they diverge, ∇×A ≠ 0 and your projectile experiences geometric forces. See Requist, Tandetzky, Gross, PRA 93, 042108 (2016); Eich & Agostini, J. Chem. Phys. 145, 054110 (2016).

---

## 5. Reading List — A Progressive Pathway

### Tier 1 — Beginner foundations (read in order, 1–4 weeks)

1. **Quantum mechanics review with scattering**: Sakurai & Napolitano, *Modern Quantum Mechanics*, 3rd ed. (2020), Ch. 6–7 on scattering; Taylor, *Scattering Theory* (1972, Dover reprint).
2. **DFT/TDDFT primer**: Engel & Dreizler, *Density Functional Theory* (Springer 2011); Ullrich, *Time-Dependent Density-Functional Theory: Concepts and Applications* (Oxford 2012) — *this is the right entry-point to RT-TDDFT*.
3. **Many-body for jellium specifically**: Giuliani & Vignale, *Quantum Theory of the Electron Liquid* (Cambridge 2005), Ch. 1, 4, 5, 8 — covers Lindhard, RPA, plasmons, FDT.
4. **Linear response (Kubo)**: Mahan, *Many-Particle Physics*, 3rd ed. (2000), Ch. 3 and 5; or Bruus & Flensberg, *Many-Body Quantum Theory in Condensed Matter Physics* (Oxford 2004), Ch. 6 (response and Kubo).

### Tier 2 — Pedagogical entry to your specific tools (1–2 weeks each)

5. **Lindhard, pedagogical**: Mihaila, "Lindhard function of a d-dimensional Fermi gas," arXiv:1111.5337 (2011).
6. **Density matrix theory**: Coleman & Yukalov, *Reduced Density Matrices*, Lecture Notes in Chemistry 72 (Springer 2000).
7. **RT-TDDFT review and overview**: Provorse & Isborn, Int. J. Quantum Chem. 116, 739 (2016); the recent review Marques, Nogueira, Maitra et al., arXiv:2509.10745 (2025) "A Snapshot of Time-Dependent Density-Functional Theory."
8. **Stopping power, ab-initio review**: Schleife & Correa et al., "Quantum dynamics simulations of the electronic stopping power," J. Chem. Phys. 155, 100901 (2021); Race, Mason, Foulkes, Horsfield, Rep. Prog. Phys. 73, 116501 (2010).
9. **Exact factorization**: Abedi, Maitra, Gross, PRL 105, 123002 (2010) [foundational]; Agostini & Gross, Adv. Phys. X 1, 463 (2016) [review]; **Nazarov & Gross, arXiv:2510.26222 (2025)** [stopping-power-specific].

### Tier 3 — Bridging to advanced topics (2–4 weeks each)

10. **Nuclear TDHF**: Simenel, Eur. Phys. J. A 48, 152 (2012) — "Nuclear quantum many-body dynamics" [the canonical review]; Negele, Rev. Mod. Phys. 54, 913 (1982); Reinhard, Maruhn, Suraud et al., PRC 85, 014614 (2012) (arXiv:1208.5805) for stopping power in TDHF.
11. **Non-equilibrium Green's functions**: Stefanucci & van Leeuwen, *Nonequilibrium Many-Body Theory of Quantum Systems* (Cambridge 2013) [textbook]; Aoki et al., RMP 86, 779 (2014) [non-equilibrium DMFT review].
12. **Open quantum systems**: Breuer & Petruccione, *The Theory of Open Quantum Systems* (Oxford 2002); Caldeira, *An Introduction to Macroscopic Quantum Phenomena and Quantum Dissipation* (Cambridge 2014).
13. **Anderson orthogonality**: Anderson PRL 18, 1049 (1967) [original]; Mahan ch. 8 [textbook]; Münder et al., PRB 85, 235104 (2012) [modern NRG treatment]; Knap et al., PRX 2, 041020 (2012) [ultracold-atom realization].
14. **Quantum work statistics**: Esposito, Harbola, Mukamel, RMP 81, 1665 (2009); Goold, Plastina, Gambassi, Silva, in *Thermodynamics in the Quantum Regime* (Springer 2018), arXiv:1804.02805.
15. **Polaron real-time**: Mishchenko et al., PRB 71, 035105 (2005); Garcia-Herrero, Caruso et al., arXiv:2601.21810 (2026).
16. **Plasma Landau damping (modern)**: Mouhot & Villani, Acta Math. 207, 29 (2011); Bedrossian summary, arXiv:1712.08498.

### Tier 4 — Specialized seminal papers (selectable)

17. **t-SURFF**: Tao & Scrinzi, NJP 14, 013021 (2012); Scrinzi, NJP 14, 085008 (2012); Mosert & Bauer, Comput. Phys. Commun. 207, 452 (2016).
18. **Stochastic mean-field, Balian–Vénéroni**: Ayik, PLB 658, 174 (2008); Simenel, PRL 106, 112502 (2011); Lacroix, Tanimura, Ayik, Yilmaz, EPJ A 52, 94 (2016).
19. **Kadanoff–Baym numerics**: Dahlen & van Leeuwen, PRL 98, 153004 (2007); Schüler et al., NESSi paper, CPC 257, 107484 (2020).
20. **Echenique–Ritchie nonlinear DFT stopping**: Echenique, Nieminen, Ritchie, Solid State Commun. 37, 779 (1981); Echenique et al., Phys. Rev. A 33, 897 (1986).
21. **Jet quenching for cross-disciplinary perspective**: Qin & Wang, IJMPE 24, 1530014 (2015); Cao & Wang, Rep. Prog. Phys. 84, 024301 (2021).
22. **Quantum-thermodynamic decomposition**: Alipour et al., PRR 4, 023034 (2022); Esposito, Lindenberg, Van den Broeck, NJP 12, 013013 (2010).
23. **Wigner-function jellium**: Bonitz, *Quantum Kinetic Theory* (Springer 2016, 2nd ed.).
24. **Free-fermion entanglement**: Peschel, J. Phys. A 36, L205 (2003); Cheong & Henley, PRB 69, 075111 (2004).

### Tier 5 — Practical implementation papers

25. **Qb@ll RT-TDDFT for stopping**: Yost, Yao, Kanai, J. Chem. Phys. 155, 100901 (2021).
26. **Octopus**: Tancogne-Dejean, Oliveira, Andrade, Marques et al., J. Chem. Phys. 152, 124119 (2020).
27. **SALMON**: Noda, Sato, Yabana et al., Comput. Phys. Commun. 235, 356 (2019).
28. **exciting RT-TDDFT (LAPW+lo)**: Pela & Draxl, arXiv:2102.02630 (2021).
29. **ML surrogate**: Ward et al., arXiv:2311.00787 (2023).

A **suggested sequenced 12-week plan** for a beginner: weeks 1–3 Tier 1 (Ullrich + Giuliani–Vignale chapters); weeks 4–5 Tier 2 (Mihaila, Provorse–Isborn, Schleife review); weeks 6–8 Tier 3 nuclear TDHF and exact factorization (Simenel review + Abedi–Maitra–Gross + Nazarov–Gross); weeks 9–10 Tier 3 NEGF + Anderson orthogonality (Stefanucci–van Leeuwen + Münder–von Delft); weeks 11–12 Tier 4/5 picking by your project needs (t-SURFF, work statistics, code-specific papers).

---

## 6. Summary of Recommended New Diagnostics (beyond your prior list)

1. **Particle-number projection P_V(N)** in a control volume (Simenel-style).
2. **Local momentum-space quadrupole tensor Q_k^{ij}(r,t)** (Reinhard et al., nuclear TDHF — *direct phase-space relaxation probe*).
3. **Anderson-orthogonality determinant** |det⟨φ_eq|φ(t)⟩| (Slater-determinant overlap; cost N³).
4. **Many-body Loschmidt echo and work distribution P(W,t)** via Fourier transform of the determinant overlap.
5. **Variance of energy loss** Var(W) and **higher cumulants** — entropy production via Crooks/Jarzynski.
6. **Reduced-system entanglement entropy** via Peschel formula on the spatially restricted γ.
7. **Mutual information** between projectile region and far-field region.
8. **GW spectral weight Z(t)** of projectile-localized natural orbital and its 1−Z transfer to satellites.
9. **Transport coefficient q̂(t)** = d⟨P_⊥²⟩/dt (jet-quenching analog).
10. **Quantum work/heat decomposition** dW = Tr[ρ Ḣ]dt, dQ = Tr[ρ̇ H]dt.
11. **t-SURFF for matter-wave packet**: scattered-projectile energy-momentum spectrum.
12. **EF vector-potential curvature** ∇×A(R,t) — purely-quantum geometric force.
13. **Plasma-echo signature** at t = k₂t₂/(k₂−k₁) for double-scattering events.
14. **Polaron-style separation time** τ_form between dressed packet and radiated-plasmon wake.
15. **Balian–Vénéroni-corrected variances** of one-body observables (beyond TDHF/TDDFT mean-field underestimate).

---

## 7. Summary and Closing Remarks

Your problem of a quantum-mechanical wavepacket projectile in a jellium under RT-TDDFT is, structurally, *the same problem* that has driven separate developments in nuclear TDHF (since the 1970s), polaron physics (since Fröhlich), plasma Landau damping (since 1946), Anderson orthogonality (since 1967), quantum impurity physics, open-quantum-systems theory (Caldeira–Leggett), and quark–gluon plasma jet quenching. Each community has contributed bookkeeping tools and observables that are directly portable. The most consequential recent development for *your specific* setting is **Nazarov & Gross, arXiv:2510.26222 (2025)**, which embeds the Exact Factorization into stopping-power theory and rigorously demonstrates new mass-dependent quantum signatures absent in classical-projectile Ehrenfest–TDDFT.

Beyond the field-, subspace-, and scattering-based observables you already have catalogued, we recommend systematically incorporating: (i) **many-body overlap-based diagnostics** (Loschmidt echo, AOC, work statistics) — these are inexpensive determinants of single-particle overlap matrices in any KS code and they connect directly to thermodynamic quantities; (ii) **nuclear-TDHF-borrowed phase-space and particle-projection observables** with Balian–Vénéroni variance corrections; (iii) **t-SURFF on the projectile coordinate** to obtain a true scattered-projectile differential cross-section; (iv) **GW post-processing for spectral-weight transfer**; (v) **Peschel-formula entanglement entropies** as basis-free, gauge-invariant complements to the projector-based metrics. Together with your existing toolkit, these provide a comprehensive observable suite that is provably consistent with sum rules, fluctuation–dissipation, and thermodynamic-energy conservation.

Two closing cautions on epistemics: (a) most RT-TDDFT energy diagnostics are **gauge- and xc-functional-dependent**; whenever feasible, double-check with a Φ-derivable (conserving) GW or 2B Kadanoff–Baym calculation as an independent check on the conservation laws; (b) **the KS orbital eigenvalues are not observables** in the dynamical regime — interpret them with care, and prefer functionals of n(r,t), j(r,t), and γ(r,r′,t).