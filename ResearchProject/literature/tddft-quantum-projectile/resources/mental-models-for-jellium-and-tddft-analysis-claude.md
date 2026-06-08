# A Literature Survey for Building Mental Models of Electronic Energy Dissipation in Matter

## Orientation for the reader

This survey is structured around twelve interlocking topics that together build the conceptual scaffolding needed to interpret real‑time TDDFT (rt‑TDDFT) simulations of quantum electron wave‑packets dissipating energy in jellium and related targets. For each topic I give (i) the central physical picture or "mental model" a beginner should commit to memory, (ii) a short, curated list of references with full bibliographic detail and DOIs where available, and (iii) explicit pointers to how the material connects to the central problem of electron wave‑packets in a homogeneous electron gas treated with the INQ code. Where useful, I have flagged limitations of the cited works (e.g., adiabatic exchange‑correlation kernels, finite‑size and trajectory artefacts) so that they are not absorbed uncritically.

A note on tools: SRIM, PSTAR, ASTAR, the IAEA stopping‑power database, and the NIST IMFP database are referenced under topic 9. The student already has the Tsubonoya, Hu and Watanabe paper (Phys. Rev. B 90, 035416 (2014), DOI 10.1103/PhysRevB.90.035416) which establishes the wave‑packet methodology; the present survey is built on that foundation and explains what is needed to extend it from elastic LEED to inelastic energy dissipation.

---

## 1. Classical and semiclassical stopping power theory

**Mental model.** A charged projectile loses energy because it polarises the medium asymmetrically: behind it, an induced cloud lags and pulls the projectile back. In the high‑velocity limit the projectile passes faster than electrons can rearrange, and energy is transferred binary‑collision by binary‑collision to bound electrons (Bethe). At low velocities the projectile is dressed by a static screening cloud and only exchanges quanta with electrons near the Fermi level, giving a friction-like S∝v law (Fermi–Teller, Lindhard–Winther). Between these limits sits the Bragg peak, where the projectile velocity matches the typical orbital velocity of valence electrons and the dielectric loss function Im[−1/ε(q,ω)] is sampled most efficiently. Nuclear stopping (elastic recoil of host nuclei) dominates only at very low velocities, while electronic stopping dominates everywhere else.

**Foundational papers and reviews.**
- H. Bethe, "Zur Theorie des Durchgangs schneller Korpuskularstrahlen durch Materie", Ann. Phys. 397, 325 (1930) — the original Bethe formula derivation, and F. Bloch, Ann. Phys. 408, 285 (1933) for the Bohr–Bloch correction.
- N. Bohr, "The penetration of atomic particles through matter", Mat.‑Fys. Medd. Dan. Vidensk. Selsk. 18, no. 8 (1948) — classical impact‑parameter treatment.
- J. Lindhard, "On the properties of a gas of charged particles", Mat.‑Fys. Medd. Dan. Vidensk. Selsk. 28, no. 8 (1954) — the dielectric formulation; gives the now‑canonical expression for stopping in terms of Im[−1/ε(q,ω)].
- J. Lindhard and A. Winther, "Stopping power of electron gas and equipartition rule", Mat.‑Fys. Medd. Dan. Vidensk. Selsk. 34, no. 4 (1964) — derives both the high‑velocity Bethe limit and the friction limit from the same dielectric expression and gives the equipartition rule that single‑particle and collective channels each contribute half.
- E. Fermi and E. Teller, "The capture of negative mesotrons in matter", Phys. Rev. 72, 399 (1947), DOI 10.1103/PhysRev.72.399 — the slow‑projectile result S∝v in an electron gas (the prototype electronic friction).
- N. D. Mermin, "Lindhard dielectric function in the relaxation‑time approximation", Phys. Rev. B 1, 2362 (1970), DOI 10.1103/PhysRevB.1.2362 — adds particle‑number conserving collisional damping to Lindhard's response.

**Authoritative textbook treatments (the clearest pedagogical entry points).**
- P. Sigmund, *Particle Penetration and Radiation Effects: General Aspects and Stopping of Swift Point Charges*, Springer Series in Solid‑State Sciences vol. 151 (Springer, 2006), ISBN 978‑3‑540‑31713‑2 — the modern standard reference, with full derivations of Bohr, Bethe, Bloch and Lindhard from a unified notation; chapters 6–9 take the reader to the research frontier.
- P. Sigmund, *Particle Penetration and Radiation Effects, Volume 2: Penetration of Atomic and Molecular Ions*, Springer Series in Solid‑State Sciences vol. 179 (Springer, 2014), DOI 10.1007/978‑3‑319‑05564‑0 — heavy ions, charge‑state effects, slow‑ion regime (chapter 8 collects most of the Echenique–Nieminen–Ritchie nonlinear‑DFT material).
- J. F. Ziegler, J. P. Biersack and U. Littmark, *The Stopping and Range of Ions in Solids*, Pergamon (1985); the SRIM/TRIM code and the underlying parametrisations are described in J. F. Ziegler, M. D. Ziegler and J. P. Biersack, "SRIM – The stopping and range of ions in matter (2010)", Nucl. Instrum. Methods B 268, 1818 (2010), DOI 10.1016/j.nimb.2010.02.091 — useful as a benchmark and as an example of the limitations of binary‑collision/local‑density Monte Carlo: it ignores band structure, channelling, threshold behaviour in insulators, and finite‑lattice effects, all of which rt‑TDDFT can capture.
- J. D. Jackson, *Classical Electrodynamics*, 3rd ed. (Wiley, 1998), chapters 13 ("Collisions, energy loss, and scattering of charged particles") and parts of 7 — the classical Bohr derivation and the relativistic Bethe formula are presented with maximum clarity.

**Key energy regimes (a quick reckoning chart for the student).**
- v ≪ vF (Fermi velocity): friction regime, S ≈ Q(rs)·v with Q the friction coefficient (Echenique–Nieminen–Ritchie nonlinear DFT, see topic 4).
- v ≈ v_Bohr ≈ Z·e²/ħ: Bragg peak, where the loss function is integrated over the entire (q,ω) plane most efficiently.
- v ≫ vF: Bethe regime, S ≈ (4π Z² e⁴ n / mv²) ln(2mv²/I), with I the mean excitation energy.

**Connection to the project.** rt‑TDDFT in INQ at constant projectile velocity directly samples the dielectric loss function — its predictions can be benchmarked against the analytic Lindhard/Mermin S(v) curves at fixed rs, and any deviation diagnoses either non‑linear screening, non‑adiabatic xc, or finite‑size effects.

---

## 2. Jellium / homogeneous electron gas as a model system

**Mental model.** Jellium replaces the discrete ions by a uniform positive background. The resulting interacting electron gas is the simplest non‑trivial many‑body system in which to test response functions and energy‑loss formulae: ground‑state properties depend on a single parameter rs (the Wigner–Seitz radius in Bohr units). Its non‑interacting density–density response χ₀(q,ω) is the Lindhard function; the corresponding RPA dielectric ε(q,ω) = 1 − v(q)χ₀(q,ω) has two famous spectral features — the particle–hole continuum (Landau damping region) and a coherent plasmon dispersion ω(q) ≈ ωp + (3vF²/10ωp)q² that lives outside it until the dispersion enters the continuum and the plasmon decays.

**Pedagogical textbooks (in order of approachability).**
- G. F. Giuliani and G. Vignale, *Quantum Theory of the Electron Liquid*, Cambridge University Press (2005), ISBN 978‑0‑521‑82112‑4 — the modern standard. Chapters 4 and 5 give a meticulous derivation of χ₀(q,ω), the loss function, plasmon dispersion, Friedel oscillations, and Thomas–Fermi screening, with all dimensional factors restored. This is the single best pedagogical entry point for the wave‑packet‑in‑jellium project.
- D. Pines and P. Nozières, *The Theory of Quantum Liquids*, Vol. I: *Normal Fermi Liquids*, W. A. Benjamin (1966), reprinted Westview Press (1999), ISBN 978‑0‑201‑40774‑7 — the classic introduction to the dielectric function, sum rules, and the f‑sum rule.
- G. D. Mahan, *Many‑Particle Physics*, 3rd ed., Plenum/Kluwer (2000), ISBN 978‑0‑306‑46338‑9 — chapters 5–6 cover RPA, beyond‑RPA local‑field corrections, and the connection to inelastic scattering cross sections.
- W. M. C. Foulkes, "Tutorial on the homogeneous electron gas" — appears as Chapter 3 of the Springer volume *Theoretical and Computational Methods in Mineral Physics* and in summer‑school proceedings; an exceptionally clear short pedagogical introduction.
- N. W. Ashcroft and N. D. Mermin, *Solid State Physics*, Saunders (1976) — the friendliest textbook treatment of static Thomas–Fermi screening, Friedel oscillations, and Lindhard's static response.

**Foundational papers.**
- J. Lindhard, Mat.‑Fys. Medd. Dan. Vidensk. Selsk. 28, no. 8 (1954) — derivation of χ₀(q,ω).
- D. Bohm and D. Pines, "A collective description of electron interactions: III. Coulomb interactions in a degenerate electron gas", Phys. Rev. 92, 609 (1953), DOI 10.1103/PhysRev.92.609 — the RPA viewpoint on plasmons.
- N. D. Mermin, Phys. Rev. B 1, 2362 (1970) — number‑conserving relaxation‑time dielectric function.
- M. Corradini, R. Del Sole, G. Onida, M. Palummo, "Analytical expressions for the local‑field factor G(q) and the exchange‑correlation kernel Kxc(q,ω) of the homogeneous electron gas", Phys. Rev. B 57, 14569 (1998), DOI 10.1103/PhysRevB.57.14569 — convenient parametrisation of beyond‑RPA static local‑field corrections.

**Key parameters at metallic densities (worth memorising).** For rs = 2 (close to Al), one has kF ≈ 0.96 a.u., vF ≈ 0.96 a.u., EF ≈ 0.46 Ha, ωp ≈ 0.58 Ha, and a Thomas–Fermi screening length λTF ≈ 0.79 a.u. These set the natural scales for any wave‑packet calculation.

**Connection to the project.** The jellium target in INQ is the cleanest possible test bed. Many predictions are analytic: Im[−1/ε(q,ω)] is a known Lindhard kernel, the plasmon dispersion is computable to all orders in q, and the TDLDA xc kernel admits a closed‑form ALDA expression. Any departure of rt‑TDDFT results from these analytics flags either an aspect of the simulation that is unconverged or a many‑body feature beyond ALDA (e.g., short‑range correlations giving the Corradini local‑field factor).

---

## 3. Wake effects and dynamic screening

**Mental model.** A swift charge in an electron gas creates a comoving polarisation cloud whose Fourier transform involves the dielectric function at the Doppler‑shifted frequency ω = q·v. For v > vF the screening becomes dynamical and oscillatory: along the trajectory behind the projectile, the induced potential exhibits damped oscillations of wavelength ≈ 2π v / ωp ("plasmon wake"). In transverse directions it forms cone‑like structures, the Cherenkov‑type wake. These wake oscillations are the origin of the asymmetric induced charge density that exerts the retarding force responsible for stopping; they can also bind trailing electrons into "wake‑riding" states.

**Foundational papers.**
- J. Neufeld and R. H. Ritchie, "Passage of charged particles through plasma", Phys. Rev. 98, 1632 (1955), DOI 10.1103/PhysRev.98.1632 — the first explicit calculation of the wake potential.
- R. H. Ritchie, "Interaction of charged particles with a degenerate Fermi–Dirac electron gas", Phys. Rev. 114, 644 (1959), DOI 10.1103/PhysRev.114.644 — extension to a quantum electron gas.
- P. M. Echenique, R. H. Ritchie and W. Brandt, "History‑ and velocity‑dependent screening of ions in solids", Phys. Rev. B 20, 2567 (1979), DOI 10.1103/PhysRevB.20.2567 — the now‑standard ERB wake formalism, and the paper that introduced the language of wake‑bound states.
- J. C. Ashley and P. M. Echenique, "Influence of damping in an electron gas on vicinage effects in ion‑cluster energy loss", Phys. Rev. B 35, 8701 (1987), DOI 10.1103/PhysRevB.35.8701 — incorporation of plasmon damping into the wake potential.
- P. M. Echenique, F. J. García de Abajo, V. H. Ponce and M. E. Uranga, "Dynamic screening of ions in solids", Nucl. Instrum. Methods B 96, 583 (1995), DOI 10.1016/0168‑583X(95)00235‑9 — broad review of the wake‑potential program with applications to charge states and resonant‑coherent excitation.

**Connection to the project.** With a projectile electron rather than an ion the wake structure is similar but the sign of the induced charge differs and indistinguishability matters; nevertheless the dynamic structure factor S(q,ω) — directly accessible in rt‑TDDFT through the time‑dependent density — encodes exactly the same physics. Visualising the wake in real time is one of the most physically transparent diagnostics available to the project; the predicted wake wavelength 2πv/ωp gives an immediate check of the simulation.

---

## 4. Electronic friction and low‑velocity stopping

**Mental model.** At v ≪ vF the projectile sees a quasi‑static screened potential and inelastically scatters from the Fermi sea via low‑energy electron–hole pair creation. The friction coefficient Q (so that S = Qv) is determined by the transport (momentum‑transfer) cross section evaluated at kF in the self‑consistently screened potential. Equivalently, by the Hellmann–Feynman / Drude‑weight argument, Q is proportional to the Fermi‑level density of states of the perturbed system. Crucially, the linear‑response result (Lindhard) underestimates Q substantially for protons in low‑rs electron gases; the nonlinear DFT result (Echenique–Nieminen–Ritchie, Puska–Nieminen) is essential for quantitative agreement.

**Foundational papers.**
- P. M. Echenique, R. M. Nieminen and R. H. Ritchie, "Density functional calculation of stopping power of an electron gas for slow ions", Solid State Commun. 37, 779 (1981), DOI 10.1016/0038‑1098(81)91173‑X — the first nonlinear DFT calculation of the friction coefficient.
- P. M. Echenique, R. M. Nieminen, J. C. Ashley and R. H. Ritchie, "Nonlinear stopping power of an electron gas for slow ions", Phys. Rev. A 33, 897 (1986), DOI 10.1103/PhysRevA.33.897 — the canonical ENRA paper, now the standard benchmark for slow‑ion friction in jellium.
- M. J. Puska and R. M. Nieminen, "Atoms embedded in an electron gas: Phase shifts and cross sections", Phys. Rev. B 27, 6121 (1983), DOI 10.1103/PhysRevB.27.6121 — the embedding‑potential viewpoint, basis for the local‑density friction approximation (LDFA).
- I. Nagy, A. Arnau and P. M. Echenique, "Low‑velocity antiproton stopping power", Phys. Rev. B 40, 11983(R) (1989), DOI 10.1103/PhysRevB.40.11983 — extension to negative projectiles, illustrating Z₁³ Barkas‑like effects.
- A. Sarasola, R. H. Ritchie, E. Zaremba and P. M. Echenique, "Density functional theory‑based stopping power for 3D and 2D systems", Adv. Quantum Chem. 46, 1 (2004), DOI 10.1016/S0065‑3276(04)46001‑8 — modern review including 2D.

**Orbital‑dependent / tensorial electronic friction (relevant to the dynamics community).**
- M. Askerka, R. J. Maurer, V. S. Batista and J. C. Tully, "Role of tensorial electronic friction in energy transfer at metal surfaces", Phys. Rev. Lett. 116, 217601 (2016), DOI 10.1103/PhysRevLett.116.217601.
- R. J. Maurer, M. Askerka, V. S. Batista and J. C. Tully, "Ab initio tensorial electronic friction for molecules on metal surfaces: Nonadiabatic vibrational relaxation", Phys. Rev. B 94, 115432 (2016), DOI 10.1103/PhysRevB.94.115432.
- S. P. Rittmeyer, J. Meyer, J. I. Juaristi and K. Reuter, "Electronic friction‑based vibrational lifetimes of molecular adsorbates: Beyond the independent‑atom approximation", Phys. Rev. Lett. 115, 046102 (2015), DOI 10.1103/PhysRevLett.115.046102.
- C. P. Race, D. R. Mason, M. W. Finnis, W. M. C. Foulkes, A. P. Horsfield and A. P. Sutton, "The treatment of electronic excitations in atomistic models of radiation damage in metals", Rep. Prog. Phys. 73, 116501 (2010), DOI 10.1088/0034‑4885/73/11/116501 — an excellent, accessible review of how electronic friction is used in MD of irradiated metals, with critical assessment of the approximations.

**Connection to the project.** The friction coefficient is the lowest non‑trivial quantity computable in rt‑TDDFT: at small v one extracts Q = lim(v→0) S(v)/v from a series of constant‑velocity runs. INQ on jellium can compare directly to the ENRA tabulated Q(rs) values; any beyond‑ALDA or finite‑size deviation can then be ascribed cleanly. The student should also note that the Maurer/Tully orbital‑friction tensor and the local‑density friction approximation are conceptually different ways of compressing the same physics that rt‑TDDFT computes from scratch.

---

## 5. The quantum kick and current‑based mental models

**Mental model.** In an extended periodic system, switching on a sudden uniform vector potential A is equivalent to giving every electron a momentum kick. In the linear‑response limit, the Fourier transform of the resulting current J(t) directly yields the dynamical conductivity σ(ω) (and, by Kubo, all of optical response). The kick deposits energy that lives partly in the density distortion and partly in the macroscopic current (Drude weight); for a metal the latter persists indefinitely without dissipation in the closed system, while in jellium the f‑sum rule fixes the total integrated weight of σ(ω) to ne²/m. When a projectile passes through, it injects a localised current; rapidly established countercurrents in the host carry away the projectile's momentum loss in the form of plasmon‑like oscillations.

**Key references (these are simultaneously the references on extracting EELS and conductivity from rt‑TDDFT).**
- K. Yabana, T. Nakatsukasa, J.‑I. Iwata and G. F. Bertsch, "Real‑time, real‑space implementation of the linear response time‑dependent density‑functional theory", phys. stat. sol. (b) 243, 1121 (2006), DOI 10.1002/pssb.200642005 — the foundational paper on the impulsive‑kick technique for extracting linear response from rt‑TDDFT.
- G. F. Bertsch, J.‑I. Iwata, A. Rubio and K. Yabana, "Real‑space, real‑time method for the dielectric function", Phys. Rev. B 62, 7998 (2000), DOI 10.1103/PhysRevB.62.7998 — explicit derivation showing that a δ(t) kick gives the full ε(ω) in one run.
- W. Kohn, "Theory of the insulating state", Phys. Rev. 133, A171 (1964), DOI 10.1103/PhysRev.133.A171 — the conceptual basis for the Drude weight as a diagnostic of metallicity.
- G. F. Giuliani and G. Vignale (textbook above) — chapter on the f‑sum rule and the Kubo formalism.
- D. Vanderbilt, *Berry Phases in Electronic Structure Theory*, Cambridge (2018), ISBN 978‑1‑107‑15765‑1 — chapter 4 gives a particularly clean modern treatment of the Drude weight and current response.

**Connection to the project.** The "quantum kick" that the projectile delivers to the jellium is, to leading order, precisely an impulsive current source localised along its trajectory. Plotting the induced current density in the simulation and decomposing the energy into "kinetic in current" versus "potential in density" parts is the most physically transparent way to visualise where the dissipated energy lives at each instant; comparison to the f‑sum rule (which must be obeyed exactly in the closed simulation) provides a stringent numerical check.

---

## 6. Electron projectiles versus ion projectiles

**Mental model.** An electron projectile differs from an ion in three ways. (i) Sign and magnitude of charge change the wake (Barkas/Z₁³ effects). (ii) The projectile is identical to the target electrons, so antisymmetrisation between projectile and target is in principle required (in practice ignored at high energy where exchange corrections are small, but important near and below the Fermi energy). (iii) An electron of energy E above the Fermi level has an inelastic mean free path (IMFP) set by the on‑shell imaginary part of its self‑energy, ‑Im Σ(E,k); near EF this scales as (E−EF)² (Fermi liquid), with corrections from plasmon emission once E−EF > ωp. The classic universal IMFP curve (Seah and Dench) shows a minimum of a few Å in the 30–100 eV range — precisely the LEED energy regime relevant to the student's existing simulations.

**Key references.**
- J. J. Quinn and R. A. Ferrell, "Electron self‑energy approach to correlation in a degenerate electron gas", Phys. Rev. 112, 812 (1958), DOI 10.1103/PhysRev.112.812 — the foundational Quinn–Ferrell result for the electron lifetime in the high‑density limit, τ ∝ rs^(−5/2)(E−EF)^(−2).
- P. M. Echenique, J. M. Pitarke, E. V. Chulkov and A. Rubio, "Theory of inelastic lifetimes of low‑energy electrons in metals", Chem. Phys. 251, 1 (2000), DOI 10.1016/S0301‑0104(99)00313‑4 — the canonical review covering the FEG, GW, and band‑structure refinements; it is the single best entry point for the student because it joins the dielectric/loss‑function picture (topic 1) to the self‑energy/IMFP picture (this topic).
- I. Campillo, J. M. Pitarke, A. Rubio, E. Zarate and P. M. Echenique, "Inelastic lifetimes of hot electrons in real metals", Phys. Rev. Lett. 83, 2230 (1999), DOI 10.1103/PhysRevLett.83.2230 — first‑principles GW lifetimes in Al, Cu and noble metals.
- W. S. M. Werner, "Electron transport in solids for quantitative surface analysis", Surf. Interface Anal. 31, 141 (2001), DOI 10.1002/sia.973 — practical compilation of IMFPs and elastic mean free paths.
- C. J. Powell and A. Jablonski, *NIST Electron Inelastic‑Mean‑Free‑Path Database*, SRD 71, version 1.2 (NIST, Gaithersburg) — implements the TPP‑2M predictive equation of S. Tanuma, C. J. Powell and D. R. Penn, Surf. Interface Anal. 21, 165 (1994), DOI 10.1002/sia.740210302; revised in H. Shinotsuka, S. Tanuma, C. J. Powell and D. R. Penn, Surf. Interface Anal. 47, 871 (2015), DOI 10.1002/sia.5789. The Penn algorithm, on which TPP‑2M is built, is in D. R. Penn, Phys. Rev. B 35, 482 (1987), DOI 10.1103/PhysRevB.35.482.
- M. P. Seah and W. A. Dench, "Quantitative electron spectroscopy of surfaces: A standard data base for electron inelastic mean free paths in solids", Surf. Interface Anal. 1, 2 (1979), DOI 10.1002/sia.740010103 — origin of the famous "universal curve".
- R. F. Egerton, *Electron Energy‑Loss Spectroscopy in the Electron Microscope*, 3rd ed. (Springer, 2011), DOI 10.1007/978‑1‑4419‑9583‑4 — connects everything to actual EELS data.

**The Bethe ridge and electron–electron scattering kinematics.** For a projectile electron of momentum k₀ scattering at (q,ω), the Bethe ridge ω = q²/2 + q·k₀ marks the locus of energy‑momentum conservation for free‑particle binary collisions; in (q,ω) plots it is the straight line that bounds where the bulk of e–e scattering weight lives. The plasmon line is below (smaller ω at given q) and the particle–hole continuum above. This picture is treated cleanly in Egerton's textbook (chapter 3) and in Mahan's textbook.

**Connection to the project.** When the projectile is an electron wave‑packet, the "stopping power" is conceptually the rate of conversion of its kinetic energy into the loss function of the medium, but with the additional features that exchange should be checked (single‑determinant TDDFT does include it via the xc functional; for a free wave‑packet entering a target this is approximate) and that the wave‑packet has a finite width in (q,ω) so it samples a finite region of the Bethe ridge rather than a δ function. The IMFP of the wave‑packet in jellium can be directly compared to TPP‑2M / Penn algorithm predictions and to the Echenique–Pitarke–Chulkov–Rubio FEG curves.

---

## 7. rt‑TDDFT methodology for stopping power

**Mental model.** Treat the projectile as a classical point charge (or pseudo‑atom) moving on a constant‑velocity Ehrenfest trajectory through the periodic target. Propagate the Kohn–Sham orbitals in real time; the work done against the projectile equals the rate of increase of the KS energy, dE/dt = F·v, where F is the instantaneous Hellmann–Feynman force on the projectile. After a transient (during which the wake is established), dE/dt becomes stationary and gives the electronic stopping S(v). The same simulation also gives the force on the projectile, which fluctuates with the local density along the trajectory and whose path‑averaged component is S(v).

**Key methodological references.**
- J. M. Pruneda, D. Sánchez‑Portal, A. Arnau, J. I. Juaristi and E. Artacho, "Electronic stopping power in LiF from first principles", Phys. Rev. Lett. 99, 235501 (2007), DOI 10.1103/PhysRevLett.99.235501 — the seminal demonstration in an insulator, establishing the threshold behaviour expected from the band gap and showing that Ehrenfest‑TDDFT reproduces it.
- A. Schleife, Y. Kanai and A. A. Correa, "Accurate atomistic first‑principles calculations of electronic stopping power in aluminum", Phys. Rev. B 91, 014306 (2015), DOI 10.1103/PhysRevB.91.014306 — the canonical aluminium benchmark; demonstrates channel/off‑channel differences and the need for explicit core electrons at high v.
- A. A. Correa, "Calculating electronic stopping power in materials from first principles", Comput. Mater. Sci. 150, 291 (2018), DOI 10.1016/j.commatsci.2018.03.064 — the most readable comprehensive review of rt‑TDDFT stopping methodology, written by the lead author of INQ; this should be top of the reading list.
- D. C. Yost, Y. Yao and Y. Kanai, "Examining real‑time TDDFT non‑equilibrium simulations for the calculation of electronic stopping power", Phys. Rev. B 96, 115134 (2017), DOI 10.1103/PhysRevB.96.115134 — exhaustive convergence study (basis, xc, finite‑size, pseudopotentials, core corrections).
- A. Kononov, T. Hentschel, S. B. Hansen and A. D. Baczewski, "Trajectory sampling and finite‑size effects in first‑principles stopping power calculations", npj Comput. Mater. 9, 205 (2023), DOI 10.1038/s41524‑023‑01157‑7 — the authoritative analysis of trajectory averaging, finite‑size scaling, and "ouroboros" self‑interaction artefacts.
- A. Kononov, A. J. White, K. A. Nichols, S. X. Hu and A. D. Baczewski, "Reproducibility of real‑time time‑dependent density functional theory calculations of electronic stopping power in warm dense matter", Phys. Plasmas 31, 043904 (2024), DOI 10.1063/5.0198008 — cross‑code (Qball, SIESTA, Octopus, GPAW, INQ) reproducibility study.
- X. Andrade, C. D. Pemmaraju, A. Kartsev, J. Xiao, A. Lindenberg, S. Rajpurohit, L. Z. Tan, T. Ogitsu and A. A. Correa, "INQ, a modern GPU‑accelerated computational framework for (time‑dependent) density functional theory", J. Chem. Theory Comput. 17, 7447 (2021), DOI 10.1021/acs.jctc.1c00562 — the INQ code paper itself; the student should read this end to end as their methods reference.
- J. Simoni, X. Andrade, W. Fang, A. C. Grieder, A. A. Correa, T. Ogitsu and Y. Ping, "Spin non‑collinear real‑time TDDFT and implementation in INQ", APL Comput. Phys. 1, 026108 (2024), arXiv:2506.21908 — for spin‑aware extensions.

**Other key implementations (for cross‑checking and historical context).**
- A. Schleife, E. W. Draeger, V. M. Anisimov, A. A. Correa and Y. Kanai, "Quantum dynamics simulation of electrons in materials on high‑performance computers", Comput. Sci. Eng. 16, 54 (2014), DOI 10.1109/MCSE.2014.55 — Qball.
- X. Andrade, D. Strubbe, U. De Giovannini et al., "Real‑space grids and the Octopus code as tools for the development of new simulation approaches for electronic systems", Phys. Chem. Chem. Phys. 17, 31371 (2015), DOI 10.1039/C5CP00351B — Octopus.
- A. Ojanperä, V. Havu, L. Lehtovaara and M. Puska, "Nonadiabatic Ehrenfest molecular dynamics within the projector augmented‑wave method", J. Chem. Phys. 136, 144103 (2012), DOI 10.1063/1.3700800 — GPAW.

**On the limitations of ALDA and the xc kernel.** The adiabatic LDA misses memory effects and mistreats the high‑frequency tail; this matters most at v ≫ vF where dE/dt is dominated by deep electron–hole excitations. The Vignale–Kohn current‑dependent functional and the Corradini parametrisation of the static local field (PRB 57, 14569 (1998), DOI 10.1103/PhysRevB.57.14569) are the standard beyond‑ALDA tools. See also G. Vignale and W. Kohn, "Current‑dependent exchange‑correlation potential for dynamical linear response theory", Phys. Rev. Lett. 77, 2037 (1996), DOI 10.1103/PhysRevLett.77.2037.

**Connection to the project.** This topic is the methodological backbone of the thesis. The Correa 2018 review and the Kononov 2023 trajectory paper together give a checklist of every source of systematic and statistical error the student must control.

---

## 8. Beyond jellium — real materials

**Mental model.** Real targets break the homogeneity of jellium in three ways with distinct consequences. (a) A band gap suppresses low‑energy electron–hole pair production — there is a velocity threshold v_th ≈ Eg/(2 kF) for direct excitation across the gap, but below threshold local defects and dynamical states may still funnel electrons across (the "electron elevator" mechanism). (b) Localised d states near EF dramatically enhance low‑velocity stopping in noble and transition metals through their large local density and high transport cross sections. (c) Crystallinity introduces channels along which the projectile sees lower average density; channelled stopping is therefore lower than off‑channel ("random") stopping, an effect controllable in TDDFT only by careful trajectory averaging.

**Key papers.**
- J. M. Pruneda et al., Phys. Rev. Lett. 99, 235501 (2007) — threshold behaviour in LiF.
- A. Lim, W. M. C. Foulkes, A. P. Horsfield, D. R. Mason, A. Schleife, E. W. Draeger and A. A. Correa, "Electron elevator: Excitations across the band gap via a dynamical gap state", Phys. Rev. Lett. 116, 043201 (2016), DOI 10.1103/PhysRevLett.116.043201 — the elevator mechanism in self‑irradiated Si.
- F. Matias, P. L. Grande, N. E. Koval, J. M. B. Shorto, T. F. Silva and N. R. Arista, "Deeper‑band electron contributions to stopping power of silicon for low‑energy ions", Phys. Rev. A 110, 022811 (2024), arXiv:2405.07794 — systematic evidence that "elevator" and "promotion" mechanisms dynamically excite deeper bands.
- M. A. Zeb, J. Kohanoff, D. Sánchez‑Portal, A. Arnau, J. I. Juaristi and E. Artacho, "Electronic stopping power in gold: The role of d electrons and the H/He anomaly", Phys. Rev. Lett. 108, 225504 (2012), DOI 10.1103/PhysRevLett.108.225504 — the H/He anomaly resolved.
- D. C. Yost and Y. Kanai, "Electronic stopping power for protons and α‑particles in silicon carbide", Phys. Rev. B 94, 115107 (2016), DOI 10.1103/PhysRevB.94.115107 — semiconductor benchmark.
- M. Caro, A. A. Correa, E. Artacho and A. Caro, "Coupled electron‑ion dynamics: Application to ion irradiation in metals", Sci. Rep. 7, 2618 (2017), DOI 10.1038/s41598‑017‑02780‑3.
- S. Bubin, B. Wang, S. Pantelides and K. Varga, "Simulation of high‑energy ion collisions with graphene fragments", Phys. Rev. B 85, 235435 (2012), DOI 10.1103/PhysRevB.85.235435.

**Warm dense matter (WDM).**
- A. Kononov, A. J. White, K. A. Nichols, S. X. Hu and A. D. Baczewski, Phys. Plasmas 31, 043904 (2024), DOI 10.1063/5.0198008 — TDDFT cross‑code WDM stopping benchmark.
- W. R. Johnson, J. Nilsen and K. T. Cheng, "Average‑atom treatment of relaxation time in X‑ray Thomson scattering from warm‑dense matter", Phys. Rev. E 86, 036410 (2012), DOI 10.1103/PhysRevE.86.036410 — Mermin DF in finite‑T plasmas.
- N. D. Mermin, "Thermal properties of the inhomogeneous electron gas", Phys. Rev. 137, A1441 (1965), DOI 10.1103/PhysRev.137.A1441 — finite‑T DFT, basis for Mermin–Kohn–Sham at WDM conditions.
- A. Kononov, T. Hentschel, S. B. Hansen and A. D. Baczewski, "Electronic stopping in warm dense matter using Ehrenfest dynamics and TDDFT", OSTI 2004179 (2022) — conference talk slides give a clear pictorial summary.

**Connection to the project.** Understanding these "real‑material" complications clarifies why jellium, despite its simplicity, remains the cleanest test bed for benchmarking the rt‑TDDFT method itself. Once the wave‑packet approach is validated in jellium, each of the above mechanisms (gaps, d electrons, channelling, elevators) constitutes a natural follow‑up question.

---

## 9. Pedagogical and review resources

**Most useful reviews and lectures (curated, in order of value to the project).**
- A. A. Correa, Comput. Mater. Sci. 150, 291 (2018), DOI 10.1016/j.commatsci.2018.03.064 — the cleanest modern overview of rt‑TDDFT for stopping; written for newcomers.
- C. P. Race, D. R. Mason, M. W. Finnis, W. M. C. Foulkes, A. P. Horsfield and A. P. Sutton, Rep. Prog. Phys. 73, 116501 (2010), DOI 10.1088/0034‑4885/73/11/116501 — the most pedagogical broad review on electronic excitations in radiation damage; covers Ehrenfest, friction, and atomistic models; written by the Imperial group whose papers form a coherent literature thread the student should master.
- P. M. Echenique, J. M. Pitarke, E. V. Chulkov and A. Rubio, Chem. Phys. 251, 1 (2000) — the classic review of inelastic lifetimes that bridges dielectric and self‑energy formalisms.
- J. M. Pitarke, V. M. Silkin, E. V. Chulkov and P. M. Echenique, "Theory of surface plasmons and surface‑plasmon polaritons", Rep. Prog. Phys. 70, 1 (2007), DOI 10.1088/0034‑4885/70/1/R01 — the standard reference on surface vs. bulk plasmons.
- ICRU Report 49, *Stopping Powers and Ranges for Protons and Alpha Particles* (1993) and ICRU Report 73, *Stopping of Ions Heavier than Helium* (2005, with errata in J. ICRU 5(1) 2005) — the standard data compilations and methodological recommendations used as benchmarks throughout the field.
- IAEA Stopping Power Database (https://nds.iaea.org/stopping/) — searchable experimental compilation maintained by the IAEA Nuclear Data Section, the most up‑to‑date public archive of experimental S(v) curves.
- NIST PSTAR / ASTAR / ESTAR online databases (https://physics.nist.gov/PhysRefData/Star/Text/intro.html) — continuous‑slowing‑down approximation S(v) for protons, α and electrons in standard materials.
- SRIM/TRIM software and book by J. F. Ziegler, J. P. Biersack and M. D. Ziegler, *SRIM – The Stopping and Range of Ions in Matter* (Lulu, 2008), ISBN 978‑0‑9654207‑1‑6 — the de facto industrial standard, with all its limitations.

**Textbooks (essentials for a thesis‑level command).**
- P. Sigmund, *Particle Penetration and Radiation Effects*, Vols. 1 and 2 (Springer, 2006/2014).
- G. F. Giuliani and G. Vignale, *Quantum Theory of the Electron Liquid* (Cambridge, 2005).
- G. D. Mahan, *Many‑Particle Physics* (Plenum, 2000).
- D. Pines and P. Nozières, *Theory of Quantum Liquids* (Westview, reprint 1999).
- C. A. Ullrich, *Time‑Dependent Density‑Functional Theory: Concepts and Applications* (Oxford, 2012), ISBN 978‑0‑19‑956302‑9 — the most pedagogical TDDFT textbook, including a chapter on real‑time methods and the f‑sum rule.
- M. A. L. Marques, N. T. Maitra, F. M. S. Nogueira, E. K. U. Gross and A. Rubio (eds.), *Fundamentals of Time‑Dependent Density Functional Theory*, Lect. Notes Phys. 837 (Springer, 2012), DOI 10.1007/978‑3‑642‑23518‑4 — collects the seminal lectures on rt‑TDDFT.
- R. F. Egerton, *EELS in the Electron Microscope*, 3rd ed. (Springer, 2011).
- W. Eckstein, *Computer Simulation of Ion–Solid Interactions*, Springer (1991), ISBN 978‑3‑642‑73513‑4 — particularly clear for the binary‑collision picture and MC framework underlying SRIM.

**Lecture notes and summer school proceedings.**
- The Erice and Trieste summer school proceedings on "Interaction of Charged Particles with Solids and Surfaces" (NATO ASI series; the 1990 Alicante volume edited by A. Gras‑Marti, H. M. Urbassek, N. R. Arista and F. Flores is especially clear on the dielectric formalism applied to stopping).
- Lectures on rt‑TDDFT by X. Andrade, A. Castro and A. Rubio at the Octopus/INQ tutorials available on the INQ project website (https://gitlab.com/npneq/inq) and at https://octopus-code.org — the most direct way for the student to gain practical fluency.

---

## 10. Quasiparticle picture and spectral functions

**Mental model.** The interacting Green's function G(k,ω) develops a complex self‑energy Σ(k,ω) = Re Σ + i Im Σ. The pole structure of G defines quasiparticles: the real part of Σ shifts the dispersion, the imaginary part broadens it (giving the inelastic lifetime τ = ‑1/(2 Im Σ)). For a Fermi liquid, Im Σ ∝ (E‑EF)² near the Fermi surface, recovering the Quinn–Ferrell result. The spectral function A(k,ω) = (1/π)|Im G| also develops plasmonic satellites at energies ω_qp ± ωp; the cumulant expansion is the natural language for these. The IMFP of an electron at energy E is precisely ℏv(E)/(2|Im Σ(E,k(E))|), tying this section back to topic 6.

**Key references.**
- L. Hedin, "New method for calculating the one‑particle Green's function with application to the electron‑gas problem", Phys. Rev. 139, A796 (1965), DOI 10.1103/PhysRev.139.A796 — birth of the GW approximation.
- L. Hedin and S. Lundqvist, in *Solid State Physics*, vol. 23, eds. F. Seitz, D. Turnbull and H. Ehrenreich (Academic, 1969) — long pedagogical review of Σ and satellites.
- F. Aryasetiawan and O. Gunnarsson, "The GW method", Rep. Prog. Phys. 61, 237 (1998), DOI 10.1088/0034‑4885/61/3/002 — the standard review.
- B. Holm and U. von Barth, "Fully self‑consistent GW self‑energy of the electron gas", Phys. Rev. B 57, 2108 (1998), DOI 10.1103/PhysRevB.57.2108.
- F. Aryasetiawan, L. Hedin and K. Karlsson, "Multiple plasmon satellites in Na and Al spectral functions from ab initio cumulant expansion", Phys. Rev. Lett. 77, 2268 (1996), DOI 10.1103/PhysRevLett.77.2268 — cumulant expansion and plasmonic satellites.
- M. Guzzo, G. Lani, F. Sottile, P. Romaniello, M. Gatti, J. J. Kas, J. J. Rehr, M. G. Silly, F. Sirotti and L. Reining, "Valence electron photoemission spectrum of semiconductors: Ab initio description of multiple satellites", Phys. Rev. Lett. 107, 166401 (2011), DOI 10.1103/PhysRevLett.107.166401 — modern application to real materials.
- P. Nozières, *Theory of Interacting Fermi Systems*, Benjamin (1964) / Westview reprint — the canonical Fermi‑liquid textbook.

**Connection to the project.** When the projectile is itself an electron, its TDDFT propagation effectively realises (within the adiabatic xc approximation) the GW‑like decay process of a hot electron in the medium, with τ extractable directly from the wave‑packet's amplitude decay. The spectral function language gives the cleanest interpretation of how rt‑TDDFT redistributes the projectile's spectral weight into plasmonic and particle‑hole channels.

---

## 11. Finite jellium systems and clusters

**Mental model.** A jellium sphere of radius R has volume plasmon at ωp and a surface (Mie) plasmon at ωp/√3 in the dipolar limit; multipolar surface plasmons appear at frequencies ωp √[ℓ/(2ℓ+1)]. Quantum size effects produce shell closures (magic numbers 2, 8, 20, 40, 58, 92, …) that distort the Mie peak through Landau damping into single particle‑hole transitions. For a wave‑packet impinging on a small jellium sphere, the dominant excitation channels are the surface plasmon for grazing trajectories and the volume plasmon for transmission through the centre.

**Key references.**
- M. Brack, "The physics of simple metal clusters: Self‑consistent jellium model and semiclassical approaches", Rev. Mod. Phys. 65, 677 (1993), DOI 10.1103/RevModPhys.65.677 — the standard pedagogical review.
- W. Ekardt, "Size‑dependent photoabsorption and photoemission of small metal particles", Phys. Rev. B 31, 6360 (1985), DOI 10.1103/PhysRevB.31.6360 — RPA in spherical jellium.
- W. A. de Heer, "The physics of simple metal clusters: Experimental aspects and simple models", Rev. Mod. Phys. 65, 611 (1993), DOI 10.1103/RevModPhys.65.611 — companion experimental review.
- G. F. Bertsch and R. A. Broglia, *Oscillations in Finite Quantum Systems*, Cambridge (1994), ISBN 978‑0‑521‑41148‑3 — collective modes in finite systems.
- C. Yannouleas and R. A. Broglia, "Landau damping and wall dissipation in large metal clusters", Ann. Phys. 217, 105 (1992), DOI 10.1016/0003‑4916(92)90340‑R.
- F. Calvayrac, P.‑G. Reinhard, E. Suraud and C. A. Ullrich, "Nonlinear electron dynamics in metal clusters", Phys. Rep. 337, 493 (2000), DOI 10.1016/S0370‑1573(00)00043‑0 — comprehensive TDDFT review for clusters; close in spirit to what the student is doing for nanoflakes.

**Connection to the project.** Because the student has already done LEED‑style scattering from finite nanoflakes, the cluster/jellium‑sphere literature is the natural bridge: replacing a graphene flake by a jellium sphere converts the simulation from elastic‑Bragg‑dominated to inelastic‑plasmon‑dominated, and the Mie formula gives an analytic check on the dominant feature in the loss spectrum.

---

## 12. Connections to experimental techniques

**Mental model.** Every theoretical loss function ‑Im[1/ε(q,ω)] is, modulo geometry factors, exactly what is measured by an EELS or X‑ray inelastic experiment. Stopping power S(v) is the first frequency moment of the loss function weighted by 1/ω, integrated over an appropriate (q,ω) window: S(v) = (2/πv²) ∫₀^(qmax) (dq/q) ∫₀^(qv) ω Im[‑1/ε(q,ω)] dω. RBS measures the energy of backscattered ions and inverts the kinematic equations to get S(E). Thin‑foil transmission and time‑of‑flight directly measure ΔE versus thickness and yield S(v) over a wide range. Mean excitation energies I are extracted from inverse‑v² fits to the Bethe formula in the high‑v regime.

**Key references.**
- R. F. Egerton, *EELS in the Electron Microscope*, 3rd ed. (Springer, 2011) — chapters 3–4 are the canonical pedagogical link between theory and EELS measurement.
- H. Raether, *Excitation of Plasmons and Interband Transitions by Electrons*, Springer Tracts in Modern Physics 88 (Springer, 1980), ISBN 978‑3‑540‑09677‑8 — historical and pedagogical; the standard early reference for the connection between EELS and Im[‑1/ε].
- F. J. García de Abajo, "Optical excitations in electron microscopy", Rev. Mod. Phys. 82, 209 (2010), DOI 10.1103/RevModPhys.82.209 — modern review of relativistic EELS theory and applications including STEM.
- W.‑K. Chu, J. W. Mayer and M.‑A. Nicolet, *Backscattering Spectrometry*, Academic Press (1978), ISBN 978‑0‑12‑173850‑1 — the standard RBS text.
- H. Paul, "A comparison of recent stopping‑power tables for light and medium‑heavy ions with experimental data, and applications to radiotherapy dosimetry", Nucl. Instrum. Methods B 247, 166 (2006), DOI 10.1016/j.nimb.2006.01.059 — assesses how well stopping tables agree with experiment.
- D. Primetzhofer and collaborators have produced a long series of high‑precision low‑velocity stopping measurements; representative paper: D. Roth, B. Bruckner, M. V. Moro, S. Gruber, D. Goebl, J. I. Juaristi, M. Alducin, R. Steinberger, J. Duchoslav, D. Primetzhofer and P. Bauer, "Electronic stopping of slow protons in transition and rare earth metals: Breakdown of the free electron gas concept", Phys. Rev. Lett. 118, 103401 (2017), DOI 10.1103/PhysRevLett.118.103401.
- For the EELS−loss‑function workflow in a TDDFT context, see V. U. Nazarov, J. M. Pitarke, C. S. Kim and Y. Takada, "Time‑dependent density functional theory of the dynamical response of an inhomogeneous electron gas", Phys. Rev. B 71, 121106(R) (2005), DOI 10.1103/PhysRevB.71.121106.

**Connection to the project.** Once the wave‑packet approach in INQ is fully validated, the loss function extracted from a uniform‑kick run on the same jellium can be directly compared to (a) the analytical Lindhard/Mermin curves and (b) tabulated EELS data for nearly‑free‑electron metals (Al, Mg, Na). The same pipeline yields S(v) via the integral above, providing a fully consistent cross‑check between dielectric‑formalism and Ehrenfest‑projectile predictions of the same simulation.

---

## A suggested reading order for a beginner

1. Start with Correa, *Comput. Mater. Sci.* 150, 291 (2018) for the rt‑TDDFT methodological picture, paired with the INQ paper (Andrade et al., JCTC 2021).
2. Read Sigmund Vol. 1 chapters 1–6 for the classical/Bethe/Bohr/Lindhard theory, with Jackson chapter 13 as a parallel text.
3. Work through chapters 4–5 of Giuliani–Vignale to internalise χ₀, ε, the loss function, and the f‑sum rule.
4. Read Echenique, Pitarke, Chulkov and Rubio, *Chem. Phys.* 251, 1 (2000) to see lifetimes and dielectric stopping unified.
5. Read Race et al., *Rep. Prog. Phys.* 73, 116501 (2010) for the radiation‑damage view that puts everything into application context.
6. Read the trio of jellium‑sphere papers (Ekardt 1985; Brack 1993; Calvayrac et al. 2000) to bridge from infinite jellium to the nanoflake/cluster geometries the student already simulates.
7. Read the canonical first‑principles stopping papers in chronological order: Pruneda 2007 (LiF), Schleife 2015 (Al), Zeb 2012 (Au), Lim 2016 (Si), Yost–Kanai 2017 (Si convergence), Kononov 2023 (trajectories), Kononov 2024 (WDM reproducibility).
8. Use Egerton's EELS textbook and the IAEA/NIST databases as continuous reference companions throughout.

This reading list should give a well‑prepared beginner the conceptual vocabulary, the historical thread, and the practical methodological background needed to connect rt‑TDDFT wave‑packet simulations in jellium to the broader and well‑developed field of charged‑particle stopping in matter, and to identify clean, publishable next steps from the LEED‑pattern reproduction the student has already achieved.