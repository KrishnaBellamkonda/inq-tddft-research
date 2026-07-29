**Where this run sits.** This is the **highest-velocity point** of the Phase-5
quantum stopping-power velocity sweep for the localised jellium slab
(r$_s$ ≈ 5.67, σ$_{\rm WP}$ = 0.5, two-sided CAP η = −0.7 Ha, L$_z$ = 25 Bohr):
**v₀ = 6.0 a.u., k₀ = 6, drift energy E = ½k₀²·27.211 ≈ 490 eV**, launched at
z₀ = −23.75 Bohr. The sweep measures the quantum (wavepacket) electronic stopping
power S(E) = [E$_{\rm total}$(t$_f$) − E$_{\rm GS}$]/L$_z$ at
v ∈ {1.3, 2.0, 3.0, 4.0, 5.0, 6.0} and compares it with the matched bulk classical
(σ$_{\rm WP}$=0.5) curve and point-charge Lindhard (see the sweep study notebook
`qsp_phase5_study.ipynb`).

**Aim of this run.** Extend the S(E) curve to the high-velocity end and test whether
the low-v quantum enhancement (S ≈ 2.4–2.6 eV/Bohr, several × the classical/Lindhard
values) persists at v₀ = 6, or turns over toward the Bethe/Lindhard high-velocity
decline.

**Hypothesis.** At v₀ = 6 the Gaussian wavepacket's momentum content
(σ$_p$ = 1/(2σ$_{\rm WP}$) = 1.0, centred at k₀ = 6) pushes a large fraction of the
packet past the real-space **grid Nyquist** wavevector k$_{\rm Nyq}$ = π/h = 6.28
a.u. (single-particle cutoff E$_{\rm cut}$ = ½(π/h)² = 537 eV). The prediction — to
be checked by the diagnostics below — is that this run is **grid-aliased**: the
supra-Nyquist tail (≈ 39 % of the norm) wraps to spurious momenta and *injects*
energy, so E$_{\rm total}$(t) should **rise** (positive late slope) and the extracted
S is an artificially large **lower bound**, not physical stopping. If confirmed, v₀ = 6
must be **excluded** from the physical S(E) curve and re-run on a finer grid
(h ≤ 0.35 Bohr, new GS). The decisive tests are the **momentum panel**
(n$_{\rm wp}$(k$_z$) at t=0 vs k$_{\rm Nyq}$) and the **energy panel** (rising
E$_{\rm total}$).

*Grounding: k₀, energy, CAP, launch z from `run_summary.txt`; the Nyquist / aliasing
figures (537 eV cutoff, ≈39 % supra-Nyquist fraction) are the sweep's grid-resolution
analysis in `qsp_phase5_study.ipynb`. This run is expected to be the aliased endpoint,
not a trusted physical point.*
