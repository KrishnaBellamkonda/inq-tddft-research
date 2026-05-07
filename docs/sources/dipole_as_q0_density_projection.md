# `dipole_x` as the q→0 longitudinal density projection

## Statement

For an impulsive ionic kick along x (the QBall-style quantum kick used in
the Li 54-atom runs), the **electronic dipole moment**

$$
d_x(t) \;=\; \int x \,\rho(\mathbf r, t)\,d^3 r
$$

is mathematically identical to the **q→0 limit of the longitudinal
density-density response function** along x. This is why, in the
plasmon-vs-electron-hole diagnostic plan for these runs, `dipole_x` (a
column already present in `observables.csv`) is the correct cheap proxy
for the bulk plasmon mode — and **no new C++ observable is needed**.

## Identity

For a translationally-invariant electronic system perturbed by an
external scalar potential $V_{\rm ext}(\mathbf r, t)$, linear response
gives

$$
\delta\rho(\mathbf q, \omega) \;=\; \chi(\mathbf q, \omega)\,
                                    V_{\rm ext}(\mathbf q, \omega).
$$

The induced dipole along x in the long-wavelength limit is

$$
\delta d_x(\omega) \;=\; \int x\,\delta\rho(\mathbf r, \omega)\,d^3 r
                  \;=\; -i\,\partial_{q_x}\,\delta\tilde\rho(\mathbf q,\omega)
                        \Big|_{\mathbf q\to 0}.
$$

So the **time-Fourier transform** of $d_x(t)$ probes
$\chi(\mathbf q\to 0, \omega)$ along the kick axis. A pole of $\chi$ at
$\omega = \omega_p$ — the bulk plasmon — appears as a peak in
$|\widehat d_x(\omega)|^2$.

This is the same identity that makes the dipole spectrum the standard
TDDFT route to the optical absorption spectrum (see Yabana & Bertsch
1996, Marques et al. *TDDFT* textbook 2012, Ullrich *TDDFT* 2012). For a
kicked metal in the low-velocity regime, the same identity says the
plasmon dispersion at q=0 (the bulk plasmon) shows up directly in the
dipole-x FFT.

## Why no new C++ observable is needed

The inq `real_time::propagate(...)` call already supports
`.observables_dipole()`, and the inqkit `ObservablesWriter` / `StepContext`
already plumb the `dipole_x, dipole_y, dipole_z` triple to
`observables.csv`. The Li run template
(`QuantumKickExtension/inq-codebase/Li/shared/cpp/run_template.hpp`)
sets `sel.dipole_x = sel.dipole_y = sel.dipole_z = true`. So the data is
already collected.

In the post-processing, the inqview `_extended_spectra` pipeline
(`inq-stack/python/inqview/postprocess/observables.py`) applies the same
plateau-detrend / Hann / 8x-zero-pad / FFT recipe to `dipole_x` (and
`dipole_z`) as it does to `energy_total`. The dipole spectrum, alongside
the energy spectrum and the gamma-transitions histogram, is what the
plasmon-vs-electron-hole diagnostic compares.

## Caveat: not the *only* signal of the plasmon

The energy-FFT peak is bilinear in the perturbation
($\Delta E \propto |\delta\rho|^2$ at leading order); the dipole-FFT
peak is linear ($\delta d \propto \delta\rho$). Both should peak at the
same $\omega_p$ for a linear-response perturbation, but at high kick
velocity (the paper's "softening" / "hardening" regime) the perturbation
is no longer linear and the two peaks can split. The Li 54-atom runs at
v = 0.0626 a.u. (top of low-v) should still be in the linear regime;
v = 0.300 a.u. is at the entrance to the softening branch and may show a
small split between $|\widehat{\Delta E}|^2$ and $|\widehat d_x|^2$ peaks.

## References

- *Time-Dependent Density Functional Theory*, Marques et al. (eds.),
  Springer Lecture Notes in Physics 837, 2012 — chapter on TDDFT in the
  linear-response regime, dipole-spectrum identity.
- C. Ullrich, *Time-Dependent Density-Functional Theory: Concepts and
  Applications*, Oxford University Press, 2012 — §4.2 derives the
  dipole-FFT / absorption-spectrum identity.
- Yabana & Bertsch, *Phys. Rev. B* 54, 4484 (1996) — the original
  real-time TDDFT impulse-perturbation formulation that this kick
  paradigm extends.
- Santervás-Arranz, Stengel, Artacho, *Phys. Rev. Research* 7, 033292
  (2025), and the BCN:1719P quantum-kick paper (Cavendish, 2026) —
  application to the metallic Li kick paradigm; Figures 4-5.
- Giuliani & Vignale, *Quantum Theory of the Electron Liquid*, Cambridge
  2005 — §7 derives $\chi(\mathbf q \to 0, \omega)$ for the homogeneous
  electron gas, including the bulk plasmon dispersion that the Li
  6.5 eV peak quantitatively matches.
