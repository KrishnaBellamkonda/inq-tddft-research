# Future TODOs

## Probability current density on PlaneScreens

The scalar LEED density pattern (∫ ρ(r,t) dt accumulated on a 2D plane) loses
directionality information. The vector **probability current density**:

  **J(r,t) = (ℏ/2mi) [ψ* ∇ψ − ψ ∇ψ*]**

could instead be measured on-the-fly at each real-time step and accumulated on
PlaneScreens analogously to the scalar density. This would capture the direction
of electron flow through the screen and give richer information than the scalar
LEED pattern — for example, distinguishing forward transmission from backward
reflection that happen to overlap spatially.

Requires evaluating the gradient of each KS orbital on the GPU at each screen
z-position. See `inq/src/operations/gradient.hpp` for the relevant operator.
Not implemented; post-hoc analysis from saved orbitals is not currently feasible
due to the large data volume.

## Post-hoc LEED analysis from saved density series

An alternative to on-the-fly LEED accumulation is to reconstruct the 2D plane
density from the saved `density_rt_total` frames. This would allow choosing
different screen positions and time windows after the run is complete, without
rerunning. The tradeoff is that the full 3D density must be saved at every
WRITE_EVERY step, which is already done for the current runs. The post-hoc
extraction would need a Python utility: load each frame, slice at the desired
z-index, integrate over the time window.

## Systematic loop-back check for wide-σ runs

run_05_wide_sigma has WP_CZ = 5×σ ≈ 18.9 bohr (close to the cell midpoint).
This means the WP starts nearly at the centre, travels only ~20 bohr before
hitting the far boundary, and the loop-back time is very tight. Future wide-σ
runs should use either a larger cell or a reduced initial position (e.g. 3×σ)
to allow more travel distance. See plan `docs/plans/jellium-wp-rt-initial-exploration.md`.

## Convergence tests

- Grid spacing: current runs use 0.50 bohr. Convergence in LEED pattern vs 0.35 bohr not checked.
- Energy cutoff equivalent for jellium: only affects WP injection quality, not GS density.
- Time step dt=0.02 a.u.: ALDA propagation. Not tested against dt=0.01.
