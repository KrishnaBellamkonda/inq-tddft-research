# Plan: draft3 free-WP figures, KS-energy stopping sweep, equation PNGs

Output folder: `docs/presentations/assets/draft3_freewp/` (own `INDEX.md`).
All Python via `/local/data/public/skcb2/tddft/venv/bin/python3`. All figures `.png`
(plus GIFs where an animation is requested). Report plot rules
(`inqview.report1._shared_style`, TufteCritic, fixed dims / no tight-bbox pitfalls).

## Source data (verified present)

- Free WP: `ResearchProject/systems/jellium/run_free_wp_L50_E100_dens`
  (sigma=5, E=100 eV, L=50 cubic, dt=0.02, write_every=2). Has density VTI series
  (`density_total`, `density_wp`, `wavefunction_wp`, 232 frames) and
  `wp_momentum_stats.csv`, `wp_real_space_stats.csv`.
- Interference-free window (IFW): t <= 3.5 a.u. -> step <= 175 -> frames
  `density_t000000` .. `density_t000174` (even steps, ~88 frames). Established in A1.
- Stopping sweep (L=50, r_s=5.69, KS-energy def S2): collected by
  `inqview.report1.stopping_power_data.collect_L50_data()`:
  - sigma=1 full velocity sweep: E20/25/50/100/200/300 (`*_sigma1_v2`).
  - sigma=5 full velocity sweep: E20/25/50/100/300/600 (unsuffixed base runs).
  - sigma=0.5/3/8: single E100 points (`*_sigma{0p5,3,8}`).

## Deliverables

1. `render_freewp_equations.py` -> 5 transparent-bg equation PNGs (extend
   `docs/presentations/drafts/render_equations.py` style; usetex->mathtext fallback):
   - `eq_spreading.png`  sigma(t)=sigma0 sqrt(1+(hbar t/2 m sigma0^2)^2)
   - `eq_bethe.png`      S(v)=(4 pi n Z^2 e^4/m_e v^2) ln(2 m_e v^2/I)
   - `eq_bloch.png`      S(v)=(4 pi n Z^2 e^4/m_e v^2)[ln(2 m_e v^2/I)+psi(1)-Re psi(1+i Z e^2/hbar v)]
   - `eq_lindhard.png`   S(v)=(2/pi v^2) int dq/q int_0^{qv} omega Im[-1/eps(q,omega)] domega
   - `eq_ks_stopping.png` S(v)=-d<T_WP>/dz, <T_WP>=<psi_WP|-1/2 grad^2|psi_WP>
   LHS uses S(v) (user). RHS explicit constants (default; switchable to a.u.).
   Ground: Bethe (1930), Bloch (1933), Lindhard (1954). Record in docs/sources.

2. `freewp_sigma_t.py` -> standalone sigma(t) over IFW (density width
   sqrt(sigma_z2); analytic sqrt(sigma0^2/2 + t^2/(2 sigma0^2)), sigma0=5).
   Known-case: sigma_r(0) ~ sigma0/sqrt2 = 3.5355.

3. `freewp_energy.py` -> standalone absolute total energy E(t)=E_kin(t) over IFW
   (Hartree=XC=0). Annotate max |Delta E| drift. Known-case: E(0)~100 eV.

4. `freewp_xz_density_gif.py` -> xz mid-slice (y=L/2) total-density heatmap GIF over
   IFW. Fixed colorbar across frames (log), per shared-colorbar rule. Known-case:
   slab-integral of one frame vs total norm ~ 1.

4b. RECREATE coronene `xz_density.gif` with propagation axis horizontal:
   edit `docs/presentations/drafts/make_coronene_anim_gifs.py` `build_xz()` so the
   slab is transposed -> z on the x-axis (length-wise, left->right), x on the y-axis
   (height-wise). xlabel "z (Bohr)", ylabel "x (Bohr)". Same orientation schema used
   for the free-WP xz gif (deliverable 4). Rationale: left-to-right reading of the
   projectile path (user, 2026-06-02). Overwrite in place (referenced by the decks).

5. `freewp_3d_gif.py` -> pvbatch volume-render GIF of free-WP density over IFW
   (reuse `inqview.paraview.ParaViewPipeline`, log scalar map, clean white bg,
   same camera idea as `render_coronene_total3d_gif.py`). density_total series.

6. `stopping_KSenergy_sigmasweep.py` -> adapt
   `docs/reports/report1/drafts/draft5/scripts/make_fig_master_stopping_vaxis.py`:
   REMOVE red "loss function (simulation)" and grey "analytical Lindhard (box q)".
   KEEP classical Ehrenfest. WP sigma=0.5/1/3/5/8 (S2_eV_per_bohr). Legend ascending
   0.5,1,3,5,8 (+ Classical). Master marker scheme. Known-case: S2 ranges finite>0.

## Validation (per testing rule)

- Equation PNGs: smoke test (size>50px, alpha, nonzero pixels) as in render harness.
- sigma(t): sigma_r(0) ~ 3.5355 assert.
- energy: E(0) ~ 100 eV; report drift.
- xz gif: per-frame norm sanity; fixed clim verified.
- 3D gif: frame count == IFW frame count; nonzero rendered frames.
- stopping: assert each WP series nonempty, S2 finite and >0; classical present.

No new simulations. Nothing destructive. User previews PNGs (no agent preview).
