# Plan: Report 1 Plots Remake

**Created:** 2026-05-27
**Spec:** `docs/reports/report1/drafts/draft3/plots_remake_resolved.md`
**Status:** In progress — grilling complete, production starting

---

## Phase 0: Style infrastructure

- [ ] Update `_shared_style.py` with tuneable `STYLE_CONFIG` dict
- [ ] Switch LaTeX preamble to `lmodern` + `mlmodern`
- [ ] Verify fonts render correctly with a test plot

## Phase 1: Coronene plots (collaborative, one at a time)

- [ ] Plot 1: `fig_coronene_target.png` — cyan crosshair
- [ ] Plot 2: `fig_leed_backscatter_centre.png`
- [ ] Plot 3: `fig_leed_transmission.png` — cyan FFT overlay + legend
- [ ] Plot 4: `fig_leed_backscatter_ccbond.png` (appendix)
- [ ] Plot 5: `fig_leed_validation.png` — 4-panel remake

## Phase 2: Jellium density plots (collaborative)

- [ ] Plot 6–9: `fig_density_diff_2d_t{1–4}.png`
- [ ] Plot 10–13: `fig_density_profile_t{1–4}.png`

## Phase 3: Jellium energetics (collaborative)

- [ ] Plot 14: `fig_energy_decomp_system.png`
- [ ] Plot 15: `fig_energy_decomp_wp.png`
- [ ] Plot 16: `fig_gs_decomposition.png`

## Phase 4: Jellium momentum (collaborative)

- [ ] Plot 17: `fig_momentum_1d.png`
- [ ] Plot 18: `fig_momentum_2d_before.png`
- [ ] Plot 19: `fig_momentum_2d_after.png`
- [ ] Plot 20: `fig_momentum_2d_diff.png`

## Phase 5: Jellium mechanisms (collaborative)

- [ ] Plot 21: `fig_plasmon_fft.png`
- [ ] Plot 22: `fig_loss_function.png`

## Phase 6: Stopping power summary

- [ ] Plot 23: `fig_master_stopping_lowdens.png`
- [ ] Plot 24: `fig_master_stopping_highdens.png`
- [ ] Plot 25–28: `fig_stopping_def_S{1–4}.png`

## Phase 7: High-density duplicates (deferred)

- Deferred. Canonical run: `run_wp_n162_L30_E50_highdens_sigma1_v2`

---

## Key data paths

| Data | Path |
|------|------|
| E25 WP run | `ResearchProject/systems/jellium/run_wp_n162_L50_E25_sigma1_v2/` |
| E25 free run | `ResearchProject/systems/jellium/run_free_wp_L50_E25_sigma1_v2/` |
| E25 classical | `ResearchProject/systems/jellium/run_classical_n162_L50_E25/` |
| Coronene centre | `ResearchProject/systems/coronene/run_propagate_paper_replica/` |
| Coronene CC bond | `ResearchProject/systems/coronene/run_cc_bond/` |
| Coronene GS | `ResearchProject/systems/coronene/run_save_gs_paper_replica/` |
| Output dir | `docs/reports/report1/drafts/draft3/figures/` |
| Style module | `inq-stack/python/inqview/report1/_shared_style.py` |
| Stopping data | `inq-stack/python/inqview/report1/stopping_power_data.py` |
