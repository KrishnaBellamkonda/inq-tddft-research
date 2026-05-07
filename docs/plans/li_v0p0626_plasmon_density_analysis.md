# Plan: identify the bulk Li plasmon mode in the v=0.0626 a.u. density VTI series

**Status:** draft — awaiting user review.
**Target run:** `QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0626_xyz/`
**Linked entry:** `docs/journals/quantumkickextension/2026-05-06_run_propagate_v0p0626_xyz.md`
**Author of plan:** generated 2026-05-06; will be edited based on user review before execution.

---

## 1. Scientific goal

Identify the bulk Li plasmon (DFT-RPA value 6.56 eV; observed near 6.5 eV
across the entire low-v family in the BCN:1719P paper / Santervás-Arranz
et al. *PRR* 7 033292) in the **spatial-temporal evolution** of the
electron density during the v=0.0626 a.u. impulsive kick run. The
existing `dipole_x` FFT gives the q→0 longitudinal-density projection of
the plasmon as a *single scalar* (see
`docs/sources/dipole_as_q0_density_projection.md`). The pixel-by-pixel
FFT analysis adds a **spatial map** — *where* in the cell the plasmon
oscillation is concentrated — and is the diagnostic the user requested
based on a manual inspection of ~70 delta-density VTI frames that
yielded a candidate peak in the 7–8 eV range.

## 2. What we have to work with

After the run completes (~25 h from launch on 2026-05-05 22:03 BST),
the run directory will contain:

| Path | Frames | Cadence | Size estimate |
|---|---:|---|---|
| `results/raw/vti/density_rt_total/` | 156 | every 100 prop steps (Δt = 0.097 fs) | ~13 MB/frame × 156 = **~2 GB total** |
| `results/raw/vti/density_rt_delta/` | 156 | same | same |
| `results/raw/vti/density_rt_delta_coarse/` | 156 | same | ~25 MB total (3 Bohr binned) |
| `results/raw/observables/observables.csv` | 15500 | every step | ~5 MB |
| `results/raw/observables/state_energies.csv` | 1550 × 808 | every 10 steps × n_states | ~100 MB |
| `results/raw/observables/eigenvalues/eigenvalues.csv` | 1 | mirrored from GS | ~50 KB |

User confirmed (Q1) that loading all 156 full-grid VTIs into memory at
once is OK.

## 3. Method overview

The analysis is staged. We start with the cheapest, most direct
approach (lab-frame on raw delta-density). If the plasmon is not
clearly visible because the rigid ion-translation component dominates,
we escalate to the comoving (RF2) frame, exposing a clean longitudinal
electronic-mode signal. At every stage we cross-check against three
independent diagnostics that the run already produces.

### Stage A — lab-frame pixel-by-pixel FFT of `density_rt_delta`

For every grid point `r` on the 54³ grid:

1. Read `δρ(r, t_i)` for i = 0..155 from
   `density_rt_delta/density_delta_t<step>.vti`.
2. Apply the **plateau-detrend** variant (subtract the second-half mean,
   matching the QBall recipe and `state_energy_spectra` postprocess).
   No interpolation — temporal samples used as written.
3. Hann window + 8× zero-pad + `|rfft|²` (same recipe as
   `inqview.postprocess.observables`).
4. Map the resulting per-pixel `|FFT|²(r, ω)` into:
   - **(ii) frequency-averaged 2D maps**: integrate over ω-bins
     {0–2 eV, 2–4 eV, 4–6 eV, **6–7 eV** (plasmon target),
     7–9 eV, 9–13 eV, 13–20 eV} and project the 3D power maps onto
     three 2D slices (xy, xz, yz). One PNG per (band × plane) =
     7 × 3 = 21 PNGs. Headline figure: the `6–7 eV` map at all three
     planes side-by-side.
   - **(i) full 3D isosurface VTI for the 5–8 eV band only** — single
     VTI per peak, opens in ParaView. Lets the user rotate the cell
     and identify the spatial mode (e.g. uniform vs nodal patterns
     that distinguish bulk plasmon from surface modes).
5. **Cell-averaged 1D spectrum**: integrate `|FFT|²(r, ω)` over r → one
   curve `S(ω)`. Look for peak at 6.5 eV. Sanity-check against the
   `dipole_x` FFT we already have.

**Expected outcome — case A1.** If the plasmon is dominant over ion
translation in `δρ`, the 6–7 eV band map shows uniform power across
the cell (consistent with a bulk longitudinal mode at q→0; the
"uniform" pattern is the spatial signature of low-q). The cell-averaged
S(ω) has a peak at ~6.5 eV. **Done — proceed to validation cross-checks
in Stage C.**

**Expected outcome — case A2.** If the rigid ion-translation dominates,
the spectrum is overwhelmed by sub-1 eV power and the 6–7 eV band map
is noisy. The diagnostic of A2 is: integrated power in [0, 1 eV] is >
10× the integrated power in [5, 8 eV]. **Escalate to Stage B.**

### Stage B — comoving frame (Galilean shift via Fourier-shift) — **gated on Stage A2**

For every frame `i = 0..155`:
1. Read the *total* density `ρ_lab(r, t_i)` (not the delta — the delta
   is by definition centred at t=0 ions; we want the whole density to
   shift cleanly with the moving ions in the lab frame).
2. Compute the cumulative ion displacement `Δr_i = v · t_i` where
   `v = (0.0626, 0, 0)` a.u.
3. **Fourier-shift** by `+Δr_i` (note the sign — we want `ρ_RF2(r) =
   ρ_lab(r + v·t)` so we Fourier-shift by `+v·t`):
   ```python
   rho_k    = np.fft.fftn(rho_lab[i])
   ramp     = np.exp(-1j * (KX*dr[0] + KY*dr[1] + KZ*dr[2]))
   rho_rf2  = np.real(np.fft.ifftn(rho_k * ramp))
   ```
   This is **exact** for band-limited fields on a periodic grid (no
   interpolation in any conventional sense — it's the unitary
   translation operator on the discrete-FFT representation of a
   continuous field, accurate to numerical round-off). User confirmed
   this is permitted (the "no interpolation" instruction was about
   temporal interpolation between VTI frames, not spatial translation
   on each existing frame).
4. Subtract `ρ_RF2(r, t=0)` from each shifted frame to form a comoving
   delta `δρ_RF2(r, t)`. (This is the analogue of `density_rt_delta`
   but in the comoving frame.)
5. Repeat the Stage A pixel-by-pixel FFT on `δρ_RF2`. The
   ion-translation component is by construction zero in the comoving
   frame (the ions are stationary), so any spectral peak is genuine
   electronic dynamics.

**Sanity check before trusting Stage B:** integrate
`δρ_RF2(r, t=0..N)` over r at each t. The integrated value should be
zero to numerical precision (no charge added by the boost — the boost
is a unitary translation, not a re-normalisation). If it's not, there's
a bug.

### Stage C — cross-validation against the existing diagnostics

Whichever Stage produces a plasmon-candidate peak, validate against:

| Cross-check | Expected at the plasmon (6.5 eV) | Expected at e-h |
|---|---|---|
| `dipole_x` FFT (already built) | peak at the same ω (q→0 longitudinal mode) | absent at 6.5 eV; present at 2.8 eV in the high-v run |
| `gamma_transitions` histogram | **gap** at 6.5 eV (paper Figure 5; max Γ-Γ transition is 5.74 eV) | cluster around the e-h ω |
| `state_energy_spectra` anti-phase pair diagnostic | **no** (n, n′) pair with opp_metric > 0.7 near 6.5 eV (collective mode) | cluster of opp_metric > 0.7 pairs near the e-h ω |

If all three cross-checks are consistent with the plasmon
interpretation, write up the result. If any cross-check disagrees,
investigate.

## 4. Implementation

### 4.1 New postprocess module

`inq-stack/python/inqview/postprocess/density_spectra.py` — phase
named `density_spectra`, runnable via the existing pipeline driver.

- Reads VTI series via the existing `inqview.vti` / `_load_vti_array`
  helpers from `density.py`. Uses the `dx_bohr` from VTI metadata —
  no hardcoded grid spacing.
- Has two modes: `frame="lab"` (Stage A) and `frame="comoving"`
  (Stage B). The comoving mode requires `v_au` (vector of length 3),
  parsed from `run_summary.txt::kick_velocity_au` × kick direction.
- Outputs go in
  `results/analysis/density_spectra/<frame>/{2d_maps,3d_isosurface,
  cell_averaged_S_omega.png, S_omega.csv}`.
- Reuses `_hann_fft` and the `plateau_detrend` recipe from
  `observables.py` for consistency.

### 4.2 Known-case test (per `.claude/rules/development-feedback-loop.md`)

Before applying to the v=0.0626 run, validate the module on a
synthetic input:
- 54³ grid, 156 frames, dt = 4 a.u.
- Synthetic δρ(r, t) = A(r) · cos(ω₀ t + φ(r)) with ω₀ = 6.5 eV and a
  spatially-uniform A. Add a low-amplitude ω = 1 eV "ion-drift"
  background.
- Run `density_spectra.run(..., frame="lab")` and verify the
  cell-averaged S(ω) peaks at 6.5 eV and the 6–7 eV map is uniform.
- Run `density_spectra.run(..., frame="comoving", v_au=(0.0626, 0, 0))`
  on a *Galilean-boosted* version of the same synthetic input. Verify
  the result matches the unboosted lab-frame result to numerical
  precision (~1e-12 max diff).

### 4.3 Dry-run while v=0.0626 is propagating (no GPU needed)

Per Q-B, we do Stage A on the existing v=0.0123 data **before** the
v=0.0626 run finishes. v=0.0123 is in the same low-v family
(0.27–1.37 Å/fs); the paper's Figure 4(a) shows the 6.5 eV peak is
*near-constant* across this family. So the v=0.0123 lab-frame Stage A
should already show a 6.5 eV peak in the cell-averaged S(ω). This is a
free, end-to-end pipeline test on real data. If v=0.0123 doesn't show
the peak, our analysis is wrong; if it does, we have high confidence
the same recipe will work on v=0.0626.

We also dry-run Stage B on v=0.0123 to check the comoving-frame
implementation — though for v=0.0123 the ion displacement per frame is
0.026 Å = 0.13 pixels, so the Galilean shift is small and Stage A and
Stage B should give nearly identical answers.

## 5. Outputs (under `results/analysis/density_spectra/`)

```
density_spectra/
├── lab/
│   ├── cell_averaged_S_omega.png      # 1D spectrum, mark expected 6.5 eV
│   ├── cell_averaged_S_omega.csv      # raw spectrum data
│   ├── 2d_maps/
│   │   ├── band_0_0_2eV_xy.png        ; one PNG per (band, plane)
│   │   ├── band_0_0_2eV_xz.png
│   │   └── ... (7 bands × 3 planes = 21 PNGs)
│   ├── 3d_isosurface/
│   │   └── plasmon_5_8eV.vti          # ParaView-friendly
│   └── diagnostics/
│       ├── stage_a1_or_a2_decision.txt   # summary of low-ω vs 5-8 eV power
│       └── ion_translation_dominance_ratio.txt
├── comoving/                          # only if Stage A2 escalated
│   └── (same layout)
└── cross_validation/
    ├── dipole_x_vs_density_S_omega.png # overlay of two FFTs
    ├── gamma_transitions_with_peak.png # histogram with peak marked
    └── state_energy_pairs_summary.txt  # any opp-phase pairs near peak
```

## 6. Risks and mitigations

| Risk | Mitigation |
|---|---|
| 156 frames × 54³ × 4 bytes ≈ 2 GB of float32 in memory | streaming-FFT version available if RAM is tight; checked with user, confirmed OK |
| Wrap-around aliasing in the Fourier-shift (Stage B) when ion drift = 2 cell lengths | None needed — the shift is exact under periodic BCs (which this cell is); the periodic image is what we want for the comoving frame |
| Plasmon peak shifted vs paper's 6.5 eV due to GS offset (see `li_gs_xyz_vs_fractional_offset_analysis.md`) | Treat as confirmed if peak is in [6.0, 7.0] eV; if outside that, investigate |
| Hann + plateau-detrend insufficient at v=0.0626 if signal is dominated by t<1 fs transient | Existing `t_skip_fs` knob in `_extended_spectra` handles this; module exposes the same knob |
| User's manual annotation suggested 7–8 eV (paper says 6.5 eV) — my method might land on the user's value if the GS offset shifts things | Document both peaks, mark which is closer to paper; the 1.5 eV gap is too large to be GS noise |

## 7. Effort and wall-time

- Stage A scaffold + known-case test: ~1 h dev (during plasmon run)
- Stage A dry-run on v=0.0123: ~5 min wall on this CPU
- Stage A on v=0.0626 (when run finishes): ~5 min wall
- Stage B (if needed): +15 min wall (Fourier-shift adds 156 FFTs)
- Cross-validation cross-checks: instant (existing CSVs)
- Journal entry update (canonical run_summary.txt + figures + writeup):
  ~30 min

Total compute: ~30 min. Total wall (including dev): ~3 h spread across
the remaining plasmon-run window (~16 h to go), so by the time the
plasmon run finishes the analysis is ready to apply with one command.

## 8. Decision gates

Before writing the module:
- **G1.** Plan signed off by user.

Before applying to the v=0.0626 results:
- **G2.** Module passes the synthetic known-case test (Stage A and
  Stage B identity check on Galilean-boosted synthetic).
- **G3.** Module dry-runs cleanly on the existing v=0.0123 data and
  reproduces the expected 6.5 eV peak (paper Figure 4a).

Before flipping the journal entry status to `complete`:
- **G4.** Stage A or B on v=0.0626 produces a peak; all three Stage C
  cross-checks are consistent with the plasmon attribution (or
  flagged as inconsistent for further investigation).

## 9. Out of scope

- q-resolved plasmon dispersion (would need separate per-q-mode FFT;
  bulk plasmon at q=0 is sufficient for the goal).
- Correlation with density-of-states near E_F (would need DOS
  postprocess wired through `eigenvalues_gs`; not needed for the
  identification question).
- Comparison against QBall density (QBall doesn't write density VTI
  series in this run set — only energy + forces are stored).

## 10. Sign-off

User to review. After review, suggested edits will be applied to this
plan, then executed under the gates above. The journal entry
(`docs/journals/quantumkickextension/2026-05-06_run_propagate_v0p0626_xyz.md`)
will be updated in place with the final results — not a new entry.
