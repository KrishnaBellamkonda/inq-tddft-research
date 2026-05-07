# Plan — Jellium WP-RT hypothesis comparisons

Mirror of the coronene `hypotheses/` framework for `ResearchProject/jellium/jellium-wp-rt/`.

---

## 1. Existing run inventory (skip `run_08_coronene_leed`)

| Run | E (Ha) | E (eV) | σ (Bohr) | k₀ (Bohr⁻¹) | direction | N_e | r_s (Bohr) |
|---|---|---|---|---|---|---|---|
| `run_01_base` | 7.35 | 200 | 1.00 | 3.834 | +z | 38 (closed) | 7.38 |
| `run_02_low_energy` | 1.84 | 50 | 1.00 | 1.917 | +z | 38 | 7.38 |
| `run_03_high_energy` | 14.70 | 400 | 1.00 | 5.422 | +z | 38 | 7.38 |
| `run_04_tilted_45` | 7.35 | 200 | 1.00 | 3.834 | 45° xz | 38 | 7.38 |
| `run_05_wide_sigma` | 7.35 | 200 | 3.78 | 3.834 | +z | 38 | 7.38 |
| `run_06_narrow_sigma` | 7.35 | 200 | 0.50 | 3.834 | +z | 38 | 7.38 |
| `run_07_open_shell` | 7.35 | 200 | 1.00 | 3.834 | +z | 40 (open) | 7.26 |

Cell L = 40 Bohr cubic, periodic. dt = 0.02 a.u. throughout.

---

## 2. New runs to add (2 GPU runs, ≈ 30–60 min each)

| Run | E (eV) | σ (Bohr) | Reason |
|---|---|---|---|
| `run_08_fast_narrow` | 400 | 0.50 | Fast projectile + narrow σ — extreme Gaussian-broadening test (high k₀, low Δx ⇒ rapid σ_z(t) growth; expected analytic spread ≈ 1.4 σ₀ at T_sim) |
| `run_09_slow_wide` | 50 | 3.78 | Slow + wide WP — analogue of coronene `run_E30_s3`; tests electron capture into jellium quasi-bound states |

(I am **not** renaming `run_08_coronene_leed`; the new runs slot in as `run_08_fast_narrow` after deleting the obsolete dir and `run_09_slow_wide`. If you'd rather keep the coronene_leed dir, I can use `run_10` / `run_11` instead — say which.)

---

## 3. Hypotheses folder structure

```
ResearchProject/jellium/jellium-wp-rt/hypotheses/
  00_base/                            run_01 only
  01_wp_energy_spread/                run_02 + run_01 + run_03 (50/200/400 eV)
  02_wp_sigma_spread/                 run_06 + run_01 + run_05 (σ=0.5/1.0/3.78)
  03_open_vs_closed_shell/            run_01 + run_07
  04_tilted_propagation/              run_01 + run_04
  05_fast_narrow_broadening/          run_03 + run_08_fast_narrow (NEW)
                                      Gaussian-broadening analytic overlay
  06_electron_capture/                run_02 + run_09_slow_wide (NEW)
                                      capture diagnostic
```

Per folder: LEED-screen overlays, current/dipole/energy spectra overlays (using existing `jellium_spectra.py` outputs), residual-norm bar at t_final, **per-config metadata table** (E_Ha / σ_Bohr / r_s / k₀ / N_e / shell / direction).

---

## 4. New analyses

### 4a. Free-particle Gaussian broadening overlay

For runs `run_02`, `run_03`, `run_06`, `run_08_fast_narrow`:

Analytic free-particle σ(t):
```
σ(t) = σ₀ √(1 + (t / (2 m σ₀²))²)
```
where m = 1 (electron in atomic units), so τ_disp = 2σ₀² is the dispersion time.

Measured σ(t) from `density_rt_wp/`: at each frame, fit a 1-D Gaussian along z to the
projected density `n_WP(z, t) = ∫∫ n_WP(x, y, z, t) dx dy` over the WP slab, extract
σ_z(t).

Output: `<hypothesis>/gaussian_broadening_overlay.png` — measured σ_z(t) (markers)
vs analytic σ_free(t) (line). Departure magnitude indicates the strength of the
jellium response (additional broadening from many-body effects).

### 4b. Electron capture diagnostic for `run_09_slow_wide`

Three quantities, each computed at every snapshot:

1. **Trapped-density fraction**:
   ```
   f_trap(t) = ∫_{slab} n_WP(r, t) d³r / ∫_{cell} n_WP(r, 0) d³r
   ```
   where `slab` = a 4 σ₀-thick window around the jellium centre. `f_trap` rising
   above ≈ 0.05 over time = capture signature.

2. **Mean WP momentum drift**:
   ```
   ⟨k_z⟩(t) = ⟨J_z⟩(t) / ⟨n_WP⟩(t) (already in observables.csv)
   ```
   Captured electrons lose their initial k₀ ⇒ ⟨k_z⟩ decays.

3. **WP overlap with occupied GS orbitals** (already saved in `results/overlap/`).
   Capture into a quasi-bound state shows up as a non-vanishing overlap with one
   or more occupied orbitals at t_final.

Output: `06_electron_capture/capture_diagnostics.png` (3-panel) plus a
`capture_summary.txt` numerical readout.

---

## 5. Metadata table generator

`inq-stack/python/inqview/postprocess/jellium_compare.py` (new): mirrors
`inqview/postprocess/compare.py` but for the flat jellium results layout. Reads each
run's `run.cpp` header comment block (E_eV, σ_bohr, r_s, k₀, direction, N_e),
emits a markdown + PNG table.

---

## 6. Cost & ordering

| Stage | Cost | Notes |
|---|---|---|
| Inventory + metadata table | 5 min CPU | reads existing run.cpp comments |
| Hypotheses 00–04 from existing 7 runs | 10–15 min CPU | reuse `jellium_spectra.py` outputs |
| 4a Gaussian broadening overlay (existing high/low-E + narrow-σ runs) | 5 min CPU | density_rt_wp already on disk |
| **NEW** `run_08_fast_narrow` build + run | ≈ 20 min build (cached) + 30 min RT on GPU | smaller σ⇒small dt cost; dispatcher slot |
| **NEW** `run_09_slow_wide` build + run | ≈ 20 min build + 60 min RT on GPU | larger σ⇒longer T_sim |
| 4b electron-capture diagnostic on run_09 | 5 min CPU | reads run_09 results |
| Hypotheses 05 + 06 with new runs | 5 min CPU |

Total wall time ≈ **2–2.5 h** if both GPUs available; ≈ 3 h sequential.

---

## 7. Open questions for you

1. **`run_08_coronene_leed`** — keep it (use `run_09`/`run_10`/`run_11` for new ones)
   or delete it? You said "ignore run_08", which I read as "skip in hypotheses"
   not "delete". Confirm.
2. **`hypotheses/04_tilted_propagation`** — worth a folder, or just a side-note in
   `00_base`? It's only 2 runs.
3. **Per-config density viz** in each hypothesis folder — yes/no? Adds ≈ 100 MB of
   GIFs per hypothesis.
4. **Are the 2 new GPU runs authorised?** Cost above. If not, I'll deliver
   hypotheses 00–05 (Gaussian broadening overlay only on existing data) and
   defer 06 + the new runs.

I'll start on the cheap parts (inventory, metadata, hypotheses 00–04, Gaussian
overlay on existing data) immediately and pause before launching the two new
GPU runs.
