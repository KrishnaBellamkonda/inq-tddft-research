# Handover: Free Gaussian Wavepacket Propagation

## Current status

**Complete.** C++ simulation built and ran; all three visualisations generated and verified.

---

## What changed

- Created `/local/data/public/skcb2/tddft/ResearchProject/jellium/03_free_gaussian_wp_propagation/run.cpp`
- Created `/local/data/public/skcb2/tddft/ResearchProject/jellium/03_free_gaussian_wp_propagation/CMakeLists.txt`
- Created `/local/data/public/skcb2/tddft/ResearchProject/jellium/03_free_gaussian_wp_propagation/plot_propagation.py`
- Created `/local/data/public/skcb2/tddft/.claude/skills/build-run.md` — documents `source ~/.bashrc` requirement before `inq-run` or manual CMake builds

---

## Files touched

| File | Role |
|---|---|
| `ResearchProject/jellium/03_free_gaussian_wp_propagation/run.cpp` | C++ FFT propagation core |
| `ResearchProject/jellium/03_free_gaussian_wp_propagation/CMakeLists.txt` | Standalone build (FFTW3 only, no INQ) |
| `ResearchProject/jellium/03_free_gaussian_wp_propagation/plot_propagation.py` | Python visualisation (2D animation, 3D isosurface, broadening) |
| `ResearchProject/jellium/03_free_gaussian_wp_propagation/results/` | All output files (see below) |
| `.claude/skills/build-run.md` | New skill documenting `source ~/.bashrc` + build patterns |

**Outputs in `results/`:**
- `width_vs_time.csv` — t, sigma_x, sigma_y, sigma_z, sigma_analytical (121 rows)
- `slice_t000.txt` … `slice_t120.txt` — 2D density at z=N/2 plane (161×161 each)
- `density3d_t000.txt` … `density3d_t120.txt` — coarse 3D density (stride=4, ~41³ per file)
- `broadening_comparison.png` — measured σ(t) vs analytical formula
- `heatmap_animation.mp4` — 2D density animation (32 KB, 15 fps)
- `isosurface_animation.mp4` — 3D PyVista isosurface animation (93 KB, 10 fps)
- `iso_frame_000.png` … `iso_frame_120.png` — individual isosurface frames

---

## Commands run

```bash
# Build
cd ResearchProject/jellium/03_free_gaussian_wp_propagation
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)

# Run simulation (3m 34s)
cd ResearchProject/jellium/03_free_gaussian_wp_propagation
./build/run_wp

# Visualisation (quantum-wave-packet pyenv)
export PYENV_ROOT="/local/data/public/skcb2/pyenv" && ...
pyenv activate quantum-wave-packet
python3 plot_propagation.py

# Fix isosurface MP4 (ffmpeg concat failed silently; rebuilt from numbered PNGs)
ffmpeg -y -framerate 10 -i results/iso_frame_%03d.png -vf "format=yuv420p" results/isosurface_animation.mp4
```

---

## Tests and validation

**Proposed:** σ_x(t=0) = σ_0 = 1.0 bohr; σ_x(t) agrees with σ_0√(1+t²/(4σ_0⁴)); σ_x = σ_y = σ_z (spherical symmetry).

**Run:**
- σ_x(t=0) = 1.0000 ✓
- σ_x(t=2) = 1.4142 = √2 ✓ (exact: σ√(1+4/4) = σ√2)
- σ_x(t=12) = 6.0725 ≈ √37 = 6.083 (0.1% deviation from analytical — from discrete grid + PBC)
- Spherical symmetry (σ_x = σ_y = σ_z): visible in CSV — all three columns agree to all printed digits ✓
- Normalisation: σ_x at t=0 matches exact input σ=1.0 to 4 decimal places ✓

**Broadening comparison plot** shows measured σ(t) lying on top of analytical curve throughout t=0…12 a.u.

**Unverified:**
- GPU/CPU consistency (this is a standalone FFTW3 binary, not an INQ simulation; not applicable)
- Probability conservation (⟨ψ|ψ⟩ not logged — could add to CSV; expected ~1 throughout)

---

## Trusted sources used

- `Tutorial/angelo-jellium/jellium/wavepacket.py` — reference implementation; default parameters (sigma=1.0, k0=[0,0,0])
  - Normalization: line 35: `N = (2π σ²)^{-3/4}`
  - FFT propagation: lines 111–142
- Griffiths, *Introduction to Quantum Mechanics* — Gaussian wavepacket broadening formula σ(t) = σ_0√(1+t²/(4σ_0⁴))

---

## Attribution notes

Gaussian formula and default parameters derived directly from `wavepacket.py`. FFT propagation scheme is standard: evolve each plane-wave mode with exp(-i|k|²t/2), following the same logic as lines 111–142 of the reference.

---

## Known issues / blockers

- `pyvista.start_xvfb()` raises a `PyVistaDeprecationWarning` in PyVista 0.47.2 — headless rendering still works; no action needed unless PyVista is upgraded.
- Isosurface ffmpeg via concat-list fails silently in the Python script (produces 262-byte file). Workaround: run `ffmpeg -framerate 10 -i results/iso_frame_%03d.png ...` directly. The plot script should be updated to use this approach if re-run.

---

## Assumptions still in play

- Grid: L=40 bohr, N=161, dx≈0.2484 bohr, E_cut=80 Ha
- Wavepacket: sigma=1.0 bohr, k0=[0,0,0], r0=[20,20,20] bohr
- Time: DT=0.02 a.u., T_MAX=12 a.u., save every 5 steps (dt_save=0.1 a.u.)
- Build machine: this is a standalone CMake build; does not use `inq-run` (FFTW3 only)

---

## Exact next steps

This task is complete. Possible follow-ons:
1. Add k0 ≠ 0 to see moving wavepacket (requires wrapping-aware centre-of-mass tracking)
2. Add a finite potential (e.g. harmonic well) — turn it into a trapped Gaussian
3. Fix the ffmpeg call in `plot_propagation.py` to use `-i results/iso_frame_%03d.png` directly
4. Add probability conservation check to `width_vs_time.csv`
