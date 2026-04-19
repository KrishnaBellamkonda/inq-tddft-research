# Checkpoint: free-gaussian-wp-propagation

**Conversation name:** free-gaussian-wp-propagation  
**Date:** 2026-04-13  
**Status:** Complete

Full handover: see `docs/handovers/free_wavepacket_propagation.md`

## Summary

Implemented and validated a free-particle 3D Gaussian wavepacket propagation in C++.

### Files created

| File | Purpose |
|---|---|
| `ResearchProject/jellium/03_free_gaussian_wp_propagation/run.cpp` | C++ simulation (FFTW3, 161³ grid) |
| `ResearchProject/jellium/03_free_gaussian_wp_propagation/CMakeLists.txt` | Standalone CMake build |
| `ResearchProject/jellium/03_free_gaussian_wp_propagation/plot_propagation.py` | Visualisation (PyVista + matplotlib) |
| `.claude/skills/build-run.md` | New skill: `source ~/.bashrc` + inq-run / manual CMake patterns |

### Outputs

- `results/broadening_comparison.png` — σ(t) measured vs analytical
- `results/heatmap_animation.mp4` — 2D density animation (15 fps)
- `results/isosurface_animation.mp4` — 3D PyVista animation (10 fps)
- `results/width_vs_time.csv` — 121 frames, t, σ_x, σ_y, σ_z, σ_analytical

### Validation

- σ_x(t=0) = 1.0000 ✓, σ_x(t=2) = 1.4142 = √2 ✓ (exact analytical match)
- Spherical symmetry: σ_x = σ_y = σ_z throughout ✓

### Rules updated

- `.claude/rules/file-placement.md` — always PNG, never PDF/SVG
- `.claude/rules/testing.md` — use GPU whenever available
