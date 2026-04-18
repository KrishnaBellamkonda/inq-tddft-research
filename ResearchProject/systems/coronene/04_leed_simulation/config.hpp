#pragma once
// ============================================================================
// 04_leed_simulation/config.hpp
//
// Core configuration for the coronene TDDFT electron WP scattering simulation.
// Reuses the shared parameter block from 02_ground_state_analysis/config.hpp.
//
// After running 03_ecut_convergence, update ECUT_HA below to the converged
// value (expected ~54 Ha for paper-quality 0.16 Å grid).
// ============================================================================

#include "../02_ground_state_analysis/config.hpp"

// Override E_cut once convergence test is done.
// Update this value after running 03_ecut_convergence.
// Paper grid: 0.16 Å = 0.302 bohr → E_cut ≈ 54 Ha.
// For a quick proof-of-concept run first, use 30 Ha; for paper comparison use 54 Ha.
namespace cfg {
    // E_cut set from 03_ecut_convergence sweep (20–60 Ha, LDA, pseudodojo_pbe).
    //
    // Convergence result: energy is non-monotonic above 40 Ha.
    //   E_cut=40 Ha gives the lowest total energy (-150.837 Ha, h=0.186 Å).
    //   Energy rises from 40→60 Ha (+90 meV), which is atypical and likely
    //   reflects aliasing between the pseudopotential projectors and the real-space
    //   grid at these cutoffs with pseudodojo_pbe.  The paper (Tsubonoya et al.
    //   PRB 90, 035416, 2014) used 54 Ha with a different code and pseudopotential.
    //
    // Choice: 40 Ha — energy minimum in our sweep; h=0.186 Å gives ~3 grid points
    // across the WP width (d=1.001 bohr), comparable to the paper's 3.3 points.
    // If paper-matching spatial resolution is required, switch to 55 Ha.
    inline constexpr double ECUT_HA_LEED = 40.0;  // Ha  (energy minimum in 03 sweep)

    // Density snapshot z-slices (for Fig. 1 replica):
    // Store 2D slices at z=0 (flake plane) and z=Z_OBS_BOHR (observation plane)
    inline constexpr double Z_FLAKE_BOHR = 0.0;    // coronene plane

    // Maximum number of snapshots to save
    inline constexpr int MAX_SNAPSHOTS = 20;
}
