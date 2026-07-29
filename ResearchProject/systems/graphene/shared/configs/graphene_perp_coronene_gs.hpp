#pragma once
// ============================================================================
// graphene_perp_coronene_gs.hpp — GS config for the PERPENDICULAR arm of the
// impact-parameter comparison.
//
// SAME finite coronene C24H12 flake as the grazing arm, but in its NATIVE x-y
// plane (flake normal = z = the beam axis) so a +z projectile hits it
// head-on / PERPENDICULAR. Identical box / grid / cutoff / electron count /
// CAP / projectile to graphene_grazing_gs.hpp so the perpendicular and grazing
// runs are DIRECTLY COMPARABLE — only the flake orientation differs (the same
// molecule rotated 90 deg). User correction 2026-06-21: same target + same
// box/grid/sim-time in both arms.
//
// Same namespace `graphene_cfg` (distinct file) so the GS run.cpp and the
// cl/wp run.cpp are reused verbatim with only the #include swapped.
// ============================================================================
namespace graphene_cfg {

// --- Box (Bohr) — IDENTICAL to the grazing arm (same grid @50 Ha) ------------
inline constexpr double LX_BOHR = 20.0;
inline constexpr double LY_BOHR = 22.0;
inline constexpr double LZ_BOHR = 60.0;   // +z traversal + two-sided CAP; flake at z=0

// --- Electronic structure (identical) ---------------------------------------
inline constexpr double CUTOFF_HA      = 50.0;
inline constexpr int    EXTRA_STATES   = 24;
inline constexpr double TEMPERATURE_EV = 0.01;   // gapped molecule

// --- SCF --------------------------------------------------------------------
inline constexpr double SCF_TOL_HA   = 1e-6;
inline constexpr int    SCF_MAX_STEPS = 300;
inline constexpr int    SCF_MIX_NDIM  = 8;
inline constexpr double SCF_MIX_ALPHA = 0.3;

// --- Bookkeeping ------------------------------------------------------------
inline constexpr int N_C_ATOMS   = 36;          // total atoms (24 C + 12 H)
inline constexpr int N_ELECTRONS = 108;

inline constexpr char GEOMETRY_XYZ[] =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared/geometry/coronene_flake_perp.xyz";
inline constexpr char GS_CHECKPOINT_DIR[] =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared_gs/gs_perp_coronene_50ha";

} // namespace graphene_cfg
