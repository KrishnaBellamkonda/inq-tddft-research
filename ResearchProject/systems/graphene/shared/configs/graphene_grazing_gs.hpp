#pragma once
// ============================================================================
// graphene_grazing_gs.hpp — GS config for the GRAZING / impact-parameter study.
//
// A FINITE 2-D graphene flake (coronene C24H12, H-passivated) reoriented into
// the y-z plane (flake normal = x = the impact-parameter axis), so a +z
// projectile travels PARALLEL to the sheet and grazes it at perpendicular
// distance b = x-offset. This is the user's intended grazing geometry (a finite
// sheet, NOT periodic bulk; 2026-06-21). Box is large with vacuum on all sides
// (the flake is isolated); z keeps the 60-Bohr traversal + two-sided z-CAP of
// the perpendicular runs.
//
// Uses namespace `graphene_cfg` (same symbol names as graphene_gs.hpp) so the
// GS run.cpp is reused verbatim with only the #include swapped.
//
// Coronene C24H12: 24*4 + 12*1 = 108 valence e (ONCV C=4, H=1), closed-shell
// molecule with a HOMO-LUMO gap -> small smearing only.
// ============================================================================
namespace graphene_cfg {

// --- Box (Bohr) — finite flake + vacuum; z = traversal + CAP -----------------
inline constexpr double LX_BOHR = 20.0;   // impact-parameter axis (flake normal); vacuum + b range
inline constexpr double LY_BOHR = 22.0;   // flake in-plane (~18 Bohr) + margin
inline constexpr double LZ_BOHR = 60.0;   // +z traversal + two-sided CAP (|z|<20 free), flake |z|<8.7

// --- Electronic structure ---------------------------------------------------
inline constexpr double CUTOFF_HA      = 50.0;   // match the perpendicular-campaign / CAP numerics
inline constexpr int    EXTRA_STATES   = 24;     // above 54 occupied (closed shell) + WP/CAP headroom
inline constexpr double TEMPERATURE_EV = 0.01;   // tiny smearing (gapped molecule; numerical safety)

// --- SCF --------------------------------------------------------------------
inline constexpr double SCF_TOL_HA   = 1e-6;
inline constexpr int    SCF_MAX_STEPS = 300;
inline constexpr int    SCF_MIX_NDIM  = 8;
inline constexpr double SCF_MIX_ALPHA = 0.3;

// --- Bookkeeping ------------------------------------------------------------
inline constexpr int N_C_ATOMS   = 36;          // TOTAL atoms (24 C + 12 H) — for the bounds/expect check
inline constexpr int N_ELECTRONS = 108;         // 24*4 (C) + 12*1 (H)

inline constexpr char GEOMETRY_XYZ[] =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared/geometry/coronene_flake_grazing.xyz";
inline constexpr char GS_CHECKPOINT_DIR[] =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared_gs/gs_grazing_coronene_50ha";

} // namespace graphene_cfg
