#pragma once
// ============================================================================
// graphene_gs.hpp — ground-state config for the graphene+CAP feasibility replica
//
// Feasibility replica of Yao & Schleife (32 C, not 112). Decisions locked in
// docs/plans/graphene-cap.md (grill 2026-06-18). NOT the paper's converged
// numbers.
//
// Cell: commensurate 4x2 rectangular graphene supercell (a = 2.46 Ang),
//   in-plane Lx x Ly = 18.5949 x 16.1037 Bohr (= 9.84 x 8.52 Ang), 32 C atoms,
//   sheet at z = 0, box z = 60 Bohr. ORTHORHOMBIC + PERIODIC (Gamma-point).
// XC: LDA/ALDA. Cutoff: 50 Ha (paper). Norm-conserving (INQ default ONCV-C).
// Smearing: electronic temperature (graphene is a semimetal; folded Gamma
//   supercell needs smearing to converge).
// ============================================================================
namespace graphene_cfg {

// --- Box (Bohr), commensurate with the 4x2 graphene supercell ---------------
inline constexpr double LX_BOHR = 13.9462;   // 3 * 2.46 Ang (nx=3 folds K->Gamma)
inline constexpr double LY_BOHR = 16.1037;   // 2 * 2.46*sqrt(3) Ang
inline constexpr double LZ_BOHR = 60.0;      // z-vacuum + two-sided CAP

// --- Electronic structure ---------------------------------------------------
inline constexpr double CUTOFF_HA      = 50.0;   // paper plane-wave cutoff
inline constexpr int    EXTRA_STATES   = 24;     // above 48 occupied (semimetal)
inline constexpr double TEMPERATURE_EV = 0.10;   // electronic smearing (metallic GS)

// --- SCF --------------------------------------------------------------------
inline constexpr double SCF_TOL_HA   = 1e-6;
inline constexpr int    SCF_MAX_STEPS = 300;
inline constexpr int    SCF_MIX_NDIM  = 8;
inline constexpr double SCF_MIX_ALPHA = 0.3;

// --- Bookkeeping ------------------------------------------------------------
inline constexpr int N_C_ATOMS   = 24;
inline constexpr int N_ELECTRONS = 96;          // 4 valence x 24 C (ONCV-C)

inline constexpr char GEOMETRY_XYZ[] =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared/geometry/graphene_3x2.xyz";
inline constexpr char GS_CHECKPOINT_DIR[] =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared_gs/gs_3x2_50ha";

} // namespace graphene_cfg
