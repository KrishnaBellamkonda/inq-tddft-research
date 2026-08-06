#pragma once
// ============================================================================
// twodef_gs.hpp — ground-state config for the TWO-DEFINITIONS stopping campaign
// (bilayer + monolayer graphene). Plan: docs/plans/real-material-stopping-comparison.md
//
// Cell: commensurate 3x2 rectangular graphene supercell (a = 2.46 Ang),
//   Lx x Ly = 13.9462 x 16.1037 Bohr; nx=3 folds the Dirac K point -> Gamma
//   (inherited from graphene_gs.hpp — verified design constraint).
//   ORTHORHOMBIC, periodicity(2): xy periodic, z FINITE (slab, corpus per2
//   convention — resolves the campaign's periodicity route (a)).
// Variants: mono (24 C, z=0) and AB bilayer (48 C, layers at z = +/-1.675 Ang,
//   d = 3.35 Ang Bernal; layer 2 shifted by one C-C bond vector).
// XC: LDA/ALDA. Pseudos: INQ default ONCV norm-conserving C (4 valence e-).
// Cutoff: 50 Ha default ("paper" value from graphene_gs.hpp); Phase 0 battery
//   sweeps {35,45,50} on the bilayer to CHECK it, not assume it.
// Smearing: 0.1 eV (semimetal; folded-Gamma supercell needs it to converge).
// ============================================================================
namespace graphene_twodef {

// --- Box (Bohr), commensurate with the 3x2 supercell ------------------------
inline constexpr double LX_BOHR = 13.9462;   // 3 * 2.46 Ang
inline constexpr double LY_BOHR = 16.1037;   // 2 * 2.46*sqrt(3) Ang
inline constexpr double LZ_BOHR = 80.0;      // default z extent (CAP scan may revise)

// --- Electronic structure ---------------------------------------------------
inline constexpr double CUTOFF_HA      = 50.0;
inline constexpr double TEMPERATURE_EV = 0.10;
inline constexpr int    EXTRA_STATES_MONO = 24;   // above 48 occupied
inline constexpr int    EXTRA_STATES_BI   = 48;   // above 96 occupied

// --- SCF (as graphene_gs.hpp) ----------------------------------------------
inline constexpr double SCF_TOL_HA    = 1e-6;
inline constexpr int    SCF_MAX_STEPS = 300;
inline constexpr int    SCF_MIX_NDIM  = 8;
inline constexpr double SCF_MIX_ALPHA = 0.3;

// --- Bookkeeping ------------------------------------------------------------
inline constexpr int N_C_MONO = 24;   inline constexpr int N_ELEC_MONO = 96;
inline constexpr int N_C_BI   = 48;   inline constexpr int N_ELEC_BI   = 192;

inline constexpr char REPO[] = "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research";
inline constexpr char GEOM_MONO[] =
  "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/graphene/shared/geometry/graphene_3x2.xyz";
inline constexpr char GEOM_BI[] =
  "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/graphene/shared/geometry/graphene_3x2_bilayer.xyz";
inline constexpr char GS_DIR_BASE[] =
  "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/graphene/shared_gs";

} // namespace graphene_twodef
