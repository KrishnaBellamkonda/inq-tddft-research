// ============================================================================
// graphene/scripts/lovelace_test/gs/run.cpp
//
// Ground-state builder for the lovelace_test stopping campaign.
// Bilayer AB graphene (default), orthorhombic 3x2 supercell.
//
// Key choices (lovelace_test):
//   spacing = 0.5 Bohr  (GR_DX_BOHR; equiv. E_cut ~ 19.7 Ha — coarser than
//                         the production 50 Ha; intentional for fast test runs)
//   Lz     = 90 Bohr    (CAPs [−45,−35] and [+35,+45], WP launch at z=−19)
//   periodicity(2): xy periodic, z finite slab
//
// Env overrides:
//   GR_VARIANT  (bi|mono, default bi)
//   GR_GEOM     (path to .xyz; defaults to shared/geometry/graphene_3x2_bilayer.xyz)
//   GR_LZ_BOHR  (default 90)
//   GR_DX_BOHR  (default 0.5)
//   GR_OUT      (default results/bi_dx0p5_Lz90)
// ============================================================================
#include <inq/inq.hpp>
#include "../../../shared/configs/twodef_gs.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
namespace Cfg = graphene_twodef;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static std::string iso_now(){
    auto t = std::time(nullptr); auto tm = *std::localtime(&t);
    char b[64]; std::strftime(b, sizeof(b), "%Y-%m-%dT%H:%M:%S", &tm); return b;
}

int main() {
    auto t0 = std::chrono::steady_clock::now();

    const std::string VARIANT = env_s("GR_VARIANT", "bi");
    const bool BI = (VARIANT == "bi");
    const int  N_ATOMS = BI ? Cfg::N_C_BI   : Cfg::N_C_MONO;
    const int  EXTRA   = env_i("GR_EXTRA", BI ? Cfg::EXTRA_STATES_BI : Cfg::EXTRA_STATES_MONO);
    const double LZ    = env_d("GR_LZ_BOHR", 90.0);
    const double DX    = env_d("GR_DX_BOHR",  0.5);

    // Local geometry paths — CSD3 /rds paths in Cfg::GEOM_* are used as fallback only
    const std::string LOCAL_GEOM_BI =
        std::string(std::getenv("TDDFT_ROOT") ? std::getenv("TDDFT_ROOT") : "/local/data/public/skcb2/tddft")
        + "/ResearchProject/systems/graphene/shared/geometry/graphene_3x2_bilayer.xyz";
    const std::string LOCAL_GEOM_MONO =
        std::string(std::getenv("TDDFT_ROOT") ? std::getenv("TDDFT_ROOT") : "/local/data/public/skcb2/tddft")
        + "/ResearchProject/systems/graphene/shared/geometry/graphene_3x2.xyz";

    const std::string GEOM = env_s("GR_GEOM", BI ? LOCAL_GEOM_BI : LOCAL_GEOM_MONO);
    const std::string OUT  = env_s("GR_OUT", "results/" + VARIANT + "_dx0p5_Lz90");

    if (!std::filesystem::exists(GEOM)) {
        std::cerr << "FATAL: geometry file not found: " << GEOM
                  << "\n  Set GR_GEOM=/path/to/graphene_3x2_bilayer.xyz\n"; return 2;
    }

    std::cout << "\n=== lovelace_test/gs [" << VARIANT << "] ===\n"
              << "  cell    = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x " << LZ
              << " Bohr  periodicity(2)\n"
              << "  geom    = " << GEOM << "  (" << N_ATOMS << " C atoms)\n"
              << "  spacing = " << DX << " Bohr  extra_states = " << EXTRA << "\n"
              << "  temp    = " << Cfg::TEMPERATURE_EV << " eV (semimetal smearing)\n"
              << "  out     = " << OUT << "\n\n";

    auto cell = systems::cell::orthorhombic(
        Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, LZ * 1.0_b).periodicity(2);

    auto ions = systems::ions::parse(GEOM, cell);
    if (int(ions.size()) != N_ATOMS) {
        std::cerr << "FATAL: parsed " << ions.size() << " atoms, expected " << N_ATOMS << "\n";
        return 2;
    }

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(DX * 1.0_b)
            .extra_states(EXTRA)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(Cfg::SCF_TOL_HA * 1.0_Ha)
            .max_steps(Cfg::SCF_MAX_STEPS)
            .broyden_mixing()
            .mixing_ndim(Cfg::SCF_MIX_NDIM)
            .mixing(Cfg::SCF_MIX_ALPHA));

    const double E_tot = gs.energy.total();
    std::cout << "  GS energy = " << E_tot << " Ha\n";

    std::filesystem::create_directories(OUT);
    electrons.save(OUT);

    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = graphene/lovelace_test/gs/" << VARIANT << "\n"
          << "date_finished = " << iso_now() << "\n"
          << "variant = " << VARIANT << "\n"
          << "geometry_xyz = " << GEOM << "\n"
          << "n_atoms = " << N_ATOMS << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << LZ << "\n"
          << "periodicity = 2\n"
          << "dx_bohr = " << DX << "\n"
          << "temperature_ev = " << Cfg::TEMPERATURE_EV << "\n"
          << "extra_states = " << EXTRA << "\n"
          << "ground_state_energy_ha = " << E_tot << "\n"
          << "energy_per_atom_ha = " << E_tot / N_ATOMS << "\n"
          << "num_states = " << electrons.states().num_states() << "\n"
          << "checkpoint_dir = " << std::filesystem::absolute(OUT).string() << "\n"
          << "wall_time_s = " << wall << "\n"
          << "run_completed = true\n";
    }
    std::cout << "Done. Wall " << wall << " s. Checkpoint: " << OUT << "\n";
    return 0;
}
