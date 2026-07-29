// Stage 1 of 2-GPU MPI Run-8 attempt (2026-05-20).
//
// Loads the dx=0.30 GS, injects the WP single-rank, saves the post-injection
// state to a new checkpoint that Stage 2 (run_wp_n162_L50_E700_mpi_propagate)
// will load under mpirun -np 2. This is the workaround for the WP injector's
// hard single-rank constraint (wavepacket.hpp:129).

#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include "../shared/configs/knudsen_sweep_L50_cubic.hpp"

#include <filesystem>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::KnudsenSweep_L50_cubic_WP_E700;

int main() {
    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p30";
    const std::string OUT_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p30_E700_wp_injected";

    std::cout << "Stage 1: inject WP at E=" << Cfg::WP_EKIN_EV
              << " eV into dx=0.30 GS, save to " << OUT_DIR << "\n";

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);

    auto report = inqkit::WavePacket{}
                      .center(Cfg::WP_CX_BOHR, Cfg::WP_CY_BOHR, Cfg::WP_CZ_BOHR)
                      .sigma(Cfg::WP_SIGMA_BOHR)
                      .k0(Cfg::WP_KX, Cfg::WP_KY, Cfg::WP_KZ)
                      .orthogonalise_against_occupied(electrons)
                      .inject_into_last_extra_state(electrons, 1.0);
    std::cout << "  WP injected, state=" << report.state_index
              << ", norm=" << report.norm_after << "\n";

    std::filesystem::create_directories(OUT_DIR);
    electrons.save(OUT_DIR);
    std::cout << "  Saved " << OUT_DIR << "\n";
    return 0;
}
