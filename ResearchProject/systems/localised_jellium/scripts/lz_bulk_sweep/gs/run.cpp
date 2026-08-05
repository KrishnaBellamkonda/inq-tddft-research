// ============================================================================
// systems/localised_jellium/scripts/lz_bulk_sweep/gs/run.cpp
//
// Ground state for the slab->bulk L_slab sweep (`lz_bulk_sweep`).
// Plan: docs/plans/jellium-slab-extend-Lz.md
//
// Clone of scripts/sigma56_sv/gs/run.cpp with the geometry moved from the
// compile-time Cfg to the RUNTIME box preset (env LZB_CFG, lzb_boxes.hpp) so one
// binary serves all four boxes of the sweep: s0p5_L15 (75 Bohr, N=60),
// s0p5_L35 (95, 140), s5p0_L15 (95, 60), s5p0_L35 (115, 140). The SCF options,
// background construction, gates and outputs are unchanged.
//
// GATES (same philosophy as sigma56): electron count and r_s are HARD — they
// only move if the background wiring is broken, and every downstream S would be
// garbage. E_GS has NO reference gate at all here: each box is a new
// calculation and its E_GS becomes the per-box deposit reference, read from
// this run's run_summary.txt by the analysis (never hard-coded).
//
// BULK-LIKENESS is checked DOWNSTREAM, not here: the GS density VTI written to
// results/<cfg>/density_gs is read by hypotheses/lz_bulk_sweep/pilot_gate.py,
// which compares n(z=0) against n0. A 15-Bohr slab holds only ~2 Friedel
// oscillation periods (lambda_F/2 = 6.9 Bohr at r_s = 4.18), so its interior
// may not plateau — that is a WARN with consequences for the 1/L fit, not a
// correctness failure of this GS (checkpoint-dont-block: hard gates are for
// correctness only).
//
// Env: LZB_CFG(REQUIRED) GS_SPACING(0.40) GS_PERIODICITY(2) GS_DIR(derived from
//      the preset's GS_TAG when unset)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/lzb_boxes.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Sh = localised_jellium::config::LzbShared;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
    const std::string REPO = "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research";

    const auto B = localised_jellium::config::lzb_box_from_env();

    const double SPACING     = env_d("GS_SPACING", Sh::SPACING_BOHR);
    const int    PERIODICITY = env_i("GS_PERIODICITY", 2);
    const std::string RES    = "results/" + B.name;
    const std::string GS_DIR = env_s("GS_DIR",
        REPO + "/ResearchProject/systems/localised_jellium/shared_gs/" + B.GS_TAG);

    if (PERIODICITY != 2 && PERIODICITY != 3) {
        std::cerr << "FATAL: GS_PERIODICITY must be 2 or 3, got " << PERIODICITY << "\n";
        return 2;
    }

    const double rs = inqkit::jellium::rs_from_n0(B.n0());

    std::cout << std::setprecision(12)
              << "\n=== lz_bulk_sweep GS [" << B.name << "] ===\n"
              << "  cell        = " << Sh::LX_BOHR << " x " << Sh::LY_BOHR
              << " x " << B.LZ_BOHR << " Bohr, periodicity(" << PERIODICITY << ")\n"
              << "  spacing     = " << SPACING << " Bohr  (nz ~ "
              << int(std::lround(B.LZ_BOHR / SPACING)) << ")\n"
              << "  slab        = half-width " << B.SLAB_HALF
              << " (thickness " << B.l_slab()
              << "), edge " << Sh::EDGE_WIDTH_BOHR << " Bohr\n"
              << "  N_e         = " << B.N_ELECTRONS
              << "  n0 = " << B.n0() << "  r_s = " << rs << "\n"
              << "  extra_states= " << B.EXTRA_STATES
              << "  T = " << Sh::TEMPERATURE_EV << " eV\n"
              << "  CAP (RT)    = bands +/-[" << B.cap_z_inner() << ", "
              << B.LZ_BOHR/2.0 << "] Bohr  (not applied to the GS)\n"
              << "  launch (RT) = z = " << B.LAUNCH_Z
              << "  (standoff " << B.standoff() << " Bohr)\n"
              << "  checkpoint  = " << GS_DIR << "\n\n";

    auto cell0 = systems::cell::orthorhombic(Sh::LX_BOHR * 1.0_b,
                                             Sh::LY_BOHR * 1.0_b,
                                             B.LZ_BOHR * 1.0_b);
    auto cell  = (PERIODICITY == 2) ? cell0.periodicity(2) : cell0.periodic();
    auto ions  = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(B.N_ELECTRONS)
            .extra_states(B.EXTRA_STATES)
            .temperature(Sh::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

    // Positive jellium background as a static perturbation (slab along z).
    inqkit::jellium::localised_background_params bg;
    bg.shape       = inqkit::jellium::background_shape::slab;
    bg.n0          = B.n0();
    bg.half_width  = B.SLAB_HALF;
    bg.slab_axis   = Sh::SLAB_AXIS;
    bg.center      = {0.0, 0.0, Sh::SLAB_CENTER};
    bg.edge_width  = Sh::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation pert(bg);

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions, electrons, options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(Sh::SCF_TOL_HA * 1.0_Ha)
            .max_steps(Sh::SCF_MAX_STEPS)
            .broyden_mixing().mixing_ndim(Sh::SCF_MIX_NDIM).mixing(Sh::SCF_MIX_ALPHA),
        pert);

    const int    n_states = electrons.states().num_states();
    const double E_GS     = gs.energy.total();
    const double n_int    = operations::integral(electrons.density());

    std::cout << "\n  --- gates ---\n" << std::setprecision(12)
              << "  E_GS        = " << E_GS << " Ha   (per-box reference; no gate)\n"
              << "  num_states  = " << n_states << "\n"
              << "  integral n  = " << n_int << "   (ref " << B.N_ELECTRONS << ")\n"
              << "  r_s         = " << rs << "   (ref 4.183)\n";

    int fails = 0;

    // HARD: the electron count. A wrong integral means the background or the
    // extra_electrons wiring is broken.
    if (std::abs(n_int - double(B.N_ELECTRONS)) > 0.01) {
        std::cerr << "  [FAIL] integral n dV = " << n_int << " != " << B.N_ELECTRONS << "\n";
        ++fails;
    } else {
        std::cout << "  [PASS] electron count.\n";
    }

    // HARD: r_s. Pure preset arithmetic — a failure means lzb_boxes.hpp was
    // edited wrongly. Catch it here, not 90 GPU-hours later.
    if (std::abs(rs - 4.183) > 0.005) {
        std::cerr << "  [FAIL] r_s = " << rs << " != 4.183 — the slab density is NOT "
                     "the campaign's. Thickness must scale N, never n0.\n";
        ++fails;
    } else {
        std::cout << "  [PASS] r_s matches every prior slab campaign (n0 unchanged).\n";
    }

    const int expect_states = B.N_ELECTRONS / 2 + B.EXTRA_STATES;
    if (n_states != expect_states)
        std::cout << "  [WARN] num_states = " << n_states << " != " << expect_states
                  << "; the RT binaries' extra_states must match this box.\n";

    if (fails > 0) {
        std::cerr << "\nFATAL: " << fails
                  << " gate(s) failed — refusing to save a GS the sweep would "
                     "silently build on.\n";
        return 3;
    }

    std::filesystem::create_directories(GS_DIR);
    std::filesystem::create_directories(RES);
    electrons.save(GS_DIR);

    // GS density VTI: (a) the n(z) bulk-likeness check in pilot_gate.py reads
    // this; (b) physical order — never fftshift a VTI (vti-coordinate-mapping).
    {
        inqkit::io::RealField3DWriter gs_wr(RES + "/density_gs",
            { .field_name = "density", .include_meta = false, .emit_raw = false,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
            { .overwrite = true });
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs");
    }

    if (electrons.root()) {
        std::ofstream s(RES + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = localised_jellium/lz_bulk_sweep/gs/" << B.name << "\n"
          << "purpose = ground state for the slab->bulk L_slab sweep\n"
          << "plan = docs/plans/jellium-slab-extend-Lz.md\n"
          << "engine = inq-study\nxc = LDA\n"
          << "box_preset = " << B.name << "\n"
          << "cell_bohr = " << Sh::LX_BOHR << "x" << Sh::LY_BOHR << "x" << B.LZ_BOHR
          << " (orthorhombic)\n"
          << "periodicity = " << PERIODICITY << "\n"
          << "slab_half_width = " << B.SLAB_HALF
          << "  slab_thickness = " << B.l_slab()
          << "  edge_width = " << Sh::EDGE_WIDTH_BOHR << "\n"
          << "n0_a0m3 = " << B.n0() << "\nr_s = " << rs << "\n"
          << "spacing_bohr = " << SPACING << "\n"
          << "extra_electrons = " << B.N_ELECTRONS
          << "\nextra_states = " << B.EXTRA_STATES << "\n"
          << "temperature_ev = " << Sh::TEMPERATURE_EV << "\n"
          << "ground_state_energy_ha = " << E_GS << "\n"
          << "electron_count_integral = " << n_int << "\n"
          << "num_states = " << n_states << "\n"
          << "launch_z = " << B.LAUNCH_Z << "  standoff_bohr = " << B.standoff() << "\n"
          << "cap_z_inner = " << B.cap_z_inner()
          << "  cap_width_bohr = " << Sh::CAP_L_BOHR << "\n"
          << "checkpoint_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. Checkpoint saved to " << GS_DIR << "\n";
    return 0;
}
