// ============================================================================
// run_01_base: free wavepacket propagation in a finite box
//
// System: empty finite cell 34.771 × 34.771 × 89.856 bohr, no ions.
//         extra_states(1) → 0 occupied + 1 WP slot (index 0).
//         GS trivially empty; max_steps=10 for API compliance.
//
// Wavepacket (base run):
//   sigma = 0.53 Å = 1.002 bohr
//   E_kin = 200 eV → k0 = 3.834 bohr⁻¹
//   center: (Lx/2, Ly/2, Lz − 5σ) → moving in −z
//
// TDDFT: dt=0.02 a.u., 10000 steps ≈ 4.83 fs
//   density written every 100 steps (100 frames)
//   observables (all) every step
//   3 screen accumulators at Lz/4, Lz/2, 3Lz/4
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;

static constexpr double ANG_TO_BOHR = 1.8897259886;
static constexpr double HA_TO_EV = 27.21138625;

static constexpr double LX_BOHR = 34.771;
static constexpr double LY_BOHR = 34.771;
static constexpr double LZ_BOHR = 89.856;

static constexpr double WP_SIGMA_ANG = 0.53;
static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;
static constexpr double WP_EKIN_EV = 200.0;
static constexpr double WP_EKIN_HA = WP_EKIN_EV / HA_TO_EV;
static const double WP_K0 = std::sqrt(2.0 * WP_EKIN_HA);

static const double WP_CX = LX_BOHR / 2.0;
static const double WP_CY = LY_BOHR / 2.0;
static const double WP_CZ = LZ_BOHR - 5.0 * WP_SIGMA_BOHR;

static constexpr int N_STEPS = 10000;
static constexpr double DT_AU = 0.02;
static constexpr int WRITE_EVERY = 100;

static void add_field_inplace(inqkit::fields::RealField3D &a,
                              inqkit::fields::RealField3D const &b) {
  for (std::size_t i = 0; i < a.values.size(); i++)
    a.values[i] += b.values[i];
}

int main() {
  std::cout << "\n=== run_01_base: free WP propagation, 200 eV, sigma=0.53A, "
               "-z ===\n";

  auto cell = systems::cell::orthorhombic(LX_BOHR * 1.0_b, LY_BOHR * 1.0_b,
                                          LZ_BOHR * 1.0_b)
                  .finite();
  auto ions = systems::ions(cell);

  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(40.0_Ha).extra_states(1));

  ground_state::initial_guess(ions, electrons);
  auto gs = ground_state::calculate(
      ions, electrons, options::theory{}.lda(),
      options::ground_state{}.energy_tolerance(1e-4_Ha).max_steps(10));

  std::cout << "  GS energy = " << gs.energy.total()
            << " Ha (trivially zero)\n";

  // WP injection
  auto wp = inqkit::WavePacket{}
                .center(WP_CX, WP_CY, WP_CZ)
                .sigma(WP_SIGMA_BOHR)
                .k0(0.0, 0.0, -WP_K0);
  // No orthogonalisation: 0 occupied states
  auto report = wp.inject_into_last_extra_state(electrons, 1.0);

  std::cout << "  state_index = " << report.state_index
            << "  norm_after = " << report.norm_after << "  (expect ≈ 1.0)\n";

  // t=0 density: total() is zero + WP orbital
  inqkit::io::RealField3DWriter density_writer(
      "results/density_rt", {.field_name = "density", .include_meta = true},
      {.overwrite = true});

  auto rho_t0 = inqkit::fields::density::total(electrons);
  auto rho_wp = inqkit::fields::density::orbital(electrons, report.state_index);
  add_field_inplace(rho_t0, rho_wp);
  density_writer.write(rho_t0, 0.0, 0);

  // Observables writer — all columns
  inqkit::io::ObservableSelection sel;
  sel.step = sel.time_au = true;
  sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc =
      true;
  sel.current_x = sel.current_y = sel.current_z = true;
  sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
  inqkit::io::ObservablesWriter obs_writer("results/observables.csv", sel);
  obs_writer.write_header();

  // Screen accumulators at Lz/4, Lz/2, 3Lz/4
  inqkit::screens::LeedPatternAccumulator sc1(
      inqkit::screens::PlaneScreen{LZ_BOHR * 0.25, "screen_z25"});
  inqkit::screens::LeedPatternAccumulator sc2(
      inqkit::screens::PlaneScreen{LZ_BOHR * 0.50, "screen_z50"});
  inqkit::screens::LeedPatternAccumulator sc3(
      inqkit::screens::PlaneScreen{LZ_BOHR * 0.75, "screen_z75"});

  // RT session
  inqkit::RealTimeSession rt(ions, electrons, WRITE_EVERY);
  rt.add([&](inqkit::StepContext const &ctx) {
    auto rho = inqkit::fields::density::total(*ctx.electrons);
    density_writer.write(rho, ctx.time_au, ctx.step);
  });

  // Observables and screens run every step
  inqkit::RealTimeSession rt_obs(ions, electrons, 1);
  rt_obs.add([&](inqkit::StepContext const &ctx) {
    obs_writer.append(ctx);
    sc1.accumulate(*ctx.electrons, DT_AU);
    sc2.accumulate(*ctx.electrons, DT_AU);
    sc3.accumulate(*ctx.electrons, DT_AU);
  });

  real_time::propagate(
      ions, electrons,
      [&](auto const &data) {
        rt.step(data);
        rt_obs.step(data);
      },
      options::theory{}.lda(),
      options::real_time{}
          .num_steps(N_STEPS)
          .dt(DT_AU * 1.0_atomictime)
          .observables_current()
          .observables_dipole());

  std::filesystem::create_directories("results/screens");
  sc1.save("results/screens/screen_z25.dat");
  sc2.save("results/screens/screen_z50.dat");
  sc3.save("results/screens/screen_z75.dat");

  std::cout << "Done. Output in results/\n";
  return 0;
}
