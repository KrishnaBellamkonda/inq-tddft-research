// ============================================================================
// twosided_cap_vs_mask: one env-driven run of a TWO-SIDED absorber on a free
// wavepacket in vacuum, emitting the FULL free-WP MINIMUM OBSERVABLE SET +
// manifest (ADR 0006). Plan: docs/plans/twosided-cap-vs-mask.md.
//
// Two absorber modes (CAP_MODE = "cap" | "mask"), same geometry & ε:
//   cap  : built-in sin² CAP on BOTH z-ends, composed as
//            absorbing(η, +mid, w) + absorbing(η, −mid, w)   (perturbations::sum)
//          integrated in H (ETRS). NO new engine code — inq/ AND inq-study are
//          untouched; this only #includes existing headers. Build vs inq-study.
//   mask : inqkit::absorbers::TwoSidedMaskAbsorber applied in the propagate
//          callback after each ETRS step (wrapper-only; no engine dependency).
//
// Geometry (two-sided, energy-scaled, quasi-monochromatic packet):
//   σ = 4√2/k0;  z_in = 6σ (inner half-width);  Lhalf = L_total/2 (per-end);
//   Lcell_z = 2·z_in + L = 12σ + L;  WP launched at z0 = −z_in + 4σ (≥4σ from the
//   near absorber ⇒ negligible density in the CAP region), moving +z toward the
//   far absorber. ε = ∫_{|z|<z_in}|ψ_wp(τ)|² / N0  (surviving un-absorbed fraction).
//
// Env: CAP_MODE(cap), CAP_K0, CAP_L(total), CAP_ETA(<0 Ha, cap only),
//      CAP_OUTDIR, CAP_NPERP(8), CAP_DT(0.01).
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/absorbers/mask_absorber.hpp>   // TwoSidedMaskAbsorber, inner_region_norm[_twosided]
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/eigenvalue_dump.hpp>
#include <inqkit/observables/minimum_observable_set.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;
namespace abs_ = inqkit::absorbers;
namespace obs_ = inqkit::observables;

static double env_d(const char *k, double d) { const char *v = std::getenv(k); return v ? std::atof(v) : d; }
static int    env_i(const char *k, int d)    { const char *v = std::getenv(k); return v ? std::atoi(v) : d; }
static std::string env_s(const char *k, const std::string &d) { const char *v = std::getenv(k); return v ? std::string(v) : d; }

int main() {
  const double HA_TO_EV = 27.211386245988;
  const std::string mode = env_s("CAP_MODE", "cap");   // "cap" | "mask"
  const bool mask_mode = (mode == "mask");
  const double k0   = env_d("CAP_K0", 1.0);
  const double Ltot = env_d("CAP_L", 20.0);            // TOTAL absorber width
  const double eta  = env_d("CAP_ETA", -0.5);          // CAP depth (cap mode only)
  const int nperp   = env_i("CAP_NPERP", 8);
  const double dt   = env_d("CAP_DT", 0.01);
  const std::string outdir = env_s("CAP_OUTDIR", "results");

  const double sigma   = 4.0 * std::sqrt(2.0) / k0;
  const double Lhalf   = Ltot / 2.0;                   // per-end absorber width
  const double z_in    = 6.0 * sigma;                  // inner half-width
  const double Lcell_z = 2.0 * z_in + Ltot;            // = 12σ + L
  const double z0      = -z_in + 4.0 * sigma;          // launch ≥4σ from near absorber
  const double tau     = 2.0 * (z_in + Ltot) / k0;     // cross inner region + absorb
  const double dx      = std::min(0.30, std::max(0.18, 0.75 / k0));
  const int N_STEPS    = std::max(1, (int)std::llround(tau / dt));
  const double Lperp   = nperp * dx;
  const double ec      = 0.5 * std::pow(M_PI / dx, 2.0);
  const double E_eV    = 0.5 * k0 * k0 * HA_TO_EV;
  const double width_frac = Lhalf / Lcell_z;                       // per-end slab, fractional
  const double mid_frac   = (z_in + Lhalf / 2.0) / Lcell_z;        // high-slab centre, fractional

  fs::create_directories(outdir + "/raw/observables/eigenvalues");
  fs::create_directories(outdir + "/raw/vti");
  const int obs_every = std::max(1, N_STEPS / 100);
  const int wf_every  = std::max(1, N_STEPS / 60);

  obs_::write_manifest(outdir + "/observables_manifest.json",
                       obs_::RunType::free_wp, obs_every, N_STEPS);

  std::printf("\n=== twosided[%s]: k0=%.4f E=%.3f eV L=%.1f eta=%.3f Ha ===\n",
              mode.c_str(), k0, E_eV, Ltot, eta);
  std::printf("  sigma=%.3f Lcell=%.2f z_in=%.3f Lhalf=%.2f z0=%.3f tau=%.2f N=%d wfrac=%.4f midfrac=%.4f\n",
              sigma, Lcell_z, z_in, Lhalf, z0, tau, N_STEPS, width_frac, mid_frac);

  auto cell = systems::cell::orthorhombic(Lperp * 1.0_b, Lperp * 1.0_b, Lcell_z * 1.0_b).periodic();
  auto ions = systems::ions(cell);
  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(ec * 1.0_Ha).extra_states(1).extra_electrons(2.0));
  ground_state::initial_guess(ions, electrons);

  auto rep = inqkit::WavePacket{}.center(0.0, 0.0, z0).sigma(sigma).k0(0.0, 0.0, k0)
                 .inject_into_last_extra_state(electrons, 1.0);
  const long wp_idx = rep.state_index;
  const double N0   = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);
  std::printf("  WP injected: state_index=%ld N0=%.6f\n", wp_idx, N0);

  // GS system density VTI (core), once at t=0.
  {
    using RLay = inqkit::io::RealField3DLayout;
    using ROpt = inqkit::io::RealField3DWriteOptions;
    const auto vti = inqkit::io::VTIWriteOptions::Format::binary;
    RLay lay{.field_name = "density_gs_system", .include_meta = true, .emit_raw = true, .emit_vti = true, .vti_format = vti};
    ROpt ow{.overwrite = true};
    inqkit::io::RealField3DWriter gw(outdir + "/raw/vti/density_gs_system", lay, ow);
    auto sys0 = inqkit::fields::density::total_excluding_orbital(electrons, wp_idx, 1.0);
    gw.write(sys0, "density_gs_system");
  }

  // Two-sided CAP = sum of two built-in absorbing slabs (composition only).
  auto cap = perturbations::absorbing(eta * 1.0_Ha,  mid_frac, width_frac)
           + perturbations::absorbing(eta * 1.0_Ha, -mid_frac, width_frac);
  // Two-sided mask (applied in the callback; constructed in both modes, used in mask mode).
  abs_::TwoSidedMaskAbsorber mask(2, z_in, Lhalf, wp_idx);

  // ---- observable writers (paths exactly as the manifest declares) -------
  inqkit::io::ObservableSelection sel;
  sel.step = sel.time_au = true;
  sel.energy_total = sel.energy_kinetic = true;
  sel.energy_hartree = sel.energy_xc = true;
  sel.density_l2 = true;
  sel.current_x = sel.current_y = sel.current_z = true;
  sel.dipole_z = true;
  inqkit::io::ObservablesWriter obs_writer(outdir + "/raw/observables/observables.csv", sel);
  obs_writer.write_header();

  obs_::WPRealSpaceStats wp_rs(outdir + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every = obs_every});
  obs_::WPMomentumStats  wp_mom(outdir + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every = obs_every});
  obs_::MomentumDistribution mom_dist(outdir + "/raw/observables/momentum_distribution.csv", wp_idx, Lcell_z,
                                      {.n_bins = 64, .k_max_bohr_inv = 0.0, .write_every = obs_every});
  obs_::DensityDelta density_delta(outdir + "/raw/vti/density_delta",
                                   outdir + "/raw/vti/density_delta_coarse",
                                   {.compute_l2 = true, .coarse_bin_bohr = 3.0});

  std::unique_ptr<inqkit::io::RealField3DWriter> wp_wr;
  {
    using RLay = inqkit::io::RealField3DLayout;
    using ROpt = inqkit::io::RealField3DWriteOptions;
    const auto vti = inqkit::io::VTIWriteOptions::Format::binary;
    RLay lay_wp{.field_name = "density_wp", .include_meta = true, .emit_raw = true, .emit_vti = true, .vti_format = vti};
    ROpt ow{.overwrite = true};
    wp_wr = std::make_unique<inqkit::io::RealField3DWriter>(outdir + "/raw/vti/density_wp", lay_wp, ow);
  }

  std::ofstream inner_csv(outdir + "/raw/observables/inner_norm_vs_time.csv");
  inner_csv << "step,time_au,total_wp_norm,inner_norm_over_N0\n";

  obs_::OccupationsWriter occ(outdir + "/raw/observables/eigenvalues/occupations.csv");

  // Shared per-step callback. In mask mode the absorber is applied AFTER the ETRS
  // step (step>0), mutating the captured non-const electrons (paper Eq. 12).
  auto step_cb = [&](auto const &data) {
    const int step = data.iter();
    if (mask_mode && step > 0) mask.apply(electrons);
    if (step % obs_every == 0) {
      inqkit::StepContext ctx;
      ctx.step = step; ctx.time_au = data.time();
      ctx.ions = &ions; ctx.electrons = &electrons;
      ctx.energy_total = data.energy().total();
      ctx.energy_kinetic = data.energy().kinetic();
      ctx.energy_hartree = data.energy().hartree();
      ctx.energy_xc = data.energy().xc();
      try { auto c = data.current(); ctx.current = {c[0], c[1], c[2]}; } catch (...) {}
      try { auto dp = data.dipole(); ctx.dipole = {dp[0], dp[1], dp[2]}; } catch (...) {}
      auto sys_f = inqkit::fields::density::total_excluding_orbital(electrons, wp_idx, 1.0);
      ctx.density_l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
      obs_writer.append(ctx);

      double tot = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);
      double inn = abs_::inner_region_norm_twosided(electrons, 2, z_in, wp_idx) / N0;
      inner_csv << step << ',' << data.time() << ',' << tot << ',' << inn << '\n';
    }
    if (step % wf_every == 0) {
      auto wp = inqkit::fields::density::orbital(electrons, wp_idx);
      wp_wr->write(wp, data.time(), step);
    }
    if (step == N_STEPS) occ.snapshot(data);
    mom_dist.maybe_accumulate(data);
    wp_mom.maybe_accumulate(data);
    wp_rs.maybe_accumulate(data);
  };

  auto theory = options::theory{}.non_interacting();
  auto rt = options::real_time{}.num_steps(N_STEPS).dt(dt * 1.0_atomictime);   // ETRS

  if (mask_mode) {
    real_time::propagate(ions, electrons, step_cb, theory, rt);                // mask in callback
  } else {
    real_time::propagate(ions, electrons, step_cb, theory, rt, cap);           // CAP in H
  }

  inner_csv.close();
  obs_::dump_eigenvalues(electrons, outdir + "/raw/observables/eigenvalues/eigenvalues.csv");

  // ---- ε + provenance -----------------------------------------------------
  const double inner_tau = abs_::inner_region_norm_twosided(electrons, 2, z_in, wp_idx);
  const double total_tau = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);
  const double eps = inner_tau / N0;

  std::ofstream eps_f(outdir + "/epsilon.txt");
  eps_f << "epsilon " << eps << "\ninner_norm_tau " << inner_tau
        << "\ntotal_wp_norm_tau " << total_tau << "\nabsorbed_fraction " << (1.0 - total_tau / N0)
        << "\nN0 " << N0 << "\nk0 " << k0 << "\nE_eV " << E_eV << "\nsigma " << sigma
        << "\nmode " << mode << "\nL_total " << Ltot << "\nLhalf " << Lhalf
        << "\neta_Ha " << (mask_mode ? 0.0 : eta) << "\nz_in " << z_in << "\nz0 " << z0
        << "\nLcell_z " << Lcell_z << "\nmid_frac " << mid_frac << "\nwidth_frac " << width_frac
        << "\ntau " << tau << "\nN_STEPS " << N_STEPS << "\ndt " << dt << "\ndx " << dx
        << "\nnperp " << nperp << "\npropagator etrs\n";
  eps_f.close();

  std::ofstream sum_f(outdir + "/run_summary.txt");
  sum_f << "run_completed = true\nmode = " << mode << "\nepsilon = " << eps
        << "\nabsorbed_fraction = " << (1.0 - total_tau / N0) << "\nN0 = " << N0
        << "\nE_eV = " << E_eV << "\nL_total = " << Ltot << "\neta_Ha = " << (mask_mode ? 0.0 : eta)
        << "\nk0 = " << k0 << "\nN_STEPS = " << N_STEPS << "\n";
  sum_f.close();

  std::printf("  [%s] epsilon=%.6f absorbed=%.4f -> %s/\n", mode.c_str(), eps, 1.0 - total_tau / N0, outdir.c_str());
  return 0;
}
