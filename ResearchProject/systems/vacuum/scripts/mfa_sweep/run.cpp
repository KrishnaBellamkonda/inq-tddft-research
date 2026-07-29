// ============================================================================
// mfa_sweep: one env-driven MFA reflection run for the ε(E,L) parameter study.
//
// Built ONCE; invoked many times by dispatch_sweep.py with different env vars.
// Implements the sin² mask absorber of De Giovannini, Larsen & Rubio (2014),
// arXiv:1409.1689, Eq. 12–13, ENTIRELY in the inq-stack wrapper (inqkit::
// absorbers::MaskAbsorber applied in the per-step propagate callback). inq/ and
// inq-study are untouched.
//
// Geometry (plan §2; propagation along z):
//   σ = 4√2/k₀                (NEVER a free parameter — pins the ε formula)
//   Lcell_z = 6σ + L,  z_abs0 = (6σ − L)/2,  z0 = −L/2 (WP launch)
//   absorber [z_abs0, z_abs0+L] (= up to the right cell edge), inner region z<z_abs0
//   τ = 2(3σ+L)/k₀,  N_STEPS = round(τ/dt)
//   ε = ∫_{z<z_abs0}|ψ_wp(τ)|² / N₀   (N₀ = full WP norm at t=0; cancels the
//       transverse-box normalisation so a minimal transverse box is exact)
//
// Run types:
//   masked  (MFA_ANCHOR=0): periodic cell, mask applied every step.
//   anchor  (MFA_ANCHOR=1): FINITE cell (hard wall), NO mask → ε≈1 reference.
//
// Env vars: MFA_K0, MFA_LABS, MFA_ANCHOR, MFA_SHOWCASE, MFA_OUTDIR,
//           MFA_NPERP(=12), MFA_DX(=0.1), MFA_DT(=0.01), MFA_BC(periodic|finite)
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/absorbers/mask_absorber.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
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

static double env_d(const char *k, double def) {
  const char *v = std::getenv(k);
  return v ? std::atof(v) : def;
}
static int env_i(const char *k, int def) {
  const char *v = std::getenv(k);
  return v ? std::atoi(v) : def;
}
static std::string env_s(const char *k, const std::string &def) {
  const char *v = std::getenv(k);
  return v ? std::string(v) : def;
}

int main() {
  const double HA_TO_EV = 27.211386245988;

  const double k0 = env_d("MFA_K0", 1.0);
  const double Labs = env_d("MFA_LABS", 20.0);
  const bool anchor = env_i("MFA_ANCHOR", 0) != 0;
  const bool showcase = env_i("MFA_SHOWCASE", 0) != 0;
  const std::string outdir = env_s("MFA_OUTDIR", "results");
  const int nperp = env_i("MFA_NPERP", 12);
  // Adaptive grid spacing: the WP density std is 4/k0, so low-k0 (wide) packets
  // are hugely over-resolved at a fixed dx — use dx≈0.75/k0 (resolves both the
  // packet ~0.8/k0 and its k-content ~0.75/k0), clamped [0.18, 0.30]. The 0.18
  // floor keeps the grid's max kinetic eigenvalue E_max=½(π/dx)²≤152 Ha so the
  // ETRS Taylor propagator stays stable (dt·E_max≤1.5; dx=0.15 aborts at ~2.2).
  // Env MFA_DX overrides (used by the gate-2 dx-convergence check).
  const double dx_auto = std::min(0.30, std::max(0.18, 0.75 / k0));
  const double dx = env_d("MFA_DX", dx_auto);
  const double dt = env_d("MFA_DT", 0.01);
  const std::string bc = env_s("MFA_BC", anchor ? "finite" : "periodic");

  // ---- geometry (σ pinned to 4√2/k₀) -------------------------------------
  const double sigma = 4.0 * std::sqrt(2.0) / k0;
  const double Lcell_z = 6.0 * sigma + Labs;
  const double z_abs0 = (6.0 * sigma - Labs) / 2.0;
  const double z0 = -Labs / 2.0;
  const double tau = 2.0 * (3.0 * sigma + Labs) / k0;
  const int N_STEPS = std::max(1, (int)std::llround(tau / dt));
  const double Lperp = nperp * dx;
  const double ec = 0.5 * std::pow(M_PI / dx, 2.0);
  const double E_eV = 0.5 * k0 * k0 * HA_TO_EV;

  fs::create_directories(outdir);
  fs::create_directories(outdir + "/raw");
  std::printf("\n=== mfa_sweep: k0=%.4f  E=%.3f eV  L=%.1f  %s  %s ===\n", k0,
              E_eV, Labs, anchor ? "ANCHOR" : "masked", bc.c_str());
  std::printf("  sigma=%.3f Lcell_z=%.2f z_abs0=%.3f z0=%.3f tau=%.2f "
              "N_STEPS=%d Lperp=%.2f\n",
              sigma, Lcell_z, z_abs0, z0, tau, N_STEPS, Lperp);

  // ---- system (free particle: empty ions, ghost occupied, non_interacting) -
  auto base = systems::cell::orthorhombic(Lperp * 1.0_b, Lperp * 1.0_b,
                                          Lcell_z * 1.0_b);
  auto cell = (bc == "finite") ? base.finite() : base.periodic();
  auto ions = systems::ions(cell);
  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(ec * 1.0_Ha).extra_states(1).extra_electrons(2.0));
  ground_state::initial_guess(ions, electrons);

  auto rep = inqkit::WavePacket{}
                 .center(0.0, 0.0, z0)
                 .sigma(sigma)
                 .k0(0.0, 0.0, k0)
                 .inject_into_last_extra_state(electrons, 1.0);
  const long wp_idx = rep.state_index;
  std::printf("  WP injected: state_index=%ld norm_after=%.6f\n", wp_idx,
              rep.norm_after);

  // N₀ = full WP norm at t=0 (normalisation that makes ε transverse-box-robust)
  const double N0 = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);

  abs_::MaskAbsorber absb(2, z_abs0, Labs, wp_idx);

  // ---- cadences ----------------------------------------------------------
  const int obs_every = std::max(1, N_STEPS / 100);
  const int inner_every = std::max(1, N_STEPS / 200);
  const int wf_every = std::max(1, N_STEPS / 60);

  // ---- observers ---------------------------------------------------------
  inqkit::io::ObservableSelection sel;
  sel.step = sel.time_au = true;
  sel.energy_total = sel.energy_kinetic = true;
  sel.current_x = sel.current_y = sel.current_z = true;
  sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
  inqkit::io::ObservablesWriter obs_writer(outdir + "/raw/observables.csv", sel);
  obs_writer.write_header();

  std::ofstream inner_csv(outdir + "/raw/inner_norm_vs_time.csv");
  inner_csv << "step,time_au,inner_norm_over_N0\n";

  inqkit::observables::WPRealSpaceStats wp_rs(
      outdir + "/raw/wp_real_space_stats.csv", wp_idx,
      {.write_every = obs_every});
  inqkit::observables::WPMomentumStats wp_mom(
      outdir + "/raw/wp_momentum_stats.csv", wp_idx, {.write_every = obs_every});
  // momentum distribution at t=0 and t=τ only (before/after excitation metric)
  inqkit::observables::MomentumDistribution mom_dist(
      outdir + "/raw/momentum_distribution.csv", wp_idx, Lcell_z,
      {.n_bins = 64, .k_max_bohr_inv = 0.0, .write_every = N_STEPS});

  // showcase-only writers
  std::unique_ptr<inqkit::io::RealField3DWriter> tot_wr, sys_wr, wp_wr;
  std::unique_ptr<inqkit::io::ComplexField3DWriter> wf_wr;
  if (showcase) {
    fs::create_directories(outdir + "/raw/vti");
    using RLay = inqkit::io::RealField3DLayout;
    using ROpt = inqkit::io::RealField3DWriteOptions;
    const auto vti = inqkit::io::VTIWriteOptions::Format::binary;
    RLay lay_tot{.field_name = "density_total", .include_meta = true, .emit_vti = true, .vti_format = vti};
    RLay lay_sys{.field_name = "density_system", .include_meta = true, .emit_vti = true, .vti_format = vti};
    RLay lay_wp{.field_name = "density_wp", .include_meta = true, .emit_vti = true, .vti_format = vti};
    ROpt ow{.overwrite = true};
    tot_wr = std::make_unique<inqkit::io::RealField3DWriter>(outdir + "/raw/vti/density_total", lay_tot, ow);
    sys_wr = std::make_unique<inqkit::io::RealField3DWriter>(outdir + "/raw/vti/density_system", lay_sys, ow);
    wp_wr = std::make_unique<inqkit::io::RealField3DWriter>(outdir + "/raw/vti/density_wp", lay_wp, ow);
    inqkit::io::ComplexField3DLayout lay_wf{.field_name = "wavefunction_wp", .include_meta = false, .emit_raw = false, .emit_vti = true, .vti_format = vti};
    inqkit::io::ComplexField3DWriteOptions cow{.overwrite = true};
    wf_wr = std::make_unique<inqkit::io::ComplexField3DWriter>(outdir + "/raw/vti/wavefunction_wp", lay_wf, cow);
  }

  // ---- propagate: mask each step (after U), then observe -----------------
  real_time::propagate(
      ions, electrons,
      [&](auto const &data) {
        const int step = data.iter();
        // Eq. 12: apply M AFTER the ETRS step (iter>0); skip the iter=0 initial
        // observer call (no U has acted yet).
        if (!anchor && step > 0) absb.apply(electrons);

        if (step % obs_every == 0) {
          inqkit::StepContext ctx;
          ctx.step = step;
          ctx.time_au = data.time();
          ctx.ions = &ions;
          ctx.electrons = &electrons;
          ctx.energy_total = data.energy().total();
          ctx.energy_kinetic = data.energy().kinetic();
          try { auto c = data.current(); ctx.current = {c[0], c[1], c[2]}; } catch (...) {}
          try { auto d = data.dipole();  ctx.dipole  = {d[0], d[1], d[2]}; } catch (...) {}
          obs_writer.append(ctx);
        }
        if (step % inner_every == 0) {
          double e = abs_::inner_region_norm(electrons, 2, z_abs0, wp_idx) / N0;
          inner_csv << step << ',' << data.time() << ',' << e << '\n';
        }
        if (showcase && step % wf_every == 0) {
          auto tot = inqkit::fields::density::total(electrons);
          auto wp = inqkit::fields::density::orbital(electrons, wp_idx);
          auto sys = inqkit::fields::density::total_excluding_orbital(
              electrons, wp_idx, 1.0);
          tot_wr->write(tot, data.time(), step);
          sys_wr->write(sys, data.time(), step);
          wp_wr->write(wp, data.time(), step);
          char nm[64];
          std::snprintf(nm, sizeof(nm), "wf_t%06d", step);
          auto wf = inqkit::fields::orbital::wavefunction(electrons, wp_idx);
          wf_wr->write(wf, std::string(nm));
        }
        mom_dist.maybe_accumulate(data);
        wp_mom.maybe_accumulate(data);
        wp_rs.maybe_accumulate(data);
      },
      options::theory{}.non_interacting(),
      // ETRS (default): it does NOT renormalise the orbital each step, so the
      // mask's norm reduction (absorption) is preserved — Crank–Nicolson
      // re-normalises the WP to unit norm every step, which silently undoes the
      // absorption. Stability is ensured by the dx≥0.16 clamp above.
      options::real_time{}.num_steps(N_STEPS).dt(dt * 1.0_atomictime));

  inner_csv.close();

  // ---- ε at τ ------------------------------------------------------------
  const double inner_tau = abs_::inner_region_norm(electrons, 2, z_abs0, wp_idx);
  const double eps = inner_tau / N0;

  std::ofstream eps_f(outdir + "/epsilon.txt");
  eps_f << "epsilon " << eps << "\n"
        << "inner_norm_tau " << inner_tau << "\n"
        << "N0 " << N0 << "\n"
        << "k0 " << k0 << "\n"
        << "E_eV " << E_eV << "\n"
        << "sigma " << sigma << "\n"
        << "L_abs " << Labs << "\n"
        << "z_abs0 " << z_abs0 << "\n"
        << "z0 " << z0 << "\n"
        << "Lcell_z " << Lcell_z << "\n"
        << "tau " << tau << "\n"
        << "N_STEPS " << N_STEPS << "\n"
        << "dt " << dt << "\n"
        << "dx " << dx << "\n"
        << "nperp " << nperp << "\n"
        << "anchor " << (anchor ? 1 : 0) << "\n"
        << "showcase " << (showcase ? 1 : 0) << "\n"
        << "bc " << bc << "\n";
  eps_f.close();

  std::ofstream sum_f(outdir + "/run_summary.txt");
  sum_f << "run_completed = true\n"
        << "epsilon = " << eps << "\n"
        << "norm_after = " << rep.norm_after << "\n"
        << "N0 = " << N0 << "\n"
        << "E_eV = " << E_eV << "\n"
        << "L_abs = " << Labs << "\n"
        << "k0 = " << k0 << "\n"
        << "N_STEPS = " << N_STEPS << "\n";
  sum_f.close();

  std::printf("  epsilon = %.6f  (inner_tau=%.6f / N0=%.6f)\n", eps, inner_tau,
              N0);
  std::printf("Done. Output in %s/\n", outdir.c_str());
  return 0;
}
