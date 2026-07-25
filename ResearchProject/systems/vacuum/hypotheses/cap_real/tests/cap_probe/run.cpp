// ============================================================================
// cap_probe: does INQ's built-in perturbations::absorbing (a sin^2 CAP) actually
// absorb, and which propagator tolerates the resulting non-Hermitian H?
//
// This is a MECHANISM probe (ADR 0007: hypotheses/<NN>/tests/), not a sweep. It
// drives the SAME free-WP geometry as the MFA study, but instead of masking in
// the inq-stack callback it passes perturbations::absorbing as the trailing
// real_time::propagate argument — i.e. the OUT-OF-THE-BOX engine CAP, inq/ and
// inq-study untouched.
//
//   perturbations::absorbing(amplitude<energy>, mid_pos_frac, width_frac)
//     adds  V = +i*amplitude*sin^2(...)  in the FRACTIONAL z-slab
//     [mid_pos-width/2, mid_pos+width/2]   (absorbing.hpp:44-46).
//   amplitude < 0 absorbs:  exp(-iVt) = exp(amplitude*sin^2*t).
//
// Geometry (matches the MFA runs so epsilon is directly comparable):
//   sigma=4sqrt2/k0, Lcell_z=6sigma+L, centred cell z in [-Lcell/2,+Lcell/2],
//   WP launched at z0=-L/2, absorber Cartesian [(6sigma-L)/2,(6sigma+L)/2] which
//   is fractional [0.5-width,0.5] with width=L/Lcell, mid_pos=0.5-width/2.
//   inner region (epsilon, Eq.7): z < z_abs0=(6sigma-L)/2.
//
// Diagnoses two distinct failure modes:
//   - NO-OP   : if even a near-full-cell CAP (CAP_WIDE=1) leaves |psi| unchanged,
//               the imaginary potential is being dropped (vscalar is real).
//   - MISPLACE: if the wide CAP absorbs but the MFA-placed slab does not, the
//               fractional placement / unit convention is wrong.
//
// Env: CAP_K0(1.0) CAP_L(20) CAP_ETA(-0.5 Ha) CAP_PROP(etrs|cn) CAP_WIDE(0)
//      CAP_NPERP(8) CAP_DT(0.01) CAP_OUTDIR(results)
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/absorbers/mask_absorber.hpp>
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

static double env_d(const char *k, double d) { const char *v = std::getenv(k); return v ? std::atof(v) : d; }
static int    env_i(const char *k, int d)    { const char *v = std::getenv(k); return v ? std::atoi(v) : d; }
static std::string env_s(const char *k, const std::string &d) { const char *v = std::getenv(k); return v ? std::string(v) : d; }

int main() {
  const double HA_TO_EV = 27.211386245988;
  const double k0   = env_d("CAP_K0", 1.0);
  const double Labs = env_d("CAP_L", 20.0);
  const double eta  = env_d("CAP_ETA", -0.5);          // Ha; <0 absorbs
  const std::string prop = env_s("CAP_PROP", "etrs");
  const bool wide   = env_i("CAP_WIDE", 0) != 0;
  const int nperp   = env_i("CAP_NPERP", 8);
  const double dt   = env_d("CAP_DT", 0.01);
  const std::string outdir = env_s("CAP_OUTDIR", "results");

  const double sigma   = 4.0 * std::sqrt(2.0) / k0;
  const double Lcell_z = 6.0 * sigma + Labs;
  const double z_abs0  = (6.0 * sigma - Labs) / 2.0;
  const double z0      = -Labs / 2.0;
  const double tau     = 2.0 * (3.0 * sigma + Labs) / k0;
  const double dx      = std::min(0.30, std::max(0.18, 0.75 / k0));
  const int N_STEPS    = std::max(1, (int)std::llround(tau / dt));
  const double Lperp   = nperp * dx;
  const double ec      = 0.5 * std::pow(M_PI / dx, 2.0);
  const double E_eV    = 0.5 * k0 * k0 * HA_TO_EV;

  // CAP slab in FRACTIONAL coords. MFA-placed: last L of the cell. Wide: ~full.
  const double width_frac = wide ? 0.90 : (Labs / Lcell_z);
  const double mid_frac   = wide ? 0.0  : (0.5 - width_frac / 2.0);

  fs::create_directories(outdir);
  std::printf("\n=== cap_probe: k0=%.3f E=%.2f eV L=%.1f eta=%.3f Ha prop=%s %s ===\n",
              k0, E_eV, Labs, eta, prop.c_str(), wide ? "WIDE" : "MFA-slab");
  std::printf("  sigma=%.3f Lcell=%.2f z_abs0=%.3f tau=%.2f N=%d width_frac=%.4f mid_frac=%.4f\n",
              sigma, Lcell_z, z_abs0, tau, N_STEPS, width_frac, mid_frac);

  auto cell = systems::cell::orthorhombic(Lperp * 1.0_b, Lperp * 1.0_b, Lcell_z * 1.0_b).periodic();
  auto ions = systems::ions(cell);
  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(ec * 1.0_Ha).extra_states(1).extra_electrons(2.0));
  ground_state::initial_guess(ions, electrons);

  auto rep = inqkit::WavePacket{}.center(0.0, 0.0, z0).sigma(sigma).k0(0.0, 0.0, k0)
                 .inject_into_last_extra_state(electrons, 1.0);
  const long wp_idx = rep.state_index;
  const double N0   = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);  // total WP norm t=0
  std::printf("  WP injected: state_index=%ld N0=%.6f\n", wp_idx, N0);

  perturbations::absorbing cap(eta * 1.0_Ha, mid_frac, width_frac);

  std::ofstream trace(outdir + "/norm_vs_time.csv");
  trace << "step,time_au,total_wp_norm,inner_wp_norm\n";
  const int every = std::max(1, N_STEPS / 100);

  auto observer = [&](auto const &data) {
    const int step = data.iter();
    if (step % every == 0) {
      double tot = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);
      double inn = abs_::inner_region_norm(electrons, 2, z_abs0, wp_idx);
      trace << step << ',' << data.time() << ',' << tot << ',' << inn << '\n';
    }
  };

  auto rt = options::real_time{}.num_steps(N_STEPS).dt(dt * 1.0_atomictime);
  if (prop == "cn") rt = rt.crank_nicolson();   // else default ETRS
  real_time::propagate(ions, electrons, observer, options::theory{}.non_interacting(), rt, cap);
  trace.close();

  const double total_tau = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);
  const double inner_tau = abs_::inner_region_norm(electrons, 2, z_abs0, wp_idx);
  const double eps = inner_tau / N0;
  const double absorbed_frac = 1.0 - total_tau / N0;

  std::ofstream f(outdir + "/probe_result.txt");
  f << "k0 " << k0 << "\nE_eV " << E_eV << "\nL_abs " << Labs << "\neta_Ha " << eta
    << "\nprop " << prop << "\nwide " << (wide ? 1 : 0)
    << "\nmid_frac " << mid_frac << "\nwidth_frac " << width_frac
    << "\nN0 " << N0 << "\ntotal_wp_norm_tau " << total_tau
    << "\nabsorbed_fraction " << absorbed_frac
    << "\ninner_wp_norm_tau " << inner_tau << "\nepsilon " << eps
    << "\ntau " << tau << "\nN_STEPS " << N_STEPS << "\ndt " << dt << "\ndx " << dx << "\n";
  f.close();

  std::printf("  total_wp_norm(tau)=%.6f  absorbed=%.4f  inner_eps=%.6f\n",
              total_tau, absorbed_frac, eps);
  std::printf("Done -> %s/\n", outdir.c_str());
  return 0;
}
