// ============================================================================
// localised_jellium / wp_highdensity_sv / inject_scan
//
// PURPOSE. Answer one question, cheaply and without propagating anything:
//
//   "If we launch the sigma_WP = 0.5 wavepacket only 1.5 Bohr from the jellium
//    slab face, does orthogonalisation against the occupied bath deform it?"
//
// The effective-sigma hypothesis (docs/plans/effective-sigma-near-launch.md)
// says the packet's ARRIVAL width, not its launch sigma, sets the interaction.
// Testing it means launching close to the slab — which puts the packet inside
// the electronic spill-out, where it acquires real overlap with the occupied
// manifold. inqkit's injector then projects that component out (two-pass
// modified Gram-Schmidt) and renormalises, so the packet that actually
// propagates is NOT the Gaussian we asked for. This program measures how big
// that discrepancy is.
//
// WHY NO PROPAGATION IS NEEDED. Orthonormalisation happens exactly ONCE, at
// injection (wavepacket.hpp:299-392). These runs use INQ's default ETRS
// propagator (wp/run.cpp sets none), and inq-study/src/real_time/etrs.hpp
// contains no orthogonalize() call — only crank_nicolson.hpp:139,162 does. So
// the whole risk is a t = 0 property and each trial costs a GS load, not a run.
//
// THE MEASUREMENT (user criterion 2026-08-01: removed weight < 3 %):
//   removed_weight = 1 - (||psi||_post-GS / ||psi||_raw-Gaussian)^2
// measured against the RAW GAUSSIAN's discrete norm, not against 1.0 — at
// dx = 0.40 with density std sigma/sqrt2 = 0.354 Bohr the grid barely resolves
// the packet, and that discretisation error would otherwise masquerade as
// orthogonalisation loss.
//
// A free closure check runs on every trial: because the KS states are mutually
// orthonormal, sum_i |<psi_i|psi_wp>|^2 must equal the drop in squared norm.
// The two are computed by different reductions, so agreement is real evidence.
//
// The t = 0 WP orbital is written as a complex VTI so the k_z marginal (which
// for a Gaussian is exactly N(k0, sigma_p^2), sigma_p = 1/(sqrt2 sigma)) can be
// fitted offline. Gaussianity is REPORTED, never a veto (user, 2026-08-01).
//
// USAGE
//   LJ_LAUNCH_Z=-14.0 LJ_K0=2.0 LJ_GS_DIR=<gs> ./run
// Appends one row to results/scan/scan.csv and writes a per-trial directory.
//
// ENGINE: builds against stock inq OR inq-study (no CAP here, nothing complex).
// Plan: docs/plans/effective-sigma-near-launch.md
// ============================================================================

#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/analytics.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include "../../../shared/configs/slab_n100_L35x35x85.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;
namespace obs_ = inqkit::observables;
using Cfg = localised_jellium::config::SlabN100_L35x35x85;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

// "-14.5" -> "zm14p5", for a filesystem-safe per-trial directory name.
static std::string ztag(double z){
    std::ostringstream o;
    o << "z" << (z < 0 ? "m" : "") << std::fixed << std::setprecision(1) << std::abs(z);
    std::string s = o.str();
    for (auto& c : s) if (c == '.') c = 'p';
    return s;
}

int main() {
    const double HA = 27.211386245988;

    const double SIGMA_WP = env_d("LJ_SIGMA", Cfg::WP_SIGMA_BOHR);   // 0.5
    const double SPACING  = env_d("LJ_SPACING", 0.40);               // production grid
    const double K0       = env_d("LJ_K0", 2.0);                     // = v (m = 1)
    const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", -14.0);
    const double TOL_PC   = env_d("LJ_REMOVED_TOL_PC", 3.0);         // user criterion
    const std::string GS_DIR = env_s("LJ_GS_DIR", "");
    const std::string OUT    = "results/scan/" + env_s("LJ_TRIAL", ztag(LAUNCH_Z)
                                                       + "_k" + std::to_string((int)std::lround(K0*10)));

    if (GS_DIR.empty() || !fs::exists(GS_DIR)) {
        std::cerr << "FATAL: LJ_GS_DIR missing or unset: '" << GS_DIR << "'\n";
        return 2;
    }

    const double face_z   = Cfg::SLAB_CENTER_BOHR - Cfg::SLAB_HALF_WIDTH;   // -12.5
    const double standoff = face_z - LAUNCH_Z;                             // Bohr outside the face
    const double sigma_d  = SIGMA_WP / std::sqrt(2.0);                      // density std
    const double sigma_p  = 1.0 / (std::sqrt(2.0) * SIGMA_WP);              // momentum std
    const double rs       = inqkit::jellium::rs_from_n0(Cfg::N0);
    const double kF       = std::cbrt(3.0 * M_PI * M_PI * Cfg::N0);

    fs::create_directories(OUT + "/raw/observables");
    fs::create_directories(OUT + "/raw/vti/wavefunction_wp");
    fs::create_directories(OUT + "/raw/vti/density_wp");

    std::cout << std::setprecision(10)
              << "\n=== inject_scan  OUT=" << OUT << " ===\n"
              << "  cell     = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x "
              << Cfg::LZ_BOHR << " Bohr, periodicity(2), dx=" << SPACING << "\n"
              << "  slab     = faces +/-" << Cfg::SLAB_HALF_WIDTH << ", edge "
              << Cfg::EDGE_WIDTH_BOHR << " Bohr, N=" << Cfg::N_ELECTRONS
              << ", n0=" << Cfg::N0 << ", r_s=" << rs << ", k_F=" << kF << "\n"
              << "  WP       = sigma_WP " << SIGMA_WP << " (density std " << sigma_d
              << ", sigma_p " << sigma_p << ")  k0=" << K0 << "\n"
              << "  launch_z = " << LAUNCH_Z << "  -> standoff " << standoff
              << " Bohr outside the face at " << face_z
              << " (= " << standoff/sigma_d << " density-std)\n"
              << "  criterion: removed weight < " << TOL_PC << " %\n"
              << "  GS       = " << GS_DIR << "\n\n";

    // ---- system (identical construction to wp/run.cpp) ---------------------
    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            Cfg::LZ_BOHR * 1.0_b).periodicity(2);
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

    electrons.load(GS_DIR);
    std::cout << "  Loaded GS (" << electrons.states().num_states() << " states)\n";

    // ---- inject with orthogonalisation, exactly as production does ---------
    auto report = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z)
                      .sigma(SIGMA_WP)
                      .k0(0.0, 0.0, K0)
                      .orthogonalise_against_occupied(electrons)
                      .inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;

    const double removed_pc = 100.0 * report.removed_weight;
    const bool   accept     = (removed_pc < TOL_PC);

    // ---- closure check: two independent routes to the same number ----------
    // sum_i |<psi_i|psi_wp>|^2 == ||psi||^2_pre - ||psi||^2_post exactly, since
    // the KS states are mutually orthonormal. Cheap, and it catches a wrong
    // bookkeeping change immediately.
    const double closure_abs = report.ortho_closure_residual();
    const double closure_rel = (report.sum_overlap_sq > 0.0)
                             ? closure_abs / report.sum_overlap_sq : 0.0;
    const bool   closure_ok  = (closure_rel < 1.0e-8);

    std::cout << std::setprecision(10)
              << "\n  --- orthogonalisation loss ---\n"
              << "  norm_pre_ortho    = " << report.norm_pre_ortho
              << "   (raw discrete Gaussian; deviation from 1 is grid, not loss)\n"
              << "  norm_pre_renorm   = " << report.norm_pre_renorm << "\n"
              << "  removed_weight    = " << report.removed_weight
              << "  =  " << removed_pc << " %\n"
              << "  sum_overlap_sq    = " << report.sum_overlap_sq << "\n"
              << "  max_overlap       = " << report.max_overlap << "\n"
              << "  norm_after        = " << report.norm_after
              << "   (post-renormalise; ~1 BY CONSTRUCTION, not a loss measure)\n"
              << "  closure |lhs-rhs| = " << closure_abs
              << "  (rel " << closure_rel << ")  "
              << (closure_ok ? "[OK]" : "[**CLOSURE FAILED**]") << "\n";

    if (!closure_ok) {
        std::cerr << "\nFATAL: the two routes to removed_weight disagree — the "
                     "bookkeeping is wrong, so no number here can be trusted.\n";
        return 5;
    }

    // ---- t=0 moments (real space + momentum) -------------------------------
    obs_::WPMomentumStats  wp_mom(OUT + "/raw/observables/wp_momentum_stats.csv",  wp_idx, {.write_every = 1});
    obs_::WPRealSpaceStats wp_pos(OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every = 1});
    auto m0 = wp_mom.compute(electrons);
    auto r0 = wp_pos.compute(electrons);

    const double T1 = m0.ekin;
    const double T2 = 0.5*(m0.px*m0.px + m0.py*m0.py + m0.pz*m0.pz);

    auto pc = [](double got, double want){
        return (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0; };

    std::cout << std::setprecision(8)
              << "\n  --- t=0 moments (post-orthogonalisation) ---\n"
              << "  norm (real space) = " << r0.N << "\n"
              << "  centroid z        = " << r0.zc << "   (launch " << LAUNCH_Z
              << ", dev " << (r0.zc - LAUNCH_Z) << " Bohr)\n"
              << "  density std z     = " << std::sqrt(r0.sz2)
              << "   (Gaussian " << sigma_d << ", dev " << pc(std::sqrt(r0.sz2), sigma_d) << " %)\n"
              << "  <p_z>             = " << m0.pz
              << "   (k0 " << K0 << ", dev " << pc(m0.pz, K0) << " %)\n"
              << "  sigma_pz^2        = " << m0.sz2
              << "   (1/(2 s^2) = " << sigma_p*sigma_p
              << ", dev " << pc(m0.sz2, sigma_p*sigma_p) << " %)\n"
              << "  T1 - T2           = " << (T1-T2)*HA << " eV"
              << "   (3/(4 s^2) = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA
              << " eV, dev " << pc(T1-T2, 3.0/(4.0*SIGMA_WP*SIGMA_WP)) << " %)\n";

    // ---- fields for the offline k_z / shape analysis -----------------------
    // Complex orbital: the k_z marginal is obtained from this by FFT offline.
    // VTIs are written in PHYSICAL order (fft_shift applied at write) — the
    // Python side must ifftshift BEFORE transforming. See
    // .claude/rules/vti-coordinate-mapping.md.
    inqkit::io::ComplexField3DWriter wf_wr(
        OUT + "/raw/vti/wavefunction_wp",
        {.field_name = "wavefunction", .include_meta = false, .emit_raw = false,
         .emit_vti = true,  // defaults to FALSE — omitting this writes nothing
         .vti_format = inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite = true});
    wf_wr.write(inqkit::fields::orbital::wavefunction(electrons, wp_idx), "wavefunction_t000000");

    inqkit::io::RealField3DLayout lay{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter wp_wr(OUT + "/raw/vti/density_wp", lay, {.overwrite = true});
    wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), "density_wp_t000000");

    // ---- per-state overlap spectrum ----------------------------------------
    {
        std::ofstream f(OUT + "/raw/observables/overlap_by_state.csv");
        f << std::setprecision(12) << "state,abs_overlap,overlap_sq\n";
        for (std::size_t i = 0; i < report.overlap_by_state.size(); ++i) {
            const double o = report.overlap_by_state[i];
            f << i << ',' << o << ',' << o*o << '\n';
        }
    }

    // ---- machine-readable trial record + the scan ledger -------------------
    {
        std::ofstream s(OUT + "/trial.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/wp_highdensity_sv/inject_scan/" << OUT << "\n"
          << "purpose = t=0 orthogonalisation-loss scan vs launch distance\n"
          << "plan = docs/plans/effective-sigma-near-launch.md\n"
          << "gs_dir = " << GS_DIR << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x " << Cfg::LZ_BOHR
          << "  periodicity = 2  spacing = " << SPACING << "\n"
          << "slab_face_z = " << face_z << "  edge_width = " << Cfg::EDGE_WIDTH_BOHR << "\n"
          << "wp_sigma_bohr = " << SIGMA_WP << "  wp_sigma_density = " << sigma_d
          << "  sigma_p = " << sigma_p << "\n"
          << "wp_k0 = " << K0 << "  launch_z = " << LAUNCH_Z
          << "  standoff_bohr = " << standoff << "  wp_state_index = " << wp_idx << "\n"
          << "norm_pre_ortho = " << report.norm_pre_ortho << "\n"
          << "norm_pre_renorm = " << report.norm_pre_renorm << "\n"
          << "norm_after = " << report.norm_after << "\n"
          << "removed_weight = " << report.removed_weight << "\n"
          << "removed_percent = " << removed_pc << "\n"
          << "sum_overlap_sq = " << report.sum_overlap_sq << "\n"
          << "max_overlap = " << report.max_overlap << "\n"
          << "closure_residual_abs = " << closure_abs << "\n"
          << "closure_residual_rel = " << closure_rel << "\n"
          << "centroid_z = " << r0.zc << "\n"
          << "density_std_z = " << std::sqrt(r0.sz2) << "\n"
          << "mean_pz = " << m0.pz << "\n"
          << "sigma_pz2 = " << m0.sz2 << "\n"
          << "t1_minus_t2_ev = " << (T1-T2)*HA << "\n"
          << "criterion_percent = " << TOL_PC << "\n"
          << "accept = " << (accept ? "true" : "false") << "\n"
          << "run_completed = true\n";
    }
    {
        // One append-only ledger for the whole scan, so the driver just reads
        // the last row rather than parsing per-trial files.
        const std::string ledger = "results/scan/scan.csv";
        const bool fresh = !fs::exists(ledger);
        std::ofstream f(ledger, std::ios::app);
        if (fresh)
            f << "launch_z,standoff_bohr,k0,sigma_wp,removed_weight,removed_percent,"
                 "sum_overlap_sq,max_overlap,closure_rel,centroid_z,density_std_z,"
                 "mean_pz,sigma_pz2,accept\n";
        f << std::setprecision(12)
          << LAUNCH_Z << ',' << standoff << ',' << K0 << ',' << SIGMA_WP << ','
          << report.removed_weight << ',' << removed_pc << ','
          << report.sum_overlap_sq << ',' << report.max_overlap << ','
          << closure_rel << ',' << r0.zc << ',' << std::sqrt(r0.sz2) << ','
          << m0.pz << ',' << m0.sz2 << ',' << (accept ? 1 : 0) << '\n';
    }

    std::cout << "\n  ================================================\n"
              << "   launch_z = " << LAUNCH_Z << "  (standoff " << standoff << " Bohr)\n"
              << "   removed  = " << removed_pc << " %   criterion < " << TOL_PC << " %\n"
              << "   VERDICT  : " << (accept ? "ACCEPT" : "REJECT — retreat 0.5 Bohr")
              << "\n  ================================================\n"
              << "  done -> " << OUT << "/\n";

    // Exit code IS the scan decision, so the driver loop needs no parsing:
    //   0 = accept, 3 = reject (retreat). Reserve non-zero-non-3 for errors.
    return accept ? 0 : 3;
}
