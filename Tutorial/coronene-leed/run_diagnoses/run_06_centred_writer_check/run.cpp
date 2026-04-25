// ============================================================================
// Diagnostic run 06: coronene with atoms truly centred at origin (Qball-parity)
//
// Hypothesis under test (H1, very high confidence):
//   The cross-shaped / split-into-quarters appearance of INQ density .vti
//   files is a writer-side indexing bug, NOT a SCF bug. INQ stores its
//   real-space grid in FFT-natural order (array index 0 -> physical position
//   0, i.e. cell centre). inqkit::fields::density::total/orbital writes the
//   array values straight through with origin metadata set to -L/2. The .vti
//   reader places array index 0 at the metadata origin (-L/2, the cell
//   corner), so a molecule that physically lives at the centre is rendered
//   at the corner -- with the four quadrants visible at the four corners of
//   any 2D slice.
//
// What this run does:
//   1. Uses Qball's coronene geometry (atoms truly centred on origin),
//      not the +L/2-shifted xyz used in run_01..run_05.
//   2. Adds a defensive pre-SCF assertion that all atoms lie inside
//      [-L/2, +L/2] (catches future xyz files written under the wrong
//      convention).
//   3. Runs the same SCF parameters as run_01 (PBE, cutoff 54 Ha,
//      energy_tolerance 1e-6 Ha, broyden mixing).
//   4. Writes:
//        results/density/                             total density
//        results/orbital_density/orbital_XXXX/        per-KS-state density
//        results/grid_diagnostics.txt                 INDEX/POSITION map +
//                                                     density probes along
//                                                     x, y, z. This is the
//                                                     key new instrumentation.
//        results/ground_state_summary.txt
//        results/checkpoint/
//
// Decision criteria (see plan):
//   - GS energy NEGATIVE and order -100 to -200 Ha => SCF is healthy
//   - sum(rho) * dV ~= 108 (electron count) => density is well-normalised
//   - Density probe along x at (y=0, z=0): peak at array index 0 (and at
//     ix near nx-1) => INQ uses FFT-natural ordering => writer bug confirmed
//   - .vti rendered with molecule at the four corners => writer bug confirmed
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace inq;
using namespace inq::magnitude;

static constexpr double LX_BOHR = 34.9222;
static constexpr double LY_BOHR = 34.9222;
static constexpr double LZ_BOHR = 59.9043;

static std::string orbital_tag(int ist) {
    std::ostringstream os;
    os << "orbital_" << std::setw(4) << std::setfill('0') << ist;
    return os.str();
}

int main() {
    std::cout << "\n=== diagnostic run 06: coronene centred at origin (writer-bug test) ===\n";

    auto cell = systems::cell::orthorhombic(
        LX_BOHR * 1.0_b, LY_BOHR * 1.0_b, LZ_BOHR * 1.0_b
    ).finite();

    auto ions = systems::ions::parse("coronene_centred.xyz", cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

    // Defensive assertion: every atom must lie inside [-L/2, +L/2] for INQ
    // to treat it without periodic wrap-around. If this fires, the xyz file
    // is using the wrong cell-origin convention (likely +L/2-shifted).
    {
        const double half_lx = 0.5 * LX_BOHR;
        const double half_ly = 0.5 * LY_BOHR;
        const double half_lz = 0.5 * LZ_BOHR;
        double max_excess = 0.0;
        int max_excess_iatom = -1;
        for (int iatom = 0; iatom < static_cast<int>(ions.size()); ++iatom) {
            auto const & p = ions.positions()[iatom];
            double ex = std::max({0.0,
                                  std::fabs(p[0]) - half_lx,
                                  std::fabs(p[1]) - half_ly,
                                  std::fabs(p[2]) - half_lz});
            if (ex > max_excess) { max_excess = ex; max_excess_iatom = iatom; }
        }
        if (max_excess > 1e-9) {
            std::cerr << "FATAL: atom " << max_excess_iatom
                      << " is outside [-L/2, +L/2] by " << max_excess
                      << " Bohr. INQ uses centred-cell convention; the xyz file\n"
                      << "must place atoms in [-L/2, +L/2]. Aborting.\n";
            return 2;
        }
        std::cout << "  All atoms are inside [-L/2, +L/2]. Max boundary slack: "
                  << (half_lx - 0.0) << " Bohr (cell half-extent x).\n";
    }

    auto electrons = systems::electrons(
        ions,
        options::electrons{}.cutoff(54.0_Ha).extra_states(8)
    );

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions,
        electrons,
        options::theory{}.pbe(),
        options::ground_state{}
            .energy_tolerance(1e-6_Ha)
            .max_steps(1000)
            .broyden_mixing()
            .mixing_ndim(8)
            .mixing(0.1)
    );

    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    std::filesystem::create_directories("results/density");
    std::filesystem::create_directories("results/orbital_density");
    std::filesystem::create_directories("results/checkpoint");

    // ---- Total density and orbital densities (existing pipeline) -----------
    {
        inqkit::io::RealField3DWriter density_writer(
            "results/density",
            {.field_name = "density", .include_meta = true},
            {.overwrite = true}
        );
        auto rho_total = inqkit::fields::density::total(electrons);
        density_writer.write(rho_total, 0.0, 0);
    }

    const int nstates = electrons.states().num_states();
    std::cout << "  Writing densities for " << nstates << " KS states\n";

    {
        std::ofstream map("results/orbital_density/orbital_index_map.csv");
        map << "state_index,directory\n";
        for (int ist = 0; ist < nstates; ++ist)
            map << ist << ',' << orbital_tag(ist) << "\n";
    }

    for (int ist = 0; ist < nstates; ++ist) {
        const auto out_dir = std::string("results/orbital_density/") + orbital_tag(ist);
        std::filesystem::create_directories(out_dir);

        inqkit::io::RealField3DWriter orbital_density_writer(
            out_dir,
            {.field_name = "density", .include_meta = true},
            {.overwrite = true}
        );
        auto rho_orbital = inqkit::fields::density::orbital(electrons, ist);
        orbital_density_writer.write(rho_orbital, 0.0, 0);
    }

    // ---- Grid diagnostics: the key new instrumentation --------------------
    // We pull the basis directly from electrons.density() and read hc[ix][iy][iz]
    // to dump the raw INQ index ordering, before any inqkit transformation.
    if (electrons.root()) {
        auto density = electrons.density();
        auto const & basis = density.basis();
        if (basis.comm().size() != 1) {
            std::cerr << "WARNING: multi-rank basis; grid_diagnostics.txt may be incomplete.\n";
        }
        const int nx = basis.sizes()[0];
        const int ny = basis.sizes()[1];
        const int nz = basis.sizes()[2];
        const auto spacing = basis.rspacing();
        const double dx = spacing[0];
        const double dy = spacing[1];
        const double dz = spacing[2];

        auto hc = density.cubic();

        auto fft_natural_pos = [](int idx, int size, double sp) -> double {
            int sym = (idx >= (size + 1) / 2) ? idx - size : idx;
            return sym * sp;
        };
        auto naive_pos = [](int idx, int size, double sp) -> double {
            return -0.5 * size * sp + idx * sp;
        };

        // Find the array index whose physical position (under FFT-natural
        // mapping) is closest to 0 along each axis.
        auto closest_to_zero = [](int size) {
            return 0; // by FFT-natural, idx=0 IS position 0
        };
        const int ix0 = closest_to_zero(nx);
        const int iy0 = closest_to_zero(ny);
        const int iz0 = closest_to_zero(nz);

        // Also locate the array index of the physical "cell corner" (-L/2)
        // under FFT-natural -> this is array index nx/2 (for even nx).
        const int ix_corner = nx / 2;
        const int iy_corner = ny / 2;
        const int iz_corner = nz / 2;

        std::ofstream g("results/grid_diagnostics.txt");
        g << std::setprecision(12);
        g << "# Grid layout diagnostics for run_06 centred coronene\n";
        g << "# Hypothesis: hc[0][0][0] = density at cell CENTRE (FFT-natural),\n";
        g << "# but inqkit writer claims origin = -L/2, so the .vti renders\n";
        g << "# the molecule at the corners.\n";
        g << "\n";
        g << "[grid]\n";
        g << "nx = " << nx << "\n";
        g << "ny = " << ny << "\n";
        g << "nz = " << nz << "\n";
        g << "dx_bohr = " << dx << "\n";
        g << "dy_bohr = " << dy << "\n";
        g << "dz_bohr = " << dz << "\n";
        g << "Lx_bohr = " << nx * dx << "\n";
        g << "Ly_bohr = " << ny * dy << "\n";
        g << "Lz_bohr = " << nz * dz << "\n";
        g << "symmetric_range_begin_x = " << basis.symmetric_range_begin(0) << "\n";
        g << "symmetric_range_begin_y = " << basis.symmetric_range_begin(1) << "\n";
        g << "symmetric_range_begin_z = " << basis.symmetric_range_begin(2) << "\n";
        g << "symmetric_range_end_x   = " << basis.symmetric_range_end(0) << "\n";
        g << "symmetric_range_end_y   = " << basis.symmetric_range_end(1) << "\n";
        g << "symmetric_range_end_z   = " << basis.symmetric_range_end(2) << "\n";

        g << "\n[corner_positions]\n";
        g << "# Physical position (Bohr) of the array index 0 under both mapping conventions.\n";
        g << "# x_fft_natural[0]   = " << fft_natural_pos(0, nx, dx) << "  (= 0, cell centre)\n";
        g << "# x_naive_origin[0]  = " << naive_pos(0, nx, dx) << "  (= -Lx/2, cell -corner)\n";
        g << "# y_fft_natural[0]   = " << fft_natural_pos(0, ny, dy) << "\n";
        g << "# y_naive_origin[0]  = " << naive_pos(0, ny, dy) << "\n";
        g << "# z_fft_natural[0]   = " << fft_natural_pos(0, nz, dz) << "\n";
        g << "# z_naive_origin[0]  = " << naive_pos(0, nz, dz) << "\n";

        g << "\n[hot_cells]\n";
        g << "# Density at array indices that should be 'centre' or 'corner' under each convention.\n";
        g << "rho[0][0][0]                     = " << hc[0][0][0]
          << "  (FFT-natural: cell centre)\n";
        g << "rho[nx/2][ny/2][nz/2]            = " << hc[ix_corner][iy_corner][iz_corner]
          << "  (FFT-natural: -L/2 corner)\n";
        g << "rho[0][0][nz/2]                  = " << hc[0][0][iz_corner]
          << "  (FFT-natural: x=0,y=0,z=-Lz/2)\n";
        g << "rho[nx-1][0][0]                  = " << hc[nx-1][0][0]
          << "  (FFT-natural: x=-dx, y=0, z=0)\n";
        g << "rho[1][0][0]                     = " << hc[1][0][0]
          << "  (FFT-natural: x=+dx, y=0, z=0)\n";

        // ---- 1D probes -----------------------------------------------------
        g << "\n[probe_x_at_y0_z0]\n";
        g << "# columns: array_ix, rho, x_fft_natural[Bohr], x_naive[Bohr]\n";
        for (int ix = 0; ix < nx; ++ix) {
            g << ix << "\t" << hc[ix][iy0][iz0]
              << "\t" << fft_natural_pos(ix, nx, dx)
              << "\t" << naive_pos(ix, nx, dx) << "\n";
        }
        g << "\n[probe_y_at_x0_z0]\n";
        g << "# columns: array_iy, rho, y_fft_natural[Bohr], y_naive[Bohr]\n";
        for (int iy = 0; iy < ny; ++iy) {
            g << iy << "\t" << hc[ix0][iy][iz0]
              << "\t" << fft_natural_pos(iy, ny, dy)
              << "\t" << naive_pos(iy, ny, dy) << "\n";
        }
        g << "\n[probe_z_at_x0_y0]\n";
        g << "# columns: array_iz, rho, z_fft_natural[Bohr], z_naive[Bohr]\n";
        for (int iz = 0; iz < nz; ++iz) {
            g << iz << "\t" << hc[ix0][iy0][iz]
              << "\t" << fft_natural_pos(iz, nz, dz)
              << "\t" << naive_pos(iz, nz, dz) << "\n";
        }

        // ---- Integrated density (independent of indexing) -----------------
        long double total_rho = 0.0L;
        for (int ix = 0; ix < nx; ++ix)
            for (int iy = 0; iy < ny; ++iy)
                for (int iz = 0; iz < nz; ++iz)
                    total_rho += static_cast<long double>(hc[ix][iy][iz]);
        long double integrated = total_rho * dx * dy * dz;
        g << "\n[integrated]\n";
        g << "sum_rho_times_dV = " << static_cast<double>(integrated)
          << "  (should equal num_electrons = "
          << electrons.states().num_electrons() << ")\n";

        std::cout << "  grid_diagnostics.txt written. Integrated rho = "
                  << static_cast<double>(integrated)
                  << " (target = " << electrons.states().num_electrons() << ")\n";
    }

    electrons.save("results/checkpoint");

    if (electrons.root()) {
        std::ofstream summary("results/ground_state_summary.txt");
        summary << std::setprecision(16);
        summary << "run = diagnostic_06_centred_writer_check\n";
        summary << "system = coronene_C24H12\n";
        summary << "geometry_file = coronene_centred.xyz\n";
        summary << "atom_convention = origin-centred (Qball-parity)\n";
        summary << "cell_bohr = " << LX_BOHR << ' ' << LY_BOHR << ' ' << LZ_BOHR << "\n";
        summary << "boundary = finite\n";
        summary << "xc = pbe\n";
        summary << "cutoff_ha = 54.0\n";
        summary << "extra_states = 8\n";
        summary << "energy_tolerance_ha = 1e-6\n";
        summary << "max_steps = 1000\n";
        summary << "mixing = broyden\n";
        summary << "mixing_ndim = 8\n";
        summary << "mixing_alpha = 0.1\n";
        summary << "num_atoms = " << ions.size() << "\n";
        summary << "num_electrons = " << electrons.states().num_electrons() << "\n";
        summary << "num_states = " << nstates << "\n";
        summary << "ground_state_energy_ha = " << gs.energy.total() << "\n";
    }

    std::cout << "Done. Output written to results/\n";
    return 0;
}
