// ============================================================================
// localised_jellium / classical_highdensity_sv / vac_exit / run.cpp
//
// Phase 1 (vacuum exit test) of the classical-highdensity-sv campaign.
//
// GOAL: prove that a moving Gaussian CHARGE, built by the REAL inqkit code path
// (inqkit::jellium::gaussian_density on a z-OPEN periodicity(2) basis), is
// CLIPPED by the finite z-grid as its center leaves the box (only the in-box
// portion survives) and does NOT wrap around to the opposite face.
//
// Pure density-construction test: EMPTY cell (no jellium background, no
// projectile potential added to any Hamiltonian, no propagation). For each
// z_center on a sweep crossing the far face +Lz/2 = +42.5, we build n_proj with
// the SAME free function the real moving-projectile perturbation uses
// (moving_gaussian_projectile_perturbation.hpp:48), write a physical-order VTI,
// and log integral(n_proj) + the near-face (z<-38) max density (the "wrap
// witness", must stay ~0).
//
// Cell:  35 x 35 x 85 Bohr orthorhombic, periodicity(2) (z OPEN), dx = 0.5.
// sigma_pot = 0.35355 (= sigma_WP/sqrt(2) for sigma_WP = 0.5).
//
// No inq/ or inq-study/ edit — wrapper-only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/real_field_3d.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace inq;
using namespace inq::magnitude;

// Convert an INQ real-space field (FFT-natural array order) into an inqkit
// RealField3D in PHYSICAL order, mirroring inqkit::fields::density::total
// (density.hpp:56-116). Same fft_shift into contiguous -L/2 .. +L/2 layout, so
// the emitted VTI is physical-order and loads via inqview.load_vti.
template <class Field>
inqkit::fields::RealField3D field_to_physical(Field const & f) {
	auto const & basis = f.basis();
	if (basis.comm().size() != 1) {
		throw std::runtime_error("field_to_physical: multi-rank basis not supported.");
	}
	auto const nx = basis.sizes()[0];
	auto const ny = basis.sizes()[1];
	auto const nz = basis.sizes()[2];
	auto const spacing = basis.rspacing();

	inqkit::fields::RealField3D field;
	field.nx = nx; field.ny = ny; field.nz = nz;
	field.origin_x_bohr = basis.symmetric_range_begin(0) * spacing[0];
	field.origin_y_bohr = basis.symmetric_range_begin(1) * spacing[1];
	field.origin_z_bohr = basis.symmetric_range_begin(2) * spacing[2];
	field.dx_bohr = spacing[0];
	field.dy_bohr = spacing[1];
	field.dz_bohr = spacing[2];
	field.values.resize(static_cast<std::size_t>(nx) * ny * nz);

	boost::multi::array<double, 3> host_hc{f.cubic()};
	for (int ix = 0; ix < nx; ++ix) {
		int sx = inqkit::detail::grid_layout::fft_shift_index(ix, nx);
		for (int iy = 0; iy < ny; ++iy) {
			int sy = inqkit::detail::grid_layout::fft_shift_index(iy, ny);
			for (int iz = 0; iz < nz; ++iz) {
				int sz = inqkit::detail::grid_layout::fft_shift_index(iz, nz);
				auto flat = inqkit::detail::grid_layout::flatten_index(ix, iy, iz, ny, nz);
				field.values[flat] = host_hc[sx][sy][sz];
			}
		}
	}
	return field;
}

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }

int main() {
	const double LX = 35.0, LY = 35.0, LZ = 85.0;   // Bohr
	const double SPACING   = env_d("VE_SPACING", 0.5);
	const double SIGMA_POT = env_d("VE_SIGMA_POT", 0.35355);   // = 0.5/sqrt(2)
	const double FAR_FACE  = LZ / 2.0;                          // +42.5
	const double NEAR_FACE_ZMAX = -38.0;   // "wrap witness" region z < -38

	std::cout << "\n=== vac_exit: EMPTY " << LX << "x" << LY << "x" << LZ
	          << " periodicity(2)  dx=" << SPACING
	          << "  sigma_pot=" << SIGMA_POT << " ===\n"
	          << "  far face (+Lz/2) = " << FAR_FACE
	          << "   wrap-witness region z < " << NEAR_FACE_ZMAX << "\n\n";

	// z-open cell + a real-space basis (via an electrons object). A small nonzero
	// extra_electrons count makes INQ build a normal populated density/basis (the
	// grid geometry is set by cell+spacing, independent of electron count — the
	// gaussian_density path only reads the basis). No background/projectile is
	// added to any Hamiltonian; there is no SCF and no propagation.
	auto cell = systems::cell::orthorhombic(LX * 1.0_b, LY * 1.0_b, LZ * 1.0_b).periodicity(2);
	auto ions = systems::ions(cell);
	auto electrons = systems::electrons(
		ions,
		options::electrons{}.spacing(SPACING * 1.0_b).extra_electrons(2),
		input::kpoints::gamma());

	auto basis = electrons.density().basis();
	std::cout << "  basis local_sizes = ("
	          << basis.local_sizes()[0] << ", "
	          << basis.local_sizes()[1] << ", "
	          << basis.local_sizes()[2] << ")\n"
	          << "  basis sizes       = ("
	          << basis.sizes()[0] << ", " << basis.sizes()[1] << ", "
	          << basis.sizes()[2] << ")\n" << std::flush;

	// Full transit sweep for the along-z animation (enter -> cross slab region ->
	// exit +42.5), with a fine crossing of the far face.
	std::vector<double> z_centers;
	for (double z = -44.0; z <= 50.0 + 1e-9; z += 2.0) z_centers.push_back(z);
	for (double z = 41.0; z <= 44.0 + 1e-9; z += 0.5) z_centers.push_back(z); // fine face
	std::sort(z_centers.begin(), z_centers.end());
	z_centers.erase(std::unique(z_centers.begin(), z_centers.end(),
		[](double a, double b){ return std::abs(a - b) < 1e-6; }), z_centers.end());

	std::filesystem::create_directories("results");
	std::ofstream csv("results/exit_scan.csv");
	csv << "z_center,integral,wrap_witness_max,phi_peak\n";
	csv << std::setprecision(12);

	for (double zc : z_centers) {
		inq::vector3<double> center{0.0, 0.0, zc};
		auto nproj = inqkit::jellium::gaussian_density(basis, center, SIGMA_POT);

		const double integ = inq::operations::integral(nproj);

		// Physical-order host copy for the VTI + wrap-witness scan.
		auto phys = field_to_physical(nproj);

		// wrap witness: max density over the near face z < NEAR_FACE_ZMAX.
		double wrap_max = 0.0;
		for (int ix = 0; ix < phys.nx; ++ix)
			for (int iy = 0; iy < phys.ny; ++iy)
				for (int iz = 0; iz < phys.nz; ++iz) {
					double zcoord = phys.origin_z_bohr + iz * phys.dz_bohr;
					if (zcoord < NEAR_FACE_ZMAX) {
						auto flat = inqkit::detail::grid_layout::flatten_index(ix, iy, iz, phys.ny, phys.nz);
						wrap_max = std::max(wrap_max, phys.values[flat]);
					}
				}

		// Perturbation POTENTIAL the bath actually feels: φ_proj = poisson(n_proj)
		// (exactly moving_gaussian_projectile_perturbation.hpp:49; a −1 projectile
		// adds +φ ⇒ a repulsive bump). z-open ⇒ slab-truncated Poisson.
		auto phi = inq::solvers::poisson::solve(nproj);
		auto phi_phys = field_to_physical(phi);
		double phi_peak = 0.0;
		for (double v : phi_phys.values) phi_peak = std::max(phi_peak, v);

		csv << zc << "," << integ << "," << wrap_max << "," << phi_peak << "\n";
		std::cout << "  z_center=" << std::setw(7) << zc
		          << "  integral=" << std::setw(12) << integ
		          << "  wrap_witness_max=" << wrap_max
		          << "  phi_peak=" << phi_peak << "\n";

		// VTIs (physical order, binary) named by z_center (n_proj + φ_proj).
		char zbuf[32];
		std::snprintf(zbuf, sizeof(zbuf), "nproj_z%+06.1f", zc);
		for (char* p = zbuf; *p; ++p) if (*p == '.') *p = 'p';
		inqkit::io::RealField3DWriter wr("results",
			{ .field_name = "nproj", .include_meta = false, .emit_raw = false,
			  .emit_vti = true,
			  .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
			{ .overwrite = true });
		wr.write(phys, std::string(zbuf));

		char pbuf[32];
		std::snprintf(pbuf, sizeof(pbuf), "phi_z%+06.1f", zc);
		for (char* p = pbuf; *p; ++p) if (*p == '.') *p = 'p';
		inqkit::io::RealField3DWriter wrp("results",
			{ .field_name = "phi_proj", .include_meta = false, .emit_raw = false,
			  .emit_vti = true,
			  .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
			{ .overwrite = true });
		wrp.write(phi_phys, std::string(pbuf));
	}

	{
		std::ofstream s("results/run_summary.txt");
		s << std::setprecision(14)
		  << "run = localised_jellium/classical_highdensity_sv/vac_exit\n"
		  << "engine = inq\n"
		  << "test = vacuum_exit_clipping\n"
		  << "cell_bohr = 35x35x85 (orthorhombic)\n"
		  << "periodicity = 2\n"
		  << "spacing_bohr = " << SPACING << "\n"
		  << "sigma_pot_bohr = " << SIGMA_POT << "\n"
		  << "far_face_bohr = " << FAR_FACE << "\n"
		  << "wrap_witness_zmax_bohr = " << NEAR_FACE_ZMAX << "\n"
		  << "n_z_centers = " << z_centers.size() << "\n"
		  << "scan_csv = results/exit_scan.csv\n"
		  << "run_completed = true\n";
	}
	std::cout << "\nDone.\n";
	return 0;
}
