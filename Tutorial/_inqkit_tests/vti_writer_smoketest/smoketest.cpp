// Smoke test for inqkit::io::VTIImageDataWriter.
//
// Builds a 4x4x4 RealField3D with values f(ix,iy,iz) = 100*ix + 10*iy + iz
// and a tiny 3x3x3 ComplexField3D with g(ix,iy,iz) = (ix+iy+iz, ix-iy+iz),
// then writes them to disk in both ASCII and binary VTI formats. The
// companion script `verify_smoketest.py` reads each file back and asserts
// that origin, spacing, dimensions, and per-point values are correct (and
// that ASCII vs binary agree element-wise).
//
// Build and run: `g++ -std=c++17 -I../../../inq-stack/include smoketest.cpp -o smoketest && ./smoketest`

#include <inqkit/fields/real_field_3d.hpp>
#include <inqkit/fields/complex_field_3d.hpp>
#include <inqkit/io/vti_image_data_writer.hpp>

#include <complex>
#include <iostream>

int main() {
    using inqkit::fields::RealField3D;
    using inqkit::fields::ComplexField3D;
    using inqkit::io::VTIImageDataWriter;
    using inqkit::io::VTIWriteOptions;

    // ─── Real field: 4x4x4, value(ix,iy,iz) = 100*ix + 10*iy + iz ──────────
    RealField3D rf;
    rf.nx = 4; rf.ny = 4; rf.nz = 4;
    rf.origin_x_bohr = -1.5; rf.origin_y_bohr = -1.5; rf.origin_z_bohr = -1.5;
    rf.dx_bohr = 1.0; rf.dy_bohr = 1.0; rf.dz_bohr = 1.0;
    rf.values.resize(static_cast<std::size_t>(rf.nx) * rf.ny * rf.nz);
    for (int ix = 0; ix < rf.nx; ++ix)
        for (int iy = 0; iy < rf.ny; ++iy)
            for (int iz = 0; iz < rf.nz; ++iz) {
                auto flat = (static_cast<std::size_t>(ix) * rf.ny + iy) * rf.nz + iz;
                rf.values[flat] = 100.0 * ix + 10.0 * iy + iz;
            }

    {
        VTIImageDataWriter w({.format = VTIWriteOptions::Format::ascii});
        w.write_real(rf, "out_real_ascii.vti", "scalar");
        std::cout << "wrote out_real_ascii.vti\n";
    }
    {
        VTIImageDataWriter w({.format = VTIWriteOptions::Format::binary});
        w.write_real(rf, "out_real_binary.vti", "scalar");
        std::cout << "wrote out_real_binary.vti\n";
    }

    // ─── Complex field: 3x3x3, value = (ix+iy+iz, ix-iy+iz) ────────────────
    ComplexField3D cf;
    cf.nx = 3; cf.ny = 3; cf.nz = 3;
    cf.origin_x_bohr = 0.0; cf.origin_y_bohr = 0.0; cf.origin_z_bohr = 0.0;
    cf.dx_bohr = 0.5; cf.dy_bohr = 0.5; cf.dz_bohr = 0.5;
    cf.values.resize(static_cast<std::size_t>(cf.nx) * cf.ny * cf.nz);
    for (int ix = 0; ix < cf.nx; ++ix)
        for (int iy = 0; iy < cf.ny; ++iy)
            for (int iz = 0; iz < cf.nz; ++iz) {
                auto flat = (static_cast<std::size_t>(ix) * cf.ny + iy) * cf.nz + iz;
                cf.values[flat] = std::complex<double>(
                    static_cast<double>(ix + iy + iz),
                    static_cast<double>(ix - iy + iz));
            }

    {
        VTIImageDataWriter w({.format = VTIWriteOptions::Format::ascii});
        w.write_complex(cf, "out_complex_ascii.vti", "psi");
        std::cout << "wrote out_complex_ascii.vti\n";
    }
    {
        VTIImageDataWriter w({.format = VTIWriteOptions::Format::binary});
        w.write_complex(cf, "out_complex_binary.vti", "psi");
        std::cout << "wrote out_complex_binary.vti\n";
    }

    return 0;
}
