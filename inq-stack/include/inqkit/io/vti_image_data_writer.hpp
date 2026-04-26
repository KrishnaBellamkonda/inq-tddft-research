// ============================================================================
// inqkit::io::VTIImageDataWriter
//
// Header-only writer that emits a single VTK XML ImageData (`.vti`) file from
// an inqkit RealField3D or ComplexField3D. Sits next to (or replaces) the
// `.raw` / `.meta.txt` outputs produced by RealField3DWriter and
// ComplexField3DWriter so that runs can publish ParaView-ready data without a
// post-processing Python step.
//
// Format support:
//   - ASCII (Float64, x-fastest point order, deterministic and diff-friendly)
//   - Inline binary (base64-encoded VTK RawBinary block; ~3x smaller, still
//     a single file). Uses a UInt64 little-endian byte-count header inside
//     the base64 payload, matching the `header_type="UInt64"` declaration on
//     <VTKFile>. See VTK XML format spec.
//
// Layout transform:
//   inqkit fields store values as flat = ((ix*ny)+iy)*nz + iz (x slowest,
//   z fastest). VTK ImageData expects PointData in x-fastest order, i.e.
//   point index = ix + nx * (iy + ny * iz). The writer reorders during the
//   ASCII stream / binary copy by walking iz outermost, iy middle, ix
//   innermost. No transpose buffer is required.
//
// Single-rank only — same restriction as the existing inqkit writers.
// ============================================================================
#pragma once

#include <inqkit/fields/complex_field_3d.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace inqkit::io {

struct VTIWriteOptions {
  // ascii: human-readable, ~3-5x larger.
  // binary: base64-encoded RawBinary inline, smaller, still single-file.
  enum class Format { ascii, binary };
  Format format = Format::ascii;

  // If false the writer throws if the target file already exists.
  bool overwrite = true;
};

class VTIImageDataWriter {
public:
  explicit VTIImageDataWriter(VTIWriteOptions options = {})
      : options_(std::move(options)) {}

  // Write a scalar grid as a single Float64 PointData array named
  // `array_name` (default "density").
  void write_real(inqkit::fields::RealField3D const &field,
                  std::string const &output_path,
                  std::string const &array_name = "density") const {
    validate_dims_real_(field);
    open_or_throw_(output_path);

    std::ofstream out(output_path, std::ios::binary);
    if (!out) {
      throw std::runtime_error(
          "VTIImageDataWriter: could not open file for writing: " +
          output_path);
    }

    write_prologue_(out, field.nx, field.ny, field.nz,
                    field.origin_x_bohr, field.origin_y_bohr, field.origin_z_bohr,
                    field.dx_bohr, field.dy_bohr, field.dz_bohr,
                    array_name);

    write_real_data_array_(out, field, array_name);

    write_epilogue_(out);

    if (!out) {
      throw std::runtime_error(
          "VTIImageDataWriter: failed while writing file: " + output_path);
    }
  }

  // Write a complex grid as two Float64 PointData arrays:
  //   <array_name>_real and <array_name>_imag
  void write_complex(inqkit::fields::ComplexField3D const &field,
                     std::string const &output_path,
                     std::string const &array_name = "psi") const {
    validate_dims_complex_(field);
    open_or_throw_(output_path);

    std::ofstream out(output_path, std::ios::binary);
    if (!out) {
      throw std::runtime_error(
          "VTIImageDataWriter: could not open file for writing: " +
          output_path);
    }

    auto const real_name = array_name + "_real";
    write_prologue_(out, field.nx, field.ny, field.nz,
                    field.origin_x_bohr, field.origin_y_bohr, field.origin_z_bohr,
                    field.dx_bohr, field.dy_bohr, field.dz_bohr,
                    real_name);

    write_complex_part_array_(out, field, array_name + "_real",
                              /*take_real=*/true);
    write_complex_part_array_(out, field, array_name + "_imag",
                              /*take_real=*/false);

    write_epilogue_(out);

    if (!out) {
      throw std::runtime_error(
          "VTIImageDataWriter: failed while writing file: " + output_path);
    }
  }

private:
  // ───────────────── Validation / fs helpers ─────────────────────────────
  static void validate_dims_real_(inqkit::fields::RealField3D const &f) {
    if (f.nx <= 0 || f.ny <= 0 || f.nz <= 0) {
      throw std::runtime_error(
          "VTIImageDataWriter: field dimensions must be positive.");
    }
    auto const expected =
        static_cast<std::size_t>(f.nx) * f.ny * f.nz;
    if (f.values.size() != expected) {
      throw std::runtime_error(
          "VTIImageDataWriter: real field values size does not match nx*ny*nz.");
    }
  }

  static void validate_dims_complex_(inqkit::fields::ComplexField3D const &f) {
    if (f.nx <= 0 || f.ny <= 0 || f.nz <= 0) {
      throw std::runtime_error(
          "VTIImageDataWriter: field dimensions must be positive.");
    }
    auto const expected =
        static_cast<std::size_t>(f.nx) * f.ny * f.nz;
    if (f.values.size() != expected) {
      throw std::runtime_error("VTIImageDataWriter: complex field values size "
                               "does not match nx*ny*nz.");
    }
  }

  void open_or_throw_(std::string const &path) const {
    auto const p = std::filesystem::path(path);
    if (p.has_parent_path()) {
      std::filesystem::create_directories(p.parent_path());
    }
    if (std::filesystem::exists(p) && !options_.overwrite) {
      throw std::runtime_error(
          "VTIImageDataWriter: file already exists and overwrite=false: " +
          path);
    }
  }

  // ───────────────── XML scaffolding ─────────────────────────────────────
  void write_prologue_(std::ostream &out,
                       int nx, int ny, int nz,
                       double ox, double oy, double oz,
                       double dx, double dy, double dz,
                       std::string const &scalars_name) const {
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"1.0\" "
           "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out << "  <ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1)
        << " 0 " << (nz - 1) << "\""
        << " Origin=\"" << format_double_(ox) << " " << format_double_(oy)
        << " " << format_double_(oz) << "\""
        << " Spacing=\"" << format_double_(dx) << " " << format_double_(dy)
        << " " << format_double_(dz) << "\">\n";
    out << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 "
        << (nz - 1) << "\">\n";
    out << "      <PointData Scalars=\"" << escape_xml_(scalars_name) << "\">\n";
  }

  void write_epilogue_(std::ostream &out) const {
    out << "      </PointData>\n";
    out << "      <CellData>\n";
    out << "      </CellData>\n";
    out << "    </Piece>\n";
    out << "  </ImageData>\n";
    out << "</VTKFile>\n";
  }

  // ───────────────── Real DataArray ──────────────────────────────────────
  void write_real_data_array_(std::ostream &out,
                              inqkit::fields::RealField3D const &f,
                              std::string const &array_name) const {
    if (options_.format == VTIWriteOptions::Format::ascii) {
      out << "        <DataArray type=\"Float64\" Name=\""
          << escape_xml_(array_name) << "\" format=\"ascii\">\n";
      stream_ascii_real_(out, f);
      out << "        </DataArray>\n";
    } else {
      // Inline RawBinary, base64-encoded.
      out << "        <DataArray type=\"Float64\" Name=\""
          << escape_xml_(array_name) << "\" format=\"binary\">\n";
      stream_binary_real_(out, f);
      out << "\n        </DataArray>\n";
    }
  }

  // x-fastest stream order: walk iz outermost, iy middle, ix innermost.
  // Reads our storage f.values[((ix*ny)+iy)*nz + iz].
  void stream_ascii_real_(std::ostream &out,
                          inqkit::fields::RealField3D const &f) const {
    auto const ny = f.ny, nz = f.nz;
    constexpr int values_per_line = 8;
    int on_line = 0;
    out << std::setprecision(17);
    out << "          ";
    for (int iz = 0; iz < f.nz; ++iz) {
      for (int iy = 0; iy < f.ny; ++iy) {
        for (int ix = 0; ix < f.nx; ++ix) {
          auto const flat =
              (static_cast<std::size_t>(ix) * ny + iy) * nz + iz;
          if (on_line == values_per_line) {
            out << "\n          ";
            on_line = 0;
          } else if (on_line > 0) {
            out << ' ';
          }
          out << f.values[flat];
          ++on_line;
        }
      }
    }
    out << '\n';
  }

  void stream_binary_real_(std::ostream &out,
                           inqkit::fields::RealField3D const &f) const {
    std::vector<double> buf;
    buf.resize(f.values.size());
    auto const ny = f.ny, nz = f.nz;
    std::size_t k = 0;
    for (int iz = 0; iz < f.nz; ++iz) {
      for (int iy = 0; iy < f.ny; ++iy) {
        for (int ix = 0; ix < f.nx; ++ix) {
          auto const flat =
              (static_cast<std::size_t>(ix) * ny + iy) * nz + iz;
          buf[k++] = f.values[flat];
        }
      }
    }
    write_base64_rawbinary_(
        out,
        reinterpret_cast<unsigned char const *>(buf.data()),
        buf.size() * sizeof(double));
  }

  // ───────────────── Complex DataArray (one component per call) ──────────
  void write_complex_part_array_(std::ostream &out,
                                 inqkit::fields::ComplexField3D const &f,
                                 std::string const &array_name,
                                 bool take_real) const {
    if (options_.format == VTIWriteOptions::Format::ascii) {
      out << "        <DataArray type=\"Float64\" Name=\""
          << escape_xml_(array_name) << "\" format=\"ascii\">\n";
      auto const ny = f.ny, nz = f.nz;
      constexpr int values_per_line = 8;
      int on_line = 0;
      out << std::setprecision(17);
      out << "          ";
      for (int iz = 0; iz < f.nz; ++iz) {
        for (int iy = 0; iy < f.ny; ++iy) {
          for (int ix = 0; ix < f.nx; ++ix) {
            auto const flat =
                (static_cast<std::size_t>(ix) * ny + iy) * nz + iz;
            if (on_line == values_per_line) {
              out << "\n          ";
              on_line = 0;
            } else if (on_line > 0) {
              out << ' ';
            }
            out << (take_real ? f.values[flat].real() : f.values[flat].imag());
            ++on_line;
          }
        }
      }
      out << "\n        </DataArray>\n";
    } else {
      out << "        <DataArray type=\"Float64\" Name=\""
          << escape_xml_(array_name) << "\" format=\"binary\">\n";
      std::vector<double> buf;
      buf.resize(f.values.size());
      auto const ny = f.ny, nz = f.nz;
      std::size_t k = 0;
      for (int iz = 0; iz < f.nz; ++iz) {
        for (int iy = 0; iy < f.ny; ++iy) {
          for (int ix = 0; ix < f.nx; ++ix) {
            auto const flat =
                (static_cast<std::size_t>(ix) * ny + iy) * nz + iz;
            buf[k++] =
                take_real ? f.values[flat].real() : f.values[flat].imag();
          }
        }
      }
      write_base64_rawbinary_(
          out,
          reinterpret_cast<unsigned char const *>(buf.data()),
          buf.size() * sizeof(double));
      out << "\n        </DataArray>\n";
    }
  }

  // ───────────────── Inline RawBinary helper ─────────────────────────────
  // VTK RawBinary inline payload = uint64 little-endian header (#bytes that
  // follow) + raw bytes, then the whole thing base64-encoded as one chunk.
  static void write_base64_rawbinary_(std::ostream &out,
                                      unsigned char const *payload,
                                      std::size_t payload_bytes) {
    std::vector<unsigned char> blob;
    blob.resize(8 + payload_bytes);

    // Little-endian uint64 header. Code is byte-order agnostic.
    std::uint64_t hdr = static_cast<std::uint64_t>(payload_bytes);
    for (int i = 0; i < 8; ++i) {
      blob[i] = static_cast<unsigned char>((hdr >> (8 * i)) & 0xFF);
    }
    std::memcpy(blob.data() + 8, payload, payload_bytes);

    out << "          ";
    base64_stream_(out, blob.data(), blob.size());
  }

  static void base64_stream_(std::ostream &out, unsigned char const *data,
                             std::size_t n) {
    static char const *tbl =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    std::size_t i = 0;
    while (i + 3 <= n) {
      std::uint32_t v = (static_cast<std::uint32_t>(data[i]) << 16) |
                        (static_cast<std::uint32_t>(data[i + 1]) << 8) |
                        (static_cast<std::uint32_t>(data[i + 2]));
      out.put(tbl[(v >> 18) & 0x3F]);
      out.put(tbl[(v >> 12) & 0x3F]);
      out.put(tbl[(v >> 6) & 0x3F]);
      out.put(tbl[v & 0x3F]);
      i += 3;
    }
    std::size_t rem = n - i;
    if (rem == 1) {
      std::uint32_t v = static_cast<std::uint32_t>(data[i]) << 16;
      out.put(tbl[(v >> 18) & 0x3F]);
      out.put(tbl[(v >> 12) & 0x3F]);
      out.put('=');
      out.put('=');
    } else if (rem == 2) {
      std::uint32_t v = (static_cast<std::uint32_t>(data[i]) << 16) |
                        (static_cast<std::uint32_t>(data[i + 1]) << 8);
      out.put(tbl[(v >> 18) & 0x3F]);
      out.put(tbl[(v >> 12) & 0x3F]);
      out.put(tbl[(v >> 6) & 0x3F]);
      out.put('=');
    }
  }

  // ───────────────── small utilities ─────────────────────────────────────
  static std::string format_double_(double v) {
    std::ostringstream ss;
    ss << std::setprecision(17) << v;
    return ss.str();
  }

  static std::string escape_xml_(std::string const &s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
      switch (c) {
      case '&': out += "&amp;"; break;
      case '<': out += "&lt;"; break;
      case '>': out += "&gt;"; break;
      case '"': out += "&quot;"; break;
      case '\'': out += "&apos;"; break;
      default: out.push_back(c);
      }
    }
    return out;
  }

private:
  VTIWriteOptions options_{};
};

} // namespace inqkit::io
