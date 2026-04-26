#pragma once

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/complex_field_3d.hpp>
#include <inqkit/io/vti_image_data_writer.hpp>

#include <complex>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace inqkit::io {

// Same opt-in model as RealField3DLayout. emit_raw governs both
// <basename>_real.raw and <basename>_imag.raw outputs; emit_vti emits one
// <basename>.vti containing two arrays named <field_name>_real and
// <field_name>_imag.
struct ComplexField3DLayout {
  std::string field_name = "field";
  bool include_meta = true;
  bool emit_raw = true;
  bool emit_vti = false;
  VTIWriteOptions::Format vti_format = VTIWriteOptions::Format::ascii;
};

struct ComplexField3DWriteOptions {
  bool overwrite = true;
};

class ComplexField3DWriter {
public:
  ComplexField3DWriter(std::string path, ComplexField3DLayout layout = {},
                       ComplexField3DWriteOptions options = {})
      : path_(std::move(path)), layout_(std::move(layout)), options_(options) {}

  void write(inqkit::fields::ComplexField3D const &field,
             std::string const &basename) const {
    if (basename.empty()) {
      throw std::runtime_error(
          "ComplexField3DWriter: basename must not be empty.");
    }

    if (field.nx <= 0 || field.ny <= 0 || field.nz <= 0) {
      throw std::runtime_error(
          "ComplexField3DWriter: field dimensions must be positive.");
    }

    auto const expected_size =
        static_cast<std::size_t>(field.nx) * field.ny * field.nz;

    if (field.values.size() != expected_size) {
      throw std::runtime_error("ComplexField3DWriter: field.values size does "
                               "not match nx * ny * nz.");
    }

    std::filesystem::create_directories(path_);

    auto const schema =
        inqkit::detail::grid_layout::complex_field_3d_raw_schema();
    auto const stem = (std::filesystem::path(path_) / basename).string();

    if (layout_.emit_raw) {
      write_binary_parts_(stem + schema.real_suffix, stem + schema.imag_suffix,
                          field);
    }

    if (layout_.emit_raw && layout_.include_meta) {
      write_meta_file_(stem + schema.meta_suffix, field, basename, schema);
    }

    if (layout_.emit_vti) {
      VTIWriteOptions vti_opts;
      vti_opts.format = layout_.vti_format;
      vti_opts.overwrite = options_.overwrite;
      VTIImageDataWriter vti_wr(vti_opts);
      vti_wr.write_complex(field, stem + ".vti", layout_.field_name);
    }
  }

  void operator()(inqkit::fields::ComplexField3D const &field,
                  std::string const &basename) const {
    write(field, basename);
  }

private:
  void write_binary_parts_(std::string const &real_filename,
                           std::string const &imag_filename,
                           inqkit::fields::ComplexField3D const &field) const {
    auto const real_path = std::filesystem::path(real_filename);
    auto const imag_path = std::filesystem::path(imag_filename);

    if (std::filesystem::exists(real_path) && !options_.overwrite) {
      throw std::runtime_error("ComplexField3DWriter: real file already exists "
                               "and overwrite=false: " +
                               real_path.string());
    }

    if (std::filesystem::exists(imag_path) && !options_.overwrite) {
      throw std::runtime_error("ComplexField3DWriter: imag file already exists "
                               "and overwrite=false: " +
                               imag_path.string());
    }

    std::vector<double> real_values(field.values.size());
    std::vector<double> imag_values(field.values.size());

    for (std::size_t i = 0; i < field.values.size(); ++i) {
      real_values[i] = std::real(field.values[i]);
      imag_values[i] = std::imag(field.values[i]);
    }

    {
      std::ofstream out(real_path, std::ios::binary);
      if (!out) {
        throw std::runtime_error(
            "ComplexField3DWriter: could not open real file for writing: " +
            real_path.string());
      }

      out.write(
          reinterpret_cast<char const *>(real_values.data()),
          static_cast<std::streamsize>(real_values.size() * sizeof(double)));

      if (!out) {
        throw std::runtime_error(
            "ComplexField3DWriter: failed while writing real file: " +
            real_path.string());
      }
    }

    {
      std::ofstream out(imag_path, std::ios::binary);
      if (!out) {
        throw std::runtime_error(
            "ComplexField3DWriter: could not open imag file for writing: " +
            imag_path.string());
      }

      out.write(
          reinterpret_cast<char const *>(imag_values.data()),
          static_cast<std::streamsize>(imag_values.size() * sizeof(double)));

      if (!out) {
        throw std::runtime_error(
            "ComplexField3DWriter: failed while writing imag file: " +
            imag_path.string());
      }
    }
  }

  void write_meta_file_(
      std::string const &filename, inqkit::fields::ComplexField3D const &field,
      std::string const &basename,
      inqkit::detail::grid_layout::ComplexField3DRawSchema const &schema)
      const {
    auto const filepath = std::filesystem::path(filename);

    if (std::filesystem::exists(filepath) && !options_.overwrite) {
      throw std::runtime_error("ComplexField3DWriter: metadata file already "
                               "exists and overwrite=false: " +
                               filepath.string());
    }

    std::ofstream out(filepath);
    if (!out) {
      throw std::runtime_error(
          "ComplexField3DWriter: could not open metadata file: " +
          filepath.string());
    }

    out << "type = " << schema.type << "\n";
    out << "dtype = " << schema.dtype << "\n";
    out << "field_name = " << layout_.field_name << "\n";

    out << "nx = " << field.nx << "\n";
    out << "ny = " << field.ny << "\n";
    out << "nz = " << field.nz << "\n";

    out << "origin_bohr = " << field.origin_x_bohr << " " << field.origin_y_bohr
        << " " << field.origin_z_bohr << "\n";

    out << "spacing_bohr = " << field.dx_bohr << " " << field.dy_bohr << " "
        << field.dz_bohr << "\n";

    out << "layout = " << schema.layout << "\n";
    out << "real_file = " << basename << schema.real_suffix << "\n";
    out << "imag_file = " << basename << schema.imag_suffix << "\n";
  }

private:
  std::string path_;
  ComplexField3DLayout layout_;
  ComplexField3DWriteOptions options_;
};

} // namespace inqkit::io
