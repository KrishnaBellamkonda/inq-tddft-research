#pragma once

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace inqkit::io {

struct RealField3DLayout {
  std::string field_name = "field";
  bool include_meta = true;
};

struct RealField3DWriteOptions {
  bool overwrite = true;
};

class RealField3DWriter {
public:
  RealField3DWriter(std::string path, RealField3DLayout layout = {},
                    RealField3DWriteOptions options = {})
      : path_(std::move(path)), layout_(std::move(layout)), options_(options) {}

  void write(inqkit::fields::RealField3D const &field,
             std::string const &basename) const {
    if (basename.empty()) {
      throw std::runtime_error(
          "RealField3DWriter: basename must not be empty.");
    }

    if (field.nx <= 0 || field.ny <= 0 || field.nz <= 0) {
      throw std::runtime_error(
          "RealField3DWriter: field dimensions must be positive.");
    }

    auto const expected_size =
        static_cast<std::size_t>(field.nx) * field.ny * field.nz;

    if (field.values.size() != expected_size) {
      throw std::runtime_error(
          "RealField3DWriter: field.values size does not match nx * ny * nz.");
    }

    std::filesystem::create_directories(path_);

    auto const schema = inqkit::detail::grid_layout::real_field_3d_raw_schema();
    auto const stem = (std::filesystem::path(path_) / basename).string();

    write_binary_file_(stem + schema.value_suffix, field);

    if (layout_.include_meta) {
      write_meta_file_(stem + schema.meta_suffix, field, basename, schema);
    }
  }

  void operator()(inqkit::fields::RealField3D const &field,
                  std::string const &basename) const {
    write(field, basename);
  }

private:
  void write_binary_file_(std::string const &filename,
                          inqkit::fields::RealField3D const &field) const {
    auto const filepath = std::filesystem::path(filename);

    if (std::filesystem::exists(filepath) && !options_.overwrite) {
      throw std::runtime_error(
          "RealField3DWriter: file already exists and overwrite=false: " +
          filepath.string());
    }

    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
      throw std::runtime_error(
          "RealField3DWriter: could not open file for writing: " +
          filepath.string());
    }

    out.write(
        reinterpret_cast<char const *>(field.values.data()),
        static_cast<std::streamsize>(field.values.size() * sizeof(double)));

    if (!out) {
      throw std::runtime_error(
          "RealField3DWriter: failed while writing file: " + filepath.string());
    }
  }

  void write_meta_file_(
      std::string const &filename, inqkit::fields::RealField3D const &field,
      std::string const &basename,
      inqkit::detail::grid_layout::RealField3DRawSchema const &schema) const {
    auto const filepath = std::filesystem::path(filename);

    if (std::filesystem::exists(filepath) && !options_.overwrite) {
      throw std::runtime_error("RealField3DWriter: metadata file already "
                               "exists and overwrite=false: " +
                               filepath.string());
    }

    std::ofstream out(filepath);
    if (!out) {
      throw std::runtime_error(
          "RealField3DWriter: could not open metadata file: " +
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
    out << "value_file = " << basename << schema.value_suffix << "\n";
  }

private:
  std::string path_;
  RealField3DLayout layout_;
  RealField3DWriteOptions options_;
};

} // namespace inqkit::io
