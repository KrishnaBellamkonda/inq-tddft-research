/*
 * This module deals with the writing of RealField3D as defined in the class
 * fields::RealField3D. It takes the schema for these files given in grid_layout.hpp
 * and writes the data files. The go to method is to write a 
 * - binary file (for 3D field data flattened as a 1D vector<>)
 * - meta data file (for recreating the grid coordinates and values easily)
 *
 * */

#pragma once

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace inqkit::io {

// This class handles the o
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

  // Ground-state overload: caller supplies the full basename.
  void write(inqkit::fields::RealField3D const &field,
             std::string const &basename) const {
    write_impl_(field, basename, std::nullopt);
  }

  // Real-time overload: generates basename = field_name + "_t{step:06d}".
  // Writes time_au into the sidecar metadata so Python readers see time
  // ordering.
  void write(inqkit::fields::RealField3D const &field, double time_au,
             int step) const {
    auto const basename =
        layout_.field_name + inqkit::detail::grid_layout::step_suffix(step);
    write_impl_(field, basename, time_au);
  }

  void operator()(inqkit::fields::RealField3D const &field,
                  std::string const &basename) const {
    write(field, basename);
  }

private:
  // This function is called multiple times in this class.
  // This is an umbrella function that runs both
  // write_binary_file and write_meta_data for a given configuration
  void write_impl_(inqkit::fields::RealField3D const &field,
                   std::string const &basename,
                   std::optional<double> time_au) const {

    // Sanity checks
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

    // Create the directory if not present
    std::filesystem::create_directories(path_);

    auto const schema = inqkit::detail::grid_layout::real_field_3d_raw_schema();
    auto const stem = (std::filesystem::path(path_) / basename).string();

    write_binary_file_(stem + schema.value_suffix, field);

    if (layout_.include_meta) {
      write_meta_file_(stem + schema.meta_suffix, field, basename, schema,
                       time_au);
    }
  }


  // Writes the 3D object as a binary file
  void write_binary_file_(std::string const &filename,
                          inqkit::fields::RealField3D const &field) const {
    auto const filepath = std::filesystem::path(filename);

    if (std::filesystem::exists(filepath) && !options_.overwrite) {
      throw std::runtime_error(
          "RealField3DWriter: file already exists and overwrite=false: " +
          filepath.string());
    }

    // Sets the output stream to binary and writes the
    // field.values() array as a binary file
    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
      throw std::runtime_error(
          "RealField3DWriter: could not open file for writing: " +
          filepath.string());
    }

    // Writes the field as a flattened array of doubles
    out.write(
        reinterpret_cast<char const *>(field.values.data()),
        static_cast<std::streamsize>(field.values.size() * sizeof(double)));

    if (!out) {
      throw std::runtime_error(
          "RealField3DWriter: failed while writing file: " + filepath.string());
    }
  }


  // This is an important function that writes a meta files using the data
  // gathered in the RealField3D and schema supplied to it. The meta data
  // must be sufficient enough to recreate the grid coordinates, the 
  // layout of the results (how does the flattened array map to the 3D vector)
  // and information about what it is (orbital density or total density) that 
  // is being stored in this file.
  // TODO: In the future, this can be expanded to include
  // 1. The entire simulation configuration information 
  //    - ion positions
  //    - no. of electrons, no. of extra states, no. of extra electrons, temp. etc. 
  //    - ground state configuration
  //    - real time propagation information 
  void write_meta_file_(
      std::string const &filename, inqkit::fields::RealField3D const &field,
      std::string const &basename,
      inqkit::detail::grid_layout::RealField3DRawSchema const &schema,
      std::optional<double> time_au) const {
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

    if (time_au.has_value()) {
      out << "time_au = " << std::fixed << *time_au << "\n";
    }

    out << "value_file = " << basename << schema.value_suffix << "\n";
  }

private:
  std::string path_;
  RealField3DLayout layout_;
  RealField3DWriteOptions options_;
};

} // namespace inqkit::io
