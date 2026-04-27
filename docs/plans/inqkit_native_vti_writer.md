# Plan: Native C++ VTI writer for inqkit

**Author:** investigation continued from `coronene-gs-diagnostics.md` and the run_08 wave-packet check.
**Status:** PROPOSED — awaiting user approval before implementation.
**Last updated:** 2026-04-26.

---

## 1. Goal

Replace the two-step "write `.raw` → run `to_vti.py`" pipeline with a single
direct path: a C++ writer that emits VTK XML ImageData (`.vti`) files in-place
during the run, alongside (or instead of) the existing `.raw` outputs.

The existing slot for this is already in the library:

```
inq-stack/include/inqkit/io/vti_image_data_writer.hpp   (currently a stub: "// TODO: Write this file")
```

## 2. Scope

In scope (this plan):

1. A new header `vti_image_data_writer.hpp` that exposes a `VTIImageDataWriter`
   class implementing the on-disk VTK XML ImageData format for a single
   `RealField3D` (scalar grid) and a single `ComplexField3D` (two scalar
   components: `<field>_real` and `<field>_imag`, in the same `.vti`).
2. Internal use from `RealField3DWriter` and `ComplexField3DWriter`: when a new
   `emit_vti = true` option is enabled the existing writers also call the new
   VTI writer, so call sites do not change.
3. A new diagnostic run `Tutorial/coronene-leed/run_diagnoses/run_09_gs_vti_writer/`
   that re-runs the same GS-only configuration as run_08 with the updated
   writers and produces `.vti` files directly. Visual comparison against
   run_08's Python-converted `.vti` files is the acceptance test.

Out of scope (not this plan):

- Time-series `.pvd` collection files. Each step still produces an independent
  `.vti`; ParaView's `Open File…` recognises a numbered series automatically.
- Compression / `appended` `RawBinary` payloads. We will support `ascii` and
  `binary` (base64-encoded `RawBinary` inline). Compression can come later.
- Multi-rank export. Same restriction as the existing writers (single-rank).

## 3. VTK XML ImageData format reference

VTI is an XML wrapper around regularly-gridded scalar data:

```xml
<?xml version="1.0"?>
<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt64">
  <ImageData WholeExtent="0 NX-1 0 NY-1 0 NZ-1"
             Origin="OX OY OZ"
             Spacing="DX DY DZ">
    <Piece Extent="0 NX-1 0 NY-1 0 NZ-1">
      <PointData Scalars="density">
        <DataArray type="Float64" Name="density" format="ascii|binary">
          ... values ...
        </DataArray>
      </PointData>
    </Piece>
  </ImageData>
</VTKFile>
```

Conventions we will commit to (chosen to match what `inqview.vti.write_vti`
already produces, so the Python-converted and C++-converted files can be
diffed):

| Property | Value |
|---|---|
| `VTKFile/type` | `ImageData` |
| `VTKFile/version` | `1.0` |
| `VTKFile/byte_order` | `LittleEndian` |
| `WholeExtent` | `0 NX-1 0 NY-1 0 NZ-1` (point data) |
| `Origin` | `(origin_x_bohr, origin_y_bohr, origin_z_bohr)` from `RealField3D` |
| `Spacing` | `(dx_bohr, dy_bohr, dz_bohr)` |
| Data location | `PointData` (matches the Python writer; voxel = grid sample) |
| Data type | `Float64` |
| Default `format` | `ascii` (deterministic, diff-friendly, slow) |
| Optional `format` | `binary` = base64-encoded VTK RawBinary block (UInt64 byte-count header + payload) |
| Point order in payload | x fastest, then y, then z (VTK requirement) |

**Layout transform.** Our in-memory `RealField3D::values` is laid out
`flat = ((ix*ny)+iy)*nz + iz` (x slowest, z fastest). VTK ImageData expects
the opposite (x fastest). The VTI writer must reorder while streaming, so the
inner loop reads `values[((ix*ny)+iy)*nz + iz]` while iz, iy, ix are walked in
the order (iz outermost, iy middle, ix innermost). No transpose buffer needed.

For complex fields we emit two arrays inside the same `<PointData>` block:

```xml
<PointData Scalars="<name>_real">
  <DataArray type="Float64" Name="<name>_real" format="ascii">…</DataArray>
  <DataArray type="Float64" Name="<name>_imag" format="ascii">…</DataArray>
</PointData>
```

## 4. Public API

```cpp
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
  explicit VTIImageDataWriter(VTIWriteOptions options = {});

  // Write one scalar grid to <output_path> (typically ending in .vti).
  // Throws on I/O errors or shape mismatch.
  void write_real(inqkit::fields::RealField3D const& field,
                  std::string const& output_path,
                  std::string const& array_name = "density") const;

  // Write a complex grid as two arrays (real_, imag_) in the same .vti.
  void write_complex(inqkit::fields::ComplexField3D const& field,
                     std::string const& output_path,
                     std::string const& array_name = "psi") const;

private:
  VTIWriteOptions options_;
};

}  // namespace inqkit::io
```

The class is **stateless apart from options**; one instance can be reused for
many fields. There is no caching of an output directory — the caller decides
the full path. This mirrors the role of `vtkXMLImageDataWriter` in VTK and
keeps lifetime concerns simple.

### Integration with the existing writers

`RealField3DLayout` and `ComplexField3DLayout` get a new boolean field:

```cpp
struct RealField3DLayout {
  std::string field_name = "field";
  bool include_meta = true;
  bool emit_raw      = true;     // new: keep .raw output (default true for back-compat)
  bool emit_vti      = false;    // new: also write .vti next to .raw
  VTIWriteOptions::Format vti_format = VTIWriteOptions::Format::ascii;
};
```

Default values keep the *current* behaviour byte-for-byte. New runs opt in
explicitly:

```cpp
inqkit::io::RealField3DWriter wr("results/density_gs", {
    .field_name = "density",
    .include_meta = true,
    .emit_raw = false,           // skip .raw if you no longer need it
    .emit_vti = true,
    .vti_format = inqkit::io::VTIWriteOptions::Format::binary,
});
wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
```

Internally `write_impl_` then:

1. (optional) writes `<basename>.raw` exactly as today,
2. (optional) writes `<basename>.meta.txt` exactly as today,
3. (optional) writes `<basename>.vti` via `VTIImageDataWriter`.

The complex writer follows the same pattern with `_real.raw` / `_imag.raw` /
`<basename>.vti` outputs.

## 5. Implementation plan

### Phase A — VTI writer header
1. Implement `VTIImageDataWriter` in `vti_image_data_writer.hpp` (header-only,
   matching the rest of inqkit). Key helpers:
   - `write_xml_prologue_` — emits the `<VTKFile><ImageData WholeExtent=…
     Origin=… Spacing=…><Piece …><PointData …>` opening.
   - `write_data_array_real_` — emits one `<DataArray …>` for a `Float64`
     scalar grid in either ASCII or base64-encoded RawBinary. ASCII path
     iterates with the (iz, iy, ix) loop order to produce x-fastest output;
     binary path needs a small staging buffer.
   - `write_xml_epilogue_` — closes the open tags.
   - `base64_encode_` — minimal base64 encoder (no dependency on a third-party
     lib). For binary `RawBinary` blocks VTK expects an 8-byte little-endian
     UInt64 header (number of payload bytes), then the raw bytes, all
     base64-encoded together.
2. Validate via `static_assert`-style sanity (nx, ny, nz > 0, vector size
   matches, sizeof(double) == 8, host endianness == little).

### Phase B — Wire into existing writers
3. Extend `RealField3DLayout` and `ComplexField3DLayout` with the new fields
   listed above. Defaults preserve existing behaviour.
4. In `RealField3DWriter::write_impl_`, after the raw/meta writes, branch on
   `emit_vti` and call `VTIImageDataWriter::write_real`. The output path is
   `<dir>/<basename>.vti`.
5. Same change in `ComplexField3DWriter::write` (passes through to
   `write_complex`).
6. Make `emit_raw = false` actually skip the binary write — useful for
   storage-tight runs where only the `.vti` is needed.

### Phase C — Smoke test (no INQ run)
7. Write a tiny `Tutorial/_inqkit_tests/vti_writer_smoketest/run.cpp` that
   builds a 4x4x4 `RealField3D` with values `f(ix,iy,iz) = 100*ix + 10*iy + iz`
   and a known origin/spacing, writes both ASCII and binary `.vti`, and reads
   them back with a Python check (`vtkXMLImageDataReader` or just XML+base64
   parsing). Confirms:
   - `Origin`, `Spacing`, `WholeExtent` round-trip correctly,
   - point order is x-fastest (the value at point (3,0,0) is 300, at (0,3,0)
     is 30, at (0,0,3) is 3),
   - binary and ASCII outputs agree element-wise.
8. Add this smoke test to the repo (does not need to run automatically yet —
   a one-shot manual verification before the GS run, per
   `.claude/rules/development-feedback-loop.md`).

### Phase D — Diagnostic GS run with native VTI
9. Create `Tutorial/coronene-leed/run_diagnoses/run_09_gs_vti_writer/` with:
   - `run.cpp` — copy of run_08's GS-only configuration but each
     `RealField3DWriter` is constructed with `.emit_vti = true`. Two passes
     are kept (`density_wp_pre_normalisation`, `density_wp_post_normalisation`)
     plus the GS density and per-orbital densities. ASCII format by default;
     binary as a follow-up if file sizes are a problem.
   - `coronene_centred.xyz` — copied from run_08.
   - `to_vti.py` — left in place but only as a fallback / sanity tool. A new
     `compare_vti.py` is added to diff the C++-emitted `.vti` against
     `inqview`'s Python-converted `.vti` from run_08, point-by-point. Tolerance
     ~1e-12 absolute (ASCII writer keeps 16 sig figs).
10. Run with `inq-run`. Tail `run.log`. Expected outputs (≈ 65 `.vti` files):
    - `results/density_gs/density_t000000.vti`
    - `results/density_gs_orbitals/orbital_NNNN/density_t000000.vti`
    - `results/density_wp_pre_normalisation/density_t000000.vti`
    - `results/density_wp_post_normalisation/density_t000000.vti`
11. Open the wave-packet pre/post `.vti`s in ParaView and check (a) WP
    centred on (0, 0, +12 Bohr), (b) flake at z=0 in the GS density,
    (c) `compare_vti.py` reports max abs diff below tolerance against the
    `inqview`-emitted reference set in run_08.

### Phase E — Commit policy
12. One commit for Phase A+B (the writer + integration). One commit for
    Phase D (the run_09 directory and any tutorial-side opt-in changes).
    No Claude attribution on either commit (per the project rule).

## 6. Validation criteria

A. **Smoke test** (Phase C) — round-trip on the 4×4×4 toy field passes for
   both ASCII and binary formats.

B. **Numerical equivalence** (Phase D step 11) — for every `(ix, iy, iz)`,
   `value_from_cpp_vti == value_from_python_vti` to machine precision.

C. **Visual equivalence** (Phase D step 11) — the WP density blob in the
   C++-emitted `.vti` is centred at `(0, 0, +12 Bohr)`; the GS total density
   has a hexagonal coronene footprint at z=0.

D. **Back-compat** — the existing tutorials that use `RealField3DWriter` /
   `ComplexField3DWriter` produce byte-identical `.raw` and `.meta.txt`
   outputs as before (no `.vti` emitted unless they opt in).

## 7. Risks / open questions

- **ASCII size.** A 120×120×200 grid is ~3M points → ~70 MB ASCII per `.vti`.
  Using `binary` (base64 RawBinary) brings that to ~31 MB per file. For
  run_09 with 65 fields this is the difference between ~4 GB and ~2 GB on
  disk. Default to ASCII for diff-ability, but expose `binary` for the
  per-orbital outputs.
- **Endianness.** Code will assert little-endian host (true on the lovelace
  GPU box). If we ever target a different host, swap bytes in the binary
  encoder.
- **Compression.** VTK supports `<DataArray … format="appended" compressor="…">`
  with zlib; deferring to a follow-up plan.
- **Cell vs. point data.** We are writing PointData. ParaView interpolates
  PointData smoothly; for voxelised iso-surfaces CellData might be preferable.
  Decision: keep PointData for parity with the existing Python pipeline,
  revisit if needed.
- **Float32 fallback.** Some users may prefer Float32 to halve disk usage.
  Not in this plan; the writer takes Float64 directly from `RealField3D` and
  emits Float64. Add as a separate option later if it becomes a real
  requirement.

## 8. References

- VTK XML file formats spec (point ordering, RawBinary header convention):
  Kitware, *VTK File Formats* — https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
- Existing Python-side reference implementation:
  `/local/data/public/skcb2/tddft/inq-stack/python/inqview/vti.py`
- Existing C++ writer contract for `.raw` / `.meta.txt`:
  `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/detail/grid_layout.hpp`
- Centred-cell coordinate convention used by the upstream field builders:
  `inq-stack/include/inqkit/fields/density.hpp` (`fft_shift_index`) and the
  `run_06_centred_writer_check` diagnostic.
- Wave-packet kernel coordinate fix (recent commit `474af1a`):
  `inq-stack/include/inqkit/wavepacket/wavepacket.hpp`.
