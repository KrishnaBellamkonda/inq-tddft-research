// ============================================================================
// Canonical paths for the Tsubonoya 2014 coronene paper-replica configuration.
//
// The numeric configuration (cell, cutoff, WP parameters, RT cadence, screen
// positions, …) lives in the existing header:
//
//     inq-stack/include/inqkit/config/tsubonoya_2014_coronene.hpp
//
// Per-system *paths* live here so multiple runs (GS-save, RT-load, future
// parameter sweeps) all reference the same on-disk geometry file and the
// same electrons.save() checkpoint directory.
//
// Layout under ResearchProject/systems/coronene/:
//
//     configurations/
//       tsubonoya_2014_paper_replica/
//         coronene_centred.xyz     <- canonical geometry
//         paths.hpp                <- this file
//     checkpoints/                 <- gitignored, holds electrons.save() dirs
//     run_<name>/                  <- per-experiment run folders
//
// All paths below are absolute so a run.cpp compiled in any working
// directory can reach the same checkpoint.
// ============================================================================
#pragma once

#include <string>

namespace coronene::paper_replica {

inline constexpr char const *GEOMETRY_XYZ =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/"
    "configurations/tsubonoya_2014_paper_replica/coronene_centred.xyz";

// Where electrons.save() writes the GS checkpoint and electrons.load() reads
// it back. The directory is gitignored.
inline constexpr char const *GS_CHECKPOINT_DIR =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/"
    "checkpoints/tsubonoya_2014_paper_replica_gs";

}  // namespace coronene::paper_replica
