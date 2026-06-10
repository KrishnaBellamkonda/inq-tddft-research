# inqkit ENGINE test tier (ADR 0001). Links INQ; GPU; runs locally.
#
# Mirrors the inq-run build path: add_subdirectory(inq) provides the `inq`
# target (and, via INQ's FetchContent, Catch2::Catch2 — version-matched, no
# extra download). Kept in its own subdir so the PURE tier stays INQ-free.
#
# Configure (GPU, matching shared/config.sh):
#   cmake -S . -B build -DENABLE_CUDA=ON \
#         -DCMAKE_CUDA_ARCHITECTURES=80 \
#         -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/nvcc \
#         -DPython_EXECUTABLE=/local/data/public/skcb2/tddft/venv/bin/python3
#   cmake --build build -j && ctest --test-dir build -L engine --output-on-failure

cmake_minimum_required(VERSION 3.21)

option(ENABLE_CUDA "Enable CUDA GPU support" ON)

if(ENABLE_CUDA)
  project(inqkit_engine_tests CXX Fortran CUDA)
  set(ENABLE_GPU ON)
  set(GPU_LANGUAGE CUDA)
else()
  project(inqkit_engine_tests CXX Fortran)
endif()

set(INQ_SOURCE "/local/data/public/skcb2/tddft/inq"
    CACHE PATH "INQ source tree (provides the `inq` and Catch2 targets)")
set(INQKIT_INCLUDE "/local/data/public/skcb2/tddft/inq-stack/include"
    CACHE PATH "inqkit header include directory")

enable_testing()

# Pull in INQ (headers + linked libs + Catch2). EXCLUDE_FROM_ALL: only build
# what the test executables actually need.
add_subdirectory(${INQ_SOURCE} inq EXCLUDE_FROM_ALL)

function(inqkit_add_engine_test name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE inq Catch2::Catch2)
  target_include_directories(${name} PRIVATE ${INQKIT_INCLUDE})
  if(ENABLE_GPU)
    set_source_files_properties(${ARGN} PROPERTIES LANGUAGE ${GPU_LANGUAGE})
  endif()
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES
    LABELS "engine"
    ENVIRONMENT
      "INQ_SHARE_PATH=${INQ_SOURCE}/install/share;PSEUDOPOD_SHARE_PATH=${INQ_SOURCE}/install/share/pseudopod")
endfunction()

inqkit_add_engine_test(test_density_total_engine
  test_engine_main.cpp
  test_density_total_engine.cpp)

inqkit_add_engine_test(test_coord_mapping_engine
  test_engine_main.cpp
  test_coord_mapping_engine.cpp)

inqkit_add_engine_test(test_density_semantics_engine
  test_engine_main.cpp
  test_density_semantics_engine.cpp)

inqkit_add_engine_test(test_orbital_engine
  test_engine_main.cpp
  test_orbital_engine.cpp)

inqkit_add_engine_test(test_observables_writer_engine
  test_engine_main.cpp
  test_observables_writer_engine.cpp)

inqkit_add_engine_test(test_plane_screen_engine
  test_engine_main.cpp
  test_plane_screen_engine.cpp)

inqkit_add_engine_test(test_eigenvalue_dump_engine
  test_engine_main.cpp
  test_eigenvalue_dump_engine.cpp)

inqkit_add_engine_test(test_wp_momentum_two_route_engine
  test_engine_main.cpp
  test_wp_momentum_two_route_engine.cpp)

inqkit_add_engine_test(test_wp_real_space_replica_engine
  test_engine_main.cpp
  test_wp_real_space_replica_engine.cpp)

inqkit_add_engine_test(test_occupations_writer_engine
  test_engine_main.cpp
  test_occupations_writer_engine.cpp)

inqkit_add_engine_test(test_orbital_overlap_engine
  test_engine_main.cpp
  test_orbital_overlap_engine.cpp)

inqkit_add_engine_test(test_leed_engine
  test_engine_main.cpp
  test_leed_engine.cpp)

# Multi-rank test (E01): built like an engine test but launched under mpirun -np 2.
add_executable(test_plane_screen_parallel_engine
  test_engine_main.cpp test_plane_screen_parallel_engine.cpp)
target_link_libraries(test_plane_screen_parallel_engine PRIVATE inq Catch2::Catch2)
target_include_directories(test_plane_screen_parallel_engine PRIVATE ${INQKIT_INCLUDE})
if(ENABLE_GPU)
  set_source_files_properties(test_engine_main.cpp test_plane_screen_parallel_engine.cpp
    PROPERTIES LANGUAGE ${GPU_LANGUAGE})
endif()
add_test(NAME test_plane_screen_parallel_engine
  COMMAND mpirun --oversubscribe -np 2 $<TARGET_FILE:test_plane_screen_parallel_engine>)
set_tests_properties(test_plane_screen_parallel_engine PROPERTIES
  LABELS "engine;parallel"
  ENVIRONMENT
    "INQ_SHARE_PATH=${INQ_SOURCE}/install/share;PSEUDOPOD_SHARE_PATH=${INQ_SOURCE}/install/share/pseudopod")

inqkit_add_engine_test(test_wp_momentum_compute_engine
  test_engine_main.cpp
  test_wp_momentum_compute_engine.cpp)

inqkit_add_engine_test(test_wp_real_space_compute_engine
  test_engine_main.cpp
  test_wp_real_space_compute_engine.cpp)

inqkit_add_engine_test(test_momentum_distribution_engine
  test_engine_main.cpp
  test_momentum_distribution_engine.cpp)

inqkit_add_engine_test(test_state_energy_engine
  test_engine_main.cpp
  test_state_energy_engine.cpp)

# Integration test (IV-M11): free Gaussian WP vs analytic free-particle evolution.
inqkit_add_engine_test(test_free_wp_engine
  test_engine_main.cpp
  test_free_wp_engine.cpp)

# Feature: per-orbital norm diagnostic (every KS orbital ∫|ψ|²≈1).
inqkit_add_engine_test(test_state_norm_engine
  test_engine_main.cpp
  test_state_norm_engine.cpp)
