# inqkit C++ test harness (ADR 0001: external tests/cpp + Catch2 + ctest).
#
# Two tiers (ADR 0001):
#   pure   - no INQ link, C++17 only, fast; runs in every CI (ADR 0002).
#   engine - links INQ, GPU; runs locally / occasionally.
#
# This file currently wires the PURE tier only. The engine tier is added
# once the first engine test (fields) proves the INQ-link CMake path.
#
# Run:   cmake -S . -B build && cmake --build build && ctest --test-dir build
#        ctest --test-dir build -L pure          # pure tier only

cmake_minimum_required(VERSION 3.20)
project(inqkit_tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

# ---------------------------------------------------------------------------
# Catch2 v3. Prefer an installed package; otherwise reuse the Catch2 source
# that the INQ build already downloaded (no network needed). Override with
# -DCATCH2_SOURCE_DIR=<path> if INQ's build tree has been cleaned.
# ---------------------------------------------------------------------------
find_package(Catch2 3 QUIET)

if(NOT Catch2_FOUND)
  if(NOT DEFINED CATCH2_SOURCE_DIR)
    set(CATCH2_SOURCE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../inq/build/_deps/catch2-src")
  endif()
  if(EXISTS "${CATCH2_SOURCE_DIR}/CMakeLists.txt")
    message(STATUS "inqkit-tests: using Catch2 source at ${CATCH2_SOURCE_DIR}")
    add_subdirectory("${CATCH2_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/catch2-build")
  else()
    # Final fallback: fetch Catch2 v3 (used by CI, where neither an installed
    # package nor INQ's build tree is present). Needs network at configure time.
    message(STATUS "inqkit-tests: fetching Catch2 v3 via FetchContent")
    include(FetchContent)
    FetchContent_Declare(Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG        v3.5.4
      GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(Catch2)
  endif()
endif()

# inqkit is header-only; the pure tier just needs the include path.
set(INQKIT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../include")

# ---------------------------------------------------------------------------
# Pure-tier tests. Each .cpp here must NOT include <inq/...>; only inqkit
# headers that are free of INQ (detail/, real_field_3d, center_of_density, ...).
# ---------------------------------------------------------------------------
add_library(inqkit_pure_iface INTERFACE)
target_include_directories(inqkit_pure_iface INTERFACE "${INQKIT_INCLUDE_DIR}")
target_link_libraries(inqkit_pure_iface INTERFACE Catch2::Catch2WithMain)

function(inqkit_add_pure_test name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE inqkit_pure_iface)
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS "pure")
endfunction()

# Sources mirror inq-stack/include/inqkit/<module>/.
inqkit_add_pure_test(test_grid_layout inqkit/detail/test_grid_layout.cpp)
inqkit_add_pure_test(test_fft_shift inqkit/detail/test_fft_shift.cpp)
inqkit_add_pure_test(test_vec3 inqkit/detail/test_vec3.cpp)
inqkit_add_pure_test(test_density_delta inqkit/observables/test_density_delta.cpp)
inqkit_add_pure_test(test_center_of_density inqkit/observables/test_center_of_density.cpp)
inqkit_add_pure_test(test_minimum_observable_set inqkit/observables/test_minimum_observable_set.cpp)
inqkit_add_pure_test(test_real_field_writer inqkit/io/test_real_field_writer.cpp)
inqkit_add_pure_test(test_complex_field_writer inqkit/io/test_complex_field_writer.cpp)
inqkit_add_pure_test(test_vti_writer inqkit/io/test_vti_writer.cpp)
inqkit_add_pure_test(test_jellium_shells inqkit/jellium/test_jellium_shells.cpp)
