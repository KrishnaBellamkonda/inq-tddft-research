cmake_minimum_required(VERSION 3.16)
project(free_gaussian_wp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3 -march=native)

# FFTW3 double-precision via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

add_executable(run_wp run.cpp)
target_include_directories(run_wp PRIVATE ${FFTW3_INCLUDE_DIRS})
target_link_libraries(run_wp PRIVATE ${FFTW3_LIBRARIES} m)
target_compile_options(run_wp PRIVATE ${FFTW3_CFLAGS_OTHER})
