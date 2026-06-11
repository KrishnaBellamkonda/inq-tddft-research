// Catch2 main for the inqkit ENGINE test tier.
//
// Mirrors INQ's own unit-test driver (inq/src/main/unit_tests_main.cpp): the
// INQ environment (MPI + FFTW) must be initialised before any INQ operation,
// so we provide our own main rather than using Catch2WithMain.

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <input/environment.hpp>
#include <fftw3.h>

int main(int argc, char *argv[]) {
  inq::input::environment::global();  // initialise MPI
  int result = Catch::Session().run(argc, argv);
  fftw_cleanup();
  return result;
}
