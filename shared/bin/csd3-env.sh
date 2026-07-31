# csd3-env.sh — the ONE toolchain environment for building and running INQ on CSD3.
# Source (do not execute) from configure-inq.sh, build-inq.slurm, and job scripts:
#   . "$REPO_ROOT/shared/bin/csd3-env.sh"
#
# Why this file exists — `module load rhel8/default-amp` pulls in cuda/11.4 as a
# hard requirement, and that leaves BOTH toolkits visible:
#   CPATH     = .../cuda/12.1/include : ... : .../cuda-11.4.../include
#   CUDA_HOME = .../cuda-11.4...
# nvcc 12.1 then implicitly includes its own crt/sm_80_rt.hpp AND picks up 11.4's
# copy off CPATH, so every sm_80 translation unit dies with
#   crt/sm_80_rt.hpp(141): error: more than one instance of overloaded function
#   "__nv_associate_access_property_impl" has "C" linkage
# That killed build job 32352675 at 92% (inq_executable, _pinq, all tests and
# benchmarks). The 11.4 modules cannot be unloaded — they are locked dependencies
# of rhel8/default-amp — so the entries are stripped from the paths instead.
# Verified 2026-07-29 by re-running the exact failing nvcc command for
# src/main/unit_tests_main.cpp: fails with 11.4 on CPATH, compiles without it.
#
# NOTE: CSD3 also has cuda/12.8.1 (module cuda/12.8.1/gcc/kdeps6ab, visible after
# `module load rhel8/ampere-env/2025-06-01`) with a clean CPATH out of the box.
# It is NOT used here because that environment ships no openmpi that INQ can use;
# revisit if an MPI becomes available there.

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-amp      # openmpi 4.1.1 / gcc-9.4.0, slurm, A100 stack

# HOST COMPILER — must be gcc >= 9. rhel8/default-amp loads no compiler module, so
# the host compiler defaults to system gcc 8.5.0, where std::filesystem lives in a
# SEPARATE library (-lstdc++fs) that INQ's CMake does not link. Every executable
# then fails to link with dozens of
#   undefined reference to `std::filesystem::create_directories(...)'
#   undefined reference to `std::filesystem::__cxx11::path::_M_split_cmpts()'
# which is what killed build job 32353243 (6 targets). gcc 9+ has <filesystem> in
# libstdc++ proper. 9.4.0 is chosen to match the gcc-9.4.0 ABI of the openmpi/ucx
# stack that rhel8/default-amp loads (zen2 view, same as those dependencies).
# Verified 2026-07-29: gcc 8.5 -> 4 undefined refs on a std::filesystem probe;
# gcc 9.4.0 -> 0, links clean.
module load gcc/9.4.0/gcc-11.2.0-tfj3hud

module load cuda/12.1

# Strip the CUDA 11.4 that rhel8/default-amp dragged in, from EVERY search path.
#
# CPATH alone is not enough. LIBRARY_PATH is what gcc/ld searches at LINK time and
# CMAKE_PREFIX_PATH is what find_package(CUDAToolkit) searches, so leaving those two
# populated links the binaries against CUDA 11.4's runtime even though nvcc is 12.1.
# That produced a build that compiled and linked cleanly but was unrunnable:
#   ./run: error while loading shared libraries: libcufft.so.10
# (CUDA 11.x cuFFT soname; 12.x ships libcufft.so.11). Both inq/install/bin/inq and
# the run binary came out NEEDING libcufft.so.10 / libcudart.so.11.0 / libcublas.so.11.
# Verified 2026-07-29 via objdump -p ... | grep NEEDED. Scrub all of them.
for _v in CPATH LD_LIBRARY_PATH LIBRARY_PATH CMAKE_PREFIX_PATH PKG_CONFIG_PATH LD_RUN_PATH; do
  _cur="$(eval printf '%s' "\"\${$_v:-}\"")"
  [ -n "$_cur" ] || continue
  export "$_v=$(printf '%s' "$_cur" | tr ':' '\n' | grep -vE 'cuda-11\.4|cuda/11\.4' | paste -sd: -)"
done
unset _v _cur

export CUDA_HOME=/usr/local/software/cuda/12.1
export CUDA_PATH="$CUDA_HOME"
export CUDAToolkit_ROOT="$CUDA_HOME"
export CMAKE_PREFIX_PATH="$CUDA_HOME${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export LIBRARY_PATH="$CUDA_HOME/lib64${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

export INQ_NVCC="$CUDA_HOME/bin/nvcc"
export INQ_CUDA_ARCH="${INQ_CUDA_ARCH:-80}"   # A100

# Fail loudly rather than silently building against the wrong toolkit.
for _v in CPATH LD_LIBRARY_PATH LIBRARY_PATH CMAKE_PREFIX_PATH; do
  if eval printf '%s' "\"\${$_v:-}\"" | grep -qE 'cuda-11\.4|cuda/11\.4'; then
    echo "csd3-env.sh: ERROR — CUDA 11.4 still on $_v; build would mis-link." >&2
    return 1 2>/dev/null || exit 1
  fi
done
unset _v
