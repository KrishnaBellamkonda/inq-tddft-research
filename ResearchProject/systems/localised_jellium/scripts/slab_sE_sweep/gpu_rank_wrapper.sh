#!/usr/bin/env bash
# gpu_rank_wrapper.sh — bind each MPI rank to a distinct physical GPU.
#
# INQ's gpurun layer always talks to CUDA device 0 (it never calls cudaSetDevice
# per MPI rank). So a bare `mpirun -np 2 ./run` would place BOTH ranks on GPU0.
# This wrapper masks the visible devices down to exactly ONE per rank, using the
# OpenMPI node-local rank, so rank 0 -> physical GPU0, rank 1 -> physical GPU1.
# The result is a genuine 2-GPU domain/state decomposition: the orbital working
# set (the memory that overflows a single 24 GB A30 at large L_z) is split across
# both cards, eliminating the UVM host<->device thrash.
#
# Usage:  mpirun.openmpi -np 2 --bind-to none ./gpu_rank_wrapper.sh ./run
set -euo pipefail
export CUDA_VISIBLE_DEVICES="${OMPI_COMM_WORLD_LOCAL_RANK:-0}"
exec "$@"
