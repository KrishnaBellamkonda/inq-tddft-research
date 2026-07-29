#!/usr/bin/env bash
# mpi_wrapper.sh — bind each OpenMPI rank to a distinct CUDA device by
# setting CUDA_VISIBLE_DEVICES from OMPI_COMM_WORLD_LOCAL_RANK.
#
# INQ's GPU init does not call cudaSetDevice and uses the default device
# 0 visible to the process. Without this wrapper, both ranks of a 2-rank
# mpirun would target device 0 and we'd lose the parallel benefit.
set -e
export CUDA_VISIBLE_DEVICES="${OMPI_COMM_WORLD_LOCAL_RANK:-0}"
echo "[wrapper] rank ${OMPI_COMM_WORLD_RANK:-?} (local ${OMPI_COMM_WORLD_LOCAL_RANK:-?}) -> CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
exec ./run
