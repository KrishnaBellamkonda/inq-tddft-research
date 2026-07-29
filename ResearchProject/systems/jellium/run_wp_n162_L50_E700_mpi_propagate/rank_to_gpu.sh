#!/bin/bash
# Bind each MPI rank to its own GPU.
export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}
exec "$@"
