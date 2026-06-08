# INQ shared configuration
# Edit this file if INQ source, CUDA version, or Python path changes.

export INQ_SOURCE="/local/data/public/skcb2/tddft/inq"
export INQ_SHARE_PATH="$INQ_SOURCE/install/share"
export PSEUDOPOD_SHARE_PATH="$INQ_SOURCE/install/share/pseudopod"

export INQ_CUDA_COMPILER="/lsc/opt/cuda-12.6.2/bin/nvcc"
export INQ_CUDA_ARCH="80"
export INQ_PYTHON_EXE="/local/data/public/skcb2/tddft/venv/bin/python3"
