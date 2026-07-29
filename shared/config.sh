# INQ shared configuration
# Edit this file if INQ source, CUDA version, or Python path changes.

# INQ_SOURCE may be overridden in the environment to build against a fork
# (e.g. inq-study). Defaults are unchanged when nothing is pre-set, so existing
# inq-run usage is byte-for-byte identical. Share paths default off INQ_SOURCE but
# can be pinned separately (a fork without its own install/share can reuse inq's).
export INQ_SOURCE="${INQ_SOURCE:-/local/data/public/skcb2/tddft/inq}"
export INQ_SHARE_PATH="${INQ_SHARE_PATH:-$INQ_SOURCE/install/share}"
export PSEUDOPOD_SHARE_PATH="${PSEUDOPOD_SHARE_PATH:-$INQ_SOURCE/install/share/pseudopod}"

export INQ_CUDA_COMPILER="/lsc/opt/cuda-12.6.2/bin/nvcc"
export INQ_CUDA_ARCH="80"
export INQ_PYTHON_EXE="/local/data/public/skcb2/tddft/venv/bin/python3"
