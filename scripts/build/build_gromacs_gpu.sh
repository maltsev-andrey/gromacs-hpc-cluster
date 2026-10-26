#!/usr/bin/env bash
# Build GROMACS with CUDA on gpu-host (local NVMe install).
# Targets sm_60 (P100) AND sm_70 (V100), so the same binary keeps working
# after the P100 -> V100 32GB swap. No rebuild needed.
#
# Prereqs (run once, as root):
#   dnf install -y gcc gcc-c++ make git hwloc-devel python3-pip
#   NVIDIA RHEL9 repo, then:  dnf install -y cuda-toolkit-12-4
#   !! Install ONLY cuda-toolkit-12-4. Do NOT install the "cuda" metapackage:
#   it pulls a new driver and can break the RAG service on the V100.
#   Driver 550 = CUDA 12.4 max, so 12.4 is the matching toolkit.
#   Do NOT use CUDA 13.x: it dropped Pascal and Volta.

set -euo pipefail

# ---------- settings ----------
GMX_VER="2026.2"                       # check ftp.gromacs.org for the latest 2026.x patch
CUDA_HOME="/usr/local/cuda-12.4"
SRC_ROOT="/data/build"                 # on local NVMe
PREFIX="/opt/gromacs/${GMX_VER}-cuda"
SIMD="AVX2_256"                        # Xeon E5-2678 v3 (Haswell)
JOBS=16                                # leave cores for the RAG stack
# ------------------------------

export PATH="${CUDA_HOME}/bin:${PATH}"

echo "== toolchain =="
gcc --version | head -1
nvcc --version | tail -1
nvidia-smi --query-gpu=index,name,pci.bus_id,memory.total,driver_version --format=csv

# CMake >= 3.28 required, RHEL 9.5 ships 3.26 -> private venv
mkdir -p "${SRC_ROOT}"
python3 -m venv "${SRC_ROOT}/cmake-venv"
"${SRC_ROOT}/cmake-venv/bin/pip" install -q --upgrade "cmake>=3.28"
CMAKE="${SRC_ROOT}/cmake-venv/bin/cmake"
${CMAKE} --version | head -1

cd "${SRC_ROOT}"
[[ -f gromacs-${GMX_VER}.tar.gz ]] || curl -fLO "https://ftp.gromacs.org/gromacs/gromacs-${GMX_VER}.tar.gz"
rm -rf "gromacs-${GMX_VER}" && tar xzf "gromacs-${GMX_VER}.tar.gz"
mkdir -p "gromacs-${GMX_VER}/build-cuda" && cd "gromacs-${GMX_VER}/build-cuda"

${CMAKE} .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DGMX_GPU=CUDA \
  -DCUDAToolkit_ROOT="${CUDA_HOME}" \
  -DGMX_CUDA_TARGET_SM="60;70" \
  -DGMX_SIMD="${SIMD}" \
  -DGMX_MPI=OFF \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_HWLOC=ON \
  -DREGRESSIONTEST_DOWNLOAD=ON

make -j"${JOBS}"
# Regression tests use GPU 0 of the visible set: pin them to the experimental card
CUDA_VISIBLE_DEVICES=0 make check -j"${JOBS}"
sudo make install

echo
echo "Installed: ${PREFIX}"
echo "Activate:  source ${PREFIX}/bin/GMXRC"
