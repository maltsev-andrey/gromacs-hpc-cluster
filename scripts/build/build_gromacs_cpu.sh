#!/usr/bin/env bash
# Build CPU-only GROMACS on head-node (head has internet) and install to NFS
# so the offline compute nodes node02..05 can run it.
#
# Builds two binaries:
#   gmx      thread-MPI, for one independent run per node (the main use case)
#   gmx_mpi  OpenMPI, ONLY for the 1/2/4-node scaling test over 1 GbE
#
# Prereqs on head (once):
#   dnf install -y gcc gcc-c++ make hwloc-devel openmpi-devel python3-pip
#
# Set SIMD from the Phase 0 survey of the COMPUTE NODES, not of head:
#   avx2 + fma visible in the VM  -> AVX2_256
#   only avx (ESXi EVC masking)    -> AVX_256
#   neither                        -> SSE4.1

set -euo pipefail

# ---------- settings ----------
GMX_VER="2026.2"
SIMD="AVX2_256"
SRC_ROOT="/root/build"
PREFIX="/nfs/shared/apps/gromacs/${GMX_VER}-cpu-${SIMD,,}"
JOBS=8
BUILD_MPI=1
RUN_TESTS=1        # set 0 if head's CPU lacks the SIMD level chosen above
# ------------------------------

mkdir -p "${SRC_ROOT}"
python3 -m venv "${SRC_ROOT}/cmake-venv"
"${SRC_ROOT}/cmake-venv/bin/pip" install -q --upgrade "cmake>=3.28"
CMAKE="${SRC_ROOT}/cmake-venv/bin/cmake"

cd "${SRC_ROOT}"
[[ -f gromacs-${GMX_VER}.tar.gz ]] || curl -fLO "https://ftp.gromacs.org/gromacs/gromacs-${GMX_VER}.tar.gz"
rm -rf "gromacs-${GMX_VER}" && tar xzf "gromacs-${GMX_VER}.tar.gz"

common_flags=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX="${PREFIX}"
  -DGMX_GPU=OFF
  -DGMX_SIMD="${SIMD}"
  -DGMX_BUILD_OWN_FFTW=ON
  -DGMX_HWLOC=ON
  -DREGRESSIONTEST_DOWNLOAD=ON
)

# 1) thread-MPI build -> gmx
mkdir -p "gromacs-${GMX_VER}/build-tmpi" && cd "gromacs-${GMX_VER}/build-tmpi"
${CMAKE} .. "${common_flags[@]}" -DGMX_MPI=OFF
make -j"${JOBS}"
[[ "${RUN_TESTS}" == 1 ]] && make check -j"${JOBS}"
make install
cd ../..

# 2) OpenMPI build -> gmx_mpi
if [[ "${BUILD_MPI}" == 1 ]]; then
  source /etc/profile.d/modules.sh
  module load mpi/openmpi-x86_64
  mkdir -p "gromacs-${GMX_VER}/build-mpi" && cd "gromacs-${GMX_VER}/build-mpi"
  ${CMAKE} .. "${common_flags[@]}" -DGMX_MPI=ON
  make -j"${JOBS}"
  make install
fi

echo
echo "Installed: ${PREFIX}"
echo "On nodes:  source ${PREFIX}/bin/GMXRC"
echo "Verify on a node: srun -N1 -p compute bash -c 'source ${PREFIX}/bin/GMXRC && gmx --version | grep -E \"SIMD|FFT\"'"
