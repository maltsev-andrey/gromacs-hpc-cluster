# Chapter 2. Build and environment

*Andrey Maltsev · last updated 25 Oct 2026*

GROMACS 2026.2 runs on this cluster in two builds from the same source: a CUDA build on the GPU host and a CPU build on NFS for the offline compute nodes. Both passed their checks, and both are pinned at 2026.2 until the P100 vs V100 comparison is done, so a patch release can never slip between the two halves.

## 2.1 Host survey

Before compiling anything I checked what each host actually exposes, because two of the decisions below depend on it.

| Check | Command | GPU host (gpu-host) | Compute node (node02) |
| --- | --- | --- | --- |
| CPU flags | `lscpu \| grep -oE 'avx2\|fma\|avx512f'` | avx2, fma | avx2, fma, no avx512 |
| Topology | `lscpu` | 12 cores, 24 threads, 1 socket | 12 vCPUs as 12 single-core sockets, no HT inside the VM |
| GPUs | `nvidia-smi --query-gpu=index,name,pci.bus_id,uuid,memory.total --format=csv` | index 0: P100 16GB, bus 03:00.0; index 1: V100 16GB, bus 04:00.0 | none |
| Driver | `nvidia-smi \| head -4` | 550.163.01, CUDA 12.4 max | n/a |
| Toolkit | `/usr/local/cuda-12.4/bin/nvcc --version` | 12.4, already installed | n/a |
| Compiler | `gcc --version` | 11.5.0 | 11.5.0 |
| CMake | `cmake --version` | 3.25 (too old, GROMACS wants 3.28) | 3.25 |
| Internet | `curl -sI https://ftp.gromacs.org \| head -1` | 200 OK | none; only the head node and the GPU host have egress |
| MPI | `ompi_info \| head -1` | n/a | RHEL OpenMPI package, module `mpi/openmpi-x86_64` |

The compute VMs all live on one ESXi host: 2 x Xeon E5-2696 v3 at 2.30 GHz, 36 physical cores, 72 logical threads, 256 GB, on a Huananzhi X99-F8D board. Four nodes at 12 vCPU therefore ask for 48 vCPUs against 36 cores, which is the confound recorded in section 4.3.

Two things this survey settled. ESXi EVC mode can mask AVX2 inside a VM, and it does not on this host, so the node build can use AVX2_256. And the CUDA toolkit on the GPU host was already 12.4, which is the last version driver 550 supports and the last major version that still compiles for Pascal and Volta.

## 2.2 CUDA build on the GPU host

Script: `build_gromacs_gpu.sh`. It installs to `/opt/gromacs/2026.2-cuda` on the GPU host's local NVMe.

| Option | Value | Why |
| --- | --- | --- |
| `GMX_GPU` | CUDA | the only backend that matters for Tesla cards |
| `CUDAToolkit_ROOT` | /usr/local/cuda-12.4 | matches driver 550; CUDA 13 dropped Pascal and Volta |
| `GMX_CUDA_TARGET_SM` | 60;70 | P100 and V100 in one binary, so a card swap needs no rebuild |
| `GMX_SIMD` | AVX2_256 | Xeon E5-2678 v3 (Haswell) |
| `GMX_MPI` | OFF | thread-MPI; one host, one rank, OpenMP threads |
| `GMX_BUILD_OWN_FFTW` | ON | avoids the distro FFTW and its SIMD level |
| `GMX_HWLOC` | ON | correct core pinning |
| `REGRESSIONTEST_DOWNLOAD` | ON | `make check` runs the official suite |

Three rules in the script that are easy to get wrong:

- Install only `cuda-toolkit-12-4`, never the `cuda` metapackage. The metapackage pulls a new driver, and a driver change takes the RAG service on the V100 down with it.
- CMake 3.28 is required and RHEL 9.5 ships 3.25. The script installs a current CMake into a private venv (`pip install cmake`) instead of touching the system package.
- Build as a normal user from a clean environment. With the `cuda_env` venv active, `deactivate` was not even defined in the child shell, so the build ran under `env -i` with an explicit PATH; the regression tests were pinned to the experimental card with `CUDA_VISIBLE_DEVICES=0`.

Build time including `make check` was under 40 minutes with 16 parallel jobs, leaving cores for the RAG stack.

## 2.3 CPU builds for the offline compute nodes

Script: `build_gromacs_cpu.sh`, run on the head node because the nodes have no internet and the build downloads FFTW and the test suite. It installs to `/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256`, so a node only needs `source .../bin/GMXRC`.

Two binaries come out of it:

| Binary | Build | Used for |
| --- | --- | --- |
| `gmx` | thread-MPI (`GMX_MPI=OFF`) | one independent simulation per node; the production build |
| `gmx_mpi` | OpenMPI (`GMX_MPI=ON`, module `mpi/openmpi-x86_64`) | the 1/2/4-node scaling test in chapter 4, nothing else |

The SIMD level comes from the survey of the compute nodes, not of the head node, because head and nodes could sit on different ESXi hosts. The build directory has to be on a local disk (`$HOME/build`), not on NFS: compiling GROMACS writes thousands of small object files, and NFS makes that slow. The install prefix on NFS was created once as root and handed to the `ansible` user, so nothing in the build needs sudo.

Verification on a node:

```
srun -N1 -p compute bash -c 'source /nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/GMXRC && gmx --version | grep -E "SIMD|MPI library|FFT"'
# MPI library: thread_mpi   SIMD instructions: AVX2_256   CPU FFT library: fftw-3.3.10-sse2-avx-avx2-avx2_128
srun -N1 -p compute bash -c 'module load mpi/openmpi-x86_64; ldd /nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/gmx_mpi | grep -c "not found"'
# 0  (grep exits 1 on a zero count, which is not an error)
```

## 2.4 Regression test results

| Build | Suite | Result |
| --- | --- | --- |
| CUDA, gpu-host | `make check`, 106 tests including the MPI-rank, PME and multi-sim suites, on the P100 | 106 of 106 passed |
| CPU thread-MPI, head-node | `make check` | passed; install completed |
| CPU OpenMPI, head-node | not run separately; validated by the scaling jobs in chapter 4 | 34 of 36 runs completed |

The script runs with `set -euo pipefail`, so reaching `make install` is itself proof that `make check` exited cleanly.

## 2.5 Analysis environment

Analysis runs on the head node in its own venv, never inside `cuda_env` on the GPU host, so the benchmark environment stays untouched:

```
python3 -m venv ~/venvs/gmx-analysis
source ~/venvs/gmx-analysis/bin/activate
pip install pandas numpy matplotlib
cd /nfs/shared/gromacs/results
python /nfs/shared/projects/gromacs/analyze_bench.py     # GPU matrices: ns/day, CPU + GPU power, Wh/ns, charts
python /nfs/shared/projects/gromacs/analyze_cpu.py       # phase 2: node sweep, scaling efficiency, chart
```

`analyze_bench.py` merges several CSVs and RAPL logs, uses each run's exact start time from the CSV, lets clean reruns replace disturbed points, and drops any run more than 5% below the median of its group, listing what it dropped. Everything configurable sits in the settings block at the top of each script; no command-line arguments.

Phase 5 will add MDAnalysis to the same venv.
