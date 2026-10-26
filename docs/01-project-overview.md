# Chapter 1. Project overview

*Andrey Maltsev · last updated 25 Oct 2026*

## 1.1 Overview

I run classical molecular dynamics (MD) with GROMACS on my home HPC cluster: real protein systems, GPU-accelerated on the GPU host, with independent ensemble runs on the CPU compute nodes. The project has two goals. First, measure honestly what each layer of hardware and offload buys me, in ns/day and in energy per simulated nanosecond. Second, use that baseline to run a real physics study: protein stability under different temperatures.

MD integrates Newton's equations of motion for every atom in the system, typically 80k to 100k atoms at a 2 fs timestep. Each step computes short-range nonbonded forces, long-range electrostatics via PME (particle mesh Ewald, an FFT-based method), bonded forces, then updates positions under constraints. That gives four distinct workloads, which GROMACS can place on the CPU or the GPU separately. This makes MD an ideal subject for measured offload experiments.

**Why MD and not something else.** I evaluated seven candidate projects: seismic full-waveform inversion, weather downscaling, SETI, MD, cosmological N-body, financial Monte Carlo, and extending my gravitational-wave search. MD won on hardware fit:

- My cluster is a throughput machine, not a capability machine: 4 small VMs on 1 GbE with about 112 µs latency. Strong-scaling one huge simulation does not work here. Independent replicas do, and MD is naturally ensemble-shaped.
- A single GPU runs a medium MD system completely, with no multi-GPU or multi-node dependency.
- Input data is a few MB from the Protein Data Bank, so my offline compute nodes need no data pipeline.
- It is a new domain with real physics (force fields, integrators, PME, thermostats), and it is not another signal search like my LIGO matched-filter work.

I sequenced it before the weather-downscaling project because MD needs almost no data engineering and does not compete with my ongoing swiss-hpc-jobs ETL for head-node resources.

## 1.2 Infrastructure and constraints

The project uses one usable GPU and four small CPU nodes, and the design follows from those limits rather than from what a datacenter would offer.

| Host | Role | Hardware | Network | Constraint that matters |
| --- | --- | --- | --- | --- |
| gpu-host | GPU host, bare metal | Xeon E5-2678 v3 (12 cores, 24 threads), 128 GB RAM, ~931 GB local NVMe, Tesla P100 16GB + Tesla V100 16GB | Home LAN 192.168.10.0/24 only | V100 is reserved for my RAG assistant (~15 GB VRAM in use), so only the P100 is available for MD. Not in Slurm yet. Also a worker node of my RKE2 cluster: the RAG runs there as a vLLM pod, and the NVIDIA GPU Operator exposes both GPUs to Kubernetes |
| head-node | Head node, NFS server, Slurm controller | ESXi VM, 32 GB RAM | Dual-homed: 10.10.10.0/24 and home LAN | Only host with internet. All builds and downloads happen here or on the GPU host |
| node02 to 05 | Compute nodes, Slurm partition compute | ESXi VMs, 12 vCPU and ~7.6 GB RAM each, AVX2 exposed | 1 GbE fabric 10.10.10.0/24, ~112 µs latency | No internet egress. About 630 MB RAM per MPI rank. Communication-heavy codes stop scaling after 2 nodes |

All of these VMs run on one ESXi host: 2 x Xeon E5-2696 v3, 36 physical cores, 72 logical threads, 256 GB. Four compute nodes at 12 vCPU each therefore ask for 48 vCPUs, which is more than the host has physical cores; section 4.3 records where that mattered.

Decisions that follow from this table:

- **One run per GPU, one run per CPU node.** I don't run multi-node GROMACS for production. PME needs all-to-all FFT communication, which TCP on 1 GbE cannot sustain.
- **RAG card isolation is enforced, not assumed.** Every GPU job pins `CUDA_VISIBLE_DEVICES` to the P100 and uses `CUDA_DEVICE_ORDER=PCI_BUS_ID`. The benchmark script also checks the card name and its idle VRAM before starting, and aborts otherwise.
- **Binaries for the compute nodes live on NFS.** I build on the head node and install to `/nfs/shared/apps/gromacs/`, so the offline nodes only source the environment file.
- **Hot data stays on local NVMe.** My earlier GPU pipeline work showed that moving data from NFS to local NVMe doubled real throughput, so all benchmark inputs and outputs sit on the GPU host's NVMe.

Since 24 Oct the P100 is card B of the host's GPU-sharing design: Kubernetes with Kueue owns it, ollama runs on it as a low-priority service, and every host batch job claims it with `gpu-run.sh <name> -- <command>`, which evicts ollama, waits for the card to drain, and releases it on exit. The bench scripts no longer stop containers or services themselves. The V100 32GB that was meant to replace the P100 was canceled by the seller; when a card does arrive, all P100 measurements become the before half of the comparison.

## 1.3 Software stack and build

I built GROMACS 2026.2 from source in two variants. The CUDA build passed all 106 regression tests, and the CPU build is installed on NFS and verified on a compute node (thread-MPI gmx and gmx_mpi, AVX2_256). I pinned 2026.2 for the entire P100 vs V100 comparison and will upgrade to a newer patch release only after it.

| Build | Host | Install path | Key options | Used for |
| --- | --- | --- | --- | --- |
| CUDA | gpu-host | `/opt/gromacs/2026.2-cuda` | CUDA 12.4, `GMX_CUDA_TARGET_SM=60;70`, AVX2_256, thread-MPI, own FFTW | All GPU runs |
| CPU, thread-MPI (`gmx`) | built on head-node | `/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256` | No GPU, AVX2_256, own FFTW | One independent run per compute node |
| CPU, OpenMPI (`gmx_mpi`) | built on head-node | same prefix | `GMX_MPI=ON` | Only the 1/2/4-node scaling test |

Why these choices:

- **CUDA 12.4, not 13.** Driver 550 supports CUDA up to 12.4, and CUDA 13 dropped Pascal and Volta entirely. I installed only `cuda-toolkit-12-4`, never the `cuda` metapackage, which would pull a new driver and could take down the RAG service.
- **Two compute capabilities in one binary.** Compiling for sm_60 and sm_70 means the same build keeps working after the swap, with no rebuild between the two halves of the comparison.
- **AVX2_256 for the nodes.** I checked the CPU flags inside the VMs first, because ESXi EVC mode can mask AVX2. It doesn't on my host: avx, avx2 and fma are all exposed.
- **Thread-MPI as the main CPU build.** One run per node needs no real MPI. The OpenMPI build exists only to measure where the 1 GbE fabric breaks.

The toolchain needed one workaround: GROMACS 2026 requires CMake 3.28, and RHEL 9.5 ships 3.25. Both build scripts install a current CMake into a private venv, so I didn't touch the system package.

The build steps, the regression results and the analysis environment are in chapter 2.

## 1.4 Benchmark methodology

I measure every offload level at every thread count, three times each, on two standard input systems, and I record GPU and CPU power alongside throughput.

**Input systems**

| System | Atoms | Source | Constraints | Electrostatics | Why I use it |
| --- | --- | --- | --- | --- | --- |
| benchMEM | 81,743 | MPI-NAT Grubmüller benchmark set (tpr built with GROMACS 4.6, 2013) | Legacy, chains of coupled constraints | PME, cutoff 1.0 nm | Widely published, so my numbers are comparable with other hardware |
| ADH dodec | 95,561 | GROMACS developers' benchmark set (Páll and Alekseenko 2024), tpr rebuilt with `gmx grompp` 2026.2 | h-bonds | PME, cutoff 0.9 nm | Modern topology that allows the fully GPU-resident update step |

I compare each system across GPUs, never the two systems against each other: different cutoffs give a different cost per atom.

**Offload matrix**

1. `cpu_only`: everything on the CPU.
2. `nb_gpu`: short-range nonbonded forces on the GPU.
3. `nb_pme_gpu`: plus PME on the GPU.
4. `nb_pme_bonded_gpu`: plus bonded forces on the GPU.
5. `full_gpu_resident`: plus update and constraints on the GPU (ADH only; benchMEM refuses it).

Each level runs at 4, 8, 12 and 16 OpenMP threads with one thread-MPI rank and `-pin on`. The CPU has 12 physical cores, so 16 threads deliberately tests hyperthreading.

**Timing rules**

- 20,000 steps per run, counters reset at step 10,000, so only the second half counts.
- GROMACS tunes the PME grid during the first steps. Tuning must finish before the reset, or mdrun aborts. On benchMEM it finished at step ~4,900, leaving about 2x margin.
- Three repeats per point. The coefficient of variation stayed below 1% in all 48 benchMEM runs. On ADH, 4 of 60 runs were disturbed by other activity on the host. I reran those points with both ollama containers stopped, and the reruns replace the originals. The analysis drops any run more than 5% below the median of its group and lists it, since clean repeats vary by less than 1%.
- Every matrix runs under gpu-run.sh, so ollama is evicted from the card and a competing host service is paused for the duration; the matrices before 24 Oct relied on checking by hand that nothing touched the P100, which is how four disturbed runs got through.

**Energy measurement**

- GPU board power: `nvidia-smi power.draw`, sampled every 500 ms during each run.
- CPU package power: the RAPL counter `/sys/class/powercap/intel-rapl:0/energy_uj`, logged every 5 s and matched to each run's time window.
- Idle baseline: 5 minutes after the matrix with nothing running. It came out at 16.6 W for the CPU package and about 25 W for the idle P100.
- Scope: CPU package plus GPU board only. The board exposes no DRAM RAPL domain, and fans, NVMe and PSU losses are not covered. Energy per simulated nanosecond = (CPU W + GPU W) × 24 / (ns/day), in Wh/ns.

## 1.5 Results so far

The measured picture, with the full data in chapters 3 and 4:

| What | Result | Where |
| --- | --- | --- |
| Best P100 config on benchMEM | 111.4 ns/day, nonbonded + PME on the GPU, 12 threads, 6.5x CPU-only | 3.1 |
| Best P100 config on ADH | 120.2 ns/day fully GPU-resident on a quiet host, 8.6x CPU-only, 4 threads already give 96% of it | 3.2, 3.6 |
| Energy | GPU-resident mode uses 7.3x less energy per simulated ns than CPU-only (34 vs 247 Wh/ns, CPU package + GPU board) | 3.3 |
| PME tuning | Worth +52% (benchMEM) and +85% (ADH) when PME runs on the CPU, nothing once PME is on the GPU | 3.4 |
| Host noise | Four of 60 ADH runs were disturbed by a host process (brief_worker); before the scheduler the noise floor was ~1.5% in the most GPU-bound config | 3.5, 3.6 |
| One compute VM | 96% to 97% of the bare-metal GPU host in CPU-only mode | 4.1 |
| Multi-node gmx_mpi | Slower than one node at every node count over 1 GbE: 27% to 40% efficiency on 2 nodes (clean), 14% to 16% on 4 (confounded by CPU oversubscription, see 4.3) | 4.2, 4.3 |
| Production rule | One replica per node with the thread-MPI build; one GPU replica in resident mode on 4 threads | 4.4 |

## 1.6 Issues and lessons learned

Most of what went wrong was environment and methodology, not GROMACS itself, and each fix is now built into the scripts.

| Issue | Cause | Fix |
| --- | --- | --- |
| `-update gpu` refused on benchMEM | The 2013 topology has chains of 3+ coupled constraints, beyond what GPU LINCS supports | Dropped that config for benchMEM; added ADH (h-bonds only) to measure the GPU-resident mode |
| mdrun aborted: "PME tuning was still active" | Counter reset at step 2,500 came before PME grid tuning finished | Reset at step 10,000 of 20,000. Kept tuning on, because it is part of real GROMACS performance |
| 4 of 60 ADH runs 30% to 50% slower than their repeats | Interference from another process on the GPU host. The journal rules out the an unrelated translate sweep, whose runs don't overlap the slow ones. CPU contention from brief_worker, a host process of my RAG stack. The per-run CPU log of the isolated rerun caught it at 1100% CPU (11 of 12 cores) for about 20 s while GROMACS held 12 pinned threads; that run lost 9%, and nothing foreign touched the bench GPU. It wakes every 60 s but only bursts when it has work, which fits 4 disturbed runs out of 60. The 4 overnight points, rerun with both ollama containers stopped, reproduced within 1% | `analyze_bench.py` drops runs more than 5% below their group median and lists them. One slow run (nb_gpu) also chose a worse PME grid during tuning. Since 24 Oct the matrix runs under gpu-run.sh, which evicts ollama and pauses a competing host service; the worker itself is now capped to 4 CPUs and pinned to the V100 by UUID (GPU Sharing chapter). Per-run GPU and CPU process logs stay on. Still open, fix the SELinux denial |
| `GMXRC: shell: unbound variable` | The script's `set -u` collides with a variable GMXRC tests | `set +u; source GMXRC; set -u`. It only worked before because my shell already had GROMACS on the PATH |
| Build failed: cannot create `/root` | Build directory pointed at root's home while running as the ansible user | Build as a normal user in `$HOME/build`, install directory owned by ansible. Never compile as root |
| `deactivate: command not found` with an active venv | The venv was activated in a parent shell; variables are inherited, shell functions are not | Run builds under `env -i` with an explicit PATH instead of deactivating |
| Scripts not executable after transfer | Execute bit lost on copy | `chmod +x` on arrival |
| `wh_per_ns` meaningless for CPU-only runs | `nvidia-smi` measures only the GPU board, which idles at ~25 W | Added RAPL logging of CPU package power |

Three lessons I will carry into the next phases:

- **CPU utilization is not CPU work.** The GPU-resident ADH run still showed 1199% CPU, because GROMACS OpenMP threads busy-wait for the GPU. Only the RAPL power log shows real CPU load.
- **Legacy inputs can lock out modern features.** A published benchmark is good for comparison, but it can't exercise code paths its topology forbids.
- **Never edit a running bash script.** Bash reads scripts progressively from disk, so I make a copy for each variant (`bench_gpu_adh.sh`) instead.

## 1.7 Roadmap and current status

Phases 0 and 1 are done on the P100, on both benchmark systems. The V100 32GB order was cancelled by the seller, so the P100 vs V100 comparison moves back by about a month and phases 2 and 4 come first. Each later phase starts only after the previous one is benchmarked, and phase 4 is deliberately left open until the measured numbers decide its design.

| Phase | Content | Status |
| --- | --- | --- |
| 0 | Survey hosts, build GROMACS for the GPU host and the CPU nodes | Done. GPU build 106/106 tests passed; CPU build installed on NFS and verified on a compute node |
| 1a | P100 offload matrix on benchMEM | Done: 48 runs, CV below 1%, plus the fixed-grid control block |
| 1b | P100 offload matrix on ADH, including fully GPU-resident | Done: 60 runs with RAPL logging, 4 disturbed points rerun clean, fixed-grid control block. GPU-resident best at 118.4 ns/day |
| 1c | Reproducibility point under the card-B scheduler (`gpu-run.sh`): ADH nb_pme_gpu and GPU-resident at 12 threads, 3 runs each, within 1% of 105.2 and 118.4 | Done 24 Oct: nb_pme_gpu 105.20 (-0.02%), GPU-resident 120.16 (+1.5%, CV 0.2%). ollama evicted and re-admitted, a competing host service paused and restarted by gpu-run.sh |
| 2 | Same inputs on one compute node (`gmx`), then `gmx_mpi` on 1/2/4 nodes | Done 25 Oct, Slurm jobs 30 to 33. One VM = 96% of bare metal; multi-node gmx_mpi is slower than one node at every count (2 nodes 27% to 40% efficiency, 4 nodes 14% to 16%). Rule: one replica per node, thread-MPI |
| 3 | Swap P100 for V100 32GB, rerun phases 1a and 1b unchanged | Postponed: the seller cancelled the order, about one month delay. P100 stays card B |
| 4 | Science run: lysozyme (PDB 1AKI) stability at several temperatures | On the P100 under `gpu-run.sh`, GPU-resident on 4 threads (~120 ns/day on ADH-sized systems), plus four CPU replicas at one per node (13 to 16 ns/day each). Next |
| 5 | Analysis in Python with MDAnalysis (RMSD, RMSF, radius of gyration, H-bonds), first-person write-up on GitHub | Not started |

Infrastructure work that runs alongside:

- **Slurm GPU integration** is no longer a prerequisite. Card B is arbitrated by Kubernetes and Kueue (see the GPU Sharing chapter). If Slurm ever runs GPU jobs on gpu-host, its prolog and epilog call `gpu-run.sh`, so there is still one queue for the card.
- **CPU-only energy reference.** Covered by the ADH matrix (117 W package at 12 threads, 247 Wh/ns). A CPU-only benchMEM rerun with the logger is optional.

**A testable prediction for phase 3.** On benchMEM, the P100 averaged about 100 W of its 250 W TDP, because the CPU-side update step limits every config. I therefore expect the V100's advantage on benchMEM to be clearly smaller than the raw hardware difference, and ADH in GPU-resident mode to show the real gap between the two cards. ADH on the P100 already supports this: in resident mode the card draws 126 to 135 W instead of about 100 W.

## 1.8 Files and repository layout

All scripts and results live on the shared NFS volume; inputs and per-run outputs stay on the GPU host's local NVMe.

| Path | Host | Content |
| --- | --- | --- |
| `/nfs/shared/projects/gromacs/build_gromacs_gpu.sh` | GPU host | CUDA build: toolchain check, CMake venv, configure, `make check`, install |
| `/nfs/shared/projects/gromacs/build_gromacs_cpu.sh` | head | CPU builds `gmx` and `gmx_mpi`, installed to NFS for the offline nodes |
| `/nfs/shared/projects/gromacs/bench_gpu.sh` | GPU host | First benchMEM matrix (22 Oct) with the RAG-card guard; superseded by bench_matrix.sh, kept for the record |
| `/nfs/shared/projects/gromacs/bench_matrix.sh` | GPU host | Replaces `bench_gpu.sh`, `bench_gpu_adh.sh` and `bench_p100_extra.sh`: both systems, all offload levels, fixed-grid control block, per-run GPU and CPU process logs. Runs under `gpu-run.sh`; settings overridable from the environment |
| `/nfs/shared/projects/gromacs/analyze_bench.py` | head, analysis venv | Joins the CSV with the RAPL log, computes ns/day, CPU + GPU power, Wh/ns, writes the summary and chart |
| `/nfs/shared/gromacs/results/` | NFS | `bench_gpu.csv`, `bench_gpu_adh.csv`, `cpu_rapl_*.log`, summaries, charts |
| `/data/gromacs/bench/` | GPU host NVMe | `benchMEM.tpr`, `adh/ADH/adh_dodec/adh_dodec.tpr`, per-run directories with `md.log` and power samples |
| `/opt/gromacs/2026.2-cuda/` | GPU host | CUDA build of GROMACS |
| `/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/` | NFS | CPU builds of GROMACS for the compute nodes |
| /nfs/shared/projects/gromacs/slurm_node_sweep.sh, slurm_scaling.sh | head, sbatch | Phase 2: single-node rank x thread sweep (thread-MPI) and 1/2/4-node gmx_mpi scaling over the 10.10.10 fabric. Chain with -d afterany so jobs never overlap on the ESXi host |
| /nfs/shared/projects/gromacs/analyze_cpu.py | head, analysis venv | Phase 2 summaries (node sweep, scaling efficiency) and the scaling chart |

To reproduce the analysis on the head node:

```
source ~/venvs/gmx-analysis/bin/activate
cd /nfs/shared/gromacs/results
python /nfs/shared/projects/gromacs/analyze_bench.py
```

For a new data set, only the settings block at the top of `analyze_bench.py` changes: CSV, RAPL log, idle start time, output names and the chart title.
