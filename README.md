# GROMACS molecular dynamics on a home HPC cluster

Measuring what each layer of hardware actually buys you in molecular dynamics, in ns/day **and** in energy per simulated nanosecond, on a 5-node RHEL cluster with one usable Tesla P100. Every number here is measured on that hardware, including the ones that contradicted my expectations.

**The short version:** moving the whole MD step onto the GPU gives 8.6x the throughput of the same host's CPUs and uses 7.3x less energy per nanosecond. Spreading one simulation across multiple nodes over 1 GbE makes it *slower* than a single node, at every node count I tried. The practical consequence is ensembles, not domain decomposition. That rule then produced the science: five 200 ns lysozyme replicas from 300 K to 450 K, one microsecond in total, showing a fold that holds at every temperature while a quarter of the chain comes loose.

---

## Hardware

| Host | Role | Spec | Network |
|---|---|---|---|
| gpu-host | GPU host, bare metal | Xeon E5-2678 v3 (12 cores), 128 GB, Tesla P100 16GB + Tesla V100 16GB | home LAN only |
| head-node | head, NFS, Slurm controller | ESXi VM, 32 GB | dual-homed |
| node02..05 | compute nodes | ESXi VMs, 12 vCPU, 7.6 GB each | 1 GbE, ~112 µs latency |

All VMs share one ESXi host: 2 × Xeon E5-2696 v3, 36 physical cores, 72 threads, 256 GB. That matters for the multi-node results below.

The V100 is reserved for a RAG assistant that runs permanently, so all GPU work here uses the P100. GROMACS 2026.2, CUDA 12.4, driver 550, AVX2_256, pinned for the whole study.

Two benchmark systems: **benchMEM** (81,743 atoms, the MPI-NAT standard) and **ADH dodec** (95,561 atoms, from the GROMACS developers' own set, rebuilt with `gmx grompp` 2026.2).

---

## Result 1: offload the whole step, not part of it

![benchMEM offload scaling on the Tesla P100](images/p100_benchmem.png)

![ADH offload scaling including GPU-resident mode](images/p100_adh.png)

| Config (ADH, 12 threads) | ns/day | vs CPU-only | CPU pkg | GPU board | Wh per ns |
|---|---|---|---|---|---|
| CPU only | 13.8 | 1.0x | 117 W | 25 W idle | 247 |
| + nonbonded on GPU | 68.8 | 5.0x | 93 W | 97 W | 66 |
| + PME on GPU | 105.2 | 7.6x | 82 W | 114 W | 45 |
| + bonded on GPU | 101.5 | 7.3x | 77 W | 118 W | 46 |
| **fully GPU-resident** | **118.4** | **8.6x** | **41 W** | 126 W | **34** |

Three things worth pulling out:

- **PME offload is the big step, not nonbonded.** At 4 threads, moving PME to the GPU more than doubles throughput (39.9 → 83.9 ns/day on benchMEM). Offloading nonbonded alone leaves the CPU doing the FFT work and bottlenecking everything.
- **Offloading bonded forces makes it slightly slower.** About 3%, consistently, on both systems. Once nonbonded and PME are on the GPU, the GPU is on the critical path and the CPU has idle time, so bonded forces on the CPU are free. Moving them just lengthens the GPU queue. More offload is not automatically faster.
- **In GPU-resident mode the thread count barely matters.** 4 threads give 96% of the peak, because the CPU only launches kernels. CPU package power drops from 117 W to 41 W. That frees 8 cores on the host for other work, which turned out to matter a lot (see "the measurement that wasn't" below).

Energy is CPU package (RAPL) plus GPU board (nvidia-smi), sampled during every run. No DRAM domain exists on this board, and fans, NVMe and PSU losses are out of scope.

---

## Result 2: multi-node MD over 1 GbE doesn't work

![gmx_mpi strong scaling across 1, 2 and 4 nodes](images/cpu_scaling.png)

| System, layout per node | 1 node | 2 nodes | 4 nodes |
|---|---|---|---|
| benchMEM, 12 ranks × 1 thread | 14.9 | 11.8 (40% eff.) | 8.2 (14%) |
| benchMEM, 4 ranks × 3 threads | 14.0 | 7.5 (27%) | 9.2 (16%) |
| ADH, 12 ranks × 1 thread | 11.5 | 7.6 (33%) | 7.3 (16%) |
| ADH, 4 ranks × 3 threads | 12.4 | 6.8 (28%) | 7.5 (15%) |

No configuration on 2 or 4 nodes beats a single node. The logs say exactly why. On 4 nodes mdrun assigns 12 of 48 ranks to PME, and the PME grid exchange is an all-to-all every single step. The tuner responds by pushing the cutoff from 0.9 to 1.226 nm to move work away from PME, and notes in its own output that this *raised the nonbonded workload by more than 50%*. It is trading compute for communication and still losing: 46% average load imbalance, 26% of the run spent waiting on it, 9.5% more lost to PP/PME imbalance.

**A confound I found afterwards, which belongs here rather than in a footnote.** All four VMs sit on one ESXi host with 36 physical cores. At 4 nodes the test wants 48 vCPUs, so about half of them share a physical core with an HT sibling and run at roughly 60% of a dedicated core. The 4-node efficiency figures therefore mix the network limit with CPU oversubscription and overstate the network penalty. The 2-node case is clean (24 vCPUs against 36 cores) and still shows 40% efficiency, so the conclusion stands on that row plus the log evidence above. One caveat on that cell: only one of its three repeats survived, so the 40% carries no spread. The three 2-node configurations that do have full repeats land at 27%, 33% and 28%, all well under a single node, and that is the result carrying the argument. Re-running the 4-node test at 6 threads per node would isolate it properly; it is on the list, and it would sharpen the number without changing the recommendation.

A detail I expected to help and didn't: using fewer, fatter ranks (4×3 instead of 12×1) made no difference. In earlier MPI work on this same fabric, message *count* dominated and that trick won. Here it doesn't, because PME's all-to-all volume is set by the grid, not the rank count.

**So the production rule for this cluster is ensembles, not decomposition.** Four independent replicas, one per node, give 4 × 16.5 = 66 ns/day of aggregate benchMEM throughput. The best single 4-node run gives 9.2.

One more measurement worth having: a compute VM delivers **96 to 97%** of the bare-metal host's CPU-only performance. ESXi costs essentially nothing for this workload.

---

## Result 3: PME tuning is worth up to 85%, or nothing

GROMACS retunes the PME grid and cutoff at the start of each run. Pinning them with `-notunepme` shows what that is worth:

| System | Config, 12 threads | Tuned | Fixed grid | Tuning effect |
|---|---|---|---|---|
| benchMEM | + nonbonded on GPU | 78.7 | 51.7 | **+52%** |
| benchMEM | + PME on GPU | 111.4 | 111.1 | none |
| ADH | + nonbonded on GPU | 68.8 | 37.2 | **+85%** |
| ADH | + PME on GPU | 105.2 | 105.3 | none |
| ADH | fully GPU-resident | 118.4 | 119.9 | none |

When PME runs on the CPU, the tuner saves you from a factor of two. When PME is on the GPU, the default grid is already right. This also explained an outlier: one run that landed at 33.7 ns/day had simply kept the default grid during tuning, landing close to the fixed-grid 37.2.

---

## Result 4: the fold holds to 450 K, but not uniformly

![per-residue C-alpha fluctuation at five temperatures](images/lyso_rmsf.png)

Five 200 ns replicas of lysozyme (PDB 1AKI) at 300, 325, 350, 400 and 450 K, one per temperature, built on the ensemble rule from Result 2.

| Measure | 300 K | 450 K | against temperature |
|---|---|---|---|
| RMSD | 0.085 nm | 0.259 nm | linear, r² = 0.97 |
| Radius of gyration | 1.415 nm | 1.399 nm | flat, 2% spread across the ladder |
| Backbone H-bonds | 43.3 | 34.8 | linear, r² = 0.98 |
| Rigid core, median RMSF | 0.052 nm | 0.090 nm | +73% |
| Mobile regions, median RMSF | 0.072 nm | 0.341 nm | +376% |
| Residues below 0.15 nm | 98% | 75% | |

**Nothing unfolds.** Rg never leaves 1.39 to 1.42 nm; unfolding would carry it past 2 nm. RMSD and H-bond loss are both linear across 150 K with no cooperative step, which is progressive loosening rather than melting.

**The loosening is local.** At 450 K three quarters of the chain is still as rigid as a folded protein, while the C-terminal tail and the loop at residues 101 to 107 take nearly all of the extra motion. Residue 129 reaches 0.726 nm, eight times the core median.

A hypothesis that did not survive: lysozyme's four disulfide bridges do not predict this. Spearman correlation between sequence distance to the nearest bridged cysteine and RMSF gain is +0.046, essentially zero, and the longest unanchored stretch is among the least mobile.

Caveats that matter: production was NVT at fixed volume, so the 450 K system sits under pressure that opposes unfolding; 200 ns is short for denaturation; and one replica per temperature is not enough for the loop statistics or for the 325 to 350 K step, which happens to sit next to lysozyme's experimental melting point. [Chapter 5](docs/05-phase5-thermal-stability.md) has the detail, including three analysis bugs that each produced confident wrong numbers first.

---

## The measurement that wasn't

Four of 60 runs in the overnight ADH matrix came in 30 to 50% below their own repeats. Same input, same settings, same card.

The cycle accounting showed both the GPU wait *and* the CPU-side steps growing, which rules out a simple "something else used the GPU" story. Logs cleared the obvious suspects one by one: a Velero backup job (no time overlap), the vLLM RAG pod and two ollama containers (no requests in the window), a translation sweep (each run under a second), and `setroubleshootd` (fires every minute, so it would have hit nearly every run, not 4 of 60).

Rerunning the four points with per-run process logging caught it: a Python worker from an unrelated service on the same host bursting to **1100% CPU, 11 of 12 cores, for about 20 seconds**, while GROMACS held 12 pinned threads. It wakes every 60 s but only does work sometimes, which fits 4 disturbed runs out of 60 exactly.

The fix was structural rather than a pre-flight check: that host now has a single arbiter for its GPU, and benchmark runs claim the card through a queue that evicts competing services for the duration. Reruns under that regime reproduced the baseline within 1%, and showed the host's noise floor beforehand had been about 1.5% in the most GPU-bound config.

Two lessons I'd carry anywhere:

- **CPU utilization is not CPU work.** GROMACS reports ~1200% CPU in GPU-resident mode, but the RAPL power log shows 41 W. The threads busy-wait cheaply. Only power told the truth.
- **Legacy inputs can lock out modern code paths.** benchMEM's 2013 topology has chains of 3+ coupled constraints, beyond what GPU LINCS supports, so it simply refuses `-update gpu`. A published benchmark is great for comparability and useless for exercising features its topology forbids. That's why ADH is here too.

---

## Documentation

The README is the summary. The chapters are the working record: method, every number, and what went wrong on the way.

| Chapter | What is in it |
|---|---|
| [1. Project overview](docs/01-project-overview.md) | The cluster, the constraints it imposes, the software stack, benchmark methodology, a script index, results to date, and the running list of issues and lessons |
| [2. Build and environment](docs/02-build-and-environment.md) | Two builds of GROMACS 2026.2 from one source: CUDA on the GPU host, CPU on NFS for the offline nodes. Host survey, SIMD and toolchain choices, 106/106 regression tests, and the analysis venv |
| [3. Phase 1: GPU offload benchmarks](docs/03-phase1-gpu-benchmarks.md) | The P100 baseline on benchMEM and ADH at every offload level and 4 to 16 threads, energy per simulated nanosecond from RAPL plus board power, the fixed-grid control series, and how four disturbed runs were traced to their cause |
| [4. Phase 2: CPU compute nodes](docs/04-phase2-cpu-nodes.md) | One VM against bare metal, then gmx_mpi across 1, 2 and 4 nodes over 1 GbE. Where PME breaks on this fabric, the ESXi oversubscription confound found afterwards, and the ensemble rule that follows |
| [5. Phase 5: thermal stability of lysozyme](docs/05-phase5-thermal-stability.md) | Five 200 ns replicas from 300 to 450 K, 1 µs in total. RMSD, per-residue RMSF, radius of gyration and backbone H-bonds; why the loosening is local rather than global; and three analysis bugs that each produced confident wrong numbers first |

Chapter 6 will cover the replicate run at 350 K.

---

## Reproducing it

```bash
# 1. build GROMACS (GPU host, then head node for the compute nodes)
scripts/build/build_gromacs_gpu.sh          # CUDA 12.4, sm_60+sm_70, AVX2_256
scripts/build/build_gromacs_cpu.sh          # gmx (thread-MPI) + gmx_mpi, installed to NFS

# 2. GPU offload matrix, both systems, all offload levels, fixed-grid control block
scripts/bench/bench_matrix.sh

# 3. CPU nodes: single-node rank×thread sweep, then 1/2/4-node scaling
sbatch scripts/bench/slurm_node_sweep.sh
for n in 1 2 4; do sbatch -N $n scripts/bench/slurm_scaling.sh; done

# 4. analysis: merges CSVs with the RAPL logs, writes summaries and charts
python analysis/analyze_bench.py            # GPU matrices
python analysis/analyze_cpu.py              # CPU nodes

# 5. phase 4+5: lysozyme production, then structural analysis
scripts/lysozyme/prep_lysozyme.sh           # PDB 1AKI -> solvated, equilibrated, one tpr per T
for T in 300 325 350 400; do sbatch -J lyso-T$T --export=T=$T scripts/lysozyme/slurm_lyso_node.sh; done
T=450 HOURS=8 scripts/lysozyme/run_lyso_gpu.sh      # GPU replica, nightly windows
scripts/lysozyme/prep_analysis.sh           # protein-only topology + shared reference
scripts/lysozyme/make_whole.sh              # gmx trjconv -pbc whole, required before analysis
python analysis/analyze_lyso.py             # RMSD, RMSF, Rg, backbone H-bonds
```

The two Python scripts need `pandas`, `matplotlib` and `numpy`, and Python 3.9 or newer:

```bash
python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python analysis/analyze_bench.py
```

They resolve `data/` and `images/` from the repository root, so they run from any working directory. Every script keeps its settings in a block at the top; nothing takes command-line arguments. `analyze_bench.py` merges several runs, matches each to its exact time window in the power log, lets clean reruns supersede disturbed points, and drops anything more than 5% below its group median while listing what it dropped.

All raw CSVs and power logs are in `data/`, so the charts and tables above can be regenerated without the cluster.

---

## Layout

```
scripts/build/      GROMACS builds: CUDA on the GPU host, CPU on NFS for offline nodes
scripts/bench/      offload matrix + Slurm sweep and scaling jobs
scripts/lysozyme/   phases 4-5: system prep from PDB 1AKI, production runners, analysis prep
analysis/           CSV + RAPL joiners, summary tables, charts
data/               every measurement, raw power logs under data/raw/, phase 5 under data/phase5/
docs/               the long-form write-up, five chapters
images/             the charts used above
```

---

## Status

| Phase | What | State |
|---|---|---|
| 0 | Build and verify on both targets | done, 106/106 regression tests |
| 1 | P100 offload matrices, both systems, energy, fixed-grid control | done, 108 runs |
| 2 | CPU nodes: single-node sweep and 1/2/4-node scaling | done, 36 + 36 runs |
| 3 | P100 vs V100 32GB comparison | postponed, card order cancelled |
| 4 | Lysozyme (PDB 1AKI) thermal stability, 5 temperatures × 200 ns | done, 1 µs total |
| 5 | MDAnalysis: RMSD, RMSF, radius of gyration, H-bonds | done |
| 6 | Three replicas at 350 K, to test the step next to the melting point | next |

Phases 4 and 5 completed 2 Nov 2026: one GPU replica at 450 K in GPU-resident mode on 4 threads at 339 ns/day, four CPU replicas at 300/325/350/400 K at about 27 ns/day each, one per node. The layout came straight out of phases 1 and 2 rather than from a guess, and 100,005 frames of trajectory came back from it.

## License

MIT, see [LICENSE](LICENSE). Benchmark inputs belong to their original authors: benchMEM from the [MPI-NAT benchmark set](https://www.mpinat.mpg.de/grubmueller/bench), ADH from the [GROMACS developers' set](https://zenodo.org/records/11234002).
