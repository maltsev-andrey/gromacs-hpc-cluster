# Chapter 3. Phase 1: GPU offload benchmarks on the P100

*Andrey Maltsev · last updated 25 Oct 2026*

The P100 baseline on two standard systems, every offload level at 4 to 16 threads, with energy per simulated nanosecond. Method and timing rules are in section 1.4.

## 3.1 benchMEM offload matrix

The best configuration, nonbonded and PME on the GPU with 12 threads, reaches 111.4 ns/day. That is 6.5x the best CPU-only run (17.1 ns/day) on the same host.

![benchMEM on the P100: ns/day vs OpenMP threads per offload level](../images/p100_benchmem.png)

*Figure 3.1. benchMEM (82k atoms), GROMACS 2026.2, Tesla P100. Source: analyze_bench.py over bench_gpu.csv, 48 runs, 22 Oct 2026.*

Energy per simulated nanosecond, CPU package plus GPU board, at the most relevant points:

| Config | Threads | CPU package (W) | GPU board (W) | Energy (Wh/ns) |
| --- | --- | --- | --- | --- |
| + nonbonded on GPU | 12 | 92 | 93 | 56.5 |
| + PME on GPU | 8 | 63 | 101 | 37.5 |
| + PME on GPU | 12 | 75 | 102 | 38.2 |
| + bonded on GPU | 12 | 72 | 103 | 38.9 |

CPU-only energy is missing for benchMEM, because the RAPL log started after those runs. The ADH matrix below covers it.

What the data shows:

- **PME offload is the big step, not nonbonded alone.** At 4 threads it more than doubles throughput, from 39.9 to 83.9 ns/day, and it cuts energy per ns by a third. With PME on the CPU, the CPU was the bottleneck.
- **Bonded forces on the GPU make things slightly worse.** At 8 to 16 threads they cost about 3%. With nonbonded and PME on the GPU, the GPU is on the critical path, while the CPU has spare time and computes bonded forces in parallel for free.
- **Hyperthreading hurts in every config.** 16 threads run 3% to 14% slower than 12. This confirms my earlier MPI finding for OpenMP on the GPU host: stay at the physical core count.
- **The P100 is underutilized.** It averages about 100 W of its 250 W TDP. benchMEM forces update and constraints onto the CPU, so the GPU waits on the CPU every step.
- **The energy optimum is at 8 threads, not 12.** The last 4 cores add 12 W of CPU power for 6% more speed. The difference is small, but it shows that the fastest setting and the most efficient one are not the same.

## 3.2 ADH dodec, including fully GPU-resident mode

With update and constraints also on the GPU, ADH reaches 118.4 ns/day at 12 threads, 8.6x the best CPU-only run (13.8 ns/day) and 13% above the best partial-offload config.

![ADH on the P100: ns/day vs OpenMP threads, including GPU-resident mode](../images/p100_adh.png)

*Figure 3.2. ADH dodec (96k atoms), GROMACS 2026.2, Tesla P100. Source: analyze_bench.py over bench_gpu_adh.csv and the clean reruns of 23 Oct; 4 disturbed runs replaced.*

Energy per simulated nanosecond on ADH, CPU package plus GPU board, at 12 threads unless noted:

| Config | CPU package (W) | GPU board (W) | Energy (Wh/ns) | Less energy than CPU-only |
| --- | --- | --- | --- | --- |
| CPU only | 117 | 25 (idle) | 247 | 1.0x |
| + nonbonded on GPU | 93 | 97 | 66 | 3.7x |
| + PME on GPU | 82 | 114 | 45 | 5.5x |
| + bonded on GPU | 77 | 118 | 46 | 5.4x |
| Fully GPU-resident | 41 | 126 | 34 | 7.3x |
| Fully GPU-resident, 4 threads | 42 | 133 | 37 | 6.7x |

What ADH adds to the benchMEM picture:

- **In resident mode the thread count barely matters.** 4 threads already deliver 96% of the peak, because the CPU only launches kernels. For production I will run resident jobs on 4 threads, which leaves 8 cores free on the host.
- **The CPU package drops to about 41 W**, against 79 W with PME offload and 117 W CPU-only. GROMACS still reports about 1200% CPU in this mode, so the threads busy-wait cheaply. Only the power log shows this.
- **The GPU works harder.** The P100 draws 126 to 135 W in resident mode, against about 100 W in the best benchMEM config, because it no longer waits on the CPU each step.
- **Bonded forces on the GPU cost about 3% again**, the same pattern as benchMEM, so this finding holds across two systems.
- **GPU-resident mode uses about 7.3x less energy per simulated ns than CPU-only.** That is the energy reference point benchMEM could not provide.

## 3.3 Energy analysis

Energy per simulated nanosecond is (CPU package W + GPU board W) × 24 / (ns/day), in Wh/ns. CPU package power comes from the RAPL counter, GPU board power from nvidia-smi; DRAM, fans, NVMe and PSU losses are outside the scope (section 1.4). The tables in 3.1 and 3.2 carry the numbers; the ADH matrix also delivers the CPU-only reference that benchMEM could not, because its RAPL log started after the CPU-only runs.

The three findings that matter for production: PME offload cuts energy per ns by a third over nonbonded-only offload; the fully GPU-resident mode is both the fastest and the most efficient config, 7.3x less energy per ns than CPU-only; and the energy optimum on benchMEM sits at 8 threads, not 12, because the last four cores add 12 W for 6% more speed.

## 3.4 Fixed-grid control series

With the PME grid fixed (`-notunepme`), throughput only drops in the config where PME runs on the CPU: tuning is worth +52% on benchMEM and +85% on ADH there, and nothing once PME is on the GPU.

| System | Config, 12 threads | Tuned (ns/day) | Fixed grid (ns/day) | Tuning effect |
| --- | --- | --- | --- | --- |
| benchMEM | + nonbonded on GPU | 78.7 | 51.7 | +52% |
| benchMEM | + PME on GPU | 111.4 | 111.1 | none |
| ADH | + nonbonded on GPU | 68.8 | 37.2 | +85% |
| ADH | + PME on GPU | 105.2 | 105.3 | none |
| ADH | Fully GPU-resident | 118.4 | 119.9 | none (within noise) |

With PME on the CPU, the tuner raises the cutoff and coarsens the grid, moving work from the CPU's PME to the GPU's nonbonded kernel. This also explains the slowest overnight run: it reached 33.7 ns/day after keeping the default grid, close to the fixed-grid 37.2. For the P100 vs V100 comparison, the fixed-grid rows isolate the hardware, while the tuned rows show what a user actually gets; the nonbonded-only config is where the two differ most.

## 3.5 Disturbed runs: tracing four slow runs to a host process

Four of the 60 overnight ADH runs came in 30% to 50% below their own repeats, with identical input and settings. The analysis script dropped them as outliers, but I wanted the cause, because the same thing could hit the V100 runs.

**What the logs said.** The cycle accounting of a slow run against a clean repeat showed both the GPU wait and the CPU-side steps growing: "Wait GPU state copy" up 79%, the update step doubled. So whatever competed with GROMACS took CPU and delayed the GPU, which is the signature of a process that steals cores from pinned threads, not necessarily one that shares the card. One slow run had also kept a worse PME grid during tuning, which on its own costs a factor of two in the nonbonded-only config (see 3.4).

**What I ruled out with logs.** The Velero backup job (one container start at 22:02, no overlap with the slow runs), the vLLM RAG pod and both ollama containers (no requests in the window), the an unrelated translate sweep (each run under one second), and setroubleshootd (fires every minute, so it would have hit almost every run, not 4 of 60).

**How I caught it.** I reran the four points with both ollama containers stopped and a per-run logger recording every GPU process and the top CPU consumers every 5 seconds. All four points reproduced within 1%, and one more run in that series lost 9% with nothing foreign on the bench GPU. The CPU log for that minute shows it: `a competing host service.service`, a Python worker of my RAG stack that wakes every 60 s, went to 1100% CPU, 11 of the 12 cores, for about 20 seconds, while GROMACS held 12 pinned threads. The package power jumped from about 44 W to over 100 W at the same moment. It only bursts when it has work, which fits 4 disturbed runs out of 60.

**What changed.** The worker is now capped to 4 CPUs and pinned to the V100 by UUID, and every benchmark matrix runs under the card-B scheduler, which pauses it for the duration (chapter 5). Two lessons stayed with me: CPU utilization as GROMACS reports it is not CPU work, only the RAPL power log shows real load; and a host where several systems reach the GPUs and the CPU without knowing about each other will disturb measurements sooner or later, so the fix has to be structural, not a pre-flight check.

## 3.6 Reproducibility under the card-B scheduler

The same two points, run through `gpu-run.sh` with ollama evicted and a competing host service paused, reproduced the baseline: nb_pme_gpu at 105.20 ns/day (-0.02%) and GPU-resident at 120.16 ns/day (+1.5%), each the mean of 3 runs with CV below 0.4%. The resident mode came out faster than the overnight baseline because that baseline still carried the host's background noise (the worker's 60 s polling and setroubleshootd each minute); the fixed-grid resident runs on a quiet host had already shown 119.9. So the host noise floor before the scheduler was about 1.5% in the most GPU-bound config and invisible in the CPU-bound ones. For the V100 comparison the whole P100 matrix will be rerun once under `gpu-run.sh`, so both halves share the same regime.

