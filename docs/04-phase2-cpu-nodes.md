# Chapter 4. Phase 2: CPU compute nodes

*Andrey Maltsev · last updated 25 Oct 2026*

The same two systems on the Slurm compute nodes: first one VM on its own, then gmx_mpi across 1, 2 and 4 nodes over the 1 GbE fabric. Scripts: slurm_node_sweep.sh and slurm_scaling.sh, chained with `-d afterany` so no two jobs share the ESXi host. No energy column here, because the VMs expose no RAPL counter. The result in one line: one VM is nearly as good as bare metal, and adding nodes over 1 GbE makes every run slower than one node alone. Section 4.3 records a confound found afterwards in the 4-node measurements.

## 4.1 Single-node baseline (thread-MPI)

One compute VM delivers 96% to 97% of the bare-metal GPU host in CPU-only mode.

Sweeping how the 12 vCPUs split into ranks and threads changes little: benchMEM lands at 16.5 ns/day with either 1 rank x 12 threads or 12 ranks x 1 thread, and the mixed layouts cost 7% to 11%. ADH prefers pure OpenMP (13.5) and loses 8% at 12 ranks. Against the GPU host's 12 Haswell cores (17.15 and 13.82) the VM gives 96% and 97%, so ESXi costs almost nothing for this workload, and the vCPUs behave like full cores.

| System | Best layout | ns/day | GPU host CPU-only, 12 cores | VM / bare metal |
| --- | --- | --- | --- | --- |
| benchMEM | 1 rank x 12 threads, or 12 x 1 | 16.5 | 17.15 | 96% |
| ADH | 1 rank x 12 threads | 13.5 | 13.82 | 97% |

Repeat variation stayed below 2% at every point, so the node was quiet. ESXi presents the 12 vCPUs as 12 single-core sockets with no hyperthreading inside the VM.

## 4.2 Multi-node scaling with gmx_mpi on 1, 2 and 4 nodes

The OpenMPI build on one node is already 10% to 15% slower than thread-MPI (14.9 vs 16.5 on benchMEM). From there it only goes down:

![gmx_mpi strong scaling over 1 GbE: every added node makes it slower](../images/cpu_scaling.png)

*Figure 4.1. Strong scaling of gmx_mpi on 1, 2 and 4 compute nodes; dotted lines are the thread-MPI build on one node. Source: analyze_cpu.py over bench_scaling.csv and bench_node_node02.csv, Slurm jobs 30 to 33, 25 Oct 2026.*

| System, layout per node | 1 node | 2 nodes (efficiency) | 4 nodes (efficiency, see 4.3) |
| --- | --- | --- | --- |
| benchMEM, 12 ranks x 1 thread | 14.9 | 11.8 (40%) | 8.2 (14%) |
| benchMEM, 4 ranks x 3 threads | 14.0 | 7.5 (27%) | 9.2 (16%) |
| ADH, 12 ranks x 1 thread | 11.5 | 7.6 (33%) | 7.3 (16%) |
| ADH, 4 ranks x 3 threads | 12.4 | 6.8 (28%) | 7.5 (15%) |

### A confound in the 4-node column

I found this after the fact, while sizing the phase 4 production runs. All four compute VMs live on one ESXi host with two Xeon E5-2696 v3, so **36 physical cores and 72 logical threads**. Each job's vCPU demand therefore depends on its node count:

| Test | vCPUs in use | Against 36 physical cores | Clean measurement? |
| --- | --- | --- | --- |
| 1 node | 12 | plenty spare | yes |
| 2 nodes | 24 | still spare | yes |
| 4 nodes | 48 | 1.33x oversubscribed | no |

At 4 nodes the vCPUs still fit inside the 72 logical threads, so ESXi never queues one, but roughly half of them share a physical core with an HT sibling. For floating-point-heavy code an HT pair delivers about 1.25x a single thread rather than 2x, so a shared vCPU runs near 60% of a dedicated one. **The 14 to 16% efficiency figures therefore mix the network limit with CPU oversubscription and overstate the network penalty.**

The conclusion survives, because the 2-node case carries it on its own: 24 vCPUs against 36 cores with nothing to contend for, and benchMEM still falls from 14.9 to 11.8 ns/day, 40% efficiency. That is the network, measured cleanly. The 4-node row shows the direction of travel but should not be quoted as a pure scaling number.

## 4.3 Where PME breaks on 1 GbE

The evidence for the network being the wall comes from the logs rather than from the timings, which is why the oversubscription above does not undermine it:

- **PME is the wall.** On 4 nodes mdrun assigns 12 of 48 ranks to PME, and the PME grid work is an all-to-all exchange every step. The tuner reacts by pushing the cut-off from 0.9 to 1.226 nm, which moves work from PME to the short-range kernel; mdrun itself notes that this raised the nonbonded workload by more than 50%. It is trading compute for communication, and still loses.
- **Imbalance on top.** The 4-node ADH log reports 46% average load imbalance, 26% of the run waiting on it, plus 9.5% lost to PP/PME imbalance. Over TCP with ~112 us latency, ranks spend more time waiting for each other than computing.
- **2 nodes is noisier than 4.** The 2-node points vary by 6% to 15% between repeats, the 4-node points by 2% to 8%. Two of the three benchMEM 12x1 runs on 2 nodes failed outright (cause still to be read from mdrun.out).
- **Fewer, larger messages don't rescue it.** The 4x3 layout, chosen because message count beat message size in my earlier MPI work, is no better here: PME's all-to-all volume scales with the grid, not with the rank count.

This is the same fabric limit I measured in my earlier MPI work (about 112 us latency, message count dominating), seen from the other side: there, fewer and larger messages won; here PME's all-to-all volume is set by the grid, so the rank layout cannot buy it back.

To separate the two effects properly, the 4-node test would need to run with 6 or 8 threads per node so that total vCPUs stay under 36. That is on the list for when the cluster is next idle; it would not change the recommendation, only the precision of the 4-node figure.

## 4.4 Production rule: ensembles, not decomposition

The practical rule for this cluster follows directly: **one simulation per node with the thread-MPI build**, four independent replicas across the four nodes. That gives 4 x 16.5 = 66 ns/day of aggregate benchMEM throughput, against 9.2 for the best 4-node run of a single simulation. Ensembles, not decomposition.

Two loose ends for the record: two of the three benchMEM 12x1 runs on 2 nodes failed outright (cause to be read from mdrun.out), and the 2-node points vary by 6% to 15% between repeats where the 4-node points vary by 2% to 8%, which looks like the 2-node case sitting right where the tuner and dynamic load balancing fight the network.
