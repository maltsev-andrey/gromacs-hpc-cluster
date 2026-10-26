# Measurements

Every number in the top-level README comes from these files. One row per `gmx mdrun`,
written at the end of the run by the benchmark script.

| File | What |
| --- | --- |
| `bench_gpu.csv` | benchMEM offload matrix on the P100, 48 runs, 22 Oct 2026 |
| `bench_gpu_adh.csv` | ADH offload matrix on the P100, 60 runs, 22 Oct 2026 |
| `bench_p100_extra.csv` | clean reruns of 4 disturbed ADH points + the fixed-grid control block, 23 Oct |
| `bench_node_node02.csv` | single compute node, rank × thread sweep, 36 runs, 25 Oct |
| `bench_scaling.csv` | `gmx_mpi` on 1, 2 and 4 nodes, 36 runs, 25 Oct |
| `summary_*.csv` | the aggregated tables the analysis scripts produce |
| `raw/cpu_rapl_*.log` | CPU package energy counter, sampled every 5 s, joined to runs by timestamp |

## Columns

GPU matrices: `date, gpu, config, ntomp, repeat, ns_per_day, avg_power_w, wh_per_ns[, start_epoch]`
where `avg_power_w` is the GPU board only. CPU package power is added during analysis
from the RAPL logs, which is why `wh_per_ns` in the raw CSV understates total energy and
is recomputed in the summaries.

Node sweep: `date, node, system, ntmpi, ntomp, repeat, ns_per_day, wall_s, npme, start_epoch`
Scaling: `date, nodes, system, ranks_per_node, ntomp, repeat, ns_per_day, wall_s, npme, start_epoch`

`ns_per_day` covers only the second half of each run; counters are reset at step 10,000 of
20,000 so that startup and PME tuning are excluded.

## Known-bad rows

`bench_gpu_adh.csv` contains the four runs disturbed by host CPU contention, left in
deliberately. `analyze_bench.py` drops any run more than 5% below its group median and
prints what it dropped; `bench_p100_extra.csv` supersedes those points.
