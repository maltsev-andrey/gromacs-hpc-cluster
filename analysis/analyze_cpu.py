#!/usr/bin/env python3
"""
Phase 2 analysis: GROMACS on the CPU compute nodes.
  bench_node_<host>.csv  single node, thread-MPI build, rank x thread sweep
  bench_scaling.csv      gmx_mpi on 1/2/4 nodes over 1 GbE
Outputs a summary table per file and one strong-scaling chart.
"""
# ---------- paths: resolved from the repo root, so this runs from anywhere ----------
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
DATA, RAW, IMG = REPO / "data", REPO / "data" / "raw", REPO / "images"

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- settings ----------
NODE_CSV = DATA / "bench_node_node02.csv"
SCALING_CSV = DATA / "bench_scaling.csv"
BARE_METAL = {"mem": 17.15, "adh": 13.82}     # GPU host CPU-only, 12 threads, same binary options
OUT_PNG = IMG / "cpu_scaling.png"
OUT_NODE = DATA / "summary_node_sweep.csv"
OUT_SCALING = DATA / "summary_scaling.csv"
# ------------------------------

pd.set_option("display.width", 160)

node = pd.read_csv(NODE_CSV)
ns = (node.groupby(["system", "ntmpi", "ntomp"], sort=False)["ns_per_day"]
          .agg(n="size", ns_day="mean", cv_pct=lambda x: 100 * x.std() / x.mean()).reset_index())
ns["vs_bare_metal_pct"] = [100 * r.ns_day / BARE_METAL[r.system] for r in ns.itertuples()]
ns.round(2).to_csv(OUT_NODE, index=False)
print("== single node, thread-MPI build\n", ns.round(2).to_string(index=False))

sc = pd.read_csv(SCALING_CSV)
ss = (sc.groupby(["system", "ranks_per_node", "ntomp", "nodes"], sort=False)
        .agg(n=("ns_per_day", "size"), ns_day=("ns_per_day", "mean"),
             cv_pct=("ns_per_day", lambda x: 100 * x.std() / x.mean()),
             npme=("npme", "first")).reset_index())
base = ss[ss.nodes == 1].set_index(["system", "ranks_per_node", "ntomp"])["ns_day"]
ss["speedup"] = [r.ns_day / base[(r.system, r.ranks_per_node, r.ntomp)] for r in ss.itertuples()]
ss["efficiency_pct"] = 100 * ss["speedup"] / ss["nodes"]
ss.round(2).to_csv(OUT_SCALING, index=False)
print("\n== gmx_mpi strong scaling\n", ss.round(2).to_string(index=False))

# ---------- chart: ns/day vs nodes, one line per system x layout ----------
fig, ax = plt.subplots(figsize=(8, 5.4), dpi=150)
styles = {("mem", 12): ("#4a7fc1", "-"), ("mem", 4): ("#4a7fc1", "--"),
          ("adh", 12): ("#d9822b", "-"), ("adh", 4): ("#d9822b", "--")}
names = {"mem": "benchMEM", "adh": "ADH"}
for (sysn, rpn, omp), grp in ss.groupby(["system", "ranks_per_node", "ntomp"], sort=False):
    c, ls = styles[(sysn, rpn)]
    ax.plot(grp.nodes, grp.ns_day, marker="o", color=c, ls=ls, lw=2,
            label=f"{names[sysn]}, {rpn} ranks x {omp} threads per node")
for sysn, c in (("mem", "#4a7fc1"), ("adh", "#d9822b")):
    best = ns[ns.system == sysn].ns_day.max()
    ax.axhline(best, color=c, lw=1, ls=":", alpha=0.7)
    ax.text(4.05, best, f"{names[sysn]} thread-MPI, 1 node", color=c, fontsize=8, va="center")
ax.set_xticks([1, 2, 4]); ax.set_xlim(0.8, 5.4); ax.set_ylim(0, 20)
ax.set_xlabel("compute nodes (12 vCPU each, 1 GbE)")
ax.set_ylabel("ns/day")
ax.set_title("Strong scaling of gmx_mpi over 1 GbE: every added node makes it slower", fontsize=10, loc="left")
ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.grid(axis="y", color="#eeeeee")
fig.tight_layout(); fig.savefig(OUT_PNG)
print(f"\nwrote {OUT_NODE}, {OUT_SCALING}, {OUT_PNG}")
