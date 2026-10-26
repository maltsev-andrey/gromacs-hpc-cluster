#!/usr/bin/env python3
"""
GROMACS offload benchmark analysis: ns/day + energy per simulated ns.

Inputs
  bench_gpu.csv       one row per mdrun, written at run END by bench_gpu.sh
  cpu_rapl_*.log      unix_ts,pkg_energy_uj[,dram_uj]  sampled every 5 s

Energy scope: CPU package (RAPL) + GPU board (nvidia-smi). No DRAM, fans, PSU.
Run window: [end of previous row, end of this row]. Includes mdrun startup and
PME tuning, same window the GPU power average in the CSV covers.
"""
import sys
import numpy as np
# ---------- paths: resolved from the repo root, so this runs from anywhere ----------
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
DATA, RAW, IMG = REPO / "data", REPO / "data" / "raw", REPO / "images"

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- settings ----------
# Each CSV: (path, default system for rows without a "sys:" prefix, supersede).
# supersede=True: rows from this file replace all earlier rows of the same
# (system, config, threads) point, e.g. clean reruns of disturbed points.
CSVS = [
    (DATA / "bench_gpu_adh.csv", "adh", False),
    (DATA / "bench_p100_extra.csv", None, True),
]
SYSTEM = "adh"                            # which system to summarize
RAPL_LOGS = [RAW / "cpu_rapl_p100_adh.log", RAW / "cpu_rapl_p100_extra.log"]   # oldest first
IDLE_FROM = "2026-10-22T22:36:30+02:00"   # idle window start, in the log that holds it
IDLE_TO = "2026-10-22T22:45:09+02:00"
OUT_SUMMARY = DATA / "summary_p100_adh.csv"
OUT_PNG = IMG / "p100_adh.png"
TITLE = "ADH dodec (96k atoms), GROMACS 2026.2, Tesla P100 + Xeon E5-2678 v3"
RAPL_WRAP_UJ = 262143328850               # max_energy_range_uj, typical Haswell value
OUTLIER_DROP = 0.05                       # drop a run >5% below its group median (clean CV is <1%)
# ------------------------------

frames = []
for order, (path, default_sys, supersede) in enumerate(CSVS):
    f = pd.read_csv(path, parse_dates=["date"])
    f["end"] = f["date"].map(lambda t: t.timestamp())
    if "start_epoch" in f.columns:
        f["start"] = f["start_epoch"].astype(float)
    else:
        f["start"] = f["end"].shift(1)    # first row has no start -> no CPU energy
    split = f["config"].str.split(":", n=1, expand=True)
    if split.shape[1] == 2:
        f["system"] = split[0].where(split[1].notna(), default_sys)
        f["config"] = split[1].where(split[1].notna(), split[0])
    else:
        f["system"] = default_sys
    f["order"] = order
    f["supersede"] = supersede
    frames.append(f)
runs = pd.concat(frames, ignore_index=True)
runs = runs[runs["system"] == SYSTEM].copy()

# superseding files replace earlier rows of the same point
key = ["system", "config", "ntomp"]
for order in sorted(runs.loc[runs["supersede"], "order"].unique()):
    newer = runs[runs["order"] == order][key].drop_duplicates()
    m = runs.merge(newer, on=key, how="left", indicator=True)["_merge"].eq("both").to_numpy()
    runs = runs[~(m & (runs["order"] < order).to_numpy())].copy()

logs = [pd.read_csv(p, header=None, names=["ts", "pkg_uj", "dram_uj"]) for p in RAPL_LOGS]
r = pd.concat(logs, ignore_index=True).sort_values("ts")
e = r["pkg_uj"].astype(float).to_numpy()
# unwrap counter overflow (gaps between logs are never interpolated across a run)
wraps = np.cumsum(np.r_[0, (np.diff(e) < 0).astype(int)])
e = e + wraps * RAPL_WRAP_UJ
ts = r["ts"].to_numpy(float)

def pkg_watts(t0, t1):
    """Mean package power over [t0, t1] by interpolating the cumulative counter."""
    if np.isnan(t0) or t0 < ts[0] or t1 > ts[-1]:
        return np.nan
    # reject windows that span a gap between logs
    inside = ts[(ts >= t0) & (ts <= t1)]
    edges = np.r_[t0, inside, t1]
    if np.max(np.diff(edges)) > 30:
        return np.nan
    return (np.interp(t1, ts, e) - np.interp(t0, ts, e)) / 1e6 / (t1 - t0)

runs["cpu_pkg_w"] = [pkg_watts(a, b) for a, b in zip(runs["start"], runs["end"])]

idle_t0 = pd.Timestamp(IDLE_FROM).timestamp()
idle_w = pkg_watts(idle_t0, pd.Timestamp(IDLE_TO).timestamp())

# ---------- outliers: interference from other processes on the host ----------
med = runs.groupby(["config", "ntomp"])["ns_per_day"].transform("median")
runs["group_median"] = med
runs["outlier"] = runs["ns_per_day"] < (1 - OUTLIER_DROP) * med
dropped = runs[runs["outlier"]]
runs = runs[~runs["outlier"]].copy()

runs["total_w"] = runs["cpu_pkg_w"] + runs["avg_power_w"]
runs["wh_per_ns_total"] = runs["total_w"] * 24 / runs["ns_per_day"]

g = runs.groupby(["config", "ntomp"], sort=False)
summary = g.agg(
    n=("ns_per_day", "size"),
    ns_day=("ns_per_day", "mean"),
    ns_day_cv_pct=("ns_per_day", lambda x: 100 * x.std() / x.mean()),
    gpu_w=("avg_power_w", "mean"),
    cpu_pkg_w=("cpu_pkg_w", "mean"),
    total_w=("total_w", "mean"),
    wh_per_ns=("wh_per_ns_total", "mean"),
).reset_index()

cfg_order = ["cpu_only", "nb_gpu", "nb_pme_gpu", "nb_pme_bonded_gpu", "full_gpu_resident"]
summary["_o"] = summary["config"].map(lambda c: cfg_order.index(c.replace("_notune", "")) + (10 if c.endswith("_notune") else 0)
                                      if c.replace("_notune", "") in cfg_order else 99)
summary = summary.sort_values(["_o", "ntomp"]).drop(columns="_o").reset_index(drop=True)

best_cpu = summary.loc[summary.config == "cpu_only", "ns_day"].max()
summary["speedup_vs_best_cpu"] = summary["ns_day"] / best_cpu
summary.round(2).to_csv(OUT_SUMMARY, index=False)

pd.set_option("display.width", 160)
print(summary.round(2).to_string(index=False))
if len(dropped):
    print(f"\nDropped {len(dropped)} outlier run(s), >{OUTLIER_DROP:.0%} below their group median:")
    for _, d in dropped.iterrows():
        print(f"  {d['date']:%H:%M:%S}  {d['config']:<20} nt={d['ntomp']:<3} r{d['repeat']}  "
              f"{d['ns_per_day']:.1f} ns/day vs median {d['group_median']:.1f}")
print(f"\nIdle CPU package baseline: {idle_w:.1f} W")
print(f"RAPL log covers {pd.Timestamp(ts[0], unit='s', tz='Europe/Zurich')} .. "
      f"{pd.Timestamp(ts[-1], unit='s', tz='Europe/Zurich')}")

# ---------- chart: ns/day vs OpenMP threads per offload level ----------
labels = {
    "cpu_only": "CPU only",
    "nb_gpu": "+ nonbonded on GPU",
    "nb_pme_gpu": "+ PME on GPU",
    "nb_pme_bonded_gpu": "+ bonded on GPU",
    "full_gpu_resident": "+ update on GPU (resident)",
}
# fixed-grid control configs (-notunepme) stay in the summary but out of the chart
colors = {"cpu_only": "#8a8f98", "nb_gpu": "#4a7fc1",
          "nb_pme_gpu": "#d9822b", "nb_pme_bonded_gpu": "#2e9e6b",
          "full_gpu_resident": "#8e5bb5"}

fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
for cfg in labels:
    s = summary[summary.config == cfg]
    if s.empty:
        continue
    ax.plot(s.ntomp, s.ns_day, marker="o", lw=2, color=colors[cfg])
    last = s.iloc[-1]
    dy = {"nb_pme_gpu": 7, "nb_pme_bonded_gpu": -7}.get(cfg, 0)
    ax.annotate(labels[cfg], (last.ntomp, last.ns_day), xytext=(8, dy),
                textcoords="offset points", va="center", color=colors[cfg], fontsize=9)
ax.axvline(12, color="#bbbbbb", ls="--", lw=1)
ax.text(12.15, 4, "12 physical cores", color="#888888", fontsize=8)
ax.set_xticks([4, 8, 12, 16])
ax.set_xlim(3, 21.5)
ax.set_ylim(0, None)
ax.set_xlabel("OpenMP threads (-ntomp)")
ax.set_ylabel("ns/day")
ax.set_title(TITLE, fontsize=10, loc="left")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.grid(axis="y", color="#eeeeee")
fig.tight_layout()
fig.savefig(OUT_PNG)
print(f"\nwrote {OUT_SUMMARY}, {OUT_PNG}")
