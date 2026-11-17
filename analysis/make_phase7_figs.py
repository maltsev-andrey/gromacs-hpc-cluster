#!/usr/bin/env python3
"""Figures for Chapter 7, built from the CSVs that analyze_lyso.py and
analyze_core.py write. No arguments: edit the settings block and run it.

  phase7_extension.png    whole-protein C-alpha RMSD of the four extended runs
  phase7_loop_states.png  which form the loop 100-117 is in, per run, over time
  phase7_core_ladder.png  three RMSD ladders against temperature
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ----------------------------------------------------------------- settings
ANALYSIS = Path("/nfs/shared/gromacs/lyso/analysis")   # where the CSVs are
OUT      = ANALYSIS                                    # where the PNGs go
EXTENDED = [("T300", "300 K, extended to 600 ns"),
            ("T325", "325 K replica 1, extended to 400 ns"),
            ("T400", "400 K, extended to 400 ns"),
            ("T450", "450 K, extended to 400 ns")]
# whole-protein means over 5-200 ns as printed in Chapter 5 (one run per temperature)
CH5 = {300: 0.0845, 325: 0.0908, 350: 0.1529, 400: 0.1913, 450: 0.2590}
PUBLISHED_NS, EQUIL_NS, LAST_NS = 200.0, 5.0, 50.0
SMOOTH_NS = 2.0                                        # smoothing before a loop state is assigned
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"  # validated categorical order
INK, MUTED, GRID, PALE = "#1f2328", "#6b7280", "#e5e7eb", "#f3f4f6"
# ---------------------------------------------------------------------------


def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5, length=0)


def runs_of(prefix):
    names = sorted(p.stem[len(prefix) + 1:] for p in ANALYSIS.glob(f"{prefix}_T*.csv"))
    return {n: pd.read_csv(ANALYSIS / f"{prefix}_{n}.csv") for n in names}


def temp(run):
    return int(run[1:4])


def fig_extension(ts):
    fig, axes = plt.subplots(len(EXTENDED), 1, figsize=(9, 2.35 * len(EXTENDED) + 0.9), sharex=True, dpi=160)
    for ax, (run, title) in zip(axes, EXTENDED):
        x = ts[run]
        t = x.time_ps / 1000.0
        per_ns = int(round(1000.0 / (x.time_ps.iloc[1] - x.time_ps.iloc[0])))
        ax.axvspan(0, PUBLISHED_NS, color=PALE, lw=0)
        ax.plot(t, x.rmsd_nm, color=BLUE, alpha=0.13, lw=0.4)
        ax.plot(t, x.rmsd_nm.rolling(5 * per_ns, center=True, min_periods=per_ns).mean(), color=BLUE, lw=2)
        pub = x[(t >= EQUIL_NS) & (t <= PUBLISHED_NS)].rmsd_nm.mean()
        ax.hlines(pub, EQUIL_NS, PUBLISHED_NS, color=INK, lw=1.2, ls=(0, (4, 3)))
        end = t.iloc[-1]
        last = x[t > end - LAST_NS].rmsd_nm.mean()
        ax.hlines(last, end - LAST_NS, end, color=INK, lw=1.2)
        ax.text(end + 6, last, f"last {LAST_NS:.0f} ns\n{last:.3f} nm", fontsize=8.5, color=INK, va="center")
        low = pub < 0.12
        ax.text(8 if low else 70, 0.040 if low else 0.085, f"dashed: mean over 5-200 ns, {pub:.3f} nm",
                fontsize=8.5, color=INK, va="bottom")
        ax.set_title(title, loc="left", fontsize=10.5, color=INK, pad=6)
        ax.set_ylim(0.03, 0.40)
        ax.set_xlim(0, 680)
        ax.set_ylabel("RMSD, nm", color=MUTED, fontsize=9)
        style(ax)
    axes[0].text(PUBLISHED_NS / 2, 0.37, "window published in Chapters 5 and 6", ha="center", fontsize=8.5, color=MUTED)
    axes[-1].set_xticks(range(0, 601, 100))
    axes[-1].set_xlabel("time, ns", color=MUTED, fontsize=9)
    fig.suptitle("Whole-protein C-alpha RMSD after extending four runs", x=0.065, ha="left", fontsize=13, color=INK, y=0.985)
    fig.text(0.065, 0.948, "Thick line: 5 ns rolling mean. Thin: every frame. Solid bar: mean of the last 50 ns.",
             fontsize=9, color=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "phase7_extension.png", facecolor="white")
    plt.close(fig)


def loop_state(x):
    """0 start, 1 form A, 2 form B, 3 open, 4 in between."""
    n = max(1, int(round(SMOOTH_NS * 1000.0 / (x.time_ps.iloc[1] - x.time_ps.iloc[0]))))
    sm = x[["rmsd_loop_nm", "d102_nm", "d105_nm"]].rolling(n, center=True, min_periods=1).mean()
    s = np.full(len(x), 4)
    s[(sm.d105_nm >= 0.30) & (sm.d102_nm < 0.50)] = 2
    s[(sm.d102_nm >= 0.50) & (sm.d105_nm < 0.30)] = 1
    s[sm.d105_nm >= 0.75] = 3
    s[sm.rmsd_loop_nm < 0.20] = 0
    return s


def fig_states(core):
    colors = ["#d1d5db", BLUE, ORANGE, GREEN, "#6b7280"]
    labels = ["starting conformation", "form A (Gly102 displaced)", "form B (residues 104-106 displaced)",
              "open (loop about 1 nm out)", "in between"]
    runs = sorted(core, key=lambda r: (temp(r), r))
    fig, ax = plt.subplots(figsize=(9, 0.46 * len(runs) + 1.9), dpi=160)
    for i, run in enumerate(runs):
        x = core[run]
        t = x.time_ps.values / 1000.0
        s = loop_state(x)
        edges = np.flatnonzero(np.diff(s)) + 1
        starts = np.r_[0, edges]
        stops = np.r_[edges, len(s) - 1]
        for a, b in zip(starts, stops):
            ax.barh(i, t[b] - t[a], left=t[a], height=0.62, color=colors[s[a]], lw=0)
        ax.text(t[-1] + 6, i, f"{t[-1]:.0f} ns", va="center", fontsize=8.5, color=MUTED)
    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels([r.replace("T", "").replace("_r", " K, replica ") + ("" if "_r" in r else " K") for r in runs],
                       fontsize=9, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 660)
    ax.set_xticks(range(0, 601, 100))
    ax.set_xlabel("time, ns", color=MUTED, fontsize=9)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(length=0, labelsize=8.5)
    ax.tick_params(axis="x", colors=MUTED)
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(colors, labels)], ncol=3, frameon=False,
              fontsize=8.5, loc="upper left", bbox_to_anchor=(-0.01, -0.22), handlelength=1.2, columnspacing=1.6)
    fig.suptitle("State of loop 100-117 in every run", x=0.065, ha="left", fontsize=13, color=INK, y=0.975)
    fig.text(0.065, 0.905, "Core-fitted C-alpha displacements, smoothed over 2 ns before a state is assigned.",
             fontsize=9, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.9))
    fig.savefig(OUT / "phase7_loop_states.png", facecolor="white", bbox_inches="tight")
    plt.close(fig)


def fig_ladder(ts, core):
    rows = []
    for run in core:
        a, c = ts[run], core[run]
        ta, tc = a.time_ps / 1000.0, c.time_ps / 1000.0
        rows.append({"T": temp(run),
                     "all_last": a[ta > ta.iloc[-1] - LAST_NS].rmsd_nm.mean(),
                     "core_last": c[tc > tc.iloc[-1] - LAST_NS].rmsd_core_nm.mean()})
    r = pd.DataFrame(rows)
    m = r.groupby("T").mean()
    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=160)
    T5 = sorted(CH5)
    series = [(T5, [CH5[t] for t in T5], None, None, BLUE, "whole protein, 5-200 ns, one run each (Chapter 5)", "o", (0, (4, 3))),
              (m.index, m.all_last, r["T"], r.all_last, ORANGE, "whole protein, last 50 ns of every run", "s", "-"),
              (m.index, m.core_last, r["T"], r.core_last, GREEN, "core only (residues 1-99), last 50 ns", "D", "-")]
    for xs, ys, px, py, col, lab, mk, ls in series:
        if px is not None:
            ax.scatter(px, py, s=22, color=col, alpha=0.35, lw=0, zorder=2)
        ax.plot(xs, ys, color=col, lw=2, ls=ls, marker=mk, ms=7, mec="white", mew=1.5, label=lab, zorder=3)
    ax.text(452, CH5[450] + 0.012, "Chapter 5", color=INK, fontsize=9, va="bottom", ha="right")
    ax.text(452, m.all_last.loc[450] - 0.012, "last 50 ns", color=INK, fontsize=9, va="top", ha="right")
    ax.text(452, m.core_last.loc[450] + 0.010, "core only", color=INK, fontsize=9, va="bottom", ha="right")
    ax.set_xticks(T5)
    ax.set_xlim(290, 460)
    ax.set_ylim(0, 0.30)
    ax.set_xlabel("temperature, K", color=MUTED, fontsize=9)
    ax.set_ylabel("C-alpha RMSD, nm", color=MUTED, fontsize=9)
    style(ax)
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    fig.suptitle("Three readings of the same temperature ladder", x=0.075, ha="left", fontsize=13, color=INK, y=0.98)
    fig.text(0.075, 0.915, "Lines join the mean per temperature. Pale dots are single runs (three each at 325 and 350 K).",
             fontsize=9, color=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(OUT / "phase7_core_ladder.png", facecolor="white")
    plt.close(fig)
    return r, m


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ts, core = runs_of("timeseries"), runs_of("core")
    if not ts or not core:
        raise SystemExit(f"no timeseries_T*.csv or core_T*.csv in {ANALYSIS}; run analyze_lyso.py and analyze_core.py first")
    fig_extension(ts)
    fig_states(core)
    _, m = fig_ladder(ts, core)
    print(m.round(4).to_string())
    print(f"wrote three figures to {OUT}")


if __name__ == "__main__":
    main()
