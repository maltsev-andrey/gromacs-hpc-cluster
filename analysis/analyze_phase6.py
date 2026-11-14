#!/usr/bin/env python3
"""
Phase 6: does the RMSD step from 325 K to 350 K survive replication?

Phase 5 measured one trajectory per temperature and found a +68% jump between
325 K and 350 K, against +7%, +25% and +35% for the other intervals. It sits
next to lysozyme's experimental melting point, which made it interesting and
made one trajectory not enough.

The error bars in the Phase 5 summary are the spread WITHIN each trajectory.
They are the wrong uncertainty for comparing temperatures: consecutive frames
are correlated, so that spread says how much the molecule wobbles, not how much
an independent 200 ns run of it would differ. This script measures the spread
BETWEEN replicas instead, which is what a difference between temperatures has
to beat.

Reads the per-run timeseries CSVs that analyze_lyso.py writes. Run that first.

    ~/venv-mdanalysis/bin/python analyze_phase6.py
"""

import re
import sys
from pathlib import Path

# ---------- settings ----------
OUT = Path("/nfs/shared/gromacs/lyso/analysis")
EQUIL_NS = 5.0          # must match what analyze_lyso.py discarded
COMPARE = (325, 350)    # the interval Phase 5 flagged
MEASURES = [("rmsd_nm", "RMSD (nm)"), ("rg_nm", "Rg (nm)"), ("hbonds", "backbone H-bonds")]
# ------------------------------

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


def runs():
    """run name -> (temperature, replica number, dataframe)"""
    out = {}
    for f in sorted(OUT.glob("timeseries_T*.csv")):
        name = f.stem.replace("timeseries_", "")
        m = re.match(r"T(\d+)(?:_r(\d+))?$", name)
        if not m:
            print(f"  skipping {f.name}: cannot read a temperature from the name")
            continue
        out[name] = (int(m.group(1)), int(m.group(2) or 1), pd.read_csv(f))
    return out


def replica_means(data):
    rows = []
    for name, (T, r, df) in data.items():
        keep = df[df["time_ps"] >= EQUIL_NS * 1000]
        if keep.empty:
            keep = df
        row = {"run": name, "T_K": T, "replica": r, "n_frames_used": len(keep)}
        for col, _ in MEASURES:
            if col in keep:
                row[col] = keep[col].mean()
                row[f"{col}_within_sd"] = keep[col].std()
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["T_K", "replica"]).reset_index(drop=True)


def by_temperature(per_run):
    rows = []
    for T, g in per_run.groupby("T_K"):
        row = {"T_K": T, "n_replicas": len(g)}
        for col, _ in MEASURES:
            if col not in g:
                continue
            v = g[col].to_numpy()
            row[f"{col}_mean"] = v.mean()
            # ddof=1 needs n>=2; with one replica there is no between-replica spread
            row[f"{col}_between_sd"] = v.std(ddof=1) if len(v) > 1 else np.nan
            row[f"{col}_sem"] = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
            row[f"{col}_within_sd"] = g[f"{col}_within_sd"].mean()
        rows.append(row)
    return pd.DataFrame(rows)


def compare(per_run, a, b):
    print(f"\n== is the {a} K to {b} K step real? ==")
    for col, label in MEASURES:
        ga = per_run.loc[per_run.T_K == a, col].dropna().to_numpy()
        gb = per_run.loc[per_run.T_K == b, col].dropna().to_numpy()
        if len(ga) < 2 or len(gb) < 2:
            print(f"  {label:18s} needs >=2 replicas at both: have {len(ga)} and {len(gb)}")
            continue
        diff = gb.mean() - ga.mean()
        t, p = stats.ttest_ind(gb, ga, equal_var=False)
        # Welch-Satterthwaite df, and a CI on the difference
        se = np.sqrt(ga.var(ddof=1) / len(ga) + gb.var(ddof=1) / len(gb))
        df = se**4 / (
            (ga.var(ddof=1) / len(ga)) ** 2 / (len(ga) - 1)
            + (gb.var(ddof=1) / len(gb)) ** 2 / (len(gb) - 1))
        crit = stats.t.ppf(0.975, df)
        lo, hi = diff - crit * se, diff + crit * se
        pooled = np.sqrt((ga.var(ddof=1) + gb.var(ddof=1)) / 2)
        print(f"  {label:18s} {ga.mean():.4f} -> {gb.mean():.4f}   "
              f"diff {diff:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]")
        d = f"   Cohen's d={diff/pooled:.1f}" if pooled > 0 else ""
        print(f"  {'':18s} Welch t={t:.2f}, df={df:.1f}, p={p:.4f}{d}")
        verdict = ("the interval is resolved: the CI excludes zero"
                   if lo * hi > 0 else
                   "not resolved: the CI spans zero, so these replicas cannot separate them")
        print(f"  {'':18s} {verdict}\n")


def error_bar_lesson(by_T):
    print("== which error bar were we using? ==")
    print("  Phase 5 quoted the within-trajectory spread. For comparing temperatures")
    print("  the between-replica spread is the relevant one. Ratio of the two:\n")
    print(f"  {'T (K)':>6} {'n':>3} {'within SD':>11} {'between SD':>11} {'ratio':>7}")
    for _, r in by_T.iterrows():
        if np.isnan(r.get("rmsd_nm_between_sd", np.nan)):
            print(f"  {int(r.T_K):>6} {int(r.n_replicas):>3} {r['rmsd_nm_within_sd']:>11.4f} "
                  f"{'-':>11} {'-':>7}   (one replica)")
        else:
            ratio = r["rmsd_nm_between_sd"] / r["rmsd_nm_within_sd"]
            print(f"  {int(r.T_K):>6} {int(r.n_replicas):>3} {r['rmsd_nm_within_sd']:>11.4f} "
                  f"{r['rmsd_nm_between_sd']:>11.4f} {ratio:>7.2f}")
    print()


def chart(per_run, by_T):
    fig, axes = plt.subplots(1, len(MEASURES), figsize=(4.2 * len(MEASURES), 3.8))
    axes = np.atleast_1d(axes)
    for ax, (col, label) in zip(axes, MEASURES):
        if f"{col}_mean" not in by_T:        # by_T holds <col>_mean, not <col>
            continue
        m = by_T.dropna(subset=[f"{col}_mean"])
        ax.errorbar(m["T_K"], m[f"{col}_mean"], yerr=m[f"{col}_sem"],
                    marker="o", capsize=4, color="#333", zorder=3,
                    label="mean ± SEM across replicas")
        single = per_run.dropna(subset=[col])
        ax.scatter(single["T_K"] + 3, single[col], s=18, color="#c2493d",
                   alpha=.75, zorder=2, label="individual replicas")
        ax.set_xlabel("temperature (K)")
        ax.set_ylabel(label)
        ax.grid(alpha=.3)
    axes[0].legend(fontsize=7, loc="best")
    fig.suptitle("Phase 6: replicated at 325 K and 350 K. Error bars are between-replica, "
                 "not within-trajectory", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "phase6_replicated.png", dpi=110)
    plt.close(fig)


def main():
    if not OUT.exists():
        sys.exit(f"no {OUT}: run analyze_lyso.py first")
    data = runs()
    if not data:
        sys.exit(f"no timeseries_T*.csv in {OUT}: run analyze_lyso.py first")

    per_run = replica_means(data)
    by_T = by_temperature(per_run)

    print("== per replica, averaged over the equilibrated portion ==")
    cols = ["run", "T_K", "replica"] + [c for c, _ in MEASURES if c in per_run]
    print(per_run[cols].round(4).to_string(index=False))

    print("\n== by temperature ==")
    show = ["T_K", "n_replicas"] + [c for c in by_T.columns
                                    if c.endswith(("_mean", "_between_sd", "_sem"))]
    print(by_T[show].round(4).to_string(index=False))

    per_run.round(5).to_csv(OUT / "phase6_per_replica.csv", index=False)
    by_T.round(5).to_csv(OUT / "phase6_by_temperature.csv", index=False)

    compare(per_run, *COMPARE)
    error_bar_lesson(by_T)
    chart(per_run, by_T)
    print(f"wrote phase6_per_replica.csv, phase6_by_temperature.csv and "
          f"phase6_replicated.png to {OUT}")


if __name__ == "__main__":
    main()
