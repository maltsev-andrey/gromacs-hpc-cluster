#!/usr/bin/env python3
"""Which residues move when a run leaves its starting basin?

For each run two time windows are compared: one before the RMSD jump and one
after it. Every frame is fitted on all C-alpha atoms to the common reference,
then the mean C-alpha position per residue is taken in each window.

  shift_nm        distance between the two window means (what moved)
  dev_before/after  distance of each window mean from the reference
  rmsf_before/after fluctuation around the window's own mean

No arguments. Edit the settings block and run:  python analyze_transition.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------- settings
LYSO   = Path("/nfs/shared/gromacs/lyso")
TOPOL  = LYSO / "common" / "ref_protein.gro"      # protein-only topology, also the fit reference
TRAJ   = "prod_whole.xtc"
OUT    = LYSO / "analysis"
STRIDE = 10                                       # every 10th frame = one per 100 ps

# run: (before_start, before_end, after_start, after_end), all in ns
WINDOWS = {
    "T300":    ( 20, 150, 400, 590),
    "T325":    ( 20, 280, 310, 400),
    "T325_r2": (  5,  60, 100, 200),
    "T325_r3": (  5,  60, 100, 200),
    "T400":    (100, 195, 250, 400),
}
TOP_N = 12                                        # residues listed per run
# ---------------------------------------------------------------------------


def kabsch(mobile, target):
    """Rotation that puts centred `mobile` onto centred `target` (both N x 3)."""
    h = mobile.T @ target
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1.0, 1.0, d]) @ vt


def window_stats(u, ca, ref0, start_ns, end_ns):
    dt_ns = u.trajectory.dt / 1000.0
    last = (u.trajectory.n_frames - 1) * dt_ns
    if end_ns > last + 1e-6:
        raise SystemExit(f"window end {end_ns} ns is beyond the trajectory ({last:.0f} ns)")
    i0, i1 = int(round(start_ns / dt_ns)), int(round(end_ns / dt_ns))
    frames = []
    for _ in u.trajectory[i0:i1 + 1:STRIDE]:
        x = ca.positions.astype(float) / 10.0            # Angstrom -> nm
        x -= x.mean(axis=0)
        x = x @ kabsch(x, ref0)
        frames.append(x)
    frames = np.asarray(frames)
    mean = frames.mean(axis=0)
    rmsf = np.sqrt(((frames - mean) ** 2).sum(axis=2).mean(axis=0))
    return mean, rmsf, len(frames)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ref_u = mda.Universe(str(TOPOL))
    ref_ca = ref_u.select_atoms("name CA")
    ref0 = ref_ca.positions.astype(float) / 10.0
    ref0 -= ref0.mean(axis=0)
    resid, resname = ref_ca.resids, ref_ca.resnames

    shifts = {}
    for run, (b0, b1, a0, a1) in WINDOWS.items():
        xtc = LYSO / run / "run" / TRAJ
        if not xtc.exists():
            print(f"{run}: {xtc} not found, skipped")
            continue
        u = mda.Universe(str(TOPOL), str(xtc))
        ca = u.select_atoms("name CA")
        if len(ca) != len(ref0):
            raise SystemExit(f"{run}: {len(ca)} C-alpha atoms, reference has {len(ref0)}")
        mb, fb, nb = window_stats(u, ca, ref0, b0, b1)
        ma, fa, na = window_stats(u, ca, ref0, a0, a1)
        df = pd.DataFrame({
            "resid": resid, "resname": resname,
            "dev_before_nm": np.linalg.norm(mb - ref0, axis=1),
            "dev_after_nm":  np.linalg.norm(ma - ref0, axis=1),
            "shift_nm":      np.linalg.norm(ma - mb, axis=1),
            "rmsf_before_nm": fb, "rmsf_after_nm": fa,
        }).round(4)
        df.to_csv(OUT / f"transition_{run}.csv", index=False)
        shifts[run] = df["shift_nm"].values

        tot = float((df.shift_nm ** 2).sum())
        top = df.sort_values("shift_nm", ascending=False).head(TOP_N)
        share = float((top.shift_nm ** 2).sum()) / tot
        print(f"\n== {run}: before {b0}-{b1} ns ({nb} frames), after {a0}-{a1} ns ({na} frames)")
        print(f"   RMS shift over all residues {np.sqrt(tot / len(df)):.4f} nm; "
              f"the top {TOP_N} residues carry {100 * share:.0f}% of the squared shift")
        print(top[["resid", "resname", "shift_nm", "dev_before_nm", "dev_after_nm",
                   "rmsf_before_nm", "rmsf_after_nm"]].to_string(index=False))

    if len(shifts) > 1:
        print("\n== do the runs move the same residues? Spearman correlation of the shift profiles")
        print(pd.DataFrame(shifts).corr(method="spearman").round(2).to_string())
        pd.DataFrame({"resid": resid, "resname": resname, **shifts}).round(4) \
            .to_csv(OUT / "transition_shift_all.csv", index=False)

    if shifts:
        fig, axes = plt.subplots(len(shifts), 1, figsize=(9, 1.9 * len(shifts) + 0.8),
                                 sharex=True, sharey=True, dpi=150, squeeze=False)
        for ax, (run, s) in zip(axes[:, 0], shifts.items()):
            b0, b1, a0, a1 = WINDOWS[run]
            ax.bar(resid, s, width=0.8, color="#2a78d6")
            ax.set_title(f"{run}: {b0}-{b1} ns vs {a0}-{a1} ns", loc="left", fontsize=10)
            ax.set_ylabel("shift, nm", fontsize=9)
            ax.grid(axis="y", color="#e5e7eb", lw=0.8)
            ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        axes[-1, 0].set_xlabel("residue", fontsize=9)
        fig.suptitle("Per-residue C-alpha shift across the RMSD transition", x=0.07, ha="left", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(OUT / "lyso_transition_shift.png", facecolor="white")
        print(f"\nwrote transition_*.csv and lyso_transition_shift.png to {OUT}")


if __name__ == "__main__":
    main()
