#!/usr/bin/env python3
"""Separate the loop rearrangement from the rest of the protein.

Every frame is fitted to the reference on the C-alpha atoms of the CORE only.
With that fit in place the script reports, per frame:

  rmsd_core_nm   RMSD of the core itself (what temperature does to the fold)
  rmsd_loop_nm   RMSD of the loop, not refitted (has the loop moved?)
  rmsd_tail_nm   RMSD of the C-terminal tail, not refitted
  d<resid>_nm    displacement of single marker residues from the reference

No arguments. Edit the settings block and run:  python analyze_core.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda

# ----------------------------------------------------------------- settings
LYSO    = Path("/nfs/shared/gromacs/lyso")
TOPOL   = LYSO / "common" / "ref_protein.gro"
TRAJ    = "prod_whole.xtc"
OUT     = LYSO / "analysis"
RUNS    = []                 # empty = every T*/run that holds TRAJ
STRIDE  = 10                 # every 10th frame = one per 100 ps
CORE    = (1, 99)            # residue ranges, inclusive
LOOP    = (100, 117)
TAIL    = (118, 129)
MARKERS = [102, 105]         # 102 leads loop form A, 104-106 lead form B
EQUIL_NS  = 5.0
COMMON_NS = 200.0            # window every run has, for like-for-like means
LAST_NS   = 50.0
# ---------------------------------------------------------------------------


def kabsch(mobile, target):
    u, _, vt = np.linalg.svd(mobile.T @ target)
    d = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1.0, 1.0, d]) @ vt


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ref = mda.Universe(str(TOPOL)).select_atoms("name CA")
    rid = ref.resids
    m_core = (rid >= CORE[0]) & (rid <= CORE[1])
    m_loop = (rid >= LOOP[0]) & (rid <= LOOP[1])
    m_tail = (rid >= TAIL[0]) & (rid <= TAIL[1])
    i_mark = [int(np.where(rid == r)[0][0]) for r in MARKERS]
    r0 = ref.positions.astype(float) / 10.0
    r0 = r0 - r0[m_core].mean(axis=0)
    print(f"core {m_core.sum()} residues, loop {m_loop.sum()}, tail {m_tail.sum()}")

    runs = RUNS or sorted(p.parent.parent.name for p in LYSO.glob(f"T*/run/{TRAJ}"))
    rows = []
    for run in runs:
        u = mda.Universe(str(TOPOL), str(LYSO / run / "run" / TRAJ))
        ca = u.select_atoms("name CA")
        if len(ca) != len(r0):
            raise SystemExit(f"{run}: {len(ca)} C-alpha atoms, reference has {len(r0)}")
        rec = []
        for ts in u.trajectory[::STRIDE]:
            x = ca.positions.astype(float) / 10.0
            x = x - x[m_core].mean(axis=0)
            x = x @ kabsch(x[m_core], r0[m_core])
            d = np.linalg.norm(x - r0, axis=1)
            rec.append([ts.time,
                        np.sqrt((d[m_core] ** 2).mean()),
                        np.sqrt((d[m_loop] ** 2).mean()),
                        np.sqrt((d[m_tail] ** 2).mean()),
                        *d[i_mark]])
        cols = ["time_ps", "rmsd_core_nm", "rmsd_loop_nm", "rmsd_tail_nm"] + [f"d{r}_nm" for r in MARKERS]
        df = pd.DataFrame(rec, columns=cols)
        df.round(4).to_csv(OUT / f"core_{run}.csv", index=False)

        t = df.time_ps / 1000.0
        length = t.iloc[-1]
        common = df[(t >= EQUIL_NS) & (t <= COMMON_NS)]
        last = df[t > length - LAST_NS]
        row = {"run": run, "length_ns": round(length)}
        for c in cols[1:]:
            key = c.replace("_nm", "")
            row[f"{key}_common"] = common[c].mean()
            row[f"{key}_last"] = last[c].mean()
        rows.append(row)
        print(f"{run}: {len(df)} frames used, {length:.0f} ns")

    s = pd.DataFrame(rows).round(4)
    s.to_csv(OUT / "core_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(f"\n== means, nm.  _common = {EQUIL_NS:.0f}-{COMMON_NS:.0f} ns,  _last = final {LAST_NS:.0f} ns of each run ==")
    print(s.to_string(index=False))
    print(f"\nwrote core_<run>.csv and core_summary.csv to {OUT}")


if __name__ == "__main__":
    main()
