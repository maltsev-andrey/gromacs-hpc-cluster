#!/usr/bin/env python3
"""
Phase 5: structural analysis of the lysozyme temperature ladder.

Per replica: RMSD against the shared equilibrated structure, per-residue RMSF,
radius of gyration, and intra-protein hydrogen-bond count.
Across the ladder: how those four move with temperature, which is the result.

Inputs  prod.xtc per temperature (protein only, 10 ps per frame) plus one
        protein-only topology and one shared reference structure, both made by
        prep_analysis.sh.
Outputs CSVs and charts under OUT. Nothing is written next to the trajectories.

All settings are in the block below; the script takes no arguments.
Run it with the analysis venv's interpreter by absolute path, e.g.
  ~/venv-mdanalysis/bin/python /nfs/shared/projects/gromacs/lyso/analyze_lyso.py
"""

import os
import sys
import warnings
from pathlib import Path

# ---------- settings ----------
LYSO = Path("/nfs/shared/gromacs/lyso")
# Protein-only topology. A GRO, not the TPR: GROMACS 2026.2 writes tpx version
# 138 and MDAnalysis 2.11's TPR parser does not support it yet ("Your tpx
# version is 138, which this parser does not support, yet"). The GRO costs
# nothing here, because MDAnalysis guesses masses from atom names (N 14.007,
# C 12.011, O 15.999, H 1.008), which is what mass-weighted Rg needs, and the
# H-bond analysis takes explicit donor/hydrogen/acceptor selections so it needs
# neither bonds nor charges. Switch back to protein.tpr if MDAnalysis gains
# tpx 138 support; prep_analysis.sh builds it either way.
TOPOL = LYSO / "common" / "ref_protein.gro"
REF = LYSO / "common" / "ref_protein.gro"      # shared equilibrated reference
TEMPS = [300, 325, 350, 400, 450]
OUT = LYSO / "analysis"

TRAJ = "prod_whole.xtc"   # written by make_whole.sh. The raw prod.xtc has the
                          # protein split across the periodic boundary for most
                          # frames, which inflates RMSD, RMSF and Rg several-fold.
EQUIL_NS = 5.0        # discarded from the start of each run before averaging
STRIDE = 1            # frame stride for RMSD / Rg (1 = every frame)
HB_STRIDE = 10        # frame stride for H-bonds: O(pairs) per frame, so coarser
SEL_FIT = "protein and name CA"   # what RMSD and RMSF are measured on
SEL_RG = "protein"                # Rg over all protein atoms, mass weighted

# Backbone hydrogen bonds: amide H donated by N, accepted by carbonyl O. These
# are the bonds secondary structure is made of, so the count is a direct read on
# whether the fold is intact. Given explicitly rather than guessed: MDAnalysis
# guesses donors and acceptors from partial charges, which only a TPR carries,
# and guessing also pulls in side-chain and solvent-facing bonds that add noise
# without telling you about the fold.
SEL_DON = "protein and name N"
SEL_H = "protein and name H"
SEL_ACC = "protein and name O"
HB_DA_CUTOFF = 3.5      # Angstrom, donor to acceptor
HB_ANGLE_CUTOFF = 150   # degrees, D-H...A
# ------------------------------

warnings.filterwarnings("ignore", category=UserWarning)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA

A_TO_NM = 0.1         # MDAnalysis works in Angstrom; GROMACS convention is nm


def load(T):
    traj = LYSO / f"T{T}" / "run" / TRAJ
    if not traj.exists():
        alt = traj.parent / "prod.xtc"
        hint = "  (prod.xtc is there: run make_whole.sh first)" if alt.exists() else ""
        return None, f"no trajectory at {traj}{hint}"
    u = mda.Universe(str(TOPOL), str(traj))
    n_top, n_traj = len(u.atoms), u.trajectory.n_atoms
    if n_top != n_traj:
        return None, (f"topology has {n_top} atoms, trajectory has {n_traj}. "
                      "TOPOL must be protein-only: see prep_analysis.sh")
    return u, None


def per_replica(T, u, ref):
    """Returns (timeseries DataFrame, rmsf DataFrame, note)."""
    sel = u.select_atoms(SEL_FIT)
    if len(sel) == 0:
        return None, None, f"selection '{SEL_FIT}' matched no atoms"

    # --- RMSD against the shared reference, with superposition
    R = rms.RMSD(u, ref, select=SEL_FIT, ref_frame=0).run(step=STRIDE)
    ts = pd.DataFrame({"time_ps": R.results.rmsd[:, 1],
                       "rmsd_nm": R.results.rmsd[:, 2] * A_TO_NM})

    # --- radius of gyration, same frames
    rg = []
    grp = u.select_atoms(SEL_RG)
    for _ in u.trajectory[::STRIDE]:
        rg.append(grp.radius_of_gyration() * A_TO_NM)
    ts["rg_nm"] = rg

    # Sanity guard. A folded protein's Rg barely moves; a protein broken across
    # the periodic boundary swings by a factor of two or more. If the median sits
    # well above the low percentile, the trajectory was not made whole.
    lo, med = ts["rg_nm"].quantile(0.05), ts["rg_nm"].median()
    if med > 1.25 * lo:
        print(f"      WARNING  Rg median {med:.2f} nm vs 5th percentile {lo:.2f} nm. "
              "That spread means the molecule is split across the periodic "
              "boundary. RMSD, RMSF and Rg below are not trustworthy. "
              "Run make_whole.sh and rerun.", flush=True)

    # --- intra-protein H-bonds, coarser stride
    for name, s in (("donors", SEL_DON), ("hydrogens", SEL_H), ("acceptors", SEL_ACC)):
        if len(u.select_atoms(s)) == 0:
            return None, None, (f"H-bond {name} selection '{s}' matched no atoms. "
                                "Check the atom naming in TOPOL.")
    hb = HBA(universe=u, donors_sel=SEL_DON, hydrogens_sel=SEL_H, acceptors_sel=SEL_ACC,
             d_a_cutoff=HB_DA_CUTOFF, d_h_a_angle_cutoff=HB_ANGLE_CUTOFF)
    hb.run(step=HB_STRIDE)
    hb_counts = hb.count_by_time()
    hb_times = np.arange(len(hb_counts)) * HB_STRIDE * dt_ps(u)
    hbs = pd.DataFrame({"time_ps": hb_times, "hbonds": hb_counts})
    ts = pd.merge_asof(ts.sort_values("time_ps"), hbs.sort_values("time_ps"),
                       on="time_ps", direction="nearest",
                       tolerance=HB_STRIDE * dt_ps(u))

    # --- RMSF per residue, on the post-equilibration portion, aligned to its
    #     own average structure (RMSF is about fluctuation, not displacement)
    start = int(EQUIL_NS * 1000 / dt_ps(u))
    avg = align.AverageStructure(u, u, select=SEL_FIT, ref_frame=0).run(start=start)
    aligner = align.AlignTraj(u, avg.results.universe, select=SEL_FIT, in_memory=True)
    aligner.run(start=start)
    # start= matters: AlignTraj only aligned frames from `start` onward, so a
    # bare .run() would mix in unaligned leading frames and report the protein's
    # diffusion through the box as fluctuation. Measured 5.3x inflation on a
    # test trajectory that drifts 2.4 nm.
    F = rms.RMSF(u.select_atoms(SEL_FIT)).run(start=start)
    rmsf = pd.DataFrame({"resid": u.select_atoms(SEL_FIT).resids,
                         "rmsf_nm": F.results.rmsf * A_TO_NM})
    return ts, rmsf, None


def dt_ps(u):
    return float(u.trajectory.dt)


def summarise(T, ts):
    """Averages over the equilibrated portion only."""
    keep = ts[ts["time_ps"] >= EQUIL_NS * 1000]
    if keep.empty:
        keep = ts
    row = {"T_K": T, "n_frames": len(ts), "n_frames_used": len(keep),
           "equil_discarded_ns": EQUIL_NS}
    for col, label in (("rmsd_nm", "rmsd"), ("rg_nm", "rg"), ("hbonds", "hbonds")):
        if col in keep:
            row[f"{label}_mean"] = keep[col].mean()
            row[f"{label}_sd"] = keep[col].std()
    return row


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in (TOPOL, REF):
        if not Path(p).exists():
            sys.exit(f"missing {p}\nRun prep_analysis.sh first: it builds the "
                     "protein-only topology and the shared reference.")
    ref = mda.Universe(str(TOPOL), str(REF))

    rows, series, rmsfs = [], {}, {}
    for T in TEMPS:
        u, err = load(T)
        if err:
            print(f"T{T}: skipped, {err}")
            continue
        prot, ca = len(u.select_atoms("protein")), len(u.select_atoms(SEL_FIT))
        print(f"T{T}: {len(u.atoms)} atoms ({prot} matched 'protein', {ca} C-alpha), "
              f"{u.trajectory.n_frames} frames, dt {dt_ps(u):.0f} ps", flush=True)
        if prot < len(u.atoms):
            print(f"      note: {len(u.atoms) - prot} atoms did not match 'protein'. "
                  "Check for non-standard residue names (HID/HIE/HIP, termini).")
        ts, rmsf, err = per_replica(T, u, ref)
        if err:
            print(f"T{T}: skipped, {err}")
            continue
        ts.to_csv(OUT / f"timeseries_T{T}.csv", index=False)
        rmsf.to_csv(OUT / f"rmsf_T{T}.csv", index=False)
        series[T], rmsfs[T] = ts, rmsf
        rows.append(summarise(T, ts))

    if not rows:
        sys.exit("nothing analysed")

    summary = pd.DataFrame(rows).round(4)
    summary.to_csv(OUT / "summary_by_temperature.csv", index=False)
    print("\n== averages over the equilibrated portion "
          f"(first {EQUIL_NS:g} ns discarded) ==")
    print(summary.to_string(index=False))

    charts(series, rmsfs, summary)
    print(f"\nwrote {len(series) * 2 + 1} CSVs and 3 charts to {OUT}")


def charts(series, rmsfs, summary):
    cmap = plt.get_cmap("coolwarm")
    temps = sorted(series)
    colour = {T: cmap((T - min(temps)) / max(1, (max(temps) - min(temps))))
              for T in temps}

    # 1. three traces against time, one row each
    panels = [("rmsd_nm", "RMSD vs equilibrated structure (nm)"),
              ("rg_nm", "Radius of gyration (nm)"),
              ("hbonds", "Intra-protein H-bonds")]
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 9), sharex=True)
    for ax, (col, label) in zip(axes, panels):
        for T in temps:
            d = series[T].dropna(subset=[col])
            if d.empty:
                continue
            ax.plot(d["time_ps"] / 1000, d[col], lw=0.7, color=colour[T], label=f"{T} K")
        ax.axvspan(0, EQUIL_NS, color="0.9", zorder=0)
        ax.set_ylabel(label, fontsize=9)
        ax.grid(alpha=.3)
    axes[0].legend(ncol=len(temps), fontsize=8, loc="upper left")
    span = max(s["time_ps"].max() for s in series.values()) / 1000
    axes[0].set_title(f"Lysozyme (PDB 1AKI), {span:.0f} ns per temperature. "
                      f"Shaded = first {EQUIL_NS:g} ns, excluded from averages", fontsize=10)
    axes[-1].set_xlabel("time (ns)")
    fig.tight_layout(); fig.savefig(OUT / "timeseries.png", dpi=110); plt.close(fig)

    # 2. per-residue RMSF, overlaid
    fig, ax = plt.subplots(figsize=(11, 4))
    for T in temps:
        ax.plot(rmsfs[T]["resid"], rmsfs[T]["rmsf_nm"], lw=1, color=colour[T], label=f"{T} K")
    ax.set_xlabel("residue"); ax.set_ylabel("C$\\alpha$ RMSF (nm)")
    ax.set_title("Per-residue fluctuation: which parts of the fold go first", fontsize=10)
    ax.legend(ncol=len(temps), fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(OUT / "rmsf.png", dpi=110); plt.close(fig)

    # 3. the result: each measure against temperature
    cols = [c for c in ("rmsd_mean", "rg_mean", "hbonds_mean") if c in summary]
    fig, axes = plt.subplots(1, len(cols), figsize=(4 * len(cols), 3.6))
    axes = np.atleast_1d(axes)
    for ax, col in zip(axes, cols):
        sd = col.replace("_mean", "_sd")
        ax.errorbar(summary["T_K"], summary[col],
                    yerr=summary[sd] if sd in summary else None,
                    marker="o", capsize=3, color="#444")
        ax.set_xlabel("temperature (K)")
        ax.set_ylabel(col.replace("_mean", "").upper())
        ax.grid(alpha=.3)
    fig.suptitle("Thermal response, averaged over the equilibrated portion", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "thermal_response.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
