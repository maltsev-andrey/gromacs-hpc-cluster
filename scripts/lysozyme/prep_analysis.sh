#!/usr/bin/env bash
# Phase 5 preparation. Builds the two inputs analyze_lyso.py needs, which the
# production run did not leave behind:
#
#   common/ref_protein.gro the shared equilibrated structure, protein only.
#                          This is BOTH the topology and the RMSD reference.
#   common/protein.tpr     protein-only topology, kept for the atom-count check
#                          and for later: MDAnalysis 2.11 CANNOT read it, because
#                          GROMACS 2026.2 writes tpx version 138 and the parser
#                          does not support that yet. The GRO loses nothing that
#                          matters: MDAnalysis guesses masses from atom names, and
#                          the H-bond step is given explicit donor/hydrogen/
#                          acceptor selections so it needs no bonds or charges.
#
# Why not simply the existing prod.gro:
#   - prod.tpr and prod.gro cover the WHOLE system, protein plus water plus
#     ions, while prod.xtc holds protein only (compressed-x-grps = Protein).
#     MDAnalysis refuses a topology whose atom count differs from the
#     trajectory's, so the full-system files cannot be used.
#   - prod.gro is also the last frame of the run, not the equilibrated start.
#
# Why one shared reference rather than each replica's own frame 0:
#   every replica got its own 1 ns NVT equilibration at its own temperature, so
#   frame 0 at 450 K is already heated. Measuring all five against the same
#   300 K NPT-equilibrated structure is what makes the comparison mean anything.
#
# Run once, on the head node, with GMXRC sourced. Read-only on the run data.

set -euo pipefail

# ---------- settings ----------
LYSO="/nfs/shared/gromacs/lyso"
COMMON="${LYSO}/common"
GMXRC="/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/GMXRC"
SRC_T=300                       # which replica's tpr to reduce; any will do
GROUP="Protein"
# ------------------------------

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export GMX_MAXBACKUP=-1

cd "${COMMON}" || { echo "no ${COMMON}: adjust COMMON to wherever prep_lysozyme.sh put the shared files"; exit 1; }

src_tpr="${LYSO}/T${SRC_T}/prod.tpr"
[[ -f "${src_tpr}" ]] || src_tpr="${LYSO}/T${SRC_T}/run/prod.tpr"
[[ -f "${src_tpr}" ]] || { echo "cannot find prod.tpr for T${SRC_T}"; exit 1; }

echo "== protein-only topology from ${src_tpr}"
echo "${GROUP}" | gmx convert-tpr -s "${src_tpr}" -o protein.tpr

# The shared equilibrated structure. prep_lysozyme.sh left the NPT output here;
# adjust NPT if it is named differently.
NPT="npt.gro"
[[ -f "${NPT}" ]] || NPT="$(ls -1 npt*.gro 2>/dev/null | head -1 || true)"
[[ -n "${NPT}" && -f "${NPT}" ]] || { echo "no npt*.gro in ${COMMON}: set NPT by hand"; ls -1; exit 1; }

echo "== protein-only reference from ${NPT}"
echo "${GROUP}" | gmx trjconv -s "${src_tpr}" -f "${NPT}" -o ref_protein.gro

echo
echo "== verify: all three atom counts must match"
printf '  protein.tpr       '; gmx check -s protein.tpr 2>&1 | awk '/^ *#atoms/{print $2}'
printf '  ref_protein.gro   '; sed -n '2p' ref_protein.gro | tr -d ' '
printf '  prod.xtc (T%-3s)   ' "${SRC_T}"
gmx check -f "${LYSO}/T${SRC_T}/run/prod.xtc" 2>&1 | awk '/^# Atoms/{print $3}'
echo
echo "If those three disagree, stop: analyze_lyso.py will refuse rather than"
echo "produce wrong numbers, but the cause is here."
