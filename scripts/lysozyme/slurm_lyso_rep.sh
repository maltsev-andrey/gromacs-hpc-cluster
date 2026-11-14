#!/usr/bin/env bash
#SBATCH --job-name=lyso-rep
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --time=23:00:00
#SBATCH --output=/nfs/shared/gromacs/results/slurm/%x_%j.out
#
# Phase 6 production replica on one compute node. Generalises
# slurm_lyso_node.sh: that one keys off a temperature, this one off a run
# directory, so replicas at the same temperature do not collide.
#
#   for RUN in T325_r2 T325_r3 T350_r2 T350_r3; do
#     sbatch -J lyso-$RUN --export=RUN=$RUN slurm_lyso_rep.sh
#   done
#
# Runs up to 22 h, checkpoints, and resubmits itself until prod.gro exists.
# --dependency=singleton keeps one job per name in the queue, so the chain for
# a given replica never overlaps itself.
#
# MAXH is 22 against a 23 h wall limit. The Phase 4 runs used 23.5 against 24,
# which left half an hour of margin; if a window is ever killed by the wall
# limit instead of stopping on -maxh, everything since the last checkpoint is
# lost. An hour of margin costs about 4% of throughput and removes that risk.

set -euo pipefail

# ---------- settings ----------
RUN="${RUN:?set RUN=T<temp>_r<n> via --export, e.g. --export=RUN=T350_r2}"
GMXRC="/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/GMXRC"
LYSO="/nfs/shared/gromacs/lyso"
MAXH=22
NTOMP=12
# ------------------------------

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export GMX_MAXBACKUP=-1

src="${LYSO}/${RUN}"
[[ -f "${src}/prod.tpr" ]] || { echo "no ${src}/prod.tpr: run prep_phase6.sh first"; exit 1; }

run="${src}/run"
mkdir -p "${run}" && cd "${run}"
[[ -f prod.tpr ]] || cp "${src}/prod.tpr" .

if [[ -f prod.gro ]]; then echo "${RUN}: complete"; exit 0; fi

cpi=""; [[ -f prod.cpt ]] && cpi="-cpi prod.cpt -append"
echo "$(date -Iseconds) ${RUN} on $(hostname -s): window ${MAXH} h ${cpi:+from checkpoint}"

status=0
# shellcheck disable=SC2086
gmx mdrun -deffnm prod ${cpi} -ntmpi 1 -ntomp "${NTOMP}" -pin on -maxh "${MAXH}" \
  >> mdrun.out 2>&1 || status=$?

done_steps=$(grep -oE 'Writing checkpoint, step [0-9]+' prod.log | tail -1 | grep -oE '[0-9]+$' || echo 0)
echo "$(date -Iseconds) ${RUN}: exit ${status}, step ${done_steps} of 100000000 ($(awk -v s="${done_steps}" 'BEGIN{printf "%.1f", s/500000}') ns)"
grep -E '^Performance:' prod.log | tail -1 || true

if [[ ! -f prod.gro && ${status} -eq 0 ]]; then
  sbatch -J "${SLURM_JOB_NAME}" --export=RUN="${RUN}" --dependency=singleton "$0"
fi
exit "${status}"
