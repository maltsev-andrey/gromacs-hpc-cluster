#!/usr/bin/env bash
#SBATCH --job-name=lyso
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=/nfs/shared/gromacs/results/slurm/lyso_T%x_%j.out
#
# Phase 4 production replica on one compute node (thread-MPI build), one
# temperature per node. Runs up to 23.5 h, checkpoints, and resubmits itself
# until prod.gro exists. The chain is bound to the temperature by job name:
#
#   for T in 300 325 350 400; do sbatch -J lyso-T$T --export=T=$T slurm_lyso_node.sh; done
#
# --dependency=singleton on the resubmit keeps one job per name in the queue.

set -euo pipefail

# ---------- settings ----------
T="${T:?set T=<K> via --export}"
GMXRC="/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/GMXRC"
LYSO="/nfs/shared/gromacs/lyso"
MAXH=23.5
NTOMP=12
# ------------------------------

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export GMX_MAXBACKUP=-1

run="${LYSO}/T${T}/run"
mkdir -p "${run}" && cd "${run}"
[[ -f prod.tpr ]] || cp ../prod.tpr .

if [[ -f prod.gro ]]; then echo "T=${T}: complete"; exit 0; fi

cpi=""; [[ -f prod.cpt ]] && cpi="-cpi prod.cpt -append"
echo "$(date -Iseconds) T=${T} on $(hostname -s): window ${MAXH} h ${cpi:+from checkpoint}"
# shellcheck disable=SC2086
gmx mdrun -deffnm prod ${cpi} -ntmpi 1 -ntomp "${NTOMP}" -pin on -maxh "${MAXH}" >> mdrun.out 2>&1
status=$?
done_steps=$(grep -oE 'Writing checkpoint, step [0-9]+' prod.log | tail -1 | grep -oE '[0-9]+$' || echo 0)
echo "$(date -Iseconds) T=${T}: exit ${status}, step ${done_steps} ($(awk -v s="${done_steps}" 'BEGIN{printf "%.1f", s/500000}') ns)"
grep -E '^Performance:' prod.log | tail -1 || true

if [[ ! -f prod.gro && ${status} -eq 0 ]]; then
  sbatch -J "${SLURM_JOB_NAME}" --export=T="${T}" --dependency=singleton "$0"
fi
exit "${status}"
