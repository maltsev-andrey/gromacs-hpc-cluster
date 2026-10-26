#!/usr/bin/env bash
# Phase 4 production on card B, in bounded windows so ollama gets the card back.
# Each call claims the card through the scheduler, continues the run from its
# checkpoint for at most HOURS, then releases. Repeat until prod.gro exists.
#
#   T=450 HOURS=8 ./run_lyso_gpu.sh          # one window, e.g. from a 22:00 cron
#
# Inside, gpu-run.sh does the claim/eviction/release; this script only runs mdrun.
# -maxh stops cleanly at 99% of HOURS, writes a checkpoint, and -cpi -append
# continues next time. Trajectory and energies are appended, not restarted.

set -euo pipefail

# ---------- settings ----------
T="${T:-450}"
HOURS="${HOURS:-8}"
GMXRC="/opt/gromacs/2026.2-cuda/bin/GMXRC"
LYSO="/nfs/shared/gromacs/lyso"
SCRATCH="/data/gromacs/lyso_runs"         # local NVMe: trajectory writes stay off NFS
NTOMP=4
GPU_UUID="GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
# ------------------------------

src="${LYSO}/T${T}"
run="${SCRATCH}/T${T}"
mkdir -p "${run}"
[[ -f "${run}/prod.tpr" ]] || cp "${src}/prod.tpr" "${run}/"

if [[ -f "${run}/prod.gro" ]]; then
  echo "T=${T}: production already complete"; exit 0
fi

if [[ -z "${GPU_RUN_NAME:-}" ]]; then
  # not inside gpu-run.sh yet: re-exec through the scheduler
  exec env GPU_RUN_STOP_SERVICES="${GPU_RUN_STOP_SERVICES:-}" gpu-run.sh "lyso-T${T}" -- \
       env GPU_RUN_NAME="lyso-T${T}" T="${T}" HOURS="${HOURS}" "$0"
fi

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export CUDA_VISIBLE_DEVICES="${GPU_UUID}"
export GMX_MAXBACKUP=-1

cd "${run}"
cpi=""; [[ -f prod.cpt ]] && cpi="-cpi prod.cpt -append"
echo "$(date -Iseconds) T=${T}: window of ${HOURS} h starting ${cpi:+from checkpoint}"
# shellcheck disable=SC2086
gmx mdrun -deffnm prod ${cpi} -ntmpi 1 -ntomp "${NTOMP}" -pin on \
  -nb gpu -pme gpu -bonded gpu -update gpu -maxh "${HOURS}" >> mdrun.out 2>&1
status=$?

done_steps=$(grep -oE 'Writing checkpoint, step [0-9]+' prod.log | tail -1 | grep -oE '[0-9]+$' || echo 0)
echo "$(date -Iseconds) T=${T}: mdrun exit ${status}, at step ${done_steps} of 100000000 ($(awk -v s="${done_steps}" 'BEGIN{printf "%.1f", s/500000}') ns)"
grep -E '^Performance:' prod.log | tail -1 || true

# sync results to NFS after every window
rsync -a --exclude 'prod.trr' "${run}/" "${src}/run/"
exit "${status}"
