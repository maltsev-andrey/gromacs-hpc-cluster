#!/usr/bin/env bash
# Phase 7: extend a finished 200 ns replica on a GPU, in bounded windows.
#
# Phase 6 showed the 200 ns runs had not reached a plateau (T300 was still
# rising when it stopped). Every run has a checkpoint, so this continues it
# instead of starting again: same trajectory, same files, more time.
#
#   RUN=T300 UNTIL_NS=600 HOURS=8 ./extend_lyso_gpu.sh            # card B (P100), via gpu-run.sh
#   RUN=T400 UNTIL_NS=400 HOURS=8 CARD=A ./extend_lyso_gpu.sh     # card A (V100), no arbiter
#
# Repeat (cron or by hand) until it prints "complete". Each call continues from
# the checkpoint for at most HOURS, then syncs back to NFS.
#
# CARD=B goes through gpu-run.sh as every other batch job on that card does.
# CARD=A runs beside the RAG stack with NO arbiter: Kubernetes does not know
# this process exists. That is a deliberate, time-boxed exception to the
# one-arbiter-per-card rule; do not leave it in cron after Phase 7.

set -euo pipefail

# ---------- settings ----------
RUN="${RUN:?set RUN, e.g. RUN=T300}"
UNTIL_NS="${UNTIL_NS:?set UNTIL_NS, the new total length, e.g. 600}"
HOURS="${HOURS:-8}"
CARD="${CARD:-B}"                           # B = P100 (scheduled), A = V100 (shared with RAG)
GMXRC="/opt/gromacs/2026.2-cuda/bin/GMXRC"
LYSO="/nfs/shared/gromacs/lyso"
SCRATCH="/data/gromacs/lyso_runs"           # local NVMe: trajectory writes stay off NFS
GPU_RUN="/usr/local/bin/gpu-run.sh"
NVSMI="/usr/bin/nvidia-smi"
NTOMP=4
PINOFFSET_B=0                               # two mdruns must not pin to the same cores
PINOFFSET_A=8                               # stride is 2 with hyperthreading: 0 -> cores 0-3, 8 -> cores 4-7
STEPS_PER_NS=500000
# ------------------------------

self="$(readlink -f "$0")"                  # absolute, so the re-exec works from cron or any cwd
src="${LYSO}/${RUN}/run"
run="${SCRATCH}/${RUN}"
target=$(( UNTIL_NS * STEPS_PER_NS ))
[[ -f "${src}/prod.cpt" ]] || { echo "no checkpoint at ${src}/prod.cpt"; exit 1; }

if [[ "${CARD}" == "B" && -z "${GPU_RUN_NAME:-}" ]]; then
  exec env GPU_RUN_STOP_SERVICES="${GPU_RUN_STOP_SERVICES:-}" "${GPU_RUN}" "lyso-ext-${RUN}" -- \
       env GPU_RUN_NAME="lyso-ext-${RUN}" RUN="${RUN}" UNTIL_NS="${UNTIL_NS}" HOURS="${HOURS}" CARD=B /bin/bash "${self}"
fi

case "${CARD}" in
  B) model="P100"; pinoffset="${PINOFFSET_B}" ;;
  A) model="V100"; pinoffset="${PINOFFSET_A}" ;;
  *) echo "CARD must be A or B"; exit 2 ;;
esac
uuid=$("${NVSMI}" --query-gpu=uuid,name --format=csv,noheader | awk -F', ' -v m="${model}" '$2 ~ m {print $1; exit}')
[[ -n "${uuid}" ]] || { echo "no ${model} found by nvidia-smi"; exit 1; }

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export CUDA_VISIBLE_DEVICES="${uuid}"
export GMX_MAXBACKUP=-1

# first call: bring the finished run to local scratch, exactly as it stands
mkdir -p "${run}"
if [[ ! -f "${run}/prod.cpt" ]]; then
  echo "$(date -Iseconds) ${RUN}: staging from ${src}"
  rsync -a "${src}/prod.tpr" "${src}/prod.cpt" "${src}/prod.log" "${src}/prod.edr" "${src}/prod.xtc" "${run}/"
fi
cd "${run}"

done_steps=$(grep -oE 'Writing checkpoint, step [0-9]+' prod.log | tail -1 | grep -oE '[0-9]+$' || echo 0)
if (( done_steps >= target )) && [[ -f prod.gro ]]; then
  echo "${RUN}: complete at ${UNTIL_NS} ns"; exit 0
fi

# a tpr that runs to the new end time; the original prod.tpr stays untouched
if [[ ! -f "prod_until${UNTIL_NS}.tpr" ]]; then
  gmx convert-tpr -s prod.tpr -until $(( UNTIL_NS * 1000 )) -o "prod_until${UNTIL_NS}.tpr" > convert_tpr.log 2>&1
fi

echo "$(date -Iseconds) ${RUN}: card ${CARD} (${model}), ${HOURS} h window, from step ${done_steps} toward ${target}"
status=0
gmx mdrun -s "prod_until${UNTIL_NS}.tpr" -deffnm prod -cpi prod.cpt -append \
  -ntmpi 1 -ntomp "${NTOMP}" -pin on -pinoffset "${pinoffset}" \
  -nb gpu -pme gpu -bonded gpu -update gpu -maxh "${HOURS}" >> mdrun_ext.out 2>&1 || status=$?

done_steps=$(grep -oE 'Writing checkpoint, step [0-9]+' prod.log | tail -1 | grep -oE '[0-9]+$' || echo 0)
echo "$(date -Iseconds) ${RUN}: mdrun exit ${status}, step ${done_steps} of ${target} ($(awk -v s="${done_steps}" -v n="${STEPS_PER_NS}" 'BEGIN{printf "%.1f", s/n}') ns)"
grep -E '^Performance:' prod.log | tail -1 || true

# sync back; the PBC-corrected copy on NFS is now stale, so remove it and let
# make_whole.sh rebuild it from the longer trajectory
rsync -a prod.cpt prod.log prod.edr prod.xtc "${src}/"
[[ -f prod.gro ]] && (( done_steps >= target )) && rsync -a prod.gro "${src}/prod.gro"
rm -f "${src}/prod_whole.xtc"
exit "${status}"
