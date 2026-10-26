#!/usr/bin/env bash
# GROMACS offload matrix on gpu-host, card B (P100 today).
# v2: runs INSIDE the card-B scheduler. gpu-run.sh queues a reservation,
# evicts ollama, waits for the card to drain and releases on exit;
# GPU_RUN_STOP_SERVICES pauses competing host services for the CPU side.
# This script no longer touches docker or systemd itself.
#
#   full matrix (both systems, ~4 h):
#     GPU_RUN_STOP_SERVICES="$MY_SERVICES" gpu-run.sh gromacs-matrix -- ./bench_matrix.sh
#   reproducibility point (6 runs, ~5 min):
#     SYSTEMS=adh THREADS=12 CONFIGS="nb_pme_gpu full_gpu_resident" CONTROL=0 TAG=repro \
#       GPU_RUN_STOP_SERVICES="$MY_SERVICES" gpu-run.sh gromacs-repro -- ./bench_matrix.sh
#
# Every setting below can be overridden from the environment as shown.
# After a card swap change only the GPU_* block.

set -euo pipefail

# ---------- settings ----------
GMXRC="/opt/gromacs/2026.2-cuda/bin/GMXRC"
TPR_MEM="/data/gromacs/bench/benchMEM.tpr"
TPR_ADH="/data/gromacs/bench/adh/ADH/adh_dodec/adh_dodec.tpr"
RESULTS="/nfs/shared/gromacs/results"

GPU_UUID="GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"   # card B, the scheduled card
GPU_LABEL="${GPU_LABEL:-P100-16GB}"
EXPECT_NAME="${EXPECT_NAME:-P100}"
EXPECT_MIN_MIB="${EXPECT_MIN_MIB:-0}"                  # V100 32GB later: 30000

SYSTEMS="${SYSTEMS:-mem adh}"            # mem | adh
THREADS="${THREADS:-4 8 12 16}"
CONFIGS="${CONFIGS:-all}"                # all | space-separated subset of the labels below
REPEATS="${REPEATS:-3}"
CONTROL="${CONTROL:-1}"                  # 1 = fixed-grid block (-notunepme, 12 threads)
TAG="${TAG:-}"                           # suffix for CSV/log names, e.g. repro
NSTEPS=20000
RESET=10000
# ------------------------------

SUFFIX="${GPU_LABEL}${TAG:+_${TAG}}"
OUT_ROOT="/data/gromacs/bench/runs_${SUFFIX}"
CSV="${RESULTS}/bench_${SUFFIX}.csv"
CPU_LOG="${RESULTS}/cpu_top_${SUFFIX}.log"
RAPL_LOG="${RESULTS}/cpu_rapl_${SUFFIX}.log"
FOREIGN_LOG="${OUT_ROOT}/foreign.log"

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u

BG_PIDS=""
cleanup() { for p in ${BG_PIDS}; do kill "${p}" 2>/dev/null || true; done; }
trap cleanup EXIT
trap 'exit 130' INT TERM

# ---------- guard the card (belt and braces; gpu-run.sh already drained it) ----------
line=$(nvidia-smi --query-gpu=uuid,index,name,memory.total,memory.used --format=csv,noheader,nounits | grep "^${GPU_UUID}" || true)
if [[ -z "${line}" ]]; then echo "ABORT: ${GPU_UUID} not present"; nvidia-smi -L; exit 1; fi
idx=$(echo "${line}"  | awk -F', ' '{print $2}')
name=$(echo "${line}" | awk -F', ' '{print $3}')
mem=$(echo "${line}"  | awk -F', ' '{print $4}')
used=$(echo "${line}" | awk -F', ' '{print $5}')
if [[ "${name}" != *"${EXPECT_NAME}"* ]] || (( mem < EXPECT_MIN_MIB )); then
  echo "ABORT: ${GPU_UUID} is '${name}' ${mem} MiB, expected ${EXPECT_NAME} >=${EXPECT_MIN_MIB} MiB"; exit 1
fi
if (( used > 1000 )); then
  echo "ABORT: ${used} MiB in use on card B. Not running under gpu-run.sh?"
  nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory --format=csv,noheader | grep "${GPU_UUID}" || true
  exit 1
fi
export CUDA_VISIBLE_DEVICES="${GPU_UUID}"
echo "card B: index ${idx}, ${name}, ${mem} MiB, ${used} MiB used -> CUDA_VISIBLE_DEVICES=${GPU_UUID}"

mkdir -p "${OUT_ROOT}" "${RESULTS}"
nvidia-smi -L > "${OUT_ROOT}/nvidia-smi-L.txt"
nvidia-smi --query-gpu=index,name,pci.bus_id,uuid,memory.total,driver_version --format=csv > "${OUT_ROOT}/gpus.csv"
[[ -f "${CSV}" ]] || echo "date,gpu,config,ntomp,repeat,ns_per_day,avg_power_w,wh_per_ns,start_epoch" > "${CSV}"

# ---------- host-wide loggers ----------
(
  while true; do
    echo "=== $(date -Iseconds)"
    top -b -n2 -d1 -o %CPU | awk '/^top -/{f++} f==2' | sed -n '8,13p'
    sleep 4
  done
) >> "${CPU_LOG}" 2>&1 &
BG_PIDS="${BG_PIDS} $!"

if sudo -n cat /sys/class/powercap/intel-rapl:0/energy_uj >/dev/null 2>&1; then
  # written as the invoking user on purpose (results dir on NFS)
  # shellcheck disable=SC2024
  sudo -n bash -c 'while true; do echo "$(date +%s),$(cat /sys/class/powercap/intel-rapl:0/energy_uj),"; sleep 5; done' \
    >> "${RAPL_LOG}" 2>/dev/null &
  BG_PIDS="${BG_PIDS} $!"
  echo "RAPL logger -> ${RAPL_LOG}"
else
  echo "NOTE: no passwordless sudo; run the RAPL logger loop yourself into ${RAPL_LOG}"
fi

# ---------- build the run list ----------
want() { [[ "${CONFIGS}" == "all" ]] || [[ " ${CONFIGS} " == *" $1 "* ]]; }
RUNS=()
for sys in ${SYSTEMS}; do
  cfgs=(
    "cpu_only|-nb cpu -pme cpu -bonded cpu -update cpu"
    "nb_gpu|-nb gpu -pme cpu -bonded cpu -update cpu"
    "nb_pme_gpu|-nb gpu -pme gpu -bonded cpu -update cpu"
    "nb_pme_bonded_gpu|-nb gpu -pme gpu -bonded gpu -update cpu"
  )
  [[ "${sys}" == "adh" ]] && cfgs+=("full_gpu_resident|-nb gpu -pme gpu -bonded gpu -update gpu")
  for cfg in "${cfgs[@]}"; do
    want "${cfg%%|*}" || continue
    for nt in ${THREADS}; do RUNS+=("${sys}|${cfg%%|*}|${nt}|${cfg#*|}"); done
  done
  if (( CONTROL == 1 )); then
    RUNS+=("${sys}|nb_gpu_notune|12|-nb gpu -pme cpu -bonded cpu -update cpu -notunepme")
    RUNS+=("${sys}|nb_pme_gpu_notune|12|-nb gpu -pme gpu -bonded cpu -update cpu -notunepme")
    [[ "${sys}" == "adh" ]] && RUNS+=("${sys}|full_gpu_resident_notune|12|-nb gpu -pme gpu -bonded gpu -update gpu -notunepme")
  fi
done
echo "${#RUNS[@]} points x ${REPEATS} repeats on ${GPU_LABEL} -> ${CSV}"

# ---------- run ----------
for spec in "${RUNS[@]}"; do
  IFS='|' read -r sys label nt flags <<< "${spec}"
  tpr="${TPR_ADH}"; [[ "${sys}" == "mem" ]] && tpr="${TPR_MEM}"

  for r in $(seq 1 "${REPEATS}"); do
    dir="${OUT_ROOT}/${sys}/${label}/nt${nt}_r${r}"
    mkdir -p "${dir}" && cd "${dir}"
    t0=$(date +%s)

    nvidia-smi --id="${GPU_UUID}" --query-gpu=power.draw --format=csv,noheader,nounits -lms 500 > power.log &
    pmon=$!
    nvidia-smi --query-compute-apps=timestamp,gpu_uuid,pid,process_name,used_memory \
      --format=csv,noheader -l 2 > gpu_procs.log 2>&1 &
    gmon=$!

    ok=1
    # shellcheck disable=SC2086
    gmx mdrun -s "${tpr}" -deffnm bench -ntmpi 1 -ntomp "${nt}" -pin on ${flags} \
      -nsteps "${NSTEPS}" -resetstep "${RESET}" -noconfout -nstlist 100 > mdrun.out 2>&1 || ok=0

    kill "${pmon}" "${gmon}" 2>/dev/null || true
    if (( ok == 0 )); then echo "FAIL ${sys}:${label} nt${nt} r${r} (see ${dir}/mdrun.out)"; continue; fi

    nsday=$(awk '/^Performance:/ {print $2}' bench.log)
    pw=$(awk '{s+=$1; n++} END {if (n) printf "%.1f", s/n; else print "NA"}' power.log)
    whns=$(awk -v p="${pw}" -v d="${nsday}" 'BEGIN {if (p!="NA" && d>0) printf "%.2f", p*24/d; else print "NA"}')

    foreign=$(grep "${GPU_UUID}" gpu_procs.log 2>/dev/null | grep -v 'gmx' \
      | awk -F', ' '{print $4}' | sort -u | tr '\n' ' ' || true)
    [[ -n "${foreign}" ]] && echo "WARN $(date -Iseconds) ${sys}:${label} nt${nt} r${r}: foreign process on card B: ${foreign}" | tee -a "${FOREIGN_LOG}"

    echo "$(date -Iseconds),${GPU_LABEL},${sys}:${label},${nt},${r},${nsday},${pw},${whns},${t0}" | tee -a "${CSV}"
  done
done

echo "Done -> ${CSV}"
