#!/usr/bin/env bash
#SBATCH --job-name=gmx-node
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/shared/gromacs/results/slurm/node_sweep_%N_%j.out
#
# Phase 2a: GROMACS on ONE compute node (thread-MPI build), both systems,
# sweeping how the 12 vCPUs are split between thread-MPI ranks (domain
# decomposition) and OpenMP threads. Compare with the GPU host CPU-only
# numbers at 12 threads: benchMEM 17.15, ADH 13.82 ns/day.
#
# Submit:  sbatch slurm_node_sweep.sh            (any free node)
#          sbatch -w node02 slurm_node_sweep.sh
# No RAPL in a VM, so no energy column here.

set -euo pipefail

# ---------- settings ----------
GMXRC="/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/GMXRC"
TPR_MEM="/nfs/shared/gromacs/bench/benchMEM.tpr"
TPR_ADH="/nfs/shared/gromacs/bench/adh_dodec.tpr"
RESULTS="/nfs/shared/gromacs/results"
SYSTEMS="${SYSTEMS:-mem adh}"
LAYOUTS="${LAYOUTS:-1x12 2x6 3x4 4x3 6x2 12x1}"   # ntmpi x ntomp
REPEATS="${REPEATS:-3}"
NSTEPS=20000
RESET=10000
# ------------------------------

node=$(hostname -s)
OUT_ROOT="/nfs/shared/gromacs/runs_cpu/${node}"
CSV="${RESULTS}/bench_node_${node}.csv"
mkdir -p "${OUT_ROOT}" "${RESULTS}/slurm"

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u

{
  echo "== ${node} $(date -Iseconds) job ${SLURM_JOB_ID}"
  lscpu | grep -E 'Model name|^CPU\(s\)|Thread\(s\) per core|Core\(s\) per socket|Flags' | sed 's/Flags:.*avx2.*/Flags: ... avx2 present/'
  free -g | head -2
  gmx --version | grep -E 'GROMACS version|SIMD|MPI library'
} | tee "${OUT_ROOT}/node_info_${SLURM_JOB_ID}.txt"

[[ -f "${CSV}" ]] || echo "date,node,system,ntmpi,ntomp,repeat,ns_per_day,wall_s,npme,start_epoch" > "${CSV}"

for sys in ${SYSTEMS}; do
  tpr="${TPR_ADH}"; [[ "${sys}" == "mem" ]] && tpr="${TPR_MEM}"
  for lay in ${LAYOUTS}; do
    ntmpi="${lay%x*}"; ntomp="${lay#*x}"
    for r in $(seq 1 "${REPEATS}"); do
      dir="${OUT_ROOT}/${sys}/tmpi${ntmpi}_omp${ntomp}_r${r}"
      mkdir -p "${dir}" && cd "${dir}"
      t0=$(date +%s)
      if gmx mdrun -s "${tpr}" -deffnm bench -ntmpi "${ntmpi}" -ntomp "${ntomp}" -pin on \
           -nsteps "${NSTEPS}" -resetstep "${RESET}" -noconfout -nstlist 100 > mdrun.out 2>&1; then
        nsday=$(awk '/^Performance:/ {print $2}' bench.log)
        wall=$(awk '/^ *Time:/ {print $3}' bench.log)
        npme=$(grep -oE 'Using [0-9]+ separate PME ranks' bench.log | grep -oE '[0-9]+' || echo 0)
        echo "$(date -Iseconds),${node},${sys},${ntmpi},${ntomp},${r},${nsday},${wall},${npme},${t0}" | tee -a "${CSV}"
      else
        echo "FAIL ${sys} ${lay} r${r} (see ${dir}/mdrun.out)"
      fi
    done
  done
done
echo "Done -> ${CSV}"
