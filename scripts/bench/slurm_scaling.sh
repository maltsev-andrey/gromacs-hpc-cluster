#!/usr/bin/env bash
#SBATCH --job-name=gmx-scale
#SBATCH --partition=compute
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=/nfs/shared/gromacs/results/slurm/scaling_N%N_%j.out
#
# Phase 2b: multi-node strong scaling with the OpenMPI build (gmx_mpi) over
# the 1 GbE fabric. Submit once per node count; the script reads the
# allocation from Slurm:
#   for n in 1 2 4; do sbatch -N $n slurm_scaling.sh; done
#
# Two rank layouts per node count, both using all 12 vCPUs of every node:
#   12x1  12 MPI ranks per node, 1 OpenMP thread  (max decomposition, max messages)
#    4x3   4 MPI ranks per node, 3 OpenMP threads (fewer, larger messages)
# GROMACS picks separate PME ranks itself; the number lands in the CSV.
#
# No RAPL in a VM, so no energy column.

set -euo pipefail

# ---------- settings ----------
GMXRC="/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/GMXRC"
TPR_MEM="/nfs/shared/gromacs/bench/benchMEM.tpr"
TPR_ADH="/nfs/shared/gromacs/bench/adh_dodec.tpr"
RESULTS="/nfs/shared/gromacs/results"
SYSTEMS="${SYSTEMS:-mem adh}"
LAYOUTS="${LAYOUTS:-12x1 4x3}"          # ranks-per-node x ntomp
REPEATS="${REPEATS:-3}"
NSTEPS=20000
RESET=10000
MPI_IF="10.10.10.0/24"                  # the MPI fabric, not the home LAN
# ------------------------------

nodes="${SLURM_JOB_NUM_NODES}"
OUT_ROOT="/nfs/shared/gromacs/runs_cpu/scaling/N${nodes}"
CSV="${RESULTS}/bench_scaling.csv"
mkdir -p "${OUT_ROOT}" "${RESULTS}/slurm"

set +u
# shellcheck source=/dev/null
source /etc/profile.d/modules.sh 2>/dev/null || true
module load mpi/openmpi-x86_64 2>/dev/null || true
# shellcheck source=/dev/null
source "${GMXRC}"
set -u

echo "== N=${nodes} nodes: ${SLURM_JOB_NODELIST}  job ${SLURM_JOB_ID}  $(date -Iseconds)"
ompi_info 2>/dev/null | head -1 || echo "ompi_info not found"
gmx_mpi --version | grep -E 'GROMACS version|SIMD|MPI library'

[[ -f "${CSV}" ]] || echo "date,nodes,system,ranks_per_node,ntomp,repeat,ns_per_day,wall_s,npme,start_epoch" > "${CSV}"

for sys in ${SYSTEMS}; do
  tpr="${TPR_ADH}"; [[ "${sys}" == "mem" ]] && tpr="${TPR_MEM}"
  for lay in ${LAYOUTS}; do
    rpn="${lay%x*}"; ntomp="${lay#*x}"
    np=$(( rpn * nodes ))
    for r in $(seq 1 "${REPEATS}"); do
      dir="${OUT_ROOT}/${sys}/rpn${rpn}_omp${ntomp}_r${r}"
      mkdir -p "${dir}" && cd "${dir}"
      t0=$(date +%s)
      if mpirun -np "${np}" --map-by "ppr:${rpn}:node:pe=${ntomp}" --bind-to core \
           --mca btl tcp,self --mca btl_tcp_if_include "${MPI_IF}" \
           -x OMP_NUM_THREADS="${ntomp}" \
           gmx_mpi mdrun -s "${tpr}" -deffnm bench -ntomp "${ntomp}" -pin on \
           -nsteps "${NSTEPS}" -resetstep "${RESET}" -noconfout -nstlist 100 > mdrun.out 2>&1; then
        nsday=$(awk '/^Performance:/ {print $2}' bench.log)
        wall=$(awk '/^ *Time:/ {print $3}' bench.log)
        npme=$(grep -oE 'Using [0-9]+ separate PME ranks' bench.log | grep -oE '[0-9]+' || echo 0)
        echo "$(date -Iseconds),${nodes},${sys},${rpn},${ntomp},${r},${nsday},${wall},${npme},${t0}" | tee -a "${CSV}"
      else
        echo "FAIL N${nodes} ${sys} ${lay} r${r} (see ${dir}/mdrun.out)"
      fi
    done
  done
done
echo "Done -> ${CSV}"
