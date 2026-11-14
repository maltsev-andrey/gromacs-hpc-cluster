#!/usr/bin/env bash
# Phase 6 preparation: independent replicas at 325 K and 350 K.
#
# The question Phase 5 left open is whether the RMSD step from 325 K to 350 K
# (+68%, against +7%, +25% and +35% for the other intervals) is real. It sits
# next to lysozyme's experimental melting point, which makes it interesting and
# makes it worth more than one trajectory.
#
# The error bars in the Phase 5 summary are the spread WITHIN one trajectory.
# They are not the right uncertainty for comparing temperatures: consecutive
# frames are correlated and one 200 ns run samples one region of conformational
# space. This phase measures the BETWEEN-replica scatter at both ends instead.
#
# Two new replicas at each temperature, giving n=3 with the Phase 5 runs as r1.
# What makes them independent is fresh velocities from the same equilibrated
# structure, with explicit seeds so the whole thing can be repeated.
#
#   gpu-run.sh lyso-prep6 -- ./prep_phase6.sh       # equilibration on the P100
#
# Equilibration is 1 ns per replica: four minutes each on the card, about an
# hour each on a node, which is why it runs here. Production then goes to the
# nodes, where four replicas run in parallel.

set -euo pipefail

# ---------- settings ----------
LYSO="/nfs/shared/gromacs/lyso"
COMMON="${LYSO}/common"
MDP_DIR="/nfs/shared/projects/gromacs/lyso/mdp"
GMXRC="/opt/gromacs/2026.2-cuda/bin/GMXRC"     # GPU host build
TEMPS="${TEMPS:-325 350}"
REPS="${REPS:-2 3}"                            # r1 is the Phase 5 run
NTOMP="${NTOMP:-4}"
# ------------------------------

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export GMX_MAXBACKUP=-1

log() { echo "$(date -Iseconds) $*"; }

for f in "${COMMON}/npt.gro" "${COMMON}/topol.top" "${MDP_DIR}/eq_nvt_T.mdp" "${MDP_DIR}/prod_T.mdp"; do
  [[ -f "${f}" ]] || { echo "missing ${f}"; exit 1; }
done

for T in ${TEMPS}; do
  for r in ${REPS}; do
    d="${LYSO}/T${T}_r${r}"
    if [[ -f "${d}/prod.tpr" ]]; then log "T${T} r${r}: prod.tpr exists, skipped"; continue; fi
    mkdir -p "${d}"

    # Explicit, distinct, reproducible seed. The Phase 5 runs used gen-seed = -1
    # (random), so r1 cannot be reproduced; from here on the seed is recorded.
    seed=$(( T * 1000 + r ))
    sed -e "s/__TEMP__/${T}/g" -e "s/^gen-seed.*/gen-seed        = ${seed}/" \
        "${MDP_DIR}/eq_nvt_T.mdp" > "${d}/eq_nvt.mdp"
    sed -e "s/__TEMP__/${T}/g" "${MDP_DIR}/prod_T.mdp" > "${d}/prod.mdp"
    grep -q "gen-seed        = ${seed}" "${d}/eq_nvt.mdp" || { echo "seed substitution failed"; exit 1; }

    log "T${T} r${r}: equilibrating 1 ns, gen-seed ${seed}"
    gmx grompp -f "${d}/eq_nvt.mdp" -c "${COMMON}/npt.gro" -r "${COMMON}/npt.gro" \
        -p "${COMMON}/topol.top" -o "${d}/eq_nvt.tpr" > "${d}/grompp_eq.log" 2>&1

    status=0
    ( cd "${d}" && gmx mdrun -deffnm eq_nvt -ntmpi 1 -ntomp "${NTOMP}" -pin on \
        -nb gpu -pme gpu -bonded gpu -update gpu > eq_nvt.out 2>&1 ) || status=$?
    (( status == 0 )) || { echo "equilibration failed for T${T} r${r}, see ${d}/eq_nvt.out"; exit 1; }

    gmx grompp -f "${d}/prod.mdp" -c "${d}/eq_nvt.gro" -t "${d}/eq_nvt.cpt" \
        -p "${COMMON}/topol.top" -o "${d}/prod.tpr" > "${d}/grompp_prod.log" 2>&1
    log "T${T} r${r}: prod.tpr ready ($(grep -c WARNING "${d}/grompp_prod.log" || true) grompp warnings)"
  done
done

echo
echo "== ready to submit =="
for T in ${TEMPS}; do for r in ${REPS}; do
  d="${LYSO}/T${T}_r${r}"
  printf '  %-14s %s\n' "T${T}_r${r}" "$([[ -f ${d}/prod.tpr ]] && echo 'prod.tpr ok' || echo 'MISSING')"
done; done
cat <<'NOTE'

Submit the four production runs, one per node:

  cd /nfs/shared/projects/gromacs/lyso
  for RUN in T325_r2 T325_r3 T350_r2 T350_r3; do
    sbatch -J lyso-$RUN --export=RUN=$RUN slurm_lyso_rep.sh
  done

Each is 200 ns at about 27 ns/day, so roughly 7.5 days, self-resubmitting in
22 h windows. Verify the seeds actually differ before you commit the week:

  grep -h '^gen-seed' /nfs/shared/gromacs/lyso/T3*_r*/eq_nvt.mdp
NOTE
