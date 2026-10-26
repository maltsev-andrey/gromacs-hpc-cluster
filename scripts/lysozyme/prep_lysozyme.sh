#!/usr/bin/env bash
# Phase 4, step 1: build the lysozyme system (PDB 1AKI) and the production
# inputs for every temperature in the ladder. Run ONCE on gpu-host
# (internet for the PDB download, GPU for the short equilibration), under
# the card-B scheduler:
#
#   GPU_RUN_STOP_SERVICES="$MY_SERVICES" gpu-run.sh lyso-prep -- ./prep_lysozyme.sh
#
# Output: /nfs/shared/gromacs/lyso/T<K>/prod.tpr for each temperature, plus
# the shared topology and equilibrated structure in /nfs/shared/gromacs/lyso/common.
# The same .tpr files run on card B and on the compute nodes.
#
# Pipeline (all at 300 K until the ladder):
#   pdb2gmx (AMBER99SB-ILDN, TIP3P) -> dodecahedron box, 1.0 nm margin
#   -> solvate -> 0.15 M NaCl -> EM -> NVT 100 ps -> NPT 500 ps
#   -> per temperature: NVT 1 ns at T (velocities regenerated at T) -> prod.tpr

set -euo pipefail

# ---------- settings ----------
GMXRC="/opt/gromacs/2026.2-cuda/bin/GMXRC"
MDP_DIR="/nfs/shared/projects/gromacs/lyso/mdp"
OUT="/nfs/shared/gromacs/lyso"
WORK="/data/gromacs/lyso_prep"              # local NVMe scratch
PDB_ID="1AKI"
FORCEFIELD="amber99sb-ildn"
WATER="tip3p"
BOX_MARGIN="1.0"                            # nm, protein to box edge
SALT="0.15"                                 # mol/l NaCl
TEMPS="${TEMPS:-300 325 350 400 450}"       # K
NTOMP=4                                     # resident mode needs no more
GPU_UUID="GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
# ------------------------------

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export CUDA_VISIBLE_DEVICES="${GPU_UUID}"
export GMX_MAXBACKUP=-1                     # no #file.1# backups

mkdir -p "${WORK}" "${OUT}/common"
cd "${WORK}"

log() { echo "$(date -Iseconds) prep: $*"; }

# ---------- 1. structure ----------
if [[ ! -f "${PDB_ID}.pdb" ]]; then
  log "downloading ${PDB_ID}"
  curl -fsSL -o "${PDB_ID}.pdb" "https://files.rcsb.org/download/${PDB_ID}.pdb"
fi
# keep only the protein: drop crystal waters and hetero groups
grep -E '^(ATOM|TER|END)' "${PDB_ID}.pdb" > protein.pdb
log "protein atoms: $(grep -c '^ATOM' protein.pdb)"

# ---------- 2. topology ----------
gmx pdb2gmx -f protein.pdb -o processed.gro -p topol.top -i posre.itp \
  -ff "${FORCEFIELD}" -water "${WATER}" -ignh > pdb2gmx.log 2>&1
grep -E 'Total charge|atoms' pdb2gmx.log | tail -3

# ---------- 3. box, solvent, ions ----------
gmx editconf -f processed.gro -o boxed.gro -c -d "${BOX_MARGIN}" -bt dodecahedron > editconf.log 2>&1
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > solvate.log 2>&1
gmx grompp -f "${MDP_DIR}/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 > grompp_ions.log 2>&1
echo SOL | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL \
  -neutral -conc "${SALT}" > genion.log 2>&1
log "system: $(tail -1 ionized.gro | xargs) box; $(sed -n 2p ionized.gro | xargs) atoms"

# ---------- 4. minimization ----------
gmx grompp -f "${MDP_DIR}/em.mdp" -c ionized.gro -p topol.top -o em.tpr > grompp_em.log 2>&1
gmx mdrun -deffnm em -ntmpi 1 -ntomp "${NTOMP}" -nb gpu > mdrun_em.log 2>&1
grep -E 'Potential Energy|Maximum force' em.log | tail -2

# ---------- 5. equilibration at 300 K ----------
gmx grompp -f "${MDP_DIR}/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr > grompp_nvt.log 2>&1
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp "${NTOMP}" -nb gpu -pme gpu -bonded gpu -update gpu -pin on > mdrun_nvt.log 2>&1
gmx grompp -f "${MDP_DIR}/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr > grompp_npt.log 2>&1
gmx mdrun -deffnm npt -ntmpi 1 -ntomp "${NTOMP}" -nb gpu -pme gpu -bonded gpu -update gpu -pin on > mdrun_npt.log 2>&1
echo Density | gmx energy -f npt.edr -o density.xvg > energy_density.log 2>&1
log "NPT density (last 100 ps): $(tail -50 density.xvg | awk '{s+=$2;n++} END {printf "%.1f kg/m3", s/n}')"

cp topol.top posre.itp npt.gro npt.cpt "${OUT}/common/"
cp "${MDP_DIR}"/*.mdp "${OUT}/common/"

# ---------- 6. temperature ladder ----------
for T in ${TEMPS}; do
  d="${OUT}/T${T}"; mkdir -p "${d}"
  sed "s/__TEMP__/${T}/g" "${MDP_DIR}/eq_nvt_T.mdp" > "${d}/eq_nvt.mdp"
  sed "s/__TEMP__/${T}/g" "${MDP_DIR}/prod_T.mdp"   > "${d}/prod.mdp"
  # 1 ns NVT at T from the 300 K NPT box, velocities regenerated at T
  gmx grompp -f "${d}/eq_nvt.mdp" -c npt.gro -r npt.gro -p topol.top -o "${d}/eq_nvt.tpr" > "${d}/grompp_eq.log" 2>&1
  ( cd "${d}" && gmx mdrun -deffnm eq_nvt -ntmpi 1 -ntomp "${NTOMP}" -nb gpu -pme gpu -bonded gpu -update gpu -pin on > mdrun_eq.log 2>&1 )
  gmx grompp -f "${d}/prod.mdp" -c "${d}/eq_nvt.gro" -t "${d}/eq_nvt.cpt" -p topol.top -o "${d}/prod.tpr" > "${d}/grompp_prod.log" 2>&1
  log "T=${T} K: prod.tpr ready ($(grep -c WARNING "${d}/grompp_prod.log" || true) warnings)"
done

log "done. inputs in ${OUT}/T*/prod.tpr"
