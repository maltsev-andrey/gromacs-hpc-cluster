#!/usr/bin/env bash
# Phase 5, step 2: make the protein whole across periodic boundaries.
#
# GROMACS wraps coordinates into the box as it writes them. Over 200 ns the
# protein diffuses across a face of the dodecahedron, and from then on its
# atoms are split between opposite sides of the box. Anything that measures
# distances between atoms of the same molecule without minimum-image
# convention then measures a protein torn in half:
#
#   before: Rg median 2.62 nm, RMSD median 2.40 nm   (nonsense)
#   after:  Rg around 1.40 nm, RMSD a few tenths     (lysozyme)
#
# H-bond counts were never affected, because those use minimum-image distances.
#
# -pbc whole uses the bond topology in protein.tpr to reassemble each molecule.
# The result may stick out of the box, which is correct and harmless: Rg is
# translation invariant and MDAnalysis superimposes before measuring RMSD and
# RMSF, so absolute position does not matter.
#
# Run on the head node with GMXRC sourced. Reads prod.xtc, writes prod_whole.xtc
# beside it, and never modifies the original.

set -euo pipefail

# ---------- settings ----------
LYSO="/nfs/shared/gromacs/lyso"
TOPOL="${LYSO}/common/protein.tpr"     # protein-only, 1960 atoms
GMXRC="/nfs/shared/apps/gromacs/2026.2-cpu-avx2_256/bin/GMXRC"
TEMPS="${TEMPS:-300 325 350 400 450}"
GROUP="System"                          # protein.tpr holds protein only
# ------------------------------

set +u
# shellcheck source=/dev/null
source "${GMXRC}"
set -u
export GMX_MAXBACKUP=-1

[[ -f "${TOPOL}" ]] || { echo "missing ${TOPOL}: run prep_analysis.sh first"; exit 1; }

for T in ${TEMPS}; do
  in="${LYSO}/T${T}/run/prod.xtc"
  out="${LYSO}/T${T}/run/prod_whole.xtc"
  if [[ ! -f "${in}" ]]; then echo "T${T}: no ${in}, skipped"; continue; fi
  if [[ -f "${out}" ]]; then echo "T${T}: ${out} exists, skipped"; continue; fi
  echo "== T${T}"
  printf '%s\n' "${GROUP}" | gmx trjconv -s "${TOPOL}" -f "${in}" -o "${out}" -pbc whole
done

echo
echo "== check: Rg of the first frame, before and after"
for T in ${TEMPS}; do
  [[ -f "${LYSO}/T${T}/run/prod_whole.xtc" ]] || continue
  printf '  T%-4s ' "${T}"
  printf 'before '; printf '%s\n' "${GROUP}" | gmx gyrate -s "${TOPOL}" \
      -f "${LYSO}/T${T}/run/prod.xtc" -o /tmp/g_before_${T}.xvg -b 0 -e 0 2>/dev/null \
      | grep -oE 'Radius of gyration[^(]*' | head -1 || true
  printf '        after  '; printf '%s\n' "${GROUP}" | gmx gyrate -s "${TOPOL}" \
      -f "${LYSO}/T${T}/run/prod_whole.xtc" -o /tmp/g_after_${T}.xvg -b 0 -e 0 2>/dev/null \
      | grep -oE 'Radius of gyration[^(]*' | head -1 || true
done

cat <<'NOTE'

Expect roughly 1.4 nm after. If a value is still far above that, the molecule
is not being reassembled: check that TOPOL really is the protein-only tpr and
that its atom count matches the trajectory.

Then rerun the analysis; it reads prod_whole.xtc now.
NOTE
