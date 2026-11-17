#!/usr/bin/env bash
# Progress and ETA for every lysozyme replica, Phases 4, 6 and 7 alike.
# Read-only, instant, no sampling.
#
#   ./lyso_eta.sh              # every run it can find
#   ACTIVE=1 ./lyso_eta.sh     # hide the ones already finished
#   RUNS="T350_r2 T350_r3" ./lyso_eta.sh
#
# Run directories are discovered, not listed, so a new replica appears here the
# moment it starts without editing this file.
#
# Rate comes from the wall-clock timestamps GROMACS writes on its checkpoint
# lines, so two of them give an exact rate over hours with nothing to wait for.
# Do NOT try to measure the rate by watching the step counter: nstlog is 50000
# steps (0.1 ns), which at ~27 ns/day is one line every ~320 s, so any sample
# under about 2 h is pure quantisation.
#
# The rate is taken from the most recent pair of checkpoints written in the same
# running stretch (no more than 20 min apart), so a pause between windows never
# leaks into it. "run hours" is how much running time is still needed. A calendar
# finish is only projected for a run whose log is growing right now; a paused run
# (a nightly GPU window, a queued resubmit) says so instead of guessing.

set -uo pipefail

# ---------- settings ----------
LYSO="/nfs/shared/gromacs/lyso"
GPU_SCRATCH="/data/gromacs/lyso_runs"      # preferred automatically when present
TOTAL_NS=200                               # default; extended runs are detected per run
STEPS_PER_NS=500000                        # dt = 0.002 ps
STALE_MIN=20
# ------------------------------

discover() {
  { find "${LYSO}" -maxdepth 2 -type d -name run 2>/dev/null
    find "${GPU_SCRATCH}" -maxdepth 1 -mindepth 1 -type d 2>/dev/null; } \
  | sed -E "s#^${LYSO}/##; s#^${GPU_SCRATCH}/##; s#/run\$##" | sort -u   # one row per run; logfor prefers scratch
}

logfor() {
  local r="${1#@}"
  if [[ "$1" == @* ]]; then echo "${GPU_SCRATCH}/${r}/prod.log"
  elif [[ -f "${GPU_SCRATCH}/${r}/prod.log" ]]; then echo "${GPU_SCRATCH}/${r}/prod.log"
  else echo "${LYSO}/${r}/run/prod.log"; fi
}

ckpts() {
  grep -oE 'Writing checkpoint, step [0-9]+ at .*[0-9]{4}' "$1" 2>/dev/null \
  | sed -E 's/Writing checkpoint, step ([0-9]+) at (.*)/\1\t\2/' \
  | while IFS=$'\t' read -r s d; do
      e=$(date -d "$d" +%s 2>/dev/null) && printf '%s\t%s\n' "$s" "$e"
    done
}

RUNS="${RUNS:-$(discover)}"
[[ -n "${RUNS}" ]] || { echo "no run directories under ${LYSO}"; exit 1; }

now=$(date +%s)
printf '%-12s %8s %6s %9s %8s %10s  %s\n' run "ns done" "pct" "ns/day" "log age" "run hours" "status"
printf '%.0s-' {1..92}; echo

worst=0
for R in ${RUNS}; do
  f=$(logfor "${R}")
  [[ -f "${f}" ]] || continue

  # Target length per run: 200 ns unless an extension tpr (prod_until<N>.tpr,
  # written by extend_lyso_gpu.sh) says otherwise. Those live on GPU scratch.
  tot="${TOTAL_NS}"
  for t in "${GPU_SCRATCH}/${R}"/prod_until*.tpr; do
    [[ -f "$t" ]] || continue
    n_ns="${t##*prod_until}"; n_ns="${n_ns%.tpr}"
    (( n_ns > tot )) && tot="${n_ns}"
  done

  mapfile -t C < <(ckpts "$f")
  n=${#C[@]}
  last_step=0; (( n > 0 )) && last_step="${C[n-1]%%$'\t'*}"
  done_ns=$(awk -v s="${last_step}" -v spn="${STEPS_PER_NS}" 'BEGIN{printf "%.1f", s/spn}')

  if [[ -f "${f%/*}/prod.gro" ]] && (( last_step >= tot * STEPS_PER_NS )); then
    [[ -n "${ACTIVE:-}" ]] && continue
    note="complete"
    (( last_step > tot * STEPS_PER_NS )) && note="extended past ${tot} ns; run this on the GPU host for the target"
    printf '%-12s %8.1f %6s %9s %8s %10s  %s\n' "${R}" "${done_ns}" "-" - - - "${note}"
    continue
  fi
  (( n > 0 )) || { printf '%-12s %8s  %s\n' "${R}" "-" "no checkpoint yet"; continue; }

  # newest pair of checkpoints from one running stretch
  rate=0
  for (( i=n-1; i>=1; i-- )); do
    IFS=$'\t' read -r s2 e2 <<< "${C[i]}"; IFS=$'\t' read -r s1 e1 <<< "${C[i-1]}"
    # a gap means an older session starts here: its rate says nothing about this one
    (( e2 - e1 > 1200 )) && break
    if (( e2 > e1 && s2 > s1 )); then
      rate=$(awk -v a="${s1}" -v b="${s2}" -v t="$(( e2 - e1 ))" -v spn="${STEPS_PER_NS}" 'BEGIN{printf "%.2f", (b-a)/spn/(t/86400)}')
      break
    fi
  done
  mtime=$(stat -c %Y "$f" 2>/dev/null || echo "${now}")

  out=$(awk -v R="${R}" -v sl="${last_step}" -v rate="${rate}" -v now="${now}" -v mt="${mtime}" \
            -v tot="${tot}" -v spn="${STEPS_PER_NS}" -v stale="${STALE_MIN}" 'BEGIN{
    ns = sl/spn; rem = tot-ns; age = (now-mt)/60; if (age < 0) age = 0;
    if (rate <= 0) { printf "0|%-12s %8.1f %5.1f%% %9s %7.0fm %10s  %s\n", R, ns, 100*ns/tot, "-", age, "-", "new session; rate after its next checkpoint"; exit }
    h = rem/rate*24;
    if (age <= stale) {
      cmd = "date -d \"+" int(h*3600) " seconds\" \"+%a %d %b %H:%M\""; cmd | getline eta; close(cmd);
      st = "running; done " eta " if it is not interrupted";
    } else st = "paused; resumes with its next window";
    printf "%.4f|%-12s %8.1f %5.1f%% %9.2f %7.0fm %10.1f  %s\n", h, R, ns, 100*ns/tot, rate, age, h, st }')

  echo "${out#*|}"
  h="${out%%|*}"
  awk -v a="${h}" -v b="${worst}" 'BEGIN{exit !(a>b)}' && worst="${h}"
done

printf '%.0s-' {1..92}; echo
awk -v h="${worst}" 'BEGIN{
  if (h <= 0) { print "Nothing in progress: every run found is either finished or has no rate yet."; exit }
  printf "Most running time still needed by one run: %.1f h\n", h }'

echo
if command -v squeue >/dev/null; then
  echo "Queue:"; squeue -u "$(id -un)" -o '%.8i %.12j %.8T %.10M %.10L %R' | sed 's/^/  /'
else
  echo "squeue not on PATH (a venv is shadowing it); run it in a plain shell"
fi
