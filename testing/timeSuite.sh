#!/usr/bin/env bash
#
# timeSuite.sh
# Equal-time comparison: BTC vs Mueller (vanilla guided_path), budgetType=seconds.
# Scenes: kitchen, mis (veach_mi), torus.  spaceship is excluded (its runs error).
# Budgets: 5 / 20 / 60 min.  One run per (scene, budget, integrator).
#
# Like sppComparison.sh: rewrites each base scene's integrator + sampler so the
# only difference between BTC and Mueller is the directional structure. No SD-tree
# dump (dumpSDTree=false); structure-size stats (spatial cells / total tiles /
# saturated cells / memory) are logged per iteration into each run's render.log.
# 3 seeds per config for equal-time error bars.
#
# Run INSIDE the btc-pg container, from /home/mitsuba. Source setpath.sh first so
# the script inherits PATH:
#   cd /home/mitsuba && source setpath.sh && \
#     nohup bash testing/timeSuite.sh >> timeSuite_progress.log 2>&1 &
#
# Resumable: a run whose .exr exists is skipped; an interrupted run (no .exr) has
# its folder wiped and re-runs from scratch. Only one render runs at a time
# (equal-time fairness), so the total wall-clock is roughly the sum of budgets.

set -u

# ----------------------------- configuration -----------------------------
# budget label -> seconds
declare -A BUDGET=( [5min]=300 [20min]=1200 [60min]=3600 )

# Per-scene budgets (only kitchen gets the 60min point).
declare -A SCENE_BUDGETS=(
  [kitchen]="5min 20min 60min"
  [mis]="5min 20min"
  [torus]="5min 20min"
)

SCENES=(kitchen mis torus)        # spaceship excluded (broken)

declare -A BASE=(
  [kitchen]=/home/scenes/kitchen/kitchen-btc-63.xml
  [mis]=/home/scenes/veach_mi/mi.xml
  [torus]=/home/scenes/torus/torus.xml
)

# optional per-scene maxDepth (empty = unbounded). Applied to both integrators.
declare -A MAXDEPTH=( )

declare -A PLUGIN=( [btc]=btc_guided_path [muller]=guided_path )
INTEGRATORS=(btc muller)

# Seeds for equal-time error bars. Each seed is an independent noise realization
# (requires the seed-aware independent sampler). NOTE: wall-clock scales linearly
# with the seed count -- 3 seeds ~= 3x the runtime below. Trim scenes/budgets or
# set SEEDS=(0) for a first pass if a full 3-seed sweep won't fit the night.
SEEDS=(0 1 2)

# BTC directional knobs. Final-pass choice: conservative t-test (P999) + a high
# cap that rarely binds (verify via "Saturated cells" in render.log).
BTC_EAGERNESS="0"
BTC_MAXSPLITS="4096"

# Per-pass sampler size. In seconds mode the integrator drives passes by time,
# so this does not affect the total sample count; it just needs to be >= 1.
SAMPLECOUNT=1024

OUTDIR=/home/mitsuba/timeSuite_out
GENPREFIX=_timegen_
# -------------------------------------------------------------------------

mkdir -p "$OUTDIR"

# setpath.sh references vars that may be unset; source it with nounset off.
if ! command -v mitsuba >/dev/null 2>&1; then
  set +u
  # shellcheck disable=SC1091
  source /home/mitsuba/setpath.sh >/dev/null 2>&1
  set -u
fi
if ! command -v mitsuba >/dev/null 2>&1; then
  echo "ERROR: mitsuba not on PATH (setpath.sh failed?)"; exit 1
fi

# Rewrite a base scene's integrator + sampler into a time-budget scene.
# args: base out plugin seconds maxdepth eagerness maxsplits samplecount seed
gen_scene() {
  python3 - "$@" <<'PY'
import re, sys
base, out, plugin, seconds, maxdepth, eagerness, maxsplits, samplecount, seed = sys.argv[1:10]
src = open(base).read()

lines = ['\t\t<boolean name="strictNormals" value="true"/>',
         '\t\t<string name="distribution" value="radiance"/>',
         '\t\t<string name="budgetType" value="seconds"/>',
         '\t\t<boolean name="dumpSDTree" value="false"/>',
         f'\t\t<float name="budget" value="{seconds}"/>']
if maxdepth:
    lines.append(f'\t\t<integer name="maxDepth" value="{maxdepth}"/>')
if plugin == 'btc_guided_path':
    if eagerness:
        lines.append(f'\t\t<integer name="eagerness" value="{eagerness}"/>')
    if maxsplits:
        lines.append(f'\t\t<integer name="maxSplits" value="{maxsplits}"/>')
integ = '\t<integrator type="%s">\n%s\n\t</integrator>' % (plugin, '\n'.join(lines))

if re.search(r'<integrator\b', src):
    src = re.sub(r'<integrator\b.*?</integrator>', integ, src, count=1, flags=re.S)
else:
    src = re.sub(r'(<scene\b[^>]*>)', r'\1\n' + integ, src, count=1)

sampler = ('\t\t<sampler type="independent">\n'
           '\t\t\t<integer name="sampleCount" value="%s"/>\n'
           '\t\t\t<integer name="seed" value="%s"/>\n'
           '\t\t</sampler>') % (samplecount, seed)
src = re.sub(r'<sampler\b.*?</sampler>', sampler, src, flags=re.S)

open(out, 'w').write(src)
PY
}

echo "=== TIME suite started: $(date) ==="
echo "scenes: ${SCENES[*]} | per-scene budgets: kitchen[${SCENE_BUDGETS[kitchen]}] mis[${SCENE_BUDGETS[mis]}] torus[${SCENE_BUDGETS[torus]}] | integrators: ${INTEGRATORS[*]}"
echo "output: $OUTDIR"
echo

for scene in "${SCENES[@]}"; do
  base="${BASE[$scene]}"
  if [[ ! -f "$base" ]]; then
    echo "[skip] $scene: base scene not found ($base)"
    continue
  fi
  scenedir=$(dirname "$base")
  md="${MAXDEPTH[$scene]:-}"

  for blabel in ${SCENE_BUDGETS[$scene]}; do
    secs="${BUDGET[$blabel]}"
    for integ in "${INTEGRATORS[@]}"; do
      plugin="${PLUGIN[$integ]}"
      for seed in "${SEEDS[@]}"; do
        rundir="$OUTDIR/${scene}-${integ}-${blabel}-s${seed}"
        out="$rundir/${scene}-${integ}-${blabel}-s${seed}.exr"
        log="$rundir/render.log"
        tag="${scene}/${integ}/${blabel}/s${seed}"

        if [[ -f "$out" ]]; then
          echo "[skip] $tag (output exists)"
          continue
        fi
        rm -rf "$rundir"; mkdir -p "$rundir"

        gen="$scenedir/${GENPREFIX}${scene}-${integ}-${blabel}-s${seed}.xml"
        gen_scene "$base" "$gen" "$plugin" "$secs" "$md" "$BTC_EAGERNESS" "$BTC_MAXSPLITS" "$SAMPLECOUNT" "$seed"

        echo "[run ] $tag (${secs}s) -> $(basename "$out")"
        start=$(date +%s)

        mitsuba "$gen" -o "$out" >"$log" 2>&1 &
        mpid=$!
        peak=0
        while kill -0 "$mpid" 2>/dev/null; do
          rss=$(ps -o rss= -C mitsuba 2>/dev/null | sort -n | tail -1)
          [[ -n "$rss" && "$rss" -gt "$peak" ]] && peak=$rss
          sleep 5
        done
        wait "$mpid"; rc=$?
        end=$(date +%s)

        if [[ $rc -eq 0 && -f "$out" ]]; then
          var=$(grep -oE 'Var: [0-9.]+' "$log" | tail -1)
          echo "       done $tag | $((end-start))s | peakRSS $((peak/1024)) MB | ${var:-Var: n/a}"
        else
          echo "       FAILED $tag (rc=$rc) -- see $log"
        fi
      done
    done
  done
done

echo
echo "=== TIME suite finished: $(date) ==="
echo "EXRs in $OUTDIR (per-run folders with .exr + render.log; stats are in render.log)"
