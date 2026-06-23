#!/usr/bin/env bash
#
# sppComparison.sh
# Equal-sample-count (SPP) comparison: BTC vs Mueller (vanilla guided_path).
# One render per (scene, spp, integrator) at 63 / 255 / 1023 spp.
#
# Scenes: kitchen, spaceship, mis (veach_mi), torus.
# For every run the script rewrites the base scene's <integrator> and <sampler>
# blocks, so the only thing that differs between BTC and Mueller is the
# directional structure (identical budget, sampler, depth, strictNormals).
#
# Run INSIDE the btc-pg container, from /home/mitsuba:
#   bash testing/sppComparison.sh 2>&1 | tee sppComparison_progress.log
# or detached from the host:
#   docker exec -d -w /home/mitsuba btc-pg bash -lc \
#     'bash testing/sppComparison.sh > sppComparison_progress.log 2>&1'
#
# Resumable: any (scene,integrator,spp) whose .exr already exists is skipped.
# Rebuild the plugin first (source setpath.sh && scons) so the run uses the
# current btc_guided_path build.

set -u

# ----------------------------- configuration -----------------------------
SPPS=(63 255 1023)
SCENES=(kitchen torus mis)

# scene label -> base scene xml (provides geometry + sensor; its integrator
# and sampler blocks are overwritten per run).
declare -A BASE=(
  [kitchen]=/home/scenes/kitchen/kitchen-btc-63.xml
  [spaceship]=/home/scenes/spaceship/spaceship.xml
  [mis]=/home/scenes/veach_mi/mi.xml
  [torus]=/home/scenes/torus/torus.xml
)

# optional per-scene maxDepth (empty = unbounded). Applied to BOTH integrators
# so the comparison stays controlled. spaceship is a deep scene; keep its bound.
declare -A MAXDEPTH=(
  [spaceship]=10
)

# integrator label -> mitsuba plugin
declare -A PLUGIN=(
  [btc]=btc_guided_path
  [muller]=guided_path
)
INTEGRATORS=(btc muller)

# BTC directional knobs. Final-pass choice: conservative t-test (P999) + a high
# cap that rarely binds (verify via "Saturated cells" in render.log).
BTC_EAGERNESS="0"     # 0 = P999 (most conservative)
BTC_MAXSPLITS="4096"  # directional tile split cap

OUTDIR=/home/mitsuba/sppComparison_out
GENPREFIX=_sppgen_    # generated scenes are written next to their base with this prefix
# -------------------------------------------------------------------------

mkdir -p "$OUTDIR"

# Make sure the mitsuba binary is on PATH. setpath.sh references variables that
# may be unset, so disable nounset while sourcing it.
if ! command -v mitsuba >/dev/null 2>&1; then
  set +u
  # shellcheck disable=SC1091
  source /home/mitsuba/setpath.sh >/dev/null 2>&1
  set -u
fi
if ! command -v mitsuba >/dev/null 2>&1; then
  echo "ERROR: mitsuba not on PATH (setpath.sh failed?)"; exit 1
fi

# Rewrite a base scene's integrator + sampler into a self-contained spp scene.
# args: base_xml out_xml plugin spp maxdepth eagerness maxsplits
gen_scene() {
  python3 - "$@" <<'PY'
import re, sys
base, out, plugin, spp, maxdepth, eagerness, maxsplits = sys.argv[1:8]
src = open(base).read()

lines = ['\t\t<boolean name="strictNormals" value="true"/>',
         '\t\t<string name="distribution" value="radiance"/>',
         '\t\t<string name="budgetType" value="spp"/>',
         '\t\t<boolean name="dumpSDTree" value="false"/>',
         f'\t\t<float name="budget" value="{spp}"/>']
if maxdepth:
    lines.append(f'\t\t<integer name="maxDepth" value="{maxdepth}"/>')
if plugin == 'btc_guided_path':
    if eagerness:
        lines.append(f'\t\t<integer name="eagerness" value="{eagerness}"/>')
    if maxsplits:
        lines.append(f'\t\t<integer name="maxSplits" value="{maxsplits}"/>')
integ = '\t<integrator type="%s">\n%s\n\t</integrator>' % (plugin, '\n'.join(lines))

# Replace the first integrator block, or insert one if the scene lacks it.
if re.search(r'<integrator\b', src):
    src = re.sub(r'<integrator\b.*?</integrator>', integ, src, count=1, flags=re.S)
else:
    src = re.sub(r'(<scene\b[^>]*>)', r'\1\n' + integ, src, count=1)

# Force an independent sampler at the requested count (mi.xml ships ldsampler,
# which asserts under PPG's pass-based rendering).
sampler = ('\t\t<sampler type="independent">\n'
           '\t\t\t<integer name="sampleCount" value="%s"/>\n'
           '\t\t\t<integer name="seed" value="0"/>\n'
           '\t\t</sampler>') % spp
src = re.sub(r'<sampler\b.*?</sampler>', sampler, src, flags=re.S)

open(out, 'w').write(src)
PY
}

echo "=== SPP comparison started: $(date) ==="
echo "scenes: ${SCENES[*]} | spp: ${SPPS[*]} | integrators: ${INTEGRATORS[*]}"
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

  for spp in "${SPPS[@]}"; do
    for integ in "${INTEGRATORS[@]}"; do
      plugin="${PLUGIN[$integ]}"
      rundir="$OUTDIR/${scene}-${integ}-${spp}spp"
      out="$rundir/${scene}-${integ}-${spp}spp.exr"
      log="$rundir/render.log"
      tag="${scene}/${integ}/${spp}spp"

      if [[ -f "$out" ]]; then
        echo "[skip] $tag (output exists)"
        continue
      fi
      mkdir -p "$rundir"

      gen="$scenedir/${GENPREFIX}${scene}-${integ}-${spp}.xml"
      gen_scene "$base" "$gen" "$plugin" "$spp" "$md" "$BTC_EAGERNESS" "$BTC_MAXSPLITS"

      echo "[run ] $tag -> $(basename "$out")"
      start=$(date +%s)

      # Render in the background and sample peak RSS while it runs.
      mitsuba "$gen" -o "$out" >"$log" 2>&1 &
      mpid=$!
      peak=0
      while kill -0 "$mpid" 2>/dev/null; do
        rss=$(ps -o rss= -C mitsuba 2>/dev/null | sort -n | tail -1)
        [[ -n "$rss" && "$rss" -gt "$peak" ]] && peak=$rss
        sleep 2
      done
      wait "$mpid"; rc=$?
      end=$(date +%s)

      if [[ $rc -eq 0 && -f "$out" ]]; then
        var=$(grep -oE 'Var: [0-9.]+' "$log" | tail -1)
        echo "       done ${tag} | $((end-start))s | peakRSS $((peak/1024)) MB | ${var:-Var: n/a}"
      else
        echo "       FAILED ${tag} (rc=$rc) -- see $log"
      fi
    done
  done
done

echo
echo "=== SPP comparison finished: $(date) ==="
echo "EXRs in $OUTDIR"
echo "Generated scenes are ${GENPREFIX}*.xml inside each scene directory (safe to delete)."
