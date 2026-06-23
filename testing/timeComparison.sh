#!/bin/bash
# ---------------------------------------------------------------------------
# Equal-time comparison, two scenes in one run (KITCHEN first, then CORNELL):
#   vanilla (Mueller PPG)  vs  BTC (eagerness 4, maxSplits 512)
#   both capped at maxDepth 512 (only the directional structure differs).
#
#   kitchen: 5 / 20 / 60 min      cornell: 5 / 20 min
#
# budgetType="seconds" gives each render that wall-clock; dumps OFF (so all the
# time goes to rendering). Run SEQUENTIALLY -- one render owns the machine at a
# time, which is what makes an equal-time comparison meaningful.
#
# RESUMABLE: any render whose .exr already exists is skipped.
#
# Outputs in timeComparison_out/:
#   kitchen-vanilla-t05min.{exr,log}
#   kitchen-btcEag4Split512-t05min.{exr,log}
#   ... kitchen t20/t60 ... cornell t05/t20 ...
#
# NOTE: requires the BTC plugin built with the eagerness/maxSplits scene props
#       (already built). Total wall-clock ~= 2*(5+20+60) + 2*(5+20) = ~3.7 h.
# ---------------------------------------------------------------------------

cd "$(dirname "$0")"; [ -f setpath.sh ] || cd ..   # run from project root even if invoked from a subdir
source setpath.sh >/dev/null 2>&1
if ! command -v mitsuba >/dev/null 2>&1; then
    echo "ERROR: 'mitsuba' not on PATH. Run from /home/mitsuba after sourcing setpath.sh." >&2
    exit 1
fi

# --- config -----------------------------------------------------------------
OUT="timeComparison_out"
DEPTH=512          # maxDepth, applied to BOTH configs (fairness)
EAG=4              # BTC eagerness
MAXSPLITS=512      # BTC maxSplits

# "tag|scene_dir|vanilla_scene|btc_scene|times_seconds"   (KITCHEN listed first)
SCENESETS=(
    "kitchen|../scenes/kitchen|kitchen-vanilla-63.xml|kitchen-btc-63.xml|300 1200 3600"
    # cornell disabled -- renders black (scene-setup issue, not the integrator). Re-enable when fixed:
    # "cornell|../scenes/cornell-box|cornell-vanilla.xml|cornell-btc.xml|300 1200"
)
# ---------------------------------------------------------------------------

mkdir -p "$OUT"

# render <scene_dir> <base_scene> <budget_seconds> <out_basename>
render() {
    local scene_dir="$1" base="$2" secs="$3" outname="$4"
    local out="${OUT}/${outname}" log="${OUT}/${outname}.log"
    local tmp="${scene_dir}/_tc_${outname}.xml"   # temp scene MUST sit in the scene dir (relative geometry)

    if [ -s "${out}.exr" ]; then
        printf '    SKIP  %-40s (already rendered)\n' "$outname"
        return
    fi

    # Force seconds budget; inject maxDepth for all; add eagerness+maxSplits only for the BTC integrator.
    awk -v secs="$secs" -v depth="$DEPTH" -v eag="$EAG" -v msplit="$MAXSPLITS" '
        /<string name="budgetType"/   { print "\t\t<string name=\"budgetType\" value=\"seconds\"/>"; next }
        /<float name="budget" value=/ { print "\t\t<float name=\"budget\" value=\"" secs "\"/>"; next }
        { print }
        /<integrator type=/ {
            print "\t\t<integer name=\"maxDepth\" value=\"" depth "\"/>"
            if ($0 ~ /btc_guided_path/) {
                print "\t\t<integer name=\"eagerness\" value=\"" eag "\"/>"
                print "\t\t<integer name=\"maxSplits\" value=\"" msplit "\"/>"
            }
        }
    ' "${scene_dir}/${base}" > "$tmp"

    echo ">>> [$(date +%T)] $outname  (budget ${secs}s = $((secs/60))m)"
    local RSTART; RSTART=$(date +%s)
    mitsuba "$tmp" -o "$out" > "$log" 2>&1 &
    local MPID=$! PEAK=0 R
    while kill -0 "$MPID" 2>/dev/null; do
        R=$(awk '/VmRSS/{print $2}' /proc/$MPID/status 2>/dev/null)
        [ -n "$R" ] && [ "$R" -gt "$PEAK" ] && PEAK=$R
        sleep 3
    done
    wait "$MPID"; local RC=$?
    rm -f "$tmp"

    local RT=$(( $(date +%s) - RSTART ))
    local VAR; VAR=$(grep -oE 'Var: [0-9.eE+-]+' "$log" | tail -1 | awk '{print $2}')
    local AP;  AP=$(grep -i "Average path length" "$log" | tail -1 | sed 's/.*Average path length *//')
    printf '    done  %-40s wall=%dm%02ds peakRSS=%dMB rc=%s finalVar=%s pathLen=%s\n' \
        "$outname" $((RT/60)) $((RT%60)) $((PEAK/1024)) "$RC" "${VAR:-n/a}" "${AP:-n/a}"
}

echo "=== time comparison START $(date) | BTC: eag=$EAG maxSplits=$MAXSPLITS | maxDepth=$DEPTH ==="

for set in "${SCENESETS[@]}"; do
    IFS='|' read -r tag scene_dir van btc times <<< "$set"
    echo "######## $tag ########"
    for secs in $times; do
        m=$(printf '%02d' $((secs/60)))
        echo "--- ${tag} ${secs}s (${m}min) ---"
        render "$scene_dir" "$van" "$secs" "${tag}-vanilla-t${m}min"
        render "$scene_dir" "$btc" "$secs" "${tag}-btcEag${EAG}Split${MAXSPLITS}-t${m}min"
    done
done

echo
echo "=== SUMMARY (finalVar; lower = less noise = better at that budget) ==="
for set in "${SCENESETS[@]}"; do
    IFS='|' read -r tag scene_dir van btc times <<< "$set"
    echo "[$tag]"
    printf '  %-26s' "config"
    for secs in $times; do printf '%10s' "$((secs/60))min"; done; echo
    for cfg in "vanilla" "btcEag${EAG}Split${MAXSPLITS}"; do
        printf '  %-26s' "$cfg"
        for secs in $times; do
            m=$(printf '%02d' $((secs/60)))
            v=$(grep -oE 'Var: [0-9.eE+-]+' "${OUT}/${tag}-${cfg}-t${m}min.log" 2>/dev/null | tail -1 | awk '{print $2}')
            printf '%10s' "${v:-n/a}"
        done; echo
    done
done
echo "=== time comparison DONE $(date) ==="
