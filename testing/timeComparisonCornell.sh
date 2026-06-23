#!/bin/bash
# ---------------------------------------------------------------------------
# Equal-time comparison on the CORNELL BOX scene:
#   vanilla (Mueller PPG) vs BTC (eagerness 4, maxDepth 512)  at  5 and 20 min.
#
# 4 renders = 2 configs x 2 time budgets. Run SEQUENTIALLY (each gets the whole
# machine for its slice -- that's what makes equal-time fair). Dumps OFF; both
# configs use the same maxDepth=512 so only the directional structure differs.
#
# Uses the 0.6-format Cornell scene (scenes/cornell-box/scene_v0.6.xml) turned
# into guided variants: scenes/cornell-box/cornell-{btc,vanilla}.xml
#
# RESUMABLE: a render whose .exr already exists is skipped.
#
# Outputs in timeComparisonCornell_out/:
#   cornell-vanilla-t05min.{exr,log}
#   cornell-btcEag4Depth512-t05min.{exr,log}
#   ... t20min ...
#
# Run INSIDE the container (project at /home/mitsuba):
#   docker exec -d -w /home/mitsuba btc-pg bash -lc 'bash timeComparisonCornell.sh > cornell_timetest.log 2>&1'
# or:  source setpath.sh && nohup bash timeComparisonCornell.sh > cornell_timetest.log 2>&1 &
#      tail -f cornell_timetest.log
#
# Total wall-clock ~= 2*(5+20) = ~50 min (+ a little overhead).
# ---------------------------------------------------------------------------

cd "$(dirname "$0")"; [ -f setpath.sh ] || cd ..   # run from project root even if invoked from a subdir
source setpath.sh >/dev/null 2>&1
if ! command -v mitsuba >/dev/null 2>&1; then
    echo "ERROR: 'mitsuba' not on PATH. Run from /home/mitsuba after sourcing setpath.sh." >&2
    exit 1
fi

# --- config -----------------------------------------------------------------
SCENE_DIR="../scenes/cornell-box"
VANILLA_SCENE="${SCENE_DIR}/cornell-vanilla.xml"
BTC_SCENE="${SCENE_DIR}/cornell-btc.xml"
OUT="timeComparisonCornell_out"
DEPTH=512                      # same path-length cap for BOTH (fairness)
TIMES_SEC="300 1200"           # 5, 20 minutes
# ---------------------------------------------------------------------------

mkdir -p "$OUT"

# run_one <label> <base_scene> <budget_seconds> <extra_integrator_props>
run_one() {
    local label="$1" base="$2" secs="$3" extra="$4"
    local mins=$(( secs / 60 ))
    local name; name="cornell-${label}-t$(printf '%02d' "$mins")min"
    local out="${OUT}/${name}"
    local log="${OUT}/${name}.log"
    local tmp="${SCENE_DIR}/_timetest_${name}.xml"   # MUST live in the scene dir (relative geometry)

    if [ -s "${out}.exr" ]; then
        printf '    SKIP  %-36s (already rendered)\n' "$name"
        return
    fi

    awk -v secs="$secs" -v depth="$DEPTH" -v extra="$extra" '
        /<string name="budgetType"/   { print "\t\t<string name=\"budgetType\" value=\"seconds\"/>"; next }
        /<float name="budget" value=/ { print "\t\t<float name=\"budget\" value=\"" secs "\"/>"; next }
        { print }
        /<integrator type=/ {
            print "\t\t<integer name=\"maxDepth\" value=\"" depth "\"/>"
            if (extra != "") print extra
        }
    ' "$base" > "$tmp"

    echo ">>> [$(date +%T)] $name  (budget ${secs}s = ${mins}m)"
    local RSTART; RSTART=$(date +%s)
    mitsuba "$tmp" -o "$out" > "$log" 2>&1 &
    local MPID=$!
    local PEAK=0 R
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
    printf '    done  %-36s wall=%dm%02ds peakRSS=%dMB rc=%s finalVar=%s pathLen=%s\n' \
        "$name" $((RT/60)) $((RT%60)) $((PEAK/1024)) "$RC" "${VAR:-n/a}" "${AP:-n/a}"
}

echo "=== Cornell time comparison START $(date) | budgets: $TIMES_SEC s | maxDepth=$DEPTH ==="

for secs in $TIMES_SEC; do
    echo "--- ${secs}s ($((secs/60)) min) ---"
    run_one "vanilla"         "$VANILLA_SCENE" "$secs" ""
    run_one "btcEag4Depth512" "$BTC_SCENE"     "$secs" "\t\t<integer name=\"eagerness\" value=\"4\"/>"
done

echo
echo "=== SUMMARY (lower finalVar = less noise = better at that budget) ==="
printf '%-24s %12s %12s\n' "config" "5min Var" "20min Var"
for label in vanilla btcEag4Depth512; do
    line=$(printf '%-24s' "$label")
    for mins in 05 20; do
        v=$(grep -oE 'Var: [0-9.eE+-]+' "${OUT}/cornell-${label}-t${mins}min.log" 2>/dev/null | tail -1 | awk '{print $2}')
        line+=$(printf ' %12s' "${v:-n/a}")
    done
    echo "$line"
done
echo "=== Cornell time comparison DONE $(date) ==="
