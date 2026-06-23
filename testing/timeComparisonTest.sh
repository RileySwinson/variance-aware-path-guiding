#!/bin/bash
# ---------------------------------------------------------------------------
# Equal-time comparison: vanilla (Mueller PPG) vs BTC (eagerness 4, maxDepth 512)
# at 5, 20 and 60 minute wall-clock budgets (budgetType="seconds").
#
# 6 renders total = 2 configs x 3 time budgets. Run SEQUENTIALLY (each gets the
# whole machine for its time slice -- that's what makes it a fair equal-time test).
# Dumps are OFF (writing .sdt would eat into the time budget). Both configs use the
# same maxDepth=512 so only the directional structure differs.
#
# RESUMABLE: a render whose .exr already exists is skipped.
#
# Outputs in timeComparison_out/:
#   kitchen-vanilla-t05min.exr / .log
#   kitchen-btcEag4Depth512-t05min.exr / .log
#   ... t20min, t60min ...
#
# Run INSIDE the container (project at /home/mitsuba):
#   docker exec -d -w /home/mitsuba btc-pg bash -lc 'bash timeComparisonTest.sh > timetest_progress.log 2>&1'
# or:  source setpath.sh && nohup bash timeComparisonTest.sh > timetest_progress.log 2>&1 &
#      tail -f timetest_progress.log
#
# Total wall-clock ~= 2*(5+20+60) = ~170 min (~3 h) plus a little per-render overhead.
# ---------------------------------------------------------------------------

cd "$(dirname "$0")"; [ -f setpath.sh ] || cd ..   # run from project root even if invoked from a subdir
source setpath.sh >/dev/null 2>&1
if ! command -v mitsuba >/dev/null 2>&1; then
    echo "ERROR: 'mitsuba' not on PATH. Run from /home/mitsuba after sourcing setpath.sh." >&2
    exit 1
fi

# --- config -----------------------------------------------------------------
SCENE_DIR="../scenes/kitchen"
VANILLA_SCENE="${SCENE_DIR}/kitchen-vanilla-63.xml"
BTC_SCENE="${SCENE_DIR}/kitchen-btc-63.xml"
OUT="timeComparison_out"
DEPTH=512                      # same path-length cap for BOTH (fairness)
TIMES_SEC="300 1200 3600"      # 5, 20, 60 minutes
# ---------------------------------------------------------------------------

mkdir -p "$OUT"

# run_one <label> <base_scene> <budget_seconds> <extra_integrator_props>
run_one() {
    local label="$1" base="$2" secs="$3" extra="$4"
    local mins=$(( secs / 60 ))
    local name; name="kitchen-${label}-t$(printf '%02d' "$mins")min"
    local out="${OUT}/${name}"
    local log="${OUT}/${name}.log"
    local tmp="${SCENE_DIR}/_timetest_${name}.xml"

    if [ -s "${out}.exr" ]; then
        printf '    SKIP  %-34s (already rendered)\n' "$name"
        return
    fi

    # Force budgetType=seconds + budget=secs, inject maxDepth (+extra), dumps stay off.
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
    printf '    done  %-34s wall=%dm%02ds peakRSS=%dMB rc=%s finalVar=%s pathLen=%s\n' \
        "$name" $((RT/60)) $((RT%60)) $((PEAK/1024)) "$RC" "${VAR:-n/a}" "${AP:-n/a}"
}

echo "=== time comparison START $(date) | budgets: $TIMES_SEC s | maxDepth=$DEPTH ==="

# Grouped by time budget so each equal-time pair finishes together.
for secs in $TIMES_SEC; do
    echo "--- ${secs}s ($((secs/60)) min) ---"
    run_one "vanilla"         "$VANILLA_SCENE" "$secs" ""
    run_one "btcEag4Depth512" "$BTC_SCENE"     "$secs" "\t\t<integer name=\"eagerness\" value=\"4\"/>"
done

echo
echo "=== SUMMARY (lower finalVar = less noise = better at that budget) ==="
printf '%-22s %12s %12s %12s\n' "config" "5min Var" "20min Var" "60min Var"
for label in vanilla btcEag4Depth512; do
    line=$(printf '%-22s' "$label")
    for mins in 05 20 60; do
        v=$(grep -oE 'Var: [0-9.eE+-]+' "${OUT}/kitchen-${label}-t${mins}min.log" 2>/dev/null | tail -1 | awk '{print $2}')
        line+=$(printf ' %12s' "${v:-n/a}")
    done
    echo "$line"
done
echo "=== time comparison DONE $(date) ==="
