#!/bin/bash
# ---------------------------------------------------------------------------
# Cross-product sweep: render the kitchen scene at a fixed spp for every
# (eagerness x maxDepth) combination, each dumping the SD-tree.
#
# Eagerness order: 4 first (current default), then 1, 2, 3.
# For each eagerness, all depths are rendered before moving to the next.
#
# RESUMABLE: a combination whose final <name>.exr already exists is SKIPPED.
# So to resume after stopping, just run this script again with the same config.
# (A combination that was interrupted has no .exr, so it re-renders cleanly.)
#
# Per combination everything lands in its own folder:
#   kitchen-spp511-btcEag<E>Depth<D>-dump/
#     kitchen-spp511-btcEag<E>Depth<D>-dump.exr        <- rendered image
#     kitchen-spp511-btcEag<E>Depth<D>-dump-00.sdt     <- SD-tree dump, iter 0
#     ... -NN.sdt (one per iteration) ... and .log
#
# Run INSIDE the mitsuba container (project bind-mounted at /home/mitsuba):
#   docker exec -d -w /home/mitsuba btc-pg \
#     bash -lc 'bash depthSweepRender.sh > sweep_progress.log 2>&1'
# or interactively:
#   source setpath.sh && nohup bash depthSweepRender.sh > sweep_progress.log 2>&1 &
#   tail -f sweep_progress.log
#
# WARNING: ~9 .sdt dumps per combination, hundreds of MB each -> several GB
# per eagerness, ~tens of GB for the full sweep. Watch free disk.
# ---------------------------------------------------------------------------

cd "$(dirname "$0")"; [ -f setpath.sh ] || cd ..   # run from project root even if invoked from a subdir
source setpath.sh >/dev/null 2>&1
if ! command -v mitsuba >/dev/null 2>&1; then
    echo "ERROR: 'mitsuba' not on PATH. Run from /home/mitsuba after sourcing setpath.sh." >&2
    exit 1
fi

# --- config (edit these) ----------------------------------------------------
BASE="../scenes/kitchen/kitchen-btc-63.xml"   # vanilla: kitchen-vanilla-63.xml
SCENE_DIR="../scenes/kitchen"
SPP=511
TAG="btc"
DEPTHS="16 32 64 124 512"                      # NOTE: 124 looks like a typo for 128 -- edit if so
EAGERNESS_ORDER="4 1 2 3"                       # 4 first, then 1, 2, 3 (per request)
# ---------------------------------------------------------------------------

NE=$(echo $EAGERNESS_ORDER | wc -w)
ND=$(echo $DEPTHS | wc -w)
TOTAL=$(( NE * ND ))
DONE=0          # combinations processed (skipped + rendered)
NR=0            # actually rendered this run (for ETA)
SUM_RT=0        # total render seconds this run (for ETA)
SWEEP_START=$(date +%s)

echo "=== sweep START $(date) | spp=$SPP | $TOTAL combos (eagerness[$EAGERNESS_ORDER] x depth[$DEPTHS]) ==="

for e in $EAGERNESS_ORDER; do
    for d in $DEPTHS; do
        NAME="kitchen-spp${SPP}-${TAG}Eag${e}Depth${d}-dump"   # folder == file base
        FOLDER="$NAME"
        OUTBASE="${FOLDER}/${NAME}"                            # -> NAME.exr , NAME-NN.sdt
        LOG="${FOLDER}/${NAME}.log"
        TMP="${SCENE_DIR}/_sweep_${NAME}.xml"

        mkdir -p "$FOLDER"
        DONE=$((DONE + 1))

        # --- resume: skip if this combination already finished ---
        if [ -s "${OUTBASE}.exr" ]; then
            printf '    SKIP   e=%s d=%-4s already rendered                              [%d/%d]\n' "$e" "$d" "$DONE" "$TOTAL"
            continue
        fi

        # Per-combo scene: set spp budget + sampleCount, inject maxDepth, eagerness, dumpSDTree.
        awk -v spp="$SPP" -v depth="$d" -v eag="$e" '
            /<float name="budget" value=/        { sub(/value="[0-9.]+"/, "value=\"" spp "\"") }
            /<integer name="sampleCount" value=/ { sub(/value="[0-9.]+"/, "value=\"" spp "\"") }
            { print }
            /<integrator type=/ {
                print "\t\t<integer name=\"maxDepth\" value=\"" depth "\"/>"
                print "\t\t<integer name=\"eagerness\" value=\"" eag "\"/>"
                print "\t\t<boolean name=\"dumpSDTree\" value=\"true\"/>"
            }
        ' "$BASE" > "$TMP"

        echo ">>> [$(date +%T)] eagerness=$e maxDepth=$d  ->  ${FOLDER}/"
        RSTART=$(date +%s)
        mitsuba "$TMP" -o "$OUTBASE" > "$LOG" 2>&1 &
        MPID=$!

        # Peak resident memory (can OOM at high depth).
        PEAK=0
        while kill -0 "$MPID" 2>/dev/null; do
            R=$(awk '/VmRSS/{print $2}' /proc/$MPID/status 2>/dev/null)
            [ -n "$R" ] && [ "$R" -gt "$PEAK" ] && PEAK=$R
            sleep 2
        done
        wait "$MPID"; RC=$?
        rm -f "$TMP"

        RT=$(( $(date +%s) - RSTART ))
        NR=$((NR + 1)); SUM_RT=$((SUM_RT + RT))
        AVG=$(( SUM_RT / NR ))
        REMAIN=$(( AVG * (TOTAL - DONE) ))
        AP=$(grep -i "Average path length" "$LOG" | tail -1 | sed 's/.*Average path length *//')
        if [ "$RC" -eq 0 ]; then STATUS="OK     "; else STATUS="FAILED "; fi

        printf '    %s e=%s d=%-4s | %dm%02ds peakRSS=%dMB pathLen=%s rc=%s | [%d/%d] est-remaining=~%dm (~%.1fh)\n' \
            "$STATUS" "$e" "$d" $((RT/60)) $((RT%60)) $((PEAK/1024)) "${AP:-n/a}" "$RC" \
            "$DONE" "$TOTAL" $((REMAIN/60)) "$(awk "BEGIN{print $REMAIN/3600}")"
    done
done

echo "=== sweep DONE $(date) | rendered $NR this run, $DONE/$TOTAL total | $(( ($(date +%s)-SWEEP_START)/60 )) min ==="
