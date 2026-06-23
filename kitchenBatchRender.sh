#!/bin/bash
set -e

echo "=== Checking all XML files exist ==="

missing=0
for spp in 63 127 255 511 1023 2047; do
    f="../scenes/kitchen/kitchen-btc-${spp}.xml"
    if [ ! -f "$f" ]; then
        echo "MISSING: $f"
        missing=1
    fi
done

for spp in 63 127 255 511 1023 2047; do
    f="../scenes/kitchen/kitchen-vanilla-${spp}.xml"
    if [ ! -f "$f" ]; then
        echo "MISSING: $f"
        missing=1
    fi
done

if [ "$missing" -eq 1 ]; then
    echo "Aborting: one or more files missing."
    exit 1
fi

echo "All files found."
echo ""
echo "=== Kitchen BTC Tests ==="

for spp in 63 127 255 511 1023 2047; do
    echo "Running kitchen-btc-${spp}..."
    mitsuba ../scenes/kitchen/kitchen-btc-${spp}.xml
    echo "Done: kitchen-btc-${spp}"
done

echo "=== Kitchen Vanilla Tests ==="

for spp in 63 127 255 511 1023 2047; do
    echo "Running kitchen-vanilla-${spp}..."
    mitsuba ../scenes/kitchen/kitchen-vanilla-${spp}.xml
    echo "Done: kitchen-vanilla-${spp}"
done

echo "=== All tests complete ==="
