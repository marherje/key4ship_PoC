#!/usr/bin/env bash
# Launch the event display.
#
#   launch.sh <ShipHits.root> <compact.xml> [window]
#
# Both paths are required: there is no default geometry. Pairing a ShipHits.root
# with a compact it was NOT simulated with snaps the hits to the wrong planes
# without any error, so the geometry is never guessed here.
#
#   ./launch.sh ../gaudi_jobs/2_mu_pipeline/ShipHits.root \
#               ../simulation/geometry/SND_compact.xml 1
#
#   ./launch.sh ../Optimization/Analysis/variants/var00007/ShipHits.root \
#               ../Optimization/Simulation/geometry/variants/var00007.xml 3
#
# `window` is optional; when omitted, event_display_eve.py's own default (0) is
# used rather than a second default duplicated here.
#
# Needs `source init_key4ship.sh` first.
set -euo pipefail

usage() {
    sed -n '2,19p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
    exit 2
}

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "error: expected 2 or 3 arguments, got $#" >&2
    echo >&2
    usage
fi

hits="$1"
geometry="$2"

for f in "$hits" "$geometry"; do
    if [ ! -f "$f" ]; then
        echo "error: file not found: $f" >&2
        exit 1
    fi
done

if [ $# -eq 3 ]; then
    python event_display_eve.py --hits "$hits" --geometry "$geometry" --window "$3"
else
    python event_display_eve.py --hits "$hits" --geometry "$geometry"
fi
