#!/usr/bin/env bash
# Launch the event display.
#
#   launch.sh [ShipHits.root] [compact.xml] [window]
#
# Defaults reproduce the previous behaviour: the 1_mu_pipeline hits on the
# committed baseline geometry, window 1. Pass a variant's ShipHits.root and its
# compact to display a geometry variant, e.g.
#
#   ./launch.sh ../Optimization/Analysis/variants/var00007/ShipHits.root \
#               ../Optimization/Simulation/geometry/variants/var00007.xml
#
# Needs `source init_key4ship.sh` first.
set -euo pipefail

hits="${1:-../gaudi_jobs/1_mu_pipeline/ShipHits.root}"
geometry="${2:-../simulation/geometry/SND_compact.xml}"
window="${3:-1}"

python event_display_eve.py --hits "$hits" --geometry "$geometry" --window "$window"
