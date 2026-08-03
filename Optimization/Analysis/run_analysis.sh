#!/usr/bin/env bash
# Manual re-run of the analysis chain for one variant.
#
#   run_analysis.sh <var_dir> <compact_xml> <sim_edm4hep>
#
# The controller does NOT use this script — it runs the same steps itself so it
# can report which one failed. This is for re-running a variant by hand.
#
# Needs `source init_key4ship.sh` first.
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "usage: $0 <var_dir> <compact_xml> <sim_edm4hep>" >&2
    exit 1
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
var_dir="$(cd "$1" && pwd)"
export SND_COMPACT="$(readlink -f "$2")"
export INPUT_FILE="$(readlink -f "$3")"

# The jobs write bare filenames (events/tracks/ShipHits.root), so they must run
# with the variant directory as cwd.
cd "$var_dir"

k4run "$HERE/job1_overlay.py"  2>&1 | tee job1.log
k4run "$HERE/job4_tracking.py" 2>&1 | tee job4.log
k4run "$HERE/job5_rntuple.py"  2>&1 | tee job5.log
python3 "$HERE/compute_fom.py" 2>&1 | tee fom.log
