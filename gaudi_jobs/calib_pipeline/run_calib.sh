#!/bin/bash
# Run overlay + PFA hit building over every calibration gun, ONE AT A TIME
# (CLAUDE.md: no parallel jobs on lxplus). Everything lives on EOS.
#
#   ./run_calib.sh              # all guns
#   ./run_calib.sh calib_e-_5GeV  # just one
set -e
here=$(cd "$(dirname "$0")" && pwd)
eos=${CALIB_EOS_DIR:-/eos/experiment/ship/user/marquezh/snd_pfa/calib}
mkdir -p "${eos}/reco"

shopt -s nullglob
for f in "${eos}"/output_*.edm4hep.root; do
    label=$(basename "$f" .edm4hep.root); label=${label#output_}
    if [ $# -gt 0 ] && [ "$1" != "$label" ]; then continue; fi
    echo "=== ${label}"
    ( cd "$here" && \
      INPUT_FILE="$f" OUTPUT_FILE="${eos}/reco/events_${label}.root" \
        k4run job1_overlay.py > "${eos}/reco/${label}.job1.log" 2>&1 && \
      INPUT_FILE="${eos}/reco/events_${label}.root" \
      OUTPUT_FILE="${eos}/reco/pfa_${label}.root" \
        k4run job6_pfa.py > "${eos}/reco/${label}.job6.log" 2>&1 )
    grep -h "written=" "${eos}/reco/${label}.job6.log" | tail -1
done
