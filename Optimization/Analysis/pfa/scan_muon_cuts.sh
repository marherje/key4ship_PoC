#!/bin/bash
# Phase 3.3.6: efficiency against fakes over the muon identification cuts.
#
# Runs job6 once per point, each writing its own file on EOS so the scan can be
# re-analysed without re-running, then prints the table. One process at a time
# (CLAUDE.md).
set -e
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "${here}/../../.." && pwd)
out=${SCAN_DIR:-/eos/experiment/ship/user/marquezh/snd_pfa/scratch/muon_scan}
nev=${NEVENTS:-500}
mkdir -p "$out"

for layers in 6 8 10; do
  for width in 10 15 25; do
    tag="L${layers}_W${width}"
    f="${out}/pfa_${tag}.root"
    if [ ! -s "$f" ]; then
      echo "=== ${tag}"
      ( cd "${repo}/gaudi_jobs/genie_pipeline" && \
        OUTPUT_FILE="$f" k4run job6_pfa.py -n "$nev" \
          --SNDMuonFinder.MinMTCLayers=${layers} \
          --SNDMuonFinder.CorridorHalfWidthMTCSciFi=${width} \
          > "${out}/${tag}.log" 2>&1 )
    fi
    grep -h "muon candidates" "${out}/${tag}.log" | tail -1
  done
done

echo
echo "=== efficiency / fakes"
for layers in 6 8 10; do
  for width in 10 15 25; do
    tag="L${layers}_W${width}"
    echo "--- MinMTCLayers=${layers} CorridorHalfWidthMTCSciFi=${width}"
    python3 "${here}/muon_perf.py" "${out}/pfa_${tag}.root" --max-events="$nev" \
      2>/dev/null | grep -E "muon had a GlobalTrack|identified among|^  NC|^  nue CC|all  "
  done
done
