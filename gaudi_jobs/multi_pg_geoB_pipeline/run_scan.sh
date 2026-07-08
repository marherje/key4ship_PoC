#!/bin/bash
#
# Reconstruction pipeline for the geoB geometry variant.
# Runs overlay + digitize + tracking over every geoB multi-PG sample
# (mu+pi and e+pi, angles 0.5/1/5/10) and keeps one tracks file per sample:
#   tracks_<pair>_angle_<a>.edm4hep.root
#
# Requires the key4hep env + install/ paths (see build.sh).

set -e

data="../../simulation/run_script/data"
angles="0.5 1.0 5.0 10.0"

for pair in mu+:mu_pi e+:e_pi; do
  IFS=: read -r p1 tag <<< "$pair"
  for a in $angles; do
    f="$data/output_geoB_QGSP_BERT_SND_${p1}_10GeV_pi+_5GeV_angle_${a}.edm4hep.root"
    if [ ! -f "$f" ]; then echo "SKIP (no existe): $f"; continue; fi
    echo "=== geoB ${tag} angle=${a}"
    INPUT_FILE="$f" k4run job1_overlay.py
    INPUT_FILE=events.edm4hep.root k4run job3_digitize.py
    k4run job4_tracking.py
    cp tracks.edm4hep.root tracks_${tag}_angle_${a}.edm4hep.root
    TRACKS_FILE=tracks_${tag}_angle_${a}.edm4hep.root \
      SHIPHITS_FILE=ShipHits_${tag}_angle_${a}.root k4run job5_rntuple.py
  done
done
echo "Done:"
ls -1 tracks_*_angle_*.edm4hep.root
