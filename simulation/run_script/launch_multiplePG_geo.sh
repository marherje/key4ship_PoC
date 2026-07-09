#!/bin/bash
#
# Multi-PG simulation pipeline for a GEOMETRY VARIANT (coexists with the
# default geometry of SND_compact.xml).
#
#   ./launch_multiplePG_geo.sh A     # geoA: SiTarget 30 + SiPad 20 + MTC 20x3
#   ./launch_multiplePG_geo.sh B     # geoB: SiTarget 20 + SiPad 20 + MTC 20x3
#   RUN_LOCAL=1 ./launch_multiplePG_geo.sh A   # run locally instead of condor
#
# The variant XML is regenerated from simulation/geometry/parameters_geo<X>.yaml
# via config.py --output (validated, SND_compact.xml untouched). Output files
# are tagged with the geometry: data/output_geo<X>_<label>.edm4hep.root


if [ $# -ne 1 ]; then
    echo "Usage: $0 <A|B>"
    exit 1
fi
GEO=$1

# ======================= CONFIGURATION =======================
nevents=100
physlist="QGSP_BERT"
angles="0.5 1.0 5.0 10.0"
pos_x=1
pos_y=1
pos_z=-1000
pairs=(
    "mu+:10:pi+:5"
    "e+:10:pi+:5"
)
# ==============================================================

mkdir -p steer data log

local=$PWD
steer_path="${local}/steer"
data_path="${local}/data"
geo_dir="${local}/../geometry"
geo_xml="${geo_dir}/SND_compact_geo${GEO}.xml"

# Regenerate the variant geometry from its parameters (validated).
python3 ${geo_dir}/config.py \
    --params ${geo_dir}/parameters_geo${GEO}.yaml \
    --output ${geo_xml} || { echo "ERROR: geometry config failed"; exit 1; }

for pair in "${pairs[@]}"; do
    IFS=: read -r p1 e1 p2 e2 <<< "${pair}"

    for angle in ${angles}; do

        label=geo${GEO}_${physlist}_SND_${p1}_${e1}GeV_${p2}_${e2}GeV_angle_${angle}
        scriptname=${label}.py

        echo "Configuring [geo${GEO}]: ${p1} ${e1} GeV + ${p2} ${e2} GeV, angle=${angle} deg"

cat > ${steer_path}/${scriptname} <<EOF
import sys
sys.path.insert(0, "${local}")

from DDSim.DD4hepSimulation import DD4hepSimulation
from multiplePG_base import configure

SIM = DD4hepSimulation()
configure(SIM, dict(
    physicsList    = "${physlist}",
    numberOfEvents = ${nevents},
    particles      = [("${p1}", ${e1}), ("${p2}", ${e2})],
    angleDeg       = ${angle},
    position       = (${pos_x}, ${pos_y}, ${pos_z}),
    compactFile    = "${geo_xml}",
    outputFile     = "${data_path}/output_${label}.edm4hep.root",
))
EOF

        ./generic_condor_multiplePG.sh ${scriptname} ${label}

    done
done
