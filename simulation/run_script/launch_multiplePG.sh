#!/bin/bash
#
# Launch a set of multiple-Particle-Gun simulations.
#
# For every configured particle pair and opening angle this script:
#   1. generates a steering file in steer/ (all the physics logic lives in
#      multiplePG_base.py — the steering file is only the configuration),
#   2. submits it with generic_condor_multiplePG.sh.
#
# Adding a new particle pair = adding ONE line to the `pairs` array below.
#
#   ./launch_multiplePG.sh              # submit everything to HTCondor
#   RUN_LOCAL=1 ./launch_multiplePG.sh  # run everything locally instead

set -e

# ======================= CONFIGURATION =======================
nevents=100
physlist="QGSP_BERT"

# Opening angle(s) between the two primaries, in degrees (space-separated).
angles="1.0 2.0 5.0 10.0"

# Common primary vertex [mm]
pos_x=1
pos_y=1
pos_z=-1000

# Particle pairs: "particle1:E1_GeV:particle2:E2_GeV"
#   mu+ (PDG -13) 10 GeV  +  pi+ (PDG 211) 5 GeV
#   e+  (PDG -11) 10 GeV  +  pi+ (PDG 211) 5 GeV
pairs=(
    "mu+:10:pi+:5"
    "e+:10:pi+:5"
)
# ==============================================================

mkdir -p steer data log

# Apply the geometry parametrization (parameters.yaml -> SND_compact.xml)
# before simulating, so the samples always match the declared geometry.
python3 ../geometry/config.py || { echo "ERROR: geometry config failed"; exit 1; }

local=$PWD
steer_path="${local}/steer"
data_path="${local}/data"

for pair in "${pairs[@]}"; do
    IFS=: read -r p1 e1 p2 e2 <<< "${pair}"

    for angle in ${angles}; do

        label=${physlist}_SND_${p1}_${e1}GeV_${p2}_${e2}GeV_angle_${angle}
        scriptname=${label}.py

        echo "Configuring: ${p1} ${e1} GeV + ${p2} ${e2} GeV, angle=${angle} deg, ${nevents} events"

        # ----------------------------------------------------------
        # STEERING FILE — thin config over multiplePG_base.configure()
        # ----------------------------------------------------------
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
    outputFile     = "${data_path}/output_${label}.edm4hep.root",
))
EOF

        ./generic_condor_multiplePG.sh ${scriptname} ${label}

    done
done
