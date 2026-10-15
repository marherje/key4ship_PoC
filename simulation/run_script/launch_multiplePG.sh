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

# ======================= CONFIGURATION =======================
nevents=${NEVENTS:-100}
physlist="QGSP_BERT"

# Opening angle(s) between the two primaries, in degrees (space-separated).
angles="1.0 2.0 5.0 10.0"

# Default primary vertex [mm]. A pair may override the z (see below); x and y
# are common to all.
pos_x=1
pos_y=1
pos_z=-1000

# Particle pairs: "particle1:E1_GeV:particle2:E2_GeV[:vertex_z_mm]"
#
# The optional fifth field moves the primary vertex along the beam for that
# pair only; omitted, it is pos_z. A pair that overrides it gets a _z<mm>
# suffix in its label, so the default samples keep the names the pipelines
# already reference.
#   mu- (PDG  13) 10 GeV  +  pi+ (PDG 211) 5 GeV
#   e-  (PDG  11) 10 GeV  +  pi+ (PDG 211) 5 GeV
#   mu- (PDG  13) 10 GeV  +  mu+ (PDG  -13) 10 GeV
#
# The two-muon pair feeds gaudi_jobs/2_mu_ang_pipeline. Note that the two
# primaries share a vertex and are separated only by the opening angle: two
# muons that are exactly PARALLEL are a degenerate configuration for the
# SiTarget, whose 1D strips then see a ghost StripX x StripY pairing that is
# itself a perfect straight line, indistinguishable from the real one by any
# fit. Physical events always carry an angle; a test sample must too.
#   mu+/pi+ at z = -2150 fires just inside the SiTarget front face instead of
#   1.2 m into it, so the pion crosses the full ~1.8 m of the detector, i.e.
#   some six nuclear interaction lengths of tungsten. It usually interacts
#   before reaching the SiPad, which is the point: it produces tracks that
#   live and die inside the SiTarget. Feeds gaudi_jobs/mu_pi_early_pipeline.
pairs=(
    "mu-:10:pi+:5"
    "e-:10:pi+:5"
    "mu-:10:mu+:10"
    "mu-:10:pi+:5:-2150"
)
# ==============================================================

mkdir -p steer log

# This script always simulates the committed baseline geometry
# (simulation/geometry/SND_compact.xml). Geometry variants are generated and
# simulated by Optimization/ instead; see Optimization/README.md.

local=$PWD
steer_path="${local}/steer"
# Simulation outputs (HepMC, edm4hep, GENIE gst/ghep) go to EOS: they do not
# fit in the 10 GB AFS quota. ./data is a symlink to the same place, which is
# what the gaudi_jobs pipelines read through. The condor submit files stay on
# AFS, since the schedd refuses /eos paths inside them.
data_path=${KEY4SHIP_DATA_DIR:-/eos/experiment/ship/user/marquezh/key4ship_PoC_data}
mkdir -p "${data_path}"

for pair in "${pairs[@]}"; do
    IFS=: read -r p1 e1 p2 e2 zvtx <<< "${pair}"
    vtx_z=${zvtx:-${pos_z}}

    for angle in ${angles}; do

        label=${physlist}_SND_${p1}_${e1}GeV_${p2}_${e2}GeV_angle_${angle}
        # Only a non-default vertex earns a suffix, so existing sample names
        # (and the pipelines pointing at them) are untouched.
        if [ "${vtx_z}" != "${pos_z}" ]; then
            label=${label}_z${vtx_z}
        fi
        scriptname=${label}.py

        echo "Configuring: ${p1} ${e1} GeV + ${p2} ${e2} GeV, angle=${angle} deg," \
             "vertex z=${vtx_z} mm, ${nevents} events"

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
    position       = (${pos_x}, ${pos_y}, ${vtx_z}),
    outputFile     = "${data_path}/output_${label}.edm4hep.root",
))
EOF

        ./generic_condor_multiplePG.sh ${scriptname} ${label}

    done
done
