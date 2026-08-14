#!/bin/bash
#
# Launch a ddsim run over GENIE-generated neutrino interactions.
#
# For every configured chunk of the GENIE "gst" summary tree this script:
#   1. converts that chunk to a HepMC file (genie_gst_to_hepmc.py) unless it
#      already exists,
#   2. generates a steering file in steer/ (all the physics logic lives in
#      genie_base.py -- the steering file is only the configuration, exactly
#      like multiplePG_base.py / launch_multiplePG.sh),
#   3. submits it with generic_condor_multiplePG.sh (generator-agnostic: it
#      just runs `ddsim --steeringFile <script>`).
#
#   ./launch_genie.sh              # submit everything to HTCondor
#   RUN_LOCAL=1 ./launch_genie.sh  # run everything locally instead
#
# Adding a new chunk = adding ONE line to the `chunks` array below.

# ======================= CONFIGURATION =======================
physlist="QGSP_BERT"

# GENIE gst.root summary tree to read from. Override with:
#   GENIE_INPUT=/path/to/gntp.0.gst.root ./launch_genie.sh
genie_input=${GENIE_INPUT:-/eos/experiment/ship/user/aiuliano/GENIE_sims/2026_07_30_key4ship_allflavours/gntp.0.gst.root}

# Chunks of the gst tree to convert+simulate: "start:nevents[:label]".
# nevents=-1 means "to the end of the tree". The label defaults to
# "start<start>_n<nevents>" when omitted.
#
# Default: one modest chunk (all flavours, CC+NC mixed, as produced by
# gevgen_fnal) so a first run stays fast; add more chunks/production-size
# ranges as needed.
chunks=(
    "0:2000:allflavours_0_2000"
)
# ==============================================================

mkdir -p steer data log

# This script always simulates the committed baseline geometry
# (simulation/geometry/SND_compact.xml), same as launch_multiplePG.sh.

local=$PWD
steer_path="${local}/steer"
data_path="${local}/data"

for chunk in "${chunks[@]}"; do
    IFS=: read -r start nevents label <<< "${chunk}"
    label=${label:-start${start}_n${nevents}}

    hepmc_file="${data_path}/genie_${label}.hepmc"
    scriptname=genie_${label}.py

    if [ -f "${hepmc_file}" ]; then
        echo "Skipping conversion for ${label}: ${hepmc_file} already exists"
    else
        echo "Converting GENIE entries [${start}, ${start}+${nevents}) -> ${hepmc_file}"
        python3 "${local}/genie_gst_to_hepmc.py" \
            --input "${genie_input}" --start "${start}" --nevents "${nevents}" \
            --output "${hepmc_file}"
    fi

    # ----------------------------------------------------------
    # STEERING FILE — thin config over genie_base.configure()
    # ----------------------------------------------------------
cat > ${steer_path}/${scriptname} <<EOF
import sys
sys.path.insert(0, "${local}")

from DDSim.DD4hepSimulation import DD4hepSimulation
from genie_base import configure

SIM = DD4hepSimulation()
configure(SIM, dict(
    physicsList = "${physlist}",
    hepmcFile   = "${hepmc_file}",
    outputFile  = "${data_path}/output_genie_${label}.edm4hep.root",
))
EOF

    ./generic_condor_multiplePG.sh ${scriptname} genie_${label}

done
