#!/bin/bash
#
# Single-particle guns for the SND-PFA energy calibration (phase 1.3 of
# docs/snd_pfa_plan.md).
#
# One job per (particle, energy): a muon to locate the MIP peak of every
# sensitive plane, electrons to fit the electromagnetic scale and its
# resolution, and pions to fit e/h and to measure how much energy leaks into
# the MTC. They all start inside the SiTarget so the shower develops in the
# detector rather than in front of it.
#
#   ./launch_calib_guns.sh              # submit every point to HTCondor
#   ./launch_calib_guns.sh _z200        # only the points whose label matches
#   RUN_LOCAL=1 ./launch_calib_guns.sh  # run in this shell (one ddsim at a
#                                       # time -- heavy, only for a smoke test)
#
# Outputs go to EOS, not AFS: these are ~10 GB of .root that the 10 GB AFS
# quota cannot hold (see docs/snd_pfa_plan.md section B). The condor submit
# files stay on AFS, since the schedd refuses /eos paths inside them.
#
# Adding a point = adding ONE line to the `guns` array.

# ======================= CONFIGURATION =======================
nevents=${NEVENTS:-500}
physlist="QGSP_BERT"
eos_out=${CALIB_EOS_DIR:-/eos/experiment/ship/user/marquezh/snd_pfa/calib}

# Vertex [mm]: 50 mm into the SiTarget (its front face is at z = -1690),
# slightly off-axis so a particle does not ride a module boundary.
pos_x=1
pos_y=1
pos_z=-1640

# "particle:E_GeV[:vertex_z_mm]"
#
# The optional third field moves the vertex for that point only, and adds a
# _zNNN suffix to its label. It exists because the two calorimeters have to be
# calibrated where particles actually start: the SiTarget set below fixes the
# electromagnetic scale, but a GENIE neutrino interacts mostly INSIDE the MTC
# -- it holds ~6300 kg of iron against the SiTarget's ~1300 kg of tungsten --
# and pions fired from the SiTarget deposit at most 2.6 % of their visible
# energy in the MTC (measured in phase 1.4), far too little to fit its scale.
# The z=200 set starts two layers into MTC40 so the shower develops there.
guns=(
    "mu-:10"
    "e-:1"    "e-:2"    "e-:5"    "e-:10"   "e-:20"
    "pi+:2"   "pi+:5"   "pi+:10"  "pi+:20"  "pi+:40"
    "mu-:10:200"
    "pi+:2:200"   "pi+:5:200"   "pi+:10:200"  "pi+:20:200"  "pi+:40:200"
    "e-:5:200"    "e-:20:200"
    # Mean muon deposit per plane vs momentum (2026-10-14): radiative losses
    # make it grow with energy, and SNDPFOBuilder subtracts the muon's mean.
    "mu-:3"   "mu-:30"  "mu-:100"
)
# ==============================================================

set -e
mkdir -p steer data log
mkdir -p "${eos_out}"

local=$PWD
steer_path="${local}/steer"

for gun in "${guns[@]}"; do
    IFS=: read -r part energy zvtx <<< "${gun}"
    vtx_z=${zvtx:-${pos_z}}

    label=calib_${part}_${energy}GeV
    if [ "${vtx_z}" != "${pos_z}" ]; then
        label=${label}_z${vtx_z}
    fi
    scriptname=${label}.py

    # With arguments, submit only the points whose label contains one of them,
    # so a set can be added without re-running (and overwriting) the rest.
    if [ $# -gt 0 ]; then
        match=0
        for pat in "$@"; do
            case "${label}" in *"${pat}"*) match=1;; esac
        done
        [ ${match} -eq 1 ] || continue
    fi

    echo "Configuring: ${part} ${energy} GeV, vertex z=${vtx_z} mm, ${nevents} events"

cat > ${steer_path}/${scriptname} <<EOF
import sys
sys.path.insert(0, "${local}")

from DDSim.DD4hepSimulation import DD4hepSimulation
from multiplePG_base import configure

# One gun: multiplePG_base sends particle 0 along +z, so angleDeg is unused.
SIM = DD4hepSimulation()
configure(SIM, dict(
    physicsList    = "${physlist}",
    numberOfEvents = ${nevents},
    particles      = [("${part}", ${energy})],
    angleDeg       = 0.0,
    position       = (${pos_x}, ${pos_y}, ${vtx_z}),
    outputFile     = "${eos_out}/output_${label}.edm4hep.root",
))
EOF

    ./generic_condor_multiplePG.sh ${scriptname} ${label}
done

echo
echo "Submitted ${#guns[@]} calibration guns; outputs land in ${eos_out}"
