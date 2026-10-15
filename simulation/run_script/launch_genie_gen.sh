#!/bin/bash
#
# GENIE generation + ddsim in ONE HTCondor job: the generator step that
# launch_genie.sh assumes has already been done on EOS.
#
# Per configured sample this script writes steer/rungenie_<label>.sh, which
# runs three steps, each in the environment it needs:
#
#   1. gevgen_fnal + gntpc   (SHiP cvmfs stack: GENIE 3.4.2 built against
#                             Pythia6, tune G21_11a)
#      -> data/genie/<label>/gntp.<run>.ghep.root, gntp.<run>.gst.root
#   2. genie_gst_to_hepmc.py (key4hep stack)
#      -> data/genie_<label>.hepmc
#   3. ddsim                 (key4hep stack, same steering as launch_genie.sh)
#      -> data/output_genie_<label>.edm4hep.root
#
# The two stacks cannot be sourced in the same shell (both carry their own
# ROOT and Python), so steps 1 and 2-3 run in separate subshells.
#
#   ./launch_genie_gen.sh              # submit to HTCondor
#   RUN_LOCAL=1 ./launch_genie_gen.sh  # run in this shell (heavy: gevgen_fnal
#                                      # parses a 470 MB spline file, then a
#                                      # full ddsim -- only for tiny tests)
#
# REFERENCE FRAME. The gsimple flux file has its rays in FairShip coordinates
# and gevgen_fnal cannot translate them (-z only slides each ray along its own
# direction), so the geometry is handed to GENIE in the FairShip frame
# (simulation/geometry/gdml/SND_compact_fairship_frame.gdml, produced by
# shift_gdml_frame.py from the key4hep export) and the vertices are moved back
# by zshift_m in step 2. Regenerate both GDMLs after any geometry change:
#
#   source ../../init_key4ship.sh
#   (cd ../geometry && geoConverter -compact2gdml -input SND_compact.xml \
#                                    -output gdml/SND_compact.gdml)
#   (cd ../geometry/gdml && ./shift_gdml_frame.py SND_compact.gdml \
#                                    SND_compact_fairship_frame.gdml)
#
# Adding a sample = adding ONE line to the `samples` array below.

# ======================= CONFIGURATION =======================
physlist="QGSP_BERT"

# --- GENIE itself lives in run_gevgen.sh (shared with mobo's payload). Its
#     inputs -- flux, splines, tune, SHiP cvmfs release -- are Antonio
#     Iuliano's and can be overridden through GENIE_FLUX / GENIE_XSEC /
#     GENIE_TUNE / GENIE_FLAVOURS / SHIP_SETUP / GENIE_CONFIG_DIR. Condor does
#     not forward the submitter's environment, so whichever of them are set
#     here are written into the job script.
run_gevgen=$(cd "$(dirname "$0")" && pwd)/run_gevgen.sh
genie_env=""
for v in GENIE_FLUX GENIE_XSEC GENIE_TUNE GENIE_FLAVOURS SHIP_SETUP GENIE_CONFIG_DIR; do
    [ -n "${!v:-}" ] && genie_env+="export ${v}=${!v}"$'\n'
done

# --- Geometry GENIE navigates, and the offset that brings its vertices back
#     to the key4ship frame (SiPad front: -37 cm here, 2539 cm in FairShip).
gdml=$(cd "$(dirname "$0")/../geometry/gdml" && pwd)/SND_compact_fairship_frame.gdml
zshift_m=-25.76

# Samples to generate+simulate: "label:nevents:run[:seed]".
# nevents is gevgen_fnal's -n (interactions written; the exposure in POT is
# in the .status file). run numbers the output files and must be unique per
# sample. seed defaults to run.
samples=(
    # Already produced (2026-09-17); listing it again would re-run its ddsim
    # and overwrite the output, so it stays commented out.
    # "gen_test_1000:1000:1"
    # SND-PFA reference sample (docs/snd_pfa_plan.md): four independent runs
    # of 2000 interactions, generated with the FairShip-frame GDML and the
    # zshift_m correction below, i.e. with the flux geometry FIXED. They
    # replace output_genie_allflavours_0_2000.edm4hep.root, whose neutrino
    # cone has its waist at z ~ 0 instead of 25.76 m upstream (measured
    # 2026-09-22: transverse RMS 44 mm at z ~ 0 growing to 105 mm at z = 2.7 m,
    # where a correctly placed source varies by ~15 % across the detector).
    "pfa_ref_a:2000:11"
    "pfa_ref_b:2000:12"
    "pfa_ref_c:2000:13"
    "pfa_ref_d:2000:14"
)
# ==============================================================

mkdir -p steer log

local=$PWD
repo=$(cd "${local}/../.." && pwd)
steer_path="${local}/steer"
# Simulation outputs (HepMC, edm4hep, GENIE gst/ghep) go to EOS: they do not
# fit in the 10 GB AFS quota. ./data is a symlink to the same place, which is
# what the gaudi_jobs pipelines read through. The condor submit files stay on
# AFS, since the schedd refuses /eos paths inside them.
data_path=${KEY4SHIP_DATA_DIR:-/eos/experiment/ship/user/marquezh/key4ship_PoC_data}
mkdir -p "${data_path}"
log_path="${local}/log"

for f in "${gdml}" "${run_gevgen}"; do
    if [ ! -e "${f}" ]; then
        echo "ERROR: missing input: ${f}"
        exit 1
    fi
done

for sample in "${samples[@]}"; do
    IFS=: read -r label nevents run seed <<< "${sample}"
    seed=${seed:-${run}}

    genie_dir="${data_path}/genie/${label}"
    gst_file="${genie_dir}/gntp.${run}.gst.root"
    hepmc_file="${data_path}/genie_${label}.hepmc"
    scriptname=genie_${label}.py
    condorsh=rungenie_${label}.sh
    condorsub=rungenie_${label}.sub
    mkdir -p "${genie_dir}"

    # ----------------------------------------------------------
    # DDSIM STEERING FILE — identical in shape to launch_genie.sh's
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

    # ----------------------------------------------------------
    # JOB SCRIPT (also usable standalone: bash steer/rungenie_<label>.sh)
    # ----------------------------------------------------------
cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
set -e
echo "Starting job on \$(hostname)"

# ---------- 1. GENIE (SHiP stack, in its own subshell) ----------
if [ -s "${gst_file}" ]; then
    echo "GENIE output ${gst_file} already exists, skipping generation"
else
    (
${genie_env}    bash ${run_gevgen} ${gdml} ${genie_dir} ${nevents} ${run} ${seed} ${genie_dir}
    )
    cp -f ${genie_dir}/gevgen_${run}.log ${log_path}/genie_${label}.gevgen.log
fi

# ---------- 2+3. HepMC conversion and ddsim (key4hep stack) ----------
(
    set -e
    if [ -z "\${KEY4HEP_STACK:-}" ]; then
        source ${repo}/init_key4hep.sh
    fi
    export LD_LIBRARY_PATH=${repo}/install/lib64:${repo}/install/lib:\$LD_LIBRARY_PATH
    export PYTHONPATH=${repo}/install/lib64:${repo}/install/lib:${repo}/install/python:\$PYTHONPATH

    python3 ${local}/genie_gst_to_hepmc.py \\
        --input ${gst_file} --start 0 --nevents -1 \\
        --z-shift-m ${zshift_m} \\
        --output ${hepmc_file}

    # ddsim occasionally aborts during ROOT/cppyy teardown AFTER the run has
    # finished and the output file is fully written. Tolerate exactly that
    # case (event loop completed) instead of failing the whole job.
    rc=0
    ddsim --steeringFile ${steer_path}/${scriptname} &> ${log_path}/genie_${label}.ddsim.log || rc=\$?
    if [ \$rc -ne 0 ]; then
        if grep -q "Finished run 0 after" ${log_path}/genie_${label}.ddsim.log; then
            echo "WARNING: ddsim exited with code \$rc after completing the run (teardown crash); output kept."
        else
            echo "ERROR: ddsim failed (code \$rc), see ${log_path}/genie_${label}.ddsim.log"
            exit \$rc
        fi
    fi
)
echo "Job finished"
EOF
    chmod +x ${steer_path}/${condorsh}

    if [ "${RUN_LOCAL:-0}" = "1" ]; then
        echo "RUN_LOCAL=1: running ${label} in this shell"
        bash ${steer_path}/${condorsh}
        continue
    fi

    # ----------------------------------------------------------
    # CONDOR SUBMIT. request_memory: gevgen_fnal parses the 470 MB spline
    # XML into memory before generating anything. SendCredential: the flux
    # and splines are read from EOS.
    # ----------------------------------------------------------
cat > ${steer_path}/${condorsub} <<EOF
executable              = ${steer_path}/${condorsh}
log                     = ${log_path}/genie_${label}.log
output                  = ${log_path}/outfile_genie_${label}.txt
error                   = ${log_path}/errors_genie_${label}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
request_memory          = 6000
MY.SendCredential       = True
+JobFlavour             = "tomorrow"
queue 1
EOF

    (cd ${steer_path} && condor_submit ${condorsub})
done
