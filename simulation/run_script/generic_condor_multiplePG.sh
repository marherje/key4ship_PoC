#!/bin/bash
#
# Submit ONE multiple-Particle-Gun ddsim job to HTCondor.
#
# Usage:
#   ./generic_condor_multiplePG.sh <steering.py> <label>
#
#   <steering.py>  steering file name inside steer/ (created by
#                  launch_multiplePG.sh, or hand-written from the template in
#                  multiplePG_base.py's docstring)
#   <label>        unique job label, used for all log file names
#
# Example:
#   ./generic_condor_multiplePG.sh \
#       QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0.py \
#       QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0
#
# Set RUN_LOCAL=1 to run ddsim in the current shell instead of submitting
# to condor (useful for quick tests):
#   RUN_LOCAL=1 ./generic_condor_multiplePG.sh <steering.py> <label>
#
# Logs (separate per configuration):
#   log/<label>.ddsim.log       ddsim output
#   log/outfile_<label>.txt     condor stdout
#   log/errors_<label>.txt      condor stderr
#   log/<label>.log             condor job log

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <steering.py> <label>"
    exit 1
fi

scriptname=$1
label=$2

# Crear directorios si no existen
mkdir -p steer data log

local=$PWD
steer_path="${local}/steer"
log_path="${local}/log"

if [ ! -f "${steer_path}/${scriptname}" ]; then
    echo "ERROR: steering file not found: ${steer_path}/${scriptname}"
    echo "Generate it with launch_multiplePG.sh or place it in steer/."
    exit 1
fi

condorsh=runddsim_${label}.sh
condorsub=runddsim_${label}.sub

# --------------------------------------------------
# JOB SCRIPT (also usable standalone: bash steer/runddsim_<label>.sh)
# --------------------------------------------------
cat > ${steer_path}/${condorsh} <<EOF
#!/bin/bash
set -e

echo "Starting job on \$(hostname)"

# Only initialize the stack if this shell does not have it yet: re-sourcing
# key4hep's setup in an initialized shell returns an error and would abort
# the job (relevant when running with RUN_LOCAL=1 in an existing session).
if [ -z "\${KEY4HEP_STACK:-}" ]; then
    source ${local}/../../init_key4hep.sh
fi
export LD_LIBRARY_PATH=${local}/../../install/lib64:${local}/../../install/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=${local}/../../install/lib64:${local}/../../install/lib:${local}/../../install/python:\$PYTHONPATH

# ddsim occasionally aborts during ROOT/cppyy teardown AFTER the run has
# finished and the output file is fully written. Tolerate exactly that case
# (event loop completed) instead of failing the whole production.
rc=0
ddsim --steeringFile ${steer_path}/${scriptname} &> ${log_path}/${label}.ddsim.log || rc=\$?
if [ \$rc -ne 0 ]; then
    if grep -q "Finished run 0 after" ${log_path}/${label}.ddsim.log; then
        echo "WARNING: ddsim exited with code \$rc after completing the run (teardown crash); output kept."
    else
        echo "ERROR: ddsim failed (code \$rc), see ${log_path}/${label}.ddsim.log"
        exit \$rc
    fi
fi

echo "Job finished"
EOF
chmod +x ${steer_path}/${condorsh}

if [ "${RUN_LOCAL:-0}" = "1" ]; then
    echo "RUN_LOCAL=1: running ddsim locally for ${label}"
    bash ${steer_path}/${condorsh}
    exit $?
fi

# --------------------------------------------------
# CONDOR SUBMIT
# --------------------------------------------------
cat > ${steer_path}/${condorsub} <<EOF
executable              = ${steer_path}/${condorsh}
log                     = ${log_path}/${label}.log
output                  = ${log_path}/outfile_${label}.txt
error                   = ${log_path}/errors_${label}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
+JobFlavour             = "tomorrow"
queue 1
EOF

cd ${steer_path}
condor_submit ${condorsub}
cd -
