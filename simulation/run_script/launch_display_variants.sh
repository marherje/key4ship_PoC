#!/bin/bash
#
# Re-run a finished MOBO trial to get its ROOT files back, for event displays.
#
# A trial deletes its intermediates when it finishes (keep_intermediate: false),
# so only the geometry, the params and the logs survive on EOS. This script
# takes a prepared workdir -- a copy of the trial directory whose trial.json has
# keep_intermediate set and, usually, a smaller nevents -- and submits the same
# payload that produced the campaign, so the events come out of the identical
# geometry rather than an approximation of it.
#
#   ./launch_display_variants.sh <workdir> [<workdir> ...]
#
# The submit file stays on AFS because the schedd refuses /eos paths inside
# one; the job itself works entirely in the EOS workdir.

set -e
if [ $# -lt 1 ]; then
    echo "Usage: $0 <workdir> [<workdir> ...]"
    exit 1
fi

local=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "${local}/../.." && pwd)
steer="${local}/steer"
logs="${local}/log"
mkdir -p "${steer}" "${logs}"

for workdir in "$@"; do
    if [ ! -f "${workdir}/trial.json" ]; then
        echo "ERROR: no trial.json in ${workdir}"
        exit 1
    fi
    label=$(basename "${workdir}")

cat > "${steer}/rundisplay_${label}.sh" <<EOF
#!/bin/bash
echo "display re-run ${label} on \$(hostname) at \$(date)"
if [ -z "\${KEY4HEP_STACK:-}" ]; then
    source ${repo}/init_key4ship.sh
fi
cd ${workdir} || exit 1
python3 ${repo}/Optimization/mobo/src/mobo/ship/payload.py --workdir ${workdir}
rc=\$?
echo "finished at \$(date) with code \$rc"
exit \$rc
EOF
    chmod +x "${steer}/rundisplay_${label}.sh"

cat > "${steer}/rundisplay_${label}.sub" <<EOF
executable              = ${steer}/rundisplay_${label}.sh
log                     = ${logs}/display_${label}.log
output                  = ${logs}/display_${label}.out
error                   = ${logs}/display_${label}.err
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
request_memory          = 6000
MY.SendCredential       = True
+JobFlavour             = "tomorrow"
queue 1
EOF

    (cd "${steer}" && condor_submit "rundisplay_${label}.sub")
done
