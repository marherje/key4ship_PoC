#!/bin/bash
#
# Generate GENIE neutrino interactions in a GDML geometry: gevgen_fnal + gntpc.
#
#   run_gevgen.sh <gdml> <outdir> <nevents> <run> <seed> <logdir>
#
# Writes <outdir>/gntp.<run>.ghep.root and gntp.<run>.gst.root, logs in
# <logdir>/gevgen_<run>.log and gntpc_<run>.log. Exit code is gevgen_fnal's
# or gntpc's. The POT normalization of the sample is the weight of the ghep
# file's gtree (gevgen_fnal stores it there); the .status file next to the
# output has the last event record and the timing.
#
# Runs on the SHiP cvmfs stack (GENIE 3.4.2, tune G21_11a, Pythia6), which
# cannot share a shell with key4hep: call this from a CLEAN environment
# (env -i, keeping HOME/USER/PATH and the kerberos cache for EOS), never from a
# shell that has sourced init_key4ship.sh. launch_genie_gen.sh does it in a
# subshell before sourcing key4hep; mobo's payload uses env -i.
#
# The GDML must be in the FairShip frame (shift_gdml_frame.py): the gsimple
# flux rays are, and gevgen_fnal cannot translate them.
#
# Inputs (Antonio Iuliano's gevgen_fnal setup) can be overridden through the
# environment: GENIE_FLUX, GENIE_XSEC, GENIE_TUNE, GENIE_FLAVOURS, SHIP_SETUP,
# GENIE_CONFIG_DIR (GXMLPATH; leave unset on this stack -- FairShip's
# genie_config selects Pythia8Decayer2023, which the cvmfs GENIE lacks).

set -e

if [ $# -ne 6 ]; then
    echo "Usage: $0 <gdml> <outdir> <nevents> <run> <seed> <logdir>" >&2
    exit 2
fi
gdml=$(readlink -f "$1"); outdir=$(readlink -f "$2"); nevents=$3; run=$4; seed=$5; logdir=$(readlink -f "$6")

genie_inputs=/eos/experiment/ship/user/aiuliano/inputs_gevgen_fnal
tune=${GENIE_TUNE:-G21_11a_00_000}
# The Shuffled file: the merged one is flavour-block-ordered and a short run
# inherits the composition of whichever block it starts in.
flux=${GENIE_FLUX:-${genie_inputs}/flux_gsimple_Hanae_makeCascade_merged_2026_06/Shuffled_Mbias2026noCharm_makeCascade_merged_nuflux_gsimple.root}
flavours=${GENIE_FLAVOURS:--12,12,-14,14,-16,16}
xsec=${GENIE_XSEC:-${genie_inputs}/splines/G21_11a_00_000/key4ship/xsec_all_key4ship.xml}
ship_setup=${SHIP_SETUP:-/cvmfs/ship.cern.ch/25.09/setUp.sh}

for f in "${gdml}" "${flux}" "${xsec}" "${ship_setup}"; do
    if [ ! -e "${f}" ]; then
        echo "ERROR: missing input: ${f}" >&2
        exit 1
    fi
done
mkdir -p "${outdir}" "${logdir}"

if [ -n "${KEY4HEP_STACK:-}" ]; then
    echo "ERROR: key4hep is sourced in this shell; run $0 from a clean environment" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${ship_setup}"
if [ -n "${GENIE_CONFIG_DIR:-}" ]; then
    export GXMLPATH=${GENIE_CONFIG_DIR}
fi

cd "${outdir}"
echo "gevgen_fnal: $(which gevgen_fnal) on $(hostname)"
# Messenger_laconic.xml (from $GENIE/config): without it GMCJDriver logs a
# NOTICE for every flux neutrino thrown, ~7500 per interaction here (only ~6%
# of the rays cross the SND), i.e. 1.2 GB of log per 1000 events.
# --event-record-print-level -1: otherwise every generated event's GHEP record
# goes to the log too, ~8 KB each (80 MB per 10k-event trial). The last one
# is still in the .status file, and the events themselves are in the ghep.
gevgen_fnal --tune "${tune}" -r "${run}" -n "${nevents}" --seed "${seed}" \
    -f "${flux},,${flavours}" \
    -g "${gdml}" -L cm -D g_cm3 \
    --cross-sections "${xsec}" \
    --message-thresholds Messenger_laconic.xml \
    --event-record-print-level -1 \
    &> "${logdir}/gevgen_${run}.log"
gntpc -i "gntp.${run}.ghep.root" -f gst -o "gntp.${run}.gst.root" \
    &> "${logdir}/gntpc_${run}.log"
echo "GENIE done: ${outdir}/gntp.${run}.gst.root"
