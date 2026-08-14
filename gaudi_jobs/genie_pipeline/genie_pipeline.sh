#!/bin/bash
#
# GENIE-driven reconstruction chain: same job3/job4/job5 as the
# particle-gun pipelines (digitization/tracking/RNTuple only care about the
# *Windowed hit collections, not about what generated the primaries).
#
# Run from this directory: all job file references are relative.
#
#   INPUT_FILE=<ddsim GENIE output> ./genie_pipeline.sh
#
# INPUT_FILE defaults to the sample produced by
# simulation/run_script/launch_genie.sh (see its --output).

k4run job1_overlay.py
sleep 1
INPUT_FILE=events.edm4hep.root k4run job3_digitize.py
sleep 1
k4run job4_tracking.py
sleep 1
k4run job5_rntuple.py
