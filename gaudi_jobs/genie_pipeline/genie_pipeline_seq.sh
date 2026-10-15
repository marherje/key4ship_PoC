#!/bin/bash
#
# GENIE-driven reconstruction chain with ACTSSequentialTracker (per-subdetector
# search + splicing) instead of the global ACTSProtoTracker used by
# genie_pipeline.sh. Writes to its own files (tracks_seq.edm4hep.root,
# ShipHits_seq.root) so both trackers can be run on the same events and
# compared.
#
# Run from this directory: all job file references are relative.
#
#   INPUT_FILE=<ddsim GENIE output> ./genie_pipeline_seq.sh
#
# INPUT_FILE defaults to the sample produced by
# simulation/run_script/launch_genie.sh (see its --output).

k4run job1_overlay.py
sleep 1
INPUT_FILE=events.edm4hep.root k4run job3_digitize.py
sleep 1
k4run job4_tracking_seq.py
sleep 1
k4run job5_rntuple_seq.py
