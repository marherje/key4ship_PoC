#!/bin/bash
#
# Same pipeline as muon_pipeline.sh, but tracking with ACTSSequentialTracker
# (per-subdetector search + splicing) instead of the global ACTSProtoTracker.
# Writes to its own files — tracks_seq.edm4hep.root and ShipHits_seq.root — so
# both trackers can be run on the same events and compared side by side.
#
# Run from this directory: all job file references are relative.

k4run job1_overlay.py
sleep 1
INPUT_FILE=events.edm4hep.root k4run job3_digitize.py
sleep 1
k4run job4_tracking_seq.py
sleep 1
k4run job5_rntuple_seq.py
