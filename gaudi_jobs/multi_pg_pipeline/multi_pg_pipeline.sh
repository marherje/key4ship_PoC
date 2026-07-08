#!/bin/bash
#
# Analyse every multiple-Particle-Gun sample found in run_script/data/.
#
# Runs job1_mcanalysis.py once per output_<physlist>_SND_*_angle_*.edm4hep.root
# file and writes one histos_<label>.root per sample in this directory.
#
# Usage (after building the repo and sourcing the key4hep environment,
# e.g. via ../../init_key4ship.sh or build.sh exports):
#   ./multi_pg_pipeline.sh

set -e

data_dir="../../simulation/run_script/data"

shopt -s nullglob
samples=("${data_dir}"/output_*GeV_*GeV_angle_*.edm4hep.root)
shopt -u nullglob

if [ ${#samples[@]} -eq 0 ]; then
    echo "No multiple-PG samples found in ${data_dir}."
    echo "Produce them with simulation/run_script/launch_multiplePG.sh first."
    exit 1
fi

for f in "${samples[@]}"; do
    echo "=== Analysing $(basename "$f")"
    INPUT_FILE="$f" k4run job1_mcanalysis.py
done

echo "Done. Histogram files:"
ls -1 histos_*.root
