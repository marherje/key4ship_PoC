#!/usr/bin/env bash

# To use only if the project is already compile and you want to play with the steering file, for example. 
# If you want to compile the project, please follow the instructions in the README.md file.
#  -r 2026-02-01
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01

# Repo root = this script's own directory, so the install paths are correct no
# matter where the script is sourced from (previously used $PWD, which broke
# when sourced from a subdirectory like event_display/ or simulation/run_script/).
_KEY4SHIP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LD_LIBRARY_PATH=$_KEY4SHIP_ROOT/install/lib64:$_KEY4SHIP_ROOT/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$_KEY4SHIP_ROOT/install/lib64:$_KEY4SHIP_ROOT/install/lib:$_KEY4SHIP_ROOT/install/python:$PYTHONPATH
unset _KEY4SHIP_ROOT
