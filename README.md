# SND_sim proof-of-concept of full key4hep sim+reco chain: dd4hep + ddsim + gaudi

SND PoC with STarget and SiPad — DD4hep v01-35 / key4hep 2026-02-01

## Target stack (key4hep latest release in Feb 2026)
- DD4hep: v01-35
- ddsim: integrated in DD4hep v01-35
- Geant4: 11.x
- Gaudi: v4x series
- ROOT: 6.36+
- Platform: lxplus.cern.ch, AlmaLinux 9
- Source: /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-02-01

## Build and run

```bash
# 1. Source key4hep
source /cvmfs/sw.hsf.org/key4hep/setup.sh

# 2. Build plugin
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" \
      -DCMAKE_INSTALL_PREFIX=../install \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make install -j4
cd ..

# 3. Expose plugins to DD4hep & Gaudi algorithms to python
export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/install/lib64:$PWD/install/lib:$PWD/install/python:$PYTHONPATH

# Example with 2 muons:
# 4. Run simulation (ddsim)
cd simulation/run_script
source launch_PG.sh
cd ../..

# 5. Run gaudi algorithm/s
cd gaudi_jobs/2_mu_pipeline
source 2_mu_pipeline.sh
# This bash script can do the whole pipeline by itself
# Different cases to test in /gaudi_jobs/ 

# For using the event display.
# In folder /eventdisplay/ do
python event_display_eve.py --hits ../gaudi_jobs/2_mu_pipeline/ShipHits.root --window 0 
# where window is the "event"

```

## Git hooks

`.githooks/commit-msg` strips AI self-attribution (`Co-Authored-By` trailers
naming an AI vendor, `Claude-Session` / `Codex-Session` trailers, "Generated
with ..." footers) from commit messages. Human co-authors are untouched.

Git does not enable hooks that come with a clone, so activate it once per
working copy:

```bash
git config core.hooksPath .githooks
```

This is repo-local and overrides a global `core.hooksPath`. It redirects
*all* hooks to that directory, so anything you keep only in `.git/hooks`
stops running — move it to `.githooks/` if you need both.