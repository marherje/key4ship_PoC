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

# 6. Event display. In folder /event_display/ do
./launch.sh ../gaudi_jobs/2_mu_pipeline/ShipHits.root \
            ../simulation/geometry/SND_compact.xml 0
# Both paths are required: the compact MUST be the one the hits were simulated
# with, or the hits get snapped to the wrong planes with no error.
# The third argument is the event index, called "window" in the data model
# (window_id, *Windowed collections). "window" and "event" are used
# interchangeably in the framework because the chain supports two modes that
# feed job3/job4/job5 the same collections:
#   - job1_overlay.py (EventOverlay): event-by-event, window i == event i.
#     This is what every committed pipeline runs.
#   - job1_shuffler.py + job2_splitter.py (EventShuffler + EventWindowSplitter):
#     merges several MC sources into one time stream with per-source delays and
#     cuts it into fixed time windows, so a window then holds several events.
# See event_display/README.md for the details, geometry variants and
# --color-by mc.

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