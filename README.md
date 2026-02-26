# SND_sim proof-of-concept of full key4hep sim+reco chain: dd4hep + ddsim + gaudi

SND PoC with STarget and SiPixel — DD4hep v01-35 / key4hep 2026-02-01

## Target stack (key4hep latest release)
- DD4hep: v01-35
- ddsim: integrated in DD4hep v01-35
- Geant4: 11.x
- Gaudi: v3x series
- ROOT: 6.36+
- Platform: lxplus.cern.ch, AlmaLinux 9
- Source: /cvmfs/sw.hsf.org/key4hep/setup.sh

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

# 4. Verify plugin symbols are exported
nm -D install/lib64/libSiPixelDetector_plugin.so | grep " T "
# Must return non-empty output

# 5. Run ddsim
cd run_scripts
ddsim --steeringFile ddsim_steering.py 

# 6. Check hits
root -l output.edm4hep.root -e \
  'cout << events->GetMaximum("ECALHits@.size()") << endl'
# Must return > 0
```

# 7. Run algorithm/s
k4run digitize.py 
