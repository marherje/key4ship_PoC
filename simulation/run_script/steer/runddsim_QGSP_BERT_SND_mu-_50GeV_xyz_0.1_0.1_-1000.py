import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV

gun_direction = (0, 0, 1)
gun_position = (0.1 * mm, 0.1 * mm, -1000 * mm)

compact_path = os.path.abspath("/afs/cern.ch/user/m/marquezh/public/key4ship/key4ship_PoC/simulation/run_script/../geometry/SND_compact.xml")

if not os.path.isfile(compact_path):
    raise RuntimeError("ERROR: geometry file not found: " + compact_path)

SIM = DD4hepSimulation()

SIM.runType        = "batch"
SIM.numberOfEvents = 1000
SIM.skipNEvents    = 0

SIM.compactFile = str(compact_path)
SIM._compactFile = SIM.compactFile
SIM.outputFile     = os.path.abspath("/afs/cern.ch/user/m/marquezh/public/key4ship/key4ship_PoC/simulation/run_script/data/output_mu-_xyz_0.1_0.1_-1000_E50.edm4hep.root")

print("COMPACT FILE =", SIM.compactFile)
print("PARTICLE =", "mu-")
print("Energy =", 50)
print("Position =", gun_position)
print("Direction =", gun_direction)

SIM.enableGun      = True
SIM.gun.particle   = "mu-"
SIM.gun.energy     = 50 * GeV
SIM.gun.position   = gun_position
SIM.gun.direction  = gun_direction

SIM.physicsList    = "QGSP_BERT"

