import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV
gun_direction = (0, 0, 1) 
gun_position = (1 * mm, -42.5 * mm, -1000 * mm)

compact_path = "simulation/geometry/SND_compact.xml"

if not os.path.isfile(compact_path):
    raise RuntimeError("ERROR: geometry file not found: " + compact_path)

SIM = DD4hepSimulation()

SIM.runType        = "run"
SIM.numberOfEvents = 10
SIM.skipNEvents    = 0

SIM.compactFile = str(compact_path)
SIM._compactFile = SIM.compactFile

particle = "mu-"
energy = 75

SIM.outputFile     = f"output_{particle}_xyz_1_-42.5_-1000_dir_0_0_1_E{energy}.edm4hep.root"


print("COMPACT FILE =", SIM.compactFile)
print("PARTICLE =", particle)
print("Energy =", energy)
print("Position =", gun_position)
print("Direction =", gun_direction)

SIM.enableGun      = True
SIM.gun.particle   = particle
SIM.gun.energy     = energy * GeV
SIM.gun.position   = gun_position
SIM.gun.direction  = gun_direction

SIM.physicsList    = "QGSP_BERT"

SIM.action.mapActions["MTC"] = "SND_SciFiAction"

