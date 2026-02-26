import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV

SIM = DD4hepSimulation()

SIM.runType         = "batch"
SIM.numberOfEvents  = 1000
SIM.skipNEvents     = 0

SIM.compactFile = "../geometry/SND_compact.xml"
SIM.outputFile      = "output.edm4hep.root"
print("COMPACT FILE =", SIM.compactFile)

SIM.enableGun       = True
SIM.gun.particle    = "mu-"
SIM.gun.energy      = 50 * GeV
SIM.gun.position    = (0, 0, -500 * mm)
SIM.gun.direction   = (0, 0, 1)

SIM.physicsList     = "QGSP_BERT"
