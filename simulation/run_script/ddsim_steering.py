import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV

particles = ["mu-", "e-"]
positions = [0, 50]  # en mm
numberOfEvents = 1000
gun_energy    = 50 * GeV
gun_direction = (0, 0, 1)

compact_path = f"../geometry/SND_compact.xml"
if not os.path.exists(compact_path):
    raise RuntimeError(f"ERROR: No se encuentra el archivo XML: {compact_path}")

# Loop sobre posiciones y partículas
for pos in positions:
    gun_position = (0, 0, -1000 * mm + pos * mm)
    for particle in particles:
        SIM = DD4hepSimulation()

        SIM.runType        = "batch"
        SIM.numberOfEvents = numberOfEvents
        SIM.skipNEvents    = 0

        SIM.compactFile    = compact_path
        SIM.outputFile     = f"output_{particle}_y{pos}.edm4hep.root"

        print(f"COMPACT FILE = {SIM.compactFile}, PARTICLE = {particle}")
        print(f"Energy = {gun_energy}, Position = {gun_position}, Direction = {gun_direction}")

        SIM.enableGun      = True
        SIM.gun.particle   = particle
        SIM.gun.energy     = gun_energy
        SIM.gun.position   = gun_position
        SIM.gun.direction  = gun_direction

        SIM.physicsList    = "QGSP_BERT"

        SIM.run()