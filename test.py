"""
DD4hep Simulation Steering Script for SND Detector
===================================================

This script configures and runs a DD4hep simulation for the SND (Scattering and Neutrino Detector).
It simulates particle beams (e.g., muons) interacting with the detector geometry.

Usage:
    source init_key4ship.sh && ddsim --steeringFile test.py

Configuration:
    - Modify the PARAMETERS section below to change particle type, energy, position, and direction
    - Output files are saved to the 'output/' directory with names encoding all simulation parameters
    - The compact geometry file path can be adjusted in the PARAMETERS section

Output:
    - ROOT file containing EDM4hep format event data
"""

import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV

# ===== PARAMETERS =====
# Gun configuration
particle = "mu-"
energy = 5  # GeV
gun_position = (5 * mm, -82.5 * mm, 750 * mm)
gun_direction = (0, 0, 1)

# Simulation configuration
compact_path = "simulation/geometry/SND_compact.xml"
numberOfEvents = 100
skipNEvents = 0

# Physics configuration
physicsList = "QGSP_BERT"
# =====================

# Validate geometry file
if not os.path.isfile(compact_path):
    raise RuntimeError("ERROR: geometry file not found: " + compact_path)

# Create output directory if it doesn't exist
os.makedirs("output", exist_ok=True)

# Initialize simulation
SIM = DD4hepSimulation()

SIM.runType        = "run"
SIM.numberOfEvents = numberOfEvents
SIM.skipNEvents    = skipNEvents
SIM.compactFile    = str(compact_path)
SIM._compactFile   = SIM.compactFile

# Generate output filename from parameters
x, y, z = gun_position[0] / mm, gun_position[1] / mm, gun_position[2] / mm
dir_x, dir_y, dir_z = gun_direction

SIM.outputFile = (
    f"output/output_{particle}_xyz_{x}_{y}_{z}_"
    f"dir_{dir_x}_{dir_y}_{dir_z}_"
    f"E{energy}_{numberOfEvents}_events.edm4hep.root"
)

# Print configuration
print("COMPACT FILE =", SIM.compactFile)
print("PARTICLE =", particle)
print("Energy =", energy, "GeV")
print("Position =", gun_position)
print("Direction =", gun_direction)
print("OUTPUT FILE =", SIM.outputFile)
print()

# Configure gun
SIM.enableGun      = True
SIM.gun.particle   = particle
SIM.gun.energy     = energy * GeV
SIM.gun.position   = gun_position
SIM.gun.direction  = gun_direction

# Configure physics
SIM.physicsList    = physicsList
SIM.action.mapActions["MTC"] = "SND_SciFiAction"