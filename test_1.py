"""
DD4hep Simulation Steering Script for SND Detector
===================================================

Usage:
    source init_key4ship.sh && ddsim --steeringFile test_1.py \
        --compactFile /absolute/path/to/SND_compact.xml

Note: compactFile is intentionally NOT set here.
      Pass it via --compactFile on the command line so ddsim's
      argument parser owns it and never resets it to empty.

Output:
    - ROOT file containing EDM4hep format event data in output/ directory
"""

import os
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV


class DD4hepSimulationWithGDML(DD4hepSimulation):
    """Subclass that exports GDML immediately after geometry is loaded."""

    def __init__(self, gdml_output_path="detector_geometry.gdml"):
        super().__init__()
        self._gdml_output_path = gdml_output_path

    def _setupGeant4(self):
        try:
            mgr = self.detectorDescription.manager()
            mgr.Export(self._gdml_output_path)
            print(f"[GDML] Geometry exported to: {self._gdml_output_path}")
        except Exception as exc:
            print(f"[GDML] WARNING: export failed — {exc}")
        super()._setupGeant4()


# ===== PARAMETERS =====
particle      = "mu-"
energy        = 5        # GeV
gun_position  = (5 * mm, -82.5 * mm, 750 * mm)
gun_direction = (0, 0, 1)

numberOfEvents = 100
skipNEvents    = 0
physicsList    = "QGSP_BERT"
# =====================

os.makedirs("output", exist_ok=True)

# ── Initialise simulation ──────────────────────────────────────────────────
# compactFile is NOT set here — it is passed via --compactFile on the CLI.
# ddsim's parseOptions() reads sys.argv AFTER exec()ing this file, so any
# SIM.compactFile assignment here gets overwritten by the CLI parser when
# --compactFile is present, or reset to "" when it is absent.
SIM = DD4hepSimulationWithGDML(gdml_output_path="detector_geometry.gdml")

SIM.runType        = "run"
SIM.numberOfEvents = numberOfEvents
SIM.skipNEvents    = skipNEvents

# ── Output filename ────────────────────────────────────────────────────────
x, y, z             = gun_position[0] / mm, gun_position[1] / mm, gun_position[2] / mm
dir_x, dir_y, dir_z = gun_direction

SIM.outputFile = (
    f"output/output_{particle}_xyz_{x}_{y}_{z}_"
    f"dir_{dir_x}_{dir_y}_{dir_z}_"
    f"E{energy}_{numberOfEvents}_events.edm4hep.root"
)

# ── Print configuration ────────────────────────────────────────────────────
print("PARTICLE     =", particle)
print("Energy       =", energy, "GeV")
print("Position     =", gun_position)
print("Direction    =", gun_direction)
print("OUTPUT FILE  =", SIM.outputFile)
print()

# ── Particle gun ───────────────────────────────────────────────────────────
SIM.enableGun     = True
SIM.gun.particle  = particle
SIM.gun.energy    = energy * GeV
SIM.gun.position  = gun_position
SIM.gun.direction = gun_direction

# ── Physics ────────────────────────────────────────────────────────────────
SIM.physicsList = physicsList
SIM.action.mapActions["MTC"] = "SND_SciFiAction"

# ── Run ────────────────────────────────────────────────────────────────────
SIM.run()