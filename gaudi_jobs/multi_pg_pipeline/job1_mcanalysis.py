from k4FWCore import ApplicationMgr, IOSvc
from Configurables import MCParticleAnalyzer
import os
import sys

# MC-truth analysis for multiple-Particle-Gun samples.
#
# Usage:
#   INPUT_FILE=<edm4hep file from ddsim> k4run job1_mcanalysis.py
#   (optional) OUTPUT_FILE=<histos.root>
#
# Reads the MCParticles collection, identifies the generator primaries and
# fills basic histograms: primary energy, angular distributions, opening
# angle, vertex position, multiplicities and primaries-vs-final-state
# comparison. Calorimeter-level studies (deposited energy, hits, shower
# shape, PID) can be added to MCParticleAnalyzer.cpp later — the sim hit
# collections are already in the input file.

infile = os.environ.get("INPUT_FILE", "")
if not infile:
    print("Use: INPUT_FILE=<input_file> [OUTPUT_FILE=<histos.root>] k4run job1_mcanalysis.py")
    sys.exit(1)

# Default output: histos_<input basename>.root in the current directory.
default_out = "histos_" + os.path.basename(infile).replace(".edm4hep.root", "") + ".root"
outfile = os.environ.get("OUTPUT_FILE", default_out)

iosvc = IOSvc()
iosvc.Input = [infile]

ana = MCParticleAnalyzer("MCParticleAnalyzer")
ana.InputCollection = "MCParticles"
ana.OutputFile      = outfile
ana.MaxEnergy       = 20.0    # GeV — fits the 10+5 GeV configurations
ana.MaxNMCParticles = 2000.0  # e+ showers produce O(1000) secondaries

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [ana],
    ExtSvc  = [iosvc]
)
