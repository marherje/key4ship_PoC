from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion, BasicDigitizer
from Gaudi.Configuration import MessageSvc, DEBUG
import os

# Take file name from environment variable
infile = os.environ.get("INPUT_FILE", "")
if not infile:
    print("Use: INPUT_FILE=<input_file> k4run job1_digitize.py")
    import sys; sys.exit(1)

# Generate output file name based on input file name
outfile = "digitized.edm4hep.root"

# IOSvc
iosvc = IOSvc()
iosvc.Input = infile
iosvc.Output = outfile

# Digi SiTarget
mip_1 = GeV2MIPConversion("GeV2MIP_SiTarget")
mip_1.InputCollection  = "SiTargetHitsWindowed"
mip_1.OutputCollection = "SiTargetHitsMIP"
mip_1.MIPValue = 0.00009
mip_1.DebugFrequency = 500
mip_1.OutputLevel = DEBUG

# Digi SiPad
mip_2 = GeV2MIPConversion("GeV2MIP_SiPad")
mip_2.InputCollection  = "SiPadHitsWindowed"
mip_2.OutputCollection = "SiPadHitsMIP"
mip_2.MIPValue = 0.0002

# MIP-cut Digitizer SiTarget
dig_1 = BasicDigitizer("BasicDigitizer_SiTarget")
dig_1.InputCollection  = "SiTargetHitsMIP"
dig_1.OutputCollection = "SiTargetHitsDigi"
dig_1.Threshold = 0.5
dig_1.DebugFrequency = 500

# MIP-cut Digitizer SiPad
dig_2 = BasicDigitizer("BasicDigitizer_SiPad")
dig_2.InputCollection  = "SiPadHitsMIP"
dig_2.OutputCollection = "SiPadHitsDigi"
dig_2.Threshold = 0.5
dig_2.DebugFrequency = 500

# Application Manager
ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [mip_1, mip_2, dig_1, dig_2],
    ExtSvc  = [iosvc]
)
