from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion, BasicDigitizer, SciFiDigitizer
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

# # MIP SciFi
mip_3 = GeV2MIPConversion("GeV2MIP_MTCSciFi")
mip_3.InputCollection  = "MTCSciFiHitsWindowed"
mip_3.OutputCollection = "MTCSciFiHitsMIP"
mip_3.MIPValue = 2e-3*1.35

# SciFi Digitizer: light attenuation + Poisson + SiPM saturation + QDC
# Geometry and physics parameters from SND_compact.xml / FairRoot MTCDetHit model.
# Geometry properties are here so they can be tuned without recompiling.
scifi_digi = SciFiDigitizer("SciFiDigitizer_MTCSciFi")
scifi_digi.InputCollection   = "MTCSciFiHitsWindowed"
scifi_digi.OutputCollection  = "MTCSciFiHitsDigi"
# Fiber geometry (from SND_compact.xml)
scifi_digi.FiberAngleDeg  = 5.0                     # MTC_fiber_angle_deg
scifi_digi.EnvHeightHalf  = [200.0, 250.0, 300.0]   # mm; MTC40/50/60 env_height / 2
scifi_digi.SiPMSide       = +1                      # SiPM at +y end
# Attenuation model (FairRoot MTCDetHit defaults)
scifi_digi.AttenuationLength = 300.0   # cm (lambda)
scifi_digi.AttenuationOffset =  20.0   # cm (x0)
scifi_digi.PhotonsPerGeV     = 1.6e5   # = energy[GeV] * 1e6 * 0.16 photons
# scifi_digi.MirrorReflectivity = 0.9
scifi_digi.MirrorReflectivity = 0.0
# SiPM + readout model
scifi_digi.SiPMNpixels     = 104.0
scifi_digi.PhotonThreshold =   3.5
scifi_digi.QDC_A           =   0.172
scifi_digi.QDC_B           =  -1.31
scifi_digi.QDC_sigmaA      =   0.006
scifi_digi.QDC_sigmaB      =   0.33
scifi_digi.DebugFrequency  =   1

# Digi MTC Scintillator (15 mm thick plastic pads — plain MIP threshold, no fiber model)
mip_4 = GeV2MIPConversion("GeV2MIP_MTCScint")
mip_4.InputCollection  = "MTCScintHitsWindowed"
mip_4.OutputCollection = "MTCScintHitsMIP"
mip_4.MIPValue = 0.003

dig_4 = BasicDigitizer("BasicDigitizer_MTCScint")
dig_4.InputCollection  = "MTCScintHitsMIP"
dig_4.OutputCollection = "MTCScintHitsDigi"
dig_4.Threshold = 0.5


# Application Manager
ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [mip_1, mip_2, dig_1, dig_2, mip_3, scifi_digi, mip_4, dig_4],
    ExtSvc  = [iosvc]
)
