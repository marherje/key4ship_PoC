from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion, BasicDigitizer
from Gaudi.Configuration import MessageSvc, DEBUG, VERBOSE, INFO

# Initializer
iosvc = IOSvc()
iosvc.Input = "output.edm4hep.root"
iosvc.Output = "snd_digi.edm4hep.root"

# To fix verbose level globally for debugging
# MessageSvc().OutputLevel = DEBUG

# Loading the digitization algorithms
mip_1 = GeV2MIPConversion("GeV2MIP_SiTarget")
mip_1.InputCollection  = "SiTargetHits"
mip_1.OutputCollection = "SiTargetHitsMIP"
mip_1.MIPValue = 0.00009
mip_1.DebugFrequency = 500 # Print debug information every 500 events, this variable is algorithm-specific
mip_1.OutputLevel = DEBUG # Set the output level to DEBUG for this algorithm exclusively

mip_2 = GeV2MIPConversion("GeV2MIP_SiPixel")
mip_2.InputCollection  = "SiPixelHits"
mip_2.OutputCollection = "SiPixelHitsMIP"
mip_2.MIPValue = 0.0002

dig_1 = BasicDigitizer("BasicDigitizer_SiTarget")
dig_1.InputCollection  = "SiTargetHitsMIP"
dig_1.OutputCollection = "SiTargetHitsDigi"
dig_1.Threshold = 0.5

dig_2 = BasicDigitizer("BasicDigitizer_SiPixel")
dig_2.InputCollection  = "SiPixelHitsMIP"
dig_2.OutputCollection = "SiPixelHitsDigi"
dig_2.Threshold = 0.5

# Application manager
ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [mip_1, mip_2, dig_1, dig_2],
    ExtSvc  = [iosvc]
)
