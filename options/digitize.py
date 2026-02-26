from k4FWCore import ApplicationMgr, IOSvc
from Configurables import GeV2MIPConversion

iosvc = IOSvc()
iosvc.Input = "run_scripts/output.edm4hep.root"
iosvc.Output = "snd_digi.edm4hep.root"

alg = GeV2MIPConversion("GeV2MIP_SiTarget")
alg.InputCollection  = "SiTargetHits"
alg.OutputCollection = "SiTargetHitsMIP"
alg.GeV2MIP = 11111.1

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 10,
    TopAlg  = [alg],
    ExtSvc  = [iosvc]
)
