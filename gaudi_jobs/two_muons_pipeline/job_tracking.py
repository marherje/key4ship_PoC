from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, SiTargetMeasConverter, \
                          SiPadMeasConverter, ACTSProtoTracker
from Gaudi.Configuration import DEBUG, INFO

iosvc = IOSvc()
iosvc.Input  = "timewindows.edm4hep.root"
iosvc.Output = "tracks.root"

geosvc = ACTSGeoSvc("ACTSGeoSvc")
geosvc.CompactFile = "../../simulation/geometry/SND_compact.xml"
geosvc.OutputLevel = INFO

sitarget_conv = SiTargetMeasConverter("SiTargetMeasConverter")
sitarget_conv.GeoSvc           = "ACTSGeoSvc"
sitarget_conv.InputCollection  = "SiTargetHitsWindowed"
sitarget_conv.OutputCollection = "SiTargetMeasurements"
sitarget_conv.BitField         = "system:8,layer:8,slice:4,plane:1,strip:14"
sitarget_conv.StripPitch       = 0.0755
sitarget_conv.OutputLevel      = DEBUG

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = "SiPadHitsWindowed"
sipad_conv.OutputCollection = "SiPadMeasurements"
sipad_conv.BitField         = "system:8,layer:8,slice:4,y:9,z:9"
sipad_conv.PixelSizeY       = 5.5
sipad_conv.PixelSizeZ       = 5.5
sipad_conv.OutputLevel      = DEBUG

proto = ACTSProtoTracker("ACTSProtoTracker")
proto.GeoSvc           = "ACTSGeoSvc"
proto.InputSiTarget    = "SiTargetMeasurements"
proto.InputSiPad       = "SiPadMeasurements"
proto.OutputCollection = "ACTSTracks"
proto.BFieldX          = 0.0
proto.BFieldY          = 0.0
proto.BFieldZ          = 0.0
proto.OutputLevel      = DEBUG

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 10,
    TopAlg  = [sitarget_conv, sipad_conv, proto],
    ExtSvc  = [iosvc, geosvc]
)
