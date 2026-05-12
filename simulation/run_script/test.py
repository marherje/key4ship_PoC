import cppyy
import ROOT
ROOT.gSystem.Load("libDDCore")

decoder = ROOT.dd4hep.DDSegmentation.BitFieldCoder("system:8,layer:8,slice:4,plane:1,strip:14")

# Try new podio API
try:
    from podio import reading
    reader = reading.get_reader("data/output_e-_xyz_1_-82.5_-1000_dir_0_0.05_1_E50.edm4hep.root")
except:
    import podio
    reader = podio.root_io.Reader("data/output_e-_xyz_1_-82.5_-1000_dir_0_0.05_1_E50.edm4hep.root")

for frame in reader.get("events"):
    hits = frame.get("SiTargetHits")
    for hit in hits:
        cid   = hit.getCellID()
        layer = decoder.get(cid, "layer")
        plane = decoder.get(cid, "plane")
        strip = decoder.get(cid, "strip")
        pos   = hit.getPosition()
        print(f"layer={layer} plane={plane} strip={strip} "
              f"x={pos.x:.1f} y={pos.y:.1f} z={pos.z:.1f}")
    break
