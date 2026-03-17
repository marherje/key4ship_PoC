import cppyy
import ROOT
ROOT.gSystem.Load("libDDCore")

desc = ROOT.dd4hep.Detector.getInstance()
desc.fromXML("../../geometry/SND_compact.xml")

mgr = ROOT.gGeoManager

def walk(node, path=""):
    vol = node.GetVolume()
    name = vol.GetName()
    current_path = path + "/" + node.GetName()
    shape = vol.GetShape()

    is_sens = False
    detector = None
    if "SiTarget" in current_path and "slice_2" in name:
        is_sens = True; detector = "SiTarget_X"
    elif "SiTarget" in current_path and "slice_4" in name:
        is_sens = True; detector = "SiTarget_Y"
    elif "SiPad" in current_path and "slice4" in name:
        is_sens = True; detector = "SiPad"

    if is_sens:
        ok = mgr.cd(current_path)
        if ok:
            t = mgr.GetCurrentMatrix().GetTranslation()
            s = shape
            print(f"{detector:<12} {name:<45} "
                  f"z={t[2]:8.2f}  "
                  f"half_xy=({s.GetDX():.1f},{s.GetDY():.1f})  "
                  f"half_z={s.GetDZ():.2f} cm")

    for i in range(node.GetNdaughters()):
        walk(node.GetDaughter(i), current_path)

walk(mgr.GetTopNode())
