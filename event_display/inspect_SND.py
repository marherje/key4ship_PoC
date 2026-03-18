import cppyy
import ROOT
ROOT.gSystem.Load("libDDCore")

desc = ROOT.dd4hep.Detector.getInstance()
desc.fromXML("../simulation/geometry/SND_compact.xml")

mgr = ROOT.gGeoManager

z_lists = {"SiPad": [], "SiTarget_Y": [], "SiTarget_Z": []}

def walk(node, path=""):
    vol = node.GetVolume()
    name = vol.GetName()
    current_path = path + "/" + node.GetName()
    shape = vol.GetShape()

    is_sens = False
    detector = None
    if "SiTarget" in current_path and "slice_2" in name:
        is_sens = True; detector = "SiTarget_Y"
    elif "SiTarget" in current_path and "slice_4" in name:
        is_sens = True; detector = "SiTarget_Z"
    elif "SiPad" in current_path and "_slice_4" in name:
        is_sens = True; detector = "SiPad"

    if is_sens:
        ok = mgr.cd(current_path)
        if ok:
            t = mgr.GetCurrentMatrix().GetTranslation()
            s = shape
            print(f"{detector:<12} {name:<45} "
                  f"z={t[2]:8.2f}  "
                  f"half_yz=({s.GetDX():.1f},{s.GetDY():.1f})  "
                  f"half_x={s.GetDZ():.2f} cm")
            z_lists[detector].append(round(t[2], 2))

    for i in range(node.GetNdaughters()):
        walk(node.GetDaughter(i), current_path)

walk(mgr.GetTopNode())

print()
var_names = {"SiPad": "SIPAD_X", "SiTarget_Y": "SITARGET_X_Y", "SiTarget_Z": "SITARGET_X_Z"}
for det, zs in z_lists.items():
    vals = ", ".join(f"{z:.2f}" for z in zs)
    # wrap at ~60 chars
    entries = [f"{z:.2f}" for z in zs]
    lines, line = [], []
    for e in entries:
        line.append(e)
        if len(", ".join(line)) > 55:
            lines.append("    " + ", ".join(line[:-1]) + ",")
            line = [line[-1]]
    if line:
        lines.append("    " + ", ".join(line) + ",")
    print(f"{var_names[det]} = [")
    print("\n".join(lines))
    print("]")
