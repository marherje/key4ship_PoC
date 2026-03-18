#!/usr/bin/env python3
"""
SND Event Display using ROOT TEve.
Geometry is fully hardcoded (no DD4hep compact XML or plugin libraries needed).
Hits are read from an RNTuple file and displayed as TEveGeoShape boxes.

Usage:
    python event_display_eve.py --hits ../ShipHits.root --window 0 \
        --source-labels "0:mixed" "1:Sig1" "2:Sig2"
"""

import argparse
import ROOT

# ---------------------------------------------------------------------------
# 1. CONSTANTS  (all Z values and sizes in cm — TGeo/TEve use cm)
# ---------------------------------------------------------------------------

# Copy here the output from inspect_SND.py                                                                                                                       
SIPAD_X = [                                                                                                                                                      
    -6.70, -5.25, -3.80, -2.35, -0.90, 0.55, 2.00, 3.45,                                                                                                         
    4.90, 6.35,                                                                                                                                                  
]                                                                                                                                                                
SITARGET_Y_X = [                                                                                                                                                 
    -18.79, -17.69, -16.59, -15.49, -14.39, -13.29, -12.19,                                                                                                      
    -11.09, -9.99, -8.89,                                                                                                                                        
]                                                                                                                                                                
SITARGET_Z_X = [                                                                                                                                                 
    -18.21, -17.11, -16.01, -14.91, -13.81, -12.71, -11.61,                                                                                                      
    -10.51, -9.41, -8.31,                                                                                                                                        
]     

# Source ID → (R, G, B) in 0-1 range
SOURCE_RGB = {
    0: (0.53, 0.53, 0.53),   # mixed  → grey
    1: (0.13, 0.59, 0.95),   # signal → blue
    2: (0.96, 0.26, 0.21),   # bkg 1  → red
    3: (0.30, 0.69, 0.31),   # bkg 2  → green
    4: (1.00, 0.60, 0.00),   # bkg 3  → orange
    5: (0.61, 0.15, 0.69),   # bkg 4  → purple
}

# Populated after ROOT is loaded (see bottom of module)
source_colors = {}

# ---------------------------------------------------------------------------
# 2. C++ NAMESPACE FOR RNTUPLE DATA
# Kept in C++ scope to avoid the PyROOT RNTuple destructor segfault
# (RColumn::~RColumn crashes when Python GC collects the RNTupleReader).
# ---------------------------------------------------------------------------
ROOT.gROOT.ProcessLine("""
namespace SNDEve {
    std::vector<float> gX, gY, gZ, gE;
    std::vector<int>   gSrc, gPlane;
}
""")


# ---------------------------------------------------------------------------
# 3. HIT DATA READING
# ---------------------------------------------------------------------------
def read_hits(hits_file, ntuple_name, window_id):
    """
    Read hits for one detector/window from an RNTuple file.
    All RNTuple objects live inside a C++ block (ProcessLine) so Python's GC
    never touches them.  Coordinates are converted from mm to cm on return.

    SiTarget has a bf_plane field (0=X strip, 1=Y strip).
    SiPad does not — gPlane is set to -1 for all SiPad entries.

    Returns: xs, ys, zs, Es, srcs, planes  (plain Python lists, coords in cm)
    """
    if ntuple_name == "SiTarget":
        plane_view_decl = 'auto vplane = reader->GetView<int>("bf_plane");'
        plane_push      = "SNDEve::gPlane.push_back(vplane(i));"
    else:
        plane_view_decl = ""
        plane_push      = "SNDEve::gPlane.push_back(-1);"

    ROOT.gROOT.ProcessLine(f"""
    {{
        SNDEve::gX.clear();
        SNDEve::gY.clear();
        SNDEve::gZ.clear();
        SNDEve::gE.clear();
        SNDEve::gSrc.clear();
        SNDEve::gPlane.clear();

        auto reader = ROOT::RNTupleReader::Open("{ntuple_name}", "{hits_file}");

        auto vwin = reader->GetView<int>("window_id");
        auto vsrc = reader->GetView<int>("source_id");
        auto vx   = reader->GetView<float>("x");
        auto vy   = reader->GetView<float>("y");
        auto vz   = reader->GetView<float>("z");
        auto vE   = reader->GetView<float>("E");
        {plane_view_decl}

        for (auto i : reader->GetEntryRange()) {{
            if (vwin(i) != {window_id}) continue;
            SNDEve::gX.push_back(vx(i));
            SNDEve::gY.push_back(vy(i));
            SNDEve::gZ.push_back(vz(i));
            SNDEve::gE.push_back(vE(i));
            SNDEve::gSrc.push_back(vsrc(i));
            {plane_push}
        }}

        printf("[reader] %s window {window_id}: %zu hits\\n",
               "{ntuple_name}", SNDEve::gX.size());
    }}
    """)

    mm2cm = 0.1
    n      = ROOT.SNDEve.gX.size()
    xs     = [ROOT.SNDEve.gX[i] * mm2cm for i in range(n)]
    ys     = [ROOT.SNDEve.gY[i] * mm2cm for i in range(n)]
    zs     = [ROOT.SNDEve.gZ[i] * mm2cm for i in range(n)]
    Es     = [ROOT.SNDEve.gE[i]         for i in range(n)]
    srcs   = [ROOT.SNDEve.gSrc[i]       for i in range(n)]
    planes = [ROOT.SNDEve.gPlane[i]     for i in range(n)]

    return xs, ys, zs, Es, srcs, planes


# ---------------------------------------------------------------------------
# 4. NEAREST LAYER SNAPPING
# ---------------------------------------------------------------------------
def nearest_z(hit_x_cm, plane_x_list):
    """Return the X of the closest sensitive plane to hit_x_cm."""
    return min(plane_x_list, key=lambda x: abs(x - hit_x_cm))


# ---------------------------------------------------------------------------
# 5. TEveGeoShape BOX HELPER
# ---------------------------------------------------------------------------
def make_box(name, dx, dy, dz, x, y, z, color, transparency=0):
    """
    Create a TEveGeoShape containing a TGeoBBox.

    dx, dy, dz : HALF-sizes of the box in cm.
    x, y, z    : center position in cm.
    color      : ROOT color index.
    transparency: 0 = opaque, 90 = nearly transparent.
    """
    shape = ROOT.TGeoBBox(dx, dy, dz)
    gs = ROOT.TEveGeoShape(name)
    gs.SetShape(shape)
    gs.SetMainColor(color)
    gs.SetMainTransparency(transparency)
    gs.RefMainTrans().SetPos(x, y, z)
    gs.ResetBBox()
    return gs


# ---------------------------------------------------------------------------
# 6. GEOMETRY BUILDER
# ---------------------------------------------------------------------------
def build_geometry(eve):
    """
    Create all detector planes as semi-transparent background shapes and add
    them to the GlobalScene.
    """
    geo_list = ROOT.TEveElementList("Geometry")

    # SiPad planes
    sipad_list = ROOT.TEveElementList("SiPad planes")
    for i, x in enumerate(SIPAD_X):
        gs = make_box(f"SiPad_{i}", 17.6, 17.6, 0.0325,
                      x, 0.0, 0.0,
                      ROOT.kGreen - 6, transparency=90)
        sipad_list.AddElement(gs)
    geo_list.AddElement(sipad_list)

    # SiTarget Y planes
    sitgt_y_list = ROOT.TEveElementList("SiTarget Y planes")
    for i, x in enumerate(SITARGET_Y_X):
        gs = make_box(f"SiTarget_Y_{i}", 20.0, 20.0, 0.015,
                      x, 0.0, 0.0,
                      ROOT.kAzure + 7, transparency=90)
        sitgt_y_list.AddElement(gs)
    geo_list.AddElement(sitgt_y_list)

    # SiTarget Z planes
    sitgt_z_list = ROOT.TEveElementList("SiTarget Z planes")
    for i, x in enumerate(SITARGET_Z_X):
        gs = make_box(f"SiTarget_Z_{i}", 20.0, 20.0, 0.015,
                      x, 0.0, 0.0,
                      ROOT.kAzure + 9, transparency=90)
        sitgt_z_list.AddElement(gs)
    geo_list.AddElement(sitgt_z_list)

    eve.GetGlobalScene().AddElement(geo_list)
    print("[Geometry] Added SiPad (10) + SiTarget X (10) + SiTarget Y (10) + SiTarget Z (10) planes")


# ---------------------------------------------------------------------------
# 7. HIT BUILDER
# ---------------------------------------------------------------------------
def build_hits(eve, hits_file, window_id, source_labels):
    """
    Read hits for SiTarget and SiPad, create one TEveGeoShape per hit, group
    by (detector, source_id) using TEveElementList, and add to EventScene.

    NOTE: For windows with thousands of hits, replace individual TEveGeoShape
    objects with batched TEveStraightLineSet or TEvePointSet per group here.
    """
    hits_root = ROOT.TEveElementList(f"Hits window {window_id}")

    for det in ["SiTarget", "SiPad"]:
        xs, ys, zs, Es, srcs, planes = read_hits(hits_file, det, window_id)
        if not xs:
            print(f"[Hits] {det}: no hits in window {window_id}")
            continue

        print(f"[Hits] {det}: {len(xs)} hits, "
              f"source_ids={sorted(set(srcs))}")

        # Group indices by source_id
        by_src = {}
        for i, sid in enumerate(srcs):
            by_src.setdefault(sid, []).append(i)

        det_list = ROOT.TEveElementList(det)

        for sid, indices in sorted(by_src.items()):
            label  = source_labels.get(sid, f"src{sid}")
            color  = source_colors.get(sid, ROOT.kWhite)
            grp    = ROOT.TEveElementList(f"{det} [{label}]")

            for i in indices:
                xc, yc, zc = xs[i], ys[i], zs[i]

                if det == "SiPad":
                    snap_x = nearest_x(xc, SIPAD_X)
                    gs = make_box(f"SiPad_hit_{i}",
                                  0.275, 0.275, 0.0325,
                                  snap_x, yc, zc,
                                  color, transparency=0)

                else:  # SiTarget
                    plane = planes[i]
                    if plane == 0:   # X strip: precise X, full Y span
                        snap_x = nearest_x(xc, SITARGET_Z_X)
                        gs = make_box(f"SiTgt_Z_hit_{i}",
                                      0.015,0.003775, 20.0,
                                      snap_x, 0.0, zc,
                                      color, transparency=0)
                    else:            # Y strip (plane==1): precise Y, full X span
                        snap_x = nearest_x(zc, SITARGET_Y_X)
                        gs = make_box(f"SiTgt_Y_hit_{i}",
                                      20.0, 0.015,0.003775,
                                      0.0, yc, snap_x,
                                      color, transparency=0)

                grp.AddElement(gs)

            det_list.AddElement(grp)
            print(f"[Hits]   {det} [{label}]: {len(indices)} hits")

        hits_root.AddElement(det_list)

    eve.GetEventScene().AddElement(hits_root)


# ---------------------------------------------------------------------------
# 8. MAIN
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="SND Event Display (TEve)")
    parser.add_argument("--hits",   default="../ShipHits.root")
    parser.add_argument("--window", type=int, default=0)
    parser.add_argument("--source-labels", nargs="*",
                        default=["0:mixed", "1:signal", "2:background"],
                        metavar="ID:LABEL",
                        help="source_id labels e.g. '0:mixed' '1:Sig1'")
    args = parser.parse_args()

    source_labels = {}
    for item in args.source_labels:
        sid, lbl = item.split(":", 1)
        source_labels[int(sid)] = lbl

    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gSystem.Load("libEve")
    ROOT.gSystem.Load("libEG")
    ROOT.gSystem.Load("libROOTNTuple")

    # Populate source_colors now that ROOT is loaded
    for sid, (r, g, b) in SOURCE_RGB.items():
        idx = ROOT.TColor.GetFreeColorIndex()
        ROOT.TColor(idx, r, g, b)
        source_colors[sid] = idx

    # Create TEveManager with full interface.
    # Add ALL elements before calling ProcessEvents() — once the GL draw
    # cycle starts it acquires DrawLocks on the scenes, and subsequent
    # AddElement calls would produce ModifyLock errors.
    eve = ROOT.TEveManager.Create(True, "FI")
    viewer = eve.SpawnNewViewer("3D View", "")
    viewer.AddScene(eve.GetGlobalScene())
    viewer.AddScene(eve.GetEventScene())

    build_geometry(eve)
    build_hits(eve, args.hits, args.window, source_labels)

    eve.Redraw3D(True)
    ROOT.gSystem.ProcessEvents()

    print(f"\n[Ready] Displaying window {args.window}.")
    print("[Ready] Close the TEve window to exit.")

    ROOT.gApplication.Run(True)


if __name__ == "__main__":
    main()
