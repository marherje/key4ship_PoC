#!/usr/bin/env python3
"""
SND Event Display using ROOT TEve.
Geometry is fully hardcoded (no DD4hep compact XML or plugin libraries needed).
Hits are read from an RNTuple file and displayed as TEveGeoShape boxes.
Reconstructed tracks are optionally read from a podio/edm4hep tracks.root
file and drawn as polylines.

Usage:
    # Hits only:
    python event_display_eve.py --hits ../gaudi_jobs/two_muons_pipeline/ShipHits.root --window 0 \
        --source-labels "0:mixed" "1:Sig1" "2:Sig2"

    # Hits + reconstructed tracks (tracks read from same ShipHits.root):
    python event_display_eve.py --hits ../gaudi_jobs/two_muons_pipeline/ShipHits.root --window 0 \
        --source-labels "0:mixed" "1:Sig1" "2:Sig2"
"""

import argparse
import ROOT

# ---------------------------------------------------------------------------
# 1. CONSTANTS  (sizes in cm — TGeo/TEve use cm)
# ---------------------------------------------------------------------------

# Copy here the output from inspect_SND.py
# SiPad stacks along Z  → SIPAD_Z holds Z positions of sensitive planes
# SiTarget stacks aloong Z with strips in X direction → SITARGET_X_Z hold Z positions of sensitive planes
SIPAD_Z = [
    -13.95, -12.50, -11.05, -9.60, -8.15, -6.70, -5.25,
    -3.80, -2.35, -0.90, 0.55, 2.00, 3.45, 4.90, 6.35, 7.80,
    9.25, 10.70, 12.15, 13.60,
]
SITARGET_X_Z = [
    -37.04, -35.94, -34.84, -33.74, -32.64, -31.54, -30.44,
    -29.34, -28.24, -27.14, -26.04, -24.94, -23.84, -22.74,
    -21.64, -20.54, -19.44, -18.34, -17.24, -16.14,
]
SITARGET_Y_Z = [
    -36.46, -35.36, -34.26, -33.16, -32.06, -30.96, -29.86,
    -28.76, -27.66, -26.56, -25.46, -24.36, -23.26, -22.16,
    -21.06, -19.96, -18.86, -17.76, -16.66, -15.56,
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

    SiTarget has a bf_plane field (0=Y strip, 1=Z strip).
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
def nearest_plane(hit_cm, plane_list):
    """Return the position of the closest sensitive plane to hit_cm."""
    return min(plane_list, key=lambda p: abs(p - hit_cm))


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

    # SiPad planes — stacks along X, each plane is a YZ face
    sipad_list = ROOT.TEveElementList("SiPad planes")
    for i, z in enumerate(SIPAD_Z):
        gs = make_box(f"SiPad_{i}", 17.6, 17.6, 0.0325, 
                      0.0, 0.0, z,
                      ROOT.kGreen - 6, transparency=90)
        sipad_list.AddElement(gs)
    geo_list.AddElement(sipad_list)

    # SiTarget planes — stacks along Z, each plane is an XY face
    sitgt_x_list = ROOT.TEveElementList("SiTarget X planes")
    for i, z in enumerate(SITARGET_X_Z):
        gs = make_box(f"SiTarget_X_{i}", 20.0, 20.0, 0.015, 
                      0.0, 0.0, z,
                      ROOT.kAzure + 7, transparency=90)
        sitgt_x_list.AddElement(gs)
    geo_list.AddElement(sitgt_x_list)

    sitgt_y_list = ROOT.TEveElementList("SiTarget Y planes")
    for i, z in enumerate(SITARGET_Y_Z):
        gs = make_box(f"SiTarget_Y_{i}", 20.0, 20.0, 0.015, 
                      0.0, 0.0, z,
                      ROOT.kAzure + 9, transparency=90)
        sitgt_y_list.AddElement(gs)
    geo_list.AddElement(sitgt_y_list)

    eve.GetGlobalScene().AddElement(geo_list)
    print("[Geometry] Added SiPad (10) + SiTarget X (10) + SiTarget Y (10) planes")


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
                    # SiPad: thin in X (stacking axis), hit at precise Y,Z
                    snap_z = nearest_plane(zc, SIPAD_Z)
                    gs = make_box(f"SiPad_hit_{i}",
                                  0.275, 0.275, 0.0325, 
                                  xc, yc, snap_z,
                                  color, transparency=0)

                else:  # SiTarget: stacks in Z, strips in XY plane
                    plane = planes[i]
                    if plane == 0:   # StripX: mide X, Y indefinido
                        snap_z = nearest_plane(zc, SITARGET_X_Z)
                        gs = make_box(f"SiTgt_X_hit_{i}",
                            0.003775, 20.0, 0.015,
                            xc, 0.0, snap_z,        # Y=0 (indefinido)
                            color, transparency=0)
                    else:            # StripY: mide Y, X indefinido
                        snap_z = nearest_plane(zc, SITARGET_Y_Z)
                        gs = make_box(f"SiTgt_Y_hit_{i}",
                            20.0, 0.003775, 0.015,
                            0.0, yc, snap_z,        # X=0 (indefinido)
                            color, transparency=0)

                grp.AddElement(gs)

            det_list.AddElement(grp)
            print(f"[Hits]   {det} [{label}]: {len(indices)} hits")

        hits_root.AddElement(det_list)

    eve.GetEventScene().AddElement(hits_root)


# ---------------------------------------------------------------------------
# 8. TRACK DATA READER (call BEFORE TEveManager is created)
# ---------------------------------------------------------------------------
def read_track_points(hits_file, window_id):
    """
    Read reconstructed tracks from the ACTSTracks RNTuple inside ShipHits.root.
    Uses seed_x, seed_y to assign hits to each track.
    Must be called BEFORE TEveManager.Create().
    Returns list of track data dicts.
    """
    import os
    if not os.path.exists(hits_file):
        print(f"[Tracks] ERROR: hits file not found: '{hits_file}'")
        raise SystemExit(1)

    import ROOT as _R
    _R.gSystem.Load("libROOTNTuple")

    # Read track seeds from ACTSTracks RNTuple
    _R.gROOT.ProcessLine(f"""
    {{
      SNDEve::gX.clear(); SNDEve::gY.clear();
      SNDEve::gZ.clear(); SNDEve::gE.clear();
      SNDEve::gSrc.clear(); SNDEve::gPlane.clear();
      try {{
        auto r = ROOT::RNTupleReader::Open("ACTSTracks", "{hits_file}");
        auto vwin   = r->GetView<int>("window_id");
        auto vtid   = r->GetView<int>("track_id");
        auto vchi2  = r->GetView<float>("chi2");
        auto vndf   = r->GetView<int>("ndf");
        auto vseedx = r->GetView<float>("seed_x");
        auto vseedy = r->GetView<float>("seed_y");
        for (auto i : r->GetEntryRange()) {{
          if ((int)vwin(i) != {window_id}) continue;
          SNDEve::gX.push_back(vseedx(i));
          SNDEve::gY.push_back(vseedy(i));
          SNDEve::gZ.push_back(vchi2(i));
          SNDEve::gE.push_back((float)vndf(i));
          SNDEve::gSrc.push_back(vtid(i));
          SNDEve::gPlane.push_back(0);
        }}
      }} catch (...) {{}}
    }}
    """)

    import ROOT
    n = ROOT.SNDEve.gX.size()
    if n == 0:
        print(f"[Tracks] No tracks found in window {window_id}.")
        return []

    seeds = []
    for i in range(n):
        seeds.append({
            'track_id': ROOT.SNDEve.gSrc[i],
            'seed_x':   ROOT.SNDEve.gX[i],
            'seed_y':   ROOT.SNDEve.gY[i],
            'chi2':     ROOT.SNDEve.gZ[i],
            'ndf':      int(ROOT.SNDEve.gE[i]),
        })
    print(f"[Tracks] Window {window_id}: {n} track(s) found.")

    # Read SiTarget hits
    ROOT.gROOT.ProcessLine(f"""
    {{
      SNDEve::gX.clear(); SNDEve::gY.clear(); SNDEve::gZ.clear();
      SNDEve::gE.clear(); SNDEve::gSrc.clear(); SNDEve::gPlane.clear();
      auto reader = ROOT::RNTupleReader::Open("SiTarget", "{hits_file}");
      auto vwin   = reader->GetView<int>("window_id");
      auto vx     = reader->GetView<float>("x");
      auto vy     = reader->GetView<float>("y");
      auto vz     = reader->GetView<float>("z");
      auto vplane = reader->GetView<int>("bf_plane");
      for (auto i : reader->GetEntryRange()) {{
        if ((int)vwin(i) != {window_id}) continue;
        SNDEve::gX.push_back(vx(i));
        SNDEve::gY.push_back(vy(i));
        SNDEve::gZ.push_back(vz(i));
        SNDEve::gE.push_back(0.f);
        SNDEve::gSrc.push_back(0);
        SNDEve::gPlane.push_back(vplane(i));
      }}
    }}
    """)
    mm2cm = 0.1
    n = ROOT.SNDEve.gX.size()
    sitarget_hits = []
    for i in range(n):
        sitarget_hits.append({
            'x': ROOT.SNDEve.gX[i], 'y': ROOT.SNDEve.gY[i],
            'z': ROOT.SNDEve.gZ[i], 'plane': ROOT.SNDEve.gPlane[i]
        })

    # Read SiPad hits
    ROOT.gROOT.ProcessLine(f"""
    {{
      SNDEve::gX.clear(); SNDEve::gY.clear(); SNDEve::gZ.clear();
      SNDEve::gE.clear(); SNDEve::gSrc.clear(); SNDEve::gPlane.clear();
      auto reader = ROOT::RNTupleReader::Open("SiPad", "{hits_file}");
      auto vwin = reader->GetView<int>("window_id");
      auto vx   = reader->GetView<float>("x");
      auto vy   = reader->GetView<float>("y");
      auto vz   = reader->GetView<float>("z");
      for (auto i : reader->GetEntryRange()) {{
        if ((int)vwin(i) != {window_id}) continue;
        SNDEve::gX.push_back(vx(i));
        SNDEve::gY.push_back(vy(i));
        SNDEve::gZ.push_back(vz(i));
        SNDEve::gE.push_back(0.f);
        SNDEve::gSrc.push_back(0);
        SNDEve::gPlane.push_back(-1);
      }}
    }}
    """)
    n = ROOT.SNDEve.gX.size()
    sipad_hits = []
    for i in range(n):
        sipad_hits.append({
            'x': ROOT.SNDEve.gX[i], 'y': ROOT.SNDEve.gY[i],
            'z': ROOT.SNDEve.gZ[i]
        })

    # Build per-track trajectories using seed-based hit selection
    track_data = []
    for seed in seeds:
        sx, sy = seed['seed_x'], seed['seed_y']

        # Group SiTarget hits by station (beam Z rounded to 1mm)
        stations = {}
        for h in sitarget_hits:
            key = round(h['z'], 0)
            if key not in stations:
                stations[key] = {'x_hits': [], 'y_hits': [], 'z': h['z']}
            if h['plane'] == 0:
                stations[key]['x_hits'].append(h['x'])
            else:
                stations[key]['y_hits'].append(h['x'])  # Y stored in pos.x

        points = []
        for key in sorted(stations):
            data  = stations[key]
            z     = data['z']
 
            points.append((sx * mm2cm, sy * mm2cm, z * mm2cm, z))

        # SiPad: best hit per layer closest to seed
        sipad_by_z = {}
        for h in sipad_hits:
            key  = round(h['z'], 0)
            dist = (h['x']-sx)**2 + (h['y']-sy)**2
            if key not in sipad_by_z or dist < sipad_by_z[key][2]:
                sipad_by_z[key] = (h['x'], h['y'], dist, h['z'])
        for key in sorted(sipad_by_z):
            x, y, _, z = sipad_by_z[key]
            points.append((sx*mm2cm, sy*mm2cm, z*mm2cm, z))

        points.sort(key=lambda p: p[3])

        if len(points) >= 2:
            track_data.append({
                'points': points,
                'chi2':   seed['chi2'],
                'ndf':    seed['ndf'],
                'seed_x': sx,
                'seed_y': sy,
            })
            print(f"[Tracks]   Track {seed['track_id']}: "
                  f"chi2={seed['chi2']:.1f} seed=({sx:.1f},{sy:.1f})mm "
                  f"{len(points)} points")

    return track_data


# ---------------------------------------------------------------------------
# 9. TRACK DRAWER (call AFTER TEveManager is created)
# ---------------------------------------------------------------------------
def draw_tracks(eve, track_data):
    """
    Draw pre-loaded track trajectories as TEve line sets.
    track_data is the output of read_track_points().
    Must be called AFTER TEveManager.Create() but BEFORE ProcessEvents().
    """
    if not track_data:
        print("[Tracks] No track data to draw.")
        return

    track_display_colors = [
        ROOT.kRed,
        ROOT.kOrange + 7,
        ROOT.kMagenta,
        ROOT.kYellow + 1,
    ]

    tracks_list = ROOT.TEveElementList("Reconstructed Tracks")

    for iTrack, td in enumerate(track_data):
        color  = track_display_colors[iTrack % len(track_display_colors)]
        points = td['points']
        chi2   = td['chi2']
        ndf    = td['ndf']

        label = f"Track {iTrack}  chi2={chi2:.1f}  ndf={ndf}"
        line  = ROOT.TEveStraightLineSet(label)
        line.SetLineColor(color)
        line.SetLineWidth(3)
        line.SetMarkerColor(color)
        line.SetMarkerSize(0.8)
        line.SetMarkerStyle(4)

        for j in range(len(points) - 1):
            x0, y0, z0, _ = points[j]
            x1, y1, z1, _ = points[j + 1]
            line.AddLine(x0, y0, z0, x1, y1, z1)

        for pt in points:
            line.AddMarker(pt[0], pt[1], pt[2])

        tracks_list.AddElement(line)

    eve.GetEventScene().AddElement(tracks_list)
    print(f"[Tracks] Added {len(track_data)} track line(s) to event scene.")


# ---------------------------------------------------------------------------
# 10. MAIN
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="SND Event Display (TEve)")
    parser.add_argument("--hits",   default="../ShipHits.root")
    parser.add_argument("--window", type=int, default=0,
                        help="Window/event index for both hits and tracks.")
    parser.add_argument("--source-labels", nargs="*",
                        default=["0:mixed", "1:signal", "2:background"],
                        metavar="ID:LABEL",
                        help="source_id labels e.g. '0:mixed' '1:Sig1'")
    args = parser.parse_args()

    source_labels = {}
    for item in args.source_labels:
        sid, lbl = item.split(":", 1)
        source_labels[int(sid)] = lbl

    # ---- Read track data BEFORE loading ROOT GL libraries -------------------
    # Tracks are read from the same --hits file (ACTSTracks RNTuple inside ShipHits.root)
    track_data = []
    if True:  # always attempt — tracks may or may not be in the hits file
        track_data = read_track_points(args.hits, args.window)

    # ---- Load ROOT and create display ---------------------------------------
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gSystem.Load("libEve")
    ROOT.gSystem.Load("libEG")
    ROOT.gSystem.Load("libROOTNTuple")

    # Populate source_colors now that ROOT is loaded
    for sid, (r, g, b) in SOURCE_RGB.items():
        idx = ROOT.TColor.GetFreeColorIndex()
        ROOT.TColor(idx, r, g, b)
        source_colors[sid] = idx

    eve = ROOT.TEveManager.Create(True, "FI")
    viewer = eve.SpawnNewViewer("3D View", "")
    viewer.AddScene(eve.GetGlobalScene())
    viewer.AddScene(eve.GetEventScene())

    build_geometry(eve)
    build_hits(eve, args.hits, args.window, source_labels)

    # Draw pre-loaded tracks (no podio loading here — already done above)
    if track_data:
        draw_tracks(eve, track_data)

    eve.Redraw3D(True)
    ROOT.gSystem.ProcessEvents()

    print(f"\n[Ready] Displaying window {args.window}.")
    print("[Ready] Close the TEve window to exit.")

    ROOT.gApplication.Run(True)


if __name__ == "__main__":
    main()
