#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

// ---------------------------------------------------------------------------
// CLUE3D on the SND's layer points.
//
// CLUE (Rovere et al., Front. Big Data 3 (2020) 591315) clusters by density:
// every point gets a local energy density rho, and the distance delta to the
// nearest point of higher density. Points with high rho and large delta are
// SEEDS (shower cores); points with low rho and large delta are OUTLIERS
// (isolated noise); everything else FOLLOWS its nearest-higher, and the chain
// of followers ends at a seed, which names the cluster. CLUE3D is the version
// CMS TICL runs on HGCAL layer clusters: the density and the nearest-higher
// search look a few layers up and down the beam, so a cluster is a 3D object
// grown along the shower rather than a per-layer blob.
//
// Two departures from TICL, both because of the SND's mixed granularity:
//   * the critical distance and the critical density are PER DETECTOR
//     (SiTarget strips of 75 um under tungsten, SiPad 5.5 mm pads, MTC iron
//     read by 1 mm SciFi and 50 mm tiles): one value cannot suit all three;
//   * the distance is transverse only, in mm, with the layer window doing the
//     longitudinal job -- the layers of the three detectors have different
//     pitches, so a mixed (mm, layer) metric would mean different things in
//     each.
//
// It is O(N^2) without the tiles TICL uses: an SND event has a few hundred
// layer points, and the tiles only pay off at HGCAL pileup occupancies.
// ---------------------------------------------------------------------------

namespace snd::clue {

struct Point {
  double x = 0.0, y = 0.0;   // transverse position [mm]
  int    layer = 0;          // index in the beam-ordered layer sequence
  int    det = 0;            // 0 = SiTarget, 1 = SiPad, 2 = MTC
  double energy = 0.0;       // [GeV]
};

struct Params {
  std::array<double, 3> dc{10.0, 15.0, 40.0};    // critical distance [mm]
  std::array<double, 3> rhoc{0.02, 0.04, 0.2};   // seed density [GeV]
  int    densityLayers = 2;        // +- layers summed into rho
  int    nearestHigherLayers = 3;  // +- layers searched for the nearest-higher
  double outlierDeltaFactor = 2.0; // delta > factor x dc and rho < rhoc: outlier
  double neighbourWeight = 0.5;    // weight of the other points in rho (TICL)
};

struct Result {
  std::vector<double> rho, delta;
  std::vector<int> nearestHigher;  // -1 if none in the window
  std::vector<int> cluster;        // -1 for outliers / unassigned
  std::vector<bool> isSeed;
  int nClusters = 0;
};

inline Result run(const std::vector<Point>& pts, const Params& par) {
  const int n = static_cast<int>(pts.size());
  Result r;
  r.rho.assign(n, 0.0);
  r.delta.assign(n, std::numeric_limits<double>::infinity());
  r.nearestHigher.assign(n, -1);
  r.cluster.assign(n, -1);
  r.isSeed.assign(n, false);

  auto dist = [&](int i, int j) {
    return std::hypot(pts[i].x - pts[j].x, pts[i].y - pts[j].y);
  };

  // 1. local density
  for (int i = 0; i < n; ++i) {
    const double dc = par.dc[pts[i].det];
    for (int j = 0; j < n; ++j) {
      if (std::abs(pts[i].layer - pts[j].layer) > par.densityLayers) continue;
      if (i != j && dist(i, j) >= dc) continue;
      r.rho[i] += (i == j ? 1.0 : par.neighbourWeight) * pts[j].energy;
    }
  }

  // 2. nearest point of higher density, within the layer window and a
  //    transverse window of outlierDeltaFactor x dc. Equal densities are
  //    ordered by index so the relation has no cycles.
  auto higher = [&](int j, int i) {
    return r.rho[j] > r.rho[i] || (r.rho[j] == r.rho[i] && j < i);
  };
  for (int i = 0; i < n; ++i) {
    const double window = par.outlierDeltaFactor * par.dc[pts[i].det];
    double best = std::numeric_limits<double>::infinity();
    int bestLayerGap = std::numeric_limits<int>::max();
    for (int j = 0; j < n; ++j) {
      if (j == i || !higher(j, i)) continue;
      const int gap = std::abs(pts[i].layer - pts[j].layer);
      if (gap > par.nearestHigherLayers) continue;
      const double d = dist(i, j);
      if (d > window) continue;
      if (d < best || (d == best && gap < bestLayerGap)) {
        best = d;
        bestLayerGap = gap;
        r.nearestHigher[i] = j;
      }
    }
    r.delta[i] = best;
  }

  // 3. seeds, outliers, followers
  std::vector<std::vector<int>> followers(n);
  std::vector<int> stack;
  for (int i = 0; i < n; ++i) {
    const double dc = par.dc[pts[i].det];
    const double rhoc = par.rhoc[pts[i].det];
    const bool seed = r.rho[i] >= rhoc && r.delta[i] > dc;
    const bool outlier = !seed && r.rho[i] < rhoc &&
                         r.delta[i] > par.outlierDeltaFactor * dc;
    if (seed) {
      r.isSeed[i] = true;
      r.cluster[i] = r.nClusters++;
      stack.push_back(i);
    } else if (!outlier && r.nearestHigher[i] >= 0) {
      followers[r.nearestHigher[i]].push_back(i);
    }
  }
  // 4. hand the seed's label down its chain of followers
  while (!stack.empty()) {
    const int i = stack.back();
    stack.pop_back();
    for (int f : followers[i]) {
      r.cluster[f] = r.cluster[i];
      stack.push_back(f);
    }
  }
  return r;
}

}  // namespace snd::clue
