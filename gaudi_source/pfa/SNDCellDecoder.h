#pragma once

#include "DDSegmentation/BitFieldCoder.h"

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// SNDCellDecoder — one place that turns a cellID into the address SND-PFA
// works in.
//
// The three subdetectors describe themselves differently: SiTarget strips are
// (layer, plane) with plane 0 = X and 1 = Y, SiPad pads are (layer, x, y) with
// no plane at all, and the MTC restarts its layer count in every station and
// uses plane 0/1 for the stereo SciFi views and 2 for the scintillator tiles.
// Every PFA stage needs the same three answers regardless — which detector,
// which layer counted along the beam, and which view — so they are computed
// here once instead of in each algorithm.
//
// `layerOrdinal` is contiguous within a detector (0 = most upstream), so
// subtracting two of them is a layer count. That is the same convention as
// ISNDGeoSvc::layerOrdinalOf, which the trackers use.
// ---------------------------------------------------------------------------

namespace snd {

enum class Det : int { SiTarget = 0, SiPad = 1, MTC = 2, Unknown = -1 };

// View of a sensitive plane: which coordinate it measures.
enum class View : int {
  X = 0,   // SiTarget X strips, MTC SciFi U
  Y = 1,   // SiTarget Y strips, MTC SciFi V
  TwoD = 2 // SiPad pads, MTC scintillator tiles
};

struct CellAddress {
  Det   det = Det::Unknown;
  int   station = -1;      // MTC only; -1 elsewhere
  int   layer = -1;        // raw layer, restarts per station in the MTC
  int   plane = -1;        // raw plane field, -1 when the readout has none
  View  view = View::TwoD;
  int   layerOrdinal = -1; // contiguous along the beam within the detector
  bool  valid = false;
};

class SNDCellDecoder {
public:
  // `readoutByCollection` maps a hit collection name to the id descriptor of
  // its readout (the same strings the steering files take from
  // parse_geometry.SNDGeometry.bitfields). `mtcLayersPerStation` folds the
  // MTC station into the layer ordinal; it comes from the DD4hep constant
  // MTC_NLayers, never from a literal.
  SNDCellDecoder(const std::map<std::string, std::string>& readoutByCollection,
                 const std::map<std::string, Det>& detByCollection,
                 int mtcLayersPerStation)
      : m_mtcLayersPerStation(mtcLayersPerStation) {
    for (const auto& [coll, desc] : readoutByCollection) {
      Readout r;
      r.coder =
          std::make_shared<dd4hep::DDSegmentation::BitFieldCoder>(desc);
      auto it = detByCollection.find(coll);
      r.det = (it != detByCollection.end()) ? it->second : Det::Unknown;
      // Resolve the field positions once: the descriptors differ between
      // detectors (no `plane` in the SiPad, no `station` outside the MTC) and
      // BitFieldCoder::get() throws on an unknown name, which is far too
      // expensive to do per hit.
      r.iLayer   = indexOf(*r.coder, "layer");
      r.iPlane   = indexOf(*r.coder, "plane");
      r.iStation = indexOf(*r.coder, "station");
      r.iStrip   = indexOf(*r.coder, "strip");
      r.iX       = indexOf(*r.coder, "x");
      r.iY       = indexOf(*r.coder, "y");
      m_readouts.emplace(coll, std::move(r));
    }
  }

  bool knows(const std::string& collection) const {
    return m_readouts.find(collection) != m_readouts.end();
  }

  CellAddress decode(const std::string& collection, uint64_t cellID) const {
    CellAddress a;
    auto itR = m_readouts.find(collection);
    if (itR == m_readouts.end()) return a;
    const auto& r = itR->second;
    a.det = r.det;

    a.layer = value(r, r.iLayer, cellID);
    a.plane = value(r, r.iPlane, cellID);

    switch (a.det) {
      case Det::SiTarget:
        // plane 0 = CartesianStripX, 1 = CartesianStripY.
        a.view = (a.plane == 1) ? View::Y : View::X;
        a.layerOrdinal = a.layer;
        break;
      case Det::SiPad:
        // No plane field: one 2D pad plane per layer.
        a.view = View::TwoD;
        a.layerOrdinal = a.layer;
        break;
      case Det::MTC:
        a.station = value(r, r.iStation, cellID);
        // plane 0/1 = stereo SciFi U/V, 2 = scintillator tiles.
        a.view = (a.plane == 0) ? View::X
               : (a.plane == 1) ? View::Y
                                : View::TwoD;
        a.layerOrdinal = (a.station >= 0 && a.layer >= 0)
            ? a.station * m_mtcLayersPerStation + a.layer
            : -1;
        break;
      default:
        return a;
    }
    a.valid = (a.layerOrdinal >= 0);
    return a;
  }

  // Channel index within a plane: the strip number for 1D readouts, and a
  // (x,y) pair folded into one integer for 2D ones. Neighbouring channels
  // differ by one in the 1D case, which is what the layer clusterer needs.
  //
  // Which of the two applies is decided by the VIEW, not by which fields the
  // descriptor happens to define: the MTC readout declares `strip`, `x` and
  // `y` for all of its planes, and the unused ones simply read back as 0. An
  // earlier version returned `strip` whenever the field existed, which mapped
  // every scintillator tile to channel 0.
  int channel(const std::string& collection, uint64_t cellID,
              int* x = nullptr, int* y = nullptr) const {
    auto itR = m_readouts.find(collection);
    if (itR == m_readouts.end()) return -1;
    const auto& r = itR->second;
    const auto a = decode(collection, cellID);
    if (a.view != View::TwoD) {
      const int strip = value(r, r.iStrip, cellID);
      if (x) *x = strip;
      if (y) *y = 0;
      return strip;
    }
    const int cx = value(r, r.iX, cellID);
    const int cy = value(r, r.iY, cellID);
    if (x) *x = cx;
    if (y) *y = cy;
    return (cx >= 0 && cy >= 0) ? (cx * 4096 + cy) : -1;
  }

private:
  struct Readout {
    std::shared_ptr<dd4hep::DDSegmentation::BitFieldCoder> coder;
    Det det = Det::Unknown;
    int iLayer = -1, iPlane = -1, iStation = -1;
    int iStrip = -1, iX = -1, iY = -1;
  };

  // Position of a field in the descriptor, or -1 when it does not define one.
  static int indexOf(const dd4hep::DDSegmentation::BitFieldCoder& coder,
                     const std::string& name) {
    const auto& fields = coder.fields();
    for (std::size_t i = 0; i < fields.size(); ++i)
      if (fields[i].name() == name) return static_cast<int>(i);
    return -1;
  }

  static int value(const Readout& r, int fieldIndex, uint64_t cellID) {
    if (fieldIndex < 0) return -1;
    return static_cast<int>(
        r.coder->get(cellID, static_cast<std::size_t>(fieldIndex)));
  }

  std::map<std::string, Readout> m_readouts;
  int m_mtcLayersPerStation = 1;
};

}  // namespace snd
