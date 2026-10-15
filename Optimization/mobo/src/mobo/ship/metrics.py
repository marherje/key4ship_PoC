"""The metrics of a trial: what is measured, and what is computed.

Split deliberately in two, because the two halves want to run in different
places:

* `analytic_metrics` is a pure function of the resolved geometry constants —
  cost, silicon area, tungsten mass, channel count. No simulation, no ROOT, no
  worker node: the driver can evaluate it in microseconds, which also makes it
  testable against numbers worked out by hand.
* `hit_metrics` reads the RNTuples the analysis chain produced, so it needs ROOT
  and runs on the worker as part of the payload.

Every run records *everything* it can, not just the objectives. Which subset is
optimized is a config decision (`ObjectiveSpec`), so adding an objective later
must never mean re-running the campaign.

Stdlib-only at import time on purpose: the payload imports this module by path
with the bare key4hep python, which has no torch and no mobo package installed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

# Density of tungsten. The absorber is not pure W in reality, but the whole
# point of a proxy cost is to rank designs, and every design pays the same
# alloy correction.
RHO_W_G_CM3 = 19.3

# RNTuple names written by job5_rntuple.py, and the metric each becomes.
HIT_CONTAINERS = {
    "SiPad": "nhits_sipad",
    "SiTarget": "nhits_sitarget",
    "MTCSciFi": "nhits_mtcscifi",
    "MTCScint": "nhits_mtcscint",
}


# Track RNTuples written by job5_rntuple_seq.py, and the metric suffix each
# becomes. GlobalTracks first: it is the one `eff_track` is defined on.
TRACK_CONTAINERS = {
    "GlobalTracks": "global",
    "SiTargetTracks": "sitarget",
    "SiPadTracks": "sipad",
    "MTCTracks": "mtc",
}


@dataclass(frozen=True)
class CostModel:
    """Coefficients of the proxy cost, in k€ per unit.

    Order-of-magnitude figures for detector-grade silicon and for tungsten;
    what matters is that they are explicit and configurable, since the trade-off
    the optimizer finds is a direct consequence of their ratio.
    """

    si_per_m2: float = 30.0  # k€ / m² of silicon sensor
    w_per_kg: float = 0.1  # k€ / kg of tungsten
    rho_w_g_cm3: float = RHO_W_G_CM3

    @classmethod
    def from_config(cls, cfg: Any) -> CostModel:
        if cfg is None:
            return cls()
        fields = set(cls.__dataclass_fields__)
        return cls(**{k: float(v) for k, v in dict(cfg).items() if k in fields})


# ── the analytic half ────────────────────────────────────────────────────────


def silicon_areas_m2(c: dict[str, float]) -> dict[str, float]:
    """Silicon *sensor* area per detector, in m².

    Sensor area, not envelope area: what is paid for is the wafers, and the
    inactive border between them is not silicon. In SiTarget every layer carries
    two sensor planes (one measuring x, one measuring y), in SiPad one.
    """
    sipad_plane = (c["SiPad_NASUsX"] * c["SiPad_NWafersX"] * c["SiPad_WaferSizeX"]) * (
        c["SiPad_NASUsY"] * c["SiPad_NWafersY"] * c["SiPad_WaferSizeY"]
    )
    sitarget_plane = (c["SiTarget_ncols"] * c["SiTarget_sensor_width"]) * (
        c["SiTarget_nrows"] * c["SiTarget_sensor_height"]
    )
    return {
        "sipad": c["SiPad_NLayers"] * sipad_plane / 1e6,
        "sitarget": 2.0 * c["SiTarget_NLayers"] * sitarget_plane / 1e6,
    }


def tungsten_masses_kg(
    c: dict[str, float], rho_g_cm3: float = RHO_W_G_CM3
) -> dict[str, float]:
    """Absorber mass per detector, in kg.

    The absorber spans the full transverse envelope of its detector (it is a
    plate, not a tiling of sensors), so this uses the envelope dimensions.
    """
    sipad_volume = (
        c["SiPad_NLayers"] * c["SiPad_WThickness"] * c["SiPad_dim_x"] * c["SiPad_dim_y"]
    )
    sitarget_volume = (
        c["SiTarget_NLayers"]
        * c["SiTarget_WThickness"]
        * c["SiTarget_env_width"]
        * c["SiTarget_env_height"]
    )
    # mm^3 -> cm^3 (/1000), g -> kg (/1000)
    return {
        "sipad": sipad_volume / 1e3 * rho_g_cm3 / 1e3,
        "sitarget": sitarget_volume / 1e3 * rho_g_cm3 / 1e3,
    }


def channel_counts(c: dict[str, float]) -> dict[str, float]:
    """Readout channels per detector.

    SiPad: one pad per cell of the segmentation, over the whole plane.
    SiTarget: strips per sensor, times the sensor grid, times two planes per
    layer — the x-measuring plane is striped across the sensor width and the
    y-measuring one across its height.
    """
    sipad = c["SiPad_NLayers"] * c["SiPad_NcellsX"] * c["SiPad_NcellsY"]
    pitch = c["SiTarget_strip_pitch"]
    strips_x = math.floor(c["SiTarget_sensor_width"] / pitch)
    strips_y = math.floor(c["SiTarget_sensor_height"] / pitch)
    sensors = c["SiTarget_ncols"] * c["SiTarget_nrows"]
    sitarget = c["SiTarget_NLayers"] * sensors * (strips_x + strips_y)
    return {"sipad": float(sipad), "sitarget": float(sitarget)}


def analytic_metrics(
    c: dict[str, float], cost: CostModel | None = None
) -> dict[str, float]:
    """Everything derivable from the geometry alone, for one variant."""
    cost = cost or CostModel()
    areas = silicon_areas_m2(c)
    masses = tungsten_masses_kg(c, cost.rho_w_g_cm3)
    channels = channel_counts(c)

    si_area = areas["sipad"] + areas["sitarget"]
    w_mass = masses["sipad"] + masses["sitarget"]
    return {
        "cost_proxy": cost.si_per_m2 * si_area + cost.w_per_kg * w_mass,
        "si_area_m2": si_area,
        "si_area_sipad_m2": areas["sipad"],
        "si_area_sitarget_m2": areas["sitarget"],
        "w_mass_kg": w_mass,
        "w_mass_sipad_kg": masses["sipad"],
        "w_mass_sitarget_kg": masses["sitarget"],
        "n_channels": channels["sipad"] + channels["sitarget"],
        "n_channels_sipad": channels["sipad"],
        "n_channels_sitarget": channels["sitarget"],
        "sipad_nlayers": float(c["SiPad_NLayers"]),
        "sitarget_nlayers": float(c["SiTarget_NLayers"]),
    }


# ── the measured half ────────────────────────────────────────────────────────


def count_entries(path: str, container: str) -> int:
    """Rows of one RNTuple. Same call `compute_fom.py` makes."""
    import gc

    import ROOT

    reader = ROOT.RNTupleReader.Open(container, path)
    n = int(reader.GetNEntries())
    # See distinct_events() for why the reader is released by hand.
    del reader
    gc.collect()
    return n


def distinct_events(path: str, container: str) -> set[int]:
    """The set of event indices in which `container` has at least one row.

    `window_id` is the event index every RNTuple written by EDM4HEP2RNTuple
    carries (the name is historical: one event is one time window). Read
    through GetView rather than RDataFrame, which this ROOT build does not
    expose for RNTuple.
    """
    import gc

    import ROOT

    reader = ROOT.RNTupleReader.Open(container, path)
    view = reader.GetView["std::int32_t"]("window_id")
    events = {int(view(i)) for i in reader.GetEntryRange()}
    # Released by hand, and this is NOT tidiness. Left to the interpreter, the
    # view and its reader are torn down after ROOT has already dismantled
    # itself, and RColumn's destructor segfaults: the payload then dies with
    # signal 129 AFTER computing correct metrics but BEFORE writing its DONE
    # sentinel, which the loop reports as a lost trial. Measured on ROOT 6.38.
    del view
    del reader
    gc.collect()
    return events


def track_metrics(
    path: str, nevents: int, containers: dict[str, str] | None = None
) -> dict[str, float]:
    """Track counts and the per-event tracking efficiency.

    `eff_track` -- the fraction of events holding at least one GlobalTrack --
    is the quantity this project has used throughout as "did we reconstruct
    anything at all", and the one an optimization over GENIE events wants:
    a neutrino interaction either yields a track or it does not, and counting
    tracks instead would reward hadronic showers that fragment into many short
    segments.

    The denominator is the REQUESTED event count, not the number observed.
    Every trial replays the same HepMC primaries, so the denominator is a
    constant across the campaign and cannot bias the ranking; using an observed
    count instead would silently rescale a trial whose simulation was truncated.
    `hit_metrics` records `n_events_seen` next to it, so a trial whose
    simulation was truncated is visible rather than looking like a bad geometry.
    """
    containers = containers or TRACK_CONTAINERS
    metrics: dict[str, float] = {}
    for ntuple, key in containers.items():
        try:
            metrics["ntracks_" + key] = float(count_entries(path, ntuple))
            metrics["nevents_with_" + key] = float(len(distinct_events(path, ntuple)))
        except Exception:  # noqa: BLE001 - an absent RNTuple is not an error here
            continue
    if "nevents_with_global" in metrics and nevents > 0:
        metrics["eff_track"] = metrics["nevents_with_global"] / float(nevents)
    return metrics


def hit_metrics(path: str, containers: dict[str, str] | None = None) -> dict[str, float]:
    """Hit counts from a ShipHits.root. Missing containers are simply absent.

    A missing SiPad is not silently zero: it would look like a legitimately bad
    geometry to the optimizer, when in fact the analysis chain broke. The caller
    (the payload) checks that the objective it needs is present.
    """
    containers = containers or HIT_CONTAINERS
    metrics: dict[str, float] = {}
    seen: set[int] = set()
    for ntuple, key in containers.items():
        try:
            metrics[key] = float(count_entries(path, ntuple))
            seen |= distinct_events(path, ntuple)
        except Exception:  # noqa: BLE001 - an absent RNTuple is not an error here
            continue
    # Truncation detector, and it has to be built on HITS rather than on tracks:
    # the last events of a run routinely yield no track at all, so the highest
    # event index among the track RNTuples undercounts even on a complete run.
    # Every simulated event deposits something somewhere, so the hit containers
    # do see the full range.
    if seen:
        metrics["n_events_seen"] = float(max(seen) + 1)
    return metrics
