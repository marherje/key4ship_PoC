"""Parse SND_compact.xml <define> constants and <readout> BitField strings.

Usage in Gaudi job files:
    from pathlib import Path, sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "simulation" / "geometry"))
    from parse_geometry import SNDGeometry
    geo = SNDGeometry()
"""

import math
import xml.etree.ElementTree as ET
from pathlib import Path

_DEFAULT_COMPACT = Path(__file__).parent / "SND_compact.xml"


class SNDGeometry:
    """Resolved constants and BitField strings from SND_compact.xml."""

    def __init__(self, compact_file=None):
        self._file = Path(compact_file) if compact_file else _DEFAULT_COMPACT
        self._constants: dict[str, float] = {}
        self._bitfields: dict[str, str] = {}
        self._parse()

    # ── internal ────────────────────────────────────────────────────────────

    def _parse(self):
        tree = ET.parse(self._file)
        root = tree.getroot()

        define = root.find("define")
        if define is None:
            raise ValueError(f"No <define> block in {self._file}")
        raw = {
            c.get("name"): c.get("value", "0")
            for c in define.findall("constant")
            if c.get("name")
        }
        self._constants = self._resolve(raw)

        readouts = root.find("readouts")
        if readouts is not None:
            for ro in readouts.findall("readout"):
                name = ro.get("name")
                id_el = ro.find("id")
                if name and id_el is not None and id_el.text:
                    self._bitfields[name] = id_el.text.strip()

    def _resolve(self, raw: dict) -> dict:
        # Unit aliases (everything stored in mm / tesla / degrees as-is)
        ctx: dict = {
            "mm": 1.0,
            "cm": 10.0,
            "m": 1000.0,
            "tesla": 1.0,
            "deg": 1.0,
            "pi": math.pi,
        }
        todo = dict(raw)
        for _ in range(len(todo) + 1):
            if not todo:
                break
            progress = False
            for name, expr in list(todo.items()):
                try:
                    val = float(eval(expr, {"__builtins__": {}}, ctx))  # noqa: S307
                    ctx[name] = val
                    del todo[name]
                    progress = True
                except Exception:
                    pass
            if not progress:
                break
        _units = {"mm", "cm", "m", "tesla", "deg", "pi"}
        return {k: v for k, v in ctx.items() if k not in _units}

    # ── raw access ──────────────────────────────────────────────────────────

    @property
    def constants(self) -> dict:
        return self._constants

    @property
    def bitfields(self) -> dict:
        return self._bitfields

    # ── SiTarget ────────────────────────────────────────────────────────────

    @property
    def sitarget_strip_pitch(self) -> float:
        return self._constants["SiTarget_strip_pitch"]

    @property
    def sitarget_sensor_width(self) -> float:
        return self._constants["SiTarget_sensor_width"]

    @property
    def sitarget_sensor_height(self) -> float:
        return self._constants["SiTarget_sensor_height"]

    @property
    def sitarget_sensor_ncols(self) -> int:
        return int(self._constants["SiTarget_ncols"])

    @property
    def sitarget_sensor_nrows(self) -> int:
        return int(self._constants["SiTarget_nrows"])

    @property
    def sitarget_sensor_gap(self) -> float:
        return self._constants["SiTarget_sensor_gap"]

    # ── ECAL (SiPad) ────────────────────────────────────────────────────────

    @property
    def ecal_cell_size_x(self) -> float:
        return self._constants["Ecal_CellSizeX"]

    @property
    def ecal_cell_size_y(self) -> float:
        return self._constants["Ecal_CellSizeY"]

    @property
    def ecal_dim_z(self) -> float:
        return self._constants["Ecal_dim_z"]

    # ── MTC ─────────────────────────────────────────────────────────────────

    @property
    def mtc_bfield_y(self) -> float:
        return self._constants["MTC_BFieldY"]

    @property
    def mtc_outer_fe_thick(self) -> float:
        return self._constants["MTC_FeThickness"]

    @property
    def mtc_layer_thick(self) -> float:
        return self._constants["MTC_layer_thickness"]

    @property
    def mtc_n_layers(self) -> int:
        return int(self._constants["MTC_NLayers"])

    @property
    def mtc_inter_gap(self) -> float:
        return self._constants["MTC_inter_station_gap"]

    @property
    def mtc_fiber_angle_deg(self) -> float:
        return self._constants["MTC_fiber_angle_deg"]

    @property
    def mtc_scifi_channel_size(self) -> float:
        return self._constants["MTC_SciFiChannelSize"]

    @property
    def mtc_station_env_half_heights(self) -> list:
        return [
            self._constants["MTC40_env_height"] / 2.0,
            self._constants["MTC50_env_height"] / 2.0,
            self._constants["MTC60_env_height"] / 2.0,
        ]

    @property
    def mtc_station_z_centers(self) -> list:
        return [
            self._constants["MTC40_z_position"],
            self._constants["MTC50_z_position"],
            self._constants["MTC60_z_position"],
        ]

    @property
    def mtc_scifi_mip_value(self) -> float:
        return 2e-3 * self._constants["MTC_SciFiThickness"]
