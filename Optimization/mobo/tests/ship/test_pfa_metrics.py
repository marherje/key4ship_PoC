"""The SND-PFA metrics, against values worked out by hand.

Twenty events, written out so every metric can be counted on one's fingers:

  * 10 nu_mu CC: the muon is identified in 6 of them;
    E_nu,true = 10 GeV and E_nu,reco = 5, 6, ..., 14 GeV;
  * 5 NC: one of them has a (fake) muon candidate;
  * 5 nu_e CC: none has a candidate.

The definitions are tested without ROOT (pfa_metrics_from_columns), and the
reading of a real EventPFA RNTuple with it (marked key4hep).

    pytest tests/ship/test_pfa_metrics.py
"""

from __future__ import annotations

import gc

import pytest

from mobo.ship.metrics import pfa_metrics, pfa_metrics_from_columns

KEYS = {
    "eff_muon_id", "fake_muon_nc", "res_pmu", "bias_ehad", "res_ehad",
    "res_enu", "frac_leak_mtc", "nevents_pfa",
}


def rows() -> list[dict[str, float]]:
    out = []
    for i in range(10):                       # nu_mu CC
        out.append({
            "flavour_true": 14, "is_cc": 1,
            "n_muon_cand": 1 if i < 6 else 0, "muon_matched": 1 if i < 6 else 0,
            "E_nu_true": 10.0, "E_nu_reco": 5.0 + i,
            "E_had_true": 4.0, "E_had_reco": 3.0,
            # 4 of the 15 nu_mu CC + NC events have > 30 % of E_had in the MTC
            "E_had_reco_mtc": 2.0 if i < 4 else 0.5,
            "p_mu_true": 10.0, "p_mu_reco": 8.0 + i if i < 6 else 0.0,
        })
    for i in range(5):                        # NC
        out.append({
            "flavour_true": 14, "is_cc": 0,
            "n_muon_cand": 1 if i == 0 else 0, "muon_matched": 0,
            "E_nu_true": 10.0, "E_nu_reco": 3.0,
            "E_had_true": 4.0, "E_had_reco": 3.0, "E_had_reco_mtc": 0.5,
            "p_mu_true": float("nan"), "p_mu_reco": 0.0,
        })
    for _ in range(5):                        # nu_e CC
        out.append({
            "flavour_true": 12, "is_cc": 1, "n_muon_cand": 0, "muon_matched": 0,
            "E_nu_true": 10.0, "E_nu_reco": 6.0,
            "E_had_true": 4.0, "E_had_reco": 3.0, "E_had_reco_mtc": 0.5,
            "p_mu_true": float("nan"), "p_mu_reco": 0.0,
        })
    return out


def columns() -> dict[str, list[float]]:
    rs = rows()
    return {k: [float(r[k]) for r in rs] for k in rs[0]}


def check(m: dict[str, float]) -> None:
    assert KEYS <= set(m)
    assert m["nevents_pfa"] == 20
    # 6 of 10 nu_mu CC muons identified; 1 of 5 NC with a candidate.
    assert m["eff_muon_id"] == pytest.approx(0.6)
    assert m["fake_muon_nc"] == pytest.approx(0.2)
    assert m["fake_muon_nue"] == pytest.approx(0.0)
    # E_nu ratios 0.5 ... 1.4 on nu_mu CC: median (0.9 + 1.0) / 2 = 0.95, and
    # sigma68 = (s[int(0.84*9)] - s[int(0.16*9)]) / 2 = (1.2 - 0.6) / 2 = 0.3.
    assert m["bias_enu"] == pytest.approx(-0.05, abs=1e-6)
    assert m["res_enu"] == pytest.approx(0.3, abs=1e-6)
    # E_had 3 / 4 everywhere: bias -0.25 and no spread. nu_e CC are not in it.
    assert m["bias_ehad"] == pytest.approx(-0.25)
    assert m["res_ehad"] == pytest.approx(0.0)
    # 4 of the 15 nu_mu CC + NC events: 2 / 3 of E_had in the MTC.
    assert m["frac_leak_mtc"] == pytest.approx(4 / 15)
    # p_mu: (8 ... 13 - 10) / 10 = -0.2 ... 0.3 on the 6 found muons;
    # sigma68 = (s[int(0.84*5)] - s[int(0.16*5)]) / 2 = (0.2 - (-0.2)) / 2.
    assert m["res_pmu"] == pytest.approx(0.2, abs=1e-6)


def test_definitions_by_hand():
    check(pfa_metrics_from_columns(columns()))


def test_missing_file_is_not_an_error(tmp_path):
    assert pfa_metrics(str(tmp_path / "ShipPFA.root")) == {}


@pytest.mark.key4hep
def test_reads_an_eventpfa_rntuple(tmp_path):
    ROOT = pytest.importorskip("ROOT")
    path = str(tmp_path / "ShipPFA.root")
    ints = ("flavour_true", "is_cc", "n_muon_cand", "muon_matched")
    model = ROOT.RNTupleModel.Create()
    for name in rows()[0]:
        model.MakeField["std::int32_t" if name in ints else "float"](name)
    writer = ROOT.RNTupleWriter.Recreate(model, "EventPFA", path)
    for r in rows():
        entry = writer.CreateEntry()
        for name, value in r.items():
            entry[name] = int(value) if name in ints else value
        writer.Fill(entry)
        del entry
    # The writer flushes and closes the file when it goes.
    del writer
    gc.collect()

    check(pfa_metrics(path))
