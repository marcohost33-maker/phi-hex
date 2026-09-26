"""Gates fuer PHY047 (Neuauswertung committeter Daten mit Upsilon pro Flaeche).

Drift-Guard: die per-Site-Spalte muss die Original-Reports reproduzieren
(gleiche Daten, gleiche Methoden). Nur dann ist die per-Flaeche-Spalte eine
reine Konsequenz der Normierung und kein Methoden-Artefakt.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import conftest

phy047 = conftest._load(
    "phy047_reanalysis", "260926 PHY047 per-area reanalysis committed v01.py")

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def rep():
    return phy047.run()


def test_per_site_column_reproduces_committed_phy042_pairs(rep):
    s = rep["sets"]["honeycomb_phy042_wl"]
    for k, v in s["committed_pairs_per_site"].items():
        assert s["per_site"]["pairs"][k] == pytest.approx(v, abs=1e-12), k


def test_per_site_column_reproduces_committed_table_pairs(rep):
    """PHY032 (12,24)=0.6014, (24,48)=0.5917; PHY033 (12,24)=0.8437,
    (12,36)=0.8410, (24,36)=0.8377 (committete Report-Werte, 4 Dezimalen)."""
    hc = rep["sets"]["honeycomb_phy032"]["per_site"]["pairs"]
    assert round(hc["12_24"], 4) == 0.6014
    assert round(hc["24_48"], 4) == 0.5917
    kg = rep["sets"]["kagome_phy033"]["per_site"]["pairs"]
    assert (round(kg["12_24"], 4), round(kg["12_36"], 4),
            round(kg["24_36"], 4)) == (0.8437, 0.8410, 0.8377)


def test_triangular_per_site_wm_matches_committed_value(rep):
    wm = rep["sets"]["triangular_phy030v02"]["per_site"]["wm_free_c"]
    assert wm["T"] == pytest.approx(1.4007, abs=2e-3)


def test_per_area_is_pure_rescaling(rep):
    """Die per-Flaeche-Auswertung unterscheidet sich NUR durch 1/a_s."""
    t, m, e = phy047.parse_table(phy047.SOURCES["kagome_phy033"])
    a_s = phy047.A_S["kagome"]
    direct = phy047.analyse_set(t, {L: v / a_s for L, v in m.items()},
                                {L: v / a_s for L, v in e.items()}, 1.0)
    ref = rep["sets"]["kagome_phy033"]["per_area"]["pairs"]
    assert set(direct["per_site"]["pairs"]) == set(ref)
    for k, v in direct["per_site"]["pairs"].items():   # (m/a)*1 vs m*(1/a)
        assert v == pytest.approx(ref[k], rel=1e-12), k


def test_direction_of_the_normalization_effect(rep):
    """Vorzeichen-Vorhersage aus O1: a_s > 1 (honeycomb, kagome) senkt die
    Crossings, a_s < 1 (triangular) hebt sie."""
    for key in ("honeycomb_phy042_wl", "kagome_phy033"):
        s = rep["sets"][key]
        for k, v in s["per_site"]["pairs"].items():
            w = s["per_area"]["pairs"][k]
            if v is not None and w is not None:
                assert w < v, (key, k)
    tri = rep["sets"]["triangular_phy030v02"]
    assert tri["per_area"]["wm_free_c"]["T"] > tri["per_site"]["wm_free_c"]["T"]


def test_parse_table_shapes():
    t, m, e = phy047.parse_table(phy047.SOURCES["honeycomb_phy032"])
    assert len(t) == 8 and sorted(m) == [12, 24, 48]
    assert all(np.all(e[L] > 0) for L in m)


def test_committed_report_if_present():
    p = ROOT / "results" / f"{phy047.STEM}.json"
    if not p.exists():
        pytest.skip("PHY047-Report noch nicht committet")
    committed = json.loads(p.read_text(encoding="utf-8"))
    regen = json.loads(json.dumps(phy047._clean(phy047.run())))
    assert committed == regen
