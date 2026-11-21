"""Gates fuer PHY053 (Freigabe/Produktion/Adjudikation W4 v03, reiner Runner).

(a) Freigabe: nur committeter, vorregistrierter PASS-Pilot (Entscheid aus den
    Zellen neu berechnet) UND nativer PHY049-Preflight.
(b) Produktion startet ohne Freigabe nicht (kein Job, kein Checkpoint).
(c) Adjudikation reicht PHY050.assess_production unveraendert durch; PHY053
    steht nicht im Provenienz-Fingerprint von PHY049/PHY050.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

import conftest
import test_phy050_correlation_fss as t50

phy053 = conftest._load(
    "phy053_w4_v03_runner",
    "261121 PHY053 w4 v03 release and adjudication runner v01.py")
p52, p49, p50 = phy053._p52, phy053._p49, phy053._p50


def _cells(passed=True):
    return [{"L": L, "T": T, "passed": passed,
             "reasons": [] if passed else ["x"]}
            for L in p52.PILOT_LS for T in p52.PILOT_TS]


def _pilot_report(passed=True):
    pilot = {"Ls": list(p52.PILOT_LS), "Ts": list(p52.PILOT_TS),
             "n_chains": 4, "n_therm": 1000, "n_meas": 4000,
             "cells": _cells(passed)}
    return {"pilot": pilot, "decision": p52.decision(pilot)}


def _preflight(ok=True):
    gates = {k: True for k in p50.W4V3_REQUIRED_PREFLIGHT_GATES}
    if not ok:
        gates["G1_geometry"] = False
    return {"gates": gates, "production_measurement_eligible": ok}


# (a) ------------------------------------------------------------------------

def test_release_only_with_pass_pilot_and_green_preflight():
    eligible = "ELIGIBLE_TO_START_PHY049_MEASUREMENT"
    assert phy053.release_decision(_pilot_report(), _preflight())["label"] \
        == eligible
    assert phy053.release_decision(_pilot_report(False), _preflight())[
        "label"] != eligible
    assert phy053.release_decision(_pilot_report(), _preflight(False))[
        "label"] != eligible
    assert phy053.release_decision(None, _preflight())["label"] != eligible
    assert phy053.release_decision({}, _preflight())["label"] != eligible


def test_a_stored_pass_label_is_not_trusted_without_the_cells():
    forged = _pilot_report(False)
    forged["decision"] = {"preregistered_design": True, "failed_cells": [],
                          "label": "MIXING_PILOT_PASS"}
    rel = phy053.release_decision(forged, _preflight())
    assert rel["label"] == "DO_NOT_START_PHY049_MEASUREMENT"
    assert "disagrees" in rel["reason"]


def test_a_shrunken_pilot_design_cannot_release():
    rep = _pilot_report()
    rep["pilot"]["n_meas"] = 400
    rep["decision"] = p52.decision(rep["pilot"])
    assert rep["decision"]["label"] == "MIXING_PILOT_FAIL"
    assert phy053.release_decision(rep, _preflight())["label"] \
        == "DO_NOT_START_PHY049_MEASUREMENT"


# (b) ------------------------------------------------------------------------

def test_production_does_not_start_without_release(monkeypatch, tmp_path):
    monkeypatch.setattr(phy053, "load_pilot_report",
                        lambda path=None: _pilot_report(False))
    monkeypatch.setattr(p49, "preflight", lambda: _preflight())
    monkeypatch.setattr(p49, "produce",
                        lambda **kw: pytest.fail("produce trotz Sperre"))
    with pytest.raises(RuntimeError, match="not released"):
        phy053.run_production(product_path=tmp_path / "p.json")
    assert not (tmp_path / "p.json").exists()


def test_production_uses_only_the_preregistered_defaults(monkeypatch,
                                                         tmp_path):
    seen = {}
    monkeypatch.setattr(phy053, "load_pilot_report",
                        lambda path=None: _pilot_report())
    monkeypatch.setattr(p49, "preflight", lambda: _preflight())

    def fake_produce(**kw):
        seen.update(kw)
        return {"checkpoint_status": "COMPLETE", "complete": True,
                "rows": []}
    monkeypatch.setattr(p49, "produce", fake_produce)
    out = phy053.run_production(product_path=tmp_path / "p.json")
    assert set(seen) == {"checkpoint_path", "resume"}   # keine Overrides
    assert out["release"]["label"] == "ELIGIBLE_TO_START_PHY049_MEASUREMENT"


# (c) ------------------------------------------------------------------------

def test_runner_is_not_in_the_production_fingerprint():
    for paths in (p49._RUNTIME_SOURCE_PATHS, p50._RUNTIME_SOURCE_PATHS):
        assert not any("PHY053" in p or "PHY052" in p for p in paths)


def test_adjudication_passes_phy050_through_unchanged(tmp_path):
    prod = t50._synthetic_product_for_validation()
    want = p50.assess_production(copy.deepcopy(prod), n_boot=20)
    rep = phy053.adjudicate(copy.deepcopy(prod), product_sha256="0" * 64,
                            n_boot=20)
    clean = phy053._clean(rep)
    assert clean["assessment"] == phy053._clean(want)
    # ein reduzierter Bootstrap kann nie freigeben (v03a Abschnitt 9)
    assert rep["assessment"]["decision"] == "INCONCLUSIVE"
    assert rep["assessment"]["physics_interpretation_enabled"] is False
    json.dumps(clean, allow_nan=False)
    phy053.write_text_report(clean, tmp_path / "a.txt")
    txt = (tmp_path / "a.txt").read_text(encoding="utf-8")
    assert "ENTSCHEID: INCONCLUSIVE" in txt


def test_adjudication_of_invalid_product_fails_closed():
    prod = t50._synthetic_product_for_validation()
    prod["rows"] = prod["rows"][:-1]
    rep = phy053.adjudicate(prod, n_boot=20)
    assert rep["assessment"]["gates"]["G0_input"] is False
    assert rep["assessment"]["decision"] == "INCONCLUSIVE"


def test_committed_adjudication_if_present():
    p = (Path(__file__).resolve().parents[1] / "results"
         / f"{phy053.ADJ_STEM}.json")
    if not p.exists():
        pytest.skip("W4-v03-Adjudikation noch nicht committet")
    rep = json.loads(p.read_text(encoding="utf-8"))
    assert rep["n_boot"] == p50.W4V3_N_BOOT
    assert rep["assessment"]["decision"] in (
        "SUPPORTED_LOW", "SUPPORTED_LITERATURE", "OVERLAP",
        "OUTSIDE_PREREGISTERED_HYPOTHESES", "INCONCLUSIVE")
