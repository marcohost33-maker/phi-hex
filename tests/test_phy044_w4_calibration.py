"""Fast gates for PHY044 W4 honeycomb calibration."""
from __future__ import annotations

import pytest

import conftest

phy044 = conftest._load(
    "phy044_honeycomb_w4_calibration",
    "260927 PHY044 honeycomb w4 calibration v01.py",
)


def test_calibration_stream_namespace_is_deterministic_and_separate():
    assert phy044._stream_for_calibration(0) == 2_000_000
    assert phy044._stream_for_calibration(1) == 2_001_000
    assert phy044._stream_for_calibration(2) == 2_002_000
    with pytest.raises(ValueError):
        phy044._stream_for_calibration(-1)


def test_calibration_production_density_contract():
    assert phy044._prod_sweeps_for(100) == 30000
    assert phy044._prod_sweeps_for(500) == 30000
    assert phy044._prod_sweeps_for(1000) == 60000
    with pytest.raises(ValueError):
        phy044._prod_sweeps_for(0)


def test_calibration_inputs_are_fail_closed():
    assert phy044._validate_inputs(64, 2, (0.565, 0.575, 0.585)) == (
        0.565, 0.575, 0.585
    )
    with pytest.raises(ValueError, match="mindestens 2"):
        phy044._validate_inputs(64, 1, (0.575,))
    with pytest.raises(ValueError, match="L muss"):
        phy044._validate_inputs(4, 2, (0.575,))
    with pytest.raises(ValueError, match="streng aufsteigend"):
        phy044._validate_inputs(64, 2, (0.575, 0.565))
    with pytest.raises(ValueError):
        phy044._validate_inputs(64, 2, (0.575, 0.575))


def test_calibration_verdict_contract():
    common = dict(
        all_1t=True,
        leak_max=1e-6,
        uncovered_max=1e-6,
        max_y2_spread=0.01,
    )
    assert phy044._calibration_verdict(**common) == (
        "CALIBRATION_PASS_FOR_COST_MODEL"
    )

    x = dict(common)
    x["all_1t"] = False
    assert phy044._calibration_verdict(**x) == "CALIBRATION_FAIL_1T_NOT_ENGAGED"

    x = dict(common)
    x["leak_max"] = phy044.MASS_LIMIT
    assert phy044._calibration_verdict(**x) == "CALIBRATION_FAIL_COVERAGE"

    x = dict(common)
    x["uncovered_max"] = phy044.MASS_LIMIT
    assert phy044._calibration_verdict(**x) == "CALIBRATION_FAIL_COVERAGE"

    x = dict(common)
    x["max_y2_spread"] = phy044.SPREAD_LIMIT
    assert phy044._calibration_verdict(**x) == "CALIBRATION_NEEDS_MORE_WL"
