"""Gates fuer PHY048 (post-hoc V&V der W4-v02-Pipeline auf square)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import conftest

phy048 = conftest._load(
    "phy048_square_validation", "260926 PHY048 square pipeline validation v01.py")

ROOT = Path(__file__).resolve().parents[1]


def test_same_ladder_and_statistics_as_w4():
    assert phy048.LADDER == phy048._p46.W4V2_LADDER
    assert (phy048.N_SEEDS, phy048.N_THERM, phy048.N_MEAS) == (
        phy048._p46.W4V2_N_SEEDS, phy048._p46.W4V2_N_THERM,
        phy048._p46.W4V2_N_MEAS)
    assert len(phy048.T_GRID) == 16
    assert phy048.T_GRID[0] < 0.8929 < phy048.T_GRID[-1]


def test_committed_report_if_present():
    p = ROOT / "results" / f"{phy048.STEM}.json"
    if not p.exists():
        pytest.skip("PHY048-Report noch nicht committet")
    rep = json.loads(p.read_text(encoding="utf-8"))
    assert rep["post_hoc"] is True
    pairs = rep["analysis"]["estimates"]["pairs"]
    assert len(pairs) == 5 and all(v is not None for v in pairs.values())
    # Befund, der das Protokoll-Lernen traegt: Paare <= 0.5 % bei 0.8929
    assert all(abs(v / 0.8929 - 1.0) <= 0.005 for v in pairs.values())
