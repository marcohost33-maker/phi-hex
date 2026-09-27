"""Fast correctness gates for PHY050 W4-v03a deterministic FSS."""
from __future__ import annotations

from pathlib import Path

import pytest

import conftest

phy050 = conftest._load(
    "phy050_corr_ratio_fss", "260927 PHY050 correlation ratio fss recovery v01.py"
)

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "spec" / "260927 PHI HEX w4 v03a fss estimator hardening.md"


def test_amendment_binds_estimator_constants():
    text = SPEC.read_text(encoding="utf-8")
    assert "[0.585, 0.610]" in text
    assert "T0 in [0.550, 0.582]" in text
    assert "80 logarithmically spaced" in text
    assert "profile width `<= 0.012`" in text
    assert "|T_hat-T0| <= 0.003" in text
    assert phy050.FSS_T_WINDOW == (0.585, 0.610)
    assert phy050.FSS_T0_BOUNDS == (0.550, 0.582)
    assert phy050.FSS_MAX_PROFILE_WIDTH == pytest.approx(0.012)
    assert phy050.G4_MAX_ABS_ERROR == pytest.approx(0.003)


def test_g4_positive_recovery_and_adversarial_nulls_pass():
    out = phy050.g4_synthetic_recovery()
    assert out["passed"] is True
    assert out["missing_input_fail_closed"] is True
    assert out["flat_null_rejected_as_unidentifiable"] is True
    assert out["noncollapse_null_rejected_as_unidentifiable"] is True
    for case in out["recovery_cases"]:
        assert case["pass"] is True
        assert case["quotable"] is True
        assert case["abs_error"] <= phy050.G4_MAX_ABS_ERROR


def test_profile_null_cannot_create_spurious_quotable_tbkt():
    fit = phy050.fit_collapse(phy050.flat_null_aggregate())
    assert fit is not None
    assert fit.identifiable is False
    assert fit.quotable is False


def test_missing_required_size_fails_closed():
    agg = phy050.synthetic_aggregate(0.573, 0.90, seed=991)
    agg["ladder"].remove(192)
    agg["curves"].pop("192")
    assert phy050.fit_collapse(agg) is None


def test_rejected_pair_does_not_pollute_score(monkeypatch):
    monkeypatch.setattr(phy050, "FSS_MIN_PAIR_POINTS", 10_000)
    agg = phy050.synthetic_aggregate(0.573, 0.90, seed=992)
    assert phy050.collapse_score(agg, 0.573, 0.90) is None
