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


def _synthetic_product_for_validation():
    rows = []
    for L in phy050.W4V3_LADDER:
        for k, T in enumerate(phy050.W4V3_T_GRID):
            for s in range(12):
                rows.append({
                    "L": L,
                    "t_idx": k,
                    "T": T,
                    "s": s,
                    "seed": phy050.W4V3_SEED_BASE + 1000 * L + 100 * k + s,
                    "g_quarter": 0.90 + 1e-4 * s,
                    "g_half": 0.80 + 2e-4 * s,
                })
    return {
        "ladder": list(phy050.W4V3_LADDER),
        "t_grid": list(phy050.W4V3_T_GRID),
        "n_seeds": phy050.W4V3_N_SEEDS,
        "n_therm": phy050.W4V3_N_THERM,
        "n_meas": phy050.W4V3_N_MEAS,
        "rows": rows,
        "complete": True,
        "unmeasured": [],
    }


def test_product_groups_validate_temperature_seed_and_complete_flag():
    prod = _synthetic_product_for_validation()
    assert phy050._product_groups(prod) is not None

    wrong_t = {**prod, "rows": [dict(r) for r in prod["rows"]]}
    wrong_t["rows"][0]["T"] += 0.001
    assert phy050._product_groups(wrong_t) is None

    wrong_seed = {**prod, "rows": [dict(r) for r in prod["rows"]]}
    wrong_seed["rows"][0]["seed"] += 1
    assert phy050._product_groups(wrong_seed) is None

    incomplete = {**prod, "complete": False}
    assert phy050._product_groups(incomplete) is None


def test_decision_label_respects_full_hypothesis_overlap():
    assert phy050.decision_label(0.565, 0.001) == "SUPPORTED_LOW"
    assert phy050.decision_label(0.578, 0.001) == "SUPPORTED_LITERATURE"
    assert phy050.decision_label(0.571, 0.0001) == "OVERLAP"
    assert phy050.decision_label(0.5715, 0.003) == "OVERLAP"
    assert phy050.decision_label(0.565, 0.011) == "INCONCLUSIVE"
    assert phy050.decision_label(0.578, 0.001, gates_passed=False) == "INCONCLUSIVE"


def test_assess_production_fails_closed_on_incomplete_input():
    prod = _synthetic_product_for_validation()
    prod["complete"] = False
    out = phy050.assess_production(prod, n_boot=20)
    assert out["gates"]["G0_input"] is False
    assert out["decision"] == "INCONCLUSIVE"
    assert out["physics_interpretation_enabled"] is False


def test_g0_requires_exact_preregistered_budget_and_row_count():
    prod = _synthetic_product_for_validation()
    assert phy050._product_groups(prod) is not None

    wrong_seeds = {**prod, "n_seeds": 2}
    assert phy050._product_groups(wrong_seeds) is None

    wrong_therm = {**prod, "n_therm": phy050.W4V3_N_THERM - 1}
    assert phy050._product_groups(wrong_therm) is None

    wrong_meas = {**prod, "n_meas": phy050.W4V3_N_MEAS - 1}
    assert phy050._product_groups(wrong_meas) is None

    extra = {**prod, "rows": [dict(r) for r in prod["rows"]]}
    extra["rows"].append(dict(extra["rows"][0]))
    assert phy050._product_groups(extra) is None

    missing_complete = dict(prod)
    missing_complete.pop("complete")
    assert phy050._product_groups(missing_complete) is None

    missing_unmeasured = dict(prod)
    missing_unmeasured.pop("unmeasured")
    assert phy050._product_groups(missing_unmeasured) is None


@pytest.mark.parametrize("value", [None, "not-a-number", float("nan"), float("inf")])
def test_malformed_seed_moments_fail_g0_without_raising(value):
    prod = _synthetic_product_for_validation()
    prod["rows"] = [dict(r) for r in prod["rows"]]
    prod["rows"][0]["g_quarter"] = value
    assert phy050._product_groups(prod) is None
    assert phy050.product_to_aggregate(prod) is None
    out = phy050.assess_production(prod)
    assert out["gates"]["G0_input"] is False
    assert out["decision"] == "INCONCLUSIVE"


def test_nonpreregistered_bootstrap_cannot_enable_production_adjudication():
    prod = _synthetic_product_for_validation()
    out = phy050.assess_production(prod, n_boot=20)
    assert out["gates"]["G0_input"] is True
    assert out["gates"]["G5_power"] is False
    assert out["gates"]["G6_robustness"] is False
    assert out["bootstrap_required"] == phy050.W4V3_N_BOOT == 1000
    assert out["decision"] == "INCONCLUSIVE"
    assert out["physics_interpretation_enabled"] is False


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("n_seeds", "12"),
        ("n_seeds", 12.0),
        ("n_therm", "1000"),
        ("n_therm", 1000.5),
        ("n_meas", "4000"),
        ("n_meas", 4000.5),
    ],
)
def test_g0_rejects_coerced_metadata_types(field, value):
    prod = _synthetic_product_for_validation()
    prod[field] = value
    assert phy050._product_groups(prod) is None


def test_g0_rejects_stringified_row_identity_and_temperature():
    prod = _synthetic_product_for_validation()
    for field in ("L", "t_idx", "s", "seed", "T"):
        bad = {**prod, "rows": [dict(r) for r in prod["rows"]]}
        bad["rows"][0][field] = str(bad["rows"][0][field])
        assert phy050._product_groups(bad) is None


def test_float_bootstrap_count_fails_closed_in_production_assessor():
    prod = _synthetic_product_for_validation()
    out = phy050.assess_production(prod, n_boot=1000.0)
    assert out["decision"] == "INCONCLUSIVE"
    assert out["physics_interpretation_enabled"] is False
