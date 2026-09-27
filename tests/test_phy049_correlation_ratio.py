"""Fast gates for PHY049 W4-v03 correlation-ratio staging."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import conftest

phy049 = conftest._load(
    "phy049_corr_ratio", "260927 PHY049 honeycomb correlation ratio v01.py")

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "spec" / (
    "260927 PHI HEX w4 honeycomb preregistration v03 correlation-ratio.md")


def test_preregistered_constants_match_spec():
    text = SPEC.read_text(encoding="utf-8")
    assert "L = {48, 72, 96, 144, 192}" in text
    assert "0.540 .. 0.610 in Schritten von 0.0025 (29 Punkte)" in text
    assert "n_seeds = 12" in text
    assert "n_therm = 1000 Sweeps" in text
    assert "n_meas = 4000 Sweeps" in text
    assert "49_000_000 + 1000*L + 100*t_idx + s" in text
    assert phy049.W4V3_LADDER == (48, 72, 96, 144, 192)
    assert len(phy049.W4V3_T_GRID) == 29
    assert phy049.W4V3_T_GRID[0] == pytest.approx(0.540)
    assert phy049.W4V3_T_GRID[-1] == pytest.approx(0.610)
    assert np.allclose(np.diff(phy049.W4V3_T_GRID), 0.0025)


def test_seed_contract_collision_free():
    seeds = {phy049.seed_for(L, k, s) for L in phy049.W4V3_LADDER
             for k in range(len(phy049.W4V3_T_GRID))
             for s in range(phy049.W4V3_N_SEEDS)}
    assert len(seeds) == (len(phy049.W4V3_LADDER)
                          * len(phy049.W4V3_T_GRID)
                          * phy049.W4V3_N_SEEDS)


def test_aligned_state_correlation_ratio_is_one_and_directions_average():
    for L in (8, 12):
        th = np.zeros(2 * L * L)
        assert phy049._py_corr(th, L, L // 4) == pytest.approx(1.0)
        assert phy049._py_corr(th, L, L // 2) == pytest.approx(1.0)


def test_ratio_is_ratio_of_means_not_mean_of_ratios():
    q = np.array([1.0, 3.0, 2.0])
    h = np.array([0.4, 2.1, 1.0])
    got, se = phy049.jackknife_ratio(q, h)
    assert got == pytest.approx(h.mean() / q.mean())
    assert got != pytest.approx(np.mean(h / q))
    assert se is not None and se > 0


def test_ratio_fails_closed_for_small_denominator():
    q = np.array([1e-9, -1e-9])
    h = np.array([1.0, 1.0])
    assert phy049.ratio_of_means(q, h) is None


def test_persistent_splay_requires_significant_persistent_tail():
    t = [0.54, 0.55, 0.56, 0.57, 0.58]
    r1 = [0.80, 0.79, 0.78, 0.75, 0.70]
    r2 = [0.80, 0.79, 0.75, 0.70, 0.64]
    e = [0.005] * 5
    assert phy049.persistent_splay(t, r1, e, r2, e) == pytest.approx(0.56)
    broken = list(r2)
    broken[4] = 0.76
    assert phy049.persistent_splay(t, r1, e, broken, e) is None


def test_preflight_allows_measurement_after_g4_but_not_physics_interpretation():
    out = phy049.preflight()
    assert out["gates"]["VAL_BIT_numba"] is True
    assert out["gates"]["G1_geometry"] is True
    assert out["gates"]["G2_aligned_limit"] is True
    assert out["gates"]["G3_seed_unique"] is True
    assert out["gates"]["G4_fss_recovery"] is True
    assert out["stage"] == "PREPRODUCTION_G4_VALIDATED"
    assert out["production_measurement_eligible"] is True
    assert out["overall_interpretation_enabled"] is False
    assert "NO_PHYSICS_INTERPRETATION" in out["claim_ceiling"]


def test_tiny_measurement_job_compiles_and_returns_finite_correlations():
    """Exercise the actual Wolff + correlation path, including Numba when present."""
    row = phy049._job((8, 0, 0.57, 0, 2, 3))
    assert row["L"] == 8
    assert row["t_idx"] == 0
    assert np.isfinite(row["g_quarter"])
    assert np.isfinite(row["g_half"])
    assert -1.0 <= row["g_quarter"] <= 1.0
    assert -1.0 <= row["g_half"] <= 1.0
    assert row["wall_s"] >= 0.0


def test_geometry_invariants_are_checked_directly():
    for L in (8, 12):
        got = phy049._geometry_invariants(L)
        assert got["translation_L_a1"] is True
        assert got["translation_L_a2"] is True
        assert got["same_sublattice_q_h"] is True


def test_numba_backend_is_bit_identical_when_available():
    assert phy049._backend_bit_identity() is True


def test_persistent_splay_rejects_nonfinite_and_wrong_sign():
    t = [0.54, 0.55, 0.56, 0.57]
    r1 = [0.80, 0.79, 0.76, 0.70]
    e = [0.005] * 4
    good = [0.80, 0.79, 0.72, 0.64]
    assert phy049.persistent_splay(t, r1, e, good, e) == pytest.approx(0.56)
    bad_nan = list(good)
    bad_nan[-1] = float("nan")
    assert phy049.persistent_splay(t, r1, e, bad_nan, e) is None
    wrong_sign = [0.80, 0.79, 0.82, 0.86]
    assert phy049.persistent_splay(t, r1, e, wrong_sign, e) is None


def test_zero_wall_budget_marks_every_job_unmeasured(monkeypatch):
    def must_not_run(_args):
        raise AssertionError("job must not start after budget stop")

    monkeypatch.setattr(phy049, "_job", must_not_run)
    out = phy049.produce(
        ladder=(8,), t_grid=(0.57,), n_seeds=2, n_therm=1, n_meas=1,
        max_workers=1, wall_budget_h=0.0,
    )
    assert out["rows"] == []
    assert out["complete"] is False
    assert len(out["unmeasured"]) == 2
    assert {r["reason"] for r in out["unmeasured"]} == {"WALL_BUDGET_STOP"}


def test_aggregate_rejects_duplicate_seed_wrong_temperature_and_wrong_rng_seed():
    base_rows = [
        {"L": 8, "t_idx": 0, "T": 0.57, "s": s,
         "seed": phy049.seed_for(8, 0, s),
         "g_quarter": 0.9 + 0.01 * s, "g_half": 0.8 + 0.01 * s}
        for s in range(3)
    ]
    prod = {"ladder": [8], "t_grid": [0.57], "n_seeds": 3, "rows": base_rows}
    assert phy049.aggregate(prod)["curves"]["8"]["R"][0] is not None

    dup = {**prod, "rows": [base_rows[0], base_rows[0], base_rows[2]]}
    assert phy049.aggregate(dup)["curves"]["8"]["R"][0] is None

    wrong_t = {**prod, "rows": [dict(r) for r in base_rows]}
    wrong_t["rows"][1]["T"] = 0.571
    assert phy049.aggregate(wrong_t)["curves"]["8"]["R"][0] is None

    wrong_seed = {**prod, "rows": [dict(r) for r in base_rows]}
    wrong_seed["rows"][1]["seed"] += 1
    assert phy049.aggregate(wrong_seed)["curves"]["8"]["R"][0] is None
