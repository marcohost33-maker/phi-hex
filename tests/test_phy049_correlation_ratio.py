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


def test_preflight_is_fail_closed_until_fss_recovery_exists():
    out = phy049.preflight()
    assert out["gates"]["G1_geometry"] is True
    assert out["gates"]["G3_seed_unique"] is True
    assert out["gates"]["G4_fss_recovery"] is False
    assert out["overall_interpretation_enabled"] is False
    assert "NO_PHYSICS_INTERPRETATION" in out["claim_ceiling"]
