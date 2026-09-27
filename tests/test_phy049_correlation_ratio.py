"""Fast gates for PHY049 W4-v03 correlation-ratio staging."""
from __future__ import annotations

import json
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
    t = [0.54, 0.55, 0.56, 0.57, 0.58]
    r1 = [0.80, 0.79, 0.76, 0.70, 0.66]
    e = [0.005] * 5
    # Drei persistente, signifikante negative D-Punkte ab T=0.56.
    good = [0.80, 0.79, 0.72, 0.64, 0.58]
    assert phy049.persistent_splay(t, r1, e, good, e) == pytest.approx(0.56)
    bad_nan = list(good)
    bad_nan[-1] = float("nan")
    assert phy049.persistent_splay(t, r1, e, bad_nan, e) is None
    wrong_sign = [0.80, 0.79, 0.82, 0.86, 0.88]
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


def test_direct_job_cannot_bypass_failed_val_bit(monkeypatch):
    """Numba availability alone must never authorize a production kernel."""
    if not phy049.HAVE_NUMBA:
        pytest.skip("Numba not installed in this environment")
    monkeypatch.setattr(phy049, "_VAL_BIT_OK", None)
    monkeypatch.setattr(phy049, "_backend_bit_identity", lambda: False)
    with pytest.raises(RuntimeError, match="VAL-BIT failed"):
        phy049._job((8, 0, 0.57, 0, 2, 3))


def test_production_backend_caches_successful_val_bit(monkeypatch):
    """One interpreter/process validates VAL-BIT once, then reuses the verdict."""
    if not phy049.HAVE_NUMBA:
        pytest.skip("Numba not installed in this environment")
    calls = {"n": 0}

    def ok():
        calls["n"] += 1
        return True

    monkeypatch.setattr(phy049, "_VAL_BIT_OK", None)
    monkeypatch.setattr(phy049, "_backend_bit_identity", ok)
    assert phy049._production_run_backend() is phy049._nb_run_corr
    assert phy049._production_run_backend() is phy049._nb_run_corr
    assert calls["n"] == 1


def test_public_produce_cannot_bypass_failed_preflight(monkeypatch):
    monkeypatch.setattr(
        phy049,
        "preflight",
        lambda: {
            "production_measurement_eligible": False,
            "gates": {
                "VAL_BIT_numba": True,
                "G1_geometry": True,
                "G2_aligned_limit": True,
                "G3_seed_unique": True,
                "G4_fss_recovery": False,
            },
        },
    )
    monkeypatch.setattr(
        phy049,
        "_job",
        lambda _args: (_ for _ in ()).throw(
            AssertionError("job must not run after failed preflight")
        ),
    )
    with pytest.raises(RuntimeError, match="G4_fss_recovery"):
        phy049.produce(
            ladder=(8,),
            t_grid=(0.57,),
            n_seeds=2,
            n_therm=1,
            n_meas=1,
            max_workers=1,
            wall_budget_h=0.0,
        )


def test_product_records_preflight_gate_evidence(monkeypatch):
    gates = {
        "VAL_BIT_numba": True,
        "G1_geometry": True,
        "G2_aligned_limit": True,
        "G3_seed_unique": True,
        "G4_fss_recovery": True,
    }
    monkeypatch.setattr(
        phy049,
        "preflight",
        lambda: {
            "production_measurement_eligible": True,
            "gates": gates,
        },
    )
    out = phy049.produce(
        ladder=(8,),
        t_grid=(0.57,),
        n_seeds=1,
        n_therm=1,
        n_meas=1,
        max_workers=1,
        wall_budget_h=0.0,
    )
    assert out["preflight_gates"] == gates



def _green_preflight_for_checkpoint_tests():
    gates = {
        "VAL_BIT_numba": True,
        "G1_geometry": True,
        "G2_aligned_limit": True,
        "G3_seed_unique": True,
        "G4_fss_recovery": True,
    }
    return {"production_measurement_eligible": True, "gates": gates}


def _fake_checkpoint_job(args):
    L, t_idx, T, s, _n_therm, _n_meas = args
    return {
        "L": L,
        "t_idx": t_idx,
        "T": T,
        "s": s,
        "seed": phy049.seed_for(L, t_idx, s),
        "g_quarter": 0.8 + 0.001 * s,
        "g_half": 0.7 + 0.001 * s,
        "wall_s": 0.01,
    }


def test_production_checkpoint_is_atomic_complete_and_resumable(
    tmp_path, monkeypatch
):
    checkpoint = tmp_path / "phy049-production.json"
    monkeypatch.setattr(
        phy049, "preflight", _green_preflight_for_checkpoint_tests
    )
    monkeypatch.setattr(phy049, "_job", _fake_checkpoint_job)

    out = phy049.produce(
        ladder=(8,),
        t_grid=(0.57,),
        n_seeds=2,
        n_therm=1,
        n_meas=1,
        max_workers=1,
        wall_budget_h=1.0,
        checkpoint_path=checkpoint,
    )
    assert out["complete"] is True
    assert out["checkpoint_status"] == "COMPLETE"
    assert checkpoint.exists()
    assert not checkpoint.with_name(checkpoint.name + ".tmp").exists()
    persisted = json.loads(checkpoint.read_text(encoding="utf-8"))
    assert persisted == out

    def must_not_recompute(_args):
        raise AssertionError("completed checkpoint must not recompute jobs")

    monkeypatch.setattr(phy049, "_job", must_not_recompute)
    resumed = phy049.produce(
        ladder=(8,),
        t_grid=(0.57,),
        n_seeds=2,
        n_therm=1,
        n_meas=1,
        max_workers=1,
        wall_budget_h=1.0,
        checkpoint_path=checkpoint,
        resume=True,
    )
    assert resumed == out


def test_resume_rejects_partial_lattice_block(tmp_path, monkeypatch):
    checkpoint = tmp_path / "partial.json"
    monkeypatch.setattr(
        phy049, "preflight", _green_preflight_for_checkpoint_tests
    )
    monkeypatch.setattr(phy049, "_job", _fake_checkpoint_job)
    out = phy049.produce(
        ladder=(8,),
        t_grid=(0.57,),
        n_seeds=2,
        n_therm=1,
        n_meas=1,
        max_workers=1,
        wall_budget_h=1.0,
        checkpoint_path=checkpoint,
    )
    out["rows"] = out["rows"][:-1]
    out["complete"] = False
    out["checkpoint_status"] = "IN_PROGRESS"
    checkpoint.write_text(json.dumps(out), encoding="utf-8")

    with pytest.raises(RuntimeError, match="partial/invalid L=8 block"):
        phy049.produce(
            ladder=(8,),
            t_grid=(0.57,),
            n_seeds=2,
            n_therm=1,
            n_meas=1,
            max_workers=1,
            wall_budget_h=1.0,
            checkpoint_path=checkpoint,
            resume=True,
        )


def test_checkpoint_requires_explicit_resume(tmp_path, monkeypatch):
    checkpoint = tmp_path / "existing.json"
    checkpoint.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(
        phy049, "preflight", _green_preflight_for_checkpoint_tests
    )
    with pytest.raises(FileExistsError, match="resume=True"):
        phy049.produce(
            ladder=(8,),
            t_grid=(0.57,),
            n_seeds=1,
            n_therm=1,
            n_meas=1,
            max_workers=1,
            wall_budget_h=1.0,
            checkpoint_path=checkpoint,
        )


def test_wall_budget_stop_checkpoint_is_terminal(tmp_path, monkeypatch):
    checkpoint = tmp_path / "budget-stop.json"
    monkeypatch.setattr(
        phy049, "preflight", _green_preflight_for_checkpoint_tests
    )
    monkeypatch.setattr(
        phy049,
        "_job",
        lambda _args: (_ for _ in ()).throw(
            AssertionError("zero-budget campaign must not start a job")
        ),
    )
    out = phy049.produce(
        ladder=(8,),
        t_grid=(0.57,),
        n_seeds=2,
        n_therm=1,
        n_meas=1,
        max_workers=1,
        wall_budget_h=0.0,
        checkpoint_path=checkpoint,
    )
    assert out["checkpoint_status"] == "WALL_BUDGET_STOP"
    assert out["complete"] is False
    assert len(out["unmeasured"]) == 2

    with pytest.raises(RuntimeError, match="cannot be resumed"):
        phy049.produce(
            ladder=(8,),
            t_grid=(0.57,),
            n_seeds=2,
            n_therm=1,
            n_meas=1,
            max_workers=1,
            wall_budget_h=0.0,
            checkpoint_path=checkpoint,
            resume=True,
        )


def test_resume_contract_mismatch_fails_closed(tmp_path, monkeypatch):
    checkpoint = tmp_path / "contract.json"
    monkeypatch.setattr(
        phy049, "preflight", _green_preflight_for_checkpoint_tests
    )
    monkeypatch.setattr(phy049, "_job", _fake_checkpoint_job)
    phy049.produce(
        ladder=(8,),
        t_grid=(0.57,),
        n_seeds=2,
        n_therm=1,
        n_meas=1,
        max_workers=1,
        wall_budget_h=1.0,
        checkpoint_path=checkpoint,
    )
    with pytest.raises(RuntimeError, match="contract does not match"):
        phy049.produce(
            ladder=(8,),
            t_grid=(0.57,),
            n_seeds=3,
            n_therm=1,
            n_meas=1,
            max_workers=1,
            wall_budget_h=1.0,
            checkpoint_path=checkpoint,
            resume=True,
        )
