"""Gates fuer PHY052 (Mixing-Pilot vor W4 v03).

(a) Konstanten an die Spec und an W4 v03 gebunden; Seed-Namensraeume disjunkt.
(b) Diagnostik = stan-dev/posterior: Golden-Werte aus arviz-stats 1.3.3
    (unabhaengige Implementierung) auf eingefrorenen RandomState-Daten,
    analytische AR(1)-Orakel, Formel-Details (S + 1/4, biased Autokovarianz).
(c) Fail-closed: NaN/Inf/konstant -> NaN -> Gate rot.
(d) Kalt/heiss-Design erkennt einen gemeinsamen Thermalisierungs-Bias, den
    reine Kaltstarts nicht sehen (Negativ- und Positivkontrolle).
(e) Entscheid und Produktionsfreigabe: nur exaktes Design + alle Zellen.
(f) Mini-Lauf mit dem echten Kernel, Report-Round-Trip.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import conftest

phy052 = conftest._load(
    "phy052_mixing_pilot",
    "261121 PHY052 honeycomb correlation-ratio mixing pilot v01.py")
p49 = phy052._p49

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / phy052.SPEC


# ---------------------------------------------------------------------------
# (a) Bindung an Spec und W4 v03
# ---------------------------------------------------------------------------

def test_constants_match_spec_and_w4_v03():
    spec = SPEC.read_text(encoding="utf-8")
    assert "L in {48, 96, 192} x T in {0.560, 0.575, 0.595}" in spec
    assert phy052.PILOT_LS == (48, 96, 192)
    assert phy052.PILOT_TS == (0.560, 0.575, 0.595)
    assert "Je Zelle 4 unabhaengige Ketten" in spec
    assert phy052.PILOT_N_CHAINS == 4
    assert "1000 Thermalisierungs- + 4000 Mess-Sweeps" in spec
    assert (phy052.PILOT_N_THERM, phy052.PILOT_N_MEAS) == (
        p49.W4V3_N_THERM, p49.W4V3_N_MEAS) == (1000, 4000)
    assert "`79_000_000 + 1000*L + 100*t_idx + s`" in spec
    assert "`89_000_000 + 1000*L + 100*t_idx + s`" in spec
    assert ("**rhat(u) <= 1.01** und **ess_bulk(u) >= 400**\nund "
            "**|z_kh| <= 3.5**") in spec
    assert (phy052.MAX_RHAT, phy052.MIN_ESS_BULK,
            phy052.MAX_COLD_HOT_Z) == (1.01, 400.0, 3.5)
    assert set(phy052.PILOT_LS) <= set(p49.W4V3_LADDER)
    grid = [round(t, 4) for t in p49.W4V3_T_GRID]
    assert all(round(T, 4) in grid for T in phy052.PILOT_TS)


def test_seed_namespaces_are_disjoint_from_production_and_each_other():
    pilot, init, prod = set(), set(), set()
    for L in p49.W4V3_LADDER:
        for k in range(len(p49.W4V3_T_GRID)):
            for s in range(max(p49.W4V3_N_SEEDS, phy052.PILOT_N_CHAINS)):
                pilot.add(phy052.pilot_seed_for(L, k, s))
                init.add(phy052.pilot_init_seed_for(L, k, s))
                prod.add(p49.seed_for(L, k, s))
    assert not (pilot & prod) and not (init & prod) and not (pilot & init)
    assert max(pilot | init) < 2 ** 32


def test_chain_starts_are_balanced_cold_and_hot():
    modes = [phy052.init_mode_for(s) for s in range(phy052.PILOT_N_CHAINS)]
    assert modes == ["ordered", "random", "ordered", "random"]
    cold = phy052.initial_theta(5000, "ordered", 1)
    hot = phy052.initial_theta(5000, "random", 2)
    assert np.all(cold == 0.0)
    assert abs(np.mean(np.cos(hot))) < 0.05 and np.std(hot) > 1.7


# ---------------------------------------------------------------------------
# (b) Diagnostik = posterior; Golden-Werte aus arviz-stats 1.3.3
# ---------------------------------------------------------------------------

def _ar1(rs, rho, m, n):
    e = rs.standard_normal((m, n))
    x = np.empty((m, n))
    x[:, 0] = e[:, 0]
    s = math.sqrt(1.0 - rho * rho)
    for i in range(1, n):
        x[:, i] = rho * x[:, i - 1] + s * e[:, i]
    return x


# RandomState (NEP 19: eingefrorener Strom) -> plattformstabile Eingaben;
# Referenzwerte in tests/data/phy052_arviz_golden.json (arviz-stats 1.3.3).
_GOLDEN_CASES = 6


def _golden_input(seed, kind, param, m, n):
    rs = np.random.RandomState(seed)
    if kind == "ar1":
        return _ar1(rs, param, m, n)
    x = rs.standard_normal((m, n))
    if kind == "shift":
        x[0] += param
    else:
        x[1] *= param
    return x


GOLDEN_PATH = ROOT / "tests" / "data" / "phy052_arviz_golden.json"


def test_rhat_and_ess_bulk_match_independent_arviz_reference():
    """Golden-Werte aus arviz-stats (unabhaengige Implementierung derselben
    Vehtari-2021-Definitionen); Toleranz: Rundung."""
    golden = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
    assert golden["generator"].startswith("arviz-stats 1.3.3")
    assert len(golden["cases"]) == _GOLDEN_CASES
    for row in golden["cases"]:
        x = _golden_input(row["seed"], row["kind"], row["param"],
                          row["chains"], row["draws"])
        assert phy052.rhat(x) == pytest.approx(row["rhat"], abs=1e-12)
        assert phy052.ess_bulk(x) == pytest.approx(row["ess_bulk"],
                                                   rel=1e-10)


def test_golden_cases_cover_both_gate_outcomes():
    golden = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
    r = [c["rhat"] for c in golden["cases"]]
    assert min(r) < 1.01 < max(r)


def test_normal_score_denominator_is_s_plus_quarter():
    """posterior backtransform_ranks: (r - c) / (S - 2c + 1), c = 3/8."""
    x = np.arange(8.0).reshape(2, 4)
    z = phy052.z_scale(x)
    S = 8
    from scipy.stats import norm
    want = norm.ppf((np.arange(1, S + 1) - 0.375) / (S + 0.25))
    assert np.allclose(z.ravel(), want, atol=0, rtol=1e-15)
    assert not np.allclose(
        z.ravel(), norm.ppf((np.arange(1, S + 1) - 0.375) / (S - 0.25)))


def test_autocovariance_is_the_biased_geyer_estimator():
    rs = np.random.RandomState(3)
    x = rs.standard_normal(257)
    y = x - x.mean()
    n = x.size
    want = np.array([np.dot(y[:n - k], y[k:]) / n for k in range(n)])
    assert np.allclose(phy052.autocovariance(x), want, rtol=1e-10,
                       atol=1e-13)


def test_split_chains_follows_posterior_for_odd_lengths():
    x = np.arange(7.0).reshape(1, 7)
    assert phy052.split_chains(x).tolist() == [[0.0, 1.0, 2.0],
                                                [4.0, 5.0, 6.0]]
    x = np.arange(8.0).reshape(1, 8)
    assert phy052.split_chains(x).tolist() == [[0.0, 1.0, 2.0, 3.0],
                                                [4.0, 5.0, 6.0, 7.0]]


@pytest.mark.parametrize("rho", [0.0, 0.5, 0.8])
def test_ess_recovers_the_ar1_truth(rho):
    """ESS = N (1 - rho) / (1 + rho) fuer AR(1); 8 Ketten x 20000."""
    x = _ar1(np.random.RandomState(int(100 * rho) + 7), rho, 8, 20000)
    truth = x.size * (1 - rho) / (1 + rho)
    assert phy052.ess_bulk(x) == pytest.approx(truth, rel=0.06)


def test_stationary_chains_pass_and_shifted_chain_fails():
    rs = np.random.RandomState(5)
    ok = [phy052.rhat(rs.standard_normal((4, 4000))) for _ in range(20)]
    assert max(ok) < 1.01
    bad = rs.standard_normal((4, 4000))
    bad[2] += 0.5
    assert phy052.rhat(bad) > 1.01
    scale = rs.standard_normal((4, 4000))
    scale[3] *= 1.5
    assert phy052.rhat(scale) > 1.01     # nur die gefaltete Variante sieht das


def test_rhat_alone_misses_a_small_start_bias_that_the_z_gate_catches():
    """Begruendung des Zweistart-Gates (Spec Abschnitt 4): ein Versatz von
    0.15 sigma in einer Kette bleibt unter rhat 1.01, der Vergleich kalt
    gegen heiss in MCSE-Einheiten sieht ihn deutlich."""
    rs = np.random.RandomState(5)
    x = rs.standard_normal((4, 4000))
    x[2] += 0.15
    assert phy052.rhat(x) < 1.01
    q, h, modes = _synthetic_cell(0.0015, 0.0, seed=21)   # 0.15 sigma(h)
    d = phy052.cell_diagnostics(q, h, modes)
    assert d["rhat"] < 1.01 and d["ess_bulk"] >= 400
    assert abs(d["cold_minus_hot_z"]) > phy052.MAX_COLD_HOT_Z
    assert d["passed"] is False
    assert any("cold_minus_hot_z" in r for r in d["reasons"])


# ---------------------------------------------------------------------------
# (c) fail-closed
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", ["nan", "inf", "constant", "short", "1d"])
def test_invalid_input_gives_nan_and_a_red_gate(bad):
    x = np.random.RandomState(1).standard_normal((4, 100))
    if bad == "nan":
        x[2, 7] = math.nan
    elif bad == "inf":
        x[0, 0] = math.inf
    elif bad == "constant":
        x = np.full((4, 100), 0.3)
    elif bad == "short":
        x = x[:, :3]
    else:
        x = x[0]
    assert math.isnan(phy052.rhat(x)) and math.isnan(phy052.ess_bulk(x))


def test_single_chain_has_no_rhat():
    x = np.random.RandomState(2).standard_normal((1, 500))
    assert math.isnan(phy052.rhat(x))
    assert math.isfinite(phy052.ess_bulk(x))


# ---------------------------------------------------------------------------
# (d) Kalt/heiss erkennt gemeinsamen Thermalisierungs-Bias
# ---------------------------------------------------------------------------

def _synthetic_cell(bias_cold, bias_hot, seed=9, n=4000):
    """q, h je Kette um 0.70/0.40; eine Startart verschiebt h (zu kurze
    Thermalisierung). Rueckgabe (q, h, modes)."""
    rs = np.random.RandomState(seed)
    modes = [phy052.init_mode_for(s) for s in range(4)]
    q = 0.70 + 0.01 * rs.standard_normal((4, n))
    h = 0.40 + 0.01 * rs.standard_normal((4, n))
    for i, m in enumerate(modes):
        h[i] += bias_cold if m == "ordered" else bias_hot
    return q, h, modes


def test_cold_hot_design_detects_a_start_bias():
    q, h, modes = _synthetic_cell(0.004, -0.004)
    d = phy052.cell_diagnostics(q, h, modes)
    assert d["passed"] is False and d["rhat"] > 1.01
    assert abs(d["cold_minus_hot_z"]) > 5


def test_all_cold_starts_would_hide_the_same_bias():
    """Negativkontrolle des Designs: sind alle Ketten gleich (falsch)
    gestartet, sehen Zwischen-Ketten-Diagnosen nichts - genau deshalb
    verlangt die Spec heisse Starts."""
    q, h, _ = _synthetic_cell(0.004, 0.004)
    d = phy052.cell_diagnostics(q, h, ["ordered"] * 4 + [])
    assert d["rhat"] < 1.01


def test_unbiased_cell_passes():
    q, h, modes = _synthetic_cell(0.0, 0.0)
    d = phy052.cell_diagnostics(q, h, modes)
    assert d["passed"] is True and d["reasons"] == []
    assert d["ess_bulk"] >= 400


# ---------------------------------------------------------------------------
# (e) Entscheid und Freigabe
# ---------------------------------------------------------------------------

def _pilot(passed=True, drop=False, **over):
    cells = [{"L": L, "T": T, "passed": passed, "reasons": [] if passed
              else ["x"]} for L in phy052.PILOT_LS for T in phy052.PILOT_TS]
    if drop:
        cells = cells[:-1]
    base = {"Ls": list(phy052.PILOT_LS), "Ts": list(phy052.PILOT_TS),
            "n_chains": 4, "n_therm": 1000, "n_meas": 4000, "cells": cells}
    base.update(over)
    return base


def _preflight(**gates):
    g = {k: True for k in p49._p50.W4V3_REQUIRED_PREFLIGHT_GATES}
    g.update(gates)
    return {"gates": g, "production_measurement_eligible": all(
        v is True for v in g.values())}


def test_pass_needs_the_exact_design_and_every_cell():
    assert phy052.decision(_pilot())["label"] == "MIXING_PILOT_PASS"
    assert phy052.decision(_pilot(False))["label"] == "MIXING_PILOT_FAIL"
    assert phy052.decision(_pilot(drop=True))["label"] == "MIXING_PILOT_FAIL"
    for over in ({"n_meas": 400}, {"n_therm": 100}, {"n_chains": 2},
                 {"Ls": [48, 96]}, {"Ts": [0.575]}):
        assert phy052.decision(_pilot(**over))["label"] == \
            "MIXING_PILOT_FAIL", over
    one_bad = _pilot()
    one_bad["cells"][4]["passed"] = None           # kein literal True
    assert phy052.decision(one_bad)["label"] == "MIXING_PILOT_FAIL"


def test_production_release_needs_pilot_and_native_preflight():
    ok, bad = phy052.decision(_pilot()), phy052.decision(_pilot(False))
    eligible = "ELIGIBLE_TO_START_PHY049_MEASUREMENT"
    assert phy052.production_release(ok, _preflight()) == eligible
    assert phy052.production_release(bad, _preflight()) != eligible
    assert phy052.production_release(ok, _preflight(G4_fss_recovery=False)) \
        != eligible
    extra = _preflight()
    extra["gates"]["G9_extra"] = True
    assert phy052.production_release(ok, extra) != eligible
    truthy = _preflight()
    truthy["gates"]["VAL_BIT_numba"] = 1                # nicht literal True
    assert phy052.production_release(ok, truthy) != eligible
    assert phy052.production_release(ok, {}) != eligible


# ---------------------------------------------------------------------------
# (f) echter Kernel, Mini-Lauf, Report-Round-Trip
# ---------------------------------------------------------------------------

def test_mini_pilot_with_the_real_kernel_round_trips(tmp_path):
    pilot = phy052.run_pilot(Ls=(8,), Ts=(0.575,), n_chains=4, n_therm=20,
                             n_meas=200, max_workers=1)
    c = pilot["cells"][0]
    assert c["init_modes"] == ["ordered", "random", "ordered", "random"]
    assert c["seeds"] == [phy052.pilot_seed_for(8, 14, s) for s in range(4)]
    assert all(math.isfinite(c[k]) for k in ("rhat", "ess_bulk",
                                             "ratio_of_means"))
    rep = phy052._clean(phy052.build_report(pilot))
    back = json.loads(json.dumps(rep, allow_nan=False))
    assert back["decision"]["label"] == "MIXING_PILOT_FAIL"   # kein Design
    assert back["physics_interpretation_enabled"] is False
    assert set(back["runtime_provenance"]["source_sha256"]) == set(
        phy052._RUNTIME_SOURCE_PATHS)
    phy052.write_text_report(back, tmp_path / "r.txt")
    assert "MIXING_PILOT_FAIL" in (tmp_path / "r.txt").read_text(
        encoding="utf-8")


def test_pilot_chain_is_deterministic_and_uses_the_gated_kernel(monkeypatch):
    a = phy052._chain_job((8, 0.575, 1, 5, 30))
    b = phy052._chain_job((8, 0.575, 1, 5, 30))
    assert np.array_equal(a["q"], b["q"]) and np.array_equal(a["h"], b["h"])
    if p49.HAVE_NUMBA:
        monkeypatch.setattr(p49, "_VAL_BIT_OK", None)
        monkeypatch.setattr(p49, "_backend_bit_identity", lambda: False)
        with pytest.raises(RuntimeError, match="VAL-BIT failed"):
            phy052._chain_job((8, 0.575, 0, 5, 30))


def test_committed_pilot_report_if_present():
    p = ROOT / "results" / f"{phy052.REPORT_STEM}.json"
    if not p.exists():
        pytest.skip("PHY052-Report noch nicht committet")
    rep = json.loads(p.read_text(encoding="utf-8"))
    assert phy052.is_preregistered_design(rep["pilot"])
    assert rep["decision"] == phy052.decision(rep["pilot"])
    assert rep["thresholds"] == {"max_rhat": 1.01, "min_ess_bulk": 400.0,
                                 "max_cold_hot_z": 3.5,
                                 "gate_observable": "ratio_influence"}
    assert rep["runtime_provenance"]["thread_env"]["OPENBLAS_NUM_THREADS"] \
        == "1"
