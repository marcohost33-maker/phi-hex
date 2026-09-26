"""Gates fuer PHY046 (W4 v02: Upsilon pro Flaeche, Wolff, HKS).

(a) Vorregistrierte Konstanten zeichengenau an die Spec v02 gebunden.
(b) Seed-Vertrag kollisionsfrei, HKS-Paare.
(c) Schaetzer-Orakel: HKS-Fit, WM-Fit mit freiem C, Jackknife-Pipeline.
(d) Mini-Produktion end-to-end (klein, schnell).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import conftest

phy046 = conftest._load(
    "phy046_w4_v02", "260926 PHY046 honeycomb w4 wolff area-helicity v01.py")

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "spec" / "260926 PHI HEX w4 honeycomb preregistration v02.md"


def test_constants_match_spec_v02():
    spec = SPEC.read_text(encoding="utf-8")
    assert ("L in {" + ", ".join(str(L) for L in phy046.W4V2_LADDER) + "}"
            in spec)
    g = phy046.W4V2_T_GRID
    assert len(g) == 16 and f"{g[0]:.3f} .. {g[-1]:.3f}, Schritt 0.005" in spec
    assert np.allclose(np.diff(g), 0.005)
    assert f"{phy046.W4V2_N_SEEDS} unabhaengige Seeds" in spec
    assert (f"je {phy046.W4V2_N_THERM} Thermalisierungs-" in spec
            and f"{phy046.W4V2_N_MEAS} Mess-Sweeps" in spec)
    assert "`47_000_000 + 1000*L + 10*t_idx + s`" in spec
    assert phy046.W4V2_SEED_BASE == 47_000_000
    assert "Untergrenze **0.003**" in spec
    assert phy046.W4V2_NORMALIZATION == "per_area"
    assert tuple(phy046.W4_BAND) == (0.560, 0.580)


def test_seed_contract_unique_and_in_range():
    seeds = {phy046.seed_for(L, k, s) for L in phy046.W4V2_LADDER
             for k in range(len(phy046.W4V2_T_GRID))
             for s in range(phy046.W4V2_N_SEEDS)}
    assert len(seeds) == (len(phy046.W4V2_LADDER) * len(phy046.W4V2_T_GRID)
                          * phy046.W4V2_N_SEEDS)
    assert max(seeds) < 2 ** 32


def test_hks_pairs_of_the_ladder():
    assert phy046.hks_pairs(phy046.W4V2_LADDER) == [
        (32, 64), (48, 96), (64, 128), (96, 192), (128, 256)]


def test_fit_hks_recovers_synthetic_parameters():
    Ls = np.array([32, 48, 64, 96, 128])
    Ts = 0.57 + 0.4 / np.log(1.3 * Ls) ** 2
    fit = phy046.fit_hks(Ls, Ts)
    assert fit["T_c"] == pytest.approx(0.57, abs=5e-4)
    assert fit["b"] == pytest.approx(1.3, rel=0.03)
    fit1 = phy046.fit_hks(Ls, 0.57 + 0.4 / np.log(Ls) ** 2, b_fixed=1.0)
    assert fit1["T_c"] == pytest.approx(0.57, abs=1e-12)
    assert phy046.fit_hks([32, 64], [0.6, 0.59]) is None


def _wm_curves(T0, C, t_grid, ladder, slope=4.0):
    """Synthetische Kurven, die bei T0 exakt der WM-Form folgen."""
    return {L: np.array([(2 * T0 / math.pi) * (1 + 1 / (2 * math.log(L) + C))
                         - slope * (T - T0) for T in t_grid]) for L in ladder}


def test_wm_free_c_fit_finds_the_wm_temperature():
    t = phy046.W4V2_T_GRID
    ladder = [32, 64, 128]
    curves = _wm_curves(0.58, 3.0, t, ladder)
    se = {L: np.full(len(t), 1e-3) for L in ladder}
    out = phy046.wm_free_c_fit(curves, se, t, ladder)
    assert out["T"] == pytest.approx(0.58, abs=0.003)
    assert out["at_grid_edge"] is False


def _synthetic_prod(ladder, T0=0.575, noise=1e-4, seed=3):
    t = list(phy046.W4V2_T_GRID)
    rng = np.random.default_rng(seed)
    base = _wm_curves(T0, 2.0, t, ladder)
    return {"t_grid": t, "ladder": list(ladder), "n_seeds": 4,
            "ups_area": {str(L): (base[L][:, None]
                                  + noise * rng.standard_normal((len(t), 4))
                                  ).tolist() for L in ladder},
            "ups_site": {str(L): (1.299 * base[L][:, None]
                                  + np.zeros((len(t), 4))).tolist()
                         for L in ladder}}


def test_analyse_stop_rule_s1_with_two_pairs():
    out = phy046.analyse(_synthetic_prod([32, 64, 128]), "ups_area")
    assert out["estimates"]["n_pairs_with_crossing"] == 2
    assert out["verdict"]["rule"] == "S1"       # < 3 Paare -> Stop


def test_analyse_recovers_exact_wm_temperature_on_full_ladder():
    """Exakte WM-Kurven: jedes Paar kreuzt bei T0; HKS, Varianten und
    Jackknife muessen T0 liefern und das Verdikt CONSISTENT (T0 in B)."""
    out = phy046.analyse(_synthetic_prod(phy046.W4V2_LADDER), "ups_area")
    assert out["estimates"]["n_pairs_with_crossing"] == 5
    for T in out["estimates"]["pairs"].values():
        assert T == pytest.approx(0.575, abs=1e-3)
    assert out["T_w4"] == pytest.approx(0.575, abs=2e-3)
    assert out["verdict"]["verdict"] == "CONSISTENT"


def test_mini_production_end_to_end(tmp_path):
    prod = phy046.produce(ladder=(4, 8), t_grid=(0.5, 0.6, 0.7), n_seeds=2,
                          n_therm=5, n_meas=20, max_workers=1)
    assert np.asarray(prod["ups_area"]["8"]).shape == (3, 2)
    ratio = (np.asarray(prod["ups_site"]["8"])
             / np.asarray(prod["ups_area"]["8"]))
    assert np.allclose(ratio, phy046.AREA_PER_SITE)
    rep = phy046._clean(phy046.build_report(prod))
    json.dumps(rep, allow_nan=False)
    phy046.write_text_report(rep, tmp_path / "r.txt")
    assert "Vorregistriertes Verdikt" in (tmp_path / "r.txt").read_text(
        encoding="utf-8")


def test_committed_w4_report_if_present():
    p = ROOT / "results" / (phy046.REPORT_STEM + ".json")
    if not p.exists():
        pytest.skip("W4-v02-Report noch nicht committet")
    rep = json.loads(p.read_text(encoding="utf-8"))
    pr = rep["production"]
    assert pr["ladder"] == list(phy046.W4V2_LADDER)
    assert pr["n_seeds"] == phy046.W4V2_N_SEEDS
    assert pr["n_meas"] == phy046.W4V2_N_MEAS
    assert rep["primary_per_area"]["verdict"]["verdict"] in (
        "CONSISTENT", "INCONSISTENT", "NEGATIVE_RESULT")
