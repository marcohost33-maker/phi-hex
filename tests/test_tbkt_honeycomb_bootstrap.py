"""PHY032 Honeycomb WM-Fit + Seed-Bootstrap: Korrektheits-Gates.

Vier Ebenen (alle fails-before-faehig, ohne teures Monte-Carlo):

(a) BOOTSTRAP-ORAKEL: ein Estimator mit BEKANNTEM Wert auf synthetischen
    per-Seed-Daten -> Bootstrap-Punktschaetzer == bekannter Wert und das CI
    MUSS ihn ueberdecken. Ein degenerierter Estimator (immer None) MUSS ein
    None-CI liefern (kein stiller 0-Balken - das Silent-Failure-Gate).
(b) JACKKNIFE-ORAKEL: leave-one-out SE schrumpft mit kleinerer Seed-Streuung;
    bei identischen Seeds ist die SE exakt 0.
(c) ESTIMATOR-VERDRAHTUNG: die Sandvik-/WM-Estimator-Fabriken rekonstruieren
    aus WM-konformen means_LT das bekannte T_true (Wiederverwendung der
    PHY030-v02-Fit-Funktionen, hier ueber die means_LT-Schnittstelle).
(d) REPRODUZIERBARKEIT: measure_honeycomb_seedwise teilt den RNG-Stream-
    Vertrag von PHY031 v01 -> per-Seed-Mittel == PHY031-v01-Mittel (slow).
"""
from __future__ import annotations

import math

import numpy as np

from phy032_honeycomb_wm_bootstrap import (
    REF_HONEYCOMB_DEDICATED,
    REF_HONEYCOMB_MULTI,
    bootstrap_tbkt_over_seeds,
    jackknife_tbkt_over_seeds,
    make_sandvik_pair_estimator,
    make_wm_fit_estimator,
    measure_honeycomb_seedwise,
    wm_form,
)


# ---------------------------------------------------------------------------
# Synthetik-Helfer: per-Seed-Wuerfel nach WM-Form (bekanntes T_true).
# ---------------------------------------------------------------------------

def _synthetic_cube(Ls, T_grid, T_true, C_true, n_seeds, noise, seed,
                    slope=0.3):
    """per_seed_cube[L][T] = Liste von n_seeds Upsilon-Werten, WM-konform.

    Bei T=T_true exakt die WM-Form; ausserhalb eine glatte L-abhaengige
    Abweichung (slope*(T-T_true)*lnL), die die Paar-/Fit-Bedingung nur bei
    T_true erfuellt. Reproduzierbares Gauss-Rauschen je Seed.
    """
    rng = np.random.default_rng(seed)
    cube = {L: {} for L in Ls}
    for L in Ls:
        for T in T_grid:
            base = wm_form(T_true, L, C_true) + slope * (T - T_true) * math.log(L)
            cube[L][T] = [float(base + rng.normal(0.0, noise))
                          for _ in range(n_seeds)]
    return cube


# ---------------------------------------------------------------------------
# (a) Bootstrap-Orakel
# ---------------------------------------------------------------------------

def test_bootstrap_recovers_constant_estimator():
    """Bootstrap eines konstanten Estimators: Punkt == Konstante, CI eng um
    sie. (Triviales aber scharfes Orakel fuer die Resampling-Maschinerie.)"""
    cube = {12: {0.57: [1.0] * 8}, 24: {0.57: [1.0] * 8}}

    def const_est(means_LT):
        return 0.5800

    out = bootstrap_tbkt_over_seeds(cube, const_est, n_seeds=8, n_boot=500)
    assert abs(out["point"] - 0.58) < 1e-12
    assert out["ci_lower"] is not None and out["ci_upper"] is not None
    assert abs(out["ci_lower"] - 0.58) < 1e-9
    assert abs(out["ci_upper"] - 0.58) < 1e-9


def test_bootstrap_ci_covers_true_value_pair_estimator():
    """Seed-Bootstrap des Sandvik-Paar-Estimators: CI ueberdeckt T_true.

    fails-before: falscher Resample-Mechanismus / kaputter Estimator verfehlt
    die Ueberdeckung oder liefert kein endliches CI.
    """
    T_true, C_true = 0.573, 2.0
    Ls = [8, 16]
    T_grid = [round(T_true - 0.02 + 0.005 * k, 4) for k in range(9)]
    cube = _synthetic_cube(Ls, T_grid, T_true, C_true, n_seeds=8,
                           noise=2e-3, seed=4242)
    est = make_sandvik_pair_estimator(8, 16, T_grid)
    out = bootstrap_tbkt_over_seeds(cube, est, n_seeds=8, n_boot=800,
                                    master_seed=42)
    assert out["point"] is not None and math.isfinite(out["point"])
    assert abs(out["point"] - T_true) < 0.01, out["point"]
    assert out["ci_lower"] is not None
    assert out["ci_lower"] <= T_true <= out["ci_upper"], out


def test_bootstrap_degenerate_estimator_yields_none_ci():
    """Silent-Failure-Gate: ein Estimator, der immer None liefert (kein
    Nulldurchgang), darf KEIN Fake-0-CI erzeugen, sondern None markieren."""
    cube = {8: {0.57: [0.5] * 8}, 16: {0.57: [0.5] * 8}}

    def none_est(means_LT):
        return None

    out = bootstrap_tbkt_over_seeds(cube, none_est, n_seeds=8, n_boot=200)
    assert out["ci_lower"] is None and out["ci_upper"] is None
    assert out["n_valid"] == 0


def test_bootstrap_is_deterministic():
    """Gleicher master_seed -> bit-identisches CI (Reproduzierbarkeit)."""
    T_true, C_true = 0.576, 1.5
    Ls = [8, 16]
    T_grid = [round(T_true - 0.02 + 0.005 * k, 4) for k in range(9)]
    cube = _synthetic_cube(Ls, T_grid, T_true, C_true, n_seeds=8,
                           noise=2e-3, seed=11)
    est = make_sandvik_pair_estimator(8, 16, T_grid)
    a = bootstrap_tbkt_over_seeds(cube, est, n_seeds=8, n_boot=400,
                                  master_seed=7)
    b = bootstrap_tbkt_over_seeds(cube, est, n_seeds=8, n_boot=400,
                                  master_seed=7)
    assert a["ci_lower"] == b["ci_lower"]
    assert a["ci_upper"] == b["ci_upper"]


# ---------------------------------------------------------------------------
# (b) Jackknife-Orakel
# ---------------------------------------------------------------------------

def test_jackknife_zero_se_for_identical_seeds():
    """Identische Seeds -> jeder leave-one-out-Schaetzer gleich -> SE = 0."""
    Ls = [8, 16]
    T_grid = [round(0.573 - 0.02 + 0.005 * k, 4) for k in range(9)]
    # noise=0 -> alle Seeds identisch.
    cube = _synthetic_cube(Ls, T_grid, 0.573, 2.0, n_seeds=8, noise=0.0,
                           seed=1)
    est = make_sandvik_pair_estimator(8, 16, T_grid)
    out = jackknife_tbkt_over_seeds(cube, est, n_seeds=8)
    assert out["jack_se"] is not None
    assert out["jack_se"] < 1e-9, out


def test_jackknife_mean_near_true_value():
    """Jackknife-Mittel des Paar-Estimators nahe T_true (WM-Synthetik)."""
    T_true = 0.573
    Ls = [8, 16]
    T_grid = [round(T_true - 0.02 + 0.005 * k, 4) for k in range(9)]
    cube = _synthetic_cube(Ls, T_grid, T_true, 2.0, n_seeds=8, noise=1.5e-3,
                           seed=99)
    est = make_sandvik_pair_estimator(8, 16, T_grid)
    out = jackknife_tbkt_over_seeds(cube, est, n_seeds=8)
    assert out["jack_mean"] is not None
    assert abs(out["jack_mean"] - T_true) < 0.01, out


# ---------------------------------------------------------------------------
# (c) Estimator-Verdrahtung (WM-Fit ueber means_LT-Schnittstelle)
# ---------------------------------------------------------------------------

def test_wm_fit_estimator_recovers_T_true():
    """make_wm_fit_estimator rekonstruiert T_true aus WM-konformen means_LT
    (Methode A, alle L)."""
    T_true, C_true = 0.576, 2.0
    Ls = [12, 24, 48]
    T_grid = [round(T_true - 0.02 + 0.005 * k, 4) for k in range(9)]
    means_LT = {L: {T: wm_form(T_true, L, C_true)
                    + 0.3 * (T - T_true) * math.log(L) for T in T_grid}
                for L in Ls}
    est = make_wm_fit_estimator(Ls, T_grid)
    tb = est(means_LT)
    assert tb is not None, "WM-Fit liefert keinen inneren Schaetzer"
    assert abs(tb - T_true) < 0.006, tb


def test_wm_fit_estimator_rejects_edge_minimum():
    """Liegt das argmin am T-Gitter-Rand, MUSS der Estimator None liefern
    (kein verlaesslicher Punktschaetzer am Rand) - Silent-Failure-Gate.

    Haertung 2026-07-10 (Code-Audit 4.2): die fruehere Assertion
    `None or im Inneren` liess fast jedes Verhalten durch (auch einen
    Estimator ohne Rand-Erkennung, der vom Rand ins Innere extrapoliert).
    Konstruktion hier: WM-konforme Daten mit T_true OBERHALB des Gitters
    -> chi^2(T) faellt monoton zum rechten Rand -> argmin garantiert am
    Rand -> None ist die einzig richtige Antwort.
    """
    Ls = [12, 24, 48]
    T_grid = [0.55, 0.56, 0.57, 0.58]
    T_true, C = 0.70, 2.0   # weit rechts ausserhalb des Gitters
    means_LT = {L: {T: wm_form(T_true, L, C)
                    + 0.3 * (T - T_true) * math.log(L) for T in T_grid}
                for L in Ls}
    est = make_wm_fit_estimator(Ls, T_grid)
    assert est(means_LT) is None
    # Spiegelfall: T_true unterhalb des Gitters -> argmin am linken Rand.
    T_true_lo = 0.45
    means_lo = {L: {T: wm_form(T_true_lo, L, C)
                    + 0.3 * (T - T_true_lo) * math.log(L) for T in T_grid}
                for L in Ls}
    assert est(means_lo) is None


def test_both_reference_values_are_distinct_and_documented():
    """Beide Literaturwerte sind vorhanden und verschieden (kein Kanon)."""
    assert abs(REF_HONEYCOMB_MULTI - 0.573) < 1e-9
    assert abs(REF_HONEYCOMB_DEDICATED - 0.576) < 1e-9
    assert REF_HONEYCOMB_MULTI != REF_HONEYCOMB_DEDICATED


# ---------------------------------------------------------------------------
# (d) Reproduzierbarkeit gegen PHY031 v01 (slow: echtes Mini-Wolff-Sampling)
# ---------------------------------------------------------------------------

# De-slow 2026-07-10 (Code-Audit 1.3): gemessen <3.5s - gehoert als
# fails-before-Gate in die schnelle CI-Suite, nicht hinter -m slow.
def test_seedwise_mean_matches_phy031_v01():
    """measure_honeycomb_seedwise teilt den RNG-Stream-Vertrag von PHY031 v01:
    bei identischen Parametern ist das per-Seed-Mittel bit-fuer-bit gleich dem
    PHY031-v01-Mittel (kein stiller Sampling-Drift)."""
    from phy031_honeycomb import measure_honeycomb

    L, T = 12, 0.62
    a = measure_honeycomb(L, T, n_seeds=3, n_measure=40, n_burn=40,
                          master_seed=42)
    b = measure_honeycomb_seedwise(L, T, n_seeds=3, n_measure=40, n_burn=40,
                                   master_seed=42)
    assert abs(a.upsilon_mean - b.upsilon_mean) < 1e-12, (
        a.upsilon_mean, b.upsilon_mean)
    assert len(b.upsilon_per_seed) == 3
    assert abs(float(np.mean(b.upsilon_per_seed)) - a.upsilon_mean) < 1e-12
