"""PHY035 dichte-Leiter-HKS-Extrapolation: Korrektheits-Gates.

PHY035 ist ein NEGATIV/NUANCIERTES Ergebnis: die echte >=3-Punkt-Extrapolation
auf L<=32 zeigt KEINE saubere Konvergenz zur Referenz (falsifiziert die
PHY034-2-Punkt-Lesart). Diese Tests sichern die ANALYSE-INTEGRITAET (nicht die
Physik-Hypothese):

(a) GEOMETRIE-ORAKEL (instant): Upsilon(T->0)=0.75 J exakt (per-Site, z=3).
(b) EXTRAPOLATIONS-WIRING (instant): synthetische reine-WM-Daten -> alle Paare
    treffen T_true -> die HKS-Geraden-Extrapolation gibt T_true exakt zurueck.
(c) PAAR-WIRING (instant): C-freie Paar-Schaetzung auf synthetischen Daten.
(d) SMOKE (slow): Mini-Wolff-Lauf -> Pipeline laeuft, extrap endlich, Gate A.
"""
from __future__ import annotations

import math

import pytest

from phy034_hks_extrapolation import wm_form
from phy035_hks_dense_ladder import (
    DOUBLING_PAIRS,
    LADDER,
    SW_LIMIT_HONEYCOMB,
    T_GRID,
    aligned_upsilon_zero,
    make_hks_extrap_estimator,
    pair_estimates,
    run_phy035,
)


# ---------------------------------------------------------------------------
# (a) Geometrie-Orakel
# ---------------------------------------------------------------------------

def test_aligned_upsilon_zero_is_three_quarters():
    """Per-Site-T=0-Spinwellen-Grenzwert des honeycomb-Torus = 0.75 J exakt."""
    assert abs(aligned_upsilon_zero(L=8) - SW_LIMIT_HONEYCOMB) < 1e-6


# ---------------------------------------------------------------------------
# (b) Extrapolations-Wiring (synthetische reine-WM-Daten)
# ---------------------------------------------------------------------------

def _synthetic_means(t_true: float, c: float, slope: float = 0.3):
    """Reine WM-Form + glatte T-Abhaengigkeit -> Paar-Bedingung trifft t_true."""
    return {
        L: {T: wm_form(t_true, L, c) + slope * (T - t_true) * math.log(L)
            for T in T_GRID}
        for L in LADDER
    }


def test_hks_estimator_recovers_synthetic_t_true():
    """Wenn alle Paare denselben t_true treffen, MUSS die Geraden-Extrapolation
    (flach) genau t_true zurueckgeben (Wiring-/Vorzeichen-Guard)."""
    t_true = 0.60
    means = _synthetic_means(t_true, c=1.5)
    est = make_hks_extrap_estimator()  # p=2, geom
    got = est(means)
    assert got is not None
    assert abs(got - t_true) < 5e-3, got


def test_pair_estimates_finite_on_synthetic():
    """Alle drei Verdopplungspaare finden den Nulldurchgang (kein None)."""
    means = _synthetic_means(0.60, c=1.0)
    pairs = pair_estimates(means)
    assert set(pairs) == set(DOUBLING_PAIRS)
    for k, v in pairs.items():
        assert v is not None and math.isfinite(v), k
        assert abs(v - 0.60) < 8e-3, (k, v)


# ---------------------------------------------------------------------------
# (d) Smoke (slow): Mini-Wolff-Lauf
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_dense_ladder_smoke_runs_and_is_integer():
    """Mini-Lauf: Pipeline laeuft, Gate A exakt, Extrapolation endlich/plausibel.
    KEINE Aussage zum Punktschaetzer (Mini-Statistik) - nur Integritaet."""
    rep = run_phy035(n_seeds=4, n_measure=50, n_burn=25)
    assert abs(rep["ups0"] - SW_LIMIT_HONEYCOMB) < 1e-6
    ex = rep["extrap"]
    assert ex is not None and math.isfinite(ex)
    assert 0.50 < ex < 0.70, ex
    assert rep["gates"]["PASS_ALIGNED_EXACT_THREE_QUARTERS"]
    assert rep["gates"]["PASS_THREE_PAIRS_FINITE"]
