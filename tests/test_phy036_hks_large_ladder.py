"""PHY036 grosse-L-Leiter: Konfig-/Integritaets-Gates (instant, ohne MC).

Der teure Wolff-Lauf (L bis 64, ~Minuten) ist KEIN CI-Test; PHY036 nutzt die
in test_phy035 bereits gegateten Bausteine (measure_cube/pair_estimates/
make_hks_extrap_estimator) + die PHY032-Bootstrap-Maschinerie. Hier nur die
billig pruefbare Konfiguration + das Geometrie-Orakel.
"""
from __future__ import annotations

from phy036_hks_large_ladder import (
    DOUBLING_PAIRS,
    LADDER,
    SW_LIMIT_HONEYCOMB,
    aligned_upsilon_zero,
)


def test_ladder_has_three_nested_doubling_pairs():
    """L=16/24/32/48/64 -> (16,32)/(24,48)/(32,64), groesstes Paar (32,64)."""
    assert LADDER == (16, 24, 32, 48, 64)
    assert DOUBLING_PAIRS == ((16, 32), (24, 48), (32, 64))
    for (l1, l2) in DOUBLING_PAIRS:
        assert l2 == 2 * l1
        assert l1 in LADDER and l2 in LADDER


def test_aligned_upsilon_zero_exact():
    """Gate A: per-Site Υ(T->0)=0.75 J exakt (Geometrie/Konvention korrekt)."""
    assert abs(aligned_upsilon_zero() - SW_LIMIT_HONEYCOMB) < 1e-6
