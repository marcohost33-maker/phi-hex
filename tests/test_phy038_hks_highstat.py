"""PHY038 Hochstatistik-HKS: Konfig-/Integritaets-Gates (instant, ohne MC).

Der teure Wolff-Lauf (L<=64, 16 Seeds, ~32 min) ist KEIN CI-Test; PHY038 nutzt
die in test_phy035/036 gegateten Bausteine. Hier nur die billig pruefbare
Hochstatistik-Konfiguration + das Geometrie-Orakel.
"""
from __future__ import annotations

from phy038_hks_highstat import (
    DOUBLING_PAIRS,
    LADDER,
    N_BURN,
    N_MEASURE,
    N_SEEDS,
    SW_LIMIT_HONEYCOMB,
    aligned_upsilon_zero,
)


def test_highstat_config_exceeds_phy036():
    """Hochstatistik: mehr Seeds + Sweeps als PHY036 (8/140/60)."""
    assert N_SEEDS >= 16 and N_MEASURE >= 240 and N_BURN >= 120
    assert LADDER == (16, 24, 32, 48, 64)
    assert DOUBLING_PAIRS == ((16, 32), (24, 48), (32, 64))


def test_aligned_upsilon_zero_exact():
    """Gate A: per-Site Υ(T->0)=0.75 J exakt."""
    assert abs(aligned_upsilon_zero() - SW_LIMIT_HONEYCOMB) < 1e-6
