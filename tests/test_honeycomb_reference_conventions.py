"""Fast-Gates fuer Honeycomb-Referenzkonventionen (PHY041/PR #18).

Diese Tests schuetzen keinen Simulationspfad, sondern den Methodenvertrag:
Literaturwerte duerfen nicht still von inverse temperature beta_BKT zu
T_BKT verwechselt werden. Der Anlass ist das Reference-Conventions-Audit
`spec/260703 PHI HEX honeycomb reference conventions audit v01.md`.

Haertung 2026-07-10 (Code-Audit 4.1): die Tests pruefen jetzt PRODUKTIONS-
Konstanten (phy042.REF_BAND, PHY041-REFs) statt nur test-lokaler Literale -
vorher konnte der Drift, den sie verhindern sollen (0.576 wieder als
direkter Anker im Code), gar nicht auffallen.
"""
from __future__ import annotations

import math

import conftest

phy041 = conftest._load(
    "phy041_honeycomb_wl",
    "260702 PHY041 honeycomb wang-landau entropic helicity v01.py")
phy042 = conftest._load(
    "phy042_honeycomb_wl_fss", "260706 PHY042 honeycomb wl fss v01.py")


def _beta_to_t(beta: float, sigma_beta: float) -> tuple[float, float]:
    """T = 1/beta mit linearer Fehlerfortpflanzung sigma_T=sigma_beta/beta^2."""
    return 1.0 / beta, sigma_beta / (beta * beta)


def test_current_reference_ledger_is_versioned_and_direct_T():
    """Issue #45: aktuelle Produktionsreferenzen stammen aus den derzeitigen
    Quellenfassungen; supersedierte beta-v1-Werte sind keine Live-Anker.
    """
    band = phy042.REF_BAND
    assert band["okabe_otsuka_2501_v1_rough"] == (0.573, None)
    assert band["jiang_ptep_nn"] == (0.572, 0.003)
    assert band["jiang_ptep_helicity"] == (0.576, 0.004)
    assert band["andrade_v4_helicity_sa"] == (0.575, 0.008)
    assert band["andrade_v4_helicity_wl"] == (0.576, 0.003)
    assert band["andrade_v4_upsilon4_sa"] == (0.551, 0.011)
    assert band["andrade_v4_upsilon4_wl"] == (0.568, 0.001)


def test_superseded_beta_values_are_not_live_reference_keys():
    """Die alte v1-Konversion bleibt reproduzierbare Historie, darf aber nicht
    mehr im aktuellen REF_BAND als beta-Kanal erscheinen.
    """
    old = {
        "upsilon_beta": _beta_to_t(1.687, 0.003),
        "upsilon4_beta": _beta_to_t(1.635, 0.011),
        "binder_beta": _beta_to_t(1.724, 0.002),
    }
    assert math.isclose(old["upsilon_beta"][0], 0.5928, abs_tol=5e-4)
    assert math.isclose(old["upsilon4_beta"][0], 0.6116, abs_tol=5e-4)
    assert math.isclose(old["binder_beta"][0], 0.5800, abs_tol=5e-4)
    for key in old:
        assert key not in phy042.REF_BAND


def test_primary_w4_reference_set_excludes_upsilon4_diagnostic():
    """Y4 bleibt Diagnostik: die aktuelle Quelle selbst warnt vor der
    FSS-/Minimum-Systematik. W4-Primärvergleich nutzt direkte T_BKT-
    Schätzungen aus Korrelations-/NN- und Y2-Helicity-Kanaelen.
    """
    assert set(phy042.REF_PRIMARY_KEYS) == {
        "okabe_otsuka_2501_v1_rough",
        "jiang_ptep_nn",
        "jiang_ptep_helicity",
        "andrade_v4_helicity_sa",
        "andrade_v4_helicity_wl",
    }
    assert all("upsilon4" not in k for k in phy042.REF_PRIMARY_KEYS)


def test_historical_phy041_anchor_remains_lineage_only():
    """PHY041 ist ein datierter historischer Lauf; sein 0.576-Anker bleibt
    fuer Reproduktion erhalten, wird aber nicht mit dem neuen Ledger verwechselt.
    """
    assert phy041.REF_MULTI == 0.573
    assert phy041.REF_DEDIC == 0.576
    assert phy042.REF_BAND["andrade_v4_helicity_wl"] == (0.576, 0.003)

