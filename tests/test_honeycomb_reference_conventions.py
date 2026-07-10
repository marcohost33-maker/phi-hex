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


def test_production_ref_band_uses_converted_beta_channels():
    """phy042.REF_BAND fuehrt die beta-Kanaele in korrekt konvertierter
    T-Form (kein stiller beta/T-Mix im Produktionscode)."""
    band = phy042.REF_BAND
    for key, (beta, sig_beta) in (("upsilon_beta", (1.687, 0.003)),
                                  ("upsilon4_beta", (1.635, 0.011)),
                                  ("binder_beta", (1.724, 0.002))):
        t, sig_t = _beta_to_t(beta, sig_beta)
        val, sig = band[key]
        assert math.isclose(val, t, abs_tol=5e-4), key
        assert math.isclose(sig, sig_t, abs_tol=2e-4), key


def test_production_refs_keep_legacy_0576_out_of_band_channels():
    """Kein REF_BAND-Kanal darf der Legacy-Direktlesart 0.576 entsprechen;
    PHY041 fuehrt 0.576 nur als explizit benannten 'dedicated'-Anker."""
    for key, v in phy042.REF_BAND.items():
        val = v[0] if isinstance(v, tuple) else v
        if key == "multi_lattice":
            continue
        assert not math.isclose(val, 0.576, abs_tol=1e-9), key
    assert phy041.REF_MULTI == 0.573
    assert phy041.REF_DEDIC == 0.576


def test_arxiv_2406_12076_beta_values_are_not_legacy_0576() -> None:
    """de Andrade/Jorge/DaSilva berichten beta_BKT, nicht direkt T=0.576.

    Die umgerechneten T-Werte liegen bei ca. 0.593, 0.612 und 0.580. Damit darf
    `0.576(3)` nicht als direkter Y2/Y4-Wert dieser Quelle gefuehrt werden.
    """
    y2_t, y2_sigma = _beta_to_t(1.687, 0.003)
    y4_t, y4_sigma = _beta_to_t(1.635, 0.011)
    binder_t, binder_sigma = _beta_to_t(1.724, 0.002)

    assert math.isclose(y2_t, 0.5928, abs_tol=5e-4)
    assert math.isclose(y2_sigma, 0.0011, abs_tol=2e-4)
    assert math.isclose(y4_t, 0.6116, abs_tol=5e-4)
    assert math.isclose(y4_sigma, 0.0041, abs_tol=5e-4)
    assert math.isclose(binder_t, 0.5800, abs_tol=5e-4)
    assert math.isclose(binder_sigma, 0.0007, abs_tol=2e-4)

    for value in (y2_t, y4_t, binder_t):
        assert abs(value - 0.576) > 0.003


def test_honeycomb_reference_band_keeps_sources_separate() -> None:
    """Die Repo-Vergleiche behalten Multi-Lattice, Jiang und beta-Werte getrennt."""
    refs = {
        "okabe_otsuka_multi_lattice": 0.573,
        "jiang_helicity_direct_T": 0.571,
        "jiang_nn_direct_T": 0.560,
        "andrade_jorge_dasilva_y2_beta_as_T": _beta_to_t(1.687, 0.003)[0],
        "andrade_jorge_dasilva_y4_beta_as_T": _beta_to_t(1.635, 0.011)[0],
        "andrade_jorge_dasilva_binder_beta_as_T": _beta_to_t(1.724, 0.002)[0],
        "legacy_dedicated_anchor_unattributed": 0.576,
    }

    assert refs["jiang_helicity_direct_T"] < refs["andrade_jorge_dasilva_y2_beta_as_T"]
    assert refs["okabe_otsuka_multi_lattice"] < refs["andrade_jorge_dasilva_y2_beta_as_T"]
    assert not math.isclose(
        refs["legacy_dedicated_anchor_unattributed"],
        refs["andrade_jorge_dasilva_y2_beta_as_T"],
        abs_tol=0.003,
    )
    assert not math.isclose(
        refs["legacy_dedicated_anchor_unattributed"],
        refs["andrade_jorge_dasilva_y4_beta_as_T"],
        abs_tol=0.003,
    )
