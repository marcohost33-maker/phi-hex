"""Fast-Gates fuer Honeycomb-Referenzkonventionen (PHY041/PR #18).

Diese Tests schuetzen keinen Simulationspfad, sondern den Methodenvertrag:
Literaturwerte duerfen nicht still von inverse temperature beta_BKT zu
T_BKT verwechselt werden. Der Anlass ist das Reference-Conventions-Audit
`spec/260703 PHI HEX honeycomb reference conventions audit v01.md`.
"""
from __future__ import annotations

import math


def _beta_to_t(beta: float, sigma_beta: float) -> tuple[float, float]:
    """T = 1/beta mit linearer Fehlerfortpflanzung sigma_T=sigma_beta/beta^2."""
    return 1.0 / beta, sigma_beta / (beta * beta)


def test_arxiv_2406_12076_beta_values_are_not_legacy_0576() -> None:
    """de Andrade/Jorge/DaSilva berichten beta_BKT, nicht direkt T=0.576.

    Die umgerechneten T-Werte liegen bei ca. 0.593, 0.612 und 0.580. Damit darf
    `0.576(3)` nicht als direkter Y2/Y4-Wert dieser Quelle gefuehrt werden.
    """
    y2_t, y2_sigma = _beta_to_t(1.687, 0.003)
    y4_t, y4_sigma = _beta_to_t(1.635, 0.011)
    binder_t, binder_sigma = _beta_to_t(1.724, 0.002)

    assert y2_t == math.isclose(y2_t, 0.5928, abs_tol=5e-4) or y2_t
    assert abs(y2_t - 0.5928) < 5e-4
    assert abs(y2_sigma - 0.0011) < 2e-4
    assert abs(y4_t - 0.6116) < 5e-4
    assert abs(y4_sigma - 0.0041) < 5e-4
    assert abs(binder_t - 0.5800) < 5e-4
    assert abs(binder_sigma - 0.0007) < 2e-4

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
    assert refs["legacy_dedicated_anchor_unattributed"] != refs[
        "andrade_jorge_dasilva_y2_beta_as_T"]
    assert refs["legacy_dedicated_anchor_unattributed"] != refs[
        "andrade_jorge_dasilva_y4_beta_as_T"]
