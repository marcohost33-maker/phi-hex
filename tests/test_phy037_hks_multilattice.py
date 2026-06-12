"""PHY037 gitter-uebergreifende HKS-Extrapolation: Korrektheits-Gates (analytisch).

(a) DRIFT-GUARD: die aus den committed Reports geparsten Paar-Schaetzer MUESSEN
    den dokumentierten Werten entsprechen (Reality-Anchor; kein Report-Drift).
(b) EXTRAPOLATIONS-REGRESSION: die L->inf-Werte je Gitter sind festgepinnt.
(c) METHODEN-VERDIKT: honeycomb/kagome SAUBER, triangular UEBERSCHIESST -
    pinnt den ehrlichen Befund (kein stilles Schoenen).
(d) WIRING: OLS-Extrapolation rekonstruiert eine bekannte Gerade exakt.
"""
from __future__ import annotations

import math

from phy037_hks_multilattice import (
    LATTICES,
    analyse_lattice,
    extrapolate,
    parse_pairs,
    run_phy037,
)


# ---------------------------------------------------------------------------
# (a) Drift-Guard
# ---------------------------------------------------------------------------

def test_parsed_pairs_match_committed_reports():
    for name, cfg in LATTICES.items():
        pairs = parse_pairs(cfg["report"])
        for k, exp in cfg["expected"].items():
            assert k in pairs, (name, k)
            assert abs(pairs[k] - exp) < 1e-4, (name, k, pairs[k], exp)


# ---------------------------------------------------------------------------
# (b) Extrapolations-Regression
# ---------------------------------------------------------------------------

def test_extrapolation_values_pinned():
    exp = {"honeycomb": 0.5747, "kagome": 0.8236, "triangular": 1.4591}
    for name, want in exp.items():
        got = analyse_lattice(name)["primary"]
        assert abs(got - want) < 1.5e-3, (name, got, want)


# ---------------------------------------------------------------------------
# (c) Methoden-Verdikt (ehrlicher Befund)
# ---------------------------------------------------------------------------

def test_verdicts_honeycomb_kagome_clean_triangular_overshoots():
    assert analyse_lattice("honeycomb")["clean"] is True
    assert analyse_lattice("kagome")["clean"] is True
    # triangular schiesst ueber 1.418 und ist NICHT naeher als das blanke Paar
    tri = analyse_lattice("triangular")
    assert tri["clean"] is False
    assert tri["dev_extrap"] > 0  # ueberschiesst (positiv)
    assert abs(tri["dev_bare"]) < abs(tri["dev_extrap"])  # Paar besser


def test_overall_gate_passes():
    assert run_phy037()["overall"] is True


# ---------------------------------------------------------------------------
# (d) Wiring
# ---------------------------------------------------------------------------

def test_extrapolate_recovers_known_line():
    """Paare exakt auf einer Geraden in u=1/(ln Lc)^2 -> Intercept exakt."""
    from phy037_hks_multilattice import char_size, log_variable
    a, b = 0.50, 0.7
    pairs = {}
    for (l1, l2) in [(8, 16), (12, 24), (16, 32)]:
        u = log_variable(char_size(l1, l2, "geom"), 2)
        pairs[(l1, l2)] = a + b * u
    got = extrapolate(pairs, power=2, size_mode="geom")
    assert got is not None and math.isfinite(got)
    assert abs(got - a) < 1e-9, got
