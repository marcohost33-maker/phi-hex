"""Regression issue #62 (Folge zu PR #64): NaN-feste min/max-Gates in phy017
und PHY029.

Invariante: Ein PASS-Gate darf nur bei einem WAHREN Vergleich auf endlichen
Werten PASS liefern; ein NaN in irgendeiner Position des Aggregats muss FAIL
ergeben. Python ``max()``/``min()`` ueber eine Liste mit NaN ist
reihenfolgeabhaengig: ``max([nan, x])`` liefert nan, ``max([x, nan])`` liefert
x -- das NaN verschwindet. Darum je Gate: NaN an Position 0, 1 und zuletzt,
+/-inf, Negativkontrolle (endliche Werte -> altes Ergebnis) und Grenze.

Alle Faelle laufen ueber die ECHTE Gate-Zeile (run_module_3_madelung_fisher
bzw. run_phy029), die Produzenten sind per monkeypatch ersetzt -- kein
Wolff-Sampling, keine Produktionslaeufe.
"""
from __future__ import annotations

import math

import pytest

import conftest

phy017 = conftest.load_all()["phy017_validation_kernel"]
phy029 = conftest.load_all()["phy029_triangular_sandvik"]

NAN = math.nan
INF = math.inf

# ---------------------------------------------------------------------------
# phy017 Modul 3: PASS_LINK_IDENTITY = max(link_err) < 1e-10
# ---------------------------------------------------------------------------

_N_LIST_LEN = 6  # N_list = [6, 12, 24, 48, 96, 192] in run_module_3


def _run_module_3_with_link_errors(monkeypatch, errs):
    """Faehrt run_module_3_madelung_fisher mit vorgegebenen link_err je N.

    link_err = abs(8 * <Q> - I_F). I_F wird auf 0 gesetzt und <Q> = e / 8
    geliefert; Division und Multiplikation mit 8 sind in IEEE-754 exakt,
    also ist link_err == abs(e) bitgenau.
    """
    assert len(errs) == _N_LIST_LEN
    it = iter(errs)
    monkeypatch.setattr(phy017, "fisher_information_link_form",
                        lambda rho, a: 0.0)
    monkeypatch.setattr(phy017, "expectation_value_Q",
                        lambda Q, rho, a: next(it) / 8.0)
    out = phy017.run_module_3_madelung_fisher()
    assert next(it, None) is None, "nicht alle Fehlerwerte verbraucht"
    return out


def _with(value, pos, base=1e-15):
    errs = [base] * _N_LIST_LEN
    errs[pos] = value
    return errs


@pytest.mark.parametrize("pos", [0, 1, 3, _N_LIST_LEN - 1])
def test_phy017_link_identity_nan_any_position_fails(monkeypatch, pos):
    out = _run_module_3_with_link_errors(monkeypatch, _with(NAN, pos))
    assert out["PASS_LINK_IDENTITY"] is False, (pos, out["max_link_identity_error"])
    # der NaN darf im Aggregat nicht verschwinden
    assert not math.isfinite(out["max_link_identity_error"])


@pytest.mark.parametrize("pos", [0, 1, _N_LIST_LEN - 1])
def test_phy017_link_identity_inf_fails(monkeypatch, pos):
    out = _run_module_3_with_link_errors(monkeypatch, _with(INF, pos))
    assert out["PASS_LINK_IDENTITY"] is False
    assert out["max_link_identity_error"] == INF


def test_phy017_link_identity_finite_negative_control(monkeypatch):
    """Endliche Werte: identisch zum alten builtin-max-Ergebnis."""
    errs = [3e-14, 7e-12, 1e-15, 5e-11, 2e-13, 4e-16]
    out = _run_module_3_with_link_errors(monkeypatch, errs)
    assert out["max_link_identity_error"] == max(errs)
    assert out["PASS_LINK_IDENTITY"] is True


def test_phy017_link_identity_boundary(monkeypatch):
    """Strikt kleiner: genau 1e-10 ist FAIL, der Nachbar darunter PASS."""
    at = _run_module_3_with_link_errors(monkeypatch, _with(1e-10, 2))
    assert at["max_link_identity_error"] == 1e-10
    assert at["PASS_LINK_IDENTITY"] is False
    below = math.nextafter(1e-10, 0.0)
    ok = _run_module_3_with_link_errors(monkeypatch, _with(below, 2))
    assert ok["max_link_identity_error"] == below
    assert ok["PASS_LINK_IDENTITY"] is True


def test_phy017_link_identity_unpatched_still_passes():
    """Echte Produzenten (Gauss-Ring, N=6..192): Gate bleibt PASS und das
    Maximum ist das der Einzelfehler."""
    out = phy017.run_module_3_madelung_fisher()
    errs = [float(e) for _, e in out["link_errors_by_N"]]
    assert all(math.isfinite(e) for e in errs)
    assert out["max_link_identity_error"] == max(errs)
    assert out["PASS_LINK_IDENTITY"] is True


# ---------------------------------------------------------------------------
# PHY029: PASS_BRACKETS_REFERENCE = min(valid) <= ref <= max(valid)
# ---------------------------------------------------------------------------

REF = 1.418  # run_phy029: ref = 1.418


def _run_phy029_with_pairs(monkeypatch, values):
    """run_phy029 mit vorgegebenen Paar-Ergebnissen (pairs (4,8), (6,12))."""
    table = {(4, 8): values[0], (6, 12): values[1]}
    monkeypatch.setattr(phy029, "sandvik_pair",
                        lambda r1, r2, temps, cache, cfg: table[(r1, r2)])
    return phy029.run_phy029()


@pytest.mark.parametrize("values", [
    (NAN, REF),   # Position 0
    (REF, NAN),   # Position 1 = letzte: frueher min = max = REF -> PASS
    (1.40, INF),  # inf als obere Klammer klammert alles
    (-INF, 1.43),
    (NAN, NAN),
])
def test_phy029_brackets_reference_nonfinite_fails(monkeypatch, values):
    out = _run_phy029_with_pairs(monkeypatch, values)
    gates = out["pass_gates"]
    assert gates["PASS_BRACKETS_REFERENCE"] is False, values
    assert gates["PASS_ALL_WITHIN_2PCT"] is False, values
    assert out["overall_pass"] is False


@pytest.mark.parametrize("valid", [
    [NAN, 1.40, 1.43],
    [1.40, NAN, 1.43],   # Mitte: frueher min 1.40, max 1.43 -> PASS
    [1.40, 1.43, NAN],   # letzte: frueher PASS
    [1.40, 1.43, INF],
    [-INF, 1.40, 1.43],
])
def test_phy029_brackets_helper_any_position_fails(valid):
    assert phy029._brackets_reference(valid, REF) is False


def test_phy029_brackets_finite_negative_control(monkeypatch):
    """Endliche Werte: identisch zum alten Ausdruck."""
    for values, expected in (((1.40, 1.43), True), ((1.43, 1.40), True),
                             ((1.42, 1.43), False), ((1.40, 1.41), False)):
        out = _run_phy029_with_pairs(monkeypatch, values)
        old = min(values) <= REF <= max(values)
        assert old is expected
        assert out["pass_gates"]["PASS_BRACKETS_REFERENCE"] is expected, values


def test_phy029_brackets_boundary_inclusive():
    """<= ist inklusiv: genau REF klammert, der Nachbar darueber nicht."""
    assert phy029._brackets_reference([REF, 1.43], REF) is True
    assert phy029._brackets_reference([1.40, REF], REF) is True
    above = math.nextafter(REF, INF)
    assert phy029._brackets_reference([above, 1.43], REF) is False


def test_phy029_within_2pct_not_vacuous(monkeypatch):
    """Kein Paar konvergiert: all() ueber [] waere True = PASS ohne einen
    einzigen Vergleich. Ein konvergiertes Paar innerhalb 2 % bleibt PASS."""
    out = _run_phy029_with_pairs(monkeypatch, (None, None))
    assert out["pass_gates"]["PASS_ALL_WITHIN_2PCT"] is False
    one = _run_phy029_with_pairs(monkeypatch, (1.42, None))
    assert one["pass_gates"]["PASS_ALL_WITHIN_2PCT"] is True


def test_phy029_brackets_needs_two_values():
    assert phy029._brackets_reference([], REF) is False
    assert phy029._brackets_reference([REF], REF) is False
