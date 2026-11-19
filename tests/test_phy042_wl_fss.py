"""Schnelle Korrektheits-Gates fuer PHY042 (honeycomb WL-FSS L=24/32/48).

Die teuren WL-Laeufe (L=24/32/48, Multi-Walker) laufen als slow/lokal
(Gate-Log in results/). Die CI-schnellen Gates pruefen MC-frei: den
Produktions-Skalierungs-Vertrag, den RNG-Stream-Vertrag (Kollisionfreiheit),
das identische T-Gitter (PHY041-Bruecke), die deterministische Bindung der
PHY041-Brueckenkonstanten an den committed Report (Muster PHY037-Drift-Guard)
und die Referenzband-Konsistenz mit der Konventions-Audit-Spec.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

import conftest

phy042 = conftest._load(
    "phy042_honeycomb_wl_fss", "260706 PHY042 honeycomb wl fss v01.py")

ROOT = Path(__file__).resolve().parents[1]
PHY041_REPORT = (ROOT / "results" /
                 "260702 PHY041 honeycomb wang-landau entropic helicity report.txt")


def test_prod_sweeps_contract():
    """L=24 exakt PHY041 (30000, deterministische Bruecke); sonst
    max(30000, 60*nbins) gegen Bin-Ausduennung bei L=32/48."""
    assert phy042._prod_sweeps_for(512, 24) == 30000
    assert phy042._prod_sweeps_for(5000, 24) == 30000
    assert phy042._prod_sweeps_for(400, 32) == 30000       # Floor greift
    assert phy042._prod_sweeps_for(900, 32) == 54000       # 60/Bin greift
    assert phy042._prod_sweeps_for(2000, 48) == 120000


def test_stream_contract_no_collision():
    """Walker 0 nutzt den PHY041-Standard-Stream (650, Bruecke); Zusatz-
    Walker liegen bei 90000+1000*w und kollidieren fuer die verwendeten
    L in {24,32,48} mit KEINEM bisherigen Stream-Vertrag (PHY040: 1/500+s/
    700+s; PHY041: 600+s+1000L, 650+L, 660+s+1000L)."""
    assert phy042._stream_for(0) == phy042._p41._STREAM_WL == 650
    assert phy042._stream_for(1) == 91000
    assert phy042._stream_for(2) == 92000
    used = set()
    for L in (24, 32, 48):
        for s in range(8):
            used |= {1, 500 + s, 700 + s, 600 + s + 1000 * L, 650 + L,
                     660 + s + 1000 * L}
    for w in (1, 2, 3):
        for L in (24, 32, 48):
            assert phy042._stream_for(w) + L not in used, (w, L)


def test_t_grid_identical_to_phy041():
    """Identisches T-Gitter wie PHY041 (0.52..0.67 Schritt 0.005) - Vertrag
    fuer die bitgleiche L=24-Bruecke und den VAL-C-Dip-Vergleich."""
    expected = np.round(np.arange(0.52, 0.6701, 0.005), 4)
    assert np.array_equal(phy042._T_GRID, expected)
    assert phy042._T_GRID[0] == 0.52 and phy042._T_GRID[-1] == 0.67


def test_phy041_bridge_constants_match_committed_report():
    """Die im Modul gepinnten PHY041-Brueckenwerte (Paar-Tabelle + Y4-Dip
    L=24) muessen zeichengenau aus dem committed PHY041-Report folgen
    (single source of truth, Muster PHY037/parse_phy032_grid)."""
    text = PHY041_REPORT.read_text(encoding="utf-8")
    pairs = {}
    for m in re.finditer(r"\|\s*\((\d+),(\d+)\)\s*\|\s*([0-9.]+)\s*\|", text):
        pairs[(int(m.group(1)), int(m.group(2)))] = float(m.group(3))
    assert pairs == phy042.PHY041_PAIRS
    m = re.search(r"L=24: Dip bei T=([0-9.]+)", text)
    assert m is not None
    assert math.isclose(float(m.group(1)), phy042.PHY041_Y4_DIP_L24)


def test_reference_band_matches_conventions_audit_spec():
    """Aktueller Evidenz-Ledger muss mit der Vertragsquelle uebereinstimmen.

    Issue #45: keine stillen arXiv-v1/beta-Altwerte mehr als aktuelle
    Vergleichswerte; jede aktuelle Quelle ist versioniert benannt.
    """
    spec = (ROOT / "spec" /
            "260703 PHI HEX honeycomb reference conventions audit v01.md"
            ).read_text(encoding="utf-8")
    band = phy042.REF_BAND

    expected = {
        "okabe_otsuka_2501_v1_rough": (0.573, None),
        "jiang_ptep_nn": (0.572, 0.003),
        "jiang_ptep_helicity": (0.576, 0.004),
        "andrade_v4_helicity_sa": (0.575, 0.008),
        "andrade_v4_helicity_wl": (0.576, 0.003),
        "andrade_v4_upsilon4_sa": (0.551, 0.011),
        "andrade_v4_upsilon4_wl": (0.568, 0.001),
    }
    assert band == expected
    for key in expected:
        assert key in spec
    for token in ("0.573", "0.572 +/- 0.003", "0.576 +/- 0.004",
                  "0.575 +/- 0.008", "0.576 +/- 0.003",
                  "0.551 +/- 0.011", "0.568 +/- 0.001"):
        assert token in spec
    assert "superseded" in spec.lower()


def test_w4_walker_plan_requires_three_per_l():
    """W4 darf eine Spread-Domaene nicht aus 1-2 Walkern ableiten."""
    phy042._w4_validate_walker_plan({64: 3, 96: 3, 128: 4})
    with pytest.raises(ValueError, match="mindestens 3 Walker"):
        phy042._w4_validate_walker_plan({64: 3, 96: 2, 128: 3})
    # Der Schwellwert selbst ist Teil des W4-Vertrags und darf nicht durch
    # einen Caller auf Kalibrierungsniveau abgesenkt werden.
    with pytest.raises(ValueError, match="nicht unter 3"):
        phy042._w4_validate_walker_plan({64: 2}, min_walkers=2)


def test_pair_inside_domains_checks_lower_and_upper_bounds():
    a = (0.56, 0.6475)  # z.B. PHY032-Fallback
    b = (0.52, 0.60)
    assert phy042._pair_inside_domains(0.575, a, b)
    assert not phy042._pair_inside_domains(0.555, a, b)  # unter Fallback-Evidenz
    assert not phy042._pair_inside_domains(0.61, a, b)
    assert not phy042._pair_inside_domains(0.575, None, b)
    assert not phy042._pair_inside_domains(0.575, a, None)
    assert not phy042._pair_inside_domains(None, a, b)

def test_pair_trend_inputs_are_sorted_ascending():
    """Die PHY041-Paar-Sequenz (Trend-Referenz) ist nach min-L sortiert und
    monoton steigend - das dokumentierte kleine-L-Negativ-Result."""
    seq = [phy042.PHY041_PAIRS[p]
           for p in sorted(phy042.PHY041_PAIRS, key=lambda p: (min(p), max(p)))]
    assert seq == sorted(seq)
    assert seq == [0.5951, 0.6029, 0.6087]


def test_validity_domain_contiguous_from_low_edge():
    """Domaene ist zusammenhaengend ab dem unteren T-Rand und endet an der
    ERSTEN Schwellen-Verletzung (auch wenn der Spread danach wieder faellt -
    ehrliche Lesart des Lauf-1-Befunds)."""
    spread = np.array([0.001, 0.01, 0.03, 0.05, 0.02, 0.01])
    dom = phy042._validity_domain(spread, 0.04)
    assert dom.tolist() == [True, True, True, False, False, False]
    assert phy042._validity_domain(spread, 0.0005).tolist() == [False] * 6
    assert phy042._validity_domain(spread, 1.0).all()


def test_uncovered_mass_two_level_oracle():
    """Coverage-Massen-Diagnostik auf einem 2-Bin-Boltzmann-Orakel: bei
    gleicher Entropie ist die Masse des unbesetzten Bins genau sein
    kanonisches Gewicht exp(-dE/T)/(1+exp(-dE/T))."""
    from types import SimpleNamespace
    res = SimpleNamespace(lng=np.array([0.0, 0.0]),
                          centers=np.array([-10.0, -9.0]),
                          mask=np.array([True, False]))
    T = 0.5
    expected = math.exp(-1.0 / T) / (1.0 + math.exp(-1.0 / T))
    assert phy042._uncovered_mass(res, T) == pytest.approx(expected, rel=1e-12)
    res_full = SimpleNamespace(lng=res.lng, centers=res.centers,
                               mask=np.array([True, True]))
    assert phy042._uncovered_mass(res_full, T) == 0.0


@pytest.mark.slow
def test_wl_job_smoke_small_lattice():
    """Slow-Smoke: ein kompletter (L=8, walker=1)-Job durch den PHY042-
    Job-Wrapper - volle Bin-Abdeckung, dichte Kurve, fallendes Y2."""
    m_lo, s_lo = phy042.wolff_anchor(8, 0.50)
    m_hi, s_hi = phy042.wolff_anchor(8, 0.70)
    e_lo, e_hi = phy042.window_from_anchors(m_lo, s_lo, m_hi, s_hi)
    L, w, res = phy042._wl_job((8, 1, e_lo, e_hi, 42))
    assert (L, w) == (8, 1)
    assert int(res.mask.sum()) == len(res.mask)          # alles abgedeckt
    for T in (0.55, 0.60, 0.65):
        assert phy042.canonical_edge_leak(res, T) < 1e-3
    c = phy042.upsilon_curves(res, phy042._T_GRID)
    assert c["y2"][0] > c["y2"][-1] + 0.1                # klar fallend


def test_partial_phy032_grid_never_installs_the_fallback_domain():
    """Teilgitter mit passenden Zeilen: grid_ok=False -> kein Fallback."""
    rows = [{"T": t, "ok": True} for t in phy042.PHY032_L24_GRID]
    assert phy042._phy032_fallback_bounds(rows, False) is None
    assert phy042._phy032_fallback_bounds(rows, True) == (0.56, 0.6475)
    assert phy042._phy032_fallback_bounds(rows[:-1] + [{"T": 0.65, "ok": False}], True) is None
    assert phy042._phy032_fallback_bounds([], True) is None


def test_fallback_needs_all_eight_rows_retained_after_the_window_filter():
    """grid_ok kam aus der rohen Parse-Zahl; faellt eine Zeile aus dem
    T-Fenster, duerfen sieben passende Zeilen keine Domaene tragen."""
    rows = [{"T": t, "ok": True} for t in phy042.PHY032_L24_GRID]
    assert phy042._phy032_fallback_bounds(rows[:7], True) is None
    assert phy042._phy032_fallback_bounds(rows[1:], True) is None


def test_fallback_needs_the_exact_contract_grid_not_just_eight_rows():
    """Eine Dublette anstelle eines fehlenden Punkts (Report-Drift) hat acht
    Zeilen mit denselben Grenzen -- traegt die Domaene trotzdem nicht."""
    rows = [{"T": t, "ok": True} for t in phy042.PHY032_L24_GRID]
    dup = rows[:3] + [rows[2]] + rows[4:]
    assert len(dup) == 8
    assert phy042._phy032_fallback_bounds(dup, True) is None
    assert phy042._phy032_fallback_bounds(rows[::-1], True) is None


def test_contract_grid_matches_the_committed_phy032_report():
    got = [round(r[0], 4) for r in phy042.parse_phy032_grid()[24]]
    assert tuple(got) == phy042.PHY032_L24_GRID


def test_pair_table_quotes_only_the_quotable_reason():
    """Issue #63 F-4/K-3: ``quotable`` wird in _assess_pairs aus dem Verdikt
    abgeleitet (frueher Zeile ``q = reason == "QUOTABLE"`` in run_phy042, ohne
    Test; der aufruferlose _pair_quotable ist entfernt)."""
    assert not hasattr(phy042, "_pair_quotable")
    t, ya, yb = _crossing_curves(0.647)
    # L=48 so, dass (32,48) bei 0.600 kreuzt; (24,48) folgt daraus bei 0.6235
    inv48 = 1.0 + 2.0 * math.log(48 / 32) - (t - 0.600)
    y48 = 2.0 * t * (1.0 + 1.0 / inv48) / math.pi
    curves = {24: {"y2": ya}, 32: {"y2": yb}, 48: {"y2": y48}}
    results = {(24, 0): "a", (32, 0): "b", (48, 0): "c"}
    clean = lambda res, T: 0.0  # noqa: E731
    leaky = lambda res, T: 0.5 if res == "c" else 0.0  # noqa: E731
    bounds = {24: (0.56, 0.65), 32: (0.52, 0.67), 48: (0.52, 0.67)}
    table = phy042._assess_pairs(t, curves, results, [24, 32, 48], bounds,
                                 leak_metric=leaky, unc_metric=clean)
    assert list(table) == [(24, 32), (24, 48), (32, 48)]
    assert phy042._size_pairs([24, 32, 48]) == list(table)
    tb, reason, q = table[(24, 32)]
    assert abs(tb - 0.647) < 1e-9 and reason == "QUOTABLE" and q is True
    tb, reason, q = table[(24, 48)]                    # Crossing vorhanden, aber Leck
    assert abs(tb - 0.6235) < 1e-9 and reason == "COVERAGE_LEAK" and q is False
    tb, reason, q = table[(32, 48)]
    assert abs(tb - 0.600) < 1e-9 and reason == "COVERAGE_LEAK" and q is False
    # ohne Leck sind alle drei Paare quotierbar (Positiv-Kontrolle)
    ok = phy042._assess_pairs(t, curves, results, [24, 32, 48], bounds,
                              leak_metric=clean, unc_metric=clean)
    assert [q for _tb, _r, q in ok.values()] == [True, True, True]
    # jede Zeile: quotable genau dann, wenn das Verdikt QUOTABLE ist
    assert all(q == (r == "QUOTABLE") for _tb, r, q in table.values())


def test_pair_reason_names_the_failed_gate():
    a, b = (0.56, 0.59), (0.56, 0.59)
    assert phy042._pair_reason(None, a, b, True, True, bracket_valid=True) == "NO_CROSSING"
    assert phy042._pair_reason(0.575, a, b, False, True, bracket_valid=True) == "COVERAGE_LEAK"
    assert phy042._pair_reason(0.575, a, b, True, False, bracket_valid=True) == "COVERAGE_UNCOVERED"
    assert phy042._pair_reason(0.61, a, b, True, True, bracket_valid=True) == "OUTSIDE_DOMAIN"
    assert phy042._pair_reason(0.575, a, b, True, True, bracket_valid=True) == "QUOTABLE"



def test_w4_walker_preflight_rejects_non_integral_values():
    nan = float("nan")
    for bad in ({64: nan}, {64: 2.5}, {64: True}):
        with pytest.raises(ValueError):
            phy042._w4_validate_walker_plan(bad)
    for thr in (nan, 3.5, True):
        with pytest.raises(ValueError):
            phy042._w4_validate_walker_plan({64: 1}, min_walkers=thr)


def _crossing_curves(t_root, La=24, Lb=32):
    """Synthetic y2 curves whose C-eliminated pair condition crosses at t_root."""
    import numpy as np
    t = phy042._T_GRID
    target = -2.0 * math.log(Lb / La)
    y2_b = 4.0 * t / math.pi                       # R_b = 1
    inv_ra = (target + 1.0) + (t - t_root)         # 1/R_a - 1/R_b - target = t - t_root
    y2_a = 2.0 * t * (1.0 + 1.0 / inv_ra) / math.pi
    return t, np.asarray(y2_a), np.asarray(y2_b)


def test_bracket_must_lie_in_the_common_domain():
    """W4-G2(4): root at 0.647 interpolated from 0.645/0.650; a domain that ends
    at 0.6475 leaves the upper bracket uncovered -> not quotable."""
    t, ya, yb = _crossing_curves(0.647)
    tb = phy042.tbkt_pair_from_curves(t, ya, 24, yb, 32)
    assert tb is not None and abs(tb - 0.647) < 1e-9
    dom = (0.56, 0.6475)
    assert phy042._pair_inside_domains(tb, dom, dom)           # point check alone passes
    assert not phy042._bracket_in_common_domain(t, ya, 24, yb, 32, tb, dom, dom)
    assert phy042._pair_reason(tb, dom, dom, True, True, bracket_valid=False) == "OUTSIDE_DOMAIN"
    # positive control: domain covering both bracket points
    wide = (0.56, 0.65)
    assert phy042._bracket_in_common_domain(t, ya, 24, yb, 32, tb, wide, wide)
    assert not phy042._bracket_in_common_domain(t, ya, 24, yb, 32, tb, None, wide)


def test_validity_domain_nan_spread_fails_closed():
    """Ein NaN-Spread ist nicht gemessen und beendet die Domaene."""
    import numpy as np
    spread = np.array([0.01, float("nan"), 0.02])
    assert phy042._validity_domain(spread, 0.04).tolist() == [True, False, False]
    # Negativ-Kontrolle: endliche Werte unter der Schwelle bleiben gueltig
    assert phy042._validity_domain(np.array([0.01, 0.02, 0.05]), 0.04).tolist() == [True, True, False]


def test_pair_coverage_is_scoped_to_the_pair_and_its_common_domain():
    """Ein fremdes L oder ein Punkt ausserhalb der Paar-Domaene darf das Paar
    nicht kippen; innerhalb zaehlt jeder Walker beider L."""
    t = [0.55, 0.56, 0.60, 0.65, 0.66]
    bad_at = {(48, 0): {0.60}, (24, 0): {0.55, 0.66}, (32, 1): set()}
    results = {k: k for k in [(24, 0), (32, 0), (32, 1), (48, 0)]}

    def metric(res, T):
        return 0.5 if T in bad_at.get(res, set()) else 0.0

    dom24, dom32 = (0.56, 0.6475), (0.52, 0.67)
    # L=48 ist fremd, L=24-Ausreisser liegen ausserhalb [0.56, 0.6475]
    assert phy042._pair_coverage_max(results, 24, 32, t, dom24, dom32, metric) == 0.0
    # Positiv-Kontrolle: Ausreisser innerhalb der gemeinsamen Domaene zaehlt
    bad_at[(32, 1)] = {0.60}
    assert phy042._pair_coverage_max(results, 24, 32, t, dom24, dom32, metric) == 0.5
    assert phy042._pair_coverage_max(results, 24, 32, t, None, dom32, metric) is None
    nan_metric = lambda res, T: float("nan")  # noqa: E731
    assert phy042._pair_coverage_max(results, 24, 32, t, dom24, dom32, nan_metric) == math.inf


def test_missing_common_domain_is_outside_domain_not_a_coverage_failure():
    """None-Abdeckung = keine gemeinsame Domaene -> OUTSIDE_DOMAIN, nicht COVERAGE_LEAK."""
    a = (0.56, 0.59)
    assert phy042._pair_verdict(0.575, None, a, None, None, bracket_valid=False) == "OUTSIDE_DOMAIN"
    assert phy042._pair_verdict(None, None, a, None, None, bracket_valid=False) == "NO_CROSSING"
    # Positiv-Kontrolle: echte Abdeckungsverletzung bleibt COVERAGE_LEAK
    assert phy042._pair_verdict(0.575, a, a, 0.5, 0.0, bracket_valid=True) == "COVERAGE_LEAK"
    assert phy042._pair_verdict(0.575, a, a, 0.0, 0.0, bracket_valid=True) == "QUOTABLE"
    # NaN/inf-Abdeckung bleibt fail-closed
    assert phy042._pair_verdict(0.575, a, a, float("inf"), 0.0, bracket_valid=True) == "COVERAGE_LEAK"


def test_global_gate_max_keeps_nan_and_inf():
    """builtin max(0.0, nan) == 0.0 liess eine NaN-Metrik verschwinden (Audit A-1)."""
    assert max(0.0, float("nan")) == 0.0  # dokumentierte Falle
    assert phy042._max_fail_closed(0.0, float("nan")) == math.inf
    assert phy042._max_fail_closed(0.0, float("inf")) == math.inf
    assert phy042._max_fail_closed(1e-4, 5e-5) == 1e-4          # Negativ-Kontrolle
    assert phy042._max_fail_closed(1e-4, 2e-4) == 2e-4


def test_assess_pair_wires_bracket_coverage_and_verdict():
    """Der Paar-Helfer _assess_pair Ende-zu-Ende (Audit A-2): Klammer,
    Paar-Coverage und Verdikt -- mit Stub-Metriken. Die Produktions-Metriken
    und die Restmassen-Haelfte deckt
    test_assess_pair_default_metrics_on_a_real_wl_result (Issue #63 F-5)."""
    t, ya, yb = _crossing_curves(0.647)
    curves = {24: {"y2": ya}, 32: {"y2": yb}}
    results = {(24, 0): "a", (32, 0): "b", (48, 0): "c"}
    clean = lambda res, T: 0.0  # noqa: E731
    wide = {24: (0.56, 0.65), 32: (0.52, 0.67)}
    tb, reason = phy042._assess_pair(t, curves, results, 24, 32, wide,
                                     leak_metric=clean, unc_metric=clean)
    assert abs(tb - 0.647) < 1e-9 and reason == "QUOTABLE"          # Positiv-Kontrolle
    narrow = {24: (0.56, 0.6475), 32: (0.52, 0.67)}                 # Klammer 0.650 ungedeckt
    assert phy042._assess_pair(t, curves, results, 24, 32, narrow,
                               leak_metric=clean, unc_metric=clean)[1] == "OUTSIDE_DOMAIN"
    leaky = lambda res, T: 0.5 if res == "b" and abs(T - 0.60) < 1e-9 else 0.0  # noqa: E731
    assert phy042._assess_pair(t, curves, results, 24, 32, wide,
                               leak_metric=leaky, unc_metric=clean)[1] == "COVERAGE_LEAK"
    foreign = lambda res, T: 0.5 if res == "c" else 0.0  # noqa: E731  (fremdes L=48)
    assert phy042._assess_pair(t, curves, results, 24, 32, wide,
                               leak_metric=foreign, unc_metric=foreign)[1] == "QUOTABLE"
    assert phy042._assess_pair(t, curves, results, 24, 32, {32: (0.52, 0.67)},
                               leak_metric=clean, unc_metric=clean)[1] == "OUTSIDE_DOMAIN"


# --- Issue #63: Regressionsschutz fuer Claim-Pfad und Domaenengrenzen --------

def test_common_domain_is_nan_safe_and_order_independent():
    """K-1: builtin max/min verlieren NaN je nach Argumentreihenfolge."""
    nan = float("nan")
    assert max(0.56, nan) == 0.56 and math.isnan(max(nan, 0.56))  # dokumentierte Falle
    a = (0.56, 0.60)
    for bad in ((nan, 0.67), (0.52, nan), (float("-inf"), 0.67), (0.52, float("inf"))):
        assert phy042._common_domain(a, bad) is None
        assert phy042._common_domain(bad, a) is None
        # die drei Paar-Gates, in beiden Reihenfolgen fail-closed
        assert not phy042._pair_inside_domains(0.58, a, bad)
        assert not phy042._pair_inside_domains(0.58, bad, a)
    assert phy042._common_domain(a, (0.61, 0.67)) is None        # leerer Schnitt
    assert phy042._common_domain(a, None) is None
    # Negativ-Kontrolle: endliche Grenzen liefern den Schnitt
    assert phy042._common_domain(a, (0.52, 0.59)) == (0.56, 0.59)
    assert phy042._pair_inside_domains(0.58, a, (0.52, 0.59))


def test_nan_domain_bound_fails_closed_through_the_pair_gates():
    """K-1 durch _assess_pair: eine NaN-Grenze ist keine gemeinsame Domaene,
    unabhaengig davon, an welchem L sie haengt."""
    nan = float("nan")
    t, ya, yb = _crossing_curves(0.647)
    curves = {24: {"y2": ya}, 32: {"y2": yb}}
    results = {(24, 0): "a", (32, 0): "b"}
    clean = lambda res, T: 0.0  # noqa: E731
    good = (0.52, 0.67)
    assert phy042._assess_pair(t, curves, results, 24, 32, {24: good, 32: good},
                               leak_metric=clean, unc_metric=clean)[1] == "QUOTABLE"
    for bad in ((nan, 0.67), (0.52, nan)):
        for bounds in ({24: bad, 32: good}, {24: good, 32: bad}):
            assert phy042._pair_coverage_max(results, 24, 32, t, bounds[24],
                                             bounds[32], clean) is None
            assert not phy042._bracket_in_common_domain(
                t, ya, 24, yb, 32, 0.647, bounds[24], bounds[32])
            assert phy042._assess_pair(t, curves, results, 24, 32, bounds,
                                       leak_metric=clean, unc_metric=clean
                                       )[1] == "OUTSIDE_DOMAIN"


def _wl_result(L, *, peak=0.0, hole=False):
    """Kleines echtes WLResult: lng so, dass das kanonische Gewicht bei allen
    T des Gitters eine schmale Glocke um ``peak`` ist (Breite 1, Bins 0.5)."""
    centers = np.linspace(-10.0, 10.0, 41)
    lng = centers / 0.6 - 0.5 * (centers - peak) ** 2
    mask = np.ones(41, dtype=bool)
    if hole:
        mask[20] = False                                  # Bin bei E = 0
    return phy042._p41.WLResult(L=L, n=2 * L * L, centers=centers, lng=lng,
                           mask=mask, micro_x={}, micro_y={}, wl_sweeps=1,
                           prod_sweeps=1)


def test_assess_pair_default_metrics_on_a_real_wl_result():
    """F-2/F-3: die Produktions-Metriken (canonical_edge_leak, _uncovered_mass)
    laufen durch _assess_pair -- ohne Stubs, auf echten WLResult-Objekten."""
    t, ya, yb = _crossing_curves(0.647)
    curves = {24: {"y2": ya}, 32: {"y2": yb}}
    wide = {24: (0.56, 0.65), 32: (0.52, 0.67)}

    def verdict(res24, res32):
        return phy042._assess_pair(t, curves, {(24, 0): res24, (32, 0): res32},
                                   24, 32, wide)[1]

    clean24, clean32 = _wl_result(24), _wl_result(32)
    assert phy042.canonical_edge_leak(clean32, 0.6) < 1e-12
    assert phy042._uncovered_mass(clean32, 0.6) == 0.0
    assert verdict(clean24, clean32) == "QUOTABLE"              # Positiv-Kontrolle
    # Restmasse: Loch im Glockenzentrum, Fensterraender bleiben leer
    holed = _wl_result(32, hole=True)
    assert phy042.canonical_edge_leak(holed, 0.6) < 1e-12
    assert phy042._uncovered_mass(holed, 0.6) > 0.1
    assert verdict(clean24, holed) == "COVERAGE_UNCOVERED"
    assert verdict(holed, clean32) == "COVERAGE_UNCOVERED"
    # Kantenleck: Glocke am oberen Fensterrand, alle Bins besetzt
    edgy = _wl_result(32, peak=9.5)
    assert phy042.canonical_edge_leak(edgy, 0.6) > 0.5
    assert phy042._uncovered_mass(edgy, 0.6) == 0.0
    assert verdict(clean24, edgy) == "COVERAGE_LEAK"


def test_global_gate_aggregation_is_nan_safe():
    """F-1: die globale Leck-/Restmassen-Aggregation von run_phy042 als reine
    Funktion. Mit builtin max ginge eine NaN-Metrik verloren (0.0 -> PASS)."""
    t = [0.55, 0.60, 0.65]
    results = {(24, 0): "a", (32, 0): "b", (32, 1): "c"}
    domains = {24: None, 32: np.array([True, True, False])}

    def metric(bad_res, bad_T, value):
        return lambda res, T: value if (res, T) == (bad_res, bad_T) else 1e-6

    nan = float("nan")
    for agg, args in ((phy042._global_leak_max, (results, t)),
                      (phy042._global_uncovered_max, (results, t, domains))):
        # NaN an erster, mittlerer und letzter Position der Iteration
        for where in (("a", 0.55), ("b", 0.60), ("c", 0.60)):
            got = agg(*args, metric(*where, nan))
            assert got == math.inf and not got < 1e-3
            assert agg(*args, metric(*where, float("inf"))) == math.inf
        # Negativ-Kontrollen: endliche Werte liefern das echte Maximum
        assert agg(*args, metric("b", 0.60, 0.5)) == 0.5
        assert agg(*args, metric("none", 0.0, 0.0)) == 1e-6
    # Restmasse bindet nur in der Domaene: L=32 bei T=0.65 liegt ausserhalb,
    assert phy042._global_uncovered_max(results, t, domains,
                                        metric("c", 0.65, 0.5)) == 1e-6
    assert phy042._global_uncovered_max(results, t, domains,
                                        metric("c", 0.65, nan)) == 1e-6
    # das Leck-Gate und ein L ohne gemessene Domaene zaehlen dagegen ueberall
    assert phy042._global_leak_max(results, t, metric("c", 0.65, nan)) == math.inf
    assert phy042._global_uncovered_max(results, t, domains,
                                        metric("a", 0.65, nan)) == math.inf


def test_global_gate_defaults_use_the_production_metrics():
    """F-1/F-3: ohne metric-Argument greifen die Produktions-Metriken."""
    t = phy042._T_GRID
    full = {32: np.ones(len(t), dtype=bool)}
    clean = {(32, 0): _wl_result(32)}
    assert phy042._global_leak_max(clean, t) < 1e-12
    assert phy042._global_uncovered_max(clean, t, full) == 0.0
    assert phy042._global_leak_max({(32, 0): _wl_result(32, peak=9.5)}, t) > 0.5
    assert phy042._global_uncovered_max({(32, 0): _wl_result(32, hole=True)},
                                        t, full) > 0.1


# --- Legacy-API fuer PHY044-WL (Supersession #50/#52/#54, 2026-11-19) -------

def _legacy_curves(spreads):
    """Zwei Walker-Kurven, deren Spread je T genau `spreads` ist."""
    base = np.linspace(1.0, 0.2, len(spreads))
    return [base, base + np.asarray(spreads, dtype=float)]


def test_legacy_walker_domain_finite_spread_is_measured():
    """Positivkontrolle: endlicher Spread < DOMAIN_THRESHOLD ueber das ganze
    Gitter -> gemessene Domaene bis zum letzten T (sonst waere der
    NaN-Test unten auch mit einer immer-leeren Domaene gruen)."""
    t = phy042._T_GRID
    d = phy042._walker_domain(_legacy_curves([0.01] * len(t)), t)
    assert d["measured"] is True and d["n_walkers"] == 2
    assert bool(d["mask"].all())
    assert d["tmax"] == float(t[-1]) and d["reason"] is None
    st = phy042._domain_status(d)
    assert st["max_spread"] == pytest.approx(0.01)
    # Schwellen-Rand: Spread == thr ist NICHT innerhalb (Annahme nur bei s < thr)
    edge = [0.01, 0.01, phy042.DOMAIN_THRESHOLD] + [0.0] * (len(t) - 3)
    d_edge = phy042._walker_domain(_legacy_curves(edge), t)
    assert d_edge["tmax"] == float(t[1])


@pytest.mark.parametrize("pos", [0, 5, -1])
def test_legacy_walker_domain_nan_spread_fails_closed(pos):
    """Negativtest (toleranz-und-nan-gate): ein NaN in EINER Walker-Kurve ist
    ein ungemessener Spread und darf die Domaene nicht verlaengern. NaN am
    unteren Rand -> leere Domaene (tmax None, Paar nicht quotierbar); NaN
    weiter oben -> Domaene endet davor. In #50 galt `s >= thr` (NaN = gueltig)."""
    t = phy042._T_GRID
    curves = _legacy_curves([0.01] * len(t))
    curves[1] = curves[1].copy()
    curves[1][pos] = math.nan
    d = phy042._walker_domain(curves, t)
    k = pos % len(t)
    assert d["measured"] is True
    assert int(d["mask"].sum()) == k
    assert not bool(d["mask"][k:].any())
    if k == 0:
        assert d["tmax"] is None and d["reason"] is not None
    else:
        assert d["tmax"] == float(t[k - 1])


def test_legacy_walker_domain_single_walker_is_not_a_measurement():
    """< 2 Walker: keine Messung (tmax None + Grund), Gate-Maske konservativ."""
    t = phy042._T_GRID
    d = phy042._walker_domain([np.ones(len(t))], t)
    assert d["measured"] is False and d["tmax"] is None
    assert bool(d["mask"].all()) and "n_walkers=1 < 2" in d["reason"]
    assert phy042._domain_status(d)["max_spread"] is None


def test_legacy_api_is_additive_and_drops_pair_domain_limit():
    """Supersession: nur die fuenf Symbole fuer PHY044-WL; der #50-Paar-Pfad
    (`_pair_domain_limit`, 'partial' quotierbar) bleibt draussen, der
    Claim-Pfad laeuft ueber _assess_pairs (#74)."""
    for name in ("PHY042_REPORT_V01", "DOMAIN_THRESHOLD", "_UNMEASURED_REASON",
                 "_walker_domain", "_domain_status"):
        assert hasattr(phy042, name), name
    assert not hasattr(phy042, "_pair_domain_limit")
    assert not hasattr(phy042, "_pair_quotable")
    assert (ROOT / phy042.PHY042_REPORT_V01).is_file()
