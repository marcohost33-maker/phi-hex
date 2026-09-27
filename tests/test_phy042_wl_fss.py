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
    """REF_BAND (T-Form) muss mit der Vertragsquelle
    spec/260703 ... reference conventions audit uebereinstimmen."""
    # Vertrag = 260703 + Provenienz-Nachtrag 260926 (Issue #45 §3)
    spec = "\n".join(
        (ROOT / "spec" / name).read_text(encoding="utf-8") for name in (
            "260703 PHI HEX honeycomb reference conventions audit v01.md",
            "260926 PHI HEX honeycomb reference provenance addendum v01.md"))
    band = phy042.REF_BAND
    assert band["multi_lattice"] == 0.573 and "0.573" in spec
    # Haertung 2026-07-10 (Code-Audit): exakte Token statt rstrip("0")-
    # Degradierung ("0.5800" -> "0.58" haette fast alles gematcht). Wert UND
    # Unsicherheit muessen zeichengenau in der Vertragsquelle stehen.
    for key, (val, sig) in ((k, v) for k, v in band.items()
                            if isinstance(v, tuple)):
        assert f"{val:.4f}" in spec or f"{val:.3f}" in spec, key
        assert f"{sig:g}" in spec, (key, sig)
    # beta-Konversions-Kanaele explizit (exakte 4-Dezimal-Token)
    for tok in ("0.5928", "0.6116", "0.5800"):
        assert tok in spec


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


# ---------------------------------------------------------------------------
# Issue #45 §2 (2026-09-26): Einzel-Walker-L ist KEINE Messung.
# ---------------------------------------------------------------------------

ERRATUM = ROOT / "results" / "260926 PHY042 domain semantics erratum.json"


def _dom(tmax, measured=True):
    """Minimal-Domaene fuer die Paar-Grenzen-Orakel."""
    return {"measured": measured, "tmax": tmax}


def test_single_walker_domain_is_null_with_reason():
    """< 2 Walker: spread/tmax None + Grund; die Gate-Maske bleibt voll
    (konservativ: Coverage-Gate bindet dann ueber das ganze Gitter)."""
    t = phy042._T_GRID
    d = phy042._walker_domain([np.linspace(1.0, 0.2, len(t))], t)
    assert d["measured"] is False and d["n_walkers"] == 1
    assert d["spread"] is None and d["tmax"] is None
    assert "Not a measurement" in d["reason"]
    assert d["mask"].all()
    st = phy042._domain_status(d)
    assert st["tmax"] is None and st["max_spread"] is None


def test_multi_walker_domain_measured_and_empty_domain_is_explicit():
    t = phy042._T_GRID
    base = np.linspace(1.0, 0.2, len(t))
    bump = np.where(t > 0.60, 0.1, 0.0)          # Spread 0.1 ab T>0.60
    d = phy042._walker_domain([base, base + bump, base], t)
    assert d["measured"] is True and d["reason"] is None
    assert d["tmax"] == pytest.approx(0.60)
    # gemessen, aber schon am unteren Rand ueber der Schwelle -> LEER
    d0 = phy042._walker_domain([base, base + 0.05], t)
    assert d0["measured"] is True and d0["tmax"] is None
    assert "empty validity domain" in d0["reason"]


def test_pair_domain_limit_fail_closed_and_order_independent():
    """Regression gegen den latenten Fail-open: frueher min(tmax_a, tmax_b)
    mit NaN fuer leere Domaenen -> min(0.6, nan) = 0.6 (quotierbar!),
    min(nan, 0.6) = nan. Jetzt symmetrisch fail-closed."""
    full, empty = _dom(0.60), _dom(None)
    unmeasured = _dom(None, measured=False)
    for a, b in ((full, empty), (empty, full)):
        assert phy042._pair_domain_limit(a, b) == (None, "empty")
        assert phy042._pair_quotable(0.55, a, b)[0] is False
    for a, b in ((full, unmeasured), (unmeasured, full)):
        assert phy042._pair_domain_limit(a, b) == (0.60, "partial")
    assert phy042._pair_domain_limit(unmeasured, unmeasured) == (
        None, "unmeasured")
    assert phy042._pair_quotable(0.55, unmeasured, unmeasured)[0] is False
    assert phy042._pair_domain_limit(_dom(0.585), full) == (
        0.585, "both_measured")
    assert phy042._pair_quotable(None, full, full)[0] is False
    assert phy042._pair_quotable(0.60, full, full)[0] is True     # Rand ok
    assert phy042._pair_quotable(0.6001, full, full)[0] is False


def test_walker_plan_min_walkers_is_fail_closed():
    """W4-Vertrag: >= 3 Walker an JEDEM L - Einzel-Walker-Bruecke L=24
    muss laut scheitern, nicht still mitlaufen."""
    assert phy042._walker_plan((24, 32, 48), 3) == {24: 1, 32: 3, 48: 3}
    with pytest.raises(ValueError, match="24"):
        phy042._walker_plan((24, 32, 48), 3, min_walkers=3)
    assert phy042._walker_plan((32, 48), 3, min_walkers=3) == {32: 3, 48: 3}
    with pytest.raises(ValueError):
        phy042._walker_plan((32, 48), 2, min_walkers=3)


def test_reanalysis_uses_two_sided_bounds_and_green_l24_fallback():
    """Reanalysis uses the same strict domain semantics as fresh runs."""
    import json
    rep = json.loads(
        (ROOT / phy042.PHY042_REPORT_V01).read_text(encoding="utf-8"))
    out = phy042.reanalyse_domains(rep)
    st = out["domain_status"]
    assert st["24"]["measured"] is False and st["24"]["tmax"] is None
    assert st["32"]["tmax"] == rep["domain_tmax_spread004"]["32"] == 0.6
    assert st["48"]["tmax"] == rep["domain_tmax_spread004"]["48"] == 0.585

    p = out["pairs"]["24_32"]
    assert p["basis_a"] == "PHY032_drift_guard"
    assert p["basis_b"] == "walker_spread"
    assert p["bounds_a"] is not None and p["bounds_b"] is not None
    assert phy042._pair_inside_bounds(
        p["tbkt"], tuple(p["bounds_a"]), tuple(p["bounds_b"])
    ) is p["quotable"]

    broken = json.loads(json.dumps(rep))
    broken["pass_gates"]["PASS_WL_Y2_MATCHES_PHY032_GRID_L24"] = False
    out_broken = phy042.reanalyse_domains(broken)
    p_broken = out_broken["pairs"]["24_32"]
    assert out_broken["effective_domain_bounds"]["24"] is None
    assert p_broken["basis_a"] == "unmeasured"
    assert p_broken["quotable"] is False


def test_committed_erratum_remains_lineage_but_runtime_reanalysis_is_current():
    """The v01 artifact stays immutable while current code is stricter."""
    import json
    committed = json.loads(ERRATUM.read_text(encoding="utf-8"))
    regen = json.loads(json.dumps(phy042._clean(phy042.domain_erratum())))
    expected_sha = (
        "19a9ce3c799401dbb55e519c09b390bee9dfbe23792f99b2e7c650c3eeefa3cf")
    assert committed["source_sha256"] == expected_sha
    assert regen["source_sha256"] == expected_sha
    assert regen["date"] == "2026-09-27"
    assert regen["pairs"]["24_32"]["basis_a"] == "PHY032_drift_guard"
    assert regen["pairs"]["24_32"]["basis_b"] == "walker_spread"


def test_run_phy042_report_path_emits_null_for_single_walker(monkeypatch):
    """Integrations-Gate durch den ECHTEN run_phy042-Pfad (MC ersetzt durch
    synthetische, deterministische Kurven): der JSON-Report fuehrt L=24
    als ungemessen (null + Grund) statt 0.67/0.0, bleibt allow_nan-frei
    serialisierbar, und die gemessenen L tragen ihre Domaene."""
    import json
    from types import SimpleNamespace

    t = phy042._T_GRID

    def fake_job(args):
        L, w, e_lo, e_hi, seed = args
        nb = 8
        return (L, w, SimpleNamespace(
            L=L, n=2 * L * L, wl_sweeps=100 + w, lng=np.zeros(nb),
            centers=np.arange(nb, dtype=float), mask=np.ones(nb, bool)))

    def fake_curves(res, t_grid):
        L = res.L
        # WM-artige Kurven: fallen, ab T>0.60 walker-abhaengig (Spread 0.1)
        y2 = 1.2 - 1.1 * (t_grid - 0.52) - 0.02 * math.log(L)
        y2 = y2 + np.where(t_grid > 0.60, 0.05 * (res.wl_sweeps - 100), 0.0)
        y4 = np.abs(t_grid - 0.62)                   # Dip bei T=0.62
        return {"T": t_grid, "y2": y2, "y4_scaled": y4,
                "E": -1.2 * res.n * np.ones(len(t_grid))}

    monkeypatch.setattr(phy042, "wolff_anchor",
                        lambda L, T, master_seed=42: (-1.3 + T, 0.01))
    monkeypatch.setattr(phy042, "_wl_job", fake_job)
    monkeypatch.setattr(phy042, "upsilon_curves", fake_curves)
    monkeypatch.setattr(phy042, "canonical_edge_leak", lambda res, T: 0.0)
    monkeypatch.setattr(phy042, "_uncovered_mass", lambda res, T: 0.0)
    monkeypatch.setattr(phy042, "wolff_reference", lambda L, T, master_seed=42:
                        {"E_ps": -1.2, "y2": 0.9, "y2_sem": 0.01})
    rep = phy042.run_phy042(max_workers=1)
    out = json.loads(json.dumps(phy042._clean(rep), allow_nan=False))
    assert out["walker_spread"]["24"] is None
    assert out["domain_tmax_spread004"]["24"] is None
    assert out["domain_status"]["24"]["measured"] is False
    assert "Not a measurement" in out["domain_status"]["24"]["reason"]
    for L in ("32", "48"):
        assert out["domain_status"][L]["measured"] is True
        assert out["domain_tmax_spread004"][L] == pytest.approx(0.60)
        assert len(out["walker_spread"][L]) == len(t)
    # Der synthetische VAL-B-Vergleich oben scheitert absichtlich deutlich:
    # ohne gruene Fallback-Evidenz bleibt L24 fail-closed ungemessen.
    assert out["effective_domain_basis"]["24"] == "unmeasured"
    assert out["effective_domain_bounds"]["24"] is None
    assert out["effective_domain_basis"]["32"] == "walker_spread"
    assert out["pair_domain_basis"]["24_32"]["basis_a"] == "unmeasured"
    assert out["pair_domain_basis"]["24_32"]["basis_b"] == "walker_spread"
    assert out["pair_domain_basis"]["32_48"]["basis_a"] == "walker_spread"
    assert out["pair_domain_basis"]["32_48"]["basis_b"] == "walker_spread"
    with pytest.raises(ValueError):
        phy042.run_phy042(max_workers=1, min_walkers=3)


def test_every_band_channel_has_provenance_row():
    """Issue #45 §3: jeder Band-Kanal fuehrt Quelle, berichtete Groesse und
    Beleg-Status; der Nachtrag-Spec listet jeden Kanal-Namen."""
    assert set(phy042.REF_PROVENANCE) == set(phy042.REF_BAND)
    allowed = {"current_method_anchor", "primary_text_verified",
               "superseded_historical"}
    addendum = (ROOT / "spec" /
                "260926 PHI HEX honeycomb reference provenance addendum v01.md"
                ).read_text(encoding="utf-8")
    for key, (src, reported, status) in phy042.REF_PROVENANCE.items():
        assert status in allowed, key
        assert src and reported, key
        assert f"| {key} |" in addendum, key
    # Historische Werte bleiben als Lineage, aktuelle Auswertung benutzt nur
    # den expliziten aktuellen Key-Satz.
    assert phy042.REF_BAND["upsilon_wl_T"] == (0.576, 0.003)
    assert "binder_beta" in phy042.REF_BAND
    assert phy042.REF_PROVENANCE["binder_beta"][2] == "superseded_historical"
    assert phy042.REF_PROVENANCE["helicity_ptep"][2] == "primary_text_verified"
    assert phy042.REF_PROVENANCE["upsilon_wl_T"][2] == "primary_text_verified"
    assert "binder_beta" not in phy042.REF_CURRENT_KEYS
    refs = phy042._reference_payload()
    assert set(refs["current"]) == set(phy042.REF_CURRENT_KEYS)
    assert "upsilon_beta" not in refs["current"]
    assert "upsilon_beta" in refs["lineage"]
    assert refs["provenance"]["upsilon_beta"]["status"] == "superseded_historical"


def test_pair_inside_effective_bounds_checks_lower_and_upper_edges():
    a = (0.56, 0.6475)
    b = (0.52, 0.60)
    assert phy042._pair_inside_bounds(0.575, a, b)
    assert not phy042._pair_inside_bounds(0.555, a, b)
    assert not phy042._pair_inside_bounds(0.61, a, b)
    assert not phy042._pair_inside_bounds(0.575, None, b)
    assert not phy042._pair_inside_bounds(None, a, b)
