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
    spec = (ROOT / "spec" /
            "260703 PHI HEX honeycomb reference conventions audit v01.md"
            ).read_text(encoding="utf-8")
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
