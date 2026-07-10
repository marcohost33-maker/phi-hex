"""Schnelle Korrektheits-Gates fuer PHY041 (WL-entropischer Helicity, honeycomb).

Die teure WL-g(E)-Konvergenz + Wolff-Quervalidierung laufen als slow/lokal
(Gate-Log in results/). Die CI-schnellen Gates pruefen MC-frei: das exakte
Aggregat-/Geometrie-Orakel (Upsilon_2(0) = 3/4 J), die lokale Energie-Delta
des skalar-optimierten Kernels gegen die volle Energie, die Auto-Fenster-
Logik, das Leak-Mass, den deterministischen PHY032-Drift-Guard-Parser und
die Paar-Verdrahtung im honeycomb-Band.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

import conftest

phy041 = conftest._load(
    "phy041_honeycomb_wl",
    "260702 PHY041 honeycomb wang-landau entropic helicity v01.py")


def test_aligned_upsilon_exact_three_quarters():
    """Mikrokanonische Aggregat-Maschinerie auf der perfekt ausgerichteten
    Konfiguration -> Upsilon_2 = 3/4 J EXAKT (per Site, beide Richtungen),
    groessen-unabhaengig. Faengt Vorzeichen-/Projektions-/Normierungsfehler."""
    for L in (4, 6, 8):
        y0 = phy041.aligned_upsilon_zero(L)
        assert abs(y0 - 0.75) < 1e-12, (L, y0)


def test_local_delta_matches_full_energy():
    """Skalar-optimierte lokale Energie-Delta == Differenz der vollen Energie
    (der Kern-Invariant des Block-RNG-Kernels; rtol ~1e-9)."""
    _, nbr_list, ei, ej, _, _, n = phy041.honeycomb_arrays(5)
    rng = np.random.default_rng(7)
    thl = list(rng.uniform(0, 2 * math.pi, n))
    E = phy041._full_energy(np.asarray(thl), ei, ej)
    for _ in range(60):
        i = int(rng.integers(n))
        new = float(rng.uniform(0, 2 * math.pi))
        d = phy041._local_delta(thl, nbr_list[i], thl[i], new)
        thl[i] = new
        E_new = phy041._full_energy(np.asarray(thl), ei, ej)
        assert abs((E + d) - E_new) < 1e-9, (i, d, E_new - E)
        E = E_new


def test_window_from_anchors_margins_and_floor():
    """Fenster-Logik: Anker-Mittel +- k*std (mit Minimal-Margen), unterer
    Rand auf den Grundzustands-Floor geclampt (unreachable-Bin-Guard)."""
    lo, hi = phy041.window_from_anchors(-1.22, 0.02, -1.04, 0.03)
    assert lo == pytest.approx(-1.22 - 0.12)          # 6*std > min_lo
    assert hi == pytest.approx(-1.04 + 0.24)          # 8*std > min_hi
    lo2, hi2 = phy041.window_from_anchors(-1.22, 0.001, -1.04, 0.001)
    assert lo2 == pytest.approx(-1.22 - 0.05)         # min_lo greift
    assert hi2 == pytest.approx(-1.04 + 0.10)         # min_hi greift
    lo3, _ = phy041.window_from_anchors(-1.46, 0.02, -1.04, 0.03)
    assert lo3 == pytest.approx(phy041.E0_PER_SITE + 0.03)  # Floor-Clamp
    assert lo < hi and lo2 < hi2


def test_canonical_edge_leak_metric():
    """Leak-Mass: kanonisches Gewicht in den Rand-Bins. Interior-konzentrierte
    Verteilung -> Leak ~ 0; Rand-konzentrierte -> Leak ~ 1."""
    nbins = 40
    centers = np.linspace(-120.0, -80.0, nbins)
    lng = np.zeros(nbins)
    micro = {k: np.zeros(nbins) for k in phy041.OBS}
    res = phy041.WLResult(L=8, n=100, centers=centers, lng=lng,
                          mask=np.ones(nbins, bool), micro_x=micro,
                          micro_y=micro, wl_sweeps=0, prod_sweeps=0)
    # T so, dass das Boltzmann-Gewicht im Innern sitzt: dE=40 ueber 40 Bins;
    # bei T=1 dominiert der unterste Bin -> Leak am unteren Rand.
    assert phy041.canonical_edge_leak(res, 1.0) > 0.9
    # Flaches g + hohes T -> Gewicht gleichverteilt: Leak ~ 6/40.
    assert phy041.canonical_edge_leak(res, 1e6) == pytest.approx(6 / 40, rel=0.01)
    # Interior-Peak in g -> Leak winzig.
    lng2 = -0.5 * ((centers - centers[nbins // 2]) / 1.5) ** 2 * 10
    res2 = phy041.WLResult(L=8, n=100, centers=centers, lng=lng2,
                           mask=np.ones(nbins, bool), micro_x=micro,
                           micro_y=micro, wl_sweeps=0, prod_sweeps=0)
    assert phy041.canonical_edge_leak(res2, 1e6) < 1e-6


def test_uncovered_canonical_mass_two_level_oracle():
    """Coverage-Massen-Gate (Backport PHY042, Code-Audit H2): 2-Bin-
    Boltzmann-Orakel - bei gleicher Entropie ist die Masse des unbesetzten
    Bins exakt exp(-dE/T)/(1+exp(-dE/T)); volle Abdeckung -> 0."""
    nbins = 2
    micro = {k: np.zeros(nbins) for k in phy041.OBS}
    res = phy041.WLResult(L=8, n=100, centers=np.array([-10.0, -9.0]),
                          lng=np.zeros(nbins),
                          mask=np.array([True, False]), micro_x=micro,
                          micro_y=micro, wl_sweeps=0, prod_sweeps=0)
    T = 0.5
    expected = math.exp(-1.0 / T) / (1.0 + math.exp(-1.0 / T))
    got = phy041.uncovered_canonical_mass(res, T)
    assert got == pytest.approx(expected, rel=1e-12)
    res_full = phy041.WLResult(L=8, n=100, centers=res.centers, lng=res.lng,
                               mask=np.array([True, True]), micro_x=micro,
                               micro_y=micro, wl_sweeps=0, prod_sweeps=0)
    assert phy041.uncovered_canonical_mass(res_full, T) == 0.0


def test_parse_phy032_grid_drift_guard():
    """Drift-Guard: der Parser liest EXAKT die committed PHY032-Evidenz
    (Stichproben zeichengenau; 8 T-Punkte je L; alle drei L vorhanden)."""
    grid = phy041.parse_phy032_grid()
    assert sorted(grid) == [12, 24, 48]
    for L in (12, 24, 48):
        assert len(grid[L]) == 8, (L, len(grid[L]))
    assert grid[12][0] == (0.56, 0.5368, 0.0035)
    assert grid[24][3] == (0.5975, 0.4726, 0.0046)
    assert grid[48][-1] == (0.6475, 0.2326, 0.0131)


def test_pair_tbkt_recovers_synthetic_wm_in_honeycomb_band():
    """Die (aus PHY040 wiederverwendete) C-eliminierte Paar-Bedingung
    rekonstruiert ein bekanntes T_true im honeycomb-Band auf dem
    PHY041-Analyse-Gitter (fails-before-faehig)."""
    T_true = 0.575
    C = -1.0
    grid = np.round(np.arange(0.52, 0.6701, 0.005), 4)

    def y2(T, L):
        q = 2 * math.log(L) + C + math.log(L) * (T - T_true)
        return (2 * T / math.pi) * (1.0 + 1.0 / q)
    for La, Lb in [(12, 16), (12, 24), (16, 24)]:
        ya = [y2(T, La) for T in grid]
        yb = [y2(T, Lb) for T in grid]
        tb = phy041.tbkt_pair_from_curves(grid, ya, La, yb, Lb)
        assert tb is not None and abs(tb - T_true) < 1e-3, (La, Lb, tb)


def test_honeycomb_arrays_geometry():
    """Bond-Projektionen: je Zelle sum(ax^2) = sum(ay^2) = 1.5 (Isotropie
    des T=0-Grenzwerts); 3 Kanten je Zelle, Einheitslaenge."""
    lat, nbr_list, ei, ej, ax, ay, n = phy041.honeycomb_arrays(6)
    assert n == 2 * 36 and len(ax) == 3 * 36
    assert all(len(t) == 3 for t in nbr_list)
    assert float(np.sum(ax ** 2)) == pytest.approx(1.5 * 36)
    assert float(np.sum(ay ** 2)) == pytest.approx(1.5 * 36)
    assert np.allclose(ax ** 2 + ay ** 2, 1.0)


@pytest.mark.slow
def test_wl_honeycomb_smoke_energy_matches_wolff():
    """Langsam: kurzer WL-Lauf (L=8) reproduziert die Wolff-Energie grob
    und liefert ein leak-freies Analyse-Fenster."""
    res, meta = phy041.wang_landau_honeycomb(
        8, lnf_final=2e-4, prod_sweeps=4000, verbose=False)
    c = phy041.upsilon_curves(res, np.array([0.55, 0.62]))
    ref = phy041.wolff_reference(8, 0.62, n_seeds=3, n_measure=150,
                                 n_burn=100)
    assert abs(c["E"][1] / res.n - ref["E_ps"]) < 0.05
    assert phy041.canonical_edge_leak(res, 0.62) < 1e-3
    assert meta["window_ps"][0] >= phy041.E0_PER_SITE
