"""Regressions-Gates fuer den PHY028-Geometrie-Fix (Code-Audit H1) und die
Pol-Guards der C-freien Paar-Root-Finder (Code-Audit M1), 2026-07-10.

H1: Die (min,max)-Kantensortierung invertierte fuer die 2L Wrap-Bonds des
    LxL-Torus das Vorzeichen von edge_disp relativ zum gespeicherten (i,j).
    T1 (cos, gerade) blieb korrekt - das T=0-Gate KONNTE den Fehler nicht
    sehen; sin_accum war systematisch falsch (verifiziert: 25% Fehler bei
    L=8 unter uniformem Twist). Orakel hier: exakter uniformer Twist.

M1: Auf dem physikalischen Weber-Minnhagen-Ast ist 1/R = 2 ln L + C > 0.
    Ohne R>0-Guard interpretierten die Root-Finder den Pol von 1/R (am
    finite-L-NK-Crossing) als Nulldurchgang und lieferten Pol-Artefakte
    als T_BKT (verifiziert: synthetisches R2-Nulldurchgangs-Szenario ergab
    faelschlich T_BKT=1.449).
"""
from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np

import conftest

phy028 = conftest._load(
    "phy028_square_validation", "260602 PHY028 square validation v01.py")
phy029 = conftest._load(
    "phy029_triangular_sandvik", "260602 PHY029 triangular sandvik v01.py")
phy030v02 = conftest._load(
    "phy030_tbkt_wm_logfit",
    "260604 PHY030 triangular tbkt per-site v02 wm-logfit.py")
phy040 = conftest._load(
    "phy040_wang_landau", "260616 PHY040 wang-landau entropic helicity v01.py")


# ---------------------------------------------------------------------------
# H1: Quadratgitter-Geometrie (Wrap-Bond-Orientierung)
# ---------------------------------------------------------------------------

def test_square_edge_disp_is_min_image_for_all_bonds():
    """edge_disp muss fuer JEDEN Bond (auch Wrap) der min-image-Vektor
    r_j - r_i des gespeicherten Paars (i, j) sein."""
    for L in (4, 8):
        lat = phy028.build_square_lattice(L)
        assert len(lat.edges) == 2 * L * L
        for k, (i, j) in enumerate(lat.edges):
            xi, yi = divmod(i, L)
            xj, yj = divmod(j, L)
            dx = (xj - xi + L // 2) % L - L // 2
            dy = (yj - yi + L // 2) % L - L // 2
            assert (float(dx), float(dy)) == lat.edge_disp[k], (L, k, i, j)


def test_square_uniform_twist_exact_helicity_terms():
    """Uniformer Twist theta = (2 pi / L) * x: exakt loesbar.

    T1 = L^2 cos(k) (x-Bonds; y-Bonds projizieren zu 0),
    sin_accum = -L^2 sin(k) mit k = 2 pi / L. Vor dem H1-Fix fehlten
    2L*sin(k) im Betrag von sin_accum (Wrap-Bonds mit falschem Vorzeichen).
    """
    for L in (4, 8, 10):
        lat = phy028.build_square_lattice(L)
        k = 2.0 * math.pi / L
        theta = np.array([k * (i // L) for i in range(L * L)])
        t1, sa = phy028.square_helicity_terms(theta, lat)
        assert abs(t1 - L * L * math.cos(k)) < 1e-9, L
        assert abs(sa - (-L * L * math.sin(k))) < 1e-9, L


def test_square_aligned_upsilon_exact_one():
    """T=0-Orakel bleibt exakt: Upsilon(aligned) = 1.0 J per Site."""
    lat = phy028.build_square_lattice(12)
    t1, sa = phy028.square_helicity_terms(np.zeros(lat.n_nodes), lat)
    ups = phy028.square_helicity_from_ensemble(
        [t1], [sa], phy028.XYConfig(J=1.0, T=0.01), lat.n_nodes)
    assert abs(ups - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# M1: Pol-Guards der Paar-Root-Finder
# ---------------------------------------------------------------------------

def _res(u):
    return SimpleNamespace(upsilon_mean=u)


def _ups_from_R(R, T):
    return 2.0 * T / math.pi * (R + 1.0)


def test_phy028_pair_rejects_pole_artifact():
    """R2-Nulldurchgang im Scan-Fenster (kein WM-Nulldurchgang): der Guard
    muss None liefern statt des Pol-Artefakts (frueher: T_BKT=1.449)."""
    Ts = [1.40, 1.42, 1.44, 1.46, 1.48]
    data_L = {T: _res(_ups_from_R(0.30, T)) for T in Ts}
    R2s = {1.40: 0.10, 1.42: 0.05, 1.44: 0.02, 1.46: -0.02, 1.48: -0.05}
    data_2L = {T: _res(_ups_from_R(R2s[T], T)) for T in Ts}
    assert phy028.sandvik_pair_tbkt(8, data_L, data_2L) is None


def test_phy029_pair_rejects_pole_artifact():
    """Gleicher Pol-Guard in PHY029 (Cache-injiziert, kein MC)."""
    r1, r2 = 4, 8
    Ts = [1.40, 1.42, 1.44, 1.46, 1.48]
    R2s = {1.40: 0.10, 1.42: 0.05, 1.44: 0.02, 1.46: -0.02, 1.48: -0.05}
    cache = {}
    for T in Ts:
        cache[(r1, T)] = _res(_ups_from_R(0.30, T))
        cache[(r2, T)] = _res(_ups_from_R(R2s[T], T))
    assert phy029.sandvik_pair(r1, r2, Ts, cache, {}) is None


def test_phy029_pair_still_finds_genuine_wm_root():
    """Positiv-Kontrolle: echte WM-Wurzel wird weiterhin gefunden.

    Konstruktion: 1/R(T,L) = 2 ln L + C + k*(T-T_true)*ln L -> die
    Paar-Bedingung ist exakt bei T_true erfuellt und crosst linear.
    """
    r1, r2 = 4, 8
    L1, L2 = 2 * r1 + 1, 2 * r2 + 1
    T_true, C, slope = 1.44, 2.0, 3.0
    Ts = [1.40, 1.42, 1.44, 1.46, 1.48]
    cache = {}
    for T in Ts:
        for r, L in ((r1, L1), (r2, L2)):
            inv_R = 2.0 * math.log(L) + C + slope * (T - T_true) * math.log(L)
            cache[(r, T)] = _res(_ups_from_R(1.0 / inv_R, T))
    tb = phy029.sandvik_pair(r1, r2, Ts, cache, {})
    assert tb is not None and abs(tb - T_true) < 1e-9, tb


def test_phy030v02_pair_rejects_pole_artifact():
    """Pol-Guard im zentralen tbkt_pair_C_eliminated (via Import auch von
    PHY032/033/034/035 genutzt)."""
    Ts = [1.40, 1.42, 1.44, 1.46, 1.48]
    R2s = [0.10, 0.05, 0.02, -0.02, -0.05]
    ups1 = [_ups_from_R(0.30, T) for T in Ts]
    ups2 = [_ups_from_R(r, T) for r, T in zip(R2s, Ts)]
    sem = [1e-4] * len(Ts)
    assert phy030v02.tbkt_pair_C_eliminated(
        Ts, ups1, sem, 9, ups2, sem, 19) is None


def test_phy040_pair_from_curves_rejects_pole_artifact():
    """Pol-Guard in tbkt_pair_from_curves (via Import auch PHY041/042)."""
    Ts = np.array([0.56, 0.58, 0.60, 0.62, 0.64])
    R2s = [0.10, 0.05, 0.02, -0.02, -0.05]
    y2a = [_ups_from_R(0.30, T) for T in Ts]
    y2b = [_ups_from_R(r, T) for r, T in zip(R2s, Ts)]
    assert phy040.tbkt_pair_from_curves(Ts, y2a, 24, y2b, 48) is None
