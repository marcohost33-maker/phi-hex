"""PHY031 Honeycomb-T_BKT: Korrektheits-Gates (per-Site, z=3).

Drei Ebenen:
(a) GEOMETRIE-ORAKEL (instant, fails-before-faehig): der per-Site-T=0-
    Spinwellen-Grenzwert des Honeycomb-Gitters ist analytisch EXAKT 3/4 J
    (N=2L^2 Knoten, 3L^2 Bonds, sum (e.delta)^2/Zelle = 3/2). Jeder Knoten
    hat exakt Koordinationszahl 3. Beides groessen-unabhaengig.
(b) SANDVIK-PAAR-ORAKEL (instant): kuenstliche Weber-Minnhagen-Daten mit
    bekanntem T_true -> die C-freie Paar-Schaetzung MUSS T_true rekonstruieren
    (unabhaengig vom Monte-Carlo-Sampling).
(c) SMOKE (slow): Mini-Wolff-Lauf -> NK-Crossing im physikalischen Band.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from phy031_honeycomb import (
    REF_HONEYCOMB,
    SW_LIMIT_HONEYCOMB,
    XYConfig,
    build_honeycomb_lattice,
    helicity_from_ensemble,
    helicity_terms,
    measure_honeycomb,
    nk_crossing,
    sandvik_pair_tbkt,
)

TOL = 1e-6


def _upsilon_T0(L: int) -> float:
    lat = build_honeycomb_lattice(L)
    t1, sa = helicity_terms(np.zeros(lat.n_nodes), lat, lat.pos)
    return helicity_from_ensemble([t1], [sa], XYConfig(J=1.0, T=1e-12),
                                  lat.n_nodes)


def test_honeycomb_T0_is_three_quarters():
    """Υ(T→0) MUSS exakt 3/4 J sein (per-Site-Spinwellen-Grenzwert, z=3)."""
    ups = _upsilon_T0(L=8)
    assert abs(ups - SW_LIMIT_HONEYCOMB) < TOL, (
        f"Upsilon(0)={ups:.8f}, erwartet {SW_LIMIT_HONEYCOMB} (per-Site)")


def test_honeycomb_T0_size_independent():
    """Per-Site-Υ(0) ist groessen-unabhaengig (= 3/4) ueber mehrere L."""
    vals = [_upsilon_T0(L) for L in (4, 6, 8, 12)]
    for v in vals:
        assert abs(v - SW_LIMIT_HONEYCOMB) < TOL, vals
    assert max(vals) - min(vals) < TOL, vals


def test_honeycomb_coordination_is_three():
    """Jeder Knoten hat exakt 3 Nachbarn (Honeycomb z=3); N=2L^2, 3L^2 Bonds."""
    L = 6
    lat = build_honeycomb_lattice(L)
    assert lat.n_nodes == 2 * L * L
    assert len(lat.edges) == 3 * L * L
    deg = [0] * lat.n_nodes
    for i, j in lat.edges:
        deg[i] += 1
        deg[j] += 1
    assert all(d == 3 for d in deg), set(deg)


def test_honeycomb_bonds_are_unit_length():
    """Alle Bond-Verschiebungen haben Einheitslaenge (NN-Abstand a=1)."""
    lat = build_honeycomb_lattice(5)
    for dx, dy in lat.edge_disp:
        assert abs(math.hypot(dx, dy) - 1.0) < 1e-12, (dx, dy)


def _wm(T_bkt: float, L: float, C: float) -> float:
    return (2.0 * T_bkt / math.pi) * (1.0 + 1.0 / (2.0 * math.log(L) + C))


class _Fake:
    def __init__(self, upsilon_mean: float):
        self.upsilon_mean = upsilon_mean


@pytest.mark.parametrize("T_true,C_true", [(0.573, 2.0), (0.55, 0.8),
                                           (0.60, 4.0)])
def test_sandvik_pair_recovers_synthetic_T_true(T_true, C_true):
    """C-freie (L,2L)-Paar-Schaetzung rekonstruiert T_true aus WM-Daten.

    Reine WM-Form erfuellt die Paar-Bedingung exakt nur bei T_true; eine
    glatte T-Abhaengigkeit (slope*(T-T_true)*lnL) ausserhalb erzeugt den
    Nulldurchgang genau bei T_true. Mutations-/Vorzeichenfehler verfehlen ihn.
    """
    # Enges Fenster + moderate slope: haelt R(T,L)=pi*Upsilon/(2T)-1 strikt
    # positiv (keine 1/R-Pole), so wie reale MC-Daten im Scan-Fenster.
    L = 8
    T_grid = [round(T_true - 0.02 + 0.005 * k, 4) for k in range(9)]
    slope = 0.3
    data_L = {T: _Fake(_wm(T_true, L, C_true) + slope * (T - T_true)
                       * math.log(L)) for T in T_grid}
    data_2L = {T: _Fake(_wm(T_true, 2 * L, C_true) + slope * (T - T_true)
                        * math.log(2 * L)) for T in T_grid}
    tc = sandvik_pair_tbkt(data_L, data_2L)
    assert tc is not None, "kein Nulldurchgang der Paar-Bedingung"
    assert abs(tc - T_true) < 0.006, f"Paar T_BKT={tc:.4f} verfehlt {T_true}"


@pytest.mark.slow
def test_honeycomb_crossing_in_physical_band():
    """Mini-Wolff-Lauf: per-L NK-Crossing im Band um die Referenz 0.573.

    Kleine L, wenige Sweeps - kein Praezisions-Anspruch, nur dass die
    Pipeline einen physikalisch sinnvollen Crossing liefert (finite-size
    von oben; klar unter dem Triangular-Wert ~1.4).
    """
    data = {}
    for T in (0.59, 0.65, 0.71):
        data[T] = measure_honeycomb(8, T, n_seeds=2, n_measure=120,
                                    n_burn=100, master_seed=42)
    tc = nk_crossing(data)
    assert tc is not None, "kein NK-Crossing im Scan-Fenster"
    assert 0.50 < tc < 0.80, f"Crossing T_x={tc:.4f} ausserhalb Band"
    assert abs(REF_HONEYCOMB - 0.573) < 1e-9
