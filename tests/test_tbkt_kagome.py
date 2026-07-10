"""PHY033 Kagome-T_BKT: Korrektheits-Gates (per-Site, z=4).

Vier Ebenen (alle fails-before-faehig, schnelle ohne teures Monte-Carlo):

(a) GEOMETRIE-ORAKEL (instant): der per-Site-T=0-Spinwellen-Grenzwert des
    Kagome-Gitters ist analytisch EXAKT 1 J (N=3L^2 Knoten, 6L^2 Bonds = 2N,
    je Zelle sum (e.delta)^2 = 3 ueber e=(1,0)). Jeder Knoten hat exakt
    Koordinationszahl z=4. Beides groessen-unabhaengig. Alle Bonds Laenge 1.
(b) SANDVIK-PAAR-ORAKEL (instant): kuenstliche WM-Daten mit bekanntem T_true
    -> die C-freie Paar-Schaetzung MUSS T_true rekonstruieren (Estimator-
    Verdrahtung gegen die wiederverwendeten PHY030-v02-Fit-Funktionen).
(c) BANDPLAUSIBILITAET (instant): die Kagome-Referenz 0.825 liegt strikt
    zwischen Honeycomb (0.573) und Triangular (1.418) - z=4 zwischen z=3/z=6.
(d) SMOKE (slow): Mini-Wolff-Lauf -> per-L NK-Crossing im physikalischen Band.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from phy033_kagome import (
    REF_KAGOME,
    SW_LIMIT_KAGOME,
    XYConfig,
    build_kagome_lattice,
    helicity_from_ensemble,
    helicity_terms,
    make_sandvik_pair_estimator,
    measure_kagome_seedwise,
    nk_crossing_kagome,
    wm_form,
)

TOL = 1e-6


# ---------------------------------------------------------------------------
# (a) Geometrie-Orakel
# ---------------------------------------------------------------------------

def _upsilon_T0(L: int) -> float:
    lat = build_kagome_lattice(L)
    t1, sa = helicity_terms(np.zeros(lat.n_nodes), lat, lat.pos)
    return helicity_from_ensemble([t1], [sa], XYConfig(J=1.0, T=1e-12),
                                  lat.n_nodes)


def test_kagome_T0_is_one():
    """Upsilon(T->0) MUSS exakt 1 J sein (per-Site-Spinwellen-Grenzwert, z=4)."""
    ups = _upsilon_T0(L=8)
    assert abs(ups - SW_LIMIT_KAGOME) < TOL, (
        f"Upsilon(0)={ups:.8f}, erwartet {SW_LIMIT_KAGOME} (per-Site)")


def test_kagome_T0_size_independent():
    """Per-Site-Upsilon(0) ist groessen-unabhaengig (= 1) ueber mehrere L."""
    vals = [_upsilon_T0(L) for L in (4, 6, 8, 12)]
    for v in vals:
        assert abs(v - SW_LIMIT_KAGOME) < TOL, vals
    assert max(vals) - min(vals) < TOL, vals


def test_kagome_coordination_is_four():
    """Jeder Knoten hat exakt 4 Nachbarn (Kagome z=4); N=3L^2, 6L^2 Bonds."""
    L = 6
    lat = build_kagome_lattice(L)
    assert lat.n_nodes == 3 * L * L
    assert len(lat.edges) == 6 * L * L  # = 2N
    deg = [0] * lat.n_nodes
    for i, j in lat.edges:
        deg[i] += 1
        deg[j] += 1
    assert all(d == 4 for d in deg), set(deg)


def test_kagome_bonds_are_unit_length():
    """Alle Bond-Verschiebungen haben Einheitslaenge (NN-Abstand a=1)."""
    lat = build_kagome_lattice(5)
    for dx, dy in lat.edge_disp:
        assert abs(math.hypot(dx, dy) - 1.0) < 1e-12, (dx, dy)


def test_kagome_no_self_loops_and_directions_match_pos():
    """Keine Selbstkanten; edge_disp ist der korrekte minimum-image-NN-Vektor.

    Auf dem Torus muss r_b - r_a (mod Box) genau edge_disp sein - sonst
    rechnet helicity_terms mit der falschen Bond-Projektion (stiller Bug).
    """
    L = 6
    lat = build_kagome_lattice(L)
    box_x = L * 2.0          # a1=(2,0) -> Breite 2L
    box_a2 = (L * 1.0, L * math.sqrt(3.0))  # a2-Periode
    for k, (a, b) in enumerate(lat.edges):
        assert a != b, (a, b)
        dx, dy = lat.edge_disp[k]
        rawx = lat.pos[b][0] - lat.pos[a][0]
        rawy = lat.pos[b][1] - lat.pos[a][1]
        # min-image gegen die beiden Gitterperioden (a1 entlang x; a2 schraeg)
        best = None
        for m in (-1, 0, 1):
            for nn in (-1, 0, 1):
                cx = rawx - m * box_x - nn * box_a2[0]
                cy = rawy - nn * box_a2[1]
                d = math.hypot(cx - dx, cy - dy)
                if best is None or d < best:
                    best = d
        assert best < 1e-9, (a, b, (dx, dy), (rawx, rawy), best)


def test_kagome_invalid_L_raises():
    """L < 1 muss hart fehlern statt still ein leeres Gitter zu liefern.

    Schliesst den von Equalita gemeldeten Silent-Edge (build_kagome_lattice(-3)
    -> 27 Knoten, 0 Bonds) gegen das Silent-Failure-Gate.
    """
    for bad in (0, -1, -3):
        with pytest.raises(ValueError):
            build_kagome_lattice(bad)


# ---------------------------------------------------------------------------
# (b) Sandvik-Paar-Orakel (Estimator-Verdrahtung, C-frei)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("T_true,C_true", [(0.825, 2.0), (0.80, 0.8),
                                           (0.85, 4.0)])
def test_sandvik_pair_recovers_synthetic_T_true(T_true, C_true):
    """C-freie (L1,L2)-Paar-Schaetzung rekonstruiert T_true aus WM-Daten.

    Reine WM-Form erfuellt die Paar-Bedingung exakt nur bei T_true; eine glatte
    T-Abhaengigkeit (slope*(T-T_true)*lnL) ausserhalb erzeugt den Nulldurchgang
    genau bei T_true. Mutations-/Vorzeichenfehler verfehlen ihn.
    """
    Ls = (12, 24)
    T_grid = [round(T_true - 0.02 + 0.005 * k, 4) for k in range(9)]
    slope = 0.3
    means_LT = {
        L: {T: wm_form(T_true, L, C_true) + slope * (T - T_true) * math.log(L)
            for T in T_grid}
        for L in Ls}
    est = make_sandvik_pair_estimator(Ls[0], Ls[1], T_grid)
    tc = est(means_LT)
    assert tc is not None, "kein Nulldurchgang der Paar-Bedingung"
    assert abs(tc - T_true) < 0.006, f"Paar T_BKT={tc:.4f} verfehlt {T_true}"


# ---------------------------------------------------------------------------
# (c) Bandplausibilitaet
# ---------------------------------------------------------------------------

def test_kagome_reference_between_honeycomb_and_triangular():
    """Kagome (z=4) liegt strikt zwischen Honeycomb (0.573) und Tri (1.418)."""
    assert abs(REF_KAGOME - 0.825) < 1e-9
    assert 0.573 < REF_KAGOME < 1.418


# ---------------------------------------------------------------------------
# (d) Smoke (slow): echtes Mini-Wolff-Sampling
# ---------------------------------------------------------------------------

# De-slow 2026-07-10 (Code-Audit 1.3): gemessen <3.5s - gehoert als
# fails-before-Gate in die schnelle CI-Suite, nicht hinter -m slow.
def test_kagome_crossing_in_physical_band():
    """Mini-Wolff-Lauf: per-L NK-Crossing im Band um die Referenz 0.825.

    Kleine L, wenige Sweeps - kein Praezisions-Anspruch, nur dass die Pipeline
    einen physikalisch sinnvollen Crossing liefert (finite-size von oben; klar
    zwischen Honeycomb ~0.57 und Triangular ~1.4).
    """
    by_T = {}
    for T in (0.80, 0.86, 0.92):
        by_T[T] = measure_kagome_seedwise(8, T, n_seeds=2, n_measure=120,
                                          n_burn=100, master_seed=42)
    tc = nk_crossing_kagome(by_T)
    assert tc is not None, "kein NK-Crossing im Scan-Fenster"
    assert 0.65 < tc < 1.05, f"Crossing T_x={tc:.4f} ausserhalb Band"
