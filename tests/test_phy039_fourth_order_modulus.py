"""Schnelle Korrektheits-Gates fuer PHY039 (4.-Ordnungs-Helicity-Modul).

Kern-Gate ist das numerische ORAKEL: der erst-prinzipiell hergeleitete
F4-Estimator muss die numerische 4. Ableitung der EXAKTEN freien Energie
F(Delta) auf einer fixen Stichprobe auf Maschinengenauigkeit reproduzieren -
das ist die Korrektheits-Garantie unabhaengig von jeder (blockierten)
Literatur-Formel. Alle Gates sind MC-frei und CI-schnell.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

import conftest

phy039 = conftest._load(
    "phy039", "260616 PHY039 fourth-order helicity modulus v01.py")


def test_f4_estimator_numeric_oracle():
    """F4 (und F2) closed-form == numerische 4.(2.) Ableitung von F(Delta)."""
    o = phy039.f4_estimator_oracle()
    assert o["F2_relerr"] < 1e-6, o
    assert o["F4_relerr"] < 1e-3, o          # 4. Abl. -> Polyfit-limitiert
    # nicht-trivial (sonst sagt rel-err nichts):
    assert abs(o["F4_closed"]) > 1e-2, o


def test_oracle_seed_robust():
    """Das Orakel haelt fuer mehrere Stichproben-Seeds (keine Zufallstreffer)."""
    for seed in (1, 11, 123):
        o = phy039.f4_estimator_oracle(seed=seed)
        assert o["F4_relerr"] < 2e-3, (seed, o)


def test_reduced_equals_full_on_symmetric_sample():
    """Auf einer theta->-theta-SYMMETRISCHEN Stichprobe (odd-Mittel = 0 exakt)
    fallen volle und reduzierte 4.-Ordnungs-Form zusammen - die Reduktion ist
    also genau die thermische-Symmetrie-Annahme, nicht eine Naeherung."""
    rng = np.random.default_rng(3)
    ei = np.array([0, 1, 2, 3])
    ej = np.array([1, 2, 3, 0])
    a = np.array([1.0, 0.7, 1.0, 0.4])
    base = rng.uniform(0, 2 * math.pi, size=(200, 4))
    configs = np.vstack([base, -base])        # symmetrisch unter theta->-theta
    agg = np.array([phy039.bond_aggregates(c, ei, ej, a) for c in configs])
    C, S, S3, C4 = agg[:, 0], agg[:, 1], agg[:, 2], agg[:, 3]
    beta = 1.3
    full = phy039.fourth_order_F(C, S, S3, C4, beta)

    def m(v):
        return float(np.mean(v))
    reduced = (-m(C4) + beta * (4 * m(S * S3) - 3 * m(C * C) + 3 * m(C) ** 2)
               + 6 * beta ** 2 * (m(C * S * S) - m(C) * m(S * S))
               - beta ** 3 * (m(S ** 4) - 3 * m(S * S) ** 2))
    assert abs(full - reduced) < 1e-9, (full, reduced)


def test_aligned_state_y2_is_one():
    """Perfekt ausgerichtet (T->0): Upsilon_2 = 1.0 fuer das Quadratgitter
    (z=4, per-Site) - Wiring der Bond-Aggregate gegen die Helicity-Konvention."""
    lat = phy039.build_square_lattice(8)
    edges = np.asarray(lat.edges)
    a = np.array([dx for dx, dy in lat.edge_disp])
    c0, s0, _, _ = phy039.bond_aggregates(
        np.zeros(lat.n_nodes), edges[:, 0], edges[:, 1], a)
    y2 = phy039.second_order_F([c0], [s0], 1.0 / 0.01, reduced=True) / lat.n_nodes
    assert abs(y2 - 1.0) < 1e-9


def test_dip_minimum_parabolic_vertex():
    """dip_minimum liefert den parabolischen Scheitel eines synthetischen Dips
    und meldet Rand-Minima ehrlich als None."""
    from dataclasses import make_dataclass
    R = make_dataclass("R", [("y4_scaled_mean", float)])
    # Dip mit Minimum bei T=0.90 (Scheitel zwischen 0.88/0.90/0.92):
    vals = {0.86: -100, 0.88: -300, 0.90: -420, 0.92: -300, 0.94: -100}
    by_T = {T: R(v) for T, v in vals.items()}
    tmin, depth = phy039.dip_minimum(by_T)
    assert tmin is not None and abs(tmin - 0.90) < 1e-6
    assert depth == -420
    # monoton fallend -> Rand-Minimum -> None
    mono = {0.86: -100, 0.88: -200, 0.90: -300, 0.92: -400, 0.94: -500}
    tmin2, _ = phy039.dip_minimum({T: R(v) for T, v in mono.items()})
    assert tmin2 is None


@pytest.mark.slow
def test_square_y2_recovers_reference_signal():
    """Langsam: ein einzelner (L=16) Messpunkt - Upsilon_2 nahe der NK-Geraden
    am Goldstandard, und N*Upsilon_4 ist negativ (Dip-Vorzeichen)."""
    r = phy039.measure_square_y2y4(16, 0.91, n_seeds=6, n_measure=500, n_burn=300)
    assert r.y2_mean > 0.3              # physikalisch (renormiert), positiv
    assert r.y4_scaled_mean < 0.0       # 4.-Ordnungs-Dip ist negativ
