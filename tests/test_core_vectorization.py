"""Aequivalenz-Gates fuer die vektorisierten Kern-Hotspots (Optimierung 2026-06).

Die Vektorisierung von `helicity_terms` und `cliff_delta` ist eine reine
Performance-Optimierung - sie darf die Physik NICHT veraendern. Diese Tests
pinnen das fest, indem sie die vektorisierten Implementierungen gegen eine
unabhaengige, in-test gehaltene SKALARE Referenz (die fruehe Schleifen-Form)
vergleichen:

- helicity_terms: identisch bis auf Gleitkomma-Summationsreihenfolge
  (rtol 1e-10) - eine echte Mutation (falsches Vorzeichen/Index) bricht das.
- cliff_delta: BIT-IDENTISCH (nur ganzzahlige Paar-Vergleiche).

Zusaetzlich: Gitter-Agnostik (Triangular- UND Honeycomb-Gitter laufen durch
denselben Kern) und Cache-Konsistenz.
"""
from __future__ import annotations

import math

import numpy as np
from numpy.testing import assert_allclose

from phi_hex_core_v2 import (
    axial_to_xy,
    build_triangular_lattice,
    cliff_delta,
    helicity_terms,
)
from phy031_honeycomb import build_honeycomb_lattice


def _helicity_terms_scalar(theta, lattice, pos, direction=(1.0, 0.0)):
    """Unabhaengige skalare Referenz (die fruehe Schleifen-Form, byte-treu)."""
    ex, ey = direction
    norm = math.hypot(ex, ey)
    ex, ey = ex / norm, ey / norm
    t1 = 0.0
    sin_accum = 0.0
    use_disp = len(lattice.edge_disp) == len(lattice.edges)
    for k, (i, j) in enumerate(lattice.edges):
        if use_disp:
            dx, dy = lattice.edge_disp[k]
        else:
            dx = pos[j][0] - pos[i][0]
            dy = pos[j][1] - pos[i][1]
        proj = ex * dx + ey * dy
        dtheta = theta[i] - theta[j]
        t1 += math.cos(dtheta) * proj * proj
        sin_accum += math.sin(dtheta) * proj
    return t1, sin_accum


def _random_theta(n, seed):
    return np.random.default_rng(seed).uniform(0.0, 2.0 * math.pi, n)


def test_helicity_terms_matches_scalar_triangular():
    """Vektorisiert == skalar (rtol 1e-10) auf dem Dreiecks-Torus, mehrere T/Seeds."""
    lat = build_triangular_lattice(radius=4, periodic=True)
    pos = [axial_to_xy(a) for a in lat.ax_of]
    for seed in (0, 1, 7, 42):
        theta = _random_theta(lat.n_nodes, seed)
        t1_v, sa_v = helicity_terms(theta, lat, pos)
        t1_s, sa_s = _helicity_terms_scalar(theta, lat, pos)
        assert_allclose(t1_v, t1_s, rtol=1e-10, atol=1e-12)
        assert_allclose(sa_v, sa_s, rtol=1e-10, atol=1e-12)


def test_helicity_terms_matches_scalar_honeycomb():
    """Gitter-agnostik: derselbe Kern liefert auf Honeycomb dasselbe wie skalar."""
    lat = build_honeycomb_lattice(6)
    for seed in (3, 11):
        theta = _random_theta(lat.n_nodes, seed)
        t1_v, sa_v = helicity_terms(theta, lat, lat.pos)
        t1_s, sa_s = _helicity_terms_scalar(theta, lat, lat.pos)
        assert_allclose(t1_v, t1_s, rtol=1e-10, atol=1e-12)
        assert_allclose(sa_v, sa_s, rtol=1e-10, atol=1e-12)


def test_helicity_terms_direction_dependence():
    """Nicht-achsenparallele Messrichtung wird korrekt projiziert (vs skalar)."""
    lat = build_triangular_lattice(radius=3, periodic=True)
    pos = [axial_to_xy(a) for a in lat.ax_of]
    theta = _random_theta(lat.n_nodes, 5)
    for direction in ((0.0, 1.0), (1.0, 1.0), (2.0, -0.5)):
        v = helicity_terms(theta, lat, pos, direction=direction)
        s = _helicity_terms_scalar(theta, lat, pos, direction=direction)
        assert_allclose(v, s, rtol=1e-10, atol=1e-12)


def test_helicity_terms_cache_consistency():
    """Wiederholte Aufrufe (Cache aktiv) liefern bit-identische Ergebnisse."""
    lat = build_triangular_lattice(radius=3, periodic=True)
    pos = [axial_to_xy(a) for a in lat.ax_of]
    theta = _random_theta(lat.n_nodes, 9)
    first = helicity_terms(theta, lat, pos)
    assert hasattr(lat, "_phx_edge_cache")  # Cache wurde gesetzt
    for _ in range(3):
        assert helicity_terms(theta, lat, pos) == first


def _cliff_delta_scalar(a, b):
    na, nb = len(a), len(b)
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (na * nb)


def test_cliff_delta_bit_identical_to_scalar():
    """cliff_delta ist BIT-identisch zur skalaren Doppelschleife (Integer-Counts)."""
    rng = np.random.default_rng(2026)
    for _ in range(20):
        a = rng.normal(0.3, 1.0, rng.integers(5, 40))
        b = rng.normal(0.0, 1.0, rng.integers(5, 40))
        assert cliff_delta(a, b)["cliff_delta"] == _cliff_delta_scalar(a, b)


def test_cliff_delta_ties_and_extremes():
    """Gleichstaende und disjunkte Mengen: Vorzeichen/Betrag korrekt."""
    assert cliff_delta(np.array([1.0, 1.0]), np.array([1.0, 1.0]))[
        "cliff_delta"] == 0.0
    assert cliff_delta(np.array([5.0, 6.0]), np.array([1.0, 2.0]))[
        "cliff_delta"] == 1.0
    assert cliff_delta(np.array([1.0, 2.0]), np.array([5.0, 6.0]))[
        "cliff_delta"] == -1.0
