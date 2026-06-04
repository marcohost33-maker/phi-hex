"""Cross-Validierung core <-> PHY028 (Audit S3).

Vor dem Fix verglich die "Cross-Validierung" zwei VERSCHIEDENE Definitionen:
- core.helicity_from_ensemble teilte durch die FLAECHE (N*sqrt(3)/2),
- PHY028.square_helicity_from_ensemble teilte durch N (per-Site).
Damit trug der Vergleich nicht (Audit S3).

Nach dem Fix nutzen BEIDE Pfade dieselbe per-Site-Definition
Upsilon = (1/N)[J<T1> - (J^2/T)<sin^2>]. Diese Tests sichern das ab:

1. Identitaets-Test: derselbe Ensemble-Input -> identische Upsilon aus
   core- und PHY028-Ensemble-Funktion (assert_allclose, rtol 1e-12).
2. Konsistenz-Test auf demselben System: T=0-Spinwellen-Grenzwerte beider
   Pfade auf IHREM Gitter (Dreieck 1.5, Quadrat 1.0) - beide per-Site.
"""
from __future__ import annotations

import math

import numpy as np
from numpy.testing import assert_allclose

from phi_hex_core_v2 import (
    XYConfig,
    axial_to_xy,
    build_triangular_lattice,
    helicity_from_ensemble,
    helicity_terms,
)
from phy028_square_validation import (
    build_square_lattice,
    square_helicity_from_ensemble,
    square_helicity_terms,
)


def test_ensemble_normalisation_identical():
    """Selber (t1, sin)-Ensemble-Input -> identische Upsilon aus beiden
    Ensemble-Funktionen. Sichert, dass beide DIESELBE per-Site-Formel nutzen."""
    rng = np.random.default_rng(20260604)
    n_nodes = 100
    cfg = XYConfig(J=1.0, T=0.9)
    t1 = list(rng.normal(150.0, 5.0, 50))
    sin = list(rng.normal(0.0, 3.0, 50))

    ups_core = helicity_from_ensemble(t1, sin, cfg, n_nodes)
    ups_sq = square_helicity_from_ensemble(t1, sin, cfg, n_nodes)
    assert_allclose(ups_core, ups_sq, rtol=1e-12, atol=0.0)


def test_T0_spinwave_limits_both_paths_per_site():
    """Auf je IHREM Gitter liefern beide Pfade den per-Site-Spinwellen-Grenzwert:
    Dreieck = 1.5J, Quadrat = 1.0J. (Selbe Konvention, verschiedene Geometrie.)"""
    cfg = XYConfig(J=1.0, T=1e-9)

    # Dreieck (core)
    lat_t = build_triangular_lattice(radius=4, periodic=True)
    pos = [axial_to_xy(a) for a in lat_t.ax_of]
    t1_t, sin_t = helicity_terms(np.zeros(lat_t.n_nodes), lat_t, pos)
    ups_tri = helicity_from_ensemble([t1_t], [sin_t], cfg, lat_t.n_nodes)
    assert abs(ups_tri - 1.5) < 0.03, ups_tri

    # Quadrat (PHY028)
    lat_s = build_square_lattice(16)
    t1_s, sin_s = square_helicity_terms(np.zeros(lat_s.n_nodes), lat_s)
    ups_sq = square_helicity_from_ensemble([t1_s], [sin_s], cfg, lat_s.n_nodes)
    assert abs(ups_sq - 1.0) < 0.03, ups_sq


def test_core_consumes_square_ensemble_equivalently():
    """End-to-End: PHY028-Quadrat-Ensemble durch die CORE-Ensemble-Funktion
    gibt denselben Upsilon wie die PHY028-eigene Funktion (gemeinsame Defn)."""
    cfg = XYConfig(J=1.0, T=0.88)
    lat_s = build_square_lattice(12)
    rng = np.random.default_rng(7)
    t1s, sins = [], []
    for _ in range(30):
        theta = rng.uniform(0, 2 * math.pi, lat_s.n_nodes)
        a, b = square_helicity_terms(theta, lat_s)
        t1s.append(a)
        sins.append(b)
    ups_core = helicity_from_ensemble(t1s, sins, cfg, lat_s.n_nodes)
    ups_sq = square_helicity_from_ensemble(t1s, sins, cfg, lat_s.n_nodes)
    assert_allclose(ups_core, ups_sq, rtol=1e-12, atol=0.0)
