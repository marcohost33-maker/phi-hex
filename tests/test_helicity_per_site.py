"""Helicity-Modulus per-Site-Normierung: Korrektheits-Gates.

Wissenschafts-Audit 2026-06-04 (S2/S3), Marco-Entscheid 2026-06-04:
Umstellung von FLAECHEN-Normierung (Upsilon(0)=J*sqrt(3)=1.732) auf
per-Site-Normierung (Upsilon(0)=1.5*J), Literatur-Standard
(arXiv:2406.12076, Nelson-Kosterlitz-Sprung Upsilon(T_BKT)=2 T_BKT/pi).

Unabhaengiges Orakel:
  T=0-Spinwellen-Grenzwert des XY-Dreiecksgitters (per-Site) = 3/2 * J.
  Bei T=0 sind alle Spins gleich -> cos(dtheta)=1, sin=0, daher
  Upsilon(0) = (J/N) * sum_<ij> (e . r_ij)^2. Auf dem Torus (z=6, a=1)
  ergibt das exakt 1.5*J, groessen-unabhaengig.
"""
from __future__ import annotations

import numpy as np

from phi_hex_core_v2 import (
    XYConfig,
    axial_to_xy,
    build_triangular_lattice,
    helicity_from_ensemble,
    helicity_terms,
)

J = 1.0
SW_LIMIT_TRIANGULAR = 1.5  # Spinwellen-Grenzwert per-Site, XY-Dreiecksgitter
TOL = 0.03


def _upsilon_T0(radius: int) -> float:
    """Upsilon(T->0) auf dem Dreiecks-Torus, exakt ausgerichteter Zustand."""
    lat = build_triangular_lattice(radius=radius, periodic=True)
    theta = np.zeros(lat.n_nodes)  # T=0: perfekt geordnet
    node_pos = [axial_to_xy(a) for a in lat.ax_of]
    t1, sin_acc = helicity_terms(theta, lat, node_pos)
    cfg = XYConfig(J=J, T=1e-9)  # T->0
    return helicity_from_ensemble([t1], [sin_acc], cfg, lat.n_nodes)


def test_helicity_T0_triangular_is_three_halves():
    """Υ(T→0) MUSS 1.5J sein (per-Site-Spinwellen-Grenzwert), tol 0.03."""
    ups = _upsilon_T0(radius=4)
    assert abs(ups - SW_LIMIT_TRIANGULAR) < TOL, (
        f"Upsilon(0)={ups:.5f}, erwartet {SW_LIMIT_TRIANGULAR} "
        f"(per-Site). Wert ~1.732 => noch FLAECHEN-normiert (Audit S2)."
    )


def test_helicity_T0_is_size_independent():
    """Bei T=0 ist Upsilon per-Site groessen-unabhaengig (= 1.5)."""
    vals = [_upsilon_T0(radius=r) for r in (3, 4, 5)]
    for v in vals:
        assert abs(v - SW_LIMIT_TRIANGULAR) < TOL, vals
    # untereinander praktisch identisch
    assert max(vals) - min(vals) < 1e-9, vals
