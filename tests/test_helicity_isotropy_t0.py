"""T=0-Helicity-Modul-ISOTROPIE: gitter-uebergreifendes Korrektheits-Gate.

Coworker Research / Coworkerz, 2026-06-19 (Inkrement 2).

MOTIVATION (Luecke vor diesem Test):
  Der per-Site-T=0-Spinwellen-Grenzwert wird bereits geprueft - aber NUR fuer
  die Mess-Richtung e=(1,0) (test_helicity_per_site / test_tbkt_honeycomb /
  test_tbkt_kagome). Die staerkere, bisher UNGEPRUEFTE Eigenschaft ist die
  ISOTROPIE des Helicity-Modul-Tensors bei T=0:

      Upsilon_{ab}(0) = (J/N) * sum_<ij> (delta_ij)_a (delta_ij)_b .

  Fuer ein Bond-Set mit >=3-zaehliger Punktsymmetrie (Dreieck z=6, Honeycomb
  z=3, Kagome z=4) ist diese 2x2-Matrix proportional zur Einheitsmatrix
  (Neumann-Prinzip): Upsilon_{ab}(0) = Upsilon0 * delta_{ab}. Damit ist der
  gemessene Helicity-Modul Upsilon(0, theta) = e(theta)^T Upsilon e(theta)
  RICHTUNGS-UNABHAENGIG (= Upsilon0 fuer jeden Winkel) und die Off-Diagonale
  Upsilon_xy verschwindet.

UNABHAENGIGES ORAKEL (kein Monte-Carlo, Maschinengenauigkeit):
  1. Konstanz ueber einen Winkel-Scan ist algebraisch AEQUIVALENT zu
     Upsilon_xx == Upsilon_yy UND Upsilon_xy == 0 (Fourier-Argument:
     Upsilon(theta) = Upsilon_xx cos^2 + Upsilon_yy sin^2 + 2 Upsilon_xy
     sin cos = const fuer alle theta <=> Diagonal & gleich).
  2. Der konstante Wert MUSS der jeweilige analytische Grenzwert sein
     (Dreieck 3/2, Honeycomb 3/4, Kagome 1).

NICHT-VAKUOS: ein anisotropes Bond-Set (Kette, Bonds nur entlang x) liefert
  Upsilon(theta) = cos^2(theta) - der Scan-Spread ist 1.0, die Off-Diagonale
  ist 0 nur zufaellig, aber Upsilon_xx != Upsilon_yy -> das Gate FAILt. Damit
  faengt das Gate eine reale Regression in der Richtungs-Projektion von
  helicity_terms (edge_disp / Normierung), die alle bestehenden e=(1,0)-Tests
  durchwinken wuerden.

Dies haertet die Kern-Observable helicity_terms gitter-agnostisch und hebt das
in PHY031/PHY033 nur INNERHALB der langsamen Wolff-Laeufe vorhandene T=0-Orakel
auf eine schnelle, richtungs-vollstaendige CI-Ebene.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from phi_hex_core_v2 import (
    XYConfig,
    axial_to_xy,
    build_triangular_lattice,
    helicity_from_ensemble,
    helicity_terms,
)
from phy031_honeycomb import build_honeycomb_lattice
from phy033_kagome import build_kagome_lattice

# Maschinengenauigkeits-Toleranzen: T=0 ist exakt (keine Statistik).
TOL_CONST = 1e-9   # Scan-Konstanz (Isotropie)
TOL_VALUE = 1e-9   # Abstand zum analytischen Grenzwert
TOL_OFFDIAG = 1e-9  # Off-Diagonale Upsilon_xy

# Winkel-Scan inkl. nicht-symmetrie-ausgerichteter Richtungen (137 deg),
# damit ein versteckter Off-Diagonal-Term sicher sichtbar wuerde.
SCAN_ANGLES_DEG = (0.0, 17.0, 30.0, 45.0, 60.0, 90.0, 137.0)


def _upsilon_T0(lattice, pos, direction: tuple[float, float]) -> float:
    """Per-Site-Helicity bei T=0 (perfekt ausgerichteter Zustand)."""
    t1, sa = helicity_terms(np.zeros(lattice.n_nodes), lattice, pos,
                            direction=direction)
    return helicity_from_ensemble([t1], [sa], XYConfig(J=1.0, T=1e-12),
                                  lattice.n_nodes)


def _upsilon_tensor_T0(lattice, pos) -> tuple[float, float, float]:
    """Liefert (Upsilon_xx, Upsilon_yy, Upsilon_xy) bei T=0.

    Bei T=0 (theta == 0) ist sin_accum == 0, also Upsilon = (J/N) <T1> und
    T1 = sum_<ij> (e.delta)^2 ist eine reine quadratische Form in e. Damit
    extrahieren wir die Tensor-Komponenten direkt aus drei Richtungen:
      e=(1,0) -> Upsilon_xx ;  e=(0,1) -> Upsilon_yy ;
      e=(1,1)/sqrt2 -> (Upsilon_xx + Upsilon_yy)/2 + Upsilon_xy .
    """
    uxx = _upsilon_T0(lattice, pos, (1.0, 0.0))
    uyy = _upsilon_T0(lattice, pos, (0.0, 1.0))
    udiag = _upsilon_T0(lattice, pos, (1.0, 1.0))  # direction wird normiert
    uxy = udiag - 0.5 * (uxx + uyy)
    return uxx, uyy, uxy


# Gitter-Fabriken + analytischer T=0-Grenzwert je Gitter.
def _triangular(radius: int = 4):
    lat = build_triangular_lattice(radius=radius, periodic=True)
    pos = [axial_to_xy(a) for a in lat.ax_of]
    return lat, pos


_LATTICES = {
    "triangular": (lambda: _triangular(4), 1.5),
    "honeycomb": (lambda: (build_honeycomb_lattice(8),
                           build_honeycomb_lattice(8).pos), 0.75),
    "kagome": (lambda: (build_kagome_lattice(8),
                        build_kagome_lattice(8).pos), 1.0),
}


@pytest.mark.parametrize("name", sorted(_LATTICES))
def test_t0_helicity_is_direction_independent(name):
    """Upsilon(0, theta) ist richtungs-unabhaengig (Isotropie-Orakel)."""
    factory, expected = _LATTICES[name]
    lat, pos = factory()
    vals = [_upsilon_T0(lat, pos,
                        (math.cos(math.radians(a)), math.sin(math.radians(a))))
            for a in SCAN_ANGLES_DEG]
    spread = max(vals) - min(vals)
    assert spread < TOL_CONST, (
        f"{name}: Upsilon(0) ist NICHT richtungs-unabhaengig "
        f"(Scan-Spread={spread:.3e} ueber {SCAN_ANGLES_DEG} deg): {vals}")
    for a, v in zip(SCAN_ANGLES_DEG, vals):
        assert abs(v - expected) < TOL_VALUE, (
            f"{name}: Upsilon(0, {a} deg)={v:.10f}, erwartet {expected}")


@pytest.mark.parametrize("name", sorted(_LATTICES))
def test_t0_helicity_tensor_is_diagonal_isotropic(name):
    """Der T=0-Helicity-Tensor ist diagonal & isotrop: Uxx==Uyy, Uxy==0."""
    factory, expected = _LATTICES[name]
    lat, pos = factory()
    uxx, uyy, uxy = _upsilon_tensor_T0(lat, pos)
    assert abs(uxx - uyy) < TOL_CONST, (
        f"{name}: anisotrop, Uxx={uxx:.10f} != Uyy={uyy:.10f}")
    assert abs(uxy) < TOL_OFFDIAG, (
        f"{name}: Off-Diagonale Uxy={uxy:.3e} != 0 (Tensor nicht diagonal)")
    assert abs(uxx - expected) < TOL_VALUE, (
        f"{name}: Uxx={uxx:.10f}, erwartet {expected}")


@pytest.mark.parametrize("name", sorted(_LATTICES))
def test_t0_isotropy_size_independent(name):
    """Isotropie + Wert sind groessen-unabhaengig (mehrere L bzw. radius)."""
    expected = _LATTICES[name][1]
    if name == "triangular":
        builders = [lambda r=r: _triangular(r) for r in (3, 4, 5)]
    elif name == "honeycomb":
        builders = [lambda L=L: (build_honeycomb_lattice(L),
                                 build_honeycomb_lattice(L).pos)
                    for L in (4, 6, 8)]
    else:
        builders = [lambda L=L: (build_kagome_lattice(L),
                                 build_kagome_lattice(L).pos)
                    for L in (3, 4, 6)]
    for build in builders:
        lat, pos = build()
        uxx, uyy, uxy = _upsilon_tensor_T0(lat, pos)
        assert abs(uxx - uyy) < TOL_CONST and abs(uxy) < TOL_OFFDIAG, (uxx, uyy, uxy)
        assert abs(uxx - expected) < TOL_VALUE, (name, uxx, expected)


# ---------------------------------------------------------------------------
# Nicht-vakuoses Negativ-Gate: ein anisotropes Bond-Set MUSS durchfallen.
# ---------------------------------------------------------------------------

class _ToyLattice:
    """Minimal-Gitter (nur edges + edge_disp + pos + n_nodes) fuer helicity_terms."""

    def __init__(self, edges, edge_disp, n_nodes):
        self.edges = list(edges)
        self.edge_disp = list(edge_disp)
        self.n_nodes = n_nodes
        self.pos = [(0.0, 0.0)] * n_nodes


def _chain_x(n: int = 10) -> _ToyLattice:
    """Anisotrope Kette: alle Bonds entlang +x (Einheitslaenge)."""
    edges = [(i, (i + 1) % n) for i in range(n)]
    disp = [(1.0, 0.0)] * n
    return _ToyLattice(edges, disp, n)


def test_anisotropic_chain_breaks_isotropy_gate():
    """Negativ-Kontrolle: die x-Kette ist anisotrop -> Isotropie-Gate FAILt.

    Beweist, dass die Konstanz-/Diagonal-Gates NICHT vakuos sind: bei einem
    anisotropen Bond-Set ist Upsilon(theta) = cos^2(theta), Scan-Spread = 1,
    Uxx=1 != Uyy=0. Ein Bug, der die Richtungs-Projektion kollabieren liesse,
    wuerde hier sichtbar - und an den echten Gittern oben einen FAIL erzeugen.
    """
    chain = _chain_x(10)
    vals = [_upsilon_T0(chain, chain.pos,
                        (math.cos(math.radians(a)), math.sin(math.radians(a))))
            for a in SCAN_ANGLES_DEG]
    spread = max(vals) - min(vals)
    # Das echte Gate ist (spread < TOL_CONST) -> hier MUSS es verletzt sein.
    assert spread > 0.5, f"anisotrope Kette unerwartet isotrop: {vals}"
    uxx, uyy, uxy = _upsilon_tensor_T0(chain, chain.pos)
    assert abs(uxx - 1.0) < 1e-12 and abs(uyy) < 1e-12, (uxx, uyy)
    assert abs(uxx - uyy) > 0.5, "Tensor faelschlich isotrop"
