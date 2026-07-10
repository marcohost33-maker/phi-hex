"""Import- + Orakel-Smokes fuer die fuenf bis dahin ungetesteten src-Module
(Coverage-Schliessung, Code-Audit 2026-07-10).

Vorher liefen phy017, PHY027, PHY029 und die beiden phy024-Referenzen nur
durch compile+lint - ein kaputter Runtime-Import oder eine stille
NumPy-Inkompatibilitaet waere gruen gemerged. Jede Datei bekommt hier
mindestens ein billiges, MC-freies fails-before-Orakel.
"""
from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np

import conftest

phy017 = conftest.load_all()["phy017_validation_kernel"]
phy027 = conftest.load_all()["phy027_precision_fss"]
phy029 = conftest.load_all()["phy029_triangular_sandvik"]
phy024_r0 = conftest.load_all()["phy024_r0_reference"]
phy024_r2 = conftest.load_all()["phy024_r2_honeycomb"]


# ---------------------------------------------------------------------------
# phy017: C6-Ring-Spektrum (exakt analytisch loesbar)
# ---------------------------------------------------------------------------

def test_phy017_c6_spectrum_matches_analytic():
    """Tight-Binding-C_N-Ring mit Peierls-Phase: numerische Eigenwerte ==
    analytisches Spektrum E_k = -2t cos(2 pi (k + phi)/N)."""
    for phi in (0.0, 0.17, 0.5):
        H = phy017.tight_binding_hamiltonian_c6(6, 1.0, phi)
        assert np.allclose(H, H.conj().T), "Hamiltonian nicht hermitesch"
        num = np.sort(np.linalg.eigvalsh(H))
        ana = np.sort(phy017.analytic_spectrum_c6(6, 1.0, phi))
        assert np.allclose(num, ana, atol=1e-12), phi


def test_phy017_persistent_current_zero_at_zero_flux():
    """I(phi=0) = 0 (Zeitumkehr-Symmetrie des halbgefuellten C6-Rings)."""
    val = phy017.persistent_current_ground_state(6, 3, 1.0, 0.0)
    assert abs(val) < 1e-6


# ---------------------------------------------------------------------------
# PHY027: NK-Crossing-Detektion (synthetisch, kein MC)
# ---------------------------------------------------------------------------

def test_phy027_find_crossing_linear_oracle():
    """Upsilon(T) linear fallend crosst die NK-Gerade 2T/pi an bekanntem T*.

    Konstruktion: Upsilon(T) = (2 T*/pi) + a*(T* - T) mit a > 2/pi
    -> Crossing exakt bei T*.
    """
    T_star, a = 1.45, 1.5
    ms = []
    for T in (1.30, 1.40, 1.50, 1.60, 1.70):
        ups = 2.0 * T_star / math.pi + a * (T_star - T)
        ms.append(SimpleNamespace(T=T, upsilon_mean=ups,
                                  nk_line=2.0 * T / math.pi))
    tc = phy027.find_crossing(ms)
    assert tc is not None and abs(tc - T_star) < 1e-9, tc


def test_phy027_find_crossing_none_without_sign_change():
    """Ohne Crossing im Raster: None (kein erfundener Schaetzer)."""
    ms = [SimpleNamespace(T=T, upsilon_mean=2.0, nk_line=2.0 * T / math.pi)
          for T in (0.5, 0.6, 0.7)]
    assert phy027.find_crossing(ms) is None


# ---------------------------------------------------------------------------
# PHY029: R-Transformation (Kern der Sandvik-Bedingung)
# ---------------------------------------------------------------------------

def test_phy029_R_of_zero_on_nk_line():
    """R = pi*Upsilon/(2T) - 1 verschwindet exakt AUF der NK-Geraden."""
    for T in (0.8, 1.4):
        assert abs(phy029.R_of(2.0 * T / math.pi, T)) < 1e-12
    assert phy029.R_of(1.5, 1.0) > 0.0   # unterhalb T_BKT positiv


# ---------------------------------------------------------------------------
# phy024 R0: endlicher Hex-Patch (offene Raender)
# ---------------------------------------------------------------------------

def test_phy024_r0_hex_patch_geometry_and_aligned_oracles():
    """Radius-Konvention (19 Knoten bei r=2), wrap_angle-Bereich, aligned:
    keine Windungen, Kohaerenz exakt 1."""
    patch = phy024_r0.build_hex_patch(2)
    assert len(patch.coords) == 19
    assert phy024_r0.wrap_angle(3 * math.pi) in (-math.pi, math.pi)
    assert abs(phy024_r0.wrap_angle(0.3) - 0.3) < 1e-12
    theta0 = np.zeros(len(patch.coords))
    w = phy024_r0.plaquette_windings(theta0, patch.plaquette_rings)
    assert (w == 0).all()
    coh = phy024_r0.helicity_analog(theta0[None, :], patch.edges)
    assert abs(coh - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# phy024 R2: periodischer Honeycomb-Torus
# ---------------------------------------------------------------------------

def test_phy024_r2_torus_geometry():
    """z=3-Torus: N = 2 lx ly Knoten, 3 lx ly Kanten, lx ly Plaquetten,
    jeder Knoten Grad 3."""
    t = phy024_r2.build_honeycomb_torus(3, 4)
    assert t.n == 2 * 3 * 4
    assert len(t.edges) == 3 * 3 * 4
    assert len(t.plaquettes) == 3 * 4
    assert all(len(nb) == 3 for nb in t.neighbors)


def test_phy024_r2_aligned_energy_and_twist_curvature():
    """Aligned: E = -J * #Kanten exakt; Twist-Kruemmung ueber die
    x-Rand-Kanten = J * #Randkanten / N (analytisch, da nur die getwisteten
    Bonds cos(phi)-abhaengig sind)."""
    t = phy024_r2.build_honeycomb_torus(3, 4)
    theta0 = np.zeros(t.n)
    e = phy024_r2.energy(theta0, t.edges, 1.0)
    assert abs(e - (-1.0 * len(t.edges))) < 1e-12
    assert len(t.x_boundary_edges) > 0
    curv = phy024_r2.twist_energy_curvature(
        theta0[None, :], t.edges, t.x_boundary_edges, 1.0)
    expected = 1.0 * len(t.x_boundary_edges) / t.n
    assert abs(curv - expected) < 1e-5, (curv, expected)
