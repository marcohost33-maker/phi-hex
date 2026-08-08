"""Tests fuer PHY043 (konventionsfreier Quercheck triangular, Audit O1).

Drei Ebenen (Vertragsquelle: spec/260808 PHI HEX phy043 triangular
convention-free crossing method v01.md):

(a) EXAKTE ORAKEL (fails-before-faehig, kein Monte-Carlo):
    - iid-U4: fuer N iid-uniforme Winkel gilt EXAKT E[|m|^2] = 1/N und
      E[|m|^4] = (2N-1)/N^3, also U4 -> 2 - 1/N.
    - Moden-Buchhaltung: uniforme Konfiguration -> S(0) = N, S(k1) = 0;
      Spin-Welle theta_j = 2 pi q_j / L -> S(0) = 0 und genau EINE der
      drei Minimalmoden traegt N (Mittel N/3).
(b) VERTRAGS-GATES der puren Finder: adjacent-Crossing (keine
    Interpolation ueber None-Luecken), Splay-Persistenz (ein einzelner
    2-sigma-Ausreisser genuegt nicht), fails-closed bei None.
(c) SMOKE: Mini-MC (L=5/7, wenige Sweeps) end-to-end inkl. Report-Struktur.
"""
from __future__ import annotations

import math

import numpy as np

import conftest

phy043 = conftest._load(
    "phy043_convention_crossing",
    "260808 PHY043 triangular convention-free crossing v01.py")


# ---------------------------------------------------------------------------
# (a) Exakte Orakel
# ---------------------------------------------------------------------------

def test_binder_u4_iid_oracle() -> None:
    """iid-uniforme Winkel: U4 = 2 - 1/N exakt im Erwartungswert.

    Herleitung: m = (1/N) sum exp(i theta_j) mit iid-uniformen theta_j.
    E|m|^2 = 1/N; E|m|^4 = (2N(N-1) + N)/N^4 = (2N-1)/N^3.
    Statistischer Test mit festem Seed und 3-sigma-artiger Toleranz.
    """
    rng = np.random.default_rng(20260808)
    N, n_samples = 64, 6000
    m2, m4 = [], []
    for _ in range(n_samples):
        theta = rng.uniform(0.0, 2.0 * math.pi, N)
        m = np.mean(np.exp(1j * theta))
        a = float(np.abs(m)) ** 2
        m2.append(a)
        m4.append(a * a)
    u4 = phy043.binder_u4(np.array(m2), np.array(m4))
    expected = 2.0 - 1.0 / N
    assert abs(u4 - expected) < 0.06, (u4, expected)


def _torus_axes(L: int):
    lattice = phy043.build_triangular_lattice(radius=(L - 1) // 2,
                                             periodic=True)
    q = np.array([a[0] for a in lattice.ax_of], dtype=float)
    r = np.array([a[1] for a in lattice.ax_of], dtype=float)
    return lattice, q, r


def test_structure_factor_uniform_config() -> None:
    """Uniforme Winkel: S(0) = N exakt, S(k1) = 0 (Geometrie-Summe)."""
    L = 9
    lattice, q, r = _torus_axes(L)
    n = lattice.n_nodes
    theta = np.full(n, 0.7321)
    s0, sk1 = phy043.structure_factors(theta, q, r, L)
    assert abs(s0 - n) < 1e-9 * n
    assert abs(sk1) < 1e-9 * n
    # xi-Schaetzer ist hier undefiniert (S(k1)=0) -> fails-closed None.
    assert phy043.xi_over_L_from_S(s0, sk1) is None


def test_structure_factor_spin_wave_mode_bookkeeping() -> None:
    """Spin-Welle theta_j = 2 pi q_j / L liegt exakt auf Mode (1,0).

    Dann ist S(0) = 0 und von den drei Minimalmoden traegt genau (1,0)
    das volle Gewicht N -> Mittel N/3. Prueft die Moden-Buchhaltung
    (Vorzeichen/Norm der torus-periodischen Wellen) exakt.
    """
    L = 9
    lattice, q, r = _torus_axes(L)
    n = lattice.n_nodes
    theta = 2.0 * math.pi * q / L
    s0, sk1 = phy043.structure_factors(theta, q, r, L)
    assert abs(s0) < 1e-9 * n
    assert abs(sk1 - n / 3.0) < 1e-9 * n
    # S(0) < S(k1): kein valider xi-Schaetzer -> fails-closed None.
    assert phy043.xi_over_L_from_S(s0, sk1) is None


def test_xi_over_L_prefactor() -> None:
    """Vorfaktor sqrt(3)/(4 pi) und Monotonie des Schaetzers."""
    val = phy043.xi_over_L_from_S(2.0, 1.0)
    assert abs(val - math.sqrt(3.0) / (4.0 * math.pi)) < 1e-12
    assert phy043.xi_over_L_from_S(5.0, 1.0) > val


# ---------------------------------------------------------------------------
# (b) Vertrags-Gates der puren Finder
# ---------------------------------------------------------------------------

def test_adjacent_crossing_linear_oracle() -> None:
    T = [1.0, 1.1, 1.2, 1.3]
    d = [-0.02, -0.01, 0.01, 0.03]
    sig = [0.001] * 4
    out = phy043.adjacent_crossings(T, d, sig)
    assert len(out) == 1
    assert abs(out[0]["T_cross"] - 1.15) < 1e-12
    assert out[0]["significance"] >= 10.0


def test_adjacent_crossing_refuses_gap_bridging() -> None:
    """Vorzeichenwechsel QUER ueber eine None-Luecke ist KEIN Crossing
    (gleicher Vertrag wie der M1/P2-gehaertete Root-Finder-Stack)."""
    T = [1.0, 1.1, 1.2, 1.3]
    d = [-0.02, None, None, 0.03]
    sig = [0.001] * 4
    assert phy043.adjacent_crossings(T, d, sig) == []


def test_splay_requires_persistence() -> None:
    """Ein einzelner signifikanter Punkt genuegt nicht; erst ab dem Punkt,
    ab dem ALLE folgenden signifikant im Splay-Vorzeichen liegen."""
    T = [1.0, 1.1, 1.2, 1.3, 1.4]
    sig = [0.01] * 5
    # Ausreisser bei 1.1 (signifikant), aber 1.2 faellt zurueck:
    d = [0.001, -0.05, -0.01, -0.06, -0.08]
    t_splay = phy043.splay_temperature(T, d, sig, sign=-1.0)
    assert t_splay is not None
    assert abs(t_splay - 1.3) < 1e-12
    # Kein persistenter Bereich -> None:
    d2 = [0.001, -0.05, 0.01, -0.06, 0.01]
    assert phy043.splay_temperature(T, d2, sig, sign=-1.0) is None


def test_splay_none_breaks_persistence() -> None:
    """None-Punkte brechen die Persistenz (fails-closed)."""
    T = [1.0, 1.1, 1.2]
    d = [-0.05, None, -0.06]
    sig = [0.001, 0.001, 0.001]
    assert phy043.splay_temperature(T, d, sig, sign=-1.0) is None


def test_pair_analysis_and_bootstrap_deterministic() -> None:
    """pair_analysis + Bootstrap auf synthetischen Kurven: Splay wird
    gefunden, Bootstrap-CI umschliesst ihn, Determinismus via Stream."""
    T = [1.0, 1.1, 1.2, 1.3, 1.4]
    o1 = [0.70, 0.69, 0.68, 0.66, 0.64]
    o2 = [0.70, 0.69, 0.62, 0.55, 0.48]   # splayt ab 1.2 nach unten
    s1 = [0.005] * 5
    s2 = [0.005] * 5
    res = phy043.pair_analysis(T, o1, s1, o2, s2, splay_sign=-1.0)
    assert res["T_splay"] is not None and abs(res["T_splay"] - 1.2) < 1e-12
    rng_a = phy043.make_rng(42, stream=phy043.BOOT_STREAM)
    rng_b = phy043.make_rng(42, stream=phy043.BOOT_STREAM)
    ci_a = phy043.bootstrap_splay_ci(T, o1, s1, o2, s2, -1.0, rng_a,
                                     n_boot=100)
    ci_b = phy043.bootstrap_splay_ci(T, o1, s1, o2, s2, -1.0, rng_b,
                                     n_boot=100)
    assert ci_a == ci_b, "Bootstrap nicht deterministisch"
    assert ci_a["none_frac"] < 0.5
    assert ci_a["ci_low"] <= res["T_splay"] <= ci_a["ci_high"] + 1e-12


# ---------------------------------------------------------------------------
# (c) Mini-MC-Smoke (end-to-end)
# ---------------------------------------------------------------------------

def test_measure_ratios_smoke() -> None:
    m = phy043.measure_ratios_wolff(radius=2, T=1.0, t_idx=0, n_measure=30,
                                    n_burn=20, n_seeds=2, master_seed=42)
    assert m.L == 5
    assert math.isfinite(m.u4_mean) and 0.9 <= m.u4_mean <= 2.1
    assert m.xi_ratio_mean is not None and m.xi_ratio_mean > 0.0
    assert m.u4_sem >= 0.0 and m.xi_ratio_sem >= 0.0


def test_run_phy043_smoke_structure(tmp_path) -> None:
    """End-to-end-Smoke: Report-Struktur + Gate-Dict vollstaendig; die
    Gates selbst duerfen bei Mini-Statistik FAIL sein (Integritaet der
    Auswertung, nicht Physik)."""
    report = phy043.run_phy043(radii=(2, 3), temps=(1.2, 1.5, 1.8),
                               n_measure=25, n_burn=15, n_seeds=2,
                               master_seed=42)
    assert report["module"].startswith("PHY043")
    assert report["lattices_L"] == [5, 7]
    assert set(report["pass_gates"]) == {
        "PASS_INPUT_COMPLETE", "PASS_U4_RANGE",
        "PASS_HIGH_T_ORDERING_XI", "PASS_HIGH_T_ORDERING_U4",
        "PASS_LOW_T_MERGE_XI", "PASS_SPLAY_XI_LARGEST_PAIR"}
    assert isinstance(report["overall_pass"], bool)
    assert "(5,7)" in report["pair_analysis"]["xi_ratio"]
    out = tmp_path / "phy043_smoke_report.txt"
    phy043.write_report(report, out)
    text = out.read_text(encoding="utf-8")
    assert "PASS-Gates" in text and "FINDING" in text
