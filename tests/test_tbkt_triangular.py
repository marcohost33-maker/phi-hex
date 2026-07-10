"""T_BKT(triangular)-Regressionstest mit per-Site-Helicity + Nelson-Kosterlitz.

Schneller Smoke-/Regressionstest (kleine L, wenige Sweeps) - die volle
Mess-Evidenz liegt in results/260604 PHY030 ... report.txt. Hier wird nur
abgesichert, dass die Pipeline mit per-Site-Normierung einen Nelson-
Kosterlitz-Crossing in einem physikalisch sinnvollen Band um T_BKT=1.418
liefert (kein 15%-Bias mehr wie bei der alten Flaechen-Normierung, die
T_BKT ~1.6 ergeben haette).
"""
from __future__ import annotations


from phy026_wolff_cluster import measure_helicity_wolff
from phi_hex_core_v2 import nelson_kosterlitz_line
from phy030_triangular_tbkt import per_L_crossing, REF_TRIANGULAR


# De-slow 2026-07-10 (Code-Audit 1.3): gemessen <3.5s - gehoert als
# fails-before-Gate in die schnelle CI-Suite, nicht hinter -m slow.
def test_nk_crossing_in_physical_band():
    """Per-L-Crossing (kleines L) liegt im Band [1.40, 1.55] um T_BKT=1.418.

    Per-Site-Konvention: Crossing leicht ueber 1.418 (Finite-Size von oben),
    aber klar UNTER dem ~1.6, das die alte Flaechen-Norm geliefert haette.
    """
    temps = (1.34, 1.42, 1.50, 1.58)
    rows = [
        measure_helicity_wolff(radius=4, T=T, n_measure=200, n_burn=150,
                               n_seeds=3, master_seed=42)
        for T in temps
    ]
    tc = per_L_crossing(rows)
    assert tc is not None, "kein NK-Crossing im Scan-Fenster"
    # Untergrenze schliesst Mess-Rauschen ab, Obergrenze schliesst den
    # alten Flaechen-Bias (~1.6) klar aus.
    assert 1.40 < tc < 1.56, f"Crossing T_x={tc:.4f} ausserhalb Band"
    # Sanity: Reference im erwarteten Bereich
    assert abs(REF_TRIANGULAR - 1.418) < 1e-9


def test_nk_line_value():
    """Nelson-Kosterlitz-Gerade ist 2T/pi (Defn-Sanity)."""
    import math
    assert abs(nelson_kosterlitz_line(1.418) - 2 * 1.418 / math.pi) < 1e-12
