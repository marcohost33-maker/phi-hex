"""Schnelle Korrektheits-Gates fuer PHY040 (Wang-Landau-entropischer Helicity).

Die teure WL-g(E)-Konvergenz + Wolff-Quervalidierung laufen als slow/lokal
(im __main__ des Moduls, Gate-Log in results/). Die CI-schnellen Gates pruefen
MC-frei die exakte Mathematik: kanonische Rueckgewichtung, die C-eliminierte
Paar-Bestimmung auf synthetischen Weber-Minnhagen-Kurven (Soll-T_BKT bekannt),
und die Upsilon-Verdrahtung gegen die getestete PHY039-Form.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

import conftest

phy040 = conftest._load(
    "phy040", "260616 PHY040 wang-landau entropic helicity v01.py")


def _make_res(centers, lng, micro_x, micro_y=None, n=64):
    centers = np.asarray(centers, float)
    return phy040.WLResult(
        L=8, n=n, centers=centers, lng=np.asarray(lng, float),
        mask=np.ones(len(centers), bool),
        micro_x={k: np.asarray(v, float) for k, v in micro_x.items()},
        micro_y={k: np.asarray(v, float) for k, v in (micro_y or micro_x).items()},
        wl_sweeps=0, prod_sweeps=0)


def test_canonical_weights_two_level():
    """g(E) exp(-E/T) Rueckgewichtung == analytisches 2-Niveau-Boltzmann."""
    centers = [-10.0, -4.0]
    lng = [math.log(1.0), math.log(3.0)]      # g0=1, g1=3 (Entartung)
    res = _make_res(centers, lng, {k: [0.0, 0.0] for k in phy040.OBS})
    for T in (0.5, 1.0, 2.0):
        w = phy040._canonical_weights(res, T)
        b = 1.0 / T
        g = np.array([1.0, 3.0])
        E = np.array(centers)
        ref = g * np.exp(-b * E)
        ref = ref / ref.sum()
        assert np.allclose(w, ref, atol=1e-12), (T, w, ref)


def test_y2_from_micro_wiring():
    """Bei auf einen Bin konzentriertem Gewicht reduziert sich Upsilon_2 auf
    (c - beta(s2 - s^2))/n - die direkte Helicity-Definition."""
    centers = [-8.0, -6.0]
    lng = [0.0, -50.0]            # Bin 0 dominiert komplett
    micro = {k: [0.0, 0.0] for k in phy040.OBS}
    micro["c"] = [1.3, 9.9]
    micro["s2"] = [0.7, 9.9]
    micro["s"] = [0.1, 9.9]
    res = _make_res(centers, lng, micro, n=64)
    T = 0.9
    w = phy040._canonical_weights(res, T)
    got = phy040._y2_from_micro(res.micro_x, w, res.n, 1.0 / T)
    exp = (1.3 - (1.0 / T) * (0.7 - 0.1 ** 2)) / 64
    assert abs(got - exp) < 1e-9, (got, exp)


def test_y4_from_micro_connected_form():
    """N*Upsilon_4 aus den (auf einen Bin konzentrierten) kanonischen Mitteln
    == die direkt ausgeschriebene verbundene 4.-Ordnungs-Form (1:1 PHY039)."""
    vals = {"c": 1.1, "s": 0.05, "s3": 0.02, "c4": 0.9, "s2": 0.6,
            "c2": 1.4, "cs2": 0.7, "s4": 0.5, "ss3": 0.3, "s3c": 0.01,
            "cs": 0.04}
    centers = [-8.0, -6.0]
    lng = [0.0, -80.0]                  # Bin 0 dominiert -> Mittel = vals
    micro = {k: [vals[k], 7.0] for k in phy040.OBS}
    res = _make_res(centers, lng, micro)
    beta = 1.0 / 0.9
    w = phy040._canonical_weights(res, 0.9)
    got = phy040._y4_from_micro(res.micro_x, w, beta)
    sa, ca, s3a, c4a = vals["s"], vals["c"], vals["s3"], vals["c4"]
    s2 = vals["s2"]
    kE3E1 = -vals["ss3"] + s3a * sa
    kE2E2 = vals["c2"] - ca * ca
    kE2E1E1 = vals["cs2"] - ca * s2 - 2 * sa * vals["cs"] + 2 * ca * sa * sa
    kE4_1 = (vals["s4"] - 4 * sa * vals["s3c"] - 3 * s2 * s2
             + 12 * s2 * sa * sa - 6 * sa ** 4)
    exp = (-c4a - beta * (4 * kE3E1 + 3 * kE2E2)
           + 6 * beta ** 2 * kE2E1E1 - beta ** 3 * kE4_1)
    assert abs(got - exp) < 1e-9, (got, exp)


def test_pair_tbkt_recovers_synthetic_wm():
    """C-eliminierte Paar-Bedingung auf synthetischen Weber-Minnhagen-Kurven
    mit eingebautem T_BKT (1/R(T,L) = 2 ln L + C + ln L*(T-T_true)): der
    Nulldurchgang muss T_true zurueckgeben."""
    T_true = 0.90
    C = -1.0
    grid = np.round(np.arange(0.80, 1.001, 0.01), 4)

    def y2(T, L):
        q = 2 * math.log(L) + C + math.log(L) * (T - T_true)
        R = 1.0 / q
        return (2 * T / math.pi) * (1.0 + R)
    for La, Lb in [(12, 24), (16, 24), (12, 16)]:
        ya = [y2(T, La) for T in grid]
        yb = [y2(T, Lb) for T in grid]
        tb = phy040.tbkt_pair_from_curves(grid, ya, La, yb, Lb)
        assert tb is not None and abs(tb - T_true) < 1e-3, (La, Lb, tb)


def test_pair_tbkt_general_ratio_not_just_doubling():
    """Die Paar-Bedingung nutzt -2 ln(Lb/La), nicht hart -2 ln 2: ein
    Nicht-Verdopplungs-Paar (16,24) muss ebenso T_true treffen."""
    T_true = 0.88
    grid = np.round(np.arange(0.80, 0.961, 0.01), 4)

    def y2(T, L):
        q = 2 * math.log(L) - 1.2 + 0.5 * math.log(L) * (T - T_true)
        return (2 * T / math.pi) * (1.0 + 1.0 / q)
    tb = phy040.tbkt_pair_from_curves(
        grid, [y2(T, 16) for T in grid], 16,
        [y2(T, 24) for T in grid], 24)
    assert tb is not None and abs(tb - T_true) < 1e-3, tb


@pytest.mark.slow
def test_wl_energy_matches_wolff_small():
    """Langsam: ein kurzer WL-Lauf (L=8) reproduziert die Wolff-Energie grob."""
    res = phy040.wang_landau_helicity(
        8, lnf_final=1e-4, prod_sweeps=4000, verbose=False)
    c = phy040.upsilon_curves(res, np.array([0.85, 0.95]))
    Ewf, _ = phy040.wolff_energy(8, 0.95, n_seeds=4, n_measure=200)
    k = 1
    assert abs(c["E"][k] / res.n - Ewf / res.n) < 0.05, (c["E"][k], Ewf)
