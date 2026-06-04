"""Tests fuer den Weber-Minnhagen-Log-Korrektur-Fit (PHY030 v02).

Zwei Ebenen:
(a) SYNTHETIK-ORAKEL (fails-before-faehig): kuenstliche Upsilon(T,L)-Daten
    werden EXAKT nach der WM-Form mit bekanntem T_true/C_true (plus
    deterministisch reproduzierbarem Rauschen) erzeugt. Der Fit MUSS T_true
    innerhalb Toleranz rekonstruieren. Das Orakel ist unabhaengig vom Monte-
    Carlo-Sampling - faellt der Fit-Algorithmus, faellt dieser Test.
(b) SMOKE: Mini-Parameter (kleine L, wenige Sweeps) sichern Pipeline-
    Integritaet end-to-end (Messung -> Fit -> endlicher Schaetzer im Band).
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from phy030_tbkt_wm_logfit import (
    wm_form,
    fit_tbkt_fixedT,
    tbkt_pair_C_eliminated,
    measure_helicity_wolff,
)


# ---------------------------------------------------------------------------
# (a) Synthetik-Orakel
# ---------------------------------------------------------------------------

def _make_synthetic(Ls, T_grid, T_true, C_true, noise_sigma, seed):
    """Erzeugt Upsilon(T,L) nach WM-Form bei T=T_true; fuer T != T_true wird
    eine glatte, physikalisch plausible Abweichung addiert (linear in (T-T_true)),
    sodass chi^2(T) bei T_true minimal wird. Plus reproduzierbares Gauss-Rauschen.

    Wichtig: die WM-Form gilt EXAKT nur bei T_BKT. Damit chi^2(T) ein klares
    Minimum bei T_true hat, modellieren wir fuer T != T_true eine zusaetzliche,
    nicht durch die 1-Parameter-WM-Form absorbierbare L-Abhaengigkeit
    (slope*(T-T_true)*ln L) - so wie reale Daten ausserhalb T_BKT von der
    reinen Log-Form abweichen.
    """
    rng = np.random.default_rng(seed)
    ups_of_T, sems_of_T = {}, {}
    slope = 0.8  # Staerke der Nicht-WM-Abweichung pro Temperatur-Distanz
    for T in T_grid:
        ups, sems = [], []
        for L in Ls:
            base = wm_form(T_true, L, C_true)
            # Abweichung von der reinen Log-Form fuer T != T_true:
            dev = slope * (T - T_true) * math.log(L)
            val = base + dev + rng.normal(0.0, noise_sigma)
            ups.append(val)
            sems.append(noise_sigma if noise_sigma > 0 else 1e-3)
        ups_of_T[T] = ups
        sems_of_T[T] = sems
    return ups_of_T, sems_of_T


def test_wm_form_definition():
    """WM-Form-Sanity: bei C->inf reduziert sich auf 2T/pi (Nelson-Kosterlitz);
    der Korrekturterm ist positiv und faellt mit L."""
    T = 1.4
    nk = 2.0 * T / math.pi
    assert abs(wm_form(T, 1000.0, 1e9) - nk) < 1e-6
    # Korrektur positiv:
    assert wm_form(T, 9.0, 2.0) > nk
    # Korrektur faellt mit L (groesseres L -> naeher an NK):
    assert wm_form(T, 19.0, 2.0) < wm_form(T, 9.0, 2.0)


@pytest.mark.parametrize("T_true,C_true", [(1.418, 2.0), (1.40, 0.5),
                                           (1.45, 5.0)])
def test_synthetic_oracle_recovers_T_true(T_true, C_true):
    """Fixed-T-Fit rekonstruiert T_true aus WM-konformen Synthetik-Daten.

    Das ist das fails-before-taugliche Orakel: ein fehlerhafter Fit (falscher
    argmin, falsche Parabel-Verfeinerung, falsche WM-Form) verfehlt T_true.
    """
    # Feines Gitter, das T_true einschliesst und es NICHT auf dem Rand hat.
    T_grid = [round(T_true - 0.02 + 0.005 * k, 3) for k in range(9)]
    Ls = [9, 13, 19, 29]
    ups_of_T, sems_of_T = _make_synthetic(
        Ls, T_grid, T_true, C_true, noise_sigma=1e-4, seed=12345)
    fit = fit_tbkt_fixedT(T_grid, Ls, ups_of_T, sems_of_T)
    # T_true innerhalb halber Gitterweite (0.005) bei sehr kleinem Rauschen.
    assert abs(fit["T_bkt"] - T_true) < 0.006, (
        f"WM-Fit T_bkt={fit['T_bkt']:.4f} verfehlt T_true={T_true}")
    # C ebenfalls plausibel rekonstruiert (grosse Toleranz: stark korreliert).
    assert abs(fit["C"] - C_true) < 2.0


def test_synthetic_oracle_with_noise_band():
    """Mit realistischem Rauschen bleibt T_true im 2-sigma-Band des Fits."""
    T_true, C_true = 1.418, 2.0
    T_grid = [round(1.40 + 0.01 * k, 3) for k in range(7)]  # 1.40..1.46
    Ls = [9, 13, 19]
    ups_of_T, sems_of_T = _make_synthetic(
        Ls, T_grid, T_true, C_true, noise_sigma=3e-3, seed=2026)
    fit = fit_tbkt_fixedT(T_grid, Ls, ups_of_T, sems_of_T)
    assert fit["sigma"] is not None, "kein inneres Minimum -> keine sigma"
    # T_true innerhalb 2-sigma (oder mind. innerhalb einer Gitterweite).
    band = max(2.0 * fit["sigma"], 0.011)
    assert abs(fit["T_bkt"] - T_true) < band, (
        f"T_bkt={fit['T_bkt']:.4f} +- {fit['sigma']:.4f} verfehlt "
        f"T_true={T_true}")


def test_pair_C_elimination_recovers_T_true():
    """C-freie Paar-Schaetzung findet T_true aus WM-Synthetik (Methode B).

    Bei reiner WM-Form ist die Paar-Bedingung exakt nur bei T_true erfuellt;
    fuer T != T_true bricht die L-unabhaengige C-Elimination, daher Nulldurchgang
    bei T_true. Quercheck unabhaengig vom Fixed-T-Least-Squares (A)."""
    T_true, C_true = 1.418, 2.0
    T_grid = [round(1.40 + 0.005 * k, 3) for k in range(13)]  # 1.40..1.46
    L1, L2 = 9, 19
    # Rauschfrei fuer das deterministische Orakel von B.
    ups1 = [wm_form(T_true, L1, C_true)
            + 0.8 * (T - T_true) * math.log(L1) for T in T_grid]
    ups2 = [wm_form(T_true, L2, C_true)
            + 0.8 * (T - T_true) * math.log(L2) for T in T_grid]
    sem = [1e-4] * len(T_grid)
    tp = tbkt_pair_C_eliminated(T_grid, ups1, sem, L1, ups2, sem, L2)
    assert tp is not None, "kein Nulldurchgang der Paar-Bedingung"
    assert abs(tp - T_true) < 0.006, f"Paar T_BKT={tp:.4f} verfehlt {T_true}"


def test_mutation_guard_sign_flip_kills_fit(monkeypatch):
    """Regression-Guard (Equalita-Auflage 2026-06-04): Sign-Flip im WM-
    Korrekturterm muss das chi^2 um Groessenordnungen verschlechtern.

    Verankert den manuellen Mutations-Beleg (chi^2-Ratio ~1e6 bei echten
    Daten) als automatisierten CI-Case: wer wm_form (Produktions-Funktion)
    still bricht, reisst diesen Test."""
    import phy030_tbkt_wm_logfit as mod

    T_true, C_true = 1.418, 2.0
    T_grid = [round(1.40 + 0.01 * k, 3) for k in range(5)]
    Ls = [9, 13, 19]
    ups_of_T, sems_of_T = _make_synthetic(
        Ls, T_grid, T_true, C_true, noise_sigma=1e-3, seed=777)
    fit_ok = mod.fit_tbkt_fixedT(T_grid, Ls, ups_of_T, sems_of_T)

    def wm_mutated(T, L, C):
        # Vorzeichen des Log-Korrekturterms gekippt (die Audit-Mutation).
        return (2.0 * T / math.pi) * (1.0 - 1.0 / (2.0 * math.log(L) + C))

    monkeypatch.setattr(mod, "wm_form", wm_mutated)
    fit_bad = mod.fit_tbkt_fixedT(T_grid, Ls, ups_of_T, sems_of_T)
    # Mutierter Fit muss um >=3 Groessenordnungen schlechter sein.
    assert fit_bad["chi2_min"] > 1e3 * max(fit_ok["chi2_min"], 1e-12), (
        f"Mutations-Guard greift nicht: chi2_ok={fit_ok['chi2_min']:.3e} "
        f"chi2_mut={fit_bad['chi2_min']:.3e}")


# ---------------------------------------------------------------------------
# (b) Smoke: Pipeline-Integritaet auf Mini-Parametern
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_smoke_pipeline_mini():
    """End-to-end Mini-Lauf: Messung -> WM-Fit liefert endlichen Schaetzer
    in einem physikalisch sinnvollen Band um T_BKT=1.418.

    Kleine L, wenige Sweeps - kein Praezisions-Anspruch, nur Integritaet:
    der Fit darf nicht crashen, muss endliches T_bkt und C liefern, und der
    Schaetzer muss grob im Uebergangsband liegen.
    """
    Ls = [9, 13]
    T_grid = [1.40, 1.42, 1.44]
    ups_of_T, sems_of_T = {}, {}
    for T in T_grid:
        ups_of_T[T], sems_of_T[T] = [], []
        for L in Ls:
            r = (L - 1) // 2
            res = measure_helicity_wolff(radius=r, T=T, n_measure=80,
                                         n_burn=60, n_seeds=2, master_seed=42)
            ups_of_T[T].append(res.upsilon_mean)
            sems_of_T[T].append(res.upsilon_sem)
    fit = fit_tbkt_fixedT(T_grid, Ls, ups_of_T, sems_of_T)
    assert math.isfinite(fit["T_bkt"])
    assert math.isfinite(fit["C"])
    # Grobes physikalisches Band (Mini-Statistik -> weit gefasst).
    assert 1.30 < fit["T_bkt"] < 1.55, (
        f"Smoke T_bkt={fit['T_bkt']:.4f} ausserhalb Band")
