"""PHY034 HKS-L->inf-Extrapolation: Korrektheits-Gates (analytisch, kein MC).

Vier Ebenen, alle instant + fails-before-faehig:

(a) REPRODUKTIONS-GUARD: die C-eliminierten Paar-Schaetzer aus dem embedded
    PHY032-Upsilon-Gitter MUESSEN die committed PHY032-Werte (0.6014/0.5917/
    0.5968) reproduzieren. Ein Tippfehler im Daten-Embed bricht das Gate.
(b) EXTRAPOLATIONS-ORAKEL: bei rein quadratischer Drift T(Lc)=t*+a/(ln Lc)^2
    MUSS die p=2-Extrapolation t* EXAKT rekonstruieren; bei zusaetzlicher
    kubischer Komponente reduziert sie den Bias ggue. dem Roh-Paar (HKS).
(c) OLS-INTERCEPT: bekannte Gerade -> exakter Achsenabschnitt.
(d) BEFUND-GATE: die Primaer-Extrapolation (p=2, geom) bewegt T_BKT vom
    Roh-Wert 0.5917 in das Referenz-Intervall [0.573, 0.576] (Konsistenz-
    Check der honeycomb-"Diskrepanz" als finite-size-Artefakt).
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from phy034_hks_extrapolation import (
    BARE_LARGEST_PAIR,
    DOUBLING_PAIRS,
    PHY032_PAIRS_EXPECTED,
    PHY032_SEM,
    PHY032_UPS,
    REF_HONEYCOMB_DEDICATED,
    REF_HONEYCOMB_MULTI,
    REF_MID,
    char_size,
    extrapolate_pairs,
    extrapolation_ci,
    ols_intercept,
    recompute_pairs,
    run_phy034,
    synthetic_pair_value,
)

_PHY032_REPORT = (
    Path(__file__).resolve().parent.parent / "results"
    / "260607 PHY032 honeycomb wm-logfit bootstrap report.txt"
)


# ---------------------------------------------------------------------------
# (a0) DRIFT-GUARD: embedded Upsilon-Gitter == committed PHY032-Report
# ---------------------------------------------------------------------------

def _parse_phy032_grid(text: str):
    """Parse 'L=12: T=0.5600:0.5368+-0.0035 ...'-Zeilen -> {L:[(ups,sem)]}."""
    out = {}
    for m in re.finditer(r"^\s*L=(\d+):(.*)$", text, re.MULTILINE):
        L = int(m.group(1))
        cells = re.findall(r"T=[\d.]+:([\d.]+)\+-([\d.]+)", m.group(2))
        out[L] = [(float(u), float(s)) for u, s in cells]
    return out


def test_embedded_grid_matches_committed_phy032_report():
    """Reality-Anchor: die im Skript eingebetteten PHY032-Messwerte MUESSEN
    Zeichen-fuer-Zeichen dem committed Report entsprechen (kein stiller Drift
    der Datenbasis)."""
    assert _PHY032_REPORT.exists(), _PHY032_REPORT
    grid = _parse_phy032_grid(_PHY032_REPORT.read_text(encoding="utf-8"))
    assert set(grid) == set(PHY032_UPS), (set(grid), set(PHY032_UPS))
    for L in PHY032_UPS:
        ups = [c[0] for c in grid[L]]
        sem = [c[1] for c in grid[L]]
        assert ups == PHY032_UPS[L], (L, ups, PHY032_UPS[L])
        assert sem == PHY032_SEM[L], (L, sem, PHY032_SEM[L])


# ---------------------------------------------------------------------------
# (a) Reproduktions-Guard
# ---------------------------------------------------------------------------

def test_pairs_reproduce_phy032():
    """C-eliminierte Paare aus dem embedded Gitter == PHY032-Report."""
    pairs = recompute_pairs()
    for k, expect in PHY032_PAIRS_EXPECTED.items():
        got = pairs[k]
        assert got is not None, k
        assert abs(got - expect) < 1.5e-3, (k, got, expect)


# ---------------------------------------------------------------------------
# (b) Extrapolations-Orakel
# ---------------------------------------------------------------------------

def _synthetic_pairs(t_star, a, b, size_mode="geom"):
    return {
        (l1, l2): synthetic_pair_value(
            t_star, a, b, char_size(l1, l2, size_mode))
        for (l1, l2) in DOUBLING_PAIRS
    }


def test_extrapolation_quadratic_drift_is_exact():
    """Rein quadratische Drift (b=0) -> p=2 rekonstruiert t* exakt."""
    t_star = 0.5745
    synth = _synthetic_pairs(t_star, a=0.20, b=0.0)
    got = extrapolate_pairs(synth, power=2, size_mode="geom")
    assert abs(got - t_star) < 1e-12, got


def test_extrapolation_reduces_subleading_bias():
    """Mit kubischer Restkorrektur ist die p=2-Extrapolation naeher an t*
    als das groesste Roh-Paar (HKS: sub-leading logs verschieben den Wert)."""
    t_star = 0.5745
    synth = _synthetic_pairs(t_star, a=0.20, b=0.05)
    extrap = extrapolate_pairs(synth, power=2, size_mode="geom")
    bare_largest = synth[DOUBLING_PAIRS[-1]]  # groesstes L, ohne Extrapolation
    assert abs(extrap - t_star) < abs(bare_largest - t_star)


# ---------------------------------------------------------------------------
# (c) OLS-Intercept
# ---------------------------------------------------------------------------

def test_ols_intercept_recovers_known_line():
    xs = [1.0, 2.0, 3.0, 4.0]
    ys = [2.0 + 3.0 * x for x in xs]  # Achsenabschnitt 2, Steigung 3
    intercept, slope = ols_intercept(xs, ys)
    assert abs(intercept - 2.0) < 1e-12
    assert abs(slope - 3.0) < 1e-12


def test_log_variable_monotone():
    """u=1/(ln L)^2 faellt mit L -> Extrapolationsrichtung wohldefiniert."""
    u_small = 1.0 / (math.log(12) ** 2)
    u_large = 1.0 / (math.log(48) ** 2)
    assert u_large < u_small


# ---------------------------------------------------------------------------
# (d) Befund-Gate (Physik-Konsistenz)
# ---------------------------------------------------------------------------

def test_primary_extrapolation_in_reference_band():
    """Primaer (HKS p=2, geom) liegt im/am Referenz-Intervall [0.573,0.576]."""
    rep = run_phy034()
    pr = rep["primary"]
    assert REF_HONEYCOMB_MULTI - 0.01 <= pr <= REF_HONEYCOMB_DEDICATED + 0.01


def test_extrapolation_moves_toward_reference():
    """Extrapolation reduziert den Offset zur Referenz ggue. dem Roh-Paar."""
    rep = run_phy034()
    assert abs(rep["primary"] - REF_MID) < abs(BARE_LARGEST_PAIR - REF_MID)


def test_systematic_band_brackets_reference():
    """Das Systematik-Band ueberlappt das Literatur-Intervall."""
    rep = run_phy034()
    lo, hi = rep["band"]
    assert lo <= REF_HONEYCOMB_DEDICATED and hi >= REF_HONEYCOMB_MULTI


def test_overall_gate_passes():
    assert run_phy034()["overall"] is True


# ---------------------------------------------------------------------------
# (e) Statistische CI-Propagation
# ---------------------------------------------------------------------------

def test_statistical_ci_is_deterministic_and_brackets_primary():
    """Fixer Seed -> reproduzierbar; CI umschliesst den Punktschaetzer und
    deckt das Referenz-Intervall (2-Punkt-Extrapolation ist breit)."""
    pairs = recompute_pairs()
    med1, lo1, hi1, std1 = extrapolation_ci(pairs, seed=42)
    med2, lo2, hi2, std2 = extrapolation_ci(pairs, seed=42)
    assert (med1, lo1, hi1, std1) == (med2, lo2, hi2, std2)  # deterministisch
    assert lo1 < med1 < hi1
    assert lo1 <= REF_HONEYCOMB_DEDICATED and hi1 >= REF_HONEYCOMB_MULTI


def test_extrapolation_amplifies_pair_error():
    """Der 2-Punkt-Hebel verstaerkt die Paar-Streuung: die CI-Breite ist
    GROESSER als die groesste per-Paar Bootstrap-Streuung (Ehrlichkeit ueber
    die Schwaeche der 2-Punkt-Extrapolation)."""
    pairs = recompute_pairs()
    _, lo, hi, std = extrapolation_ci(pairs, seed=7)
    assert std > 0.0034328  # > max per-Paar boot_std (24,48)
