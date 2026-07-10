"""PHY030 v02: T_BKT(triangular) via Weber-Minnhagen-Log-Korrektur-Fit

Coworker Research / Coworkerz, 04. Juni 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

KONTEXT (SOTA-Verbesserung gegenueber v01):
v01 (`260604 PHY030 triangular tbkt per-site v01.py`) bestimmt T_BKT via
per-L-Crossing Upsilon(T,L)=2T/pi + 1/ln(L)-Extrapolation und liefert
konservative SCHRANKEN: groesstes-L-Crossing (L=19) = 1.456 (+2.69%) als
obere, 1/ln(L)-Extrapolation = 1.382 (-2.56%) als untere Schranke; die
Referenz 1.418 liegt zwischen beiden. v01 bleibt gueltig als konservative
Schranken-Methode (Gate gegen groesstes L).

v02 liefert einen PRAEZISEN PUNKTSCHAETZER mittels der Weber-Minnhagen-
Finite-Size-Form fuer den Helicity-Modulus exakt am BKT-Punkt:

    Upsilon(T_BKT, L) = (2 T_BKT / pi) * (1 + 1 / (2 ln L + C))

mit einem einzigen unbekannten Konstanten-Parameter C (Weber & Minnhagen,
Phys. Rev. B 37, 5986(R) (1988); single-parameter log correction).

METHODIK (Standard, arXiv:1302.2900 Hsieh/Kao/Sandvik; arXiv:2406.12076):
Zwei komplementaere Auswertungen derselben Wolff-Daten:

  (A) FIXED-T-LEAST-SQUARES (Punktschaetzer):
      Fuer jede Kandidaten-Temperatur T wird die gemessene Kurve
      Upsilon(T,L) ueber alle L gegen die WM-Form gefittet, mit C als
      EINZIGEM freien Parameter (1-Parameter-Least-Squares, gewichtet mit
      1/sem^2). Bei der WAHREN T_BKT beschreibt die WM-Form die L-Abhaengigkeit
      exakt -> minimales chi^2. Daher:  T_BKT = argmin_T chi^2(T).
      Der Vorfaktor 2T/pi ist bei festem T bekannt; der Fit ist linear in
      der Variablen u = 1/(2 ln L + C) bei festem C, also exakt loesbar je C
      und ueber C eindimensional minimiert (robust, keine Startwert-Sorge).
      Unsicherheit: chi^2(T)-Parabel, Band chi^2_min + 1 (1-sigma).

  (B) PAAR-C-ELIMINATION (Querschecks, C-frei, arXiv:1302.2900):
      Fuer ein Paar (L1,L2) eliminiert man C aus zwei WM-Gleichungen und
      definiert eine groessen-abhaengige T_BKT(L1,L2) ueber die Bedingung,
      dass beide gemessenen Upsilon mit EINEM gemeinsamen C konsistent sind.
      Konkret: R(T,L) := pi*Upsilon(T,L)/(2T) - 1 = 1/(2 ln L + C), also
      1/R(T,L1) - 1/R(T,L2) = 2 ln(L2/L1) (von C unabhaengig). Die Temperatur
      mit Nulldurchgang dieser Differenz ist T_BKT(L1,L2). Mittel ueber alle
      Paare als robuster Quercheck.

Beide Methoden nutzen DIESELBE Wolff-Pipeline und DIESELBEN Parameter wie v01
(radii 4/6/9, n_measure=400, n_burn=300, n_seeds=4, master_seed=42), nur auf
einem feineren T-Gitter um den Uebergang (Delta T = 0.01 statt 0.04-0.06).

Referenz (web-validiert, arXiv:2501.07388): T_BKT(triangular) = 1.418 J/k_B.

EVIDENZ: deterministischer Gate-Report nach results/.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parent

REF_TRIANGULAR = 1.418  # arXiv:2501.07388
# v01-Schranken (konservativ; aus dem v01-Gate-Report, gleiche Pipeline):
V01_UPPER_BOUND = 1.456  # groesstes-L-Crossing (L=19)
V01_LOWER_BOUND = 1.382  # 1/ln(L)-Extrapolation


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


core = _load("phi_hex_core_v2", "260602 PHI HEX core v2 2 hardened.py")
wolff = _load("phy026_wolff_cluster", "260602 PHY026 wolff cluster v01.py")

measure_helicity_wolff = wolff.measure_helicity_wolff


# ============================================================================
# Weber-Minnhagen-Form + Fit-Kern (pure, testbar ohne Monte-Carlo)
# ============================================================================

def wm_form(T_bkt: float, L: float, C: float) -> float:
    """Weber-Minnhagen-Finite-Size-Form des Helicity-Modulus am BKT-Punkt.

    Upsilon(T_bkt, L) = (2 T_bkt / pi) * (1 + 1 / (2 ln L + C)).
    """
    return (2.0 * T_bkt / math.pi) * (1.0 + 1.0 / (2.0 * math.log(L) + C))


def _best_C_and_chi2(T: float, Ls, ups, sems):
    """Fixed-T-Fit: minimiere ueber C das gewichtete Residuum der WM-Form.

    Bei festem T und festem C ist die WM-Form keine freie lineare Regression
    mehr (kein freier Vorfaktor) - der Vorfaktor 2T/pi ist physikalisch
    festgelegt. Wir minimieren daher chi^2(C) direkt eindimensional. chi^2 ist
    glatt und (im physikalisch relevanten Bereich C > -2 ln L_min) unimodal;
    ein robustes Golden-Section-Minimum genuegt. Rueckgabe: (C_opt, chi2_opt).
    """
    lnL = [math.log(L) for L in Ls]
    w = [1.0 / max(s, 1e-9) ** 2 for s in sems]

    def chi2(C: float) -> float:
        # Definitionsluecke vermeiden: 2 ln L + C muss > 0 fuer alle L bleiben.
        denom_min = 2.0 * min(lnL) + C
        if denom_min <= 1e-6:
            return 1e18
        s = 0.0
        for k, L in enumerate(Ls):
            # Single source of truth: IMMER die oeffentliche wm_form
            # (keine Inline-Duplikation -- Mutations-Guard im Test haengt
            # genau hier; Equalita-Auflage 2026-06-04).
            model = wm_form(T, L, C)
            s += w[k] * (ups[k] - model) ** 2
        return s

    # Suchintervall fuer C: physikalisch typ. C in [-2, 20]; untere Schranke
    # strikt oberhalb der Definitionsluecke -2 ln(L_min).
    lo = -2.0 * min(lnL) + 0.05
    hi = 50.0
    # Golden-Section-Search.
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c1 = b - gr * (b - a)
    c2 = a + gr * (b - a)
    f1, f2 = chi2(c1), chi2(c2)
    for _ in range(200):
        if f1 < f2:
            b, c2, f2 = c2, c1, f1
            c1 = b - gr * (b - a)
            f1 = chi2(c1)
        else:
            a, c1, f1 = c1, c2, f2
            c2 = a + gr * (b - a)
            f2 = chi2(c2)
        if abs(b - a) < 1e-8:
            break
    C_opt = 0.5 * (a + b)
    return C_opt, chi2(C_opt)


def fit_tbkt_fixedT(T_grid, Ls, ups_of_T, sems_of_T):
    """Methode (A): T_BKT = argmin_T chi^2(T) des fixed-T WM-Fits.

    ups_of_T[T] / sems_of_T[T] sind Listen ueber Ls. Rueckgabe-Dict mit
    pro-T (C, chi2), T_BKT (parabolisch verfeinertes Minimum), C am Minimum,
    Residuen je L bei T_BKT, und 1-sigma-Unsicherheit aus chi^2_min + 1.
    """
    per_T = []
    for T in T_grid:
        C, c2 = _best_C_and_chi2(T, Ls, ups_of_T[T], sems_of_T[T])
        per_T.append((T, C, c2))

    chi2_arr = np.array([row[2] for row in per_T])
    i_min = int(np.argmin(chi2_arr))
    T_min = per_T[i_min][0]
    C_min = per_T[i_min][1]
    chi2_min = per_T[i_min][2]

    # Parabolische Verfeinerung des Minimums (3-Punkt) + 1-sigma-Band.
    T_bkt = T_min
    sigma = None
    if 0 < i_min < len(per_T) - 1:
        x0, x1 = per_T[i_min - 1][0], per_T[i_min][0]
        y0, y1, y2 = (chi2_arr[i_min - 1], chi2_arr[i_min],
                      chi2_arr[i_min + 1])
        # Parabel a*(x-x1)^2 + ... ueber 3 Punkte; Scheitel-Verschiebung.
        denom = (y0 - 2.0 * y1 + y2)
        if abs(denom) > 1e-12:
            dx = x1 - x0  # gleichmaessiges Gitter angenommen
            T_bkt = float(x1 - 0.5 * dx * (y2 - y0) / denom)
            a = denom / (2.0 * dx * dx)  # Kruemmung
            if a > 0:
                # chi^2(T) ~ a*(T - T_bkt)^2 + chi2_min; chi^2_min+1 -> 1-sigma
                sigma = float(math.sqrt(1.0 / a))

    # Residuen je L am Gitter-argmin T_min: Daten UND Modell bei DERSELBEN
    # Temperatur (T_min, mit dem dort optimierten C_min) ausgewertet.
    # P2-Fix 2026-07-10 (Code-Audit L4): frueher wurden Daten bei T_min mit
    # dem Modell beim parabolisch verfeinerten T_bkt verglichen - das mischte
    # zwei Temperaturen und addierte einen systematischen Offset bis
    # ~(2/pi)*Gitterweite/2 in jedes Residuum.
    residuals = {}
    for k, L in enumerate(Ls):
        model = wm_form(T_min, L, C_min)
        meas = ups_of_T[T_min][k]
        residuals[L] = {"measured": meas, "wm_model": model,
                        "residual": meas - model}

    # Haertung 2026-07-10 (Code-Audit L5): C am Rand des Golden-Section-
    # Suchintervalls [-2 ln L_min + 0.05, 50] ist ein Diagnose-Signal
    # (Optimum ausserhalb des Fensters -> chi^2(T)-Profil dort verzerrt).
    c_lo = -2.0 * min(math.log(L) for L in Ls) + 0.05
    c_at_bound = (C_min - c_lo) < 1e-3 or (50.0 - C_min) < 1e-3

    return {
        "per_T": [{"T": t, "C": c, "chi2": x} for (t, c, x) in per_T],
        "T_bkt": T_bkt,
        "C": C_min,
        "C_at_bound": c_at_bound,
        "chi2_min": chi2_min,
        "sigma": sigma,
        "residuals_per_L": residuals,
        "T_grid_argmin": T_min,
    }


def tbkt_pair_C_eliminated(T_grid, ups_L1, sems_L1, L1, ups_L2, sems_L2, L2):
    """Methode (B): C-freie Paar-Schaetzung T_BKT(L1,L2).

    1/R(T,L) = 2 ln L + C mit R = pi*Upsilon/(2T) - 1; also
    1/R(T,L1) - 1/R(T,L2) = 2 ln(L1) - 2 ln(L2) = -2 ln(L2/L1).
    Nulldurchgang von g(T) = [1/R(T,L1) - 1/R(T,L2)] - (-2 ln(L2/L1)).
    """
    target = -2.0 * math.log(L2 / L1)
    g = []
    for i, T in enumerate(T_grid):
        R1 = math.pi * ups_L1[i] / (2.0 * T) - 1.0
        R2 = math.pi * ups_L2[i] / (2.0 * T) - 1.0
        # P1-Fix 2026-07-10 (Code-Audit M1): auf dem physikalischen WM-Ast
        # ist 1/R = 2 ln L + C > 0, also R > 0. Punkte mit R <= 0 (jenseits
        # des finite-L-NK-Crossings) werden uebersprungen; sonst liefert der
        # Pol von 1/R einen Schein-Nulldurchgang, der in Bootstrap-Resamples
        # als "valider" T_BKT-Wert zaehlte (Pol-Artefakt statt WM-Wurzel).
        if R1 < 1e-6 or R2 < 1e-6:
            g.append((T, None))
        else:
            g.append((T, (1.0 / R1 - 1.0 / R2) - target))
    # P2-Fix 2026-07-10b (Codex-Review PR#23): Nulldurchgang NUR zwischen
    # BENACHBARTEN Gitterpunkten akzeptieren. Das fruehere Kompaktieren der
    # None-Punkte interpolierte ueber Pol-Luecken hinweg - ein Vorzeichen-
    # wechsel QUER ueber eine R<=0-Luecke ist kein physikalischer WM-
    # Nulldurchgang (disconnected branch), sondern das Pol-Artefakt in
    # neuer Form (relevant fuer verrauschte/Bootstrap-Eingaben).
    for i in range(len(g) - 1):
        (t0, d0), (t1, d1) = g[i], g[i + 1]
        if d0 is None or d1 is None:
            continue
        if d0 == 0.0:
            return t0
        if d0 * d1 < 0.0:
            return float(t0 + (t1 - t0) * (-d0) / (d1 - d0))
    return None


# ============================================================================
# Hauptlauf
# ============================================================================

@dataclass
class Measurement:
    T: float
    L: int
    upsilon_mean: float
    upsilon_sem: float


def run_phy030_v02(radii=(4, 6, 9),
                   temps=(1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43,
                          1.44, 1.45, 1.46),
                   n_measure=400, n_burn=300, n_seeds=4, master_seed=42):
    t_start = time.time()
    print("PHY030 v02: T_BKT(triangular) via Weber-Minnhagen-Log-Fit")
    print("Coworker Research / Coworkerz, 04. Juni 2026")
    print("=" * 70)
    print(f"Referenz T_BKT(triangular) = {REF_TRIANGULAR} J/k_B "
          f"(arXiv:2501.07388)")
    print("WM-Form: Upsilon(T_BKT,L) = (2 T_BKT/pi)*(1 + 1/(2 ln L + C))")
    print("Methodik: arXiv:1302.2900 / Weber-Minnhagen PRB 37, 5986(R) "
          "(1988)\n")

    Ls = [2 * r + 1 for r in radii]
    print(f"Gitter L = {Ls}  T-Gitter (Delta=0.01) = {list(temps)}")
    print(f"Wolff: n_measure={n_measure}, n_burn={n_burn}, "
          f"n_seeds={n_seeds}, master_seed={master_seed}\n")

    # --- Messung (gleiche Pipeline/Params wie v01) ---
    meas = {L: {} for L in Ls}
    print("Per-Site-Helicity Upsilon(T,L):")
    for r in radii:
        L = 2 * r + 1
        cells = []
        for T in temps:
            res = measure_helicity_wolff(radius=r, T=T, n_measure=n_measure,
                                         n_burn=n_burn, n_seeds=n_seeds,
                                         master_seed=master_seed)
            meas[L][T] = Measurement(T=T, L=L,
                                     upsilon_mean=res.upsilon_mean,
                                     upsilon_sem=res.upsilon_sem)
            cells.append(f"T={T:.2f}:{res.upsilon_mean:.4f}"
                         f"+-{res.upsilon_sem:.4f}")
        print(f"  L={L:2d}: " + "  ".join(cells))

    # Datenstrukturen fuer Fits.
    ups_of_T = {T: [meas[L][T].upsilon_mean for L in Ls] for T in temps}
    sems_of_T = {T: [meas[L][T].upsilon_sem for L in Ls] for T in temps}

    # --- (A) Fixed-T Least-Squares-WM-Fit ---
    print("\n(A) Fixed-T WM-Fit: chi^2(T) (1-Parameter C, gewichtet 1/sem^2):")
    fitA = fit_tbkt_fixedT(list(temps), Ls, ups_of_T, sems_of_T)
    for row in fitA["per_T"]:
        marker = "  <- argmin" if row["T"] == fitA["T_grid_argmin"] else ""
        print(f"  T={row['T']:.2f}: chi^2={row['chi2']:.4e}  "
              f"C={row['C']:+.3f}{marker}")
    tA = fitA["T_bkt"]
    sigA = fitA["sigma"]
    sig_str = f"{sigA:.4f}" if sigA is not None else "n/a (Rand)"
    print(f"\n  T_BKT(WM-Fit) = {tA:.4f} +- {sig_str}  "
          f"(C={fitA['C']:+.3f}"
          f"{', C AM SUCHRAND - Diagnose pruefen' if fitA['C_at_bound'] else ''})")
    print(f"  Abw. vs Referenz {REF_TRIANGULAR}: "
          f"{(tA - REF_TRIANGULAR) / REF_TRIANGULAR * 100:+.2f}%")
    print("  Residuen je L bei T_BKT (gemessen - WM-Modell):")
    for L in Ls:
        rr = fitA["residuals_per_L"][L]
        print(f"    L={L:2d}: gemessen={rr['measured']:.4f}  "
              f"WM={rr['wm_model']:.4f}  Residuum={rr['residual']:+.4f}")

    # --- (B) Paar-C-Eliminations-Querscheck ---
    print("\n(B) Paar-C-Eliminations-Querscheck T_BKT(L1,L2) "
          "(C-frei, arXiv:1302.2900):")
    pair_estimates = {}
    for a in range(len(Ls)):
        for b in range(a + 1, len(Ls)):
            L1, L2 = Ls[a], Ls[b]
            ups1 = [meas[L1][T].upsilon_mean for T in temps]
            ups2 = [meas[L2][T].upsilon_mean for T in temps]
            sem1 = [meas[L1][T].upsilon_sem for T in temps]
            sem2 = [meas[L2][T].upsilon_sem for T in temps]
            tp = tbkt_pair_C_eliminated(list(temps), ups1, sem1, L1,
                                        ups2, sem2, L2)
            pair_estimates[f"{L1},{L2}"] = tp
            tp_str = f"{tp:.4f}" if tp is not None else "kein Nulldurchgang"
            print(f"  T_BKT(L={L1:2d},L={L2:2d}) = {tp_str}")
    valid_pairs = [v for v in pair_estimates.values() if v is not None]
    pair_mean = float(np.mean(valid_pairs)) if valid_pairs else None
    if pair_mean is not None:
        print(f"  Mittel der Paar-Schaetzer = {pair_mean:.4f}  "
              f"(Abw. {(pair_mean - REF_TRIANGULAR) / REF_TRIANGULAR * 100:+.2f}%)")

    # --- Vergleich mit v01-Schranken ---
    print("\nVERGLEICH:")
    print(f"  v01 untere Schranke (1/ln L-Extrap.) = {V01_LOWER_BOUND}  "
          f"({(V01_LOWER_BOUND - REF_TRIANGULAR) / REF_TRIANGULAR * 100:+.2f}%)")
    print(f"  v01 obere Schranke  (L=19-Crossing)  = {V01_UPPER_BOUND}  "
          f"({(V01_UPPER_BOUND - REF_TRIANGULAR) / REF_TRIANGULAR * 100:+.2f}%)")
    print(f"  v02 WM-Punktschaetzer (A)            = {tA:.4f}  "
          f"({(tA - REF_TRIANGULAR) / REF_TRIANGULAR * 100:+.2f}%)")
    if pair_mean is not None:
        print(f"  v02 Paar-Mittel (B)                  = {pair_mean:.4f}  "
              f"({(pair_mean - REF_TRIANGULAR) / REF_TRIANGULAR * 100:+.2f}%)")
    print(f"  Referenz                             = {REF_TRIANGULAR}")

    # WM naeher an Referenz als beide v01-Schranken?
    rel_A = abs(tA - REF_TRIANGULAR) / REF_TRIANGULAR
    rel_v01_upper = abs(V01_UPPER_BOUND - REF_TRIANGULAR) / REF_TRIANGULAR
    rel_v01_lower = abs(V01_LOWER_BOUND - REF_TRIANGULAR) / REF_TRIANGULAR
    wm_closer = rel_A < min(rel_v01_upper, rel_v01_lower)

    # --- Gates ---
    pass_gates = {
        "PASS_WM_FIT_WITHIN_3PCT": rel_A < 0.03,
        "PASS_WM_FIT_WITHIN_REF_BAND": V01_LOWER_BOUND <= tA <= V01_UPPER_BOUND,
        "PASS_WM_CLOSER_THAN_V01_BOUNDS": bool(wm_closer),
        "PASS_PAIR_CROSSCHECK_EXISTS": pair_mean is not None,
        "PASS_RESIDUALS_SMALL": all(
            abs(fitA["residuals_per_L"][L]["residual"]) < 0.05 for L in Ls),
    }
    print("\nPASS-Gates:")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values())
    dt = time.time() - t_start
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}   (Laufzeit {dt:.0f}s)")

    print("\nINTERPRETATION:")
    if wm_closer and rel_A < 0.03:
        print("Der Weber-Minnhagen-Log-Fit liefert einen Punktschaetzer, der")
        print("naeher an der Referenz 1.418 liegt als beide konservativen")
        print("v01-Schranken - die explizite Log-Korrektur korrigiert den")
        print("Finite-Size-Bias der kleinen Gitter (L=9..19). v01 bleibt als")
        print("konservative Schranken-Methode gueltig.")
    else:
        print("Der WM-Fit liegt NICHT naeher an 1.418 als beide v01-Schranken")
        print("oder ausserhalb 3%. Bei nur drei kleinen Gittern (L=9..19) ist")
        print("die Log-Korrektur eine reale Grenze (sub-leading Korrekturen,")
        print("arXiv:1302.2900). Ehrlich ausgewiesen, nicht ueberinterpretiert.")

    return {
        "module": "PHY030_triangular_tbkt_per_site_v02_wm_logfit",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-06-04",
        "method": "Weber-Minnhagen single-parameter log correction; "
                  "fixed-T weighted least-squares (A) + pair C-elimination (B)",
        "wm_form": "Upsilon(T_BKT,L) = (2 T_BKT/pi)*(1 + 1/(2 ln L + C))",
        "reference_T_bkt_triangular": REF_TRIANGULAR,
        "reference_source": "arXiv:2501.07388",
        "method_source": "arXiv:1302.2900; Weber & Minnhagen PRB 37, 5986(R) (1988)",
        "lattices_L": Ls,
        "temperatures": list(temps),
        "wolff": {"n_measure": n_measure, "n_burn": n_burn,
                  "n_seeds": n_seeds, "master_seed": master_seed},
        "measurements": {
            str(L): {f"{T:.2f}": asdict(meas[L][T]) for T in temps}
            for L in Ls},
        "fit_A_fixedT": fitA,
        "pair_estimates_B": pair_estimates,
        "pair_mean_B": pair_mean,
        "v01_bounds": {"lower": V01_LOWER_BOUND, "upper": V01_UPPER_BOUND},
        "wm_closer_than_v01_bounds": bool(wm_closer),
        "pass_gates": pass_gates,
        "overall_pass": overall,
        "runtime_s": dt,
    }


def _clean(o):
    if o is None:
        return None
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, np.ndarray):
        return [_clean(x) for x in o.tolist()]
    if hasattr(o, "__dataclass_fields__"):
        return _clean(asdict(o))
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    return o


def write_report(report: dict, path: Path) -> None:
    """Deterministischer Text-Report (keine Laufzeit/Timestamps im Body)."""
    fitA = report["fit_A_fixedT"]
    Ls = report["lattices_L"]
    tA = fitA["T_bkt"]
    sigA = fitA["sigma"]
    sig_str = f"{sigA:.4f}" if sigA is not None else "n/a (Randminimum)"
    ref = report["reference_T_bkt_triangular"]
    lines = []
    lines.append("PHY030 v02 - T_BKT(triangular) via Weber-Minnhagen-Log-Fit")
    lines.append("Coworker Research / Coworkerz, 2026-06-04")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"WM-Form: {report['wm_form']}")
    lines.append(f"Methode: {report['method']}")
    lines.append(f"Quellen: {report['method_source']}")
    lines.append(f"Referenz T_BKT = {ref} J/k_B ({report['reference_source']})")
    lines.append("")
    lines.append(f"Gitter L = {Ls}")
    lines.append(f"T-Gitter = {report['temperatures']}")
    w = report["wolff"]
    lines.append(f"Wolff: n_measure={w['n_measure']}, n_burn={w['n_burn']}, "
                 f"n_seeds={w['n_seeds']}, master_seed={w['master_seed']}")
    lines.append("")
    lines.append("--- Messung Upsilon(T,L) ---")
    for L in Ls:
        cells = []
        for T in report["temperatures"]:
            m = report["measurements"][str(L)][f"{T:.2f}"]
            cells.append(f"T={T:.2f}:{m['upsilon_mean']:.4f}"
                         f"+-{m['upsilon_sem']:.4f}")
        lines.append(f"  L={L:2d}: " + "  ".join(cells))
    lines.append("")
    lines.append("--- (A) Fixed-T WM-Fit chi^2(T) ---")
    for row in fitA["per_T"]:
        mark = "  <- argmin" if row["T"] == fitA["T_grid_argmin"] else ""
        lines.append(f"  T={row['T']:.2f}: chi^2={row['chi2']:.4e}  "
                     f"C={row['C']:+.4f}{mark}")
    lines.append("")
    lines.append(f"  T_BKT(WM-Fit) = {tA:.4f} +- {sig_str}")
    lines.append(f"  C = {fitA['C']:+.4f}")
    lines.append(f"  Abweichung vs Referenz {ref}: "
                 f"{(tA - ref) / ref * 100:+.2f}%")
    lines.append("")
    lines.append("  Residuen je L bei T_BKT (gemessen - WM-Modell):")
    for L in Ls:
        rr = fitA["residuals_per_L"][L]
        lines.append(f"    L={L:2d}: gemessen={rr['measured']:.4f}  "
                     f"WM={rr['wm_model']:.4f}  Residuum={rr['residual']:+.4f}")
    lines.append("")
    lines.append("--- (B) Paar-C-Eliminations-Querscheck (C-frei) ---")
    for k, v in report["pair_estimates_B"].items():
        vs = f"{v:.4f}" if v is not None else "kein Nulldurchgang"
        lines.append(f"  T_BKT(L={k}) = {vs}")
    if report["pair_mean_B"] is not None:
        pm = report["pair_mean_B"]
        lines.append(f"  Mittel = {pm:.4f}  "
                     f"(Abw. {(pm - ref) / ref * 100:+.2f}%)")
    lines.append("")
    lines.append("--- Vergleich v01-Schranken / Referenz ---")
    lo = report["v01_bounds"]["lower"]
    hi = report["v01_bounds"]["upper"]
    lines.append(f"  v01 untere Schranke (1/ln L-Extrap.) = {lo}  "
                 f"({(lo - ref) / ref * 100:+.2f}%)")
    lines.append(f"  v01 obere Schranke  (L=19-Crossing)  = {hi}  "
                 f"({(hi - ref) / ref * 100:+.2f}%)")
    lines.append(f"  v02 WM-Punktschaetzer (A)            = {tA:.4f}  "
                 f"({(tA - ref) / ref * 100:+.2f}%)")
    if report["pair_mean_B"] is not None:
        pm = report["pair_mean_B"]
        lines.append(f"  v02 Paar-Mittel (B)                  = {pm:.4f}  "
                     f"({(pm - ref) / ref * 100:+.2f}%)")
    lines.append(f"  Referenz                             = {ref}")
    lines.append(f"  WM-Fit naeher an Referenz als beide v01-Schranken: "
                 f"{report['wm_closer_than_v01_bounds']}")
    lines.append("")
    lines.append("--- PASS-Gates ---")
    for k, v in report["pass_gates"].items():
        lines.append(f"  [{'PASS' if v else 'FAIL'}] {k}")
    lines.append(f"  OVERALL: {'PASS' if report['overall_pass'] else 'FAIL'}")
    lines.append("")
    lines.append("--- Caveats ---")
    lines.append("  - Nur drei kleine Gitter (L=9,13,19); sub-leading Log-")
    lines.append("    Korrekturen (arXiv:1302.2900) sind bei kleinem L eine")
    lines.append("    reale Grenze der Praezision.")
    lines.append("  - v01 bleibt gueltig als konservative Schranken-Methode.")
    lines.append("  - Statistik: n_seeds=4 (bewusst kleines Budget, Minuten).")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    report = run_phy030_v02()
    out = (_SRC.parent / "results" /
           "260604 PHY030 v02 wm-logfit report.txt")
    out.parent.mkdir(exist_ok=True)
    write_report(report, out)
    print(f"\nReport geschrieben: {out}")
    print("\n--- JSON-Report ---")
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
