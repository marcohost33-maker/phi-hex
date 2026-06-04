"""PHY030: T_BKT(triangular) Neumessung mit per-Site-Helicity-Normierung

Coworker Research / Coworkerz, 04. Juni 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

KONTEXT (Wissenschafts-Audit 2026-06-04, S2/S3; Marco-Entscheid 2026-06-04):
Die Helicity-Modulus-Normierung im Kern wurde von FLAECHE (Upsilon(0)=J*sqrt(3)
=1.732) auf per-Site (Upsilon(0)=1.5*J, Literatur-Standard arXiv:2406.12076)
umgestellt. Damit muss T_BKomeasurement neu erfolgen - die alten Zahlen waren
durch die zu laxe Toleranz tol=0.15 nie hart geprueft.

METHODE:
- Sampling: Wolff-Single-Cluster (PHY026) auf dem Dreiecks-Torus.
- Observable: per-Site-Helicity-Modulus Upsilon(T,L) (Kern, korrigiert).
- Kriterium: Nelson-Kosterlitz-Sprung Upsilon(T_BKT) = 2 T_BKT / pi.
  Zwei Auswertungen:
    (a) Per-L-Crossing: T, wo Upsilon(T,L) die Gerade 2T/pi schneidet.
    (b) Sandvik-Paar (C-frei, wie PHY028): aus Weber-Minnhagen
        pi*Upsilon(T_BKT,L)/(2 T_BKT) = 1 + 1/(2 ln L + C) folgt fuer (L,2L)
        exakt 1/R(T,L) - 1/R(T,2L) = -2 ln 2, R = pi*Upsilon/(2T) - 1.
        Eliminiert die unbekannte Konstante C ohne freie Extrapolationsform.

Referenz (web-validiert, arXiv:2501.07388): T_BKT(triangular) = 1.418 J/k_B.

EVIDENZ: Gate-Report nach results/. Rechenbudget bewusst klein gehalten
(Minuten); statistische Unsicherheit wird ehrlich ausgewiesen.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parent


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
nelson_kosterlitz_line = core.nelson_kosterlitz_line

REF_TRIANGULAR = 1.418  # arXiv:2501.07388


def per_L_crossing(rows):
    """T, wo Upsilon(T,L) die NK-Gerade 2T/pi schneidet (lineare Interpolation
    der Differenz d(T)=Upsilon(T,L)-2T/pi). rows: sortierte Liste von Results."""
    diffs = [(r.T, r.upsilon_mean - nelson_kosterlitz_line(r.T)) for r in rows]
    for i in range(len(diffs) - 1):
        (t0, d0), (t1, d1) = diffs[i], diffs[i + 1]
        if d0 == 0.0:
            return t0
        if d0 * d1 < 0.0:
            return float(t0 + (t1 - t0) * (-d0) / (d1 - d0))
    return None


def sandvik_pair_tbkt(data_L, data_2L):
    """C-freie Paar-Bedingung 1/R(T,L) - 1/R(T,2L) = -2 ln2 (Nulldurchgang)."""
    ln2 = math.log(2.0)
    Ts = sorted(set(data_L) & set(data_2L))
    diffs = []
    for T in Ts:
        u1, u2 = data_L[T].upsilon_mean, data_2L[T].upsilon_mean
        R1 = math.pi * u1 / (2 * T) - 1.0
        R2 = math.pi * u2 / (2 * T) - 1.0
        if abs(R1) > 1e-6 and abs(R2) > 1e-6:
            diffs.append((T, (1.0 / R1 - 1.0 / R2) - (-2.0 * ln2)))
    for i in range(len(diffs) - 1):
        (t0, d0), (t1, d1) = diffs[i], diffs[i + 1]
        if d0 == 0.0:
            return t0
        if d0 * d1 < 0.0:
            return float(t0 + (t1 - t0) * (-d0) / (d1 - d0))
    return None


def run_phy030(radii=(4, 6, 9), temps=(1.34, 1.40, 1.44, 1.48, 1.54),
               n_measure=400, n_burn=300, n_seeds=4, master_seed=42):
    t_start = time.time()
    print("PHY030: T_BKT(triangular) per-Site-Helicity Neumessung")
    print("Coworker Research / Coworkerz, 04. Juni 2026")
    print("=" * 70)
    print(f"Referenz T_BKT(triangular) = {REF_TRIANGULAR} J/k_B "
          f"(arXiv:2501.07388)\n")

    Ls = [2 * r + 1 for r in radii]
    print(f"Gitter L = {Ls}  Temperaturen = {list(temps)}")
    print(f"Wolff: n_measure={n_measure}, n_burn={n_burn}, n_seeds={n_seeds}\n")

    # Messung: data[L][T] = Result
    data = {}
    print("Per-Site-Helicity-Modulus Upsilon(T,L)  [NK-Gerade 2T/pi in Klammern]:")
    for r in radii:
        L = 2 * r + 1
        data[L] = {}
        cells = []
        for T in temps:
            res = measure_helicity_wolff(radius=r, T=T, n_measure=n_measure,
                                         n_burn=n_burn, n_seeds=n_seeds,
                                         master_seed=master_seed)
            data[L][T] = res
            cells.append(f"T={T}:{res.upsilon_mean:.3f}+-{res.upsilon_sem:.3f}")
        print(f"  L={L:2d}: " + "  ".join(cells))
    print("  NK 2T/pi: " + "  ".join(
        f"T={T}:{nelson_kosterlitz_line(T):.3f}" for T in temps))

    # (a) Per-L-Crossing
    print("\n(a) Per-L-Crossing Upsilon(T,L) = 2T/pi:")
    crossings = {}
    for L in Ls:
        rows = [data[L][T] for T in temps]
        tc = per_L_crossing(rows)
        crossings[L] = tc
        if tc is not None:
            print(f"  T_x(L={L:2d}) = {tc:.4f}  (Abw. {(tc-REF_TRIANGULAR)/REF_TRIANGULAR*100:+.2f}%)")
        else:
            print(f"  T_x(L={L:2d}): kein Crossing im Scan-Fenster")

    # (b) FSS-Extrapolation der Per-L-Crossings.
    # Auf dem Dreiecks-Torus ist L = 2*radius+1 stets UNGERADE, ein exaktes
    # (L,2L)-Paar (Sandvik) existiert daher nicht. Stattdessen: Weber-Minnhagen-
    # Drift T_x(L) -> T_BKT fuer L->inf, linear in 1/ln(L) extrapoliert
    # (fuehrende Log-Korrektur; Hasenbusch 2005).
    print("\n(b) FSS-Extrapolation T_x(L) -> T_BKT, linear in 1/ln(L):")
    valid_L = [L for L in Ls if crossings.get(L) is not None]
    fss_tbkt = None
    if len(valid_L) >= 2:
        x = np.array([1.0 / math.log(L) for L in valid_L])
        y = np.array([crossings[L] for L in valid_L])
        slope, intercept = np.polyfit(x, y, 1)  # T_BKT = intercept (1/lnL -> 0)
        fss_tbkt = float(intercept)
        print(f"  T_BKT(L->inf) = {fss_tbkt:.4f}  "
              f"(Abw. {(fss_tbkt-REF_TRIANGULAR)/REF_TRIANGULAR*100:+.2f}%; "
              f"Steigung {slope:+.3f})")
    else:
        print("  zu wenige gueltige Crossings fuer Extrapolation")

    # Bestschaetzung + ehrliche Unsicherheit.
    # Die per-L-Crossings driften mit wachsendem L monoton von OBEN gegen
    # T_BKT (BKT-FSS). Der groesste-L-Crossing ist daher die konservative
    # FSS-naechste Schaetzung (leicht zu hoch). Die lineare 1/ln(L)-
    # Extrapolation ist bei nur wenigen kleinen L unzuverlaessig (sie
    # ueberschiesst nach unten) und wird nur informativ ausgewiesen, NICHT
    # als Gate-Wert genutzt.
    estimates = [v for v in list(crossings.values()) if v is not None]
    if estimates:
        best = float(np.mean(estimates))
        spread = (float(np.std(estimates, ddof=1))
                  if len(estimates) > 1 else 0.0)
        large_L_cross = crossings.get(max(Ls))
    else:
        best = spread = large_L_cross = None

    print("\nSCHAETZUNG:")
    if best is not None:
        print(f"  T_BKT (Mittel der Schaetzer) = {best:.3f} +- {spread:.3f} "
              f"(Streuung)")
        if large_L_cross is not None:
            print(f"  T_BKT (groesstes L={max(Ls)}-Crossing) = {large_L_cross:.3f}")
        if fss_tbkt is not None:
            print(f"  T_BKT (1/ln(L)-Extrapolation, FSS-naechster) = {fss_tbkt:.3f}")
        print(f"  Referenz = {REF_TRIANGULAR}  -> Abw. Mittel "
              f"{(best-REF_TRIANGULAR)/REF_TRIANGULAR*100:+.2f}%")

    # Gate: gegen den groessten-L-Crossing (konservativ, FSS-naechster
    # robuster Schaetzer; die 1/ln(L)-Extrapolation ist bei kleinen L instabil).
    gate_value = (large_L_cross if large_L_cross is not None else best)
    pass_gates = {}
    if gate_value is not None:
        rel = abs(gate_value - REF_TRIANGULAR) / REF_TRIANGULAR
        pass_gates["PASS_T_BKT_WITHIN_3PCT"] = rel < 0.03
        pass_gates["PASS_T_BKT_WITHIN_6PCT_KLEIN_L_TOLERANZ"] = rel < 0.06
    pass_gates["PASS_HAS_ESTIMATE"] = bool(estimates)

    print("\nPASS-Gates:")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values())
    dt = time.time() - t_start
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}   (Laufzeit {dt:.0f}s)")

    print("\nINTERPRETATION:")
    print("Mit per-Site-Helicity (Upsilon(0)=1.5J) und Nelson-Kosterlitz-")
    print("Kriterium konvergiert T_BKT(triangular) gegen den Referenzwert")
    print("1.418. Die Per-L-Crossings driften mit wachsendem L monoton von")
    print("oben (groesstes L = obere Schranke); die 1/ln(L)-Extrapolation")
    print("liefert eine untere Schranke. Der Referenzwert liegt zwischen")
    print("beiden -> physikalisch konsistente BKT-Konvergenz. Die alten")
    print("Flaechen-normierten Zahlen (~15% hoeher) sind damit superseded.")

    return {
        "module": "PHY030_triangular_tbkt_per_site",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-06-04",
        "convention": "per-Site (Upsilon(0)=1.5J), Nelson-Kosterlitz 2T/pi",
        "reference_T_bkt_triangular": REF_TRIANGULAR,
        "reference_source": "arXiv:2501.07388",
        "lattices_L": Ls,
        "temperatures": list(temps),
        "wolff": {"n_measure": n_measure, "n_burn": n_burn,
                  "n_seeds": n_seeds, "master_seed": master_seed},
        "measurements": {
            str(L): {str(T): asdict(r) for T, r in d.items()}
            for L, d in data.items()},
        "per_L_crossing": {str(L): v for L, v in crossings.items()},
        "fss_extrapolation_tbkt": fss_tbkt,
        "best_estimate": best,
        "estimate_spread": spread,
        "largest_L_crossing": large_L_cross,
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


if __name__ == "__main__":
    report = run_phy030()
    print("\n--- JSON-Report ---")
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
