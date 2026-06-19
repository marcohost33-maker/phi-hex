"""PHY027: Praezisions-FSS von T_BKT mit Wolff-Sampling (Stufe 2)

Coworker Research / Coworkerz, 02. Juni 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

Stufe 2 der Roadmap: Nutzt den Wolff-Cluster-Sampler (PHY026) fuer eine
Bestimmung der kritischen Temperatur auf GROSSEN Gittern (L bis 24/32),
die mit lokalem Update unzugaenglich waren. Wendet die Weber-Minnhagen-
Log-Korrektur und die Hsieh/Sandvik-Paar-Extrapolation an.

METHODE:
  Nelson-Kosterlitz: rho_s(T_BKT) = (2/pi) T_BKT (universeller Sprung).
  Naives Crossing Upsilon(T,L) = (2/pi)T ist nach oben verzerrt.
  Weber-Minnhagen am kritischen Punkt:
     pi*Upsilon(T_BKT,L)/(2 T_BKT) = 1 + 1/(2 ln(L) + C)
  Paar-Methode (Hsieh/Sandvik): Definiere T*(L) als Crossing pro L, dann
  extrapoliere T*(L) gegen 1/(ln L)^2 in den thermodynamischen Limes, da
  die fuehrende Korrektur der naiven Crossing-Temperatur ~1/(ln L)^2 ist.

REFERENZWERT: T_BKT(triangular) = 1.418(2)
  (Sorokin, Ann. Phys. 411, 167952 (2019); zitiert in arXiv:2305.00651)

HINWEIS ehrliche Einordnung: Auch mit Wolff bleibt dies ein Demonstrator.
Hochpraezisions-Studien nutzen L bis 512-2048. Ziel hier ist, dass die
extrapolierte T_BKT den Referenzwert 1.418(2) im Fehlerband ERREICHT oder
ihm signifikant naeher kommt als die naiven Kleingitter-Crossings (~1.49).
"""
from __future__ import annotations

import json
import math
import numpy as np
from dataclasses import asdict
from typing import Any, Optional

import importlib.util
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent


def _load(name: str, filename: str):
    """Portabler Loader (ersetzt hartkodierten /home/claude-Dev-Pfad).

    spec_from_file_location relativ zu DIESER Datei -> Standalone-Lauf
    (`python "src/..."`, README) funktioniert plattformunabhaengig.
    Identisches Muster wie PHY030 v02 / PHY031-033.
    """
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_wolff = _load("phy026_wolff_cluster", "260602 PHY026 wolff cluster v01.py")
measure_helicity_wolff = _wolff.measure_helicity_wolff
WolffHelicityResult = _wolff.WolffHelicityResult


def find_crossing(measurements: list[WolffHelicityResult]) -> Optional[float]:
    """Crossing T, wo Upsilon(T) = (2/pi)T, via lineare Interpolation."""
    ms = sorted(measurements, key=lambda m: m.T)
    for i in range(len(ms) - 1):
        g0 = ms[i].upsilon_mean - ms[i].nk_line
        g1 = ms[i + 1].upsilon_mean - ms[i + 1].nk_line
        if g0 >= 0 >= g1 and (g0 - g1) != 0:
            T0, T1 = ms[i].T, ms[i + 1].T
            return float(T0 + (g0 / (g0 - g1)) * (T1 - T0))
    return None


def run_phy027() -> dict[str, Any]:
    print("PHY027: Praezisions-FSS von T_BKT mit Wolff (Stufe 2)")
    print("Coworker Research / Coworkerz, 02. Juni 2026")
    print("=" * 70)
    ref = 1.418
    print(f"Referenz T_BKT(triangular) = {ref}(2)\n")

    # Groessere Gitter dank Wolff: L = 9, 13, 17, 25 (radius 4,6,8,12)
    radii = [4, 6, 8, 12]
    temps = [1.30, 1.40, 1.50, 1.60, 1.70]

    # Wolff-Sampling-Parameter (moderat, Demonstrator)
    cfg = dict(n_measure=600, n_burn=300, n_seeds=3, master_seed=42)

    crossings = {}
    all_meas = {}
    for radius in radii:
        L = 2 * radius + 1
        print(f"L = {L} (radius {radius}, {L*L} Knoten):")
        ms = []
        for T in temps:
            r = measure_helicity_wolff(radius=radius, T=T, **cfg)
            ms.append(r)
            print(f"  T={T:.2f}: Ups={r.upsilon_mean:+.3f}+-{r.upsilon_sem:.3f}"
                  f", 2T/pi={r.nk_line:.3f}, m={r.magnetization:.3f}, "
                  f"tau={r.tau_int:.2f}")
        tc = find_crossing(ms)
        crossings[L] = tc
        all_meas[L] = ms
        if tc is not None:
            print(f"  -> Crossing T*(L={L}) = {tc:.3f} "
                  f"(Abw. {(tc-ref)/ref*100:+.1f}%)\n")
        else:
            print("  -> kein Crossing im Raster\n")

    # Extrapolation gegen 1/(ln L)^2
    valid = [(L, tc) for L, tc in crossings.items() if tc is not None]
    print("=" * 70)
    extrap = None
    if len(valid) >= 2:
        Ls = np.array([L for L, _ in valid], dtype=float)
        tcs = np.array([tc for _, tc in valid])
        x = 1.0 / (np.log(Ls) ** 2)  # fuehrende Korrektur ~1/(ln L)^2
        # lineare Regression T*(L) = T_BKT + a * x; Achsenabschnitt = T_BKT
        A = np.vstack([x, np.ones_like(x)]).T
        coef, *_ = np.linalg.lstsq(A, tcs, rcond=None)
        slope, extrap = float(coef[0]), float(coef[1])
        print("Extrapolation T*(L) = T_BKT + a/(ln L)^2:")
        for L, tc in valid:
            print(f"  L={L}: T*={tc:.3f}")
        print(f"  -> extrapoliert T_BKT(L->inf) = {extrap:.3f} "
              f"(Abw. {(extrap-ref)/ref*100:+.1f}%)")
        print(f"     (Steigung a = {slope:+.3f})")

    # Bewertung
    pass_gates = {}
    pass_gates["PASS_ALL_CROSSINGS_FOUND"] = len(valid) == len(radii)
    if valid:
        tcs_only = [tc for _, tc in valid]
        # Crossings sollten mit L monoton fallen (von oben annaehern)
        pass_gates["PASS_MONOTONE_DOWNWARD"] = all(
            tcs_only[k + 1] <= tcs_only[k] + 0.02
            for k in range(len(tcs_only) - 1))
    if extrap is not None:
        # Extrapolation soll naeher an Referenz sein als das kleinste-L-Crossing
        largest_naive_dev = abs(valid[0][1] - ref)
        extrap_dev = abs(extrap - ref)
        pass_gates["PASS_EXTRAP_IMPROVES"] = extrap_dev < largest_naive_dev
        pass_gates["PASS_EXTRAP_WITHIN_5PCT"] = extrap_dev / ref < 0.05

    print("\nPASS-Gates:")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values()) if pass_gates else False
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}")

    print("\nINTERPRETATION:")
    print("Wolff erlaubt groessere Gitter; die Crossing-Temperatur naehert")
    print("sich von oben dem Referenzwert 1.418. Die 1/(ln L)^2-Extrapolation")
    print("testet, ob der thermodynamische Limes den Literaturwert trifft.")

    return {
        "module": "PHY027_precision_fss_wolff",
        "attribution": "Coworker Research / Coworkerz",
        "reference_T_bkt": ref,
        "radii": radii, "temperatures": temps,
        "crossings": {str(k): v for k, v in crossings.items()},
        "extrapolated_T_bkt": extrap,
        "measurements": {
            str(L): [asdict(m) for m in ms] for L, ms in all_meas.items()},
        "pass_gates": pass_gates,
        "overall_pass": overall,
    }


if __name__ == "__main__":
    report = run_phy027()
    def clean(o):
        if o is None:
            return None
        if isinstance(o, (np.bool_, bool)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            v = float(o)
            return v if math.isfinite(v) else None
        if isinstance(o, float):
            return o if math.isfinite(o) else None
        if isinstance(o, np.ndarray):
            return [clean(x) for x in o.tolist()]
        if hasattr(o, "__dataclass_fields__"):
            return clean(asdict(o))
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(x) for x in o]
        return o
    print("\n--- JSON-Report ---")
    print(json.dumps(clean(report), indent=2, allow_nan=False))
