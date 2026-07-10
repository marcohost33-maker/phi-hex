"""PHY038 - honeycomb T_BKT: HOCHSTATISTIK-Versuch der HKS-Aufloesung
(L bis 64, deutlich mehr Seeds/Sweeps als PHY036).

PHY036 zeigte: das (32,64)-Paar ist bei 8 Seeds rausch-dominiert (CI[0.598,
0.635]) und treibt die Extrapolation auf 0.639. PHY038 wiederholt EXAKT die
Leiter L=16/24/32/48/64, aber mit n_seeds=16 und mehr Sweeps - der ehrliche
Versuch, ob hoehere Statistik den Rausch-Sprung bei L=64 zaehmt und den
Extrapolations-Limes stabilisiert.

Erwartung NICHT vorweggenommen: entweder (i) der (32,64)-Punkt beruhigt sich
Richtung ~0.59 und die Extrapolation landet wieder ~0.57, ODER (ii) der Limes
wandert weiter -> dann ist honeycomb-T_BKT auch mit erhoehtem (aber lokalem)
Budget nicht robust extrapolierbar (das waere die endgueltige ehrliche Aussage).

ALLE Bausteine 1:1 aus PHY035 (parametrisiert) + PHY032 (Bootstrap) + PHY034
(Extrapolations-Primitive). Reiner Statistik-Hebel, KEINE neue Methode.

EVIDENZ: deterministischer Gate-Report nach results/ (seed=42 dokumentiert).
"""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

_SRC = Path(__file__).resolve().parent


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


phy032 = _load("phy032_honeycomb_wm_bootstrap",
               "260607 PHY032 honeycomb wm-logfit bootstrap v01.py")
phy035 = _load("phy035_hks_dense_ladder",
               "260611 PHY035 honeycomb hks dense ladder v01.py")

measure_cube = phy035.measure_cube
pair_estimates = phy035.pair_estimates
make_hks_extrap_estimator = phy035.make_hks_extrap_estimator
aligned_upsilon_zero = phy035.aligned_upsilon_zero
char_size = phy035.char_size
SW_LIMIT_HONEYCOMB = phy035.SW_LIMIT_HONEYCOMB
bootstrap_tbkt_over_seeds = phy032.bootstrap_tbkt_over_seeds
make_sandvik_pair_estimator = phy032.make_sandvik_pair_estimator

REF_MULTI = 0.573
REF_DEDICATED = 0.576
REF_MID = 0.5 * (REF_MULTI + REF_DEDICATED)

LADDER = (16, 24, 32, 48, 64)
DOUBLING_PAIRS = ((16, 32), (24, 48), (32, 64))
T_GRID = phy035.T_GRID

# Hochstatistik (vs PHY036: n_seeds=8, n_measure=140, n_burn=60).
N_SEEDS = 16
N_MEASURE = 240
N_BURN = 120

# Referenzwerte aus PHY036 (8 Seeds) fuer den direkten Statistik-Vergleich.
PHY036_PAIR_32_64 = 0.6178
PHY036_EXTRAP = 0.6388


def run_phy038(n_seeds=N_SEEDS, n_measure=N_MEASURE, n_burn=N_BURN,
               master_seed=42):
    t0 = time.time()
    cube, means, sems = measure_cube(
        n_seeds=n_seeds, n_measure=n_measure, n_burn=n_burn,
        master_seed=master_seed, ladder=LADDER, t_grid=T_GRID)
    pairs = pair_estimates(means, pairs=DOUBLING_PAIRS, t_grid=T_GRID)
    pair_list = [pairs[p] for p in DOUBLING_PAIRS]
    all_finite = all(p is not None for p in pair_list)

    est = make_hks_extrap_estimator(pairs=DOUBLING_PAIRS, t_grid=T_GRID)
    boot = bootstrap_tbkt_over_seeds(cube, est, n_seeds, master_seed=master_seed)
    extrap = boot["point"]

    pair_boot = {}
    for (l1, l2) in DOUBLING_PAIRS:
        pe = make_sandvik_pair_estimator(l1, l2, T_GRID)
        pair_boot[(l1, l2)] = bootstrap_tbkt_over_seeds(
            cube, pe, n_seeds, master_seed=master_seed)

    ups0 = aligned_upsilon_zero()
    gate_a = abs(ups0 - SW_LIMIT_HONEYCOMB) < 1e-6
    ci_ok = boot["ci_lower"] is not None
    # Hat hoehere Statistik das (32,64)-CI ggue. PHY036 verengt?
    pb = pair_boot[(32, 64)]
    width_38 = ((pb["ci_upper"] - pb["ci_lower"])
                if pb["ci_lower"] is not None else None)

    gates = {
        "PASS_ALIGNED_EXACT_THREE_QUARTERS": gate_a,
        "PASS_THREE_PAIRS_FINITE": all_finite,
        "PASS_BOOTSTRAP_CI_COMPUTED": ci_ok,
    }
    findings = {
        "extrap_in_reference_band": (
            ci_ok and boot["ci_lower"] <= REF_DEDICATED
            and boot["ci_upper"] >= REF_MULTI),
        "extrap_closer_than_phy036": (
            extrap is not None
            and abs(extrap - REF_MID) < abs(PHY036_EXTRAP - REF_MID)),
        "pair_32_64_ci_width": width_38,
    }
    return {"elapsed_s": time.time() - t0, "n_seeds": n_seeds,
            "n_measure": n_measure, "n_burn": n_burn, "means": means,
            "sems": sems, "pairs": pairs, "pair_boot": pair_boot,
            "extrap": extrap, "boot": boot, "ups0": ups0, "gates": gates,
            "findings": findings, "overall": all(gates.values())}


def _pct(x, ref):
    if x is None:
        return "n/a"
    return f"{(x - ref) / ref * 100:+.2f}%"


def _fmt4(x):
    """None-sichere Formatierung (Code-Audit M3, 2026-07-10): der FAIL-Pfad
    muss einen Report schreiben koennen statt mit TypeError zu crashen."""
    return "None" if x is None else f"{x:.4f}"


def write_report(rep, path):
    L = []
    L.append("PHY038 v01 - honeycomb T_BKT: HOCHSTATISTIK-HKS (L bis 64,")
    L.append(f"            n_seeds={rep['n_seeds']}, n_measure={rep['n_measure']},"
             f" n_burn={rep['n_burn']}; vs PHY036: 8/140/60)")
    L.append("Selbst-Audit / Coworkerz, 2026-06-12")
    L.append("=" * 72)
    L.append("")
    L.append(f"Leiter L = {list(LADDER)}  Paare {list(DOUBLING_PAIRS)}")
    L.append(f"T-Gitter = {list(T_GRID)}  master_seed=42")
    L.append(f"Gate A (T->0 exakt 0.75 J): {rep['ups0']:.6f}")
    L.append(f"Laufzeit: {rep['elapsed_s']:.1f}s")
    L.append("")
    L.append("--- Messung Upsilon(T,L) (Seed-Mittel +- SEM) ---")
    for Ll in LADDER:
        cells = "  ".join(
            f"T={T:.4f}:{rep['means'][Ll][T]:.4f}+-{rep['sems'][Ll][T]:.4f}"
            for T in T_GRID)
        L.append(f"  L={Ll}: {cells}")
    L.append("")
    L.append("--- C-eliminierte Verdopplungspaare (+ Seed-Bootstrap-CI) ---")
    for (l1, l2) in DOUBLING_PAIRS:
        tp = rep["pairs"][(l1, l2)]
        b = rep["pair_boot"][(l1, l2)]
        ci = (f"CI[{b['ci_lower']:.4f},{b['ci_upper']:.4f}]"
              if b["ci_lower"] is not None else "CI[n/a]")
        tps = "None" if tp is None else f"{tp:.4f}"
        L.append(f"  T_BKT({l1},{l2}) = {tps}  {ci}  [Lc={char_size(l1,l2,'geom'):.1f}]")
    L.append(f"  (PHY036-Vergleich (32,64): {PHY036_PAIR_32_64:.4f} bei 8 Seeds)")
    L.append("")
    ex = rep["extrap"]
    b = rep["boot"]
    L.append("--- HKS-L->inf-Extrapolation (3 Punkte, u=1/(lnLc)^2) ---")
    L.append(f"  T_BKT(L->inf) = {_fmt4(ex)}  {_pct(ex, REF_MULTI)}/"
             f"{_pct(ex, REF_DEDICATED)}")
    if b["ci_lower"] is not None:
        L.append(f"    Seed-Bootstrap: CI[{b['ci_lower']:.4f},{b['ci_upper']:.4f}]"
                 f" boot_std={b['boot_std']:.4f}")
    L.append(f"  (PHY036 8-Seed-Extrapolation war {PHY036_EXTRAP:.4f})")
    L.append("")
    f = rep["findings"]
    L.append("--- BEFUND ---")
    L.append(f"  Extrapolation im Referenzband [0.573,0.576]? "
             f"{'JA' if f['extrap_in_reference_band'] else 'NEIN'}")
    L.append(f"  Naeher an Referenz als PHY036 (8 Seeds)?      "
             f"{'JA' if f['extrap_closer_than_phy036'] else 'NEIN'}")
    if f["pair_32_64_ci_width"] is not None:
        L.append(f"  (32,64)-CI-Breite: {f['pair_32_64_ci_width']:.4f} "
                 f"(PHY036: ~0.0375)")
    L.append("  Interpretation im PR/README. OVERALL=PASS = Analyse-Integritaet")
    L.append("  (Gate A, Paare endlich, CI berechnet), NICHT Physik-Hypothese.")
    L.append("")
    L.append("--- PASS-Gates ---")
    for k, v in rep["gates"].items():
        L.append(f"  [{'PASS' if v else 'FAIL'}] {k}")
    L.append(f"  OVERALL: {'PASS' if rep['overall'] else 'FAIL'}")
    L.append("")
    L.append("--- Caveats ---")
    L.append("  - 16 Seeds / L=64 sind weiterhin WEIT unter Publikationsniveau")
    L.append("    (dort L=48..192, tausende Messungen). Lokaler Best-Effort.")
    L.append("  - Reiner Statistik-Hebel ueber PHY036; gleiche Methode/Leiter.")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    rep = run_phy038()
    out = (_SRC.parent / "results"
           / "260612 PHY038 honeycomb hks highstat report.txt")
    write_report(rep, out)
    print(f"PHY038 Report -> {out}")
    b = rep["boot"]
    ci = (f"CI[{b['ci_lower']:.4f},{b['ci_upper']:.4f}]"
          if b["ci_lower"] is not None else "CI[n/a]")
    print(f"  T_BKT(L->inf) = {_fmt4(rep['extrap'])}  {ci}  "
          f"(PHY036 war {PHY036_EXTRAP})")
    print(f"  Laufzeit {rep['elapsed_s']:.1f}s  "
          f"OVERALL: {'PASS' if rep['overall'] else 'FAIL'}")
    return 0 if rep["overall"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
