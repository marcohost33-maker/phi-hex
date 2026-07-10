"""PHY036 - honeycomb T_BKT: HKS-Extrapolation auf GROESSERER L-Leiter
(L bis 64; entscheidet, ob die Konvergenz bei L>=48 sichtbar wird).

PHY035 zeigte: auf L<=32 liegen die Paare flach (~0.60), keine saubere
Konvergenz gegen 0.573. HKS brauchen L=48..192. PHY036 schiebt die Leiter auf
L = 16/24/32/48/64 (drei verschachtelte Verdopplungspaare (16,32)/(24,48)/
(32,64), groesstes Paar jetzt (32,64) - groesser als PHY032/PHY035). Erst durch
die Kern-Vektorisierung (2026-06-09) ist L=64 (N=8192) lokal vertretbar.

Frage: bewegt das groessere L den Extrapolations-Limes Richtung 0.573/0.576?
Offenes, ehrliches Experiment - Ergebnis NICHT vorweggenommen.

ALLE Bausteine 1:1 aus PHY035 (parametrisiert) + PHY032 (Bootstrap) + PHY034
(Extrapolations-Primitive) wiederverwendet - single source of truth.

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
jackknife_tbkt_over_seeds = phy032.jackknife_tbkt_over_seeds
make_sandvik_pair_estimator = phy032.make_sandvik_pair_estimator

REF_MULTI = 0.573
REF_DEDICATED = 0.576
REF_MID = 0.5 * (REF_MULTI + REF_DEDICATED)

LADDER = (16, 24, 32, 48, 64)
DOUBLING_PAIRS = ((16, 32), (24, 48), (32, 64))
T_GRID = phy035.T_GRID


def run_phy036(n_seeds=8, n_measure=140, n_burn=60, master_seed=42):
    t0 = time.time()
    cube, means, sems = measure_cube(
        n_seeds=n_seeds, n_measure=n_measure, n_burn=n_burn,
        master_seed=master_seed, ladder=LADDER, t_grid=T_GRID)

    pairs = pair_estimates(means, pairs=DOUBLING_PAIRS, t_grid=T_GRID)
    pair_list = [pairs[p] for p in DOUBLING_PAIRS]
    all_finite = all(p is not None for p in pair_list)

    est = make_hks_extrap_estimator(pairs=DOUBLING_PAIRS, t_grid=T_GRID)
    boot = bootstrap_tbkt_over_seeds(cube, est, n_seeds, master_seed=master_seed)
    jack = jackknife_tbkt_over_seeds(cube, est, n_seeds)
    extrap = boot["point"]

    pair_boot = {}
    for (l1, l2) in DOUBLING_PAIRS:
        pe = make_sandvik_pair_estimator(l1, l2, T_GRID)
        pair_boot[(l1, l2)] = bootstrap_tbkt_over_seeds(
            cube, pe, n_seeds, master_seed=master_seed)

    ups0 = aligned_upsilon_zero()
    gate_a = abs(ups0 - SW_LIMIT_HONEYCOMB) < 1e-6
    # Cross-Check: L=24/48 teilen RNG-Vertrag mit PHY032 -> Paar (24,48)=0.5917.
    repro = (pairs[(24, 48)] is not None
             and abs(pairs[(24, 48)] - 0.5917) < 2e-3)
    ci_ok = boot["ci_lower"] is not None

    # FINDING (Physik, kein Build-Breaker): naehert L bis 64 die Referenz?
    bias_decreases = (all_finite
                      and abs(pair_list[-1] - REF_MID)
                      < abs(pair_list[0] - REF_MID))
    reduces = (extrap is not None and all_finite
               and abs(extrap - REF_MID) < abs(pair_list[-1] - REF_MID))
    ci_covers = (ci_ok and boot["ci_lower"] <= REF_DEDICATED
                 and boot["ci_upper"] >= REF_MULTI)

    gates = {
        "PASS_ALIGNED_EXACT_THREE_QUARTERS": gate_a,
        "PASS_THREE_PAIRS_FINITE": all_finite,
        "PASS_REPRODUCES_PHY032_PAIR_24_48": repro,
        "PASS_BOOTSTRAP_CI_COMPUTED": ci_ok,
    }
    overall = all(gates.values())
    findings = {"bias_decreases_with_L": bias_decreases,
                "extrap_reduces_offset": reduces,
                "ci_covers_reference": ci_covers}

    return {"elapsed_s": time.time() - t0, "n_seeds": n_seeds,
            "n_measure": n_measure, "n_burn": n_burn, "means": means,
            "sems": sems, "pairs": pairs, "pair_boot": pair_boot,
            "extrap": extrap, "boot": boot, "jack": jack, "ups0": ups0,
            "gates": gates, "findings": findings, "overall": overall}


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
    L.append("PHY036 v01 - honeycomb T_BKT: HKS-Extrapolation auf GROESSERER")
    L.append("            L-Leiter (L bis 64; Konvergenz-Test ueber PHY035)")
    L.append("Selbst-Audit / Coworkerz, 2026-06-11")
    L.append("=" * 72)
    L.append("")
    L.append(f"Leiter L = {list(LADDER)}  (N = 2 L^2, groesstes N = 8192)")
    L.append(f"Verschachtelte Verdopplungspaare: {list(DOUBLING_PAIRS)}")
    L.append(f"T-Gitter = {list(T_GRID)}")
    L.append(f"Wolff: n_seeds={rep['n_seeds']}, n_measure={rep['n_measure']}, "
             f"n_burn={rep['n_burn']}, master_seed=42")
    L.append(f"Gate A (T->0 exakt 0.75 J): gemessen {rep['ups0']:.6f}")
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
        dev = "" if tp is None else f"  {_pct(tp, REF_MULTI)}/{_pct(tp, REF_DEDICATED)}"
        lc = char_size(l1, l2, "geom")
        L.append(f"  T_BKT({l1},{l2}) = {tps}  {ci}{dev}  [Lc(geom)={lc:.2f}]")
    L.append("")
    L.append("--- HKS-L->inf-Extrapolation (3 Punkte, OLS in u=1/(lnLc)^2) ---")
    ex = rep["extrap"]
    b = rep["boot"]
    L.append(f"  T_BKT(L->inf) = {_fmt4(ex)}  {_pct(ex, REF_MULTI)}/"
             f"{_pct(ex, REF_DEDICATED)}")
    if b["ci_lower"] is not None:
        L.append(f"    Seed-Bootstrap: CI[{b['ci_lower']:.4f},"
                 f"{b['ci_upper']:.4f}] boot_std={b['boot_std']:.4f}")
    f = rep["findings"]
    L.append("")
    L.append("--- BEFUND ---")
    L.append(f"  Paar-Folge faellt mit L Richtung Referenz? "
             f"{'JA' if f['bias_decreases_with_L'] else 'NEIN'}")
    L.append(f"  Extrapolation im Referenzband [0.573,0.576]?  "
             f"{'JA' if f['ci_covers_reference'] else 'NEIN'}")
    L.append("  (Interpretation im PR/README; PHY036 erweitert PHY035 auf L=64.")
    L.append("   OVERALL=PASS = Analyse-Integritaet, NICHT die Physik-Hypothese.)")
    L.append("")
    L.append("--- PASS-Gates (Analyse-Integritaet) ---")
    for k, v in rep["gates"].items():
        L.append(f"  [{'PASS' if v else 'FAIL'}] {k}")
    L.append(f"  OVERALL: {'PASS' if rep['overall'] else 'FAIL'}")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    rep = run_phy036()
    out = (_SRC.parent / "results"
           / "260611 PHY036 honeycomb hks large ladder report.txt")
    write_report(rep, out)
    print(f"PHY036 Report -> {out}")
    b = rep["boot"]
    ci = (f"CI[{b['ci_lower']:.4f},{b['ci_upper']:.4f}]"
          if b["ci_lower"] is not None else "CI[n/a]")
    print(f"  T_BKT(L->inf) = {_fmt4(rep['extrap'])}  {ci}  (Refs 0.573/0.576)")
    print(f"  Laufzeit {rep['elapsed_s']:.1f}s  "
          f"OVERALL: {'PASS' if rep['overall'] else 'FAIL'}")
    return 0 if rep["overall"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
