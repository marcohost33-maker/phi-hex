"""PHY042 - honeycomb WL-FSS L=24/32/48: die in PHY041 (Spec + README
"Naechste Stufe nach PR #18") vertraglich benannten Punkte 2-4.

Coworker Research / Coworkerz, 6. Juli 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

================================ ZWECK =================================
PHY041 hat den gitter-agnostischen Wang-Landau/1-t-Apparat auf honeycomb
bei L=12/16/24 validiert und als FINDING festgehalten: die C-eliminierten
WM-Paare STEIGEN mit L (0.5951/0.6029/0.6087). PHY042 beantwortet die
vertraglich folgende Frage: was machen die Paar-Schaetzer bei L=32/48,
und liefern Y2- und Y4-Kanal getrennt konsistente Aussagen?

Vertrag: `spec/260706 PHI HEX phy042 honeycomb wl fss method v01.md`.
Referenz-Einordnung NUR ueber das Band aus
`spec/260703 PHI HEX honeycomb reference conventions audit v01.md`.

================================ METHODE ===============================
Kernel + Analyse 1:1 aus PHY041/PHY040 importiert (single source of truth,
kein Reinvent). Delta zu PHY041:

1. GROESSERE L: (24, 32, 48). L=24 laeuft mit exakt den PHY041-Parametern
   (master_seed=42, Streams 600/650, lnf_final=1e-5, prod_sweeps=30000,
   identisches T-Gitter) -> bitgleiche Reproduktion des committed
   PHY041-Stands als deterministische Bruecke (VAL-B/VAL-C).
2. SKALIERTE PRODUKTION: prod_sweeps = max(30000, 60*nbins) fuer L=32/48,
   damit die Bin-Belegungsdichte (~60 Samples/Bin, PHY041-Niveau bei L=24)
   nicht ausduennt. Die Fenster werden im Parent EINMAL je L bestimmt und
   an alle Walker desselben L uebergeben (identisches Binning).
3. MULTI-WALKER-SYSTEMATIK: bei L=32 laufen 3 unabhaengige g(E)-Walker
   (Walker 0: Standard-Stream 650+L; Walker w>=1: Stream 90000+1000*w+L,
   kollisionsfrei zu 400+/500+/600+/650+/660+/700+ inkl. +1000L-Termen).
   Der Y2-Spread ueber Walker ist der Sampler-Systematik-Schaetzer und
   wird per Paar-Neuberechnung je Walker auf T_BKT propagiert.
4. GETRENNTE KANAELE: Y2 (WM-Paare) und Y4 (Dip-Lage je L) werden separat
   ausgewertet und separat gegen das Referenzband eingeordnet.

Unabhaengige Jobs (L, Walker) laufen parallel (Prozess-Pool, fork);
jede Job-RNG ist vollstaendig durch (master_seed, stream) bestimmt -
Determinismus unabhaengig vom Scheduling.

EHRLICHKEIT: PASS = Pipeline-Integritaet (Orakel, Leak, Quervalidierung,
Glaette, Walker-Spread). Der Physik-Befund (Paar-Trend, Band-Lage,
Y4-Dip-Trend) ist FINDING, kein Build-Breaker (AGENTS.md: Negativ-Results
sind Buerger erster Klasse). Die Walker-Systematik misst NUR den Sampler,
nicht den finite-size-Bias. Kein "neuer T_BKT-Bestwert" ohne Review.

Provenance: Wang & Landau PRL 86, 2050 (2001); Belardinelli & Pereyra JCP
127, 184105 (2007); WM-Paare PRB 37, 5986(R) (1988), arXiv:1302.2900;
Methodik-Vorbild honeycomb (WL + Y2/Y4): arXiv:2406.12076.
"""
from __future__ import annotations

import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import importlib.util
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent


def _load(name: str, filename: str):
    """Portabler Loader (Muster PHY030 v02 / PHY031-033 / PHY040/041)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_p41 = _load("phy041_honeycomb_wl",
             "260702 PHY041 honeycomb wang-landau entropic helicity v01.py")

# Kernel + Analyse 1:1 aus PHY041/PHY040 (single source of truth).
honeycomb_arrays = _p41.honeycomb_arrays
wolff_anchor = _p41.wolff_anchor
window_from_anchors = _p41.window_from_anchors
wl_entropic_lattice = _p41.wl_entropic_lattice
canonical_edge_leak = _p41.canonical_edge_leak
parse_phy032_grid = _p41.parse_phy032_grid
wolff_reference = _p41.wolff_reference
aligned_upsilon_zero = _p41.aligned_upsilon_zero
upsilon_curves = _p41.upsilon_curves
tbkt_pair_from_curves = _p41.tbkt_pair_from_curves
REF_BAND = {  # T-Form, aus spec/260703 (beta-Kanaele konvertiert)
    "multi_lattice": 0.573,
    "helicity_direct": (0.571, 0.008),
    "nn_mc": (0.560, 0.009),
    "upsilon_beta": (0.5928, 0.0011),
    "upsilon4_beta": (0.6116, 0.0041),
    "binder_beta": (0.5800, 0.0007),
}

# PHY041-Bruecke: committed Reportwerte (results/260702 PHY041 ... report.txt)
PHY041_PAIRS = {(12, 16): 0.5951, (12, 24): 0.6029, (16, 24): 0.6087}
PHY041_Y4_DIP_L24 = 0.6500

# RNG-Stream-Vertrag PHY042: Walker 0 nutzt die PHY041-Streams (600er Anker,
# 650+L WL) fuer die deterministische L=24-Bruecke; Zusatz-Walker w>=1 nutzen
# 90000+1000*w (+L in wl_entropic_lattice) - kollisionsfrei zu allen
# bisherigen Vertraegen (400+/500+/600+/650+/660+/700+ inkl. +1000L-Terme).
_STREAM_WALKER_BASE = 90000

_T_GRID = np.round(np.arange(0.52, 0.6701, 0.005), 4)  # identisch PHY041
_CORE_WIN = (0.54, 0.66)


def _prod_sweeps_for(nbins: int, L: int) -> int:
    """Skalierte Produktion (~60 Samples/Bin); L=24 exakt PHY041 (30000)."""
    if L == 24:
        return 30000
    return max(30000, 60 * nbins)


def _stream_for(walker: int) -> int:
    """RNG-Stream je Walker: 0 -> PHY041-Standard (650, +L im Kernel);
    w>=1 -> 90000+1000*w (kollisionsfrei, siehe Stream-Vertrag oben)."""
    if walker == 0:
        return _p41._STREAM_WL
    return _STREAM_WALKER_BASE + 1000 * walker


def _wl_job(args: tuple) -> tuple:
    """Ein (L, walker)-Job: WL + Produktion im vorbestimmten Fenster.
    Top-level (picklebar); RNG rein (seed, stream)-bestimmt."""
    L, walker, e_lo, e_hi, master_seed = args
    _, nbr_list, ei, ej, ax, ay, n = honeycomb_arrays(L)
    nbins = int(round((e_hi - e_lo) * n))
    stream = _stream_for(walker)
    res = wl_entropic_lattice(nbr_list, ei, ej, ax, ay, n, L, e_lo, e_hi,
                              lnf_final=1e-5,
                              prod_sweeps=_prod_sweeps_for(nbins, L),
                              seed=master_seed, stream=stream, verbose=False)
    return (L, walker, res)


def run_phy042(Ls=(24, 32, 48), n_walkers_L32=3, master_seed=42,
               max_workers=4) -> dict:
    t_grid = _T_GRID
    print("PHY042: honeycomb WL-FSS L=24/32/48 + Y2/Y4-Kanaele + Multi-Walker")
    print("Coworker Research / Coworkerz, 6. Juli 2026")
    print("=" * 72)
    t_start = time.time()

    # Gate A: Geometrie-/Aggregat-Orakel (exakt, arXiv-unabhaengig)
    y0 = aligned_upsilon_zero()
    aligned_ok = abs(y0 - 0.75) < 1e-9
    print(f"\n[GATE A] aligned Upsilon_2(0) = {y0:.6f} (exakt 0.75) "
          f"-> {'PASS' if aligned_ok else 'FAIL'}")

    # --- Fenster je L EINMAL im Parent (identisches Binning je L) ---------
    windows: dict[int, dict] = {}
    for L in Ls:
        m_lo, s_lo = wolff_anchor(L, 0.50, master_seed=master_seed)
        m_hi, s_hi = wolff_anchor(L, 0.70, master_seed=master_seed)
        e_lo, e_hi = window_from_anchors(m_lo, s_lo, m_hi, s_hi)
        n = 2 * L * L
        nbins = int(round((e_hi - e_lo) * n))
        windows[L] = {"window_ps": [e_lo, e_hi], "nbins": nbins,
                      "prod_sweeps": _prod_sweeps_for(nbins, L),
                      "anchor_lo": {"T": 0.50, "mean_ps": m_lo, "std_ps": s_lo},
                      "anchor_hi": {"T": 0.70, "mean_ps": m_hi, "std_ps": s_hi}}
        print(f"[WIN] L={L} (N={n}): [{e_lo:.4f}, {e_hi:.4f}] "
              f"bins={nbins} prod={windows[L]['prod_sweeps']}")

    # --- WL-Jobs: (L, walker) parallel, deterministisch je (seed, stream) --
    jobs = []
    for L in Ls:
        e_lo, e_hi = windows[L]["window_ps"]
        n_walk = n_walkers_L32 if L == 32 else 1
        for w in range(n_walk):
            jobs.append((L, w, e_lo, e_hi, master_seed))
    # laengste Jobs zuerst (L=48 dominiert die Wall-Zeit)
    jobs.sort(key=lambda j: (-j[0], j[1]))
    print(f"\n[WL] {len(jobs)} Jobs (L, walker): "
          f"{[(j[0], j[1]) for j in jobs]} (max_workers={max_workers})")
    results: dict[tuple[int, int], object] = {}
    if max_workers > 1:
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            for L, w, res in ex.map(_wl_job, jobs):
                results[(L, w)] = res
                print(f"    fertig: L={L} walker={w} wl_sweeps={res.wl_sweeps} "
                      f"covered={int(res.mask.sum())}/{len(res.mask)} "
                      f"t={time.time() - t_start:.0f}s")
    else:
        for job in jobs:
            L, w, res = _wl_job(job)
            results[(L, w)] = res
            print(f"    fertig: L={L} walker={w} wl_sweeps={res.wl_sweeps} "
                  f"covered={int(res.mask.sum())}/{len(res.mask)} "
                  f"t={time.time() - t_start:.0f}s")

    curves = {(L, w): upsilon_curves(results[(L, w)], t_grid)
              for (L, w) in results}
    main_curve = {L: curves[(L, 0)] for L in Ls}

    # --- Leak-Gate: kanonisches Gewicht an den Fensterraendern ------------
    leak_max = 0.0
    for key, res in results.items():
        for T in t_grid:
            leak_max = max(leak_max, canonical_edge_leak(res, float(T)))
    leak_ok = leak_max < 1e-3
    print(f"\n[LEAK] max. Randgewicht ueber alle (L, walker, T): "
          f"{leak_max:.2e} (Gate < 1e-3) -> {'PASS' if leak_ok else 'FAIL'}")

    # --- VAL-A: frische Wolff-Referenz bei NEUEM L=32 ----------------------
    print("\n[VAL-A] WL vs frische Wolff-Referenz (L=32):")
    e_ok, y2_ok = True, True
    val_rows = []
    for T in (0.55, 0.60, 0.65):
        c = main_curve[32]
        k = int(np.argmin(np.abs(c["T"] - T)))
        Ewl = c["E"][k] / results[(32, 0)].n
        y2wl = c["y2"][k]
        ref = wolff_reference(32, T, master_seed=master_seed)
        de = abs(Ewl - ref["E_ps"])
        dy = abs(y2wl - ref["y2"])
        e_ok = e_ok and de < 0.03
        y2_ok = y2_ok and dy < 0.04
        val_rows.append({"T": T, "E_wl": Ewl, "E_wolff": ref["E_ps"],
                         "y2_wl": y2wl, "y2_wolff": ref["y2"],
                         "y2_wolff_sem": ref["y2_sem"]})
        print(f"      T={T:.3f}  <E>/N WL={Ewl:+.4f} Wolff={ref['E_ps']:+.4f} "
              f"(d={de:.4f})  Y2 WL={y2wl:.4f} "
              f"Wolff={ref['y2']:.4f}+-{ref['y2_sem']:.4f} (d={dy:.4f})")

    # --- VAL-B: Drift-Guard gegen committed PHY032-Gitter (L=24) ----------
    print("\n[VAL-B] WL vs committed PHY032-Messgitter (L=24, Drift-Guard):")
    grid032 = parse_phy032_grid()
    grid_ok = True
    grid_rows = []
    for (T, u032, sem032) in grid032.get(24, []):
        if not (t_grid[0] <= T <= t_grid[-1]):
            continue
        y2wl = float(np.interp(T, main_curve[24]["T"], main_curve[24]["y2"]))
        tol = max(4.0 * sem032, 0.025)
        ok = abs(y2wl - u032) <= tol
        grid_ok = grid_ok and ok
        grid_rows.append({"L": 24, "T": T, "y2_wl": y2wl, "y2_phy032": u032,
                          "sem_phy032": sem032, "ok": ok})
        print(f"      L=24 T={T:.4f}  WL={y2wl:.4f} "
              f"PHY032={u032:.4f}+-{sem032:.4f} "
              f"(d={abs(y2wl - u032):.4f}, tol={tol:.4f}) "
              f"{'ok' if ok else 'ABWEICHUNG'}")

    # --- VAL-C: PHY041-Bruecke (bitgleiche L=24-Reproduktion) --------------
    c24 = main_curve[24]
    k24 = int(np.argmin(c24["y4_scaled"]))
    y4_dip_24 = float(c24["T"][k24])
    bridge_ok = abs(y4_dip_24 - PHY041_Y4_DIP_L24) < 1e-9
    print(f"\n[VAL-C] PHY041-Bruecke: Y4-Dip L=24 = {y4_dip_24:.4f} "
          f"(committed {PHY041_Y4_DIP_L24:.4f}) "
          f"-> {'PASS' if bridge_ok else 'FAIL'}")

    # --- Y2-Kanal: C-eliminierte Paare auf glatten Kurven ------------------
    print("\n[Y2-FSS] C-eliminierte WM-Paare (Hauptwalker):")
    pairs = [(Ls[i], Ls[j]) for i in range(len(Ls)) for j in range(len(Ls))
             if Ls[j] > Ls[i]]
    pair_tbkt = {}
    for La, Lb in pairs:
        tb = tbkt_pair_from_curves(t_grid, main_curve[La]["y2"], La,
                                   main_curve[Lb]["y2"], Lb)
        pair_tbkt[(La, Lb)] = tb
        if tb is not None:
            print(f"      T_BKT({La},{Lb}) = {tb:.4f}  "
                  f"(vs Band: 0.560..0.6116; Multi-Lattice 0.573: "
                  f"{(tb - 0.573) / 0.573 * 100:+.2f}%)")
        else:
            print(f"      T_BKT({La},{Lb}): kein Nulldurchgang im T-Fenster")

    # Trend gegen PHY041-Paare (steigend = Negativ-Result fuer Paar-Lesart)
    phy41_seq = [PHY041_PAIRS[p] for p in sorted(PHY041_PAIRS, key=min)]
    valid_pairs = {p: v for p, v in pair_tbkt.items() if v is not None}
    pair_seq = [valid_pairs[p]
                for p in sorted(valid_pairs, key=lambda p: (min(p), max(p)))]
    trend_up = (len(pair_seq) >= 2
                and all(b > a for a, b in zip(pair_seq, pair_seq[1:]))
                and (not pair_seq or pair_seq[0] > phy41_seq[-1] - 0.005))
    print(f"      PHY041-Paare (L<=24): "
          f"{', '.join(f'{v:.4f}' for v in phy41_seq)}")

    # --- Y4-Kanal: Dip-Lage je L -------------------------------------------
    print("\n[Y4-FSS] N*Upsilon_4-Dip je L (deterministisch, rauschfrei):")
    y4_dips = {}
    for L in Ls:
        c = main_curve[L]
        k = int(np.argmin(c["y4_scaled"]))
        y4_dips[L] = {"T_dip": float(c["T"][k]),
                      "depth": float(c["y4_scaled"][k]),
                      "at_edge": bool(k == 0 or k == len(c["T"]) - 1)}
        print(f"      L={L}: Dip bei T={y4_dips[L]['T_dip']:.4f} "
              f"(Tiefe {y4_dips[L]['depth']:.0f}"
              f"{', am Gitterrand' if y4_dips[L]['at_edge'] else ''})")

    # --- Multi-Walker-Systematik (L=32) -------------------------------------
    print(f"\n[WALKER] {n_walkers_L32} unabhaengige g(E)-Walker (L=32):")
    core = (t_grid >= _CORE_WIN[0]) & (t_grid <= _CORE_WIN[1])
    y2_stack = np.stack([curves[(32, w)]["y2"] for w in range(n_walkers_L32)])
    spread = float((y2_stack[:, core].max(axis=0)
                    - y2_stack[:, core].min(axis=0)).max())
    walker_ok = spread < 0.02
    print(f"      max. Y2-Spread im Kernfenster {list(_CORE_WIN)}: "
          f"{spread:.4f} (Gate < 0.02) -> {'PASS' if walker_ok else 'FAIL'}")
    # Propagation auf Paar-T_BKT: Paare je Walker neu berechnen
    walker_pairs = {}
    for (La, Lb) in [(24, 32), (32, 48)]:
        vals = []
        for w in range(n_walkers_L32):
            ca = curves[(La, w)] if La == 32 else main_curve[La]
            cb = curves[(Lb, w)] if Lb == 32 else main_curve[Lb]
            tb = tbkt_pair_from_curves(t_grid, ca["y2"], La, cb["y2"], Lb)
            if tb is not None:
                vals.append(tb)
        if vals:
            walker_pairs[(La, Lb)] = {"values": vals,
                                      "spread": max(vals) - min(vals)}
            print(f"      T_BKT({La},{Lb}) ueber Walker: "
                  f"{', '.join(f'{v:.4f}' for v in vals)} "
                  f"(Spread {max(vals) - min(vals):.4f})")

    # --- Glaette-Gate (L=48, Kernfenster) -----------------------------------
    c48 = main_curve[48]["y2"][core]
    y2_smooth = bool(np.all(np.diff(c48) < 0))
    print(f"\n[Y2] Upsilon_2(T) L=48 Kernfenster {list(_CORE_WIN)}: "
          f"streng monoton fallend={y2_smooth} "
          f"({c48[0]:.3f} .. {c48[-1]:.3f})")

    gates = {
        "PASS_ALIGNED_EXACT_THREE_QUARTERS": aligned_ok,
        "PASS_NO_CANONICAL_EDGE_LEAK": leak_ok,
        "PASS_WL_ENERGY_MATCHES_WOLFF_L32": e_ok,
        "PASS_WL_Y2_MATCHES_WOLFF_L32": y2_ok,
        "PASS_WL_Y2_MATCHES_PHY032_GRID_L24": grid_ok,
        "PASS_PHY041_BRIDGE_Y4_DIP_L24": bridge_ok,
        "PASS_Y2_SMOOTH_NOISE_FREE_L48": y2_smooth,
        "PASS_WALKER_SPREAD_BELOW_GATE": walker_ok,
    }
    print("\n[GATES]")
    for k, v in gates.items():
        print(f"      [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(gates.values())
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}  "
          f"(= Pipeline-Integritaet; Physik-Befund ist FINDING, "
          f"kein Build-Breaker)")

    # --- ehrlicher Physik-Befund -------------------------------------------
    print("\n[FINDING] (ehrlich, NICHT getunt)")
    if trend_up:
        print("    - NEGATIV-RESULT (Y2-Kanal): der C-eliminierte "
              "Paar-Schaetzer steigt AUCH bei L=24..48 weiter mit L -")
        print("      die kleine-L-Paar-Lesart extrapoliert nicht von selbst "
              "ins Referenzband (deckt sich mit PHY034-038).")
    elif pair_seq:
        print("    - Y2-Kanal: Paar-Trend dreht/stagniert bei L=24..48 -"
              " Lage relativ zum Referenzband siehe Tabelle.")
    dips = [y4_dips[L]["T_dip"] for L in Ls if not y4_dips[L]["at_edge"]]
    if len(dips) >= 2 and all(b <= a for a, b in zip(dips, dips[1:])):
        print("    - Y4-Kanal: Dip wandert mit L nach unten "
              f"({', '.join(f'{d:.4f}' for d in dips)}) - konsistent mit "
              "Annaeherung an den Y4-beta-Kanal des Bands (0.6116).")
    else:
        print(f"    - Y4-Kanal: Dips {[y4_dips[L]['T_dip'] for L in Ls]} "
              f"(Rand-Flags {[y4_dips[L]['at_edge'] for L in Ls]}).")
    print(f"    - Sampler-Systematik (Walker, L=32): Y2-Spread {spread:.4f}; "
          f"T_BKT-Spread "
          f"{ {f'{a},{b}': round(v['spread'], 4) for (a, b), v in walker_pairs.items()} }")
    print("    - Grenzen: L=48 bleibt endlich; Walker-Systematik misst nur "
          "den Sampler, nicht den finite-size-Bias. Kein neuer Bestwert.")
    print(f"    - Laufzeit gesamt: {time.time() - t_start:.0f}s "
          f"({len(jobs)} WL-Jobs, max_workers={max_workers})")

    return {
        "module": "PHY042_honeycomb_wl_fss",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-07-06",
        "reference_band": {k: (list(v) if isinstance(v, tuple) else v)
                           for k, v in REF_BAND.items()},
        "Ls": list(Ls),
        "n_walkers_L32": n_walkers_L32,
        "master_seed": master_seed,
        "lnf_final": 1e-5,
        "windows": {str(L): windows[L] for L in Ls},
        "wl_sweeps": {f"{L}_{w}": results[(L, w)].wl_sweeps
                      for (L, w) in results},
        "covered": {f"{L}_{w}": [int(results[(L, w)].mask.sum()),
                                 int(len(results[(L, w)].mask))]
                    for (L, w) in results},
        "leak_max": leak_max,
        "validation_vs_wolff_L32": val_rows,
        "validation_vs_phy032_grid": grid_rows,
        "phy041_bridge": {"y4_dip_L24": y4_dip_24,
                          "committed": PHY041_Y4_DIP_L24, "ok": bridge_ok},
        "pair_tbkt": {f"{a}_{b}": v for (a, b), v in pair_tbkt.items()},
        "pair_trend_still_rising": bool(trend_up),
        "phy041_pairs": {f"{a}_{b}": v for (a, b), v in PHY041_PAIRS.items()},
        "y4_dips": {str(L): y4_dips[L] for L in Ls},
        "walker_y2_spread_L32": spread,
        "walker_pair_tbkt": {f"{a}_{b}": v
                             for (a, b), v in walker_pairs.items()},
        "curves": {f"{L}_{w}": {"T": curves[(L, w)]["T"].tolist(),
                                "y2": curves[(L, w)]["y2"].tolist(),
                                "y4_scaled": curves[(L, w)]["y4_scaled"].tolist()}
                   for (L, w) in results},
        "pass_gates": gates,
        "overall_pass": overall,
    }


def _clean(o):
    if o is None:
        return None
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, np.ndarray):
        return [_clean(x) for x in o.tolist()]
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    return o


if __name__ == "__main__":
    report = run_phy042()
    print("\n--- JSON-Report ---")
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
