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
3. MULTI-WALKER-SYSTEMATIK AN BEIDEN GROSSEN L (Haertung nach Lauf 1,
   2026-07-06): je 3 unabhaengige g(E)-Walker bei L=32 UND L=48
   (Walker 0: Standard-Stream 650+L; Walker w>=1: Stream 90000+1000*w+L,
   kollisionsfrei zu 400+/500+/600+/650+/660+/700+ inkl. +1000L-Termen).
   LAUF-1-BEFUND (Negativ-Result, Gate-Log results/): mit nur einem
   Walker je L>=32 sind die Y2-Kurven oberhalb T~0.60 sampler-limitiert
   (Walker-Spread bis 0.14, FSS-Ordnungsverletzung L32 vs L48, Y4-Dip
   nicht walker-robust). Konsequenz-Haertungen:
   a) Hauptanalyse-Kurve je L>=32 = WALKER-MITTEL (Bias ~1/sqrt(3));
   b) VALIDITAETS-DOMAENE je L aus dem Walker-Spread (< 0.04 =
      VAL-A-Y2-Toleranz; 0.02/0.01 als strengere Stufen ausgewiesen);
   c) Paare sind nur QUOTIERBAR, wenn ihr Crossing in beiden Domaenen
      liegt; sonst explizites NR;
   d) COVERAGE-MASSEN-GATE: kanonische Masse (volle lng) auf
      produktions-unbesetzten Bins < 1e-3 fuer jedes in-Domaene-T.
4. GETRENNTE KANAELE: Y2 (WM-Paare) und Y4 (Dip-Lage je L) werden separat
   ausgewertet und separat gegen das Referenzband eingeordnet; die
   Walker-Robustheit des Y4-Dips wird je L explizit geprueft.

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

import hashlib
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
    # Nachtrag 2026-09-26 (Issue #45 §3, spec/260926 ... provenance addendum):
    # DIREKT berichtete T-Kanaele. 0.576(3) war bis dahin gar nicht im Band
    # (PR #18 hatte es als "unattribuiert" gefuehrt); es ist der WL-Upsilon-
    # T-Wert derselben Quelle wie die beta-Kanaele. Die Journal-Fassung von
    # arXiv:2406.14812 (PTEP) berichtet andere Werte als arXiv v1 - beide
    # Fassungen bleiben getrennt stehen (keine stille Ersetzung).
    "upsilon_wl_T": (0.576, 0.003),
    "upsilon_sa_T": (0.575, 0.008),
    "upsilon4_wl_T": (0.568, 0.001),
    "upsilon4_sa_T": (0.551, 0.011),
    "helicity_ptep": (0.576, 0.004),
    "nn_ptep": (0.572, 0.003),
}

# Provenienz je Kanal. Historische Kanaele bleiben fuer Lineage im Dictionary,
# duerfen aber nicht als aktuelle Referenzfamilie interpretiert werden.
# Status 2026-09-27 wird durch spec/260926 ... provenance addendum supersediert.
REF_PROVENANCE = {
    "multi_lattice": ("arXiv:2501.07388 v1 = J.Phys.A 58 065003 (2025)",
                      "T, ohne Fehlerbalken", "current_method_anchor"),
    "helicity_direct": ("arXiv:2406.14812 v1", "T, Helicity",
                        "superseded_historical"),
    "nn_mc": ("arXiv:2406.14812 v1", "T, NN", "superseded_historical"),
    "upsilon_beta": ("arXiv:2406.12076 v3/v4 = Phys.Scr. 100 065953 (2025)",
                     "beta=1.687(3), Upsilon", "superseded_historical"),
    "upsilon4_beta": ("arXiv:2406.12076 v3/v4 = Phys.Scr. 100 065953 (2025)",
                      "beta=1.635(11), Upsilon_4", "superseded_historical"),
    "binder_beta": ("arXiv:2406.12076 v2..v4 Abstract",
                    "beta=1.724(2), Binder", "superseded_historical"),
    "upsilon_wl_T": ("arXiv:2406.12076 Abstract (v4) = Phys.Scr. 100 065953",
                     "T=0.576(3), Upsilon, Wang-Landau", "primary_text_verified"),
    "upsilon_sa_T": ("arXiv:2406.12076 Abstract (v4) = Phys.Scr. 100 065953",
                     "T=0.575(8), Upsilon, sim. annealing",
                     "primary_text_verified"),
    "upsilon4_wl_T": ("arXiv:2406.12076 Haupttext", "T=0.568(1), Upsilon_4, WL",
                      "primary_text_verified"),
    "upsilon4_sa_T": ("arXiv:2406.12076 Haupttext",
                      "T=0.551(11), Upsilon_4, SA", "primary_text_verified"),
    "helicity_ptep": ("PTEP 2024(10) 103A02 (Journal-Fassung arXiv:2406.14812)",
                      "T=0.576(4), Helicity", "primary_text_verified"),
    "nn_ptep": ("PTEP 2024(10) 103A02 (Journal-Fassung arXiv:2406.14812)",
                "T=0.572(3), NN", "primary_text_verified"),
}


REF_CURRENT_KEYS = (
    "multi_lattice",
    "upsilon_wl_T",
    "upsilon_sa_T",
    "upsilon4_wl_T",
    "upsilon4_sa_T",
    "helicity_ptep",
    "nn_ptep",
)

# PHY041-Bruecke: committed Reportwerte (results/260702 PHY041 ... report.txt)
PHY041_PAIRS = {(12, 16): 0.5951, (12, 24): 0.6029, (16, 24): 0.6087}
PHY041_Y4_DIP_L24 = 0.6500

# RNG-Stream-Vertrag PHY042: Walker 0 nutzt die PHY041-Streams (600er Anker,
# 650+L WL) fuer die deterministische L=24-Bruecke; Zusatz-Walker w>=1 nutzen
# 90000+1000*w (+L in wl_entropic_lattice) - kollisionsfrei zu allen
# bisherigen Vertraegen (300+/400+/500+/600+/650+/660+/700+ inkl. +1000L-
# Terme; Hinweis Code-Audit L1 2026-07-10: 300+ ist von PHY028 UND PHY039
# doppelt gebucht, dort nur durch disjunkte T-Gitter getrennt).
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


def _validity_domain(spread: np.ndarray, thr: float) -> np.ndarray:
    """Zusammenhaengende Validitaets-Domaene ab dem unteren T-Rand:
    True solange der Walker-Spread < thr bleibt; ab der ersten Verletzung
    False (Lauf-1-Befund: der Spread waechst monoton mit T - eine
    zusammenhaengende Domaene ist die ehrliche Lesart)."""
    mask = np.zeros(len(spread), dtype=bool)
    for k, s in enumerate(spread):
        if s >= thr:
            break
        mask[k] = True
    return mask


# Domaenen-Semantik (Issue #45 §2, 2026-09-26). Bis dahin setzte der
# Einzel-Walker-Pfad spread=0, Domaene="alles gueltig" und T_max=Gitterende.
# Der JSON-Report fuehrte L=24 dadurch als "domain_tmax_spread004": 0.67 in
# derselben Form wie die GEMESSENEN 0.60 (L32) / 0.585 (L48) - eine
# Nicht-Messung, die wie eine Messung aussah. Jetzt: None + Grund. Die
# Gate-Maske bleibt fuer ungemessene L "alles True" - das macht das
# Coverage-Gate STRENGER (bindet ueber das ganze Gitter), nicht laxer.
DOMAIN_THRESHOLD = 0.04  # VAL-A-Y2-Toleranz (Lauf-1-Haertung)
_UNMEASURED_REASON = (
    "n_walkers={n} < 2: walker spread undefined, no sampler validity "
    "domain measured (single-walker L; PHY041 bridge). Not a measurement.")


def _walker_domain(y2_curves: list, t_grid: np.ndarray,
                   thr: float = DOMAIN_THRESHOLD) -> dict:
    """Sampler-Validitaets-Domaene EINES L aus den Y2-Kurven seiner Walker.

    >= 2 Walker -> gemessen: spread = max-min je T, mask = _validity_domain,
    tmax = groesstes T der Domaene bzw. None, wenn die gemessene Domaene
    LEER ist (dann ist jedes Paar mit diesem L fail-closed NICHT quotierbar).
    < 2 Walker -> NICHT gemessen: spread = tmax = None, reason gesetzt;
    mask = alles True nur als konservative Gate-Maske.
    """
    n = len(y2_curves)
    if n < 2:
        return {"measured": False, "n_walkers": n, "spread": None,
                "mask": np.ones(len(t_grid), dtype=bool), "tmax": None,
                "reason": _UNMEASURED_REASON.format(n=n)}
    stack = np.stack([np.asarray(c, dtype=float) for c in y2_curves])
    spread = stack.max(axis=0) - stack.min(axis=0)
    mask = _validity_domain(spread, thr)
    tmax = float(t_grid[mask].max()) if mask.any() else None
    return {"measured": True, "n_walkers": n, "spread": spread, "mask": mask,
            "tmax": tmax,
            "reason": None if tmax is not None else
            f"measured: walker spread >= {thr} already at the lowest T "
            "(empty validity domain)"}


def _pair_domain_limit(dom_a: dict, dom_b: dict) -> tuple:
    """(T-Grenze, Basis) fuer die Quotierbarkeit eines Paar-Crossings.

    Fail-closed und reihenfolge-unabhaengig. Vorher `min(tmax_a, tmax_b)`
    mit NaN fuer eine leere Domaene: `min(0.6, nan)` = 0.6, aber
    `min(nan, 0.6)` = nan - je nach Paar-Reihenfolge haette eine LEERE
    Domaene das Paar quotierbar gelassen (latenter Fail-open, im Lauf
    2026-07-06 nicht bindend, weil keine Domaene leer war).

    Basis: "both_measured" | "partial" (ein L ungemessen - die Grenze
    stammt nur vom gemessenen L; die ungemessene Seite traegt keine
    Sampler-Evidenz) | "empty" (gemessene leere Domaene -> None) |
    "unmeasured" (keine Seite gemessen -> None: ohne Sampler-Evidenz
    kein quotierbares Paar).
    """
    for d in (dom_a, dom_b):
        if d["measured"] and d["tmax"] is None:
            return None, "empty"
    measured = [d["tmax"] for d in (dom_a, dom_b) if d["measured"]]
    if not measured:
        return None, "unmeasured"
    return min(measured), ("both_measured" if len(measured) == 2
                           else "partial")


def _domain_status(d: dict) -> dict:
    """JSON-taugliche Sicht auf _walker_domain (ohne Arrays)."""
    return {"measured": d["measured"], "n_walkers": d["n_walkers"],
            "tmax": d["tmax"],
            "max_spread": (None if d["spread"] is None
                           else float(np.max(d["spread"]))),
            "reason": d["reason"]}


def reanalyse_domains(report: dict, thr: float = DOMAIN_THRESHOLD) -> dict:
    """Domaenen + Paar-Quotierbarkeit aus einem COMMITTETEN PHY042-Report
    neu ableiten (MC-frei: nutzt nur die gespeicherten Walker-Kurven).

    Zweck (Issue #45 §2): der Report 260707 ist per voller SHA-256 gepinnt
    und wird NICHT umgeschrieben (AGENTS.md: Lineage ehrlich). Diese
    Funktion liefert die korrigierte Domaenen-Semantik als separates,
    reproduzierbares Erratum und prueft zugleich, dass die gemessenen
    Domaenen (L32/L48) und die Quotierbarkeits-Urteile aus den Kurven
    folgen.
    """
    Ls = [int(x) for x in report["Ls"]]
    t_grid = None
    dom = {}
    for L in Ls:
        nw = int(report["n_walkers"][str(L)])
        cs = [report["curves"][f"{L}_{w}"] for w in range(nw)]
        t_grid = np.asarray(cs[0]["T"], dtype=float)
        dom[L] = _walker_domain([c["y2"] for c in cs], t_grid, thr)
    pairs = {}
    for key, tb in report["pair_tbkt_mean_curves"].items():
        a, b = (int(x) for x in key.split("_"))
        q, limit, basis = _pair_quotable(tb, dom[a], dom[b])
        pairs[key] = {"tbkt": tb, "quotable": q, "limit_tmax": limit,
                      "basis": basis,
                      "quotable_as_committed": bool(
                          report["pair_quotable"][key])}
    return {"threshold": thr,
            "domain_status": {str(L): _domain_status(dom[L]) for L in Ls},
            "committed_domain_tmax_spread004":
                report["domain_tmax_spread004"],
            "pairs": pairs}


def _pair_quotable(tb, dom_a: dict, dom_b: dict) -> tuple:
    """(quotierbar, T-Grenze, Basis): Crossing existiert UND liegt in der
    fail-closed Paar-Grenze (_pair_domain_limit)."""
    limit, basis = _pair_domain_limit(dom_a, dom_b)
    q = tb is not None and limit is not None and tb <= limit
    return bool(q), limit, basis


def _pair_inside_bounds(
    tb: float | None,
    bounds_a: tuple[float, float] | None,
    bounds_b: tuple[float, float] | None,
) -> bool:
    """Require a crossing to lie inside BOTH evidenced [T_min,T_max] domains."""
    if tb is None or bounds_a is None or bounds_b is None:
        return False
    lo = max(float(bounds_a[0]), float(bounds_b[0]))
    hi = min(float(bounds_a[1]), float(bounds_b[1]))
    return bool(lo <= float(tb) <= hi)


def _walker_plan(Ls, n_walkers: int, min_walkers: int | None = None) -> dict:
    """Walker je L: n_walkers fuer L>=32, 1 fuer L<32 (PHY041-Bruecke).

    min_walkers (W4-Vertrag, Issue #45 §2: >= 3 Walker an JEDEM L eines
    W4-Laufs) prueft fail-closed VOR jeder Rechnung - ein Einzel-Walker-L
    darf in einem W4-Lauf nicht still mitlaufen.
    """
    walkers = {L: (n_walkers if L >= 32 else 1) for L in Ls}
    if min_walkers is not None:
        short = {L: w for L, w in walkers.items() if w < min_walkers}
        if short:
            raise ValueError(
                f"Walker-Vertrag verletzt: min_walkers={min_walkers}, "
                f"aber {short} (L -> Walker). L<32 laeuft in PHY042 als "
                "Einzel-Walker-Bruecke; ein W4-Lauf braucht einen eigenen "
                "Walker-Plan.")
    return walkers


def _uncovered_mass(res, T: float) -> float:
    """Kanonische Gewichtsmasse auf produktions-UNBESETZTEN Bins, mit der
    VOLLEN lng (alle Bins; die WL-Phase hat alle Bins besucht). Verallgemeinert
    das Rand-Leak-Gate: misst je T, wieviel kanonisches Gewicht auf Bins
    ohne mikrokanonische Aggregate faellt (dort wuerde die maskierte
    Renormierung die Kurve verfaelschen)."""
    beta = 1.0 / T
    x = res.lng - beta * res.centers
    x = x - x.max()
    w = np.exp(x)
    w /= w.sum()
    return float(w[~res.mask].sum())


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


def run_phy042(Ls=(24, 32, 48), n_walkers=3, master_seed=42,
               max_workers=4, min_walkers: int | None = None) -> dict:
    t_grid = _T_GRID
    # Walker-Plan zuerst: ein verletzter min_walkers-Vertrag muss VOR den
    # teuren Anker-/WL-Laeufen scheitern (fail-closed, Issue #45 §2).
    walkers = _walker_plan(Ls, n_walkers, min_walkers)
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
    # Multi-Walker an BEIDEN grossen L (Lauf-1-Haertung): der Sampler-
    # Systematik-Check muss jede Gitter-Groesse abdecken, deren Kurven in
    # Paar-Schaetzer eingehen. L=24 bleibt Einzel-Walker (PHY041-Bruecke);
    # Plan oben via _walker_plan (fail-closed gegen min_walkers).
    jobs = []
    for L in Ls:
        e_lo, e_hi = windows[L]["window_ps"]
        for w in range(walkers[L]):
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
    # Hauptanalyse-Kurve je L: Walker-MITTEL (L>=32) bzw. Walker 0 (L=24).
    # Das Mittel ueber unabhaengige g(E)-Walker reduziert den Lauf-1-Befund
    # (walker-abhaengiger Sampler-Bias) um ~1/sqrt(n_walkers); der Spread
    # bleibt als ehrlicher Systematik-Schaetzer ausgewiesen.
    main_curve = {}
    for L in Ls:
        ws = range(walkers[L])
        main_curve[L] = {
            "T": t_grid,
            "y2": np.mean([curves[(L, w)]["y2"] for w in ws], axis=0),
            "y4_scaled": np.mean([curves[(L, w)]["y4_scaled"] for w in ws],
                                 axis=0),
            "E": np.mean([curves[(L, w)]["E"] for w in ws], axis=0),
        }

    # --- Leak-Gate: kanonisches Gewicht an den Fensterraendern ------------
    leak_max = 0.0
    for key, res in results.items():
        for T in t_grid:
            leak_max = max(leak_max, canonical_edge_leak(res, float(T)))
    leak_ok = leak_max < 1e-3
    print(f"\n[LEAK] max. Randgewicht ueber alle (L, walker, T): "
          f"{leak_max:.2e} (Gate < 1e-3) -> {'PASS' if leak_ok else 'FAIL'}")

    # --- Validitaets-Domaenen aus Walker-Spread (Lauf-1-Haertung) ----------
    # Der Spread unabhaengiger Walker diagnostiziert, ab welchem T die
    # Y2-Kurven sampler-limitiert sind. Schwelle 0.04 = VAL-A-Y2-Toleranz;
    # 0.02/0.01 werden als strengere Domaenen mit ausgewiesen.
    print("\n[DOMAIN] Validitaets-Domaenen aus Walker-Spread (L>=32):")
    dom: dict[int, dict] = {}
    domains: dict[int, np.ndarray] = {}   # Gate-Maske (ungemessen: alles)
    for L in Ls:
        dom[L] = _walker_domain([curves[(L, w)]["y2"]
                                 for w in range(walkers[L])], t_grid)
        domains[L] = dom[L]["mask"]
        if not dom[L]["measured"]:
            print(f"      L={L}: NICHT GEMESSEN ({walkers[L]} Walker) - "
                  "Report fuehrt T_max=null + Grund; Sampler-Guard fuer "
                  "L=24 ist VAL-B (PHY032-Gitter bis T=0.6475).")
            continue
        sp = dom[L]["spread"]
        strict = {thr: (float(t_grid[_validity_domain(sp, thr)].max())
                        if _validity_domain(sp, thr).any() else None)
                  for thr in (0.01, 0.02)}
        tmax_txt = ("LEER" if dom[L]["tmax"] is None
                    else f"{dom[L]['tmax']:.4f}")
        print(f"      L={L}: T<= {tmax_txt} (Spread<{DOMAIN_THRESHOLD}); "
              f"strenger: <0.02 -> T<={strict[0.02]}, <0.01 -> T<={strict[0.01]}; "
              f"max Spread {sp.max():.4f}")

    # --- Coverage-Massen-Gate (nur innerhalb der Domaene bindend) ----------
    unc_max = 0.0
    for (L, w), res in results.items():
        for k, T in enumerate(t_grid):
            if domains[L][k]:
                unc_max = max(unc_max, _uncovered_mass(res, float(T)))
    unc_ok = unc_max < 1e-3
    print(f"\n[COVER] max. kanonische Masse auf unbesetzten Bins "
          f"(in-Domaene): {unc_max:.2e} (Gate < 1e-3) "
          f"-> {'PASS' if unc_ok else 'FAIL'}")

    # --- VAL-A: frische Wolff-Referenz L=32; Y2-Gate nur in-Domaene --------
    print("\n[VAL-A] WL(Mittel-Kurve) vs frische Wolff-Referenz (L=32):")
    e_ok, y2_ok = True, True
    y2_any_in_domain = False  # Haertung 2026-07-10 (Code-Audit M2)
    val_rows = []
    for T in (0.55, 0.60, 0.65):
        c = main_curve[32]
        k = int(np.argmin(np.abs(c["T"] - T)))
        in_dom = bool(domains[32][k])
        Ewl = c["E"][k] / results[(32, 0)].n
        y2wl = c["y2"][k]
        ref = wolff_reference(32, T, master_seed=master_seed)
        de = abs(Ewl - ref["E_ps"])
        dy = abs(y2wl - ref["y2"])
        e_ok = e_ok and de < 0.03
        if in_dom:
            y2_ok = y2_ok and dy < 0.04
            y2_any_in_domain = True
        val_rows.append({"T": T, "in_domain": in_dom, "E_wl": Ewl,
                         "E_wolff": ref["E_ps"], "y2_wl": y2wl,
                         "y2_wolff": ref["y2"], "y2_wolff_sem": ref["y2_sem"]})
        tag = "" if in_dom else "  [AUSSERHALB DOMAENE - nur Evidenz, kein Gate]"
        print(f"      T={T:.3f}  <E>/N WL={Ewl:+.4f} Wolff={ref['E_ps']:+.4f} "
              f"(d={de:.4f})  Y2 WL={y2wl:.4f} "
              f"Wolff={ref['y2']:.4f}+-{ref['y2_sem']:.4f} (d={dy:.4f}){tag}")

    # --- VAL-B: Drift-Guard gegen committed PHY032-Gitter (L=24) ----------
    print("\n[VAL-B] WL vs committed PHY032-Messgitter (L=24, Drift-Guard):")
    grid032 = parse_phy032_grid()
    # Haertung 2026-07-10 (Code-Audit M1): der Drift-Guard darf nicht
    # vacuous passen, wenn der Parser nichts findet (Vertrag: 8 T-Punkte
    # fuer L=24 im committeten PHY032-Report).
    grid_ok = len(grid032.get(24, [])) == 8
    if not grid_ok:
        print("      PARSE-VERTRAG VERLETZT: erwarte 8 Gitterpunkte fuer "
              f"L=24, gefunden {len(grid032.get(24, []))}")
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

    # Effektive Evidenzdomaene fuer NEUE Paarurteile: beide Grenzen binden.
    # Gemessene Walker-Domaenen beginnen am unteren T-Gitterrand. Fuer die
    # historische L=24-Einzel-Walker-Bruecke darf nur der separat validierte
    # PHY032-Drift-Guard als Fallback dienen; er wird nicht als Walker-Spread
    # umetikettiert.
    effective_domain_bounds: dict[int, tuple[float, float] | None] = {}
    effective_domain_basis: dict[int, str] = {}
    for L in Ls:
        if dom[L]["measured"] and np.any(dom[L]["mask"]):
            vt = t_grid[dom[L]["mask"]]
            effective_domain_bounds[L] = (float(vt.min()), float(vt.max()))
            effective_domain_basis[L] = "walker_spread"
        else:
            effective_domain_bounds[L] = None
            effective_domain_basis[L] = (
                "measured_empty" if dom[L]["measured"] else "unmeasured"
            )

    if 24 in Ls and effective_domain_bounds.get(24) is None and grid_rows:
        if grid_ok and all(bool(row["ok"]) for row in grid_rows):
            effective_domain_bounds[24] = (
                min(float(row["T"]) for row in grid_rows),
                max(float(row["T"]) for row in grid_rows),
            )
            effective_domain_basis[24] = "PHY032_drift_guard"

    # --- VAL-C: PHY041-Bruecke (bitgleiche L=24-Reproduktion) --------------
    c24 = main_curve[24]
    k24 = int(np.argmin(c24["y4_scaled"]))
    y4_dip_24 = float(c24["T"][k24])
    bridge_ok = abs(y4_dip_24 - PHY041_Y4_DIP_L24) < 1e-9
    print(f"\n[VAL-C] PHY041-Bruecke: Y4-Dip L=24 = {y4_dip_24:.4f} "
          f"(committed {PHY041_Y4_DIP_L24:.4f}) "
          f"-> {'PASS' if bridge_ok else 'FAIL'}")

    # --- Y2-Kanal: C-eliminierte Paare auf Mittel-Kurven + Quotierbarkeit --
    print("\n[Y2-FSS] C-eliminierte WM-Paare (Walker-Mittel-Kurven):")
    pairs = [(Ls[i], Ls[j]) for i in range(len(Ls)) for j in range(len(Ls))
             if Ls[j] > Ls[i]]
    pair_tbkt = {}
    pair_quotable = {}
    pair_basis = {}
    for La, Lb in pairs:
        tb = tbkt_pair_from_curves(t_grid, main_curve[La]["y2"], La,
                                   main_curve[Lb]["y2"], Lb)
        pair_tbkt[(La, Lb)] = tb
        bounds_a = effective_domain_bounds.get(La)
        bounds_b = effective_domain_bounds.get(Lb)
        q = _pair_inside_bounds(tb, bounds_a, bounds_b)
        pair_quotable[(La, Lb)] = q
        pair_basis[(La, Lb)] = {
            "bounds_a": bounds_a,
            "bounds_b": bounds_b,
            "basis_a": effective_domain_basis.get(La),
            "basis_b": effective_domain_basis.get(Lb),
        }
        if tb is None:
            print(f"      T_BKT({La},{Lb}): kein Nulldurchgang im T-Fenster")
        else:
            tag = (
                "QUOTIERBAR"
                if q
                else "NR: Crossing ausserhalb oder ohne zwei belegte "
                     "[T_min,T_max]-Domaenen"
            )
            tag += (
                f", Basis {effective_domain_basis.get(La)}/"
                f"{effective_domain_basis.get(Lb)}"
            )
            print(f"      T_BKT({La},{Lb}) = {tb:.4f}  "
                  f"(vs Multi-Lattice 0.573: "
                  f"{(tb - 0.573) / 0.573 * 100:+.2f}%)  [{tag}]")
    print(f"      PHY041-Paare (L<=24, Einzel-Walker): "
          f"{', '.join(f'{v:.4f}' for v in PHY041_PAIRS.values())}")

    # Sampler-Systematik der Paare: alle Walker-Kombinationen
    print("\n[WALKER] Paar-T_BKT ueber Walker-Kombinationen:")
    walker_pairs = {}
    for (La, Lb) in pairs:
        vals = []
        for wa in range(walkers[La]):
            for wb in range(walkers[Lb]):
                tb = tbkt_pair_from_curves(
                    t_grid, curves[(La, wa)]["y2"], La,
                    curves[(Lb, wb)]["y2"], Lb)
                if tb is not None:
                    vals.append(tb)
        if vals:
            walker_pairs[(La, Lb)] = {
                "values": vals, "spread": max(vals) - min(vals),
                "n_combos": walkers[La] * walkers[Lb]}
            print(f"      T_BKT({La},{Lb}): {len(vals)} Kombos, "
                  f"[{min(vals):.4f}, {max(vals):.4f}] "
                  f"(Spread {max(vals) - min(vals):.4f})")

    # --- Y4-Kanal: Dip-Lage je L (Mittel-Kurve + Walker-Robustheit) --------
    print("\n[Y4-FSS] N*Upsilon_4-Dip je L:")
    y4_dips = {}
    for L in Ls:
        c = main_curve[L]
        k = int(np.argmin(c["y4_scaled"]))
        wd = []
        for w in range(walkers[L]):
            cw = curves[(L, w)]
            wd.append(float(cw["T"][int(np.argmin(cw["y4_scaled"]))]))
        robust = (max(wd) - min(wd)) <= 0.005 + 1e-12
        y4_dips[L] = {"T_dip": float(c["T"][k]),
                      "depth": float(c["y4_scaled"][k]),
                      "at_edge": bool(k == 0 or k == len(c["T"]) - 1),
                      "walker_dips": wd, "walker_robust": bool(robust)}
        print(f"      L={L}: Mittel-Kurven-Dip T={y4_dips[L]['T_dip']:.4f}"
              f"{', am Gitterrand' if y4_dips[L]['at_edge'] else ''}; "
              f"Walker-Dips {wd} -> "
              f"{'robust' if robust else 'NICHT walker-robust'}")

    # --- Glaette-Gate (L=48 Mittel-Kurve, Kernfenster) ----------------------
    core = (t_grid >= _CORE_WIN[0]) & (t_grid <= _CORE_WIN[1])
    c48 = main_curve[48]["y2"][core]
    y2_smooth = bool(np.all(np.diff(c48) < 0))
    print(f"\n[Y2] Upsilon_2(T) L=48 (Mittel) Kernfenster {list(_CORE_WIN)}: "
          f"streng monoton fallend={y2_smooth} "
          f"({c48[0]:.3f} .. {c48[-1]:.3f})")

    gates = {
        "PASS_ALIGNED_EXACT_THREE_QUARTERS": aligned_ok,
        "PASS_NO_CANONICAL_EDGE_LEAK": leak_ok,
        "PASS_NO_UNCOVERED_MASS_IN_DOMAIN": unc_ok,
        "PASS_WL_ENERGY_MATCHES_WOLFF_L32": e_ok,
        # Haertung 2026-07-10 (Code-Audit M2): das Gate darf nicht vacuous
        # passen, wenn KEIN Validierungs-T in der Domaene liegt - genau dann
        # waere der Sampler am kaputtesten.
        "PASS_WL_Y2_MATCHES_WOLFF_L32_IN_DOMAIN": y2_ok and y2_any_in_domain,
        "PASS_WL_Y2_MATCHES_PHY032_GRID_L24": grid_ok,
        "PASS_PHY041_BRIDGE_Y4_DIP_L24": bridge_ok,
        "PASS_Y2_SMOOTH_NOISE_FREE_L48": y2_smooth,
        "PASS_PAIR_2432_WITHIN_VALIDITY": bool(pair_quotable.get((24, 32))),
    }
    print("\n[GATES]")
    for k, v in gates.items():
        print(f"      [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(gates.values())
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}  "
          f"(= Pipeline-Integritaet innerhalb der selbst-diagnostizierten "
          f"Validitaets-Domaene; Physik-Befund ist FINDING)")

    # --- ehrlicher Physik-Befund + Negativ-Results --------------------------
    print("\n[FINDING] (ehrlich, NICHT getunt)")
    tb2432 = pair_tbkt.get((24, 32))
    if tb2432 is not None and pair_quotable.get((24, 32)):
        sy = walker_pairs.get((24, 32), {}).get("spread", float("nan"))
        print(f"    - QUOTIERBAR: T_BKT(24,32) = {tb2432:.4f} "
              f"+- {sy / 2:.4f} (Sampler-Systematik = halber Walker-Spread); "
              f"Lage im Referenzband: nahe Upsilon-beta-Kanal 0.5928, "
              f"unterhalb PHY041 (16,24)=0.6087.")
    # Haertung 2026-07-10 (Code-Audit L2): default schuetzt den
    # n_walkers=1-Pfad (leerer Generator -> ValueError VOR dem Report).
    max_spread = max((dom[L]["spread"].max() for L in Ls
                      if dom[L]["measured"]),
                     default=float("nan"))
    print("    - NR-PHY042-01: Y2 ist bei L>=32 oberhalb T~0.60 "
          "sampler-limitiert (Walker-Spread bis "
          f"{max_spread:.2f}); "
          "Validitaets-Domaenen oben ausgewiesen.")
    nr02 = [f"({a},{b})" for (a, b) in pairs
            if pair_tbkt.get((a, b)) is not None
            and not pair_quotable.get((a, b))]
    if nr02:
        print(f"    - NR-PHY042-02: Paare {', '.join(nr02)} crossen "
              "ausserhalb der Domaene - bei diesem Statistik-Budget "
              "NICHT belastbar (kein Ersatz durch mehr Walker; braucht "
              "laengere Produktion/bessere Dekorrelation).")
    nr03 = [str(L) for L in Ls if walkers[L] > 1
            and not y4_dips[L]["walker_robust"]]
    if nr03:
        print(f"    - NR-PHY042-03: Y4-Dip bei L in {{{', '.join(nr03)}}} "
              "nicht walker-robust - Y4-FSS oberhalb L=24 braucht mehr "
              "Produktion; die Einzel-Walker-Dip-Wanderung aus Lauf 1 "
              "(0.65/0.64/0.615) ist als Trend NICHT quotierbar.")
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
        "n_walkers": {str(L): walkers[L] for L in Ls},
        "master_seed": master_seed,
        "lnf_final": 1e-5,
        "windows": {str(L): windows[L] for L in Ls},
        "wl_sweeps": {f"{L}_{w}": results[(L, w)].wl_sweeps
                      for (L, w) in results},
        "covered": {f"{L}_{w}": [int(results[(L, w)].mask.sum()),
                                 int(len(results[(L, w)].mask))]
                    for (L, w) in results},
        "leak_max": leak_max,
        "uncovered_mass_max_in_domain": unc_max,
        # Issue #45 §2: ungemessene L -> null + Grund (nicht 0 / Gitterende)
        "walker_spread": {str(L): (None if dom[L]["spread"] is None
                                   else dom[L]["spread"].tolist())
                          for L in Ls},
        "domain_tmax_spread004": {str(L): dom[L]["tmax"] for L in Ls},
        "domain_status": {str(L): _domain_status(dom[L]) for L in Ls},
        "effective_domain_bounds": {
            str(L): (
                None if effective_domain_bounds.get(L) is None
                else list(effective_domain_bounds[L])
            )
            for L in Ls
        },
        "effective_domain_basis": {
            str(L): effective_domain_basis.get(L) for L in Ls
        },
        "pair_domain_basis": {f"{a}_{b}": v
                              for (a, b), v in pair_basis.items()},
        "validation_vs_wolff_L32": val_rows,
        "validation_vs_phy032_grid": grid_rows,
        "phy041_bridge": {"y4_dip_L24": y4_dip_24,
                          "committed": PHY041_Y4_DIP_L24, "ok": bridge_ok},
        "pair_tbkt_mean_curves": {f"{a}_{b}": v
                                  for (a, b), v in pair_tbkt.items()},
        "pair_quotable": {f"{a}_{b}": bool(v)
                          for (a, b), v in pair_quotable.items()},
        "walker_pair_tbkt": {f"{a}_{b}": v
                             for (a, b), v in walker_pairs.items()},
        "phy041_pairs": {f"{a}_{b}": v for (a, b), v in PHY041_PAIRS.items()},
        "y4_dips": {str(L): y4_dips[L] for L in Ls},
        "curves": {f"{L}_{w}": {"T": curves[(L, w)]["T"].tolist(),
                                "y2": curves[(L, w)]["y2"].tolist(),
                                "y4_scaled": curves[(L, w)]["y4_scaled"].tolist()}
                   for (L, w) in results},
        "main_curves": {str(L): {"T": main_curve[L]["T"].tolist(),
                                 "y2": main_curve[L]["y2"].tolist(),
                                 "y4_scaled": main_curve[L]["y4_scaled"].tolist()}
                        for L in Ls},
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


PHY042_REPORT_V01 = ("results/260707 PHY042 honeycomb wl-fss L24-32-48 "
                     "gate report.json")


def domain_erratum(report_path: Path = _ROOT / PHY042_REPORT_V01) -> dict:
    """Maschinenlesbares Erratum zum committeten PHY042-Report (Issue #45
    §2): korrigierte Domaenen-Semantik, abgeleitet via reanalyse_domains.
    Der Quell-Report bleibt unveraendert; seine SHA-256 (ueber LF-normierte
    Bytes = committete Bytes, plattform-portabel) bindet das Erratum."""
    raw = Path(report_path).read_bytes().replace(b"\r\n", b"\n")
    rep = json.loads(raw)
    return {
        "module": "PHY042_domain_semantics_erratum",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-09-26",
        "issue": "#45 section 2",
        "source_report": PHY042_REPORT_V01,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "note": ("Committed report lists domain_tmax_spread004['24'] = 0.67 "
                 "and walker_spread['24'] = 0.0 x31. Both are construction "
                 "artefacts of the single-walker path (grid end / zero), "
                 "not measurements. Corrected semantics: null + reason. "
                 "Measured domains (L32, L48) and all pair-quotability "
                 "verdicts are unchanged (re-derived from stored curves)."),
        **reanalyse_domains(rep),
    }


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--reanalyse":
        src = Path(sys.argv[2]) if len(sys.argv) > 2 else (
            _ROOT / PHY042_REPORT_V01)
        print(json.dumps(_clean(domain_erratum(src)), indent=2,
                         allow_nan=False))
    else:
        report = run_phy042()
        print("\n--- JSON-Report ---")
        print(json.dumps(_clean(report), indent=2, allow_nan=False))
